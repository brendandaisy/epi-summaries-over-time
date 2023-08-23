#=
pi-spread-obs.jl
Produce the results from Figure 2 of the main text, where a fixed number of observations
are placed evenly over different time spans, and both MC and asymptotic approx. are used
=#

using MarginalDivergence
using ConditionalTransform
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using Optim
using CSV, DataFrames

include("../src/transform-funcs.jl")

function sir_from_samples(α, β, S₀)
    SIRModel{Float32}(
        stop=Float32(maximum(tmax_set)),
        S₀=Particles(Float32.(S₀)), 
        β=Particles(Float32.(β)), 
        α=Particles(Float32.(α))
    )
end

"""
For the given `θtrue`, run the SIR simiulations while holding u = f(`θtrue`) fixed for each u

Returns a `Dict` of infection curves for each u
"""
function solve_inf_cond(θtrue, θprior, lat_mod)
    nsamples = nparticles(lat_mod.β)
    inf_conds = Dict()

    inf_conds["beta"] = solve(lat_mod, (β=θtrue.β,); save_idxs=2, dense=true, abstol=1f-12) # simulate x∣θᵢ

    #= Basic Reproductive Number =#
    Rtrue = rnot(θtrue.α, θtrue.β)
    αcond, βcond = sample_cond_f(NamedTupleTools.select(θprior, (:α, :β)), Rtrue, inv_rnot, d_inv_rnot; pivot=:β, nsamples)
    sir_cond = sir_from_samples(αcond, βcond, rand(θprior.S₀, nsamples))
    inf_conds["rep-number"] = solve(sir_cond; save_idxs=2, dense=true, abstol=1f-12)

    #= Initial growth rate =#
    grate = growth_rate(θtrue...)
    αcond, βcond, Scond = sample_cond_f(
        θprior, grate, inv_growth_rate, d_inv_growth_rate; 
        pivot=:α, nsamples
    )
    sir_cond = sir_from_samples(αcond, βcond, Scond)
    inf_conds["growth-rate"] = solve(sir_cond; save_idxs=2, dense=true, abstol=1f-12)
    
    return inf_conds
end

function ident_spread_obs(;N=100, M=60_000)
    # simulate from the full prior
    lat_mod = SIRModel{Float32}(
        stop=Float32(maximum(tmax_set)),
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    sol_pri = solve(lat_mod; save_idxs=2, dense=true, abstol=1f-12, reltol=1f-9)
    # simulate true and conditionals
    sol_true = solve(lat_mod, θtrue; save_idxs=2, dense=true)
    sol_conds = solve_inf_cond(θtrue, θprior, lat_mod)

    ret = []
    for T in tmax_set
        obs_t = Float32.(range(0, T, num_obs))
        y = observe_dist(obs_mod, sol_true(obs_t).u)
        inf_pri = sol_pri(obs_t).u
        
        for lab in keys(sol_conds)
            inf_cond = sol_conds[lab](obs_t).u
            δ = marginal_divergence(y, inf_cond, inf_pri, obs_mod; N=N)
            @show (T, δ)
            push!(ret, (lab=lab, max_t=T, num_obs=num_obs, md=δ))
        end
    end
    return ret
end

function solve_adj(θ, obs_t)
    _prob = remake(prob, p=[θ[2], θ[1]], u0=[θ[3], 0.01])
    solve(_prob, Tsit5(); save_idxs=2, saveat=obs_t, reltol=1e-6, abstol=1e-10).u
end

"""
Return the "JOJ" portion of eq. 16
"""
function info_mat(obs_t)
    inf = solve(lm_true; save_idxs=2, saveat=obs_t).u
    O = Diagonal(1000 ./ inf)
    J = ForwardDiff.jacobian(p->solve_adj(p, obs_t), vcat(values(θtrue)...))
    return J' * O * J
end

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
num_obs = 40
tmax_set = 1:2:41
obs_mod = PoissonRate(1000)

res_monte_carlo = ident_spread_obs(N=2000, M=70_000)

Frnot = [1 0 0; rnot(θtrue.α, θtrue.β) 0.2 0; 0 0 1] # F matrix for eq. 16
Fgrate = [-1 0.6 1.25; 0 1 0; 0 0 1]


lm_true = SIRModel{Float32}(;stop=Float32(maximum(tmax_set)), θtrue...)
prob = de_problem(lm_true)

res_approx = []
for T in tmax_set
    obs_t = range(0, T, num_obs)
    JOJ = info_mat(obs_t)
    δβ = -0.5 * (log(inv(JOJ)[2, 2]) + log(2*π)) - logpdf(θprior.β, θtrue.β)
    δrnot = -0.5 * (log(inv(Frnot' * JOJ * Frnot)[2, 2]) + log(2*π)) - log(0.03)
    δgrate = -0.5 * (log(inv(Fgrate' * JOJ * Fgrate)[1, 1]) + log(2*π)) - log(0.5)
    append!(res_approx, [
        (lab="beta", max_t=T, num_obs=num_obs, md=δβ),
        (lab="rep-number", max_t=T, num_obs=num_obs, md=δrnot),
        (lab="growth-rate", max_t=T, num_obs=num_obs, md=δgrate)
    ])
end
        