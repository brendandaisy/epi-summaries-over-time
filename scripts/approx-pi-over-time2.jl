using MarginalDivergence
using Distributions, MonteCarloMeasurements
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using NamedTupleTools
using Optim
using CSV, DataFrames

include("../src/transform-funcs.jl")

function load_sir_from_samples(varname, M)
    df = CSV.read("results/cond-samples/samp-$(varname).csv", DataFrame)

    ret = SIRModel{Float32}(
        S₀=Particles(Float32.(collect(df.S0))), 
        β=Particles(Float32.(collect(df.beta))), 
        α=Particles(Float32.(collect(df.alpha)))
    )
    nparticles(ret.S₀) != M ? resample(ret, M) : ret
end

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
lm_true = SIRModel{Float32}(;θtrue...)
prob = de_problem(lm_true) # for the approximation methods

nrep = 1
T = 24
tss = [0.25, 0.2]
Ns = [3000]
Ms = [240_000]
# tss = [0.75, 1/3, 0.25]
# Ns = [2000, 2000, 2000]
# Ms = [150_000, 200_000, 200_000]

obs_mod = PoissonRate(1000)

function solve_adj(θ, ts)
    _prob = remake(prob, p=[θ[2], θ[1]], u0=[θ[3], 0.01])
    solve(_prob, Tsit5(); save_idxs=2, saveat=0:ts:T, reltol = 1e-6, abstol = 1e-6).u
end

function info_mat(ts)
    inf = solve(lm_true; save_idxs=2, saveat=0:ts:T).u
    O = Diagonal(1000 ./ inf)
    J = ForwardDiff.jacobian(p->solve_adj(p, ts), vcat(values(θtrue)...))
    return J' * O * J
end

F = [0.6 -1 1.25; 1 0 0; 0 0 1]

map(0.2:0.1:2) do ts
    JOJ = info_mat(ts)
    mdS0_approx = 0.5 * (log(nrep) - log(inv(JOJ)[3, 3]) - log(2*π)) - logpdf(θprior.S₀, θtrue.S₀)

    Iu = F' * JOJ * F
    # TODO wouldnt index be 1???
    mdgrate_approx = 0.5 * (log(nrep) - log(inv(Iu)[2, 2]) - log(2*π)) - log(0.5)
    (;mdS0_approx, mdgrate_approx)
end

ret = []
for i in eachindex(tss)
    ts = tss[i]
    N = Ns[i]
    M = Ms[i]

    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_true = solve(lm_true; save_idxs=2, saveat=0:ts:T).u

    # "Exact" calculation using Monte Carlo
    inf_pri = solve(lat_mod; save_idxs=2, saveat=0:ts:T).u
    y = observe_dist(obs_mod, inf_true, nrep)

    # inf_cond_S0 = solve(lat_mod, (S₀=θtrue.S₀,); save_idxs=2, saveat=0:ts:T).u
    # global mdS0 = marginal_divergence(y, inf_cond_S0, inf_pri, obs_mod; N=N)

    sir_cond_grate = load_sir_from_samples("growth-rate", M)
    inf_cond_grate = solve(sir_cond_grate; save_idxs=2, saveat=0:ts:T).u
    global mdgrate = marginal_divergence(y, inf_cond_grate, inf_pri, obs_mod; N=N)

    res = (;mdS0, mdgrate)
    println(res)
    push!(ret, res)
end