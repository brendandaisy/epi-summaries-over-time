#====
pi-over-time.jl
Run the experiment shown in main panel of Figure 1B, where PI is assessed with an increasing observation window
====#

using MarginalDivergence
using ConditionalTransform
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using UnPack
using CSV, DataFrames

include("../src/transform-funcs.jl")

"""
Load previously sampled θ from P(θ|`varname`) and construct a corresponding `SIRModel`

Assumes the previously calculated samples are correct for this setting of θtrue!!!
"""
function sir_from_samples(varname, M)
    df = CSV.read("results/cond-samples/samp-$(varname).csv", DataFrame)

    ret = SIRModel{Float32}(
        S₀=Particles(Float32.(collect(df.S0))), 
        β=Particles(Float32.(collect(df.beta))), 
        α=Particles(Float32.(collect(df.alpha)))
    )
    nparticles(ret.S₀) > M ? resample(ret, M) : ret
end

"""
Container to hold the conditions for this experiment
"""
struct Exper1Config
    θtrue
    θprior
    obs_mod
    inf_pri
    inf_conds
    inf_true
    y
end

"""
Create an `Exper1Config` by simulating necessary SIR processes based on the input `θtrue`, `θprior`, and `obs_mod`.

The result is to be based to `ident_over_tspan`.
"""
function exper1_setup(θtrue, θprior, obs_mod; M=60_000)
    # setup model and simulate infections from prior:
    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_pri = solve(lat_mod; save_idxs=2, saveat=1).u 
    # solve true infection curve and draw samples y|θtrue:
    inf_true = solve(lat_mod, θtrue; save_idxs=2, saveat=1).u
    y = observe_dist(obs_mod, inf_true)

    # save a large bank of infection curves, conditional on each variable u:
    inf_conds = Dict()
    for θᵢ in keys(θtrue) # setup for ind params
        var = string(θᵢ)
        θcond = NamedTupleTools.select(θtrue, (θᵢ,))
        inf_conds[var] = solve(lat_mod, θcond; save_idxs=2, saveat=1).u # simulate x∣θᵢ
    end
    for var in keys(get_var_transforms()) # setup for param transformations
        sir_cond = sir_from_samples(var, M)
        inf_conds[var] = solve(sir_cond; save_idxs=2, saveat=1).u # simulate x|u
    end

    (;θtrue, θprior, obs_mod, inf_pri, inf_conds, inf_true, y)
end

"""
Run experiment 1: practical identifiability over time

Returns a `Vector{NamedTuple}` containing results for each tspan
"""
function ident_over_tspan(config; N=100)      
    @unpack θtrue, θprior, obs_mod, inf_pri, inf_conds, inf_true, y = config

    ret = []
    for imax in 2:length(inf_true) # for each timespan (0, t)
        println(imax)
        ts_idx = 1:imax

        for var in keys(inf_conds)
            δ = marginal_divergence(y[ts_idx], inf_conds[var][ts_idx], inf_pri[ts_idx], obs_mod)
            push!(ret, (var=var, t=imax-1, md=δ))
        end
    end
    return ret
end

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
obs_mod = PoissonRate(1000)
exper_config = exper1_setup(θtrue, θprior, obs_mod; M=60_000)

res = ident_over_tspan(exper_config; N=3000)
# CSV.write("data/res-pi-over-time.csv", DataFrame(res))