using DrWatson
@quickactivate "optimal-test-design"
using DiffEqInformationTheory
using ConditionalTransform
using Distributions, MonteCarloMeasurements
using NamedTupleTools

include(srcdir("watson-tools.jl"))
include(srcdir("transform-funcs.jl"))

function sir_from_samples(varname)
    df = CSV.read(datadir("sims", "cond-samples", "samp-$(varname).csv"), DataFrame)

    SIRModel{Float32}(
        S₀=Particles(Float32.(collect(df.S0))), 
        β=Particles(Float32.(collect(df.beta))), 
        α=Particles(Float32.(collect(df.alpha)))
    )
end

function inc_tspan_exper!(config; N=100, M=60_000)
    @unpack θtrue, θprior, obs_mod, var_transforms = config
    
    # particle vectors that are constant wrt ϕ:
    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_pri = solve(lat_mod; save_idxs=2, saveat=1).u 
    inf_true = solve(lat_mod, θtrue; save_idxs=2, saveat=1).u
    y = Particles(40_000, observe_dist(obs_mod; observe_params(obs_mod, inf_true)...))

    inf_conds = Dict()
    for θᵢ in keys(θtrue) # setup for ind params
        lab = "md-"*string(θᵢ)
        config[lab] = []
        θcond = NamedTupleTools.select(θtrue, (θᵢ,))
        inf_conds[lab] = solve(lat_mod, θcond; save_idxs=2, saveat=1).u # simulate x∣θᵢ
    end
    for vt in keys(var_transforms) # setup for param transformations
        lab = "md-"*vt
        config[lab] = []
        sir_cond = sir_from_samples(vt)
        inf_conds[lab] = solve(sir_cond; save_idxs=2, saveat=1).u
    end
        
    for imax in 2:length(inf_true) # for each timespan (0, t)
        println(imax)
        ts_idx = 1:imax
        ynew = bootstrap(y[ts_idx], N)

        for lab in keys(inf_conds)
            push!(config[lab], marginal_divergence(ynew, inf_conds[lab][ts_idx], inf_pri[ts_idx], obs_mod))
        end
    end
end

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
obs_mod = PoissonTests(1000)
var_transforms = get_var_transforms()

fdir = datadir("sims", "increasing-tspan-marg")

exper = @strdict θtrue θprior obs_mod var_transforms
fname = mysavename(exper; ignores=[:var_transforms])
inc_tspan_exper!(exper; N=3000)
tagsave("$fdir/$fname", exper; safe=true)


