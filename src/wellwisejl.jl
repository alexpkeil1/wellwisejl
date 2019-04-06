module wellwisejl

export sendto, addmoreprocs, ccdat, rbern, expit, summarygibbs, runmod, traceplot, densplot

using Distributed
using DataFrames
import DataFrames: completecases
using Statistics
using Makie
using KernelDensity, Interpolations




#struct GibbsSampler{I <: Int64, B <: Int64, D <: DataFrame} <: DataFrame
#    iter::I
#    burnin::B
#    dat::D
#    @doc """
#    $(SIGNATURES)
#
#    """
#    function GibbsSampler(iter::I, burnin::B, dat::D) where {I <: Int64, B <: Int64, D <: DataFrame}
#    end
#end


function addmoreprocs()
  if(nprocs()<2) # need to run this before anything!
    addprocs(2-nprocs())
  end
end 

function sendto(p::Int; args...)
  for (nm, val) in args
    @spawnat(p, eval(Main, Expr(:(=), nm, val)))
  end
end 

function sendto(ps::Vector{Int}; args...)
  for p in ps
    sendto(p; args...)
  end
end 

function ccdat(data, i::Int=1)
  cols = [:Arsenic, :Cadmium, :Lead, :Manganese, :Copper, :y];
  cc = completecases(data[:,cols]);
  dat = data[(cc .& (data[:iter] .== i)) ,cols];
  return dat
end 

function rbern(mu::Float64)
  x::Int8  = (rand()<mu)*1
end 

function rbern(mu::Array{Float64,1})
  N::Int64 = length(mu)
  x::Array{Int8,1}  = (rand(N) .< mu) .* 1
end 

function expit(mu::Float64)
  x::Float64  = 1/(1+exp(-mu))
end 

function expit(mu::Array{Float64,1})
  x::Array{Float64,1}  = 1. ./ (1. .+ exp.( .- mu))
end
      
function expit(mu::Array{Union{Missing, Float64},1})
  x::Array{Float64,1}  = 1. ./ (1. .+ exp.( .- mu))
end 

function flat(arr::Array)
   mapreduce(x -> isa(x, Array) ? flat(x) : x, append!, arr,init=[])
end

function runmod(sampler::Function, rdat::DataFrame, niter::NI, burnin=0, chains=4) where {NI<:Integer}
  futureres, res = Dict(), Dict()
  sendto([i for i in procs()[1:chains]], dat=rdat)
  for i in procs()[1:chains]
    @spawnat i global dat = rdat
    push!(futureres,  i => @spawnat i sampler(niter, burnin, dat));
  end
  for i in procs()[1:chains]
    push!(res, i => fetch(futureres[i]))
  end
  return res
end

function summarygibbs(results::DataFrame)
 sets, means, medians, pl, pu, stds, lens = Array[], Array[], Array[], Array[], Array[], Array[], Array[]
 nm = names(results)
 for i in 1:size(results, 2)
   col = results[:,i]
   means = vcat(means, mean(col))
   medians = vcat(medians, median(col))
   pl = vcat(pl, quantile(col, 0.025)[1])
   pu = vcat(pu, quantile(col,  0.975)[1])
   stds = vcat(stds, std(col))
   lens = vcat(lens, length(col))
 end
 res = convert(DataFrame, hcat(nm, means, stds, medians, pl, pu, lens))
 names!(res, [:nm, :mean, :std, :median, :lower2_5, :upper97_5, :length])
 return res
end

function summarygibbs(results::Dict{Any,Any})
 sets, means, medians, pl, pu, stds, lens = Array[], Array[], Array[], Array[], Array[], Array[], Array[]
 nm = names(results[1])
 for i in 1:size(results[1], 2)
   col = flat([vcat(r[2][:,i]) for r in results])
   means = vcat(means, mean(col))
   medians = vcat(medians, median(col))
   pl = vcat(pl, quantile(col, 0.025)[1])
   pu = vcat(pu, quantile(col,  0.975)[1])
   stds = vcat(stds, std(col))
   lens = vcat(lens, length(col))
 end
 res = convert(DataFrame, hcat(nm, means, stds, medians, pl, pu, lens))
 names!(res, [:nm, :mean, :std, :median, :lower2_5, :upper97_5, :length])
 return res
end

function traceplot(res::Dict{Any,Any}, pos::Integer, burnin::Integer)
    Iteration = [i for i in 1:(size(res[1], 1) - burnin+1)]
    cols = [:black, :blue, :red, :green, :yellow]
    # make plot
    parname = string(names(res[1])[pos])
	scene = Scene()
    for chain in res
      lines!(scene, Iteration, chain[2][(burnin+1):end,pos], 
        color = cols[chain[1]],
        axis = (
           names = (axisnames = ("Iteration", parname),),
           grid = (linewidth = (0, 0),),
           )
           )
    end
    return scene
end


function traceplot(res::Dict{Any,Any}, colnm::Symbol, burnin::Integer)
    Iteration = [i for i in 1:(size(res[1], 1) - burnin+1)]
    cols = [:black, :blue, :red, :green, :yellow]
    # make plot
    parname = string(colnm)
	scene = Scene()
    for chain in res
      lines!(scene, Iteration, chain[2][colnm][(burnin+1):end], 
        color = cols[chain[1]],
        axis = (
           names = (axisnames = ("Iteration", parname),),
           grid = (linewidth = (0, 0),),
           )
           )
    end
    return scene
end
traceplot(res::Dict{Any,Any}, pos::Integer) = traceplot(res, pos, 0)
traceplot(res::Dict{Any,Any}, colnm::Symbol) = traceplot(res, colnm, 0)


function densplot(res::Dict{Any,Any}, pos::Integer, burnin::Integer)
    Iteration = [i for i in 1:(size(res[1], 1) - burnin+1)]
    cols = [:black, :blue, :red, :green, :yellow]
    # make plot
    parname = string(names(res[1])[pos])
	scene = Scene()
    for chain in res
      vals = chain[2][(burnin+1):end,pos]
      rng = (minimum(vals), maximum(vals))
      x = range(rng[1], stop=rng[2], length=101)
      kdens = kde(vals)
      dens = pdf(kdens, x)
      lines!(scene, x, dens, 
        color = cols[chain[1]],
        axis = (
           names = (axisnames = (parname, "Density"),),
           grid = (linewidth = (0, 0),),
           )
           )
    end
    return scene
end

function densplot(res::Dict{Any,Any}, colnm::Symbol, burnin::Integer)
    Iteration = [i for i in 1:size(res[1], 1)]
    cols = [:black, :blue, :red, :green, :yellow]
    # make plot
    parname = string(colnm)
	scene = Scene()
    for chain in res
      vals = chain[2][colnm][(burnin+1):end]
      rng = (minimum(vals), maximum(vals))
      x = range(rng[1], stop=rng[2], length=101)
      kdens = kde(vals)
      dens = pdf(kdens, x)
      lines!(scene, x, dens, 
        color = cols[chain[1]],
        axis = (
           names = (axisnames = (parname, "Density"),),
           grid = (linewidth = (0, 0),),
           )
           )
    end
    return scene
end
densplot(res::Dict{Any,Any}, pos::Integer) = densplot(res, pos, 0)
densplot(res::Dict{Any,Any}, colnm::Symbol) = densplot(res, colnm, 0)


end # module
