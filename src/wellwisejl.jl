module wellwisejl

export sendto, addmoreprocs, ccdat, rbern, expit, summarygibbs, runmod

using Distributed
using DataFrames
import DataFrames: completecases
using Statistics
using StatsBase
using MCMCDiagnostics




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


function addmoreprocs(chains::Int=2)
  if(nprocs()<chains) # need to run this before anything!
    addprocs(chains-nprocs())
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
   res = mapreduce(x -> isa(x, Array) ? flat(x) : x, append!, arr,init=[])
   convert(Array{Float64}, res)
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
 sets, means, medians, pl, pu, stds, ac1, ac5, ess, lens = Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[]
 nm = names(results)
 for i in 1:size(results, 2)
   col = results[:,i]
   means = vcat(means, mean(col))
   medians = vcat(medians, median(col))
   pl = vcat(pl, quantile(col, 0.025)[1])
   pu = vcat(pu, quantile(col,  0.975)[1])
   stds = vcat(stds, std(col))
   ac = autocor(col, [1,5])
   ac1 = vcat(ac1, ac[1])
   ac5 = vcat(ac5, ac[2])
   ess = vcat(ess, effective_sample_size(col))
   lens = vcat(lens, length(col))
 end
 res = convert(DataFrame, hcat(nm, means, stds, medians, pl, pu, ess, ac1, ac5, lens))
 rename!(res, [:nm, :mean, :std, :median, :lower2_5, :upper97_5, :ess, :autocor_1, :autocor_5, :length])
 return res
end

function summarygibbs(results::Dict{Any,Any})
 sets, means, medians, pl, pu, stds, ac1, ac5, ess, rhat, lens = Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[], Array[]
 nm = names(results[1])
 for i in 1:size(results[1], 2)
   col = flat([vcat(r[2][:,i]) for r in results])
   means = vcat(means, mean(col))
   medians = vcat(medians, median(col))
   pl = vcat(pl, quantile(col, 0.025)[1])
   pu = vcat(pu, quantile(col,  0.975)[1])
   stds = vcat(stds, std(col))
   ac = autocor(col, [1,5])
   ac1 = vcat(ac1, ac[1])
   ac5 = vcat(ac5, ac[2])
   ess = vcat(ess, effective_sample_size(col))
   rhat = vcat(rhat,  potential_scale_reduction([vcat(r[2][:,i]) for r in results]...))
   lens = vcat(lens, length(col))
 end
 res = convert(DataFrame, hcat(nm, means, stds, medians, pl, pu, ess, rhat, ac1, ac5, lens))
 rename!(res, [:nm, :mean, :std, :median, :lower2_5, :upper97_5, :ess, :rhat, :autocor_1, :autocor_5, :length])
 return res
end


##########
# dpp related functions
##########



end # module
