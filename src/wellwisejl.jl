module wellwisejl

export sendto, addmoreprocs, ccdat, rbern, expit, summary, runmod

using Pkg, Distributed, LinearAlgebra, Statistics, Distributions,
       DelimitedFiles, RData, DataFrames, CSV, PolyaGammaDistribution, AltDistributions 


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

function addmoreprocs()
       if(nprocs()<2) # need to run this before anything!
        addprocs(2-nprocs())
        end
      end 


function ccdat(data)
        cols = [:Arsenic, :Cadmium, :Lead, :Manganese, :Copper, :y];
        cc = completecases(data[:,cols]);
        dat = data[(cc .& (data[:iter] .== 1)) ,cols];
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
        x::Float64  = 1/(1-exp(mu))
     end 

function expit(mu::Array{Float64,1})
        x::Array{Float64,1}  = 1. ./ (1. .+ exp.( .- mu))
     end
      

function expit(mu::Array{Union{Missing, Float64},1})
        x::Array{Float64,1}  = 1. ./ (1. .+ exp.( .- mu))
     end 

function summary(results, burn=0)
 sets, means, stds, lens = Array[], Array[], Array[], Array[]
 s = 1
 for set in results
   cset = [c for c in eachcol(set)]
   sets = vcat(sets, map(mean, [s for c in eachcol(set)]))
   means = vcat(means, map(mean, cset))
   stds = vcat(stds, map(std, cset))
   lens = vcat(lens, map(x -> size(x)[1], cset))
   s = s+1
 end
 return hcat(sets, means, stds, lens)
end

function runmod(niter::Int64, burnin::Int64=0, chains::Int64=4)
      futureres, res = Dict(), Dict()
      sendto([i for i in procs()[1:",chains,"]], dat=rdat)
      for i in procs()
       #@spawnat i global dat = rdat
       push!(futureres,  i => @spawnat i gibbs(niter, burnin, dat));
      end
      for i in procs()
        push!(res, i => fetch(futureres[i]))
      end
      return res
    end

end # module
