    # wellwisejl
    using Pkg
    Pkg.add(PackageSpec(url="https://github.com/alexpkeil1/wellwisejl.git", rev="master"))

    # parallel
    using Distributed # implicitly loads when calling julia with: julia -p <ncores>
    # custom
    
    addprocs(3)
    
    
    @everywhere using PolyaGammaDistribution
    # math
    @everywhere using LinearAlgebra, Statistics, Distributions, AltDistributions
    # utility
    @everywhere using RData, DataFrames, CSV, Distributed
    
    @everywhere using wellwisejl
    
    
    ###### end header ######
    
    
    @everywhere function gibbs(iter::Int64=10, burnin::Int64=0, data=data)
      y = data[:y]
      kappa = y .- 0.5
      N = length(y)
      intercept = ones(N)
      Arsenic = data[:Arsenic]
      Manganese = data[:Manganese]
      Lead = data[:Lead]
      Cadmium = data[:Cadmium]
      Copper = data[:Copper]
      intprop = vcat(
                   [.98 .98 .98 .98 0],
                   [0 .985 0 .965 0],
                   [0 .918 0 0 0],
                   [.994 .98 .98 .987 0],
                   [0 0 0 .99 0]
                 )
       nints = size(intprop)[1]
     
      X = hcat(intercept, Arsenic, Manganese, Lead, Cadmium, Copper, Arsenic .* Arsenic,
         Arsenic .* Manganese, Arsenic .* Lead, Arsenic .* Cadmium, Arsenic .* Copper, 
         Manganese .* Manganese, Lead .* Manganese, Cadmium .* Manganese, Copper .* Manganese,
         Lead .* Lead, Cadmium .* Lead, Copper .* Lead, Cadmium .* Cadmium, Cadmium .* Copper,
         Copper .* Copper)
       Xil = Dict()
       for j in 1:nints
         push!(Xil,  j => hcat(intercept, (1-intprop[j,1]) .* Arsenic, (1-intprop[j,2]) .* Manganese, (1-intprop[j,
             3]) .* Lead, (1-intprop[j,4]) .* Cadmium, (1-intprop[j,5]) .* Copper, (1-intprop[j,
             1]) .* Arsenic .* (1-intprop[j,1]) .* Arsenic, (1-intprop[j,1]) .* Arsenic .* (1-intprop[j,
             2]) .* Manganese, (1-intprop[j,1]) .* Arsenic .* (1-intprop[j,3]) .* Lead, (1-intprop[j,
             1]) .* Arsenic .* (1-intprop[j,4]) .* Cadmium, (1-intprop[j,1]) .* Arsenic .* (1-intprop[j,
             5]) .* Copper, (1-intprop[j,2]) .* Manganese .* (1-intprop[j,2]) .* Manganese, (1-intprop[j,
             3]) .* Lead .* (1-intprop[j,2]) .* Manganese, (1-intprop[j,4]) .* Cadmium .* (1-intprop[j,
             2]) .* Manganese, (1-intprop[j,5]) .* Copper .* (1-intprop[j,2]) .* Manganese, (1-intprop[j,
             3]) .* Lead .* (1-intprop[j,3]) .* Lead, (1-intprop[j,4]) .* Cadmium .* (1-intprop[j,3]) .* Lead,
             (1-intprop[j,5]) .* Copper .* (1-intprop[j,3]) .* Lead, (1-intprop[j,4]) .* Cadmium .* (1-
             intprop[j,4]) .* Cadmium, (1-intprop[j,4]) .* Cadmium .* (1-intprop[j,5]) .* Copper, (1-
             intprop[j,5]) .* Copper .* (1-intprop[j,5]) .* Copper))
       end
      # initialize storage vectors
      p = size(X)[2]
      beta = Array{Float64}(undef, iter, p)
      rd = Array{Float64}(undef, iter, nints)
      b = rand(p)
      # priors
      Bmu = zeros(p)
      prior_sig = 0.8 .* ones(p-1)
      Bsig = zeros(p,p)
      Bsig[diagind(Bsig)] = [100; prior_sig]
      invBsig = pinv(Bsig)
      tX = transpose(X)
      beta[1,:] = b
      for i in 2:iter
        #print('.')
        mu = X * b 
        for j in 1:nints 
          mui = Xil[j] * b
          rd[i-1,j] = mean(expit(mui)) - mean(expit(mu))
        end
        while true
          global om = [rand(PolyaGamma(1, mu[m])) for m in 1:N]
          txdox = (tX * Diagonal(om) * X) + invBsig
          global V = pinv( txdox  )
          if isposdef(Symmetric(V))
          	break
          	print(".")
          end
        end
        m = V * (tX * kappa + invBsig * Bmu)
        b = rand(MvNormal(m,Symmetric(V))) # needed Symmetric hack to hide float errors?
        beta[i,:] = b
      end
      mu = X * b 
      for j in 1:nints 
        mui = Xil[j] * b
        rd[iter,j] = mean(expit(mui)) - mean(expit(mu))
      end
      allres = hcat(beta, rd)
      res = convert(DataFrame, allres)
      names!(res, vcat([Symbol("b" * "$i") for i in 0:(p-1)], [Symbol("rd" * "$i") for i in 1:nints]))
      return res[(burnin+1):iter,:]
    end
    
    
    # read in data
    dt = CSV.read("testdata20190404.csv", missingstring="NA");
    data = ccdat(dt, 1)
    #single run of four chains, 10,000 iterations each with zero burnin
    res = runmod(gibbs, data, 10000, 0, 4)
    
    summarygibbs(res)

    


