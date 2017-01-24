

using DataFrames
using Gadfly
using Mamba
using Distributions
using StatsBase

dat = readtable("test.txt.gz", separator = '\t', header=false, names=[:Length, :DeltaG])

hex=plot(dat, x=:Length, y=:DeltaG, Geom.hexbin(xbincount=35), 
         Coord.cartesian(xmin=0, ymax=15),
         Guide.xlabel("Duplex Length"), 
         Guide.ylabel("Gibbs Free Energy Change (ΔG)"))

dens=plot(dat[dat[:Length] .>= 5,:], color=:Length, x=:DeltaG, Geom.density, Guide.colorkey("Duplex Length"), 
          Guide.xlabel("Gibbs Free Energy Change (ΔG)"))

dd = DuplexDist( dat )
model=plotcdf( dd, dat[:Length], dat[:DeltaG] )

draw(PDF("deltaG.pdf", 8inch, 8inch), vstack(hstack(hex,dens), model))

type DuplexDist
   length::Distributions.Poisson
   deltag::Vector{Distributions.Normal}
end

# G[L] ~ Normal(mu[L], s[L])
# L ~ Poisson(lambda)
function DuplexDist( df::DataFrame )
   # Fit Model
   dgdat = Vector{DataArray}()
   for i in 1:maximum(df[:Length])
      push!( dgdat, dat[dat[:,1] .== i, 2] )
   end
   lendat = df[:Length]

   len = fit(Poisson, lendat)
   dg  = Vector{Distributions.Normal}(length(dgdat))
   for i in 1:N
      dg[i] = fit(Normal, dgdat[i])
   end

   DuplexDist( len, dg )
end

function cdf( dd::DuplexDist, len::Int, dg::Float64 )
   prob = 0.0
   for i in len:length(dd.deltag)
      prob += pdf( dd.length, i ) * cdf( dd.deltag[i], dg )     
   end
   prob += ccdf( dd.length, length(dd.deltag)+1 )
   prob
end

# P(G,L) = P(G|L)P(L)
pdf( dd::DuplexDist, len::Int, dg::Float64 ) = pdf( dd.length, len ) * pdf( dd.deltag[len], dg )

function plotcdf( dd::DuplexDist, lenvec, dgvec )
   l = Vector{Int}()
   d = Vector{Float64}()
   c = Vector{Float64}()
   p = Vector{Float64}()
   irange = 1:maximum(lenvec)
   jrange = minimum(dgvec):1.0:maximum(dgvec)
   for i in irange
      for j in jrange
         push!( l, i )
         push!( d, j )
         push!( c, cdf( dd, i, j ) )
         push!( p, pdf( dd, i, j ) )
      end
   end
   z = reshape( p, (length(jrange), length(irange)))
 
   pval = plot(layer(x=l, y=d, color=c, Geom.rectbin),
         Coord.cartesian(xmin=0, ymax=15),
         Guide.xlabel("Duplex Length"),
         Guide.ylabel("Gibbs Free Energy Change (ΔG)"),
         Guide.colorkey("p-value"), Scale.color_log10)

   prob =plot(layer(x=l, y=d, color=p, Geom.rectbin),
         Coord.cartesian(xmin=0, ymax=15),
         Guide.xlabel("Duplex Length"),
         Guide.ylabel("Gibbs Free Energy Change (ΔG)"),
         Guide.colorkey("Pdf(x)")) 

   hstack(prob, pval)
end

#= MAMBA

fitdat = Dict{Symbol, Any}(
  :len => len,
  :dg => dg
)
N = length(dg)

model = Model(
   dg = Stochastic(1, (len, mu, s2) -> UnivariateDistribution[ Normal(mu[i], sqrt(s2[i])) for i in 1:N ], false),

   mu = Stochastic(1, () -> Normal(0, 1)),
   s2 = Stochastic(1, () -> InverseGamma(0.01, 0.01)),

   len = Stochastic(alpha -> begin
                               lambda = max(0, 1/alpha)
                               Poisson(lambda)
                             end, false),
   alpha = Stochastic(() -> Exponential(1/10)),
)
   
inits = [
  Dict{Symbol, Any}(
     :alpha => rand(Exponential(1/10)),
     :len => mean(fitdat[:len]),
     :s2 => [rand(InverseGamma(0.01, 0.01)) for i in 1:N],
     :mu => [mean(dg[i]) for i in 1:N],
     :dg => [mean(dg[i]) for i in 1:N]
  )
for j in 1:5
]

scheme = [Slice(:alpha, 1.0, Univariate),
          Slice(:len, 1.0, Univariate),
          NUTS([:mu, :s2])]

setinputs!(model, fitdat)
setinits!(model, inits[1])
setsamplers!(model, scheme)

sim = mcmc(model, fitdat, inits, 10000, burnin=250, thin=2, chains=1)

=#
