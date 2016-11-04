

# RNADuplex objects represent the path of a duplex as a
# Vector of NucleotidePair's.  Type of the NucleotidePair
# thus defines the operations: Pair, Mismatch, Bulge
# and interior loops are combinations of Mismatches and Bulges
type RNADuplex
   path::Vector{NucleotidePair} # sequence and structure path
   energy::Vector{Float64} # gibbs free energy from nearest-neighbor
   length::Int # shortest length

   RNADuplex() = new(NucleotidePair[], Float64[4.09], 0)
end

function Base.push!( duplex::RNADuplex, pair::RNAPair )
   cur_energy = duplex.energy[end]
   if length(duplex.path) >= 1 && isa( duplex.path[end], RNAPair )
      # add stack energy
      cur_energy += stack_energy( pair, duplex.path[end] )
   else
      # add bulge/loop motif from lookup table   
   end
   push!( duplex.energy, cur_energy )
   push!( duplex.path, pair )
end

function Base.push!{NP <: Union{RNAMismatch,RNABulge}}( duplex::RNADuplex, pair::NP )
   push!( duplex.energy, duplex.energy[end] ) # don't calculate energy yet
   push!( duplex.path, pair )
end

function Base.pop!{NP <: NucleotidePair}( duplex::RNADuplex )
   pop!( duplex.path )
   pop!( duplex.energy )
end


