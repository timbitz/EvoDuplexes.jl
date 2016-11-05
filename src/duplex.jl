
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
      cur_energy += stack_energy( pair, duplex.path[end] ) # TODO CHECK!
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

# Figure out what kind of bulge/internal loop we have prior to
# the curpair closing RNAPair then score it appropriately
function motif_energy( duplex::RNADuplex, curpair::Int )

   bulge    = false
   mismatch = false
   range    = 1:0

   energy   = 0.0 # return value
   # trace back size of internal loop/bulge motif
   for i in (curpair-1):-1:1
      if     isa( duplex.path[i], RNAPair )
         range = i:curpair
         break
      elseif isa( duplex.path[i], RNAMismatch )
         mismatch = true
      elseif isa( duplex.path[i], RNABulge )
         bulge    = true
      else
         error("$( duplex.path[i] ) is not a valid NucleotidePair!")
      end
   end

   # assign proper energy function to motif
   if     bulge && !mismatch
      # bulge_energy()
   elseif mismatch && !bulge
      # mismatch_energy() 1x1 or 2x2 or ...
   else # mismatch && bulge
      if length(range) == 4 # 1x2 Internal Loop
         # mismatch_1x2
      else
         # complex_energy()
      end
   end

   energy
end

function bulge_energy( duplex::RNADuplex, range::UnitRange )
   energy = 0.0
   energy += bulge_exception( duplex, range ) # bulge C exception
   if length(range) > 3 # multiple nt bulge
      
   else # single nt bulge
      energy += stack_energy( duplex.path[range.start], duplex.path[range.stop] )
      energy += # states?
   end

   energy
end

function internal_1x1_energy( duplex::RNADuplex, range::UnitRange )
   (x,y) = split(duplex.path[range.start+1])
   TURNER_2004_INTERNAL_TWO[ index(duplex.path[range.start]),
                             index(duplex.path[range.stop]) ][x,y]
end

function internal_2x2_energy( duplex::RNADuplex, range::UnitRange )
   TURNER_2004_INTERNAL_FOUR[ index(duplex.path[range.start]),
                              index(duplex.path[range.stop]) ][index(range.start+1),index(range.stop-1)]
end

function internal_1x2_energy( duplex::RNADuplex, range::UnitRange )
   # flip?
end

function internal_complex_energy( duplex::RNADuplex, range::UnitRange )
   
end

function bulge_exception( duplex::RNADuplex, range::UnitRange )
   if (duplex.path[range.start] == CG_PAIR && duplex.path[range.start+1] == CB_BULGE) ||
      (duplex.path[range.start] == GC_PAIR && duplex.path[range.start+1] == BC_BULGE) ||
      (duplex.path[range.stop]  == CG_PAIR && duplex.path[range.stop-1]  == CB_BULGE) ||
      (duplex.path[range.stop]  == GC_PAIR && duplex.path[range.stop-1]  == BC_BULGE)
      return TURNER_2004_BULGE_EXCEPTION
   else
      return 0.0
   end
end
