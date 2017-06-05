
# Each base as 'one hot' nucleotide encoding
# 
# A   = 0b0001
# C   = 0b0010
# G   = 0b0100
# U   = 0b1000
# Gap = 0b0000

# Pairs encoded as adjacent sets of 'one hot' values
# so AA would be 0b00010001, and AU  0b00011000 and so forth.

abstract NucleotidePair

bitstype 8 RNAPair <: NucleotidePair

const NA_PAIR = reinterpret(RNAPair, 0b00000000)

const AU_PAIR = reinterpret(RNAPair, 0b00011000)
const UA_PAIR = reinterpret(RNAPair, 0b10000001)
const CG_PAIR = reinterpret(RNAPair, 0b00100100)
const GC_PAIR = reinterpret(RNAPair, 0b01000010)
const GU_PAIR = reinterpret(RNAPair, 0b01001000)
const UG_PAIR = reinterpret(RNAPair, 0b10000100)

bitstype 8 RNAMismatch <: NucleotidePair

const AA_MISMATCH = reinterpret(RNAMismatch, 0b00010001)
const AC_MISMATCH = reinterpret(RNAMismatch, 0b00010010)
const AG_MISMATCH = reinterpret(RNAMismatch, 0b00010100)
const CA_MISMATCH = reinterpret(RNAMismatch, 0b00100001)
const CC_MISMATCH = reinterpret(RNAMismatch, 0b00100010)
const CU_MISMATCH = reinterpret(RNAMismatch, 0b00101000)
const GA_MISMATCH = reinterpret(RNAMismatch, 0b01000001)
const GG_MISMATCH = reinterpret(RNAMismatch, 0b01000100)
const UC_MISMATCH = reinterpret(RNAMismatch, 0b10000010)
const UU_MISMATCH = reinterpret(RNAMismatch, 0b10001000)

bitstype 8 RNABulge <: NucleotidePair

const AB_BULGE = reinterpret(RNABulge, 0b00010000)
const CB_BULGE = reinterpret(RNABulge, 0b00100000)
const GB_BULGE = reinterpret(RNABulge, 0b01000000)
const UB_BULGE = reinterpret(RNABulge, 0b10000000)

const BA_BULGE = reinterpret(RNABulge, 0b00000001)
const BC_BULGE = reinterpret(RNABulge, 0b00000010)
const BG_BULGE = reinterpret(RNABulge, 0b00000100)
const BU_BULGE = reinterpret(RNABulge, 0b00001000)

Base.convert{NP <: NucleotidePair}(::Type{NP}, x::UInt8) = reinterpret(NP, x)
Base.convert(::Type{UInt8}, x::NucleotidePair)           = reinterpret(UInt8, x)

function Base.convert{NP <: NucleotidePair}(::Type{NP}, first::Bio.Seq.Nucleotide,
                                                        last::Bio.Seq.Nucleotide)
   enc  = reinterpret(UInt8, first) << 4
   enc |= reinterpret(UInt8, last)
   reinterpret(NP, enc)
end

function Base.convert{NP <: NucleotidePair}(::Type{NP}, first::UInt8, last::UInt8 )
   enc  = first << 4
   enc |= last
   reinterpret(NP, enc)
end

# convert two adjacent 'one hot' nucleotide 4-bit codes 
# into two adjacent 2-bit encoding 0b00 = A, 0b01 = C, 0b10 = G, 0b11 = U
# for both left and right nucleotides, so 0b00010001 becomes 0b0000
# and 0b10001000 becomes 0b1111, then add 1 to bring to 1:16 range
function index(x::NucleotidePair)
   idx   = trailing_zeros(convert(UInt8, x) >> 4) << 2
   right = trailing_zeros(convert(UInt8, x) & 0x0F)
   (idx | right) + 1
end

function index(x::RNABulge)
   trailing_zeros(reinterpret(UInt8, x)) + 1
end

isbulge(x::RNAPair)     = false
isbulge(x::RNAMismatch) = false
isbulge(x::RNABulge)    = true

ispair(x::RNAPair)     = true
ispair(x::RNAMismatch) = false
ispair(x::RNABulge)    = false

ismismatch(x::RNAPair)     = false
ismismatch(x::RNAMismatch) = true
ismismatch(x::RNABulge)    = false

function type_cnt{NP <: NucleotidePair}( v::Vector{NP}, func::Function )
   cnt = 0
   @inbounds for i in 1:length(v)
      func(v[i]) && (cnt += 1)
   end
   cnt
end

npairs{NP <: NucleotidePair}( v::Vector{NP} )      = type_cnt( v, ispair )
nmismatches{NP <: NucleotidePair}( v::Vector{NP} ) = type_cnt( v, ismismatch )
nbulges{NP <: NucleotidePair}( v::Vector{NP} )     = type_cnt( v, isbulge )

isfiveprime(x::RNABulge) = reinterpret(UInt8, x) & 0x0F == 0 ? true : false

isbulgefive(x::RNABulge)    = isfiveprime(x)
isbulgefive(x::RNAPair)     = false
isbulgefive(x::RNAMismatch) = false

isbulgethree(x::RNABulge)    = !isfiveprime(x)
isbulgethree(x::RNAPair)     = false
isbulgethree(x::RNAMismatch) = false

function flip{NP <: NucleotidePair}(x::NP)
   left   = convert(UInt8, x) >> 4
   retval = (convert(UInt8, x) & 0x0F) << 4
   reinterpret( NP, retval | left )
end

function Base.split(x::NucleotidePair)
   left  = trailing_zeros(convert(UInt8, x) >> 4) + 1
   right = trailing_zeros(convert(UInt8, x) & 0x0F) + 1
   (left, right)
end

function Base.split(x::RNABulge)
   x_int = reinterpret(UInt8, x)
   isfiveprime(x) ? trailing_zeros(x_int >> 4) + 1 : trailing_zeros(x_int) + 1
end

# this function returns the first, or last if flag is true
# five prime or three prime indexes (1:4) in a path of NucleotidePairs
function five_three{NP <: NucleotidePair}( path::Vector{NP}, range::UnitRange; 
                                           last::Bool=false )
   five  = 0
   three = 0
   for i in (last ? reverse(range) : range)
      if isa( path[i], RNABulge )
         if isfiveprime( path[i] )
            five = split( path[i] )
         else
            three = split( path[i] )
         end
      else
         x,y   = split( path[i] )
         five  = five > 0 ? five : x
         three = three > 0 ? three : y
      end
      if five > 0 && three > 0
         break
      end
   end
   five, three
end

function is_purine_pyrimidine(x::RNAPair)
   if x == AU_PAIR || x == GC_PAIR || x == GU_PAIR
      return true
   else
      return false
   end
end

