
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

isfiveprime(x::RNABulge) = reinterpret(UInt8, x) & 0x0F == 0 ? true : false

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
