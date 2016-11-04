
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

const TURNER_2004_STACK = fill( 0.0, (16,16) )

      TURNER_2004_STACK[index(AU_PAIR),index(AU_PAIR)] = -0.9
      TURNER_2004_STACK[index(AU_PAIR),index(UA_PAIR)] = -1.1
      TURNER_2004_STACK[index(AU_PAIR),index(CG_PAIR)] = -2.2
      TURNER_2004_STACK[index(AU_PAIR),index(GC_PAIR)] = -2.1
      TURNER_2004_STACK[index(AU_PAIR),index(GU_PAIR)] = -0.6
      TURNER_2004_STACK[index(AU_PAIR),index(UG_PAIR)] = -1.4

      TURNER_2004_STACK[index(CG_PAIR),index(AU_PAIR)] = -2.1
      TURNER_2004_STACK[index(CG_PAIR),index(UA_PAIR)] = -2.1
      TURNER_2004_STACK[index(CG_PAIR),index(CG_PAIR)] = -3.3
      TURNER_2004_STACK[index(CG_PAIR),index(GC_PAIR)] = -2.4
      TURNER_2004_STACK[index(CG_PAIR),index(GU_PAIR)] = -1.4
      TURNER_2004_STACK[index(CG_PAIR),index(UG_PAIR)] = -2.1

      TURNER_2004_STACK[index(GC_PAIR),index(AU_PAIR)] = -2.4
      TURNER_2004_STACK[index(GC_PAIR),index(UA_PAIR)] = -2.2
      TURNER_2004_STACK[index(GC_PAIR),index(CG_PAIR)] = -3.4
      TURNER_2004_STACK[index(GC_PAIR),index(GC_PAIR)] = -3.3
      TURNER_2004_STACK[index(GC_PAIR),index(GU_PAIR)] = -1.5
      TURNER_2004_STACK[index(GC_PAIR),index(UG_PAIR)] = -2.5

      TURNER_2004_STACK[index(GU_PAIR),index(AU_PAIR)] = -1.3
      TURNER_2004_STACK[index(GU_PAIR),index(UA_PAIR)] = -1.4
      TURNER_2004_STACK[index(GU_PAIR),index(CG_PAIR)] = -2.5
      TURNER_2004_STACK[index(GU_PAIR),index(GC_PAIR)] = -2.1
      TURNER_2004_STACK[index(GU_PAIR),index(GU_PAIR)] = -0.5 # +0.5
      TURNER_2004_STACK[index(GU_PAIR),index(UG_PAIR)] = +1.3

      TURNER_2004_STACK[index(UA_PAIR),index(AU_PAIR)] = -1.3
      TURNER_2004_STACK[index(UA_PAIR),index(UA_PAIR)] = -0.9
      TURNER_2004_STACK[index(UA_PAIR),index(CG_PAIR)] = -2.4
      TURNER_2004_STACK[index(UA_PAIR),index(GC_PAIR)] = -2.1
      TURNER_2004_STACK[index(UA_PAIR),index(GU_PAIR)] = -1.0
      TURNER_2004_STACK[index(UA_PAIR),index(UG_PAIR)] = -1.3

      TURNER_2004_STACK[index(UG_PAIR),index(AU_PAIR)] = -1.0
      TURNER_2004_STACK[index(UG_PAIR),index(UA_PAIR)] = -0.6
      TURNER_2004_STACK[index(UG_PAIR),index(CG_PAIR)] = -1.5
      TURNER_2004_STACK[index(UG_PAIR),index(GC_PAIR)] = -1.4
      TURNER_2004_STACK[index(UG_PAIR),index(GU_PAIR)] = +0.3
      TURNER_2004_STACK[index(UG_PAIR),index(UG_PAIR)] = -0.5


# duplex object
type RNADuplex
   
end
