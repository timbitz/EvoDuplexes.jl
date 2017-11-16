#  Functions and Constants required for traversal

const DNAPAIRS = [(DNA_A, DNA_T),   (DNA_T, DNA_A),
                  (DNA_G, DNA_C),   (DNA_C, DNA_G),
                  (DNA_G, DNA_T),   (DNA_T, DNA_G)]

const DNAGAPS  = [(DNA_A, DNA_Gap), (DNA_Gap, DNA_A),
                  (DNA_C, DNA_Gap), (DNA_Gap, DNA_C),
                  (DNA_G, DNA_Gap), (DNA_Gap, DNA_G),
                  (DNA_T, DNA_Gap), (DNA_Gap, DNA_T)]

const DNAMISMATCH = [(DNA_A, DNA_A), (DNA_A, DNA_C),
                     (DNA_A, DNA_G), (DNA_C, DNA_A),
                     (DNA_C, DNA_C), (DNA_C, DNA_T),
                     (DNA_G, DNA_A), (DNA_G, DNA_G),
                     (DNA_T, DNA_C), (DNA_T, DNA_T)]

const RNAPAIRS = [(RNA_A, RNA_U), (RNA_U, RNA_A),
                  (RNA_G, RNA_C), (RNA_C, RNA_G),
                  (RNA_G, RNA_U), (RNA_U, RNA_G)]

const RNAGAPS  = [(RNA_A, RNA_Gap), (RNA_Gap, RNA_A),
                  (RNA_C, RNA_Gap), (RNA_Gap, RNA_C),
                  (RNA_G, RNA_Gap), (RNA_Gap, RNA_G),
                  (RNA_U, RNA_Gap), (RNA_Gap, RNA_U)]

const RNAMISMATCH = [(RNA_A, RNA_A), (RNA_A, RNA_C),
                     (RNA_A, RNA_G), (RNA_C, RNA_A),
                     (RNA_C, RNA_C), (RNA_C, RNA_U),
                     (RNA_G, RNA_A), (RNA_G, RNA_G),
                     (RNA_U, RNA_C), (RNA_U, RNA_U)]

const PairsType = Vector{Tuple{Bio.Seq.Nucleotide, Bio.Seq.Nucleotide}}

function twobitalpha{S <: Bio.Seq.Sequence}( seq::S )
   t = eltype(seq)
   if t == DNANucleotide
      return DNAAlphabet{2}
   elseif t == RNANucleotide
      return RNAAlphabet{2}
   else
      error("RNASuffixArrays are only built for DNA/RNA input sequences!")
   end
end

pairs{n}(::Type{Bio.Seq.DNAAlphabet{n}})      = DNAPAIRS
pairs{n}(::Type{Bio.Seq.RNAAlphabet{n}})      = RNAPAIRS
gaps{n}(::Type{Bio.Seq.DNAAlphabet{n}})       = DNAGAPS
gaps{n}(::Type{Bio.Seq.RNAAlphabet{n}})       = RNAGAPS
mismatches{n}(::Type{Bio.Seq.DNAAlphabet{n}}) = DNAMISMATCH
mismatches{n}(::Type{Bio.Seq.RNAAlphabet{n}}) = RNAMISMATCH

gap{n}(::Type{DNAAlphabet{n}})                = DNA_Gap
gap{n}(::Type{RNAAlphabet{n}})                = RNA_Gap

onehot{I <: Integer}(x::I) = 0x01 << (x-1)

encodeindex{A <: Alphabet}(::Type{A}, x::Bio.Seq.Nucleotide) = (x == DNA_Gap || x == RNA_Gap) ? 0 : Bio.Seq.encode(A, x)+1
index{A <: Alphabet}(::Type{A}, func=pairs) = map( x->map(y->encodeindex(A, y), x), func(A) )

revoffset( x::Int, seq::BioSequence ) = revoffset( x, length(seq) )
revoffset( x, len ) = len - x + 1

Base.reverse( seq::Bio.Seq.ReferenceSequence ) = Bio.Seq.ReferenceSequence( reverse( String( seq ) ) )

