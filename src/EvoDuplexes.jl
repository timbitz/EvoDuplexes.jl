#__precompile__()

module EvoDuplexes.jl

using Bio
using Bio.Seq
using Bio.Intervals
using SuffixArrays

importall Bio.Seq
importall Bio.Intervals

include("pairs.jl")
include("energy.jl")
include("duplex.jl")
include("intervals.jl")
include("traverse.jl")
include("trie.jl")
include("suffix.jl")

export RNATrie,
       DuplexTrie,
       RNASuffixArray,
       RNADuplexArray,
       DuplexCollection,
       traverse,
       stitch

end
