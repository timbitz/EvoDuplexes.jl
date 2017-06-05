#__precompile__()

module RNATries

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
       RNASuffix,
       DuplexTrie,
       DuplexSuffix,
       DuplexCollection,
       traverse,
       stitch

end
