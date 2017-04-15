#__precompile__()

module RNATries

using Bio
using Bio.Seq
using Bio.Intervals

importall Bio.Seq
importall Bio.Intervals

include("pairs.jl")
include("energy.jl")
include("duplex.jl")
include("intervals.jl")
include("trie.jl")

export RNATrie,
       DuplexTrie,
       DuplexCollection,
       traverse,
       stitch

end
