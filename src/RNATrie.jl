#__precompile__()

module RNATrie

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

end
