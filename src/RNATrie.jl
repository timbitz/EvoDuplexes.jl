#__precompile__()

module RNATrie

using Bio.Seq
using Bio.Intervals

include("pairs.jl")
include("energy.jl")
include("duplex.jl")
include("regions.jl")
include("trie.jl")

end
