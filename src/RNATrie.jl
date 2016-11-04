#__precompile__()

module RNATrie

using Bio.Seq

include("pairs.jl")
include("energy.jl")
include("trie.jl")

end
