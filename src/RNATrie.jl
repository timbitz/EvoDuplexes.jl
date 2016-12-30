#__precompile__()

module RNATrie

using Bio.Seq
using Bio.Intervals

import Bio.Intervals.IntervalCollection

include("pairs.jl")
include("energy.jl")
include("duplex.jl")
include("intervals.jl")
include("trie.jl")

end
