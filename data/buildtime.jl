using Base.Test
using Base.Threads
using Bio.Seq
using Bio.Intervals
using Libz
using DataFrames
using Gadfly

importall Bio.Intervals

include("pairs.jl")
include("energy.jl")
include("duplex.jl")
include("intervals.jl")
include("io.jl")
include("trie.jl")

function trie_profile( seq, range, num_steps; raw_output="" )
   stepsize = div(length(seq), num_steps)
   steps    = collect(50:stepsize:length(seq))
   build    = Vector{Float64}(length(steps))
   trav     = Vector{Float64}(length(steps))
   dep      = Vector{Int}(length(steps))
   cnt      = Vector{Int}(length(steps))
   df = DataFrame(length=[], buildtime=[], travtime=[], depth=[])

   if raw_output != ""
      const io = open( raw_output, "w" )
      const stream = ZlibDeflateOutputStream( io )
      const output = true
   else
      const output = false
   end

   for i in range
      for j in 1:length(steps)
         trie, build[j] = @timed DuplexTrie{DNAAlphabet{2},UInt8}( seq[1:steps[j]], i:i )
         cnt[j] = nodecount( trie.fwd )
         dep[j] = i
         res, trav[j]   = @timed traverse( trie, 1:10000 )
      end
      if output
         for j in 1:length(steps)
            write( stream, string(i) * "\t" )
            write( stream, string(j) * "\t" )
            write( stream, string(steps[j]) * "\t" )
            write( stream, string(cnt[j]) * "\t" )
            write( stream, string(build[j]) * "\t" )
            write( stream, string(trav[j]) * "\n" )
         end
      else
         df = join(df, DataFrame(length=steps, buildtime=build, travtime=trav, depth=dep), on=[:length, :buildtime, :travtime, :depth], kind = :outer)
      end
   end
   if !output
      buildplt = plot(df, x=:length, y=:buildtime, color=:depth, Geom.smooth)
      travplt  = plot(df, x=:length, y=:travtime,  color=:depth, Geom.smooth)

      return df,hstack(buildplt, travplt)
   end
end
