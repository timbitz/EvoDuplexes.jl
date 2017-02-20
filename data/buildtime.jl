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
   bulge    = Vector{Float64}(length(steps))
   mismatch = Vector{Float64}(length(steps))
   bulgemat = Vector{Float64}(length(steps))
   dep      = Vector{Int}(length(steps))
   cnt      = Vector{Int}(length(steps))
   int      = Vector{Int}(length(steps))

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
         trie, build[j]   = @timed DuplexTrie{DNAAlphabet{2},UInt8}( seq[1:steps[j]], min(i, 5):i )
         cnt[j] = nodecount( trie.fwd )
         dep[j] = i
         res, trav[j]     = @timed traverse( trie, 1:10000 )
         int[j] = length(res)
         res, bulge[j]    = @timed traverse( trie, 1:10000, bulge_max=1 )
         res, mismatch[j] = @timed traverse( trie, 1:10000, mismatch_max=1 )
         res, bulgemat[j] = @timed traverse( trie, 1:10000, bulge_max=1, mismatch_max=1 )
      end
      if output
         for j in 1:length(steps)
            write( stream, string(i) * "\t" )
            write( stream, string(j) * "\t" )
            write( stream, string(steps[j]) * "\t" )
            write( stream, string(cnt[j]) * "\t" )
            write( stream, string(build[j]) * "\t" )
            write( stream, string(trav[j]) * "\t" )
            write( stream, string(bulge[j]) * "\t" )
            write( stream, string(mismatch[j]) * "\t" )
            write( stream, string(bulgemat[j]) * "\t" )
            write( stream, string(int[j]) * "\n" )
         end
      else
         df = join(df, DataFrame(length=steps, buildtime=build, travtime=trav, depth=dep), on=[:length, :buildtime, :travtime, :depth], kind = :outer)
      end
   end
   if !output
      buildplt = plot(df, x=:length, y=:buildtime, color=:depth, Geom.smooth)
      travplt  = plot(df, x=:length, y=:travtime,  color=:depth, Geom.smooth)

      return df,hstack(buildplt, travplt)
   else
      close(stream)
      close(io)
   end
end

function plot_profile( filename::String )
   df = readtable(filename, separator = '\t', header=false, names=[:depth, :step, :length, :nodecount, :buildtime, 
                                                                   :travtime, :travbulgetime, :travmistime, :travbulgemistime, :duplexes])

   buildplt = plot(df, x=:length, y=:buildtime, color=:depth, Geom.smooth,
         Coord.cartesian(ymin=0),
         Guide.xlabel("Sequence Length"),
         Guide.ylabel("Trie Build Time (sec)"))
   nodeplt  = plot(df, x=:length, y=:nodecount,  color=:depth, Geom.line,
         Coord.cartesian(ymin=0),
         Guide.xlabel("Sequence Length"),
         Guide.ylabel("Node Count"))
   travplt  = plot(df, x=:length, y=:travtime,  color=:depth, Geom.smooth,
         Coord.cartesian(ymin=0),
         Guide.xlabel("Sequence Length"),
         Guide.ylabel("Traversal Time (sec)"))
  duplextravplt = plot(df, x=:duplexes, y=:travtime, color=:depth, Geom.line,
         Coord.cartesian(ymin=0),
         Guide.xlabel("Returned Duplexes"),
         Guide.ylabel("Traversal Time (sec)"))

   draw(PDF("profile.pdf", 9inch, 8inch), vstack(hstack(buildplt, nodeplt), hstack(travplt, duplextravplt)))

   res = stack(df, [:buildtime, :travtime, :travbulgetime, :travmistime, :travbulgemistime])
   timeplot = plot(res, x=:length, y=:value, color=:variable, Geom.smooth, Coord.cartesian(ymin=0, fixed=1.0),
                   Guide.xlabel("Sequence Length"),
                   Guide.ylabel("Traversal Time (sec)"))
   draw(PDF("profile_time.pdf", 4inch, 4inch), timeplot)
end


