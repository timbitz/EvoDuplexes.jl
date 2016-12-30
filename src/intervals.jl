

type DuplexInterval{T}
   first::Interval{T}  # fwd helical region 
   last::Interval{T}   # rev helical region
   duplex::RNADuplex   # rna duplex formed
end

# a DuplexCollection contains DuplexIntervals indexed by first interval
typealias DuplexCollection{T} IntervalCollection{Vector{DuplexInterval{T}}}

isoverlapping(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) && isoverlapping(x.last, y.last) ? true : false
precedes(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) ? Bio.Intervals.precedes(x.last, y.last) : Bio.Intervals.precedes(x.first, y.first)

function Base.push!(col::DuplexCollection, int::DuplexInterval)

   hasintersect = false
   added_duplex = false
   
   if haskey( col.trees, int.first.seqname )
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         hasintersect = true
         for j in 1:length(dupvector)
            if added_duplex && isoverlapping( int, dupvector[j] )
               splice!( dupvector, j )
            elseif isoverlapping( int, dupvector[j] )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  i.first = i.first < int.first.first ? i.first : int.first.first
                  i.last  = i.last  > int.first.last  ? i.last  : int.first.last
                  dupvector[j] = int
                  added_duplex = true
                  break
               end
            elseif precedes( int, dupvector[j] )
               i.first = i.first < int.first.first ? i.first : int.first.first
               i.last  = i.last  > int.first.last  ? i.last  : int.first.last
               insert!( dupvector, j, int )
               added_duplex = true
               break
            end
         end
         if !added_duplex
            i.first = i.first < int.first.first ? i.first : int.first.first
            i.last  = i.last  > int.first.last  ? i.last  : int.first.last
            push!( dupvector, int )
            added_duplex = true
         end

         if length(dupvector) == 0
            IntervalTrees.deletefirst!( col.trees[i.seqname], i.first, i.last )
         end
      end
   end

   if !hasintersect
      ival = Interval(int.first.seqname, int.first.first, int.first.last, int.first.strand, DuplexInterval[ int ])
      push!(col, ival)
   end

end


#function Base.

#IntervalTrees.deletefirst!(ic.trees["c"], 15,19)
