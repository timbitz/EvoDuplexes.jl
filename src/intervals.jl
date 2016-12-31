

type DuplexInterval{T}
   first::Interval{T}  # fwd helical region 
   last::Interval{T}   # rev helical region
   duplex::RNADuplex   # rna duplex formed
end

# a DuplexCollection contains DuplexIntervals indexed by first interval
typealias DuplexCollection{T} IntervalCollection{Vector{DuplexInterval{T}}}


function isoverlapping{S,T}(x::Interval{S}, y::Interval{T}, identity::Float64)
   if isoverlapping( x, y )
      overlap  = min( x.last, y.last ) - max( x.first, y.first ) + 1
      smallest = min( x.last - x.first, y.last - y.first ) + 1
      @assert 0 <= overlap / smallest <= 1
      if overlap / smallest > identity
         return true
      end
   end
   return false
end
isoverlapping(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) && isoverlapping(x.last, y.last) ? true : false
isoverlapping(x::DuplexInterval, y::DuplexInterval, ident::Float64) = isoverlapping(x.first, y.first, ident) && isoverlapping(x.last, y.last, ident) ? true : false
precedes(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) ? Bio.Intervals.precedes(x.last, y.last) : Bio.Intervals.precedes(x.first, y.first)

function isordered{T}( col::DuplexCollection{T} )
   for i in col
      const dupvector = i.metadata
      for j in 2:length(dupvector)
         if !precedes( dupvector[j-1], dupvector[j] )
            return false
         end
      end
   end
   return true
end

function Base.delete!{T}(col::DuplexCollection{T}, int::DuplexInterval{T})
   fordel = Vector{Tuple{AbstractString,Int64,Int64}}()
   if haskey( col.trees, int.first.seqname ) && length(col.trees[int.first.seqname]) > 0
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         for j in 1:length(dupvector)
            if dupvector[j] == int
               splice!( dupvector, j )
            end
         end
         if length(dupvector) == 0
            push!( fordel, (i.seqname, i.first, i.last) )
         end
      end
      if length(fordel) > 0
         for (n,f,l) in fordel
            IntervalTrees.deletefirst!( col.trees[n], f, l )
            col.length -= 1
         end
      end
   end
end

function Base.push!{T}(col::DuplexCollection{T}, int::DuplexInterval{T})
   hasintersect = false
   added_duplex = false
   remove_added = false
   
   if haskey( col.trees, int.first.seqname )
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         hasintersect = true
         for j in 1:length(dupvector)
            if added_duplex && isoverlapping( int, dupvector[j], 0.5 )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  splice!( dupvector, j )
               else
                  remove_added = true
               end
            elseif isoverlapping( int, dupvector[j], 0.5 )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  i.first = i.first < int.first.first ? i.first : int.first.first
                  i.last  = i.last  > int.first.last  ? i.last  : int.first.last
                  dupvector[j] = int
                  added_duplex = true
                  break
               else
                  return
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
            col.length -= 1
         end
      end
   
      if remove_added
         delete!(col, int)
         return
      end
   end

   if !hasintersect
      ival = Interval(int.first.seqname, int.first.first, int.first.last, int.first.strand, [ int ])
      push!(col, ival)
   end

   return
end


#function Base.

#IntervalTrees.deletefirst!(ic.trees["c"], 15,19)
