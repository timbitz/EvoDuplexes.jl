

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
   for_deletion = Nullable{Vector{Tuple{AbstractString,Int64,Int64}}}()
   if haskey( col.trees, int.first.seqname ) && length(col.trees[int.first.seqname]) > 0
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         for j in 1:length(dupvector)
            if dupvector[j] == int
               splice!( dupvector, j )
            end
         end
         if length(dupvector) == 0
            if isnull(for_deletion)
               for_deletion = Nullable(Vector{Tuple{AbstractString,Int64,Int64}})()
            end
            push!( for_deletion.value, (i.seqname, i.first, i.last) )
         end
      end
      if !isnull(for_deletion) && length(for_deletion.value) > 0
         for (n,f,l) in for_deletion.value
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

   for_deletion = Nullable{Vector{Tuple{AbstractString,Int64,Int64}}}()
   reset_first  = Nullable{Tuple{Interval,Int64}}()

   if haskey( col.trees, int.first.seqname )
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         hasintersect = true
         for j in 1:length(dupvector)
            if added_duplex && isoverlapping( int, dupvector[j], 0.75 )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  splice!( dupvector, j )
               else
                  remove_added = true
               end
            elseif isoverlapping( int, dupvector[j], 0.75 )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  i.last  = i.last  > int.first.last  ? i.last  : int.first.last
                  if int.first.first < i.first
                     reset_first = Nullable((i, int.first.first))
                  end
                  dupvector[j] = int
                  added_duplex = true
                  break
               else
                  return
               end
            elseif precedes( int, dupvector[j] )
               i.last  = i.last  > int.first.last  ? i.last  : int.first.last
               if int.first.first < i.first
                  reset_first = Nullable((i, int.first.first))
               end
               insert!( dupvector, j, int )
               added_duplex = true
               break
            end
         end
         if !added_duplex
            i.last  = i.last  > int.first.last  ? i.last  : int.first.last
            if int.first.first < i.first
               reset_first = Nullable((i, int.first.first))
            end
            push!( dupvector, int )
            added_duplex = true
         end

         # push for deletion
         if length(dupvector) == 0
            if isnull(for_deletion)
               for_deletion = Nullable(Vector{Tuple{AbstractString,Int64,Int64}}())
            end
            push!( for_deletion.value, (i.seqname, i.first, i.last) )
         end
      end
   
      if remove_added
         delete!(col, int)
         return
      end
   end
   # delete empty interval entries
   if !isnull(for_deletion)
      for (n,f,l) in for_deletion.value
         IntervalTrees.deletefirst!( col.trees[n], f, l )
         col.length -= 1
      end
   end

   if !isnull(reset_first) && 
      length(reset_first.value[1].metadata) > 0 &&
      !remove_added
      inter,first = reset_first.value
      IntervalTrees.deletefirst!( col.trees[inter.seqname], inter.first, inter.last )
      col.length -= 1
      inter.first = first
      push!( col, inter )
   end

   # push new interval and duplex
   if !hasintersect
      ival = Interval(int.first.seqname, int.first.first, int.first.last, int.first.strand, [ int ])
      push!(col, ival)
   end

   return
end



