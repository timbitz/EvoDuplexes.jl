

type DuplexInterval{T}
   first::Interval{T}  # fwd helical region 
   last::Interval{T}   # rev helical region
   duplex::RNADuplex   # rna duplex formed
end

# a DuplexCollection contains DuplexIntervals indexed by first interval
typealias DuplexCollection{T} IntervalCollection{Vector{DuplexInterval{T}}}

fixed_interval{I <: Integer}( val::I, size::Int ) = (v = div(val, size)*size; return (v+1, v+size))

function isoverlapping{S,T}(x::Interval{S}, y::Interval{T}, identity::Float64)
   if isoverlapping( x, y )
      overlap  = min( x.last, y.last ) - max( x.first, y.first ) + 1
      smallest = min( x.last - x.first, y.last - y.first ) + 1
      if overlap / smallest >= identity
         return true
      end
   end
   return false
end

isoverlapping(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) && 
                                                      isoverlapping(x.last, y.last) ? true : false
isoverlapping(x::DuplexInterval, y::DuplexInterval, ident::Float64) = isoverlapping(x.first, y.first, ident) && 
                                                                      isoverlapping(x.last, y.last, ident) ? true : false


Base.isless(x::DuplexInterval, y::DuplexInterval) = (x.first.first == y.first.first && x.first.last == y.first.last) ? 
                                                    Bio.Intervals.isless(x.last, y.last) : Bio.Intervals.isless(x.first, y.first)

Base.:<(x::DuplexInterval, y::DuplexInterval) = isless( x, y )

# return boolean of whether the duplex collection is
# properly ordered or not.
function isordered{T}( col::DuplexCollection{T} )
   for i in col
      const dupvector = i.metadata
      for j in 2:length(dupvector)
         if !isless( dupvector[j-1], dupvector[j] )
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
         j = 1
         while j <= length(dupvector)
            if dupvector[j] == int
               splice!( dupvector, j )
            else
               j += 1
            end
         end
         if length(dupvector) == 0
            if isnull(for_deletion)
               for_deletion = Nullable(Vector{Tuple{AbstractString,Int64,Int64}}())
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

function Base.push!{T}(col::DuplexCollection{T}, int::DuplexInterval{T}; rate::Float64=0.75, size::Int=25)
   hasintersect = false
   added_duplex = false
   remove_added = false

   for_deletion = Nullable{Vector{Tuple{AbstractString,Int64,Int64}}}()

   if haskey( col.trees, int.first.seqname )
      for i in intersect(col, int.first)
         const dupvector = i.metadata
         hasintersect = true
         j = 1
         while j <= length(dupvector)
            if added_duplex && isoverlapping( int, dupvector[j], rate )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  splice!( dupvector, j )
                  j -= 1
               else
                  remove_added = true
               end
            elseif isoverlapping( int, dupvector[j], rate )
               if energy(int.duplex) < energy(dupvector[j].duplex)
                  splice!( dupvector, j )
                  push!( col, int, rate=rate, size=size )
                  added_duplex = true
                  break
               else
                  return
               end
            elseif !added_duplex && isless( int, dupvector[j] )
               insert!( dupvector, j, int )
               added_duplex = true
            end
            j += 1
         end
         if !added_duplex
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

   # push new interval and duplex
   if !hasintersect
      bounds = fixed_interval( int.first.first, size )
      ival = Interval(int.first.seqname, bounds[1], bounds[2], int.first.strand, [ int ])
      push!( col, ival )
   end

   return
end



