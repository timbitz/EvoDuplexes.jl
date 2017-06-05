

type DuplexInterval{T}
   first::Interval{T}  # fwd helical region 
   last::Interval{T}   # rev helical region
   duplex::RNADuplex   # rna duplex formed
end

const ZERO_DUPLEX_INTERVAL = DuplexInterval(Interval("", 0, 0), Interval("", 0, 0), RNADuplex())

# a DuplexCollection contains DuplexIntervals indexed by first interval
type DuplexCollection{T}
   names::Dict{String,Dict{UInt64,Vector{DuplexInterval{T}}}}
   binsize::Int64
   length::Int64

   DuplexCollection()               = new(Dict{String,Dict{UInt64,Vector{DuplexInterval{T}}}}(), 50, 0 )
   DuplexCollection( binsize::Int ) = new(Dict{String,Dict{UInt64,Vector{DuplexInterval{T}}}}(), binsize, 0)
end

function Base.minimum{T}( col::DuplexCollection{T}; default::Float64=0.0 )
   minval = default
   minentry = Nullable{DuplexInterval{T}}()
   for n in keys(col.names)
      for b in keys(col.names[n])
         idx    = indmin(map(x->energy(x.duplex), col.names[n][b]))
         idxdg  = energy(col.names[n][b][idx].duplex)  
         if idxdg < minval
            minval   = idxdg
            minentry = Nullable(col.names[n][b][idx])
         end
      end
   end
   minentry
end

function Base.collect{T}( col::DuplexCollection{T}; minlength::Int=1, minenergy::Float64=0.0 )
   res = Vector{DuplexInterval{T}}()
   for n in keys(col.names)
      for b in sort(collect(keys(col.names[n])))
         for d in col.names[n][b]
            if npairs( d.duplex.path ) >= minlength &&
               energy( d.duplex ) <= minenergy
                push!( res, d )
            end
         end
      end
   end
   res
end


Base.length{T}( col::DuplexCollection{T}; keep::Bool=false ) = col.length > 0 && keep ? col.length : calculate_length!( col )

function calculate_length!{T}( col::DuplexCollection{T} )
   for n in keys(col.names)
      for b in keys(col.names[n])
         col.length += length(col.names[n][b])
      end
   end
   col.length
end

# hash by identity for speed
Base.hash(x::UInt64) = x

@inline bin{I <: Integer}( val::I, size::Int ) = (v = div(val, size)*size; return UInt64(v+1))

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

@inline isoverlapping(x::DuplexInterval, y::DuplexInterval) = isoverlapping(x.first, y.first) && 
                                                              isoverlapping(x.last, y.last) ? true : false
@inline isoverlapping(x::DuplexInterval, y::DuplexInterval, ident::Float64) = isoverlapping(x.first, y.first, ident) && 
                                                                              isoverlapping(x.last, y.last, ident) ? true : false

@inline precedes(x::DuplexInterval, y::DuplexInterval) = Bio.Intervals.precedes(x.first, y.first)
@inline Base.isless(x::DuplexInterval, y::DuplexInterval) = (x.first.first == y.first.first && x.first.last == y.first.last) ? 
                                                             Bio.Intervals.isless(x.last, y.last) : Bio.Intervals.isless(x.first, y.first)

Base.:<(x::DuplexInterval, y::DuplexInterval) = isless( x, y )

# return boolean of whether the duplex collection is
# properly ordered or not.
function isordered{T}( col::DuplexCollection{T} )
   for n in keys(col.names)
      const dict = col.names[n]
      prev = ZERO_DUPLEX_INTERVAL
      for b in sort(collect(keys(dict)))
         const vect = dict[b]
         for i in 1:length(vect)
            if !isless( prev, vect[i] )
               return false
            end
            prev = vect[i]
         end
      end
   end
   return true
end

function Base.delete!{T}(col::DuplexCollection{T}, int::DuplexInterval{T})
   const key = bin( int.first.first, col.binsize )
   if haskey( col.names, int.first.seqname ) && length(col.names[int.first.seqname][key]) > 0
      const vect = col.names[int.first.seqname][key]
      delete!( vect, int )
      if length(vect) == 0
         delete!( col.names[int.first.seqname], key )
      end
   end
end

function Base.delete!{T}( vect::Vector{DuplexInterval{T}}, int::DuplexInterval{T} )
   j = searchsortedfirst( vect, int )
   while j <= length(vect)
      if vect[j] == int
         deleteat!( vect, j )
      else
         j += 1
      end
   end
end

@inline function Base.push!{T}(col::DuplexCollection{T}, int::DuplexInterval{T})
   if !haskey(col.names, int.first.seqname)
      col.names[int.first.seqname] = Dict{UInt64,Vector{DuplexInterval{T}}}()
   end
   key = bin( int.first.first, col.binsize )
   # check bins before and after for better matches
   if !findbetter(col.names[int.first.seqname], key-0x01, int, true) &&
      !findbetter(col.names[int.first.seqname], key+0x01, int, false)
        pushinterval!(col.names[int.first.seqname], key, int)
   end
end

@inline function findbetter{T, K}(dict::Dict{K,Vector{DuplexInterval{T}}}, key::K, int::DuplexInterval{T}, left::Bool=true)
   !haskey( dict, key ) && (return false)
   const vect = dict[key]
   i = left ? searchsortedfirst( vect, int, lt=precedes ) : 1
   while i <= length(vect)
      if left && precedes( int, vect[j] )
         return false
      elseif isoverlapping( vect[i], int ) 
         if energy(vect[i].duplex) < energy(int.duplex)
            return true
         else
            deleteat!( vect, i )
            continue
         end
      end
      i += 1
   end
   if length(vect) == 0
      delete!( dict, key )
   end
   return false
end

function pushinterval!{T}(vect::Vector{DuplexInterval{T}}, int::DuplexInterval{T}; rate::Float64=0.99)
   j = searchsortedfirst( vect, int, lt=precedes )
   while j <= length(vect)
      if isoverlapping( int, vect[j], rate )
         if energy(int.duplex) < energy(vect[j].duplex)
            splice!( vect, j )
            j -= 1
         else
            return
         end
      end
      j += 1
   end
   insert!( vect, searchsortedfirst( vect, int ), int ) 
end

@inline function pushinterval!{T, K}(dict::Dict{K,Vector{DuplexInterval{T}}}, key::K, int::DuplexInterval{T}; rate::Float64=0.99)
   if !haskey( dict, key )
      dict[key] = Vector{DuplexInterval{T}}()
   end
   const vect = dict[key]
   pushinterval!( vect, int, rate=rate )
end


# This function 'stitches' two overlapping and compatible
# duplexes together into one.
@inline function _stitch{T}( a::DuplexInterval{T}, b::DuplexInterval{T}, max_bulge::Int, max_mis::Int )
   ret = Nullable{DuplexInterval{T}}()
   if isoverlapping( a, b )
      npairs_first = a.first.last - b.first.first + 1
      npairs_last  = b.last.last  - a.last.first  + 1
      afirst,alast = strings(a.duplex.path)
      bfirst,blast = strings(b.duplex.path)
      npairs = index_npairs(b.duplex.path, npairs_first, npairs_last)
      if abs(npairs_first - npairs_last) <= max_bulge &&
         length(bfirst) - npairs_first > 0 &&
         length(blast)  - npairs_last  > 0 &&
         npairs_first <= length(afirst) && npairs_first <= length(bfirst) &&
         npairs_last  <= length(alast)  && npairs_last  <= length(blast) &&
         afirst[end-npairs_first+1:end] == bfirst[1:npairs_first] && 
         alast[end-npairs_last+1:end]   == blast[1:npairs_last]

         spliced = deepcopy(a)
         push!( spliced.duplex, b.duplex.path[npairs+1:end] )
         spliced.first.last = b.first.last
         spliced.last.first = b.last.first
         aspliced, bspliced = strings(spliced.duplex)
              #println("$a\n$b\n")
              #println("$afirst : $alast")
              #println("$bfirst : $blast")
              #println("$npairs_first : $npairs_last : $npairs")
              #println("$(b.duplex.path[npairs+1:end])")
              #println("RESULT: $spliced")
                                                                  
         if spliced.first.last - spliced.first.first + 1 == length(aspliced) &&
            spliced.last.last  - spliced.last.first  + 1 == length(bspliced) &&
            energy(spliced.duplex) < energy(a.duplex) && 
            energy(spliced.duplex) < energy(b.duplex) &&
            nbulges(spliced.duplex.path) <= max_bulge && 
            nmismatches(spliced.duplex.path) <= max_mis
             ret = Nullable(spliced)
         end
      end
   end
   ret
end

stitch{T}( a::DuplexInterval{T}, b::DuplexInterval{T}, max_bulge::Int, max_mis::Int ) = a.first < b.first ? 
                                                                                             _stitch( a, b, max_bulge, max_mis ) : 
                                                                                             _stitch( b, a, max_bulge, max_mis )
#perform stitching on duplex collection
function stitch{T}(col::DuplexCollection{T})
   colret = deepcopy(col)
   for n in keys(col.names)
      res  = Vector{DuplexInterval{T}}()
      bins = sort(collect(keys(col.names[n])))
      for i in 1:length(bins)
         if i < length(bins) && bins[i+1] == bins[i]+1
            rec_stitch!( res, [ col.names[n][bins[i]]; col.names[n][bins[i+1]] ] )
         else
            rec_stitch!( res, col.names[n][bins[i]] )
         end
      end
      for r in res
         push!( colret, r )
      end
   end
   colret
end

function stitch!{T}( res::Vector{DuplexInterval{T}}, vect::Vector{DuplexInterval{T}} )
   for i in 1:length(vect)
      for j in i:searchsortedlast( vect, vect[i], lt=precedes )
         j == i && continue
         sval = stitch( vect[i], vect[j], 5, 5 )
         if !isnull(sval)
            pushinterval!( res, sval.value )
         end
      end
   end
end

function rec_stitch!{T}( res::Vector{DuplexInterval{T}}, vect::Vector{DuplexInterval{T}} )
   newvec = Vector{DuplexInterval{T}}()
   stitch!( newvec, vect )
   if length(newvec) > 0
      rec_stitch!( res, newvec )
      for r in newvec
         pushinterval!( res, r )
      end
   end
   res
end

