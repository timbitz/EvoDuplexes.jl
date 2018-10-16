
isgzipped( filename::String ) = hasextension( filename, "gz" )

function hasextension( filename::String, ext::String )
   restr = "\.$ext\$"
   re = Base.match(Regex(restr), filename)
   return re == nothing ? false : true
end

plug_write{S <: AbstractString}( io, str::S; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))
plug_write( io, str::Char; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))

tab_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\t' )
tab_write( io, str::Char ) = plug_write( io, str, plug='\t' )

end_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\n' )
end_write( io, str::Char ) = plug_write( io, str, plug='\n' )

function tab_write( io, dat::Tuple )
   for i in 1:length(dat)
      write( io, string(dat[i]) )
      if i < length(dat)
         write( io, '\t' )
      end
   end
   write( io, '\n' )
end

default_colors( n::Int ) = [Scale.color_continuous().f(p) for p in linspace(0, 1, n)]
default_colors( n::Int, opacity::Float64 ) = [(t = convert(Gadfly.RGBA, c);
                                              return Gadfly.RGBA(t.r, t.g, t.b, opacity))
                                              for c in default_colors(n)]

const DEF_COLORS = default_colors( 100, 1.0 )

random_rgb() = (rand(0:255), rand(0:255), rand(0:255))
random_rgb_str() = join(random_rgb(), ",")

rgba_to_int( col::Gadfly.RGBA ) = Int(floor(col.r*255)),Int(floor(col.g*255)),Int(floor(col.b*255))

function weighted_rgb( dg::Float64, maxenergy::Float64 )
   if dg < 0.0
      ind = min(Int(floor((abs(dg) / abs(maxenergy)) * 100)), 100)
      col = DEF_COLORS[ind]
   else 
      col = DEF_COLORS[1]
   end
   join(rgba_to_int(col), ",")
end

function writebed{T}( io, dup::DuplexInterval{T}, name::String; maxdistance::Int=100_000, maxenergy::Float64=-35.0 )
   const dist = dup.last.first - dup.first.last
   if dup.first.seqname == dup.last.seqname &&
      dist <= maxdistance
      tab_write( io, string(dup.first.seqname) )
      tab_write( io, string(dup.first.first-1) )
      tab_write( io, string(dup.last.last) )
      tab_write( io, "EvoDuplexes.jl:" * name )
      tab_write( io, string(min(1000, energy(dup.duplex) * -250)) )
      tab_write( io, dup.first.strand == STRAND_POS ? '+' : '-' )
      tab_write( io, string(dup.first.first-1) )
      tab_write( io, string(dup.last.last) )
      tab_write( io, weighted_rgb(energy(dup.duplex), maxenergy) )
      tab_write( io, "2" )
      tab_write( io, string(length(dup.first.first:dup.first.last)) * "," * string(length(dup.last.first:dup.last.last)) )
      end_write( io, string(0) * "," * string(dup.last.first-dup.first.first) )
   end
end

function writebed{T}( io, v::Vector{DuplexInterval{T}}, name::String; maxdistance::Int=100_000, maxenergy::Float64=-35.0 )
   for i in v
      writebed( io, i, name, maxdistance=maxdistance, maxenergy=maxenergy )
   end
end

