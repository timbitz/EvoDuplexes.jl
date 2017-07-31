# Code based on Automa.jl multi-line 'FASTA' example:
# https://github.com/BioJulia/Automa.jl/blob/master/example/fasta.jl 

const re = Automa.RegExp

# Create a machine of MAF
const record_machine = (function ()

    lf          = re"\n"
    newline     = re"\r?" * lf
    lf.actions[:enter]          = [:countline]

    anything    = re.rep(re.diff( re".", newline ))
    header      = re.cat(re" ?#", anything, newline)

    float       = re.cat(re"[-+]?[0-9]*", re.opt("\."), re.opt(re"[0-9]*"))
    float.actions[:enter]       = [:mark]
    float.actions[:exit]        = [:score]

    sln         = re.rep1(re.space())
    emptyline   = re.cat(re" *", newline)

    scoreline   = re.cat(re" ?a score=", float, re.alt(newline, re.cat(re" ", anything, newline)))
    #scoreline.actions[:enter]   = [:setstart]

    organism    = re"[0-9A-Za-z\.\_\-]+"
    organism.actions[:enter]    = [:mark]
    organism.actions[:exit]     = [:organism]

    position    = re"[0-9]*"
    position.actions[:enter]    = [:mark]
    position.actions[:exit]     = [:position]

    size        = re"[0-9]*"
    size.actions[:enter]        = [:mark]
    size.actions[:exit]         = [:size]
 
    strand      = re"[+-]"
    strand.actions[:enter]      = [:mark]
    strand.actions[:exit]       = [:strand]

    sequence    = re"[ATGCNatgcn-]+"
    sequence.actions[:enter]    = [:mark]
    sequence.actions[:exit]     = [:sequence]

    entry       = re.cat(re" ?s", sln, organism, sln, position, sln, size, sln, strand, sln, re"[0-9]*", sln, sequence, newline)
    entry.actions[:enter]       = [:setpos]
    entry.actions[:exit]        = [:entry]

    fullheader  = re.cat(re.rep(header), re.rep(emptyline))
    fullheader.actions[:exit]   = [:anchor, :mark]

    record      = re.cat(re.opt(fullheader), re.opt(scoreline), entry, re.rep(entry), re.rep1(emptyline))
    record.actions[:enter]      = [:anchor]
    record.actions[:exit]       = [:escape]

    return Automa.compile(re.rep(record))
end)()

# It is useful to visualize the state machine for debugging.
# write("maf.dot", Automa.machine2dot(record_machine))
# run(`dot -Tsvg -o maf.svg maf.dot`)

# Bind Julia code to each action name (see the `parse_fasta` function defined below)

maf_actions_stream = Dict(
            :countline   => :(linenum += 1),
            :mark        => :(mark = p),
            :score => quote
                if !bad_record
                   range = (mark:p-1)
                   score = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                end
                mark = 0
            end,
            :organism => quote
                if !bad_record
                   range = (mark:p-1)
                   organism = mark == 0 ? "" : convert(String, data[range])
                end
                mark = 0
            end,
            :position => quote
                if !bad_record
                   range = (mark:p-1)
                   position = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                end
                mark = 0
            end,
            :size => quote
                if !bad_record
                   range = (mark:p-1)
                   size = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                end
                mark = 0
            end,
            :strand => quote
                if !bad_record
                   range = (mark:p-1)
                   strand = mark == 0 ? true : Char(data[last(range)]) == '+'
                end
                mark = 0
            end,
            :sequence => quote
                if !bad_record
                   range = (mark:p-1)
                   sequence = mark == 0 ? "" : chomp(String(data[range]))
                end
                mark = 0
            end,
            :setpos => :(stream.position = p; mark_cs = cs), #println("Set stream position $p at cs $cs")),
            :entry  => quote
                try
                   if first
                      const len    = parse(Int, size)
                      const pscore = parse(Float64, score)
                      const pos    = parse(Int, position)
                      org,chr = parsename( organism, reader.keepname )
                      if (len < reader.minlen || pscore < reader.minscore) ||
                         (regionbool && !IntervalTrees.hasintersection( regions.trees[chr], pos ) )
                         bad_record = true
                      end
                   end
                   first = false
                   if !bad_record
                      org,chr = parsename( organism, reader.keepname )
                      push!(record.species, MAFSpecies(org, chr, parse(Int, position), strand, DNASequence(sequence)))
                   end
                catch
                   println(STDERR, "size: $size and score: $score")
                   println(STDERR, "Bad Record at $p, setting p to $(stream.position)")
                   p = stream.position
                   @escape
                end
            end,
            :countline   => :(linenum += 1),
            :anchor      => :(stream.anchor = p), #; println("Anchoring at $p")),
            :escape      => :(found_record=true; @escape))

# Define a type to store a FASTA record.
type MAFSpecies
    name::String
    chr::String
    position::Int
    strand::Bool
    sequence::BioSequence{DNAAlphabet{4}}
end

Base.:(==)( a::MAFSpecies, b::MAFSpecies ) = a.name == b.name && a.position == b.position && 
                                             a.strand == b.strand && a.sequence == b.sequence ? true : false

type MAFRecord
   species::Vector{MAFSpecies}

   MAFRecord() = new(Vector{MAFSpecies}())
end

Base.endof( maf::MAFRecord ) = endof(maf.species)
Base.getindex( maf::MAFRecord, ind ) = maf.species[ind]

function Base.reverse!( maf::MAFRecord )
   for s in maf.species
      reverse!( s.sequence )
   end
end

function reverse_complement!( mblock::MAFRecord, strandflip::Bool=true )
   for s in mblock.species
      Bio.Seq.reverse_complement!( s.sequence )
      if strandflip
         s.strand = s.strand ? false : true
      end
   end
end

Base.length( maf::MAFRecord ) = length(maf.species)


function deletegaps!( mblock::MAFRecord )
   str = convert(String, mblock.species[1].sequence)
   ind = search(str, r"-+")
   if length(ind) > 0
      for i in mblock.species
         deleteat!( i.sequence, ind )
      end
      deletegaps!( mblock )
      return
   else
      return
   end
end

type MAFReader{T<:BufferedInputStream}
    stream::T
    cs::Int
    p_eof::Int
    linenum::Int
    minlen::Int
    minscore::Float64
    keepname::Bool
end

function MAFReader{BS <: BufferedInputStream}(input::BS, keepname::Bool=false)
    hits_eof = fillbuffer!(input) == 0
    p_eof = hits_eof ? input.available : -1
    return MAFReader{BS}(input, record_machine.start_state, p_eof, 1, 0, typemin(Float64), keepname)
end

Base.done(reader::MAFReader) = reader.cs <= 0

iscomment(line::String) = (length(line) >= 1 && line[1] == '#') || (length(line) >= 2 && line[2] == '#')
isheader(line::String)  = iscomment(line)

function discard_header!(reader::MAFReader)
    for l in eachline(reader.stream)
       reader.linenum += 1
       if isheader(l)
           continue
       end
       anchor!(reader.stream)
       hits_eof = fillbuffer!(reader.stream) == 0
       if hits_eof
          BufferedStreams.shiftdata!(reader.stream)
          reader.p_eof = reader.stream.available
       end
       reader.linenum -= 1
       return reader
    end
    reader.linenum -= 1
    reader
end

function parsename( name::String, keepname::Bool )
   if keepname
      return name,name
   else
      spl = split(name, '.')
      length(spl) < 2 && error("Must have a . in the name if MAFReader( ..., keepname=true ): $name")
      return convert(String, spl[1]),convert(String, spl[2])
   end 
end

const EMPTY_INTCOL = IntervalCollection{Void}()

@eval function Base.read!(reader::MAFReader, record::MAFRecord; regionbool::Bool=false, regions::IntervalCollection=EMPTY_INTCOL)
    cs      = reader.cs
    stream  = reader.stream
    data    = stream.buffer
    p       = stream.position
    p_end   = stream.available
    p_eof   = reader.p_eof
    linenum = reader.linenum
    mark = 0
    mark_cs = 0

    scoreline = organism = sequence = ""
    size = position = score = ""
    bad_record   = false
    found_record = false
    rebuffered   = false
    retry        = false
    first        = true

    if cs < 0
        error("the reader is in error state")
    elseif cs == 0
        error("the reader is finished")
    end

    while true

        p = p < 1 ? stream.anchor : p
        #error("p:$p p_end:$p_end p_eof:$p_eof cs:$cs anchor:$(stream.anchor) avail:$(stream.available) \"$(map(Char, data[1:10]))\"")

        $(Automa.generate_exec_code(record_machine, actions=maf_actions_stream, code=:goto, check=false))

        reader.cs              = cs
        reader.p_eof           = p_eof
        reader.linenum = linenum
      
        if found_record && !bad_record
           reader.stream.position = p
           return record
        elseif found_record && bad_record
           bad_record = false
           found_record = false
           first = true
           empty!(record.species)
           continue
        elseif (p == p_eof+1 && data[p_eof] == 10)
           #println("p==p_eof+1 p:$p p_eof:$p_eof cs:$cs \"$(Char(data[p-1]))\" $record")
           if !bad_record && organism != ""
              org,chr = parsename( organism, reader.keepname )
              spec = MAFSpecies(org, chr, parse(Int, position), strand, DNASequence(sequence))
              length(record.species) > 0 && !(record.species[end] == spec) && push!(record.species, spec)
           end
           reader.cs = 0
           return record
        elseif p > p_eof ≥ 0
         #  println("p:$p p_eof:$p_eof cs:$cs \"$(Char(data[p-1]))\" $record")
           error("incomplete MAF input on line ", linenum)
        elseif p == p_end+1
            #println(STDERR, "Rebuffer p=p_end+1 p:$p p_eof:$p_eof p_end:$p_end cs:$cs")
            hits_eof = BufferedStreams.fillbuffer!(reader.stream) == 0
            if hits_eof
               reader.p_eof = p_eof = p_end
            else
               p     = stream.position
               p_end = stream.available
            end
            #println("Rebuffering to $(map(Char, data[p:p+10])) at cs $mark_cs")
            cs = record_machine.start_state
        elseif cs < 0
           if !retry && p <= 4
              retry = true
              p = 1
              cs = record_machine.start_state
           else
              println(String(data))
              println("cs:$cs p:$p p_end:$p_end $record anchor:$(stream.anchor) avail:$(stream.available) \"$(map(Char, data[1:p+10]))\"")
              error("MAF file format error on line ", linenum)
           end
        else
           error("Shouldn't happen cs:$cs p:$p p_eof:$p_eof $record")
        end

    end
end


function stitch!( a::MAFRecord, b::MAFRecord, index::Dict{String,Int} )
   const aref = a.species[1]
   const bref = b.species[1]
   sort!( a.species, by=x->index[x.name] )
   sort!( b.species, by=x->index[x.name] )
   aref.name != bref.name && error("Invalid stitching of alignment blocks with two different references $(aref.name) & $(bref.name)!!")
   if aref.position+length(aref.sequence) == bref.position
       i,j = 1,1
       #println(" $( map(x->a.species[x].name, 1:length(a)) ) && \n $( map(x->b.species[x].name, 1:length(b)) )")
       while i <= length(a.species) && j <= length(b.species)
          #println("$( a.species[i].name ) vs. $( b.species[j].name ) for $i and $j")
          if a.species[i].name == b.species[j].name
             a.species[i].sequence *= b.species[j].sequence
             i += 1
             j += 1
          elseif index[a.species[i].name] < index[b.species[j].name]
             a.species[i].sequence *= dna"-" ^ length(bref.sequence)
             i += 1
          else
             b.species[j].sequence = dna"-" ^ length(aref.sequence) * b.species[j].sequence
             insert!( a.species, i, b.species[j] )
             i += 1
             j += 1
          end
       end
       return true
   else
      return false
   end
end

const MAFInterval = Interval{MAFRecord}
const MAFCollection = IntervalCollection{MAFRecord}

function readmaf!( reader::MAFReader, order::Dict{String,Int}; minspecies::Int=1, minlength::Int=1, minscore::Float64=-Inf )
   maf = MAFCollection()
   readmaf!( maf, reader, order, minspecies=minspecies, minlength=minlength, minscore=minscore )
end

function readmaf!( maf::MAFCollection, reader::MAFReader, order::Dict{String,Int}; minspecies::Int=1, minlength::Int=1, minscore::Float64=-Inf )
   #reader.minlen   = minlength
   #reader.minscore = minscore
   prev = MAFRecord()
   while !done( reader )
      rec = MAFRecord()
      read!( reader, rec )
      if length(rec.species) >= minspecies
         try
            deletegaps!(rec)
         catch
            println(rec)
            error()
         end
         if length(prev) > 0
            if !(stitch!( prev, rec, order ))
               const ref = prev.species[1]
               if length(ref.sequence) > minlength
                  push!( maf, MAFInterval(ref.chr, ref.position, ref.position+length(ref.sequence)-1, ref.strand ? '+' : '-', prev) )
               end
               prev = rec
            end
         else
            prev = rec
         end
      end
   end
   if length(prev) > 0
      const ref = prev.species[1]
      push!( maf, MAFInterval(ref.chr, ref.position, ref.position+length(ref.sequence)-1, ref.strand ? '+' : '-', prev) )
   end
   maf
end

