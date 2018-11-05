
struct DistanceForest
   forests::Vector{PyCall.PyObject}
   intervals::Vector{Tuple{Int,Int}}
   midpoints::Vector{Int}
   midvalues::Vector{Vector{Float64}}
   heuristic::Vector{Function}

   function DistanceForest( outratio::Float64=0.1, heuristic=[>, <, >, >, >] )
      inter, mids = distances()
      forests = [IsolationForest(contamination=outratio) for i in 1:length(mids)]
      new(forests, inter, mids, Vector{Vector{Float64}}(length(mids)), heuristic)
   end 
end

function train!( df::DistanceForest, array::Array{Float64,2}, dist_col::Int=1 )
   data = sortrows( array )
   training = data[:, 2:end]
   dist = data[:, 1]
   for i in 1:length(df.intervals)
      curint = df.intervals[i]
      start = searchsortedfirst(dist, curint[1])
      stop  = searchsortedlast(dist, curint[2])
      df.midvalues[i] = Float64[]
      for j in 1:size(training, 2)
         push!( df.midvalues[i], mean(training[start:stop, j]) )
      end
      if length(start:stop) > 1
         fit!( df.forests[i], training[start:stop, :] )
      end
   end
end

predict( df::DistanceForest, array::Vector{Float64} ) = predict( df, array' )

function predict( df::DistanceForest, array::Array{Float64,2} )
   dist = array[:, 1]
   assigned = map(x->searchsortedfirst( df.midpoints, x ), dist)
   outliers = Int[]

   function validate( row::Int, meanvalues::Vector{Float64} )
      for k in 2:size(array, 2)
         if !(df.heuristic[k-1]( array[row, k], meanvalues[k-1] )) 
            return false
         end
      end
      return true
   end

   for i in 1:length(df.midpoints)
      range = searchsorted( assigned, i )
      if length(range) >= 1
         res = ScikitLearn.predict( df.forests[i], array[range, 2:end] )
         for j in 1:length(res)
            if res[j] == -1 && validate( j + first(range) - 1, df.midvalues[i] )
               push!( outliers, j + first(range) - 1 )
            end
         end
      end
   end
   outliers
end

#function distances( b::Int=100, range::UnitRange=1:1 ) 
function distances( b::Int=2, range::UnitRange=2:12 )
   first = [(b^i, b^(i+1)) for i in range]
   last  = [(b^i+b^(i-1), b^(i+1)+b^i) for i in range]
   first[1] = (1, first[1][2])
   last[1]  = (1, last[1][2])
   ret = Vector{Tuple{Int,Int}}()
   retmean = Int[]
   midpoint( left, right ) = ((right - left) >> 1) + left
   for i in 1:length(first)
      push!(ret, first[i])
      push!(retmean, midpoint(first[i][1], first[i][2]))
      push!(ret, last[i])
      push!(retmean, midpoint(last[i][1], last[i][2]))
   end
   ret, retmean
end

