
# Format forced into left/right children, and additional markov probability matrix functionality added.
type PhyloNode
    label::String
    left::Nullable{PhyloNode}
    right::Nullable{PhyloNode}
    length::Real
    prob::Array{Float64,2}
end

const EMPTY_NODE = PhyloNode("", Nullable{PhyloNode}(), Nullable{PhyloNode}(), 0.0, Array{Float64}(0,0))

immutable PhyloTree
   root::PhyloNode
   order::Vector{String}
   index::Dict{String,Int}
end

const EMPTY_TREE = PhyloTree(EMPTY_NODE, String[], Dict{String,Int}())

### 'parsenewick' function based on: https://gist.github.com/porterjamesj/7672080
### License is CC-BY-NC
function parsenewick(newick::String)
    newick = rstrip(newick,';')
    newick = replace(newick,":","+")
    newick_expr = parse(newick)
    root  = parsenewick(newick_expr)
    order = collect(root)
    index = index_dict(order)
    PhyloTree(root, order, index)
end

parsenewick(newick::Symbol) = PhyloNode(string(newick), Nullable{PhyloNode}(), Nullable{PhyloNode}(), -1, Array{Float64}(0,0))

function parsenewick(newick::Expr)
    if newick.head == :tuple
        children = [parsenewick(child) for child in newick.args]
        left  = Nullable(children[1])
        right = Nullable(children[2])
        name = ""
        length = -1
    elseif newick.head == :call
        if newick.args[1] == :+
            # + indicates length
            length = newick.args[3]
            if typeof(newick.args[2]) == Expr
                if newick.args[2].head == :tuple
                    children = [parsenewick(child) for child in newick.args[2].args]
                    left  = Nullable(children[1])
                    right = Nullable(children[2])
                    name = ""
                elseif newick.args[2].head == :call && newick.args[2].args[1] == :*
                    # * indicates naming
                    name = string(newick.args[2].args[3])
                    children = [parsenewick(child) for child in newick.args[2].args[2].args]
                    left  = Nullable(children[1])
                    right = Nullable(children[2])
                end
            elseif typeof(newick.args[2]) == Symbol || typeof(newick.args[2]) == Int
                # tip node
                name = string(newick.args[2])
                left  = Nullable{PhyloNode}()
                right = Nullable{PhyloNode}()
            end
        elseif newick.args[1] == :*
            # bare * indicates a node with name but no length
            name = string(newick.args[3])
            children = [parsenewick(child) for child in newick.args[2].args]
            left  = Nullable(children[1])
            right = Nullable(children[2])
            length = -1
        end
    end
    PhyloNode(name,left,right,length, Array{Float64}(0,0))
end

# Tree collection
function Base.collect( root::PhyloNode )
   res = Vector{String}()
   function dfs_names!( vec::Vector{String}, node::PhyloNode )
      if !isnull(node.left) && !isnull(node.right)
         dfs_names!( vec, node.left.value )
         dfs_names!( vec, node.right.value )
      else
         push!( vec, node.label )
      end
      vec
   end
   dfs_names!( res, root )
end

function index_dict( labels::Vector{String} )
   ret = Dict{String,Int}()
   for i in 1:length(labels)
      ret[labels[i]] = i
   end
   ret
end

# Matrix setting and branch length modification functions
# Set P = exp(Q*b) for each node    
set_prob_mat!(tree::PhyloTree, Q::Array{Float64,2}) = set_prob_mat!(tree.root, Q)

function set_prob_mat!(node::PhyloNode, Q::Array{Float64,2})
   if !isnull(node.left) && !isnull(node.right) # internal node
      left  = set_prob_mat!(node.left.value,  Q)
      right = set_prob_mat!(node.right.value, Q)
   end
   if node.length > 0
      node.prob = expm(Q*node.length)
   end
   node
end

extend_branches!(tree::PhyloTree, multiplier::Float64) = extend_branches!(tree.root, multiplier)

function extend_branches!(node::PhyloNode, multiplier::Float64)
   if !isnull(node.left) && !isnull(node.right) # internal node
      left  = extend_branches!(node.left.value, multiplier)
      right = extend_branches!(node.right.value, multiplier)
   end
   if node.length > 0
      node.length *= multiplier
   end
   node
end

# Implementation of Felsenstein's pruning algorithm
# to calculate the L(D|T) = L(k) = Sigma_x(pi(x)*L(x))
# from the tree topology, branch lengths, instantaneous matrix (Q),
# and the data at site i
function likelihood( tree::PhyloTree, single::Vector{DNA}; 
                     background=GTR_SINGLE_PI, gapdenom::Float64=1.0 )
   function _likelihood( node::PhyloNode, single::Vector{DNA}, branchlength::Bool )
      if !isnull(node.left) && !isnull(node.right)
          if branchlength
             const left_res  = node.prob * _likelihood( node.left.value,  single, true )
             const right_res = node.prob * _likelihood( node.right.value, single, true )
          else 
             const left_res  = _likelihood( node.left.value,  single, true )
             const right_res = _likelihood( node.right.value, single, true )
          end
          #println("LEFT: $left_res .* RIGHT $right_res")
          return left_res .* right_res
      else
         #println("node: $(node.label)")
         const ind = tree.index[node.label]
         const v = isgap( single[ind] ) ? (return ones(4) / gapdenom) : zeros(4)
         v[ trailing_zeros(reinterpret(UInt8, single[ind])) + 1 ] = 1.0
         return node.prob * v
      end
   end
   sum(_likelihood( tree.root, single, false ) .* background)
end

function likelihood( tree::PhyloTree, first::Vector{DNA}, last::Vector{DNA};
                     background=GTR_PAIRED_PI )
   index( l::DNA, r::DNA ) = index( reinterpret(UInt8, l), reinterpret(UInt8, r) )
   function index( l::UInt8, r::UInt8 )
      const lidx = trailing_zeros(l) << 2
      const ridx = trailing_zeros(r)
      (lidx | ridx) + 1
   end
   function _likelihood( node::PhyloNode, first::Vector{DNA}, last::Vector{DNA}, branchlength::Bool )
      if !isnull(node.left) && !isnull(node.right)
          if branchlength
             const left_res  = node.prob * _likelihood( node.left.value,  first, last, true )
             const right_res = node.prob * _likelihood( node.right.value, first, last, true )
          else
             const left_res  = _likelihood( node.left.value,  first, last, true )
             const right_res = _likelihood( node.right.value, first, last, true )
          end
          return left_res .* right_res
      else
         const ind = tree.index[node.label]
         const v = isgap(first[ind]) || isgap(last[ind])  ? (return ones(16)) : zeros(16)
         v[ index(first[ind], last[ind]) ] = 1.0
         return node.prob * v
      end
   end
   sum(_likelihood( tree.root, first, last, false ) .* background)   
end
