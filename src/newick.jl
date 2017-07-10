### Code from: https://gist.github.com/porterjamesj/7672080
### License is CC-BY-NC

type PhyloNode
    label::String
    children::Vector{PhyloNode}
    length::Real
end

function parsenewick(newick::String)
    newick = rstrip(newick,';')
    newick = replace(newick,":","+")
    newick_expr = parse(newick)
    return parsenewick(newick_expr)
end

parsenewick(newick::Symbol) = PhyloNode(string(newick), PhyloNode[], -1)

function parsenewick(newick::Expr)
    if newick.head == :tuple
        children = [parsenewick(child) for child in newick.args]
        name = ""
        length = -1
    elseif newick.head == :call
        if newick.args[1] == :+
            # + indicates length
            length = newick.args[3]
            if typeof(newick.args[2]) == Expr
                if newick.args[2].head == :tuple
                    children = [parsenewick(child) for child in newick.args[2].args]
                    name = ""
                elseif newick.args[2].head == :call && newick.args[2].args[1] == :*
                    # * indicates naming
                    name = string(newick.args[2].args[3])
                    children = [parsenewick(child) for child in newick.args[2].args[2].args]
                end
            elseif typeof(newick.args[2]) == Symbol || typeof(newick.args[2]) == Int
                # tip node
                name = string(newick.args[2])
                children = PhyloNode[]
            end
        elseif newick.args[1] == :*
            # bare * indicates a node with name but no length
            name = string(newick.args[3])
            children = [parsenewick(child) for child in newick.args[2].args]
            length = -1
        end
    end
    PhyloNode(name,convert(Vector{PhyloNode},children),length)
end
