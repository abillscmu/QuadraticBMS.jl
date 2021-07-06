function generateDFNData(cell,iapp::Array,N)
    # Generate Dfn Data using the cell defined by pybamm
    #First, want to get the OCV Curves
    x = collect(range(0,stop=1,length=N))
    #Convert the OCP to Julia readable
    cellDict = py"dict($cell)"
    nocv = cellDict["Negative electrode OCP [V]"]
    pocv = cellDict["Positive electrode OCP [V]"]
    #TODO: Deal with array OCV's
    @assert length(pocv)==2 "POCV must be a tuple of length 2"
    @assert length(nocv)==2 "NOCV must be a tuple of length 2"

    nocv_arr = nocv[2]
    pocv_arr = pocv[2]
    #Make sure we're dealing with arrays
    @assert nocv_arr isa AbstractArray "OCV must be an array"
    @assert pocv_arr isa AbstractArray "OCV must be an array"

    #_x present for possiblity of splines
    nocv_x = nocv_arr[:,1]
    pocv_x = pocv_arr[:,1]
    nocv_y = nocv_arr[:,2]
    pocv_y = pocv_arr[:,2]

    cathodeOCV = repeat(pocv_y,inner=[1,length(iapp)])
    anodeOCV = repeat(nocv_y,inner=[1,length(iapp)])

    α_pos = cellDict["Positive electrode charge transfer coefficient"] 
    α_neg = cellDict["Negative electrode charge transfer coefficient"]

    @assert α_pos==0.5 "α Must be 0.5 for this model to be applicable"
    @assert α_neg==0.5 "α Must be 0.5 for this model to be applicable"

    #Get Value of Exchange Current Density
    cathodeOP =  similar(cathodeOCV)
    anodeOP = similar(anodeOCV)

    n,m=size(cathodeOCV)
    cathode_css = pocv_x.*cellDict["Maximum concentration in positive electrode [mol.m-3]"]
    anode_css = nocv_x.*cellDict["Maximum concentration in negative electrode [mol.m-3]"]
    aFRT_pos = α_pos*FARADAY/(R*298)
    aFRT_neg = α_neg*FARADAY/(R*298)

    #Get Overpotentials
    for i in 1:n; for j in 1:m
        cathode_tree =  cellDict["Positive electrode exchange-current density [A.m-2]"](cathode_css[i],cathode_css[i],298)
        anode_tree = cellDict["Negative electrode exchange-current density [A.m-2]"](anode_css[i],anode_css[i],298)
        i_0p = recursive_traverse_tree(cathode_tree,cellDict)
        i_0n = recursive_traverse_tree(anode_tree,cellDict)
        I_ = iapp[j]
        cathodeOP[i,j]=(1/aFRT_pos).*asinh(I_./(2*i_0p));
        anodeOP[i,j]=(1/aFRT_neg).*asinh(I_./(2*i_0n));
    end;end
    return cathodeOCV,anodeOCV,cathodeOP,anodeOP,nocv_x,pocv_x
    
    
end


function recursive_traverse_tree(tree::PyObject,cell::Dict)
    #Function to recursively traverse a pybamm expression tree and evaluate the results
    
    #First, Figure out if it's constant
    iscons = tree.is_constant()
    if iscons
        #Figure out what to do here later lol
        val = tree.value
        return val
    else
        children = tree.children
        type = tree.name
        #Catch because julia=/= python
        if type=="**"
            type="^"
        elseif type in keys(cell)
            #It's a param, which may be a tree or may be a constant
            return recursive_traverse_tree(type,cell)
        end
        operator = Symbol(type)
        @assert length(children)==2
        a = children[1]
        b = children[2]
        child1 = recursive_traverse_tree(a,cell)
        child2 = recursive_traverse_tree(b,cell)
        return @eval $operator($child1,$child2)
    end
end

function recursive_traverse_tree(tree::String,cell::Dict)
    return cell[tree]
end



