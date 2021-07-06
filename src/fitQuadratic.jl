#Convenient function to fit a quadratic: xTQx + qTx
function fitQuadraticPSD(data::Array,y_true::Array)
    #Num Points, Num Features = size(data)
    np,nf = size(data)
    #Create Model
    model = JuMP.Model(Mosek.Optimizer)
    #Create Quadratic variable
    JuMP.@variable(model,Q[1:nf,1:nf])
    JuMP.@SDconstraint(model,Q>=0)
    JuMP.@variable(model,q[1:nf])
    #Create Variable for fitting
    JuMP.@variable(model,e[1:np])
    for i in 1:np
        JuMP.@constraint(model,e[i]==data[i,:]'*Q*data[i,:]+data[i,:]'*q-y_true[i])
    end
    #Minimize Error
    JuMP.@objective(model,Min,e'*e)
    return model
end

function fitQuadraticNSD(data::Array,y_true::Array)
    #Num Points, Num Features = size(data)
    np,nf = size(data)
    #Create Model
    model = JuMP.Model(Mosek.Optimizer)
    #Create Quadratic variable
    JuMP.@variable(model,Q[1:nf,1:nf])
    JuMP.@SDconstraint(model,Q<=0)
    JuMP.@variable(model,q[1:nf])
    #Create Variable for fitting
    JuMP.@variable(model,e[1:np])
    for i in 1:np
        JuMP.@constraint(model,e[i]==data[i,:]'*Q*data[i,:]+data[i,:]'*q-y_true[i])
    end
    #Minimize Error
    JuMP.@objective(model,Min,e'*e)
    return model
end
