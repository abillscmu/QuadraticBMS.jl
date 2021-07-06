"""
    buildController()

Build a controller using the QuadraticBMS algorithm using IPOPT

Q Matrices should be fitted in the DFN Data Generation codes

Keyword arguments of max_T and pl_tol correspond to the constraints for maximum temperature and minimum anode potential, respectively


# Examples
```julia-repl
julia> buildController()
```
"""
function buildController(ic,capacity,dt,N,Q_OCV,q_OCV,Q_POCV,q_POCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,ρ,m,τ_ohm;max_T=333,pl_tol=0)
    #Construct IPOPT controller
    model = Model(Ipopt.Optimizer)
    #Build model matrices
    @variable(model,x[i=1:N,j=1:3])
    @variable(model,V[i=1:N])
    @variable(model,OCV[i=1:N])
    @variable(model,COCV[i=1:N])
    @variable(model,POP[i=1:N])
    @variable(model,NOP[i=1:N])
    @variable(model,OOP[i=1:N])
    @variable(model,I[i=1:N])
    @variable(model,ipl[i=1:N])
    @variable(model,T[i=1:N+1])
    @variable(model,z[i=1:N+1])
    @variable(model,DT[i=1:N])
    @variable(model,e[i=1:N])
    T_amb = 300;
    t_vec = collect(range(0,step=dt,length=N))
    #Initial conditions
    fix(z[1],ic[1])
    fix(OOP[1],ic[2])
    fix(T[1],ic[3])

    #Establish model constraints
    for i=1:N
        ##MODEL CONSTRAINTS: ALGEBRAIC AND DIFFERENTIAL
        #Algebraic Constraints
        @constraint(model,V[i]==OCV[i]-POP[i]-NOP[i]-OOP[i])
        @constraint(model,OCV[i]==x[i,:]'*Q_OCV*x[i,:]+q_OCV*x[i,:])
        @constraint(model,POCV[i]==x[i,:]'*Q_POCV*x[i,:]+q_POCV*x[i,:])
        @constraint(model,NOP[i]==x[i,:]'*Q_NOP*x[i,:]+q_NOP*x[i,:])
        @constraint(model,POP[i]==x[i,:]'*Q_POP*x[i,:]+q_POP*x[i,:])
        @constraint(model,ipl[i]==POCV[i]-OCV[i]+NOP[i])
        #Differential Constraints
        @constraint(model,OOP[i+1]==(-OOP[i]+x[i,:]'*Q_OOP*x[i,:]+q_OOP*x[i,:])*dt/τ_ohm+OOP[i])
        @constraint(model,z[i+1]==(I[i]/capacity)*dt+z[i])
        @constraint(model,T[i+1]==(h*(T[i]-T_amb)-I[i]*(V[i]-OCV[i]))*dt/ρ*m*cp+T[i])
        #Trivial Constraints (definitions)
        @constraint(model,x[i,3]==I[i])
        @constraint(model,x[i,2]==z[i])
        ##CONSTRAINTS ON OPERATION
        @constraint(model,ipl[i]>=pl_tol)
        @constraint(model,T[i+1]<=max_T)
        @constraint(model,e[i]==(z[i+1]-1)^2)
        
    end
    #Objective Function
    @objective(model,Min,sum(e))

end