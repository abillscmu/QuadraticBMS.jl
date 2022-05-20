"""
    buildController()

Build a controller using the QuadraticBMS algorithm using IPOPT

Q Matrices should be fitted in the DFN Data Generation codes

Keyword arguments of max_T and pl_tol correspond to the constraints for maximum temperature and minimum anode potential, respectively. This is intended to be an internal function for building a controller


"""
function buildController(ic,capacity,dt,N,Q_OCV,q_OCV,Q_POCV,q_POCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,c,τ_ohm;max_T=333,pl_tol=0)
    #THIS IS REF-TRAJ
    #Construct IPOPT controller
    model = Model(Ipopt.Optimizer)
    #Build model matrices
    @variable(model,x[i=1:N,j=1:3])
    @variable(model,V[i=1:N])
    @variable(model,OCV[i=1:N])
    @variable(model,POCV[i=1:N])
    @variable(model,POP[i=1:N])
    @variable(model,NOP[i=1:N])
    @variable(model,OOP[i=1:N+1])
    @variable(model,I[i=1:N])
    @variable(model,ipl[i=1:N])
    @variable(model,T[i=1:N+1])
    @variable(model,z[i=1:N+1])
    @variable(model,DT[i=1:N])
    @variable(model,e[i=1:N])
    @variable(model,e2[i=1:N-1])
    T_amb = 300;
    t_vec = collect(range(0,step=dt,length=N))
    oneC = capacity/3600
    #Initial conditions
    fix(z[1],ic[1])
    fix(OOP[1],ic[2])
    fix(T[1],ic[3])
    

    #Establish model constraints
    for i=1:N
        ##MODEL CONSTRAINTS: ALGEBRAIC AND DIFFERENTIAL
        #Algebraic Constraints
        @constraint(model,V[i]==OCV[i]-POP[i]-NOP[i]-OOP[i])
        @constraint(model,OCV[i]==x[i,:]'*Q_OCV*x[i,:]+q_OCV'*x[i,:])
        @constraint(model,POCV[i]==x[i,:]'*Q_POCV*x[i,:]+q_POCV'*x[i,:])
        @constraint(model,NOP[i]==x[i,:]'*Q_NOP*x[i,:]+q_NOP'*x[i,:])
        @constraint(model,POP[i]==x[i,:]'*Q_POP*x[i,:]+q_POP'*x[i,:])
        @constraint(model,ipl[i]==POCV[i]-OCV[i]+NOP[i])
        #Differential Constraints
        @constraint(model,OOP[i+1]==(-OOP[i]+x[i,:]'*Q_OOP*x[i,:]+q_OOP'*x[i,:])*dt/τ_ohm+OOP[i])
        @constraint(model,z[i+1]==(-I[i]/capacity)*dt+z[i])
        @constraint(model,T[i+1]==(-h*(T[i]-T_amb)-I[i]*(V[i]-OCV[i]))*dt/(c)+T[i])
        #Trivial Constraints (definitions)
        @constraint(model,x[i,3]==I[i])
        @constraint(model,x[i,2]==z[i])
        ##CONSTRAINTS ON OPERATION
        @constraint(model,ipl[i]>=pl_tol)
        @constraint(model,T[i+1]<=max_T)
        @constraint(model,e[i]==(z[i+1]-1)^2)
        #Restrict to 10C
        @constraint(model,I[i]<=10*oneC)
        @constraint(model,I[i]>=-10*oneC)
        @constraint(model,z[i]<=1)
        fix(x[i,1],1)
        
    end
    for i=2:N
        @constraint(model,e2[i-1]==I[i]-I[i-1])
    end
    #Objective Function
    @objective(model,Min,sum(e)+0.01*e2'*e2)

    return model

end


#THIS IS DT
function buildController(ic,capacity,N,Q_OCV,q_OCV,Q_POCV,q_POCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,c,τ_ohm;max_T=333,pl_tol=0,top_SOC=0.9)
    #Construct IPOPT controller
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 10000)
    #Build model matrices
    @variable(model,x[i=1:N,j=1:3])
    @variable(model,V[i=1:N])
    @variable(model,OCV[i=1:N])
    @variable(model,POCV[i=1:N])
    @variable(model,POP[i=1:N])
    @variable(model,NOP[i=1:N])
    @variable(model,OOP[i=1:N+1])
    @variable(model,I[i=1:N])
    @variable(model,ipl[i=1:N])
    @variable(model,T[i=1:N+1])
    @variable(model,z[i=1:N+1])
    @variable(model,DT[i=1:N])
    @variable(model,DOOP[i=1:N])
    #@variable(model,e[i=1:N])
    @variable(model,e2[i=1:N-1])
    T_amb = 300;
    #t_vec = collect(range(0,step=dt,length=N))
    #Initial conditions
    fix(z[1],ic[1])
    fix(OOP[1],ic[2])
    fix(T[1],ic[3])
    oneC = capacity/3600
    #Establish model constraints
    for i=1:N
        ##MODEL CONSTRAINTS: ALGEBRAIC AND DIFFERENTIAL
        #Algebraic Constraints
        @constraint(model,V[i]==OCV[i]-POP[i]-NOP[i]-OOP[i])
        @constraint(model,OCV[i]==x[i,:]'*Q_OCV*x[i,:]+q_OCV'*x[i,:])
        @constraint(model,POCV[i]==x[i,:]'*Q_POCV*x[i,:]+q_POCV'*x[i,:])
        @constraint(model,NOP[i]==x[i,:]'*Q_NOP*x[i,:]+q_NOP'*x[i,:])
        @constraint(model,POP[i]==x[i,:]'*Q_POP*x[i,:]+q_POP'*x[i,:])
        @constraint(model,ipl[i]==POCV[i]-OCV[i]+NOP[i])
        #Differential Constraints
        @NLconstraint(model,(OOP[i+1]-OOP[i])==DOOP[i]*DT[i])
        @constraint(model,DOOP[i]==(-OOP[i]+x[i,:]'*Q_OOP*x[i,:]+q_OOP'*x[i,:])/τ_ohm)
        @NLconstraint(model,(z[i+1]-z[i])==(-I[i]/capacity)*DT[i])
        @NLconstraint(model,(T[i+1]-T[i])==(-h*(T[i]-T_amb)-I[i]*(V[i]-OCV[i]))*DT[i]/(c))
        #Trivial Constraints (definitions)
        @constraint(model,x[i,3]==I[i])
        @constraint(model,x[i,2]==z[i])
        ##CONSTRAINTS ON OPERATION
        #@constraint(model,ipl[i]>=pl_tol)
        @constraint(model,T[i+1]<=max_T)
        #@constraint(model,e[i]==(z[i+1]-1)^2)
        #Restrict to 10C
        @constraint(model,I[i]<=0)
        @constraint(model,I[i]>=-10*oneC)
        @constraint(model,DT[i]<=1)
        @constraint(model,DT[i]>=0.1)
        fix(x[i,1],1)
        
    end
    for i=2:N-1
        @constraint(model,z[i]<=1)
    end
    for i=2:N
        @constraint(model,e2[i-1]==(I[i]-I[i-1])^2)
    end
    #Objective Function
    @constraint(model,z[end]>=top_SOC)
    @NLobjective(model,Min,sum(DT[i] for i in 1:N)/1000+sum(e2[i] for i in 1:N-1)/10000)

    return model

end

function buildControllerRungeKutta(ic,capacity,N,Q_OCV,q_OCV,Q_POCV,q_POCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,c,τ_ohm;max_T=333,pl_tol=0)
    #Construct IPOPT controller
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 10000)
    #Build model matrices
    @variable(model,x[i=1:N,j=1:3])
    @variable(model,V[i=1:N])
    @variable(model,OCV[i=1:N])
    @variable(model,POCV[i=1:N])
    @variable(model,POP[i=1:N])
    @variable(model,NOP[i=1:N])
    @variable(model,OOP[i=1:N+1])
    @variable(model,I[i=1:N])
    @variable(model,ipl[i=1:N])
    @variable(model,T[i=1:N+1])
    @variable(model,z[i=1:N+1])
    @variable(model,DT[i=1:N])
    @variable(model,DOOP[i=1:N])
    #@variable(model,e[i=1:N])
    @variable(model,e2[i=1:N-1])
    T_amb = 300;
    #t_vec = collect(range(0,step=dt,length=N))
    #Initial conditions
    fix(z[1],ic[1])
    fix(OOP[1],ic[2])
    fix(T[1],ic[3])
    oneC = capacity/3600
    #Establish model constraints
    for i=1:N
        ##MODEL CONSTRAINTS: ALGEBRAIC AND DIFFERENTIAL
        #Algebraic Constraints
        @constraint(model,V[i]==OCV[i]-POP[i]-NOP[i]-OOP[i])
        @constraint(model,OCV[i]==x[i,:]'*Q_OCV*x[i,:]+q_OCV'*x[i,:])
        @constraint(model,POCV[i]==x[i,:]'*Q_POCV*x[i,:]+q_POCV'*x[i,:])
        @constraint(model,NOP[i]==x[i,:]'*Q_NOP*x[i,:]+q_NOP'*x[i,:])
        @constraint(model,POP[i]==x[i,:]'*Q_POP*x[i,:]+q_POP'*x[i,:])
        @constraint(model,ipl[i]==POCV[i]-OCV[i]+NOP[i])
        #Differential Constraints
            @NLconstraint(model,(OOP[i+1]-OOP[i])==DOOP[i]*DT[i])
            @constraint(model,DOOP[i]==(-OOP[i]+x[i,:]'*Q_OOP*x[i,:]+q_OOP'*x[i,:])/τ_ohm)
            @NLconstraint(model,(z[i+1]-z[i])==(-I[i]/capacity)*DT[i])
            @NLconstraint(model,(T[i+1]-T[i])==(-h*(T[i]-T_amb)-I[i]*(V[i]-OCV[i]))*DT[i]/(c))
        #Trivial Constraints (definitions)
        @constraint(model,x[i,3]==I[i])
        @constraint(model,x[i,2]==z[i])
        ##CONSTRAINTS ON OPERATION
        @constraint(model,ipl[i]>=pl_tol)
        @constraint(model,T[i+1]<=max_T)
        #@constraint(model,e[i]==(z[i+1]-1)^2)
        #Restrict to 10C
        @constraint(model,I[i]<=0)
        @constraint(model,I[i]>=-20*oneC)
        @constraint(model,DT[i]<=0.25)
        @constraint(model,DT[i]>=0.025)
        fix(x[i,1],1)
        
    end
    for i=2:N-1
        @constraint(model,z[i]<=1)
    end
    for i=2:N
        @constraint(model,e2[i-1]==(I[i]-I[i-1])^2)
    end
    #Objective Function
    @constraint(model,z[end]>=0.7)
    @NLobjective(model,Min,sum(DT[i] for i in 1:N)/1000)

    return model

end

function reinitialize!(model,new_ic)
    z = model[:z]
    OOP = model[:OOP]
    T = model[:T]
    fix(z[1],new_ic[1])
    fix(OOP[1],new_ic[2])
    fix(T[1],new_ic[3])
    return model
end
