using JuMP
import Ipopt
using MATLABPlots


A = [-5 4;-3 2]./2.5
b = [1,2]
N=250
ys = zeros(N+1,2)
Δt = 1
ys[1,:] .= [1,2]
for n=1:N
     ys[n+1,:] .= ys[n,:]'*A*ys[n,:].+b
end
figure(1)
clf()
plot(ys)


#Control

model = Model(Ipopt.Optimizer)
@variable(model,u[i=1:N,1:2])
@variable(model,y[i=1:N+1,1:2])
fix(y[1,1],1)
fix(y[1,2],2)
for i in 2:N+1
    @constraint(model,y[i,:].==y[i-1,:]'*A*y[i-1,:].+b.+u[i-1,:])
end

@variable(model,e[i=1:N,1:2])
for i in 2:N+1
    @constraint(model,e[i-1,:].==(y[i,:].-[1,3]).^2)
end

@objective(model,Min,sum(e))

status = optimize!(model)
yv = value.(y)
hold_on()
plot(yv)
legend(["y1 uncontrolled","y2 uncontrolled","y1 controlled","y2 controlled"])

fix(y[1,1],5)
fix(y[1,2],6)
optimize!(model)
figure(2)
clf()
yv = value.(y)
plot(yv)
