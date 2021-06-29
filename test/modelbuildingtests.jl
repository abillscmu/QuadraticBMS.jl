#Test Fitting PSD Matrix
using QuadraticBMS
using LinearAlgebra
using JuMP

#Generate a dataset
x = rand(100,3)
A = rand(3,3)
APD = A'*A #Random PSD Matrix
a = rand(3)
y = zeros(100)
for (i,dy) in enumerate(y)
    y[i] = x[i,:]'*APD*x[i,:]+x[i,:]'*a
end

mod = QuadraticBMS.fitQuadraticPSD(x,y)
optimize!(mod)
Afit = value.(mod[:Q])
afit = value.(mod[:q])
#Check that the fitting Works (Afit should equal Apd)
@test Afit≈APD atol=1e-4
@test afit≈a atol=1e-4

x = rand(100,3)
ANSD = -APD #Random Negative Semidefinite Matrix
y = zeros(100)
for (i,dy) in enumerate(y)
    y[i] = x[i,:]'*ANSD*x[i,:]+x[i,:]'*a
end
mod = QuadraticBMS.fitQuadraticPSD(x,y)
optimize!(mod)
Afit = value.(mod[:Q])
ev = eigvals(Afit)
@test all(ev.>=0) #Test that the thing is PSD (should NOT be a good fit)
mod = QuadraticBMS.fitQuadraticNSD(x,y)
optimize!(mod)
Afit = value.(mod[:Q])
@test Afit≈ANSD atol=1e-4 #Should be a good fit
@test afit≈a atol=1e-4 
