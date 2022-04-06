#includes
using QuadraticBMS
using PyCall
using MATLABPlots
using JuMP
using MAT
using DataFrames
using CSV

dt=1
file = matopen("data/try3.mat")
Q_OCV = read(file,"QOCVfeasible")
q_OCV = vec(read(file,"qOCVfeasible"))
Q_COCV = read(file,"QPOCVfeasible")
q_COCV = vec(read(file,"qPOCVfeasible"))
Q_NOP = read(file,"QNOPfeasible")
q_NOP = vec(read(file,"qNOPfeasible"))
Q_POP = read(file,"QPOPfeasible")
q_POP = vec(read(file,"qPOPfeasible"))
Q_OOP = read(file,"QOOPfeasible")
q_OOP = vec(read(file,"qOOPfeasible"))
close(file)

capacity = 141550
cap_cell = capacity
oneC = 141550/3600
h = 0.36    # Cooling Coefficient
T_amb = 298    # Ambient Temperature
cpm = 1587    # Thermal Capacity
ic = [0.05,0,300]
N=4000
τ_ohm = 1

#Build Controller
controller_t = QuadraticBMS.buildController(ic,capacity,N,Q_OCV,q_OCV,Q_COCV,q_COCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,cpm,τ_ohm;max_T=333,pl_tol=0)

optimize!(controller_t)

I = value.(controller_t[:I])
DT_val = value.(controller_t[:DT])
T = value.(controller_t[:T])
t_val = cumsum(DT_val)
t_val = t_val.-t_val[1]



#df = DataFrame("# time [s]"=>t_val,"current [A]"=>I)
#CSV.write("optimal_input.csv",df)


