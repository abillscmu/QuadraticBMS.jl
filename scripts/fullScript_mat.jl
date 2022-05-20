#includes
using QuadraticBMS
using MAT
using DataFrames
using CSV
using JuMP

dt=1
file = matopen("../matlab_out.mat")
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
N=1000
τ_ohm = 1

N_reftraj = 10
dt_reftraj=1.0

#Build Controller
controller_t = QuadraticBMS.buildController(ic,capacity,N,Q_OCV,q_OCV,Q_COCV,q_COCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,cpm,τ_ohm;max_T=333,pl_tol=0,top_SOC=0.2)

controller_reftraj = QuadraticBMS.buildController(ic,capacity,dt_reftraj,N_reftraj,Q_OCV,q_OCV,Q_POCV,q_POCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,cpm,τ_ohm)


function eval_controller_t()
    optimize!(controller_t)
    return value.(controller_t[:I])
end

function eval_controller_reftraj(new_ic)
    QuadraticBMS.reinitialize!(controller_reftraj,new_ic)
    optimize!(controller_reftraj)
    return value.(controller_refraj[:I])
end


#df = DataFrame("# time [s]"=>t_val,"current [A]"=>I)
#CSV.write("optimal_input.csv",df)


