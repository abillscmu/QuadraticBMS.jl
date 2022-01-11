#includes
using QuadraticBMS
using PyCall
using MATLABPlots
using JuMP
pybamm = pyimport("pybamm")
pandas = pyimport("pandas")

iapp = collect(range(-30,stop=0,length=100))

#dfn data generation
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Ramadass2004)



cellDict = py"dict($param)"
println("Generating Data...")
cathodeOCV,anodeOCV,cathodeOP,anodeOP,OOP = QuadraticBMS.generateDFNData(cellDict,iapp,25)
println("Data Generation Complete!")
n,m = size(cathodeOCV)
soc = collect(range(0,stop=1,length=n))
reverse!(cathodeOCV,dims=1)
reverse!(cathodeOP,dims=1)
cocv = cathodeOCV[:,1]
aocv = anodeOCV[:,1]
using MATLABPlots
plotOptions = Dict("LineWidth"=>3)

#Plot OCV Curves
figure(1)
clf()
plot(soc,cocv,"r",options=plotOptions)
hold_on()
plot(soc,aocv,"b",options=plotOptions)
plot(soc,cocv.-aocv,"k",options=plotOptions)
xlabel("Stoichiometry")
ylabel("Voltage")
graphicsOptions = Dict(
    "FontName"=>"Palatino",
    "FontWeight"=>"Bold",
    "FontSize"=>24
)
legend(["Positive OCV","Negative OCV","Total OCV"])
setgca(graphicsOptions)

IAPP = repeat(iapp',inner=[n,1])
SOC = repeat(soc,inner=[1,m])
#Plot OverPotentials
figure(2)
subplot(2,1,1)
surf(IAPP[1:end-1,:],SOC[1:end-1,:],cathodeOP[1:end-1,:])
xlabel("Current")
ylabel("SOC")
zlabel("Overpotential (positive)")
figure(2)
subplot(2,1,2)
surf(IAPP[1:end-1,:],SOC[1:end-1,:],anodeOP[1:end-1,:])
xlabel("Current")
ylabel("SOC")
zlabel("Overpotential (negative)")

makevec(mat) = reshape(mat,(length(mat),1))

iapp = makevec(IAPP)
soc = makevec(SOC)
cathodeOCV = makevec(cathodeOCV)
anodeOCV = makevec(anodeOCV)

x = hcat(ones(length(iapp)),soc,iapp)
OCV = cathodeOCV.-anodeOCV

#Get Everything
get_val(sym::Symbol,mod::Model) = value.(mod[sym])

#OCV
mod_QOCV = QuadraticBMS.fitQuadraticNSD(x,OCV)
optimize!(mod_QOCV)
Q_OCV = value.(mod_QOCV[:Q])
q_OCV = value.(mod_QOCV[:q])

mod_COCV = QuadraticBMS.fitQuadraticNSD(x,cathodeOCV)
optimize!(mod_COCV)
Q_COCV = value.(mod_COCV[:Q])
q_COCV = value.(mod_COCV[:q])

isnotinf = !(isinf)
isnotnan = !(isnan)
cathodeOP = makevec(cathodeOP)
anodeOP = makevec(anodeOP)

#Filter Inf's
cathodeOPmask = isnotinf.(cathodeOP)
anodeOPmask = isnotinf.(anodeOP)
cathodeOPmask = findall(cathodeOPmask)
anodeOPmask = findall(anodeOPmask)
cathodeOPmask = [cathodeOPmask[i][1] for i in 1:length(cathodeOPmask)]
anodeOPmask = [anodeOPmask[i][1] for i in 1:length(anodeOPmask)]
cathodeX = x[cathodeOPmask,:]
anodeX = x[anodeOPmask,:]
cathodeOPFiltered = filter(!isinf,cathodeOP)
anodeOPFiltered = filter(!isinf,anodeOP)

#Filter Nans
cathodeOPmask = isnotnan.(cathodeOPFiltered)
anodeOPmask = isnotnan.(anodeOPFiltered)
cathodeOPmask = findall(cathodeOPmask)
anodeOPmask = findall(anodeOPmask)
cathodeOPmask = [cathodeOPmask[i][1] for i in 1:length(cathodeOPmask)]
anodeOPmask = [anodeOPmask[i][1] for i in 1:length(anodeOPmask)]
cathodeX = x[cathodeOPmask,:]
anodeX = x[anodeOPmask,:]
cathodeOPFiltered = filter(!isnan,cathodeOPFiltered)
anodeOPFiltered = filter(!isnan,anodeOPFiltered)

#Fit Overpotentials
mod_POP = QuadraticBMS.fitQuadraticNSD(cathodeX,cathodeOPFiltered)
mod_NOP = QuadraticBMS.fitQuadraticNSD(anodeX,anodeOPFiltered)
mod_OOP = QuadraticBMS.fitQuadraticNSD(x,OOP)

optimize!(mod_POP)
optimize!(mod_NOP)
optimize!(mod_OOP)

#Extract value from overpotentials
Q_POP = get_val(:Q,mod_POP)
q_POP = get_val(:q,mod_POP)
Q_NOP = get_val(:Q,mod_NOP)
q_NOP = get_val(:q,mod_NOP)
Q_OOP = get_val(:Q,mod_OOP)
q_OOP = get_val(:q,mod_OOP)


#Show Fits
fit_cathode_ocv = similar(cathodeOCV)
fit_total_ocv = similar(cathodeOCV)
fit_anode_op = similar(cathodeOCV)
fit_cathode_op = similar(cathodeOCV)
for i=1:length(fit_cathode_ocv)
    fit_cathode_ocv[i] = x[i,:]'*Q_COCV*x[i,:]+q_COCV'*x[i,:]
    fit_total_ocv[i] = x[i,:]'*Q_OCV*x[i,:]+q_OCV'*x[i,:]
    fit_anode_op[i] = x[i,:]'*Q_NOP*x[i,:]+q_NOP'*x[i,:]
    fit_cathode_op[i] = x[i,:]'*Q_POP*x[i,:]+q_POP'*x[i,:]
end

figure(1)
hold_on()
plot(x[:,2],fit_cathode_ocv,"r+",options=plotOptions)
plot(x[:,2],fit_total_ocv,"k+",options=plotOptions)


figure(2)
subplot(2,1,1)
hold_on()
plot3(x[:,3],x[:,2],fit_cathode_op,"r+")
subplot(2,1,2)
hold_on()
plot3(x[:,3],x[:,2],fit_anode_op,"r+")

#Set Initial conditions and build the controller
ic = [0.05,0,300]
capacity = 10800
dt = 1
N = 1000
h = cellDict["Cell cooling surface area [m2]"]*cellDict["Total heat transfer coefficient [W.m-2.K-1]"]
m = 0.5
τ_ohm = 1
specific_heat = 500

#ACCOUNTING FOR PARAMETERS
neg_elec_thickness = cellDict["Negative electrode thickness [m]"]
neg_cc_thickness = cellDict["Negative current collector thickness [m]"]
pos_elec_thickness = cellDict["Positive electrode thickness [m]"]
pos_cc_thickness = cellDict["Positive current collector thickness [m]"]
sep_thickness = cellDict["Separator thickness [m]"]

cell_volume = cellDict["Cell volume [m3]"]

neg_elec_density = cellDict["Negative electrode density [kg.m-3]"]
neg_cc_density = cellDict["Negative current collector density [kg.m-3]"]
pos_elec_density = cellDict["Positive electrode density [kg.m-3]"]
pos_cc_density = cellDict["Positive current collector density [kg.m-3]"]
sep_density = cellDict["Separator density [kg.m-3]"]

total_cell_thickness = neg_elec_thickness+pos_elec_thickness+pos_cc_thickness+neg_elec_thickness+sep_thickness

frac_vol_neg_elec = neg_elec_thickness/total_cell_thickness
frac_vol_pos_elec = pos_elec_thickness/total_cell_thickness
frac_vol_neg_cc = neg_cc_thickness/total_cell_thickness
frac_vol_pos_cc = pos_cc_thickness/total_cell_thickness
frac_vol_sep = sep_thickness/total_cell_thickness

mass_neg_elec = frac_vol_neg_elec*cell_volume*neg_elec_density
mass_neg_cc = frac_vol_neg_cc*cell_volume*neg_cc_density
mass_pos_elec = frac_vol_pos_elec*cell_volume*pos_elec_density
mass_pos_cc = frac_vol_pos_cc*cell_volume*pos_cc_density
mass_sep = frac_vol_sep*cell_volume*sep_density

cp_neg_cc = cellDict["Negative current collector specific heat capacity [J.kg-1.K-1]"]
cp_neg_elec = cellDict["Negative electrode specific heat capacity [J.kg-1.K-1]"]
cp_pos_cc = cellDict["Positive current collector specific heat capacity [J.kg-1.K-1]"]
cp_pos_elec = cellDict["Positive electrode specific heat capacity [J.kg-1.K-1]"]
cp_sep = cellDict["Separator specific heat capacity [J.kg-1.K-1]"]

cap_cell = specific_heat = (cp_neg_cc*mass_neg_cc)+
(cp_neg_elec*mass_neg_elec)+
(cp_pos_elec*mass_pos_elec)+
(cp_pos_cc*mass_pos_cc)+
(cp_sep*mass_sep)



#BUILD CONTROLLER
controller = QuadraticBMS.buildController(ic,capacity,dt,N,Q_OCV,q_OCV,Q_COCV,q_COCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,cap_cell,τ_ohm;max_T=333,pl_tol=0)
controller_t = QuadraticBMS.buildController(ic,capacity,N,Q_OCV,q_OCV,Q_COCV,q_COCV,Q_NOP,q_NOP,Q_POP,q_POP,Q_OOP,q_OOP,h,cap_cell,τ_ohm;max_T=333,pl_tol=0)
#optimize!(controller)
optimize!(controller_t)

I = get_val(:I,controller_t)
DT_val = get_val(:DT,controller_t)
t_val = cumsum(DT_val)
t_val = t_val.-t_val[1]

# Create interpolation function
function current_function(t,t_val,I)
    println(t)
    t_num = t.value*11346.612775644178
    ind = findfirst(t_val->t_val>=t_num,t_val)
    return I[ind]
end

current_function_specific(t) = current_function(t,t_val,I)

timescale=param.evaluate(model.timescale)
current_interpolant = pybamm.Interpolant(t_val, I_val, timescale * pybamm.t)

current_function_py = PyObject(current_function_specific)
param.update(PyDict(Dict("Current function [A]"=>current_interpolant)))

model = pybamm.lithium_ion.DFN()
simulation = pybamm.Simulation(model, parameter_values=param)
t_max = t_val[end]

simulation.solve(t_val,initial_soc=0)
simulation.plot()

