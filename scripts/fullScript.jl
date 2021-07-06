#includes
using QuadraticBMS
using PyCall
using MATLABPlots
using JuMP
pybamm = pyimport("pybamm")

iapp = collect(range(-10,stop=10,length=100))

#dfn data generation
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)
cellDict = py"dict($param)"
cathodeOCV,anodeOCV,cathodeOP,anodeOP = QuadraticBMS.generateDFNData(cellDict,iapp,25)
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

x = hcat(iapp,soc,ones(length(IAPP)))
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
cathodeOP = makevec(cathodeOP)
anodeOP = makevec(anodeOP)

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

mod_POP = QuadraticBMS.fitQuadraticNSD(cathodeX,cathodeOPFiltered)
mod_NOP = QuadraticBMS.fitQuadraticNSD(anodeX,anodeOPFiltered)

optimize!(mod_POP)
optimize!(mod_NOP)

Q_POP = get_val(:Q,mod_POP)
q_POP = get_val(:q,mod_POP)
Q_NOP = get_val(:Q,mod_NOP)
q_NOP = get_val(:q,mod_NOP)





