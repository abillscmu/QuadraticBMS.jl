#includes
using QuadraticBMS
using PyCall
using MATLABPlots
using JuMP

pybamm = pyimport("pybamm")
pandas = pyimport("pandas")
np = pyimport("numpy")

model = pybamm.lithium_ion.DFN()
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)
drive_cycle = pandas.read_csv("optimal_input.csv",comment="#").to_numpy()

timescale = param.evaluate(model.timescale)
current_interpolant = pybamm.Interpolant(drive_cycle[:, 1], drive_cycle[:, 2],timescale * pybamm.t)



# set up simulation
simulation = pybamm.Simulation(model, parameter_values=param)


#simulation.plot(output_variables, labels=label)