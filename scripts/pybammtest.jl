#includes
using QuadraticBMS
using PyCall
using MATLABPlots
using JuMP

pybamm = pyimport("pybamm")
pandas = pyimport("pandas")
np = pyimport("numpy")

model = pybamm.lithium_ion.DFN()
#parameter_values = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)
parameter_values = model.default_parameter_values
f(t) = t
f2 = PyObject(f)

parameter_values.update(PyDict(Dict("Current function [A]"=>f2)))

sim = pybamm.Simulation(model, parameter_values=parameter_values)
sim.solve([0, 15],initial_soc=0.1)
sim.plot()