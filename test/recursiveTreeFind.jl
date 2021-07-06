using Revise
using Pkg
Pkg.activate(".")
using PyCall
pybamm = pyimport("pybamm")
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Marquis2019)
cell = Dict(py"dict($param)")
exc = cell["Positive electrode exchange-current density [A.m-2]"]
tree = exc(1000,1000,298)
a,b=tree.children
c,d = b.children
e,f = c.children
