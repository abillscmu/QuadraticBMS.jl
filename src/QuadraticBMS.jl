module QuadraticBMS

using JuMP
import Ipopt
using MosekTools
using LinearAlgebra
import Dierckx

#pybamm = PyCall.pyimport("pybamm")
FARADAY = 96485.33212
R = 8.314

#Generate DFN Data
#include("dfnDataGen.jl")

#Fit QuadraticBMS
#include("fitQuadratic.jl")

#Build the quadratic controller
include("controller.jl")

end # module
