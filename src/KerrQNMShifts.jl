__precompile__()
module KerrQNMShifts

using KerrQuasinormalModes
using PyCall
using CSV
using DataFrames
using Plots
using Cubature
using ContourIntegrals
# Write your package code here.

include("KerrQuasinormalModeAdditions.jl")
include("ImportExpressionsFromFile.jl")
include("OperatorShiftType.jl")
include("OperatorSandwichType.jl")


export qnm, qnmfunctionnew, importqnm

end
