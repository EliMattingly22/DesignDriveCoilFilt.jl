module DesignDriveCoilFilt
using FFTW,PyPlot,Optim
using LinearAlgebra
using PyPlot
using CSV

# using Interpolations

include("Load_LTSpice_Net.jl")
include("RunAnalysis.jl")
include("SPICE2Matrix.jl")
include("ToroidOptimizer.jl")
include("Filter_Designer.jl")
include("ThermalModeling.jl")
include("WriteLTSPICEFile.jl")
export DesignDriveFilter_OptimDrift,
       DesignDriveFilter,
       ToroidOptimizer,
       PipeFlow,
       findResPair
       
# Write your package code here.

end
