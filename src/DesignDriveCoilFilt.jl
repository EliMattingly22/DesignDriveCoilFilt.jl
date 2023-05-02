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
include("ImpedanceTransformations.jl")
include("SaveToroidSVG.jl")
export DesignDriveFilter,
       ToroidOptimizer,
       PipeFlow,
       findResPair,
       findEquivLC,
       Z_Cap,
       Z_Ind,
       RunACAnalysis,
       SPICE2Matrix,
       SPICE_DF2Matrix_ω
       
# Write your package code here.

end
