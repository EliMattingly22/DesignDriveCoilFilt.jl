module DesignDriveCoilFilt
using FFTW,Gtk,PyPlot,Optim
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
include("Filter_Designer_GUI.jl")
export DesignDriveFilter_OptimDrift,
       DesignDriveFilter,
       ToroidOptimizer,
       PipeFlow,
       findResPair,
       useFilterDesignGUI
       
# Write your package code here.

end
