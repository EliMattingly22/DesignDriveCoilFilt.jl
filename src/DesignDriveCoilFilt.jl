module DesignDriveCoilFilt
using FFTW,Gtk,PyPlot,Optim
using LinearAlgebra
using PyPlot

using Interpolations

include("Load_LTSpice_Net.jl")
include("RunAnalysis.jl")
include("SPICE2Matrix.jl")
include("ToroidOptimizer.jl")
include("Filter_Designer.jl")
include("ThermalModeling.jl")
export DesignDriveFilter_OptimDrift,
       DesignDriveFilter,
       ToroidOptimizer,
       PipeFlow
       
# Write your package code here.

end
