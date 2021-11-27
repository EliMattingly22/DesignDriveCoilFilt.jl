# Overview
This tool is to facilitate designing a drive coil filter for MPI. It is coded with [Julia](https://julialang.org/), an open-source language.

It works by taking in the coil parameters (L, C, R) and then picking values for the matching network and filter. The tool then designs the ideal inductors (to include the ESR), simulates the full circuit, and outputs the values and transfer function.
## Including the package:
 In the Julia REPL run:
 
 1) Pkg> add https://github.com/EliMattingly22/DesignDriveCoilFilt.jl.git
 
 (to get from the "normal" command REPL, which shows "julia>", press "]"
 
 2) julia> using DesignDriveCoilFilt

This will install and precompile the package. These steps may take some time on the first run, due to all of the dependencies, but will be faster for future runs. 

# Main functions:
 

## DesignDriveFilter
[See the wiki for more details on function use](https://github.com/EliMattingly22/DesignDriveCoilFilt.jl/wiki/Use-of-tool)

The main function is: *DesignDriveFilter* and the required inputs are:

* LDrive, The drive coil inductance in Henries
    
* RDrive, The drive coil resistance in Ohms
    
* TargetZ, The target impedance you want to match to (I've only tested with real values)
    
* DriveFreq; The target drive frequency. Optionally this can be determined (assuming you have a value for CDrive that is reasonable) by setting "DetermineFreq" to be true 
* 
**keyword arguments (kwargs) are**

* CDrive = 1e6, The drive coil capacitance (Farads) it is very large by default to effectively bypass it.

* NumDriveElements = 1, This is an option if you have multiple repreated drive coil sections (e.g. 10 sub-coils that are wired in series, but you enter the value of a single section in LDrive, you can set this to 10).

* WireDiam = 2e-3, This is for the inductor designer. It makes sure the ESR of filter inductors are reasonable

* WireFillFac = 0.75, This is also for the inductor designer, it is a fill factor for the wires to account for using Litz. It simply scales the ESR

* PlotSrcR = TargetZ, This is not used. 

* PlotFTs = true, Optionally plotting

* VSrc = 2, This scales the input for simulation

* DetermineFreq = false, If you want the solver to determine the most appropriate drive frequency set this to true. Must include a CDrive then

* AddNotchFreq = nothing, If you want a notch frequency e.g. 2f0, or 3f0, set this to the frequency you want to notch in Hz. If this is a vector (length must be 1 or 2) it will make a notch at both values.

* FilterZ = TargetZ, The characteristic impednace of the filter

* RDampVal = FilterZ, The value for the damping resistor. 

* PerturbTxReactance = nothing, if you want to perturb the matching section you can set this to some value (Ohms). This allows for fine-tuning the frequency response and therefore optimize the filter for stability. 


## Optimized Drive Filter

By using *DesignDriveFilter_OptimDrift* (same inputs, with some extra kwargs for optimization ) the main function will be repeatedly called and the matching impedance is varied (changes the location of the "trough" *) in the transfer function (Vin->IDrive) such that there will be minimum drift given perturbation to the capactors in series with the drive coil (or distance to local minima "trough"). This optimization was picked, because in practice, these capacitors have the most current through them, and are most challenging to cool.

*The impedance matching only changes the location of drive relative to trough in the case that a notch is included. If no notch is included, then you can choose to use the PerturbTxReactance option. 
 
 The extra kwargs are: 
 BruteForceOpt = false,
 MinimizeDistToTrough = true
 
 BruteForceOpt is optimizing simply by running a large list of perturbations to the reactance and then picking the one with minimum cost (dist to trough or drift coeff). 
 
 MinimizeDistToTrough is the boolean option to minimize the distance ofthe drive freqeuncy to the trough. If false, it will minimize the drift coeff. in the drive capacitor. 
 

## Optimized air-core toroid
Also included are other useful functions such as the toroid optimizer, which creates a toroid core shape that has minimal ESR (See "D-Shaped toroidal cage inductors" P.N. Murgatroyd & D. Belahrache,1989)

The function is: *ToroidOptimizer* and the required inputs are:

* Dia, wire diameter in meters
* LTarget, target inductance in Henries
    
    
* NumLayers = 2, the number of layers on the toroid. Having this be above 2 may cause issues for heat removal, or parasitic effects

* CoreMu = 1, In case you want to use a ferrite core. Typically it should be 1 (air)

* Alpha = 2, Ratio of OD/ID. In theory, higher numbers are more efficient, but have more leakage. Higher numbers are also more challenging to wind. Should be between 2-5

* CuFillFactor = 1, The packing factor of the wire. Solid magnet wire is 1, Litz will be <1

### DCore_DetermineIdealInduct

This function determines the (ideal) inductance of a D shaped toroid given geometry and turn count. 

## Rogowski_Calc

This function calculates the voltage per amp in a Rogowski coil being used as a current sensor. There are methods for rectangular coils, and circular cross-sections.



