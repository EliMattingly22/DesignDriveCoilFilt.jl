

# include("Load_LTSpice_Net.jl")
# include("RunAnalysis.jl")
# include("SPICE2Matrix.jl")
# include("ToroidOptimizer.jl")


"""
This function takes in similar inputs to DesignDriveFilter. It will iterate that function with varying matching impedances to minimize the expected drift in the current. 
"""
function DesignDriveFilter_OptimDrift(
    LDrive,
    RDrive,
    TargetZ,
    DriveFreq;
    CDrive = 1e6,
    NumDriveElements = 1,
    WireDiam = 2e-3,
    WireFillFac = 0.75,
    PlotSrcR = TargetZ,
    PlotFTs = true,
    VSrc = 2,
    DetermineFreq = false,
    AddNotchFreq = nothing,
    FilterZ = TargetZ,
    RDampVal = FilterZ,
    PerturbTxReactance = nothing
)

    
        function WrapperFunkZin(Zin)
            DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList = DesignDriveFilter(
            LDrive,RDrive, Zin,DriveFreq;
            CDrive = CDrive,
            NumDriveElements = NumDriveElements,
            WireDiam = WireDiam,
            DetermineFreq=DetermineFreq,
            AddNotchFreq = AddNotchFreq,
            FilterZ = Zin)
        
            
            DetermineComponentsTempCoeffs(SPICE_DF,InputList,1,DriveFreq,"LDrive")
            if CDrive>1 #Check the drift in the drive capacitor if it is used, otherwise use the CSer 
                Drift = SPICE_DF.DriftCoeff[findfirst(isequal("CDrive"),SPICE_DF.Name)]
            else
                Drift = SPICE_DF.DriftCoeff[findfirst(isequal("CSer"),SPICE_DF.Name)]
            end

            return Drift
        end
    
        
        function WrapperFunk_Per(PerturbationX)
            DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList = DesignDriveFilter(
            LDrive,RDrive, TargetZ,DriveFreq;
            CDrive = CDrive,
            NumDriveElements = NumDriveElements,
            WireDiam = WireDiam,
            DetermineFreq=DetermineFreq,
            AddNotchFreq = AddNotchFreq,
            FilterZ = TargetZ,
            PerturbTxReactance = PerturbationX)
        
            
            DetermineComponentsTempCoeffs(SPICE_DF,InputList,1,DriveFreq,"LDrive")
            if CDrive>1 #Check the drift in the drive capacitor if it is used, otherwise use the CSer 
                Drift = SPICE_DF.DriftCoeff[findfirst(isequal("CDrive"),SPICE_DF.Name)]
            else
                Drift = SPICE_DF.DriftCoeff[findfirst(isequal("CSer"),SPICE_DF.Name)]
            end

            return Drift
        end
    
    

    if  (PerturbTxReactance === nothing)
        MinZ = optimize(WrapperFunkZin,5,25)
        DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList = DesignDriveFilter(
        LDrive,RDrive, MinZ.minimizer,DriveFreq;
        CDrive = CDrive,
        WireDiam = WireDiam,
        DetermineFreq=DetermineFreq,
        AddNotchFreq = AddNotchFreq,
        FilterZ = MinZ.minimizer)
        MinZVal = MinZ.minimizer
        DetermineComponentsTempCoeffs(SPICE_DF,InputList,1,DriveFreq,"LDrive")
        return MinZVal, DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList
    else
        println("Perturbation given")
        PerVec = -2:.1:2
        DriftVec = zeros(length(PerVec),1)
        for i in 1:length(PerVec)
            println(i)
            DriftVec[i] =WrapperFunk_Per(PerVec[i])
        end
        MinIndex = findfirst(x-> x==minimum(DriftVec),DriftVec)
        MinXVal = PerVec[MinIndex]

        DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList = DesignDriveFilter(
        LDrive,RDrive, TargetZ,DriveFreq;
        CDrive = CDrive,
        WireDiam = WireDiam,
        DetermineFreq=DetermineFreq,
        AddNotchFreq = AddNotchFreq,
        FilterZ = TargetZ,
        PerturbTxReactance = MinXVal)
        
        DetermineComponentsTempCoeffs(SPICE_DF,InputList,1,DriveFreq,"LDrive")
        return MinXVal, DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList
    end

end



"""
This function takes in  the following parameters:

    LDrive, drive coil inductance (H)

    RDrive, Drive coil resistance (Ω)

    TargetZ, Amplifier/target impedance to match to  (Ω)

    DriveFreq; drive frequency (Hz)
The following are keyword arguments and should be entered after the first four (required) as;
DesignDriveFilter(LDrive,RDrive,Amp_Z,DriveFreq; [kwarg] = value). See below for an example

    CDrive = 1e6, Any series capacitance with drive coil

    NumDriveElements = 1, The number of series copies of the L,C,R for example if there are 6 pancakes of LDrive, RDrive, Cdrive, then this would be 6

    WireDiam = 2e-3, Wire diameter for toroids

    WireFillFac = 0.75, If litz wire, this is the ratio of copper to air

    PlotFTs = true, Plotting be default is true
    VSrc = 2, The source voltage for plotting
    DetermineFreq = false, IF you have a fixed capacitor (CDrive has been set) you may want to search for an ideal drive freq. in which case, set this to true.
    AddNotchFreq = nothing, This can be a single value of a notch frequency, or an array of length 2, the two values are the notch frequencies

    FilterZ = TargetZ, The characteristic impedance of the filtering section
    RDampVal = 100, [Ω] The damping resistor added in around the series L-C sections
    PerturbTxReactance = nothing, an optional element to add some purturbation to the reactance of the drive coil, which effectively enables the user to tune where the operating point ends up on the transfer function (V in → I drive)


    Example (Copy/paste):
    RDrive = 400e-3
    LDrive = 400e-6
    Amp_Z = 4
    DriveFreq = 25e3
    AddNotchFreqList = [DriveFreq*2, DriveFreq*3]

    DesignDriveFilter(LDrive,RDrive,Amp_Z,DriveFreq; AddNotchFreq = AddNotchFreqList)




"""
function DesignDriveFilter(
    LDrive,
    RDrive,
    TargetZ,
    DriveFreq;
    CDrive = 1e6,
    NumDriveElements = 1,
    WireDiam = 2e-3,
    WireFillFac = 0.75,
    PlotSrcR = TargetZ,
    PlotFTs = true,
    VSrc = 2,
    DetermineFreq = false,
    AddNotchFreq = nothing,
    FilterZ = TargetZ,
    RDampVal = 100,
    PerturbTxReactance = nothing
)


    ZeroVal = 1e-9 #There are issues with zero-valued components
    ωDr = 2 * π * DriveFreq
    ZDrive =
        RDrive * NumDriveElements +
        im * ωDr * LDrive * NumDriveElements +
        NumDriveElements ./ (im * 2 * pi * DriveFreq * CDrive) #Equivalent complex impednace of the load

    Reactance_Load =
        ωDr * LDrive * NumDriveElements - NumDriveElements ./ (ωDr * CDrive)
    Reactance_Load = round(Reactance_Load; sigdigits = 3)
    println("The reactance of the load is: $(round(Reactance_Load;sigdigits=3)) Ω")
    if ~DetermineFreq
        if (Reactance_Load > 0)

            SerCap,CParAct = ImpMatch_LLoad(TargetZ, ZDrive, Reactance_Load,ωDr)

            LTee_2 = ZeroVal
            LTee_2_ESR = ZeroVal
            LTee_1 = ZeroVal
            LTee_1_ESR = ZeroVal


        elseif (Reactance_Load < 0)

            LTee_2,LTee_2_ESR,CParAct = ImpMatch_CLoad(TargetZ, ZDrive, Reactance_Load,ωDr)

            LTee_1 = ZeroVal
            LTee_1_ESR = ZeroVal
            SerCap = 1e6

        elseif (Reactance_Load == 0)
            matchRatio = real(TargetZ) / real(ZDrive)
            Q = sqrt(matchRatio - 1)
            Xs = Q * real(ZDrive)
            LSer2 = Xs / (ωDr)
            LTee_2 = LSer2
            LTee_1 = ZeroVal
            LTee_1_ESR = ZeroVal
            CParAct = findResPair((1 + Q^(-2)) * LSer2, DriveFreq)
            LTee_2_Geom =
                ToroidOptimizer(WireDiam, LTee_2; CuFillFactor = WireFillFac)
            LTee_2_ESR = LTee_2_Geom.DCore.Resistance
            SerCap = 1e6
        end
    else

        matchRatio = real(TargetZ) / real(ZDrive)
        Q = √(matchRatio - 1)
        Xs = Q * real(ZDrive)

        DriveFreq = findFilterFreq(Q,real(ZDrive), LDrive * NumDriveElements, CDrive / NumDriveElements)
        ωDr = 2*π*DriveFreq
        println("Determined Q to be $(round(Q))")
        println("Determined drive freq. to be $(round(DriveFreq)) Hz")
        ZSer = 1*im * ωDr*LDrive * NumDriveElements -  1*im /(ωDr* CDrive / NumDriveElements) +RDrive* NumDriveElements
        YSer = 1/ZSer

        CParAct = abs.(imag(YSer))/ωDr
        LTee_2 = ZeroVal
        LTee_2_ESR = ZeroVal
        LTee_1 = ZeroVal
        LTee_1_ESR = ZeroVal
        SerCap = 1e6
    end

    if ~(PerturbTxReactance === nothing)
        PerturbInductance = PerturbTxReactance/ωDr
        LTee_2 = LTee_2+PerturbInductance
        ZDrive = RDrive * NumDriveElements +
        im * ωDr * LDrive * NumDriveElements +
        NumDriveElements ./ (im * 2 * pi * DriveFreq * CDrive)
        Z_SerMatchingSect = LTee_2_ESR+im * ωDr * LTee_2+ 1/(im*ωDr*SerCap)
        ZSerTot = Z_SerMatchingSect+ZDrive
        Z_TotMatchSect = Par([Z_Cap(CParAct,DriveFreq),ZSerTot])
        LTee_1 = abs.(-1*imag(Z_TotMatchSect)/ωDr)

    end
    println("L Tee 1 =  $(round(LTee_1*1e6;sigdigits=3))μH ")
    println("L Tee 2 =  $(round(LTee_2*1e6;sigdigits=3))μH ")
    println("CParAct =  $(round(CParAct*1e6;sigdigits=3))μF ")
    println("SerCap =  $(round(SerCap*1e6;sigdigits=3))μF ")


    (LFilt, CFilt) = Butterworth_2(FilterZ, DriveFreq)
    LFiltMatch_C = findResPair(LFilt, DriveFreq)
    CFiltMatch_L = findResPair(CFilt, DriveFreq)

    LFilt1_Geom =
            ToroidOptimizer(WireDiam, LFilt; CuFillFactor = WireFillFac)
    LFilt1_ESR = LFilt1_Geom.DCore.Resistance

    LFilt2_Geom =
            ToroidOptimizer(WireDiam, CFiltMatch_L; CuFillFactor = WireFillFac)
    LFilt2_ESR = LFilt2_Geom.DCore.Resistance
    println(
        "LFilt =  $(round(LFilt*1e6;sigdigits=3))μH matched with: LFiltMatch_C =  $(round(LFiltMatch_C*1e6;sigdigits=3))μF ",
    )
    println(
        "CFilt =  $(round(CFilt*1e6;sigdigits=3))μF matched with: CFiltMatch_L = $(round(CFiltMatch_L*1e6;sigdigits=3))μH ",
    )

    

    if ~(AddNotchFreq===nothing)
    
        LNotch, CNotch, LNotch_Tune = makeNotchSection(AddNotchFreq[1], DriveFreq, TargetZ;LVal = 100e-6)
        LNotch_Geom =
                ToroidOptimizer(WireDiam, LNotch; CuFillFactor = WireFillFac)
        LNotch_ESR = LNotch_Geom.DCore.Resistance
        LNotch_Tune_Geom =
                ToroidOptimizer(WireDiam, LNotch_Tune; CuFillFactor = WireFillFac)
        LNotch_Tune_ESR = LNotch_Tune_Geom.DCore.Resistance
        println("Adding notch")
        if length(AddNotchFreq)==2
            LNotch2, CNotch2, LNotch_Tune2 = makeNotchSection(AddNotchFreq[2], DriveFreq, TargetZ;LVal = 100e-6)
            LNotch2_Geom =
                    ToroidOptimizer(WireDiam, LNotch2; CuFillFactor = WireFillFac)
            LNotch2_ESR = LNotch2_Geom.DCore.Resistance
            LNotch_Tune = LNotch_Tune+LNotch_Tune2
            LNotch_Tune_Geom =
                    ToroidOptimizer(WireDiam, LNotch_Tune; CuFillFactor = WireFillFac)
            LNotch_Tune_ESR = LNotch_Tune_Geom.DCore.Resistance
            println("Adding second notch")
        else
            LNotch2 = 10
            CNotch2 = 1e-12
            LNotch2_ESR = 1e-3
           
            
        end
            
    else
        LNotch = 10
        CNotch = 1e-12
        LNotch_Tune = 1e-12
        LNotch_ESR = 1e-3
        LNotch_Tune_ESR = 1e-3
        LNotch2 = 10 #Set default values is no notch is added
        CNotch2 = 1e-12 
        LNotch2_ESR = 1e-3
    end

    CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList = CircModel_SPICE(DriveFreq, VSrc,
    1e-3,
    RDrive,
    NumDriveElements,
    LDrive,
    CDrive,
    SerCap,
    LTee_2,
    LTee_2_ESR,
    LTee_1,
    LTee_1_ESR,
    CParAct,
    LFilt,
    LFilt2_ESR,
    LFilt1_ESR,
    LFiltMatch_C,
    CFilt,
    CFiltMatch_L,
    LNotch,
    LNotch_ESR,
    CNotch,
    LNotch2,
    LNotch2_ESR,
    CNotch2,
    LNotch_Tune,
    LNotch_Tune_ESR;
    PlotOn = PlotFTs,
    RDampVal = RDampVal)

    return DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList

end



function CircModel_SPICE(DriveFreq, VSrc,
    PlotSrcR,
    RDrive,
    NumDriveElements,
    LDrive,
    CDrive,
    SerCap,
    LTee_2,
    LTee_2_ESR,
    LTee_1,
    LTee_1_ESR,
    CParAct,
    LFilt,
    LFilt2_ESR,
    LFilt1_ESR,
    LFiltMatch_C,
    CFilt,
    CFiltMatch_L,
    LNotch,
    LNotch_ESR,
    CNotch,
    LNotch2,
    LNotch2_ESR,
    CNotch2,
    LNotch_Tune,
    LNotch_Tune_ESR;
    PlotOn = true,
    FreqList = 1000:10:100e3,
    ArchetypeNetFileName = nothing,
    RDampVal = nothing)




    if ArchetypeNetFileName === nothing

        ArchetypeNetFileName = joinpath(dirname(pathof(DesignDriveCoilFilt)),"Filter_Archetype_Damped_2.net")

    end
    SPICE_DF,NodeList,InputList,NumVSources = SPICE2Matrix(ArchetypeNetFileName)
    UpdateElementVal!(SPICE_DF,"RSrc",PlotSrcR)

    UpdateElementVal!(SPICE_DF,"LNotch_Tune",LNotch_Tune)
    UpdateElementVal!(SPICE_DF,"LNotch",LNotch)
    UpdateElementESR!(SPICE_DF,"LNotch_Tune",LNotch_Tune_ESR)
    UpdateElementESR!(SPICE_DF,"LNotch",LNotch_ESR)
    UpdateElementVal!(SPICE_DF,"CNotch",CNotch)
    UpdateElementVal!(SPICE_DF,"LNotch2",LNotch2)
    UpdateElementESR!(SPICE_DF,"LNotch2",LNotch2_ESR)
    UpdateElementVal!(SPICE_DF,"CNotch2",CNotch2)

    UpdateElementVal!(SPICE_DF,"LFilt1",LFilt)
    UpdateElementVal!(SPICE_DF,"LFilt2",LFilt)
    UpdateElementVal!(SPICE_DF,"CLFilt_C1",LFiltMatch_C)
    UpdateElementVal!(SPICE_DF,"CLFilt_C2",LFiltMatch_C)
    UpdateElementESR!(SPICE_DF,"LFilt1",LFilt1_ESR)

    UpdateElementVal!(SPICE_DF,"CFilt",CFilt)
    UpdateElementVal!(SPICE_DF,"LCFilt_L",CFiltMatch_L)
    UpdateElementESR!(SPICE_DF,"LCFilt_L",LFilt2_ESR)

    UpdateElementVal!(SPICE_DF,"RDrive",RDrive)
    UpdateElementVal!(SPICE_DF,"LDrive",NumDriveElements*LDrive)
    UpdateElementVal!(SPICE_DF,"CDrive",CDrive/NumDriveElements)

    UpdateElementVal!(SPICE_DF,"CSer",SerCap)
    UpdateElementVal!(SPICE_DF,"CPar",CParAct)

    UpdateElementVal!(SPICE_DF,"LTee1",LTee_1)
    UpdateElementESR!(SPICE_DF,"LTee1",LTee_1_ESR)
    UpdateElementVal!(SPICE_DF,"LTee2",LTee_2)
    UpdateElementESR!(SPICE_DF,"LTee2",LTee_2_ESR)


    if ~(RDampVal===nothing)

        UpdateElementVal!(SPICE_DF,"LDamp",LFilt)
        UpdateElementESR!(SPICE_DF,"LDamp",LFilt1_ESR)
        UpdateElementVal!(SPICE_DF,"CDamp",LFiltMatch_C)
        UpdateElementVal!(SPICE_DF,"RDamp2",RDampVal)
        UpdateElementVal!(SPICE_DF,"RDamp1",RDampVal)

    end


    inputs = zeros(length(InputList))
    inputs[end-(NumVSources-1):end] .= 1

    ResultNodeNames = vcat("V(".*NodeList.*")", "I(".*(InputList[end-(NumVSources-1):end]).*")")


    Results = [ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results = hcat(Results...)
    ResDict = Dict(ResultNodeNames[1] => Results[1,:])
    for i in 2:length(ResultNodeNames)

         merge!(ResDict,Dict(ResultNodeNames[i] => Results[i,:]))

    end

    CurrentVec = VSrc.* plotACElCurrent(SPICE_DF,FreqList,Results,"LDrive")
    
    if PlotOn
        pygui(true)
        semilogy(FreqList,abs.(CurrentVec[:]))
        
        xlabel("Frequency, Hz")
        ylabel("Current in LDrive")
    end

    getElementCurrents(SPICE_DF,Results,DriveFreq)
    return CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList
end


function ImpMatch_LLoad(TargetZ, ZDrive, Reactance_Load,ωDr)
    println("Load is inductive")
    EquivSerL = Reactance_Load / ωDr
    println(
        "Load appears to be a inductor with a value of $(round(EquivSerL*1e6;sigdigits=3))μH and a ESR of $(real(ZDrive))Ω ",
    )

    matchRatio = real(TargetZ) / real(ZDrive)
    Q = sqrt(matchRatio - 1)
    Xs = Q * real(ZDrive) #Target reactance

    X_SerCap = Reactance_Load - Xs #Reactance of series capacitor
    SerCap = 1 / (ωDr * X_SerCap) #Series Capacitor value
    CParAct = Q / (ωDr * TargetZ)


    return SerCap,CParAct

end

function ImpMatch_CLoad(TargetZ, ZDrive, Reactance_Load,ωDr; WireDiam = 2e-3,
    WireFillFac = 0.75)
    println("Load is capacitive")
        EquivSerC = abs.(1 / (Reactance_Load * ωDr))
        println(
            "Load appears to be a capacitor with a value of $(round(EquivSerC*1e6;sigdigits=3))μF and a ESR of $(real(ZDrive))Ω ",
        )
        LSer = abs.(Reactance_Load) / (ωDr)
        matchRatio = real(TargetZ) / real(ZDrive)
        Q = sqrt(matchRatio - 1)
        Xs = Q * real(ZDrive)
        LSer2 = Xs / (ωDr)
        LTee_2 = LSer2 + LSer

        CParAct = findResPair((1 + Q^(-2)) * LSer2, ωDr*2*π)
        LTee_2_Geom =
            ToroidOptimizer(WireDiam, LTee_2; CuFillFactor = WireFillFac)
        LTee_2_ESR = LTee_2_Geom.DCore.Resistance
        println("Q = $Q")

        return LTee_2,LTee_2_ESR,CParAct

end

"""
This function determines the ideal frequency of a resonant RLC such that it can be matched to a different impedance. It assumes there is already an L, C, and R, so at some frequency it the real part of the admittance will be the matched admittance. This admittance is made fully real with a shunt reactive element

"""
function findFilterFreq(Q,RLoad,L,C)
        f = ((Q*RLoad + √(Q^2*RLoad^2+4*L/C)) / (4*π*L) , (Q*RLoad - √(Q^2*RLoad^2+4*L/C)) / (4*π*L) )

        return maximum(f)
end


"""
This function makes a filter section that notches out a frequency (NotchFreq), but is apparently real at another (DriveFreq)
    the inputs are:
        * NotchFreq: the freq that will be shunted out
        * DriveFreq: The freq that should remain real
        * RMatch: The resistance the drive coil is tuned to
        * LVal is a kwarg because either the inductor or capacitor value should be pre-defined
The outputs are:
    LVal(Henries)
    CVal(Farads)
    LTune(Henries)
"""
function makeNotchSection(NotchFreq, DriveFreq, RMatch;LVal = 100e-6)
    CVal = findResPair(LVal, NotchFreq)
    Resid_Reactance = -1*(imag(Par(Z_Cap(CVal,DriveFreq)+Z_Ind(LVal,DriveFreq),RMatch)))
    LTune = Resid_Reactance/(2*π*DriveFreq)
    return LVal,CVal,LTune
end

"""

Calculates parallel impedances.
If given two inputs, it assumes both are inputs
If given an array, it takes the parallel combo of it all

Eg.
ZList = [1 + 1*im, 1]
Par(ZList[1], ZList[2])
Par(ZList)
"""
function Par(Z1, Z2)
    return 1 / (1 / Z1 + 1 / Z2)
end
function Par(Z::Array)
    Y = 1 ./Z
    YTotal = sum(Y)
    Zeff = 1/YTotal
    return Zeff
end

"""
Parallel component lists
Input is Array{Tuple{Number,String},1}
e.g.
Freq = 25e3
CompList = [(3, "R"),(30e-6, "L"),(200e-9,"C")]
Par(CompList,Freq)

"""
function Par(C::Array{Tuple{Real,String},1},freq)

    Z = Complex.(zeros(length(C)))
    for i in 1:length(C)
        if (lowercase(C[i][2])=="r") | (lowercase(C[i][2])=="resistor")
            Z[i] = C[i][1]
        elseif (lowercase(C[i][2])=="l") | (lowercase(C[i][2])=="inductor" )
            Z[i] = Z_Ind(C[i][1],freq)
        elseif (lowercase(C[i][2])=="c") | (lowercase(C[i][2])=="capacitor" )
            Z[i] = Z_Cap(C[i][1],freq)
        elseif lowercase(C[i][2])=="z"
            Z[i] = C[i][1]
        else
            error("Unknown component")
        end

    end


    Y = 1 ./Z
    YTotal = sum(Y)
    Zeff = 1/YTotal
    return Zeff
end

"""
Converts from a polar coordinate (Mag, phase) to cartesian
if you enter a tuple or array it assumes first is mag.

"""
function Polar2Cart(Mag, Phase)
    Real = real.(Mag * exp(-im * Phase / 360 * 2 * pi))
    Imag = imag.(Mag * exp(-im * Phase / 360 * 2 * pi))
    return (Real, Imag)
end

function Polar2Cart(MagPhaseTup::Union{Tuple, Array})
    Polar2Cart(MagPhaseTup[1], MagPhaseTup[2])
end

function Z_Cap(C, f)
    return (1 ./ (im * 2 * pi .* f .* C))
end
function Z_Ind(L, f)
    return im * 2 * pi .* f .* L
end


"""
This function takes in some reactance element (E) and a frequency, and finds the pair


"""
function findResPair(E, f)
    ##Finds the resonant pair for a reactive element, E
    ## ωL = 1/(ωC) → L = 1/(ω^2C),and vice versa
    return 1 ./ ((2 * pi * f) .^ 2 * E)
end

function Butterworth_2(Z, f)
    #Sqrt(L/C) = Z ; f = 1/(2π*sqrt(LC))
    #L = Z^2*C
    #f = 1/(2π*Z C) → C = 1/(2*pi*f*Z)
    #Butterworth 2nd order = [LX = sqrt(2)Z, CX = sqrt(2)Z]
    C = 1 / (2 * pi * f) / Z * sqrt(2)
    L = 1 / (2 * pi * f) * Z * sqrt(2)
    return (L, C)
end




function PlotPhasor(
    Impedance::Tuple,
    PrevImpedance;
    Color = nothing,
    NumSteps = 1,
    Admittance = false,
    HeadWidth = 0.01,
)

    ## This function takes in an impedance, and plots it using 'phasor()'
    #optionally, it takes in a previous impedance which can be added in series or parallel
    #It then plots the line using PrevImpedance as a base
    #If Admittance is true, it means the inputs are Admittance not impedance, so parallel and series are switched
    NewImp = 0 + 0im
    dx = 0
    dy = 0
    for I = 1:NumSteps
        if I > 1
            #    Color = "gray"
        end
        if lowercase(Impedance[2]) == "series"
            if !Admittance
                NewImp = PrevImpedance + Impedance[1] / NumSteps
            else
                NewImp = Par(PrevImpedance, Impedance[1] * NumSteps)
            end
        elseif lowercase(Impedance[2]) == "parallel"
            if !Admittance
                NewImp = Par(PrevImpedance, Impedance[1] * NumSteps)
            else
                NewImp = PrevImpedance + Impedance[1] / NumSteps
            end
        else
            error("must be series or parallel")
        end

        RealParts = [real(NewImp), real(PrevImpedance)]
        ImagParts = [imag(NewImp), imag(PrevImpedance)]

        dx = -(RealParts[2] - RealParts[1])
        dy = -(ImagParts[2] - ImagParts[1])
        if Color !== nothing



            PyPlot.plot(RealParts, ImagParts, color = "gray")

            if I == 1
                phasor(PrevImpedance, color = Color)
            end
        else
            PyPlot.plot(RealParts, ImagParts)
            if I == 1
                phasor(PrevImpedance)
            end
        end
        PrevImpedance = NewImp
    end
    # arrow(real(NewImp)-dx,imag(NewImp)-dy,dx,dy,head_width=HeadWidth,fc="gray",ec="k",alpha=1.,width = 0.00001,head_length= 2*HeadWidth)

    xlabel("Real")
    ylabel("Imag.")
    return NewImp
end


"""
This function takes in an array of tuples structures such as:

    [(<Complex impedance #1> , <"series" or "Parallel">),
    (<Complex impedance #2> , <"series" or "Parallel">)
    .
    .
    .
    (<Complex impedance #n> , <"series" or "Parallel">)]

    e.g.

    ZList = [(1+1*im, "series"),
             (1+5*im, "series"),
             (5+1*im, "parallel"),]
    PlotImpedanceTransformList(ZList)
"""
function PlotImpedanceTransformList(ImpList; InitialImp = nothing, ArrHeadWidth =1)
    ColorList = ["r", "g", "b", "magenta", "teal"]
    L = zeros(length(ImpList))
    for I = 1:length(ImpList)
        L[I] = abs(ImpList[I][1])
    end
    LongestVector = maximum(L)
    HeadWidth = LongestVector / ArrHeadWidth
    ##Impedance Plotting
    if InitialImp !== nothing
        PrevImp = InitialImp
        StartVal = 1
    else
        PrevImp = ImpList[1][1]
        StartVal = 2
    end
    for I = StartVal:length(ImpList)
        PrevImp = PlotPhasor(
            ImpList[I],
            PrevImp;
            NumSteps = 500,
            Color = ColorList[I],
            HeadWidth = HeadWidth,
        )
    end
    phasor(PrevImp, color = ColorList[length(ImpList)+1])
    title("Impedance")
    ##Admittance Plotting
    figure()

    L = zeros(length(ImpList))
    for I = 1:length(ImpList)
        L[I] = abs(1 / ImpList[I][1])
    end
    LongestVector = maximum(L)
    HeadWidth = LongestVector / ArrHeadWidth
    if InitialImp !=0 nothing
        PrevImp = 1 / InitialImp
        StartVal = 1
    else
        PrevImp = 1 / ImpList[1][1]
        StartVal = 2
    end
    for I = StartVal:length(ImpList)
        PrevImp = PlotPhasor(
            (1 / ImpList[I][1], ImpList[I][2]),
            PrevImp;
            NumSteps = 500,
            Color = ColorList[I],
            Admittance = true,
            HeadWidth = HeadWidth,
        )
    end
    phasor(PrevImp, color = ColorList[length(ImpList)+1])
    title("Admittance")
end



function findMinima_Bounds(x,y;StartVal=nothing,StopVal=nothing)



    if !(StartVal===nothing)
        StartIndex = minimum(findall(x .>= StartVal))
    else
        StartIndex = 1
    end

    if !(StopVal===nothing)
        StopIndex = maximum(findall(x .<= StopVal))
    else
        StopIndex = length(x)
    end

    
    xSubset = x[StartIndex:StopIndex]
    ySubset = y[StartIndex:StopIndex]
    FirstDiff = diff(ySubset)
    SecondDiff = diff(FirstDiff)
    Tmp = SecondDiff .< 0
    push!(Tmp,1)
    FirstDiff[Tmp] .= 1000 #Only look at values where the second derivative is positive

    minIndex = findfirst(isequal(minimum(abs.(FirstDiff))),abs.(FirstDiff))

    Min_x = xSubset[minIndex]
    Min_y = ySubset[minIndex]
    return Min_x, Min_y

end