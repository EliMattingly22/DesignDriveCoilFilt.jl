
"""
This function runs an AC analysis (freq. sweep) for an LTSPICE netlist file.

    if no inputs are given, it will prompt for the user to pick a file

    Inputs:
    FileName: Path to netlist (LTSPICE: View→SPICE netlist, and it will show a text file as well as populate in PWD)
    FreqList: List of frequencies to test (Hz), Default is 100 Hz to 100kHz in steps of 10
    inputs: vector of independent nodal current inputs (1: NumNodes) and independent V inputs (NumNodes+1:End)
            All voltage inputs default to 1
"""
function RunACAnalysis(FileName=nothing, FreqList = 100:10:100e3, inputs = nothing)
    # if FileName===nothing
    #         FileName = open_dialog("Pick a file")
    # end

    SPICE_DF,NodeList,InputList,NumVSources = SPICE2Matrix(FileName)
    # println(InputList)
    if inputs===nothing
        # println("No inputs given")
        inputs = zeros(length(InputList))
        inputs[end-(NumVSources-1):end] .= 1
    end

    ResultNodeNames = vcat("V(".*NodeList.*")", "I(".*(InputList[end-(NumVSources-1):end]).*")")


    Results = [ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results = hcat(Results...)
    ResDict = Dict(ResultNodeNames[1] => Results[1,:])
    for i in 2:length(ResultNodeNames)

         merge!(ResDict,Dict(ResultNodeNames[i] => Results[i,:]))

    end




    return FreqList,Results, ResDict

end

"""
This function takes in:
    Filename to LTSPICE netlist (View→ SPICE netlist)
    DriveFreq, primary frequency of interest in Hz
    Θᵣ = 1, Θc = 1 are the lumped thermal resistance of inductors and resistors ( Θᵣ) and capacitors(Θc) in ᵒC per Watt
    TCᵣ = 0.004, TCc = 0.0003 are the temperature coefficients of the components-- TCᵣis the resistance tempco, and TCc is the capacitor tempco
        Units for TC is in UL per deg C

    InputScaling is the scaling of the input vector, for instance if you use the default input vec it will default to 1V at voltage sources, but this scalar can be used to modify that value

"""
function DetermineTempCo(FileName=nothing, DriveFreq = 25e3, ComponentName = "Ldrive";
    FreqList = 100:10:100e3,
    inputs = nothing,
    Θᵣ = 1, Θc = 1,
    TCᵣ = 0.004, TCc = -0.0003,
    InputScaling = 1)



    # if FileName===nothing
    #     FileName = open_dialog("Pick a file")
    # end

    SPICE_DF,NodeList,InputList,NumVSources = SPICE2Matrix(FileName)
    # println(InputList)
    if inputs===nothing
        # println("No inputs given")
        inputs = zeros(length(InputList))
        inputs[end-(NumVSources-1):end] .= (1*InputScaling)
    end

    ResultNodeNames = vcat("V(".*NodeList.*")", "I(".*(InputList[end-(NumVSources-1):end]).*")")


    Results = [ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results = hcat(Results...)
    ResDict = Dict(ResultNodeNames[1] => Results[1,:])
    for i in 2:length(ResultNodeNames)

        merge!(ResDict,Dict(ResultNodeNames[i] => Results[i,:]))

    end
    CurVec = plotACElCurrent(SPICE_DF,FreqList,Results,ComponentName)
    plot(FreqList,abs.(CurVec))
    CurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    CurrResults = hcat(CurrResults...)

    CurrentDF = getElementCurrents(SPICE_DF,CurrResults,DriveFreq)
    CurrentElIndex = findfirst(isequal(ComponentName),CurrentDF.Name)

    BaselineCurrent = CurrentDF.Current[CurrentElIndex]


    SPICE_DF.Drift = zeros(length(SPICE_DF.Name))
    SPICE_DF.TempRise = zeros(length(SPICE_DF.Name))
    SPICE_DF.Dissipation = zeros(length(SPICE_DF.Name))
    for i in 1:length(SPICE_DF.Name)

        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²R
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ)
            SPICE_DF.Value[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) * SPICE_DF.Value[i]

        elseif SPICE_DF.Type[i] =='L'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ)
            SPICE_DF.ESR[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) * SPICE_DF.ESR[i]# for inductors, only the ESR will drift, not inductance
        elseif SPICE_DF.Type[i] =='C'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
            SPICE_DF.TempRise[i]    = Θc*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θc*SPICE_DF.Dissipation[i]*TCc)
            SPICE_DF.Value[i]       = (1 + Θc*SPICE_DF.Dissipation[i]*TCc) * SPICE_DF.Value[i]

        end
    end

    NewCurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    NewCurrResults = hcat(NewCurrResults...)

    NewCurrentDF = getElementCurrents(SPICE_DF,NewCurrResults,DriveFreq)
    PostHeatCurrent = NewCurrentDF.Current[CurrentElIndex]

    MagDriftPercent = abs((PostHeatCurrent - BaselineCurrent) / PostHeatCurrent)*100

    Results2 =[ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results2 = hcat(Results2...)
    ResDict2 = Dict(ResultNodeNames[1] => Results2[1,:])
    for i in 2:length(ResultNodeNames)

        merge!(ResDict2,Dict(ResultNodeNames[i] => Results2[i,:]))

    end
    CurVec2 = plotACElCurrent(SPICE_DF,FreqList,Results2,ComponentName)
    plot(FreqList,abs.(CurVec2))


    return SPICE_DF,MagDriftPercent,NewCurrentDF
end



"""
This function takes in the SPICE dataframe, results vector(Nx1), and the frequency to test it at.
        This is needed because the Results vector is in terms of node voltages (and currents in voltage sources) but often the current in a single component is useful to know.

        This gets the voltage across two nodes and divides by the impednace of a component


"""
function getElementCurrents(SPICE_DF,Results,Freq)
    NumPasives = sum(SPICE_DF.Type .!= 'V')
    CurrentDF = DataFrame((Name = Any[], ΔV = Any[], Z = Any[],  Current = Any[]))
    for i in 1:length(SPICE_DF.Name)
        if SPICE_DF.Type[i] !='V'
            N1 = SPICE_DF.Node1[i]
            N2 = SPICE_DF.Node2[i]
            if (N1+N2)<999 #if both are non-GND
                ΔV  = Results[N1] - Results[N2]
            elseif N1==999
                ΔV  = Results[N2]
            else
                ΔV  = Results[N1]
            end

        end
        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        elseif SPICE_DF.Type[i] =='L'
            Z = 2 .* π .* Freq.*im.*SPICE_DF.Value[i] .+ SPICE_DF.ESR[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        elseif SPICE_DF.Type[i] =='C'
            Z = 1 ./( 2 .* π .* Freq.*im.*SPICE_DF.Value[i]) .+ SPICE_DF.ESR[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        end


    end
    return CurrentDF
end


"""
Plot AC Element Current
similar to getElementCurrents, except the freq. is a vector, Results in Nx(FreqList Length), and only looks at a single component's current
"""
function plotACElCurrent(SPICE_DF,FreqList,Results,ComponentName)
    ElIndex = findfirst(isequal(ComponentName),SPICE_DF.Name)




            N1 = SPICE_DF.Node1[ElIndex]
            N2 = SPICE_DF.Node2[ElIndex]
            if (N1+N2)<999 #if both are non-GND
                ΔV  = Results[N1,:] - Results[N2,:]
            elseif N1==999
                ΔV  = Results[N2,:]
            else
                ΔV  = Results[N1,:]
            end


        if SPICE_DF.Type[ElIndex] =='R'
            Z = SPICE_DF.Value[ElIndex]
            CurrentVec = ΔV ./ Z


        elseif SPICE_DF.Type[ElIndex] =='L'
            Z = [2 .* π .* FreqList[i].*im.*SPICE_DF.Value[ElIndex] .+ SPICE_DF.ESR[ElIndex] for i in 1:length(FreqList)]
            CurrentVec = ΔV ./ Z


        elseif SPICE_DF.Type[ElIndex] =='C'
            Z = [1 ./( 2 .* π .* FreqList[i].*im.*SPICE_DF.Value[ElIndex]) .+ SPICE_DF.ESR[ElIndex] for i in 1:length(FreqList)]
            CurrentVec = ΔV ./ Z


        end




    return CurrentVec
end


"""
This function modifies the spice dataframe. Specifically, it updates a value of a component
    DF should be a SPICE Dataframe
    ComponentName should be a string of the name (e.g. ` "LDrive" `)
    NewVal should be a number, and the new value (e.g 0.005 for a 5mH inductor)

"""
function UpdateElementVal!(DF, ComponentName::String,NewVal)
    ElIndex = findfirst(isequal(ComponentName),DF.Name)
    DF.Value[ElIndex] = NewVal
end

"""
This function modifies the spice dataframe. Specifically, it updates a value of a component's ESR
    DF should be a SPICE Dataframe
    ComponentName should be a string of the name (e.g. ` "LDrive" `)
    NewVal should be a number, and the new value (e.g 0.01 for a 10mΩ ESR)

"""
function UpdateElementESR!(DF, ComponentName::String,NewVal)
    ElIndex = findfirst(isequal(ComponentName),DF.Name)
    DF.ESR[ElIndex] = NewVal
end

function GetElementVal(SPICE_DF,Name)
    SPICE_DF.Value[findfirst(isequal(Name),SPICE_DF.Name)]
end

function addPowerDiss_DF!(SPICE_DF,inputs,InputList,DriveFreq;TargetRMSCurrent=35,ComponentName="LDrive")


    CurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    CurrResults = hcat(CurrResults...)

    CurrentDF = getElementCurrents(SPICE_DF,CurrResults,DriveFreq)
    CurrentElIndex = findfirst(isequal(ComponentName),CurrentDF.Name)

    BaselineCurrent = CurrentDF.Current[CurrentElIndex]

    
    UpdatedInputVoltage = TargetRMSCurrent/abs.(BaselineCurrent)

    CurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs*UpdatedInputVoltage]
    CurrResults = hcat(CurrResults...)

    CurrentDF = getElementCurrents(SPICE_DF,CurrResults,DriveFreq)
    

    
    
    SPICE_DF.Dissipation = zeros(length(SPICE_DF.Name))
    for i in 1:length(SPICE_DF.Name)

        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²R
            
        elseif SPICE_DF.Type[i] =='L'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
           
        elseif SPICE_DF.Type[i] =='C'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
           

        end

    end
    return CurrentDF
end




"""
This function takes in:

    DriveFreq, primary frequency of interest in Hz
    Θᵣ = 1, Θc = 1 are the lumped thermal resistance of inductors and resistors ( Θᵣ) and capacitors(Θc) in ᵒC per Watt
    TCᵣ = 0.004, TCc = 0.0003 are the temperature coefficients of the components-- TCᵣis the resistance tempco, and TCc is the capacitor tempco
        Units for TC is in UL per deg C

    InputScaling is the scaling of the input vector, for instance if you use the default input vec it will default to 1V at voltage sources, but this scalar can be used to modify that value

"""
function DetermineComponentsTempCoeffs(SPICE_DF,InputList,NumVSources,DriveFreq = 25e3, ComponentName = "LDrive";
    Θᵣ = 1, Θc = 1,
    TCᵣ = 0.004, TCc = -0.0003,
    InputScaling = 1,
    δppm = 900)

    δ = δppm*1e-6



    
    inputs = zeros(length(InputList))
    inputs[end-(NumVSources-1):end] .= (1*InputScaling)
   

  
   
    
    CurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    CurrResults = hcat(CurrResults...)

    CurrentDF = getElementCurrents(SPICE_DF,CurrResults,DriveFreq)
    CurrentElIndex = findfirst(isequal(ComponentName),CurrentDF.Name)
    
    BaselineCurrent = CurrentDF.Current[CurrentElIndex]


    SPICE_DF.Drift = zeros(length(SPICE_DF.Name))
    SPICE_DF.TempRise = zeros(length(SPICE_DF.Name))
    SPICE_DF.Dissipation = zeros(length(SPICE_DF.Name))
    SPICE_DF.DriftCoeff = zeros(length(SPICE_DF.Name))
    SPICE_DF.DriftCoeffPhase = zeros(length(SPICE_DF.Name))
    for i in 1:length(SPICE_DF.Name)

        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²R
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ)
            dVal = SPICE_DF.Value[i]*δ
            SPICE_DF.Value[i]       = dVal +  SPICE_DF.Value[i]

            NewCurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
            NewCurrResults = hcat(NewCurrResults...)
        
            NewCurrentDF = getElementCurrents(SPICE_DF,NewCurrResults,DriveFreq)
            PostHeatCurrent = NewCurrentDF.Current[CurrentElIndex]
        
            SPICE_DF.DriftCoeff[i] = (abs(PostHeatCurrent) - abs(BaselineCurrent)) / abs(BaselineCurrent)/δ
            SPICE_DF.DriftCoeffPhase[i] = (angle(PostHeatCurrent) - angle(BaselineCurrent))/δ
            SPICE_DF.Value[i]       =  SPICE_DF.Value[i]-dVal

        elseif SPICE_DF.Type[i] =='L'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ)
            dVal = SPICE_DF.ESR[i]*δ
            SPICE_DF.ESR[i]       = dVal +  SPICE_DF.ESR[i]# for inductors, only the ESR will drift, not inductance
            NewCurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
            NewCurrResults = hcat(NewCurrResults...)
        
            NewCurrentDF = getElementCurrents(SPICE_DF,NewCurrResults,DriveFreq)
            PostHeatCurrent = NewCurrentDF.Current[CurrentElIndex]
            SPICE_DF.DriftCoeffPhase[i] = (angle(PostHeatCurrent) - angle(BaselineCurrent))/δ
            SPICE_DF.DriftCoeff[i] = (abs(PostHeatCurrent) - abs(BaselineCurrent)) / abs(BaselineCurrent)/δ

            SPICE_DF.ESR[i]       =  SPICE_DF.ESR[i]-dVal


        elseif SPICE_DF.Type[i] =='C'
            # TempCo of celem PP caps is -300ppm/⁰
            lossTangent = 2e-4
            
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])*abs(CurrentDF.ΔV[TmpIndex])*lossTangent #Power = Total apparent power*tandδ
            SPICE_DF.TempRise[i]    = Θc*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θc*SPICE_DF.Dissipation[i]*TCc)
            dVal = SPICE_DF.Value[i]*δ
            SPICE_DF.Value[i]       = SPICE_DF.Value[i]-dVal #Subtract val because they lower capacitance as heat is applied
            NewCurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
            NewCurrResults = hcat(NewCurrResults...)
        
            NewCurrentDF = getElementCurrents(SPICE_DF,NewCurrResults,DriveFreq)
            PostHeatCurrent = NewCurrentDF.Current[CurrentElIndex]
        
            SPICE_DF.DriftCoeff[i] = (abs(PostHeatCurrent) - abs(BaselineCurrent)) / abs(BaselineCurrent)/δ
            SPICE_DF.DriftCoeffPhase[i] = (angle(PostHeatCurrent) - angle(BaselineCurrent))/δ
            SPICE_DF.Value[i]       =  SPICE_DF.Value[i]+dVal

        end
    end

    
   



    return SPICE_DF
end