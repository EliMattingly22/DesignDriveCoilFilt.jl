pygui(true)




function useFilterDesignGUI()
win = GtkWindow("Filter Designer")
g = GtkGrid()

LDriveDesc = GtkEntry()  # a widget for entering text
set_gtk_property!(LDriveDesc, :text, "Load Inductance [H]")
makeDescBox(LDriveDesc)
LDriveVal =GtkEntry()
makeValBox(LDriveVal,100e-6)
g[1,1] = LDriveDesc
g[2,1] = LDriveVal

RDriveDesc = GtkEntry()
set_gtk_property!(RDriveDesc, :text, "Load Resistance [Ω]")
makeDescBox(RDriveDesc)
RDriveVal =GtkEntry()
makeValBox(RDriveVal,0.1)
g[1,2] = RDriveDesc
g[2,2] = RDriveVal

CDriveDesc = GtkEntry()
set_gtk_property!(CDriveDesc, :text, "Load Capacitance [F]")
makeDescBox(CDriveDesc)
CDriveVal =GtkEntry()
makeValBox(CDriveVal,1)
g[1,3] = CDriveDesc
g[2,3] = CDriveVal

TargetZDesc = GtkEntry()
set_gtk_property!(TargetZDesc, :text, "Target Impedance [Ω]")
makeDescBox(TargetZDesc)
TargetZVal =GtkEntry()
makeValBox(TargetZVal,20)
g[1,4] = TargetZDesc
g[2,4] = TargetZVal

DriveFreqDesc = GtkEntry()
set_gtk_property!(DriveFreqDesc, :text, "Drive Freq [Hz]")
makeDescBox(DriveFreqDesc)
DriveFreqVal =GtkEntry()
makeValBox(DriveFreqVal,25000)
g[1,5] = DriveFreqDesc
g[2,5] = DriveFreqVal

FilterZDesc = GtkEntry()
set_gtk_property!(FilterZDesc, :text, "Filter Characteristic Impedance [Ω]")
makeDescBox(FilterZDesc)
FilterZVal =GtkEntry()
makeValBox(FilterZVal,25)
g[1,6] = FilterZDesc
g[2,6] = FilterZVal

DampRValDesc = GtkEntry()
set_gtk_property!(DampRValDesc, :text, "Damping Resistor [Ω]")
makeDescBox(DampRValDesc)
DampRVal =GtkEntry()
makeValBox(DampRVal,50)
g[1,7] = DampRValDesc
g[2,7] = DampRVal

OptimButton = GtkCheckButton("Run Optimization?")
Optim_cb = GtkComboBoxText()
choices = ["Dist to Trough", "Magnitude Drift"]
for choice in choices
  push!(Optim_cb,choice)
end
g[1,8] = OptimButton
g[2,8] = Optim_cb

SaveButton = GtkCheckButton("Save LTSPICE File?")

g[1,9] = SaveButton
SavePathBox = GtkEntry()
set_gtk_property!(SavePathBox, :text, "Not saving")
makeDescBox(SavePathBox)
g[2,9] = SavePathBox

id = signal_connect(SaveButton, "button-press-event") do widget, event
    if ~(get_gtk_property(widget, :active, Bool))
        
        
        set_gtk_property!(SavePathBox, :has_frame,true)
        set_gtk_property!(SavePathBox, :sensitive,true)
        set_gtk_property!(SavePathBox, :editable,true)
        set_gtk_property!(SavePathBox, :text, "FilterDesignOutput.asc")
    else
        set_gtk_property!(SavePathBox, :has_frame,false)
        set_gtk_property!(SavePathBox, :sensitive,false)
        
        set_gtk_property!(SavePathBox, :text, "Not saving") 
    end
end




Notch1Button = GtkCheckButton("Add notch #1")

g[1,10] = Notch1Button
Notch1Val = GtkEntry()
set_gtk_property!(Notch1Val, :text, "No notch")
makeDescBox(Notch1Val)
g[2,10] = Notch1Val

signal_connect(Notch1Button, "button-press-event") do widget, event
    if ~(get_gtk_property(widget, :active, Bool))
        
        
        set_gtk_property!(Notch1Val, :has_frame,true)
        set_gtk_property!(Notch1Val, :sensitive,true)
        set_gtk_property!(Notch1Val, :editable,true)
        set_gtk_property!(Notch1Val, :text, "$(2*getVal(DriveFreqVal))")
    else
        set_gtk_property!(Notch1Val, :has_frame,false)
        set_gtk_property!(Notch1Val, :sensitive,false)
        
        set_gtk_property!(Notch1Val, :text, "Not notching 1") 
    end
end


Notch2Button = GtkCheckButton("Add notch #2")

g[1,11] = Notch2Button
Notch2Val = GtkEntry()
set_gtk_property!(Notch2Val, :text, "No notch")
makeDescBox(Notch2Val)
g[2,11] = Notch2Val

signal_connect(Notch2Button, "button-press-event") do widget, event
    if ~(get_gtk_property(widget, :active, Bool))
        
        
        set_gtk_property!(Notch2Val, :has_frame,true)
        set_gtk_property!(Notch2Val, :sensitive,true)
        set_gtk_property!(Notch2Val, :editable,true)
        set_gtk_property!(Notch2Val, :text, "$(3*getVal(DriveFreqVal))")
    else
        set_gtk_property!(Notch2Val, :has_frame,false)
        set_gtk_property!(Notch2Val, :sensitive,false)
        
        set_gtk_property!(Notch2Val, :text, "Not notching 2") 
    end
end


FindFreqButton = GtkCheckButton("Determine ideal freq?")
g[1,12] = FindFreqButton


RunScriptButton = GtkButton()
set_gtk_property!(RunScriptButton, :label, "RUN DESIGNER")
g[1:2,13] = RunScriptButton


set_gtk_property!(Optim_cb,:active,1)

id = signal_connect(RunScriptButton, "button-press-event") do widget, event
    println("Starting filter designer")

    if  (get_gtk_property(SaveButton, :active, Bool))
            Savepath = get_gtk_property(SavePathBox, :text,String)
            println("The SPICE file will be saved at location:"*get_gtk_property(SavePathBox, :text,String))
            if ~(Savepath[end-3:end]==".asc")
                Savepath = Savepath*".asc"
            end
            set_gtk_property!(SavePathBox, :text, Savepath) 
    else
        Savepath = nothing
        println("Will not save schematic")
    end

    NotchArray = [1e9,1e9]
    if  (get_gtk_property(Notch1Button, :active, Bool))
        NotchArray[1] = getVal(Notch1Val)
    end

    if  (get_gtk_property(Notch2Button, :active, Bool))
        NotchArray[2] = getVal(Notch2Val)
    end


    if  (get_gtk_property(OptimButton, :active, Bool))
        
        MinimizeDistToTrough = parse(Int,get_gtk_property(Optim_cb,:active,String))==0
        
        MinZ,DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList,AxesArray = DesignDriveFilter_OptimDrift(getVal(LDriveVal),getVal(RDriveVal),getVal(TargetZVal),getVal(DriveFreqVal); MinimizeDistToTrough=MinimizeDistToTrough,
        CDrive = getVal(CDriveVal),
        NumDriveElements = 1,
        WireDiam = 2e-3,
        WireFillFac = 0.75,
        PlotSrcR = 0.0001,
        PlotFTs = true,
        VSrc = 2,
        DetermineFreq = (get_gtk_property(FindFreqButton, :active, Bool)),
        AddNotchFreq = NotchArray,
        FilterZ = getVal(FilterZVal),
        RDampVal = getVal(DampRVal),
        PerturbTxReactance = nothing,
        WriteFileName=Savepath
        )
    else
    DriveFreq, CurrentVec, Results, SPICE_DF,inputs,InputList,FreqList,AxesArray = DesignDriveFilter(getVal(LDriveVal),getVal(RDriveVal),getVal(TargetZVal),getVal(DriveFreqVal);
        CDrive = getVal(CDriveVal),
        NumDriveElements = 1,
        WireDiam = 2e-3,
        WireFillFac = 0.75,
        PlotSrcR = 0.0001,
        PlotFTs = true,
        VSrc = 2,
        DetermineFreq = (get_gtk_property(FindFreqButton, :active, Bool)),
        AddNotchFreq = NotchArray,
        FilterZ = getVal(FilterZVal),
        RDampVal = getVal(DampRVal),
        PerturbTxReactance = nothing,
        AxesArray = nothing,
        WriteFileName=Savepath
        )
    end
        println(SPICE_DF)

        if  (get_gtk_property(SaveButton, :active, Bool))
            
            
            Savepath_CSV = Savepath[1:end-4]*".csv"
            println(Savepath_CSV)
            CSV.write(Savepath_CSV, SPICE_DF) 
    
        end
end


set_gtk_property!(g, :column_homogeneous, true)
set_gtk_property!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
push!(win, g)
showall(win)


end

function makeDescBox(Box)
    set_gtk_property!(Box, :has_frame,false)
    set_gtk_property!(Box, :sensitive,false)
    set_gtk_property!(Box, :editable,false)
end

function makeValBox(Box,DefaultVal)
    set_gtk_property!(Box, :text, "$DefaultVal")
end

function getVal(Box)
    parse(Float64,get_gtk_property(Box, :text, String))
end

