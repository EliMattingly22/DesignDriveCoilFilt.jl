include("FieldCalc.jl")
# BiotSav(PointPath,r;Current=1,MinThreshold = 0.01)
# DCore_PointPath(r₁, r₂; dr = (r₂-r₁)/100)
# MakeEllip(r₁,r₂;Center = [0,0,0],NPts=100)
# VecDist(X::Array) = [ X[2,1]-X[1,1], X[2,2]-X[1,2] , X[2,3]-X[1,3]]




"""
This function takes in the following (required) variables:
    IRad: Inner radius of the D-shaped toroid in meters
    ORad: Outer radius of the D-shaped toroid in meters (α*IRad)
    N: Number of (evenly distributed) turns

    kwargs:
    DownSampleFac: this is the amount the D-core point path will be down-sampled  to make the ``Test points''
        It is more accurate to have a highly sampled point path and less finely sampled interior. 
    PlotsOn is a boolean that, if true, plots the wires as they are simulated
    NPtsPath is the number of sections on the point path
    NLayers is the number of interior test point layers
    WireRadius is passed to Mutual_L_TwoLoops and (slightly) changes the way test points are found

    
"""
function eval_Induct_DToroid(IRad,ORad,N;DownSampleFac = 1,PlotOn=false,NPtsPath = 100,NLayers = 20,WireRadius=nothing)
    ΘperTurn = 360/N
    if PlotOn
        pygui(true)
        gcf()
    end
    PointPath_SingleLoop =  DCore_PointPath(IRad, ORad;NPts = NPtsPath)
    SliceArea = SymmetricPP2SliceArea(PointPath_SingleLoop)
    LMat = zeros(N,N)
    if PlotOn
        scatter3D(PointPath_SingleLoop[:,1],PointPath_SingleLoop[:,2],PointPath_SingleLoop[:,3])
    end
    NewPP = Rot_PointPath_Y(PointPath_SingleLoop, ΘperTurn)
    if PlotOn
        scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3])
    end
    Mut,Self,SavedBSelfArr = Mutual_L_TwoLoops(PointPath_SingleLoop,NewPP;DownSampleFac = DownSampleFac,MeasLayers = NLayers,MinThreshold=1e-10,WireRadius=WireRadius,IncludeWireInduct = true,SaveΦ₁ = true)
    LMat .+= Self.*OffDiagOnes(N,1)
    LMat .+= Mut.*OffDiagOnes(N,2)
    LMat .+= Mut.*OffDiagOnes(N,N)
    for i in 2:Int(ceil((N-1)/2))
        NewPP = Rot_PointPath_Y(PointPath_SingleLoop, i*ΘperTurn)
        if PlotOn
            scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3])
        end
        Mut = Mutual_L_TwoLoops(PointPath_SingleLoop,NewPP,SavedBSelfArr;DownSampleFac = DownSampleFac,MeasLayers = NLayers,MinThreshold=1e-10,WireRadius=WireRadius)
        println(Mut)
        
        LMat .+= Mut.*OffDiagOnes(N,i+1)
        LMat .+= Mut.*OffDiagOnes(N,N-(i-1))
        
        println(i)
    end
    LMat[N,N] = LMat[1,1]
    return LMat,sum(LMat)
end
"""
This function takes in the following (required) variables:
    IRad: Inner radius of the circular-core toroid in meters
    ORad: Outer radius of the circular-core toroid in meters (α*IRad)
    N: Number of (evenly distributed) turns

    kwargs:
    DownSampleFac: this is the amount the D-core point path will be down-sampled  to make the ``Test points''
        It is more accurate to have a highly sampled point path and less finely sampled interior. 
    PlotsOn is a boolean that, if true, plots the wires as they are simulated
    NPtsPath is the number of sections on the point path
    NLayers is the number of interior test point layers
    WireRadius is passed to Mutual_L_TwoLoops and (slightly) changes the way test points are found
    

"""
function eval_Induct_CircToroid(IRad,ORad,N;DownSampleFac = 1,
    PlotOn=false,
    NPtsPath = 100,
    NLayers = 20,WireRadius=nothing)
    ΘperTurn = 360/N
    if PlotOn
        pygui(true)
        gcf()
    end
    CoreRad = (ORad-IRad)/2
    CoreOffset = (ORad+IRad)/2
    PointPath_SingleLoop =  MakeEllip(CoreRad,CoreRad;Center =[CoreOffset,0,0], NPts = NPtsPath)
    SliceArea = SymmetricPP2SliceArea(PointPath_SingleLoop)
    LMat = zeros(N,N)
    if PlotOn
        scatter3D(PointPath_SingleLoop[:,1],PointPath_SingleLoop[:,2],PointPath_SingleLoop[:,3])
    end

    for i in 1:(N-1)
        NewPP = Rot_PointPath_Y(PointPath_SingleLoop, i*ΘperTurn)
        if PlotOn
            scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3])
        end
        Mut,Self = Mutual_L_TwoLoops(PointPath_SingleLoop,NewPP;DownSampleFac = DownSampleFac,MeasLayers = NLayers,MinThreshold=1e-10,WireRadius=WireRadius,IncludeWireInduct = true)
        println(Mut)
        LMat[i,i] = Self
        LMat .+= Mut.*OffDiagOnes(N,i+1)
        
        println(i)
    end
    LMat[N,N] = LMat[1,1]
    return LMat,sum(LMat)
end

"""
This function takes in the following (required) variables:
    IRad: Inner radius of the D-shaped toroid in meters
    ORad: Outer radius of the D-shaped toroid in meters (α*IRad)
    N₁: Number of (evenly distributed) turns on the primary
    N₂: Number of (evenly distributed) turns on the secondary

    kwargs:
    DownSampleFac: this is the amount the D-core point path will be down-sampled  to make the ``Test points''
        It is more accurate to have a highly sampled point path and less finely sampled interior. 
    PlotsOn is a boolean that, if true, plots the wires as they are simulated
    NPtsPath is the number of sections on the point path
    NLayers is the number of interior test point layers
    WireRadius is passed to Mutual_L_TwoLoops and (slightly) changes the way test points are found
    SelfInductScaleFactor is, as the name implies, a small ``fudge factor'' if you want to artificially scale up the self-inductance. This may be useful if there is a known offset to account for. 

"""
function eval_Induct_DToroidTransformer(IRad,ORad,N₁,N₂;DownSampleFac = 1,PlotOn=false,NPtsPath = 100,NLayers = 20,SelfInductScaleFactor = 1,UseAnalyticalSelfL=false)
    if PlotOn
        pygui(true)
        gcf()
    end

    N = N₁+N₂

    N₁List = 0:(N/N₁):(N-N/N₁)
    N₂List = 0.0001:(N/N₂):(N-N/N₂+0.0001)
    
    NIndexingArray = hcat(vcat(N₁List,N₂List), vcat(ones(length(N₁List),1),zeros(length(N₂List),1)))
    NIndexingArray= sortslices(NIndexingArray;dims=1) #This creates an array that will allow for determining which turn is associated with primary or secondary
    println(NIndexingArray)
    ΘperTurn = 360/N
    
    PointPath_SingleLoop =  DCore_PointPath(IRad, ORad;NPts = NPtsPath)
    SliceArea = SymmetricPP2SliceArea(PointPath_SingleLoop)
    LMat = zeros(N,N)
    if PlotOn
        scatter3D(PointPath_SingleLoop[:,1],PointPath_SingleLoop[:,2],PointPath_SingleLoop[:,3];color="b")
    end
    NewPP = Rot_PointPath_Y(PointPath_SingleLoop, ΘperTurn)
    Mut,Self,SavedBSelfArr = Mutual_L_TwoLoops(PointPath_SingleLoop,NewPP;DownSampleFac = DownSampleFac,MeasLayers = NLayers,MinThreshold=1e-10,IncludeWireInduct = true,SaveΦ₁ = true)
    LMat .+= Self.*OffDiagOnes(N,1)
    LMat .+= Mut.*OffDiagOnes(N,2)
    LMat .+= Mut.*OffDiagOnes(N,N)
    for i in 1:Int(ceil((N-1)/2))
        isPrimary = (NIndexingArray[i+1,2]==1)
        NewPP = Rot_PointPath_Y(PointPath_SingleLoop, i*ΘperTurn)
        if PlotOn
            if isPrimary 
                scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3];color="b",s=2)
            else
                scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3];color="r",s=2)
            end
        end
        Mut = Mutual_L_TwoLoops(PointPath_SingleLoop,NewPP,SavedBSelfArr;DownSampleFac = DownSampleFac,MeasLayers = NLayers,MinThreshold=1e-10)
        println(Mut)
        
        LMat .+= Mut.*OffDiagOnes(N,i+1)
        LMat .+= Mut.*OffDiagOnes(N,N-(i-1))
        
        println(i)
    end
    LMat[N,N] = LMat[1,1]
    Primary_Self, Secondary_Self,MutualInduct = LMatMutualIndexing(N,NIndexingArray[:,2])
    L₁₁ = LMat[Primary_Self] .*SelfInductScaleFactor
    L₂₂ = LMat[Secondary_Self] .*SelfInductScaleFactor
    LMut = LMat[MutualInduct]
    if UseAnalyticalSelfL
        println("Using analyical formula for self inductance. Make sure ORad/IRad is an integer")
        L₁₁ = DCore_DetermineIdealInduct(IRad*2,N₁,Int(round(ORad/IRad)))
        L₂₂ = DCore_DetermineIdealInduct(IRad*2,N₂,Int(round(ORad/IRad)))
    end
    Lμ = sum(LMut)*N₁/N₂ / 2 #THe division by 2 is due to the double-counting of M in a symmetric matrix
    LLeak1 = sum(L₁₁) - Lμ
    LLeak2 = sum(L₂₂) - Lμ*(N₂/N₁)^2
    
    K = (sum(LMut)/2) / √(sum(L₁₁)*sum(L₂₂))
    println("L₁₁ = $(round(sum(L₁₁);sigdigits=3))")
    println("L₂₂ = $(round(sum(L₂₂);sigdigits=3))")
    println("Lₘ = $(round(sum(LMut)/2;sigdigits=3))")
    println("K = $(round(sum(K);sigdigits=3))")
    println("Lμ = $(round(Lμ;sigdigits=3))")
    println("LLeak1 = $(round(LLeak1;sigdigits=3))")
    println("LLeak2 = $(round(LLeak2;sigdigits=3))")
    return LMat, NIndexingArray,L₁₁,L₂₂,LMut,Lμ,LLeak1,LLeak2,K
end


"""
This builds a rotation matrix, to rotate around the z axis
Radians input.
"""
function BuildRotMatrix_deg_Z(Θ)
    Θ = Θ*π/180
   Mat =  [cos(Θ) -1*sin(Θ) 0
           sin(Θ)  cos(Θ)   0
           0       0        1]
end
"""
This builds a rotation matrix, to rotate around the Y axis
Radians input.
"""
function BuildRotMatrix_deg_Y(Θ)
    Θ = Θ*π/180
   Mat =  [cos(Θ)  0   -1*sin(Θ) 
           0       1   0
           sin(Θ)       0   cos(Θ)]
end


function Rot_PointPath_Y(Point_Path, Θ)
    Point_PathTrans = transpose(Point_Path)
    return transpose(BuildRotMatrix_deg_Y(Θ)*Point_PathTrans)
end

function Rot_PointPath_Z(Point_Path, Θ)
    Point_PathTrans = transpose(Point_Path)
    return transpose(BuildRotMatrix_deg_Z(Θ)*Point_PathTrans)
end

"""
This function takes in two point paths and simulates the self inductance of the first, and mutual inductance between them.

    kwargs:
    DownSampleFac: this is the amount the point path will be down-sampled  to make the ``Test points''
        It is more accurate to have a highly sampled point path and less finely sampled interior. 
    
    
    MeasLayers is the number of interior test point layers
    WireRadius (slightly) changes the way test points are found. It mostly just adds a set of test points near the wire if the distance is large between the outer layer of test points and wire edge
    IncludeWireInduct is a boolean. If true, adds the wire's internal inductance to the self-inductance of the first
    MinThreshold should not be modified to be larger than the wire radius. It essentially makes the bior savart function ignore test points under that treshold, which can be useful to prevent divisions by zero in edge cases.

    

"""
function Mutual_L_TwoLoops(PP₁,PP₂;DownSampleFac = 1,
    MeasLayers = 55,#Number of concentric layers to calc field
    MinThreshold = 0.000001,
    WireRadius = 0.001,
    IncludeWireInduct = false,
    SaveΦ₁ = false
    )

    TestPoints, Weights = PP_2_TestPoints_drdΘ(PP₁,MeasLayers,WireRadius;DownSampleFac = DownSampleFac)
    PtsPerLayer = length(PP₁[:,1])
    NPtsPerSlice = MeasLayers*PtsPerLayer  
    SliceArea = sum(Weights) 
        Weights = Weights / sum(Weights)*NPtsPerSlice  #normalizing
        if SaveΦ₁
            B_Self = zeros(length(TestPoints[:,1]) ,3)
            Φ₁₁ = 0.
        else
            Φ₁₁ = 0.
        end
    
    Φ₂₁ = 0.
    for i in 1:length(TestPoints[:,1]) 
         
        if SaveΦ₁
            B_Self[i,:] =  BiotSav(PP₁,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice
            Φ₁₁ += √( sum( B_Self[i,:].^2 ) )
            Φ₂₁ += sum(B_Self[i,:]  .* BiotSav(PP₂,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice) /    √(sum(B_Self[i,:] .^2))
        else
            B_Self =  BiotSav(PP₁,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice
            Φ₁₁ += √(sum(B_Self.^2))
            
            Φ₂₁ += sum(B_Self .* BiotSav(PP₂,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice) /    √(sum(B_Self.^2))
        end
        
       
    end

    CumulativeDist=0
    if IncludeWireInduct
        for i in 1:(length(PP₁[:,1])-1)
            CumulativeDist += sqrt(sum((PP₁[i,:] .- PP₁[i+1,:]).^2))
        end
        println(CumulativeDist)
        WireInductance = internalInductance(CumulativeDist)
        Φ₁₁ += WireInductance
    end

    return Φ₂₁, Φ₁₁,B_Self

end


function Mutual_L_TwoLoops(PP₁,PP₂,PriorBSelf;DownSampleFac = 1,
    MeasLayers = 55,#Number of concentric layers to calc field
    MinThreshold = 0.000001,
    WireRadius = 0.001,
   
    )

    TestPoints, Weights = PP_2_TestPoints_drdΘ(PP₁,MeasLayers,WireRadius;DownSampleFac = DownSampleFac)
    PtsPerLayer = length(PP₁[:,1])
    NPtsPerSlice = MeasLayers*PtsPerLayer  
    SliceArea = sum(Weights) 
        Weights = Weights / sum(Weights)*NPtsPerSlice  #normalizing
    
    Φ₂₁ = 0.
    for i in 1:length(TestPoints[:,1]) 
        Φ₂₁ += sum(PriorBSelf[i,:] .* BiotSav(PP₂,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice) /    √(sum(PriorBSelf[i,:].^2))
       
    end

    

    return Φ₂₁

end




function SymmetricPP2SliceArea(PP) 
Total = 0
    for i in 1:(Int(round(length(PP[:,1])/2)))
        dX = PP[i+1,1] - PP[i,1]
        MeanY = (PP[i+1,2] + PP[i,2])/2
        Total += dX*MeanY
    end

    return Total*2 #integrating over half, then multiplying by 2 to account for symmetry 
end


function OffDiagOnes(FullSize,Dist)
    M = zeros(FullSize,FullSize)
    j=0
    for i in Dist:(FullSize)
        j +=1
        M[i,j] = 1
        M[j,i] = 1
    end
    return M
end

function LMatMutualIndexing(NTotal,IndexingVector_Primary)

    Primary_Self = zeros(NTotal,NTotal)
    Secondary_Self = zeros(NTotal,NTotal)
    MutualInduct = zeros(NTotal,NTotal)

    for i in 1:NTotal, j in 1:NTotal
        if ((IndexingVector_Primary[i]==1)&(IndexingVector_Primary[j]==1))
            Primary_Self[i,j] = 1
        elseif ((IndexingVector_Primary[i]==0)&(IndexingVector_Primary[j]==0))
            Secondary_Self[i,j] = 1
        end
    end
    MutualInduct = (Primary_Self.+Secondary_Self).==0
    Primary_Self = Primary_Self.==1
    Secondary_Self = Secondary_Self.==1
    return Primary_Self, Secondary_Self,MutualInduct

end





function Mutual_L_CircTwoLoops_Separation(D,SeparationVec  ;DownSampleFac=1,
    MeasLayers = 55,#Number of concentric layers to calc field
    MinThreshold = 1e-15, NPts = 100,PlotOn = false
    )

    PP₁ = MakeEllip(D/2,D/2;NPts = NPts)
    if PlotOn
        scatter3D(PP₁[:,1],PP₁[:,2],PP₁[:,3];color="r")
    end
    SliceArea = π/4*D^2

    MutualLMat = zeros(length(SeparationVec))
    Counter = 0
    for j in SeparationVec
        Counter+=1
        PP₂ = MakeEllip(D/2,D/2;Center = [0,0,j])
        if PlotOn
            scatter3D(PP₂[:,1],PP₂[:,2],PP₂[:,3];color="b")
        end
        TestPoints, Weights = PP_2_TestPoints_drdΘ(PP₁,MeasLayers;DownSampleFac = DownSampleFac)
        PtsPerLayer = length(PP₁[:,1])
        NPtsPerSlice = MeasLayers*PtsPerLayer   
        Weights = Weights / sum(Weights)*NPtsPerSlice  #normalizing
    
        Φ₁₁ = 0.  
        Φ₂₁ = 0.
        for i in 1:length(TestPoints[:,1]) 
            B_Self =  BiotSav(PP₁,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice 
            Φ₁₁ += √(sum(B_Self.^2))
            Φ₂₁ += sum(B_Self .* BiotSav(PP₂,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice) /    √(sum(B_Self.^2))
        
        end
        MutualLMat[Counter] = Φ₂₁
    end

    return MutualLMat

end



function Self_L_CircLoop(D  ; DownSampleFac=1,
    MeasLayers = 55,#Number of concentric layers to calc field
    MinThreshold = 1e-15, NPts = 100,PlotOn = false,WireRadius = 0.0005,IncludeWireInduct = true
    )

    PP₁ = MakeEllip(D/2,D/2;NPts = NPts)
    if PlotOn
        scatter3D(PP₁[:,1],PP₁[:,2],PP₁[:,3];color="r")
    end
    SliceArea = π/4*D^2

    
        PP₂ = MakeEllip(D/2,D/2)
        if PlotOn
            scatter3D(PP₂[:,1],PP₂[:,2],PP₂[:,3];color="b")
        end
        DistBetweenLayers = D/2 / MeasLayers
        if DistBetweenLayers>WireRadius
            TestPoints, Weights = PP_2_TestPoints_drdΘ(PP₁,MeasLayers,WireRadius;DownSampleFac = DownSampleFac)
        else
            TestPoints, Weights = PP_2_TestPoints_drdΘ(PP₁,MeasLayers,WireRadius;DownSampleFac = DownSampleFac)
        end
        PtsPerLayer = length(PP₁[:,1])
        NPtsPerSlice = MeasLayers*PtsPerLayer   
        Weights = Weights / sum(Weights)*NPtsPerSlice  #normalizing
    
        Φ₁₁ = 0.  
        Φ₂₁ = 0.
        for i in 1:length(TestPoints[:,1]) 
            B_Self =  BiotSav(PP₁,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice 
            Φ₁₁ += √(sum(B_Self.^2))
            #Φ₂₁ += sum(B_Self .* BiotSav(PP₂,TestPoints[i,:];MinThreshold) .* Weights[i].*SliceArea ./NPtsPerSlice) /    √(sum(B_Self.^2))
        
        end
       
        
        if IncludeWireInduct
            if DistBetweenLayers<WireRadius
                WireInductance = 0
            else
                WireInductance = internalInductance(π*D)
            end
            println("Adding wire internal inductance")
            println("Distance between layers = : $(round(DistBetweenLayers;sigdigits=3))")
            Φ₁₁ += WireInductance
        end
        
    return Φ₁₁

end

function internalInductance(L)
    μ₀ = 4*pi*1e-7  
    return μ₀*L/(8*π)
end


function makeIndexinVec(BasisVec,Scale)
    NewVec = []
    N = length(BasisVec)
    for i in 1:length(BasisVec), j in 1:Scale
        
        push!(NewVec,BasisVec[i])
    
    end
    NewVec = NewVec[1:N]
    return NewVec
end



function ScatterPlotToroid(IRad,ORad,N₁,N₂,BundleScale;NPtsPath = 100,NLayers = 20)
    
        pygui(true)
        gcf()
    

    N = N₁+N₂

    N₁List = 0:(N/N₁):(N-N/N₁)
    N₂List = 0.0001:(N/N₂):(N-N/N₂+0.0001)
    
    NIndexingArray = hcat(vcat(N₁List,N₂List), vcat(ones(length(N₁List),1),zeros(length(N₂List),1)))
    NIndexingArray= sortslices(NIndexingArray;dims=1) #This creates an array that will allow for determining which turn is associated with primary or secondary
    PrimaryIndex = makeIndexinVec(NIndexingArray[:,2],BundleScale)
    ΘperTurn = 360/N
    
    PointPath_SingleLoop =  DCore_PointPath(IRad, ORad;NPts = NPtsPath)
    
    
        scatter3D(PointPath_SingleLoop[:,1],PointPath_SingleLoop[:,2],PointPath_SingleLoop[:,3];color="b",s=3)
   
    for i in 1:(N-1)
        isPrimary = (PrimaryIndex[i+1]==1)
        NewPP = Rot_PointPath_Y(PointPath_SingleLoop, i*ΘperTurn)
        
            if isPrimary 
                scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3];color="b",s=3)
            else
                scatter3D(NewPP[:,1],NewPP[:,2],NewPP[:,3];color="r",s=3)
            end
        
        
    end
    
    return 
end

"""
This function takes in a point path and produces test points and weights for each. 

    Inputs are a point path and the number of test point layers. A ``layer'' is just a concentric set of points to the point path with a smaller radius. Using this sampling approach gives an uneven distribution of points, but ensures they are all on the interior. 

"""
function PP_2_TestPoints_drdΘ(PP,Layers = 10, WireRadius=0.001; DownSampleFac::Int=1)
    
    PP = PP[1:DownSampleFac:end,:]
    PP_Mean = [sum(PP[:,i]) / length(PP[:,i]) for i in 1:3]
    PP_Cent = hcat([PP[:,1] .- PP_Mean[1] , PP[:,2] .- PP_Mean[2], PP[:,3] .- PP_Mean[3]]...)

    TestPoints_Cent = vcat([PP_Cent .* √(NLayer / (Layers+1))  for NLayer in 1:Layers]...)
    TestPoints_Cent = reverse(TestPoints_Cent;dims=1)
    MaxDist = maximum(sqrt.(sum(PP_Cent.^2;dims=2)))
    TestPoints_Cent = vcat(reverse(PP_Cent .* ((MaxDist-WireRadius) / MaxDist);dims=1),TestPoints_Cent)
    Layers = Layers+1 #Because adding an additional layer near wire
   
    TestPoints_Polar = Cart2Polar(TestPoints_Cent)
    PP_Cent_Polar = Cart2Polar(PP_Cent)
    dΘ = zeros(length(TestPoints_Polar[:,1]),1)
    dr = similar(dΘ)

    # dΘ = abs.(TestPoints_Polar[:,2] .- circshift(TestPoints_Polar[:,2],1))
    # dΘ = minimum(hcat(dΘ,  abs.(TestPoints_Polar[:,2] .- circshift(TestPoints_Polar[:,2],-1)));dims=2)
    L = length(PP_Cent[:,1])
    
    Counter = 0
    for j in 1:Layers
        
        for k in 1:L
            Counter+=1

            if k==1
                dΘ[Counter] = AbsAngBetw(TestPoints_Polar[Counter,2],TestPoints_Polar[Counter+1,2])
                
            elseif k==L
                dΘ[Counter] = AbsAngBetw(TestPoints_Polar[Counter,2],TestPoints_Polar[Counter-1,2])
            else
                dΘ[Counter] = AbsAngBetw(TestPoints_Polar[Counter-1,2],TestPoints_Polar[Counter+1,2])/2
            end

        end
        StartInd = ((j-1)*L+1)
        StopInd = j*L

        StartInd2 = (j*L+1)
        StopInd2 = (j+1)*L

        
        if j !== Layers
            dr[StartInd : StopInd] = TestPoints_Polar[StartInd : StopInd,1] .- TestPoints_Polar[StartInd2 : StopInd2,1]
        else
            dr[StartInd : StopInd] = TestPoints_Polar[StartInd : StopInd,1]
        end
    end

    TestPoints = hcat([TestPoints_Cent[:,1] .+ PP_Mean[1] , TestPoints_Cent[:,2] .+ PP_Mean[2], TestPoints_Cent[:,3] .+ PP_Mean[3]]...)

    Weights = abs.(TestPoints_Polar[:,1] .* dr[:] .* dΘ[:])
    r = TestPoints_Polar[:,1] 
    return TestPoints, Weights,dr,dΘ,r
end

function AbsAngBetw(θ₁,θ₂)
    if θ₁==θ₂
        return 0.
    elseif (abs(θ₁-θ₂)>π)
        return (2*π - mod2pi(abs(θ₁-θ₂)))
    else
        return abs(θ₁-θ₂)
    end
end



function PP_2_TestPoints(PP,Layers = 10;  DownSampleFac::Int=1,WireRadius = nothing)
    PP = PP[1:DownSampleFac:end,:]
    PP_Mean = [sum(PP[:,i]) / length(PP[:,i]) for i in 1:3]
    PP_Cent = hcat([PP[:,1] .- PP_Mean[1] , PP[:,2] .- PP_Mean[2], PP[:,3] .- PP_Mean[3]]...)

    TestPoints_Cent = vcat([PP_Cent .* (NLayer / (Layers+1))  for NLayer in 1:Layers]...)

    if ~(WireRadius===nothing)
        println("Adding inner test layer")
        MaxDist = maximum(sqrt.(sum(PP_Cent.^2;dims=2)))
        TestPoints_Cent = vcat(TestPoints_Cent, PP_Cent .* ((MaxDist-WireRadius) / MaxDist))

    end
    TestPoints = hcat([TestPoints_Cent[:,1] .+ PP_Mean[1] , TestPoints_Cent[:,2] .+ PP_Mean[2], TestPoints_Cent[:,3] .+ PP_Mean[3]]...)


    Weights = [√(TestPoints_Cent[i,1]^2 + TestPoints_Cent[i,2]^2 + TestPoints_Cent[i,3]^2) for i in 1:length(TestPoints_Cent[:,1]) ]

    return TestPoints, Weights
end

function Cart2Polar(Arr)
    PolarArr = similar(Arr)   
    for i in 1:length(Arr[:,1])
        PolarArr[i,1] = √(Arr[i,1]^2 + Arr[i,2]^2)
        PolarArr[i,2] = atan(Arr[i,2],Arr[i,1])
        PolarArr[i,3] = Arr[i,3]
    end

    return PolarArr
end
