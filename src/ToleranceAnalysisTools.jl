
function InitializeTolerance!(DF,Val=0.0)
    DF.Tol .= Val
end

function  UpdateElementTolerance!(DF, ComponentName::String,NewVal)
    ElIndex = findfirst(isequal(ComponentName),DF.Name)
    try
        DF.Tol[ElIndex] = NewVal
    catch
        InitializeTolerance!(DF)
        DF.Tol[ElIndex] = NewVal
    end


end

function  UpdateTypeTolerance!(DF, ComponentType::Char,NewVal)
    ElIndex = findall(isequal(ComponentType),DF.Type)
   
    try
        DF.Tol[ElIndex] .= NewVal
    catch
        InitializeTolerance!(DF)
        DF.Tol[ElIndex] .= NewVal
    end
end

"""
Assigns a 0 or 1 depending on the run number and element index number (should be 0 indexed, so if there are 10 components with tolerances, each component should be numbered 0:9)
...Used in the worst-case analysis
"""
function BinaryVal(currentRun,index)  
    floor(currentRun/(2^index))-2*floor(currentRun/(2^(index+1)))
end

BinaryValTable(MaxRuns, NumInds) = [BinaryVal(c,i) for  i in 0:(NumInds-1),c in 1:MaxRuns]


"""
Assigns a value to a given component for a worst-case tolerance analysis. 

    Inputs

    nom = nominal component value
    tol = tolerance where 0.1 = ±10%
    index = component index. (should be 0 indexed, so if there are 10 components with tolerances, each component should be numbered 0:9)
    numruns = max number of runs 
    currentrun = run that the sovler is currently on
"""
function WorstCaseTol(nom,tol,index,numruns,currentRun)
    if(currentRun==numruns)
        return nom
    elseif BinaryVal(currentRun,index)==1
        return nom*(1+tol)
    else 
        return nom*(1-tol)
    end
end


"""
Gives an updated value perturbed by some component tolerance with a pseudo-Gauss distribution. The distribution is cropped at a maximum spread of the input Tol.

Val = nominal component value
Tol = maximum tolerance, so 0.1 = ±10% maximum deviation

kwarg:
Stdevs = 1, this is the standard deviation division factor of the gaussian distribution in terms of maximum tolerance.  Val*(1+Tol/Stdevs*randn(1)[1])
High values of Stdevs will yield narrow distribution, small values will yield flat distributions
"""
function GaussTol(Val, Tol;Stdevs = 1)
    MaxAllowable = Val*(1+Tol)
    MinAllowable = Val*(1-Tol)
    OutputVal = MaxAllowable
    while ((OutputVal<=MinAllowable)||(OutputVal>=MaxAllowable))
    OutputVal = Val*(1+Tol/Stdevs*randn(1)[1])
    end
    return OutputVal
end

