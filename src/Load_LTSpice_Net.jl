using Gtk, DataFrames

"""
This function takes in a file path/name and exports a dataframe consisting of the circuit parameters
Specifically the file should be from LTSPICE (or it wont work)

If no input is given, the user will be prompted for a FileName
"""
function LTSpiceLoad(FileName=nothing)
        if FileName===nothing
                FileName = open_dialog("Pick a file")
        end
        f =open(FileName)
         Lines =  readlines(f)

         close(f)

         SPICE_DF = DataFrame((Type = Any[], Value=Any[], Node1=Any[], Node2=Any[], Name=Any[],ESR = Any[]))
        for i in 1:length(Lines)

                Spaces = vcat(findall(" ",Lines[i])...)
                
                if ~(Lines[i][1]=='*')
                        # println(Lines[i])
                end
                if (Lines[i][1]=='R')|(Lines[i][1]=='L')|(Lines[i][1]=='C')
                        if length(Spaces)==3
                                Value = Lines[i][Spaces[3]+1:end]
                        else
                                Value = Lines[i][Spaces[3]+1:Spaces[4]-1]
                        end
                        Node1 = Lines[i][Spaces[1]+1:Spaces[2]-1]
                        Node2 = Lines[i][Spaces[2]+1:Spaces[3]-1]
                        Name = Lines[i][1:Spaces[1]-1]
                        RSerInds = findall("Rser=", Lines[i])
                        if length(RSerInds)>0
                                RSerIndEnd =RSerInds[1][end]+1
                                if RSerInds[1][end]>maximum(Spaces)
                                        RSerVal = MakeNumericalVals(Lines[i][RSerIndEnd:end])
                                else
                                        Tmp = findfirst(x-> x>RSerIndEnd,Spaces)
                                        RSerVal = MakeNumericalVals(Lines[i][RSerIndEnd:Spaces[Tmp]])
                                end
                        else
                                RSerVal=0
                        end
                        push!(SPICE_DF,[Lines[i][1] MakeNumericalVals(Value) Node1 Node2 Name RSerVal])

                elseif (Lines[i][1]=='V')
                        Node1 = Lines[i][Spaces[1]+1:Spaces[2]-1]
                        Node2 = Lines[i][Spaces[2]+1:Spaces[3]-1]
                        Name = Lines[i][1:Spaces[1]-1]
                        push!(SPICE_DF,[Lines[i][1] 1 Node1 Node2 Name 0])

                end


        end

        return SPICE_DF
end


"""
A simple string to numeric conversion for metric prefixes
Input: String 
Output Float64

e.g. 
        MakeNumericalVals("1n")
         > 1.0e-9
        MakeNumericalVals.(["1n" , "2m"])
                2-element Array{Float64,1}:
                1.0e-9
                0.002
"""
function MakeNumericalVals(ValString::String)
        NewString = replace(ValString,"Meg"=>"e6")
        NewString = replace(NewString,"k"=>"e3")
        NewString = replace(NewString,"m"=>"e-3")
        NewString = replace(NewString,"u"=>"e-6")
        NewString = replace(NewString, "n"=>"e-9")
        NewString = replace(NewString, "p"=>"e-12")
        return parse(Float64,NewString)
end
