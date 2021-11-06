



function SPICE2Matrix(FileName::String)
    SPICE_DF = LTSpiceLoad(FileName)
    NodeList,InputList,NumVSources = ProcessSPICE_DF(SPICE_DF)
    return SPICE_DF,NodeList,InputList,NumVSources
end






function ProcessSPICE_DF(DF::DataFrame)

    NodeList = []
    NodeList = reverse(sort(unique(vcat(DF[:,3],DF[:,4])))) #picking out unique nodes, and making 0 at the end
    

    InputList = NodeList;
    InputList = convert(Array{String,1}, InputList)
    InputList = "Iin_".*InputList

    NumVSources=0
    for k in 1:length(DF.Type)
           if DF.Type[k]=='V'
               NumVSources = NumVSources+1
               push!(InputList, "Vin_".*DF.Name[k])
           end
    end
    for j in 1:length(NodeList)
        
            DF[:,3] = replace(DF[:,3],NodeList[j]=>j)
            DF[:,4] = replace(DF[:,4],NodeList[j]=>j)
        
    end
    DF[:,3] = replace(DF[:,3],length(NodeList)=>999) #Making GND node 999
    DF[:,4] = replace(DF[:,4],length(NodeList)=>999)

    filter!(x-> x!="0",NodeList) #Removing GND node
    filter!(x-> x!="Iin_0",InputList)

return NodeList,InputList,NumVSources
end

function addG!(G,Value,N1,N2)
    
    Z=Value
    if N1<999
        G[N1,N1] = G[N1,N1]+ 1 ./Z
    end
    if N2<999
        G[N2,N2] = G[N2,N2]+ 1 ./Z
    end
    if ((N1<999)&(N2<999))
        G[N1,N2] = G[N1,N2]- 1 ./Z
        G[N2,N1] = G[N2,N1]- 1 ./Z
    end
    
    return G
end

function SPICE_DF2Matrix_ω(DF,ω,InputList)
      
        

    Y = complex(zeros(length(InputList),length(InputList)))
    
    SrcMat = zeros(size(Y))


 for i in 1:size(DF,1)
     N1 = DF.Node1[i]
     N2 = DF.Node2[i]
     if (DF.Type[i]=='R')|(DF.Type[i]=='L')|(DF.Type[i]=='C')
         #println(i)


         if DF.Type[i]=='R'
             Z = DF.Value[i]
            
             addG!(Y,Z,N1,N2)
         elseif DF.Type[i]=='L'
             Z = im*ω*DF.Value[i]
             if DF.ESR[i]>0
                 Z += DF.ESR[i]
             end
             addG!(Y,Z,N1,N2)
         elseif DF.Type[i]=='C'
            Z = 1 /(im*ω*DF.Value[i])
            if DF.ESR[i]>0
                Z += DF.ESR[i]
            end
            addG!(Y,Z,N1,N2)
         end



     elseif DF.Type[i]=='V'

         VecLoc = findfirst( x -> x=="Vin_"*DF.Name[i], InputList)

         if N1<999
             SrcMat[VecLoc,N1] = 1
             SrcMat[N1,VecLoc] = 1

         end
         if N2<999

             SrcMat[VecLoc,N2] = -1
             SrcMat[N2,VecLoc] = -1
         end

     end

 end
 FullMat = Y .+ SrcMat
 
 return FullMat

end