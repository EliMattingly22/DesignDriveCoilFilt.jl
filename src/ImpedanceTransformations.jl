"""
This function takes in some reactance (real value, not complex) at frequency fₒ and determines the values for a series LC pair that resonates at fᵣ and has the same reactance,X at fₒ

    For future reference, see Everitt, Communication Engineering, P. 267

    To find the equivalent parallel LC use "findEquivLC_Par"
"""
function findEquivLC(X,fₒ,fᵣ)
    ωₒ = 2*π*fₒ
    ωᵣ = 2*π*fᵣ
    L = ωₒ/(ωₒ^2 - ωᵣ^2)*X
    C = (ωₒ^2 - ωᵣ^2) / (ωₒ*ωᵣ^2*X)
    return L,C
end
"""
This function takes in some reactance (real value, not complex) at frequency fₒ and determines the values for a parallel LC pair that resonates at fᵣ and has the same reactance,X at fₒ

    For future reference, see Everitt, Communication Engineering, P. 267 (it is just a impedance/admittance transformation)

"""
function findEquivLC_Par(X,fₒ,fᵣ)

    C,L = abs.(findEquivLC(-1/(X),fₒ,fᵣ))

    return L,C
end






"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

It outputs a T network to match to that load, although only two of the three elements of the 'T' will be populated
\n

\t   ---CS-------------- *\n
\t           |           |\n
\t          CP           Zₗ\n
\t           |           |\n
\t ----------------------*\n

returns CP,CS
"""
function lumpedElementMatch_CapCap(Zₗ,Z₀,f)

    Rₛ = real(Z₀)
    Xₛ = imag(Z₀)
    Rₗ = real(Zₗ)
    Xₗ = imag(Zₗ)
    ω = 2*pi*f

    a = Rₛ*ω^2*(Rₗ^2+Xₗ^2)
    b = -2*Xₗ*Rₛ*ω
    c = -Rₗ+Rₛ


#    x1 = (Rₛ*Xₗ*ω - sqrt(-1*(Rₛ*Rₗ*ω)^2+Rₛ*Xₗ^2*Rₗ*ω^2+Rₛ*Rₗ^3*ω^2))/(Rₛ*ω^2*(Rₗ^2+Xₗ^2))

    CP = (-b-sqrt(b^2-4*a*c))/(2*a)
    
    X = imag(Par(Z_Cap(CP,f),Zₗ))+Xₛ
    CS = 1/(ω*X)

    return CP,CS

end
    













"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

It outputs a T network to match to that load, although only two of the three elements of the 'T' will be populated
\n

\t jX[First] ----------- jX₂[second]--- *\n
\t                 |                    |\n
\t                jB                    Zₗ\n
\t                 |                    |\n
\t -------------------------------------*\n


"""
function lumpedElementMatch(Zₗ,Z₀)

    Xₗ = imag(Zₗ)
    Rₗ = real(Zₗ)
    if Rₗ>real(Z₀)
        B1 = (Xₗ + sqrt(Rₗ/Z₀)*sqrt(Rₗ^2 + Xₗ^2 - Z₀*Rₗ)) /(Rₗ^2 + Xₗ^2)
        B2 = (Xₗ - sqrt(Rₗ/Z₀)*sqrt(Rₗ^2 + Xₗ^2 - Z₀*Rₗ)) /(Rₗ^2 + Xₗ^2)
        
        X1 = (B1*Z₀*Rₗ-Xₗ) / (1-B1*Xₗ)
        X2 = (B2*Z₀*Rₗ-Xₗ) / (1-B2*Xₗ)
        SerElLoc = "First"

    else
    
        B1 = 1/Z₀ * sqrt((Z₀-Rₗ)/Rₗ)
        B2 = -1/Z₀ * sqrt((Z₀-Rₗ)/Rₗ)
        
        X1 = sqrt(Rₗ*(Z₀-Rₗ)) - Xₗ
        X2 = -1*sqrt(Rₗ*(Z₀-Rₗ)) - Xₗ
        SerElLoc = "Second"
    end


    return B1,B2,X1,X2,SerElLoc

end

"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

It outputs a T network to match to that load, although only two of the three elements of the 'T' will be populated
\n

\t   ---SerEl------ SerEl---- *\n
\t                |           |\n
\t             ShuntEl        Zₗ\n
\t                |           |\n
\t ---------------------------*\n

returns ShuntEl, SerEl

"""
function lumpedElementMatch(Zₗ,Z₀,f;ShuntType = 'C')

    B1,B2,X1,X2,SerElLoc = lumpedElementMatch(Zₗ,Z₀)
    ω = 2*pi*f
    ShuntEl = 0.0

    if B1>=0
        if ShuntType=='C'
            ShuntEl = B1/(ω)
            SerEl = X1/(ω)
        elseif ShuntType=='L'
            ShuntEl = -1/(B2*ω)
            SerEl = -1/(X2*ω)
        else
            error("ShuntType must be L or C")
        end
    else
        if ShuntType=='C'
            ShuntEl = B2/(ω)
            SerEl = X2/(ω)
        elseif ShuntType=='L'
            ShuntEl = -1/(B1*ω)
            SerEl = -1/(X1*ω)
        else
                error("ShuntType must be L or C")
        end
        
    end


    return ShuntEl,SerEl,SerElLoc

end



# The following functions use equations from: https://home.sandiego.edu/~ekim/e194rfs01/jwmatcher/matcher2.html
# though they are found elsewhere as well. 

"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

It outputs a T network to match to that load, although only two of the three elements of the 'T' will be populated
\n

\t   -----Cs------ LL---- *\n
\t    |       |           |\n
\t    LS      CL          Zₗ\n
\t    |       |           |\n
\t -----------------------*\n

returns Ls,Cs,CL,LL,q

"""
function lumpedElementMatch_BP1(Zₗ,Z₀,f)
    Rₛ = real(Z₀)
    Xₛ = imag(Z₀)
    Rₗ = real(Zₗ)
    Xₗ = imag(Zₗ)
    
    q=Xₛ/Rₛ
    ω = 2*pi*f

    Rₚ = (1+q^2)*Rₛ
    Rᵥ = sqrt(Rₚ*Rₗ)
            
 
            if Rₚ>Rᵥ 
                qₛ=sqrt(Rₚ /Rᵥ -1)
                ql=sqrt(Rᵥ /Rₗ -1)
                lp=Rₚ / ω / q
                Cs=1 / ω / Rᵥ / qₛ
                
                Ls = Rₚ/ω/qₛ
                if Xₛ!=0
                        if lp==Ls
                            LS = inf  
                        else  
                            LS=1*lp/(lp-l)
                        end
                    end
                
                LL=ql*Rₗ/ω-Xₗ/ω
                
                CL=ql/ω/Rᵥ
                
            else
                println("Impedance is stepping wrong direction")
                Ls = NaN
                Cs = NaN
                CL = NaN
                LL = NaN
            end
return Ls,Cs,CL,LL,q
end        


"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

\n

\t --LS------CL---------- *\n
\t        |       |       |\n
\t       CS       LL      Zₗ\n
\t        |       |       |\n
\t -----------------------*\n

returns Ls,Cs,CL,LL,q

"""
function lumpedElementMatch_BP2(Zₗ,Z₀,f)
    Rₛ = real(Zₗ)
    Xₛ = imag(Zₗ)
    Rₗ = real(Z₀)
    Xₗ = imag(Z₀)
    q=Xₛ/Rₛ
    ω = 2*pi*f
    Rₚ = (1+q^2)*Rₛ
    Rᵥ = sqrt(Rₚ*Rₗ)
            
 
            if Rₚ>Rᵥ 
                qₛ=sqrt(Rₚ /Rᵥ -1)
                ql=sqrt(Rᵥ /Rₗ -1)
                lp=Rₚ / ω / q
                Cs=1 / ω / Rᵥ / qₛ
                
                Ls = Rₚ/ω/qₛ
                if Xₛ!=0
                        if lp==Ls
                            LS = inf  
                        else  
                            LS=lp/(lp-l)
                        end
                    end
                
                LL=ql*Rₗ/ω-Xₗ/ω
                
                CL=ql/ω/Rᵥ
                
            else
                
                println("Impedance is stepping wrong direction")
                
                Ls = NaN
                Cs = NaN
                CL = NaN
                LL = NaN
            end
return Ls,Cs,CL,LL,q
end        

"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

\n

\t   -----Ls------ CL---- *\n
\t    |       |           |\n
\t    CS      LL          Zₗ\n
\t    |       |           |\n
\t -----------------------*\n

returns Ls,Cs,CL,LL,q

"""
function lumpedElementMatch_BP3(Zₗ,Z₀,f)
    Rₛ = real(Z₀)
    Xₛ = imag(Z₀)
    Rₗ = real(Zₗ)
    Xₗ = imag(Zₗ)

    q=-1*Xₛ/Rₛ
    ω = 2*pi*f
    Rₚ = (1+q^2)*Rₛ
    Rᵥ = sqrt(Rₚ*Rₗ)
            
 
            if Rₚ>Rᵥ 
                qₛ=sqrt(Rₚ  /Rᵥ -1)
                ql=sqrt(Rᵥ /Rₗ  -1)
                Cₚ=q / ω / Rₚ
                Cs=qₛ / ω / Rₚ - Cₚ
                
                Ls = qₛ*Rᵥ/ω
                
                
                LL=Rᵥ/ω/ql
                Cₛ = -1/ω/Xₗ
                CL = 1/ω/Rₗ/ql
                    if Xₗ==0
                        CL = Inf
                    else
                        CL=CL*Cₛ/(Cₛ-CL)
                    end
            else
                
                println("Impedance is stepping wrong direction")
                
                Ls = NaN
                Cs = NaN
                CL = NaN
                LL = NaN
            end
return Ls,Cs,CL,LL,q
end        


"""
Function takes in a load impedance, Zₗ, matching impedance,Z₀, and a frequency

\n

\t --CS------LL---------- *\n
\t        |       |       |\n
\t       LS       CL      Zₗ\n
\t        |       |       |\n
\t -----------------------*\n

returns Ls,Cs,CL,LL,q

"""
function lumpedElementMatch_BP4(Zₗ,Z₀,f)
    Rₛ = real(Zₗ)
    Xₛ = imag(Zₗ)
    Rₗ = real(Z₀)
    Xₗ = imag(Z₀)
    q=Xₛ/Rₛ
    ω = 2*pi*f
    Rₚ = (1+q^2)*Rₛ
    Rᵥ = sqrt(Rₚ*Rₗ)
    if Rₚ>Rᵥ 
        qₛ=sqrt(Rₚ  /Rᵥ -1)
        ql=sqrt(Rᵥ /Rₗ  -1)
        Cₚ=q / ω / Rₚ
        Cs=qₛ / ω / Rₚ - Cₚ
        
        Ls = qₛ*Rᵥ/ω
        
        
        LL=Rᵥ/ω/ql
        Cₛ = -1/ω/Xₗ
        CL = 1/ω/Rₗ/ql
            if Xₗ==0
                CL = Inf
            else
                CL=CL*Cₛ/(Cₛ-CL)
            end
    else
        
        println("Impedance is stepping wrong direction")
        
        Ls = NaN
        Cs = NaN
        CL = NaN
        LL = NaN
    end
return Ls,Cs,CL,LL,q
end        
