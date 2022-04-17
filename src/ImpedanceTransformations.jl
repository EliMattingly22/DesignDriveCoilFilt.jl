"""
This function takes in some reactance (at frequency fₒ) and determines the values for a series LC pair that resonates at fᵣ and has the same reactance,X at fₒ

    For future reference, see Everitt, Communication Engineering, P. 267
"""
function findEquivLC(X,fₒ,fᵣ)
    ωₒ = 2*π*fₒ
    ωᵣ = 2*π*fᵣ
    L = ωₒ/(ωₒ^2 - ωᵣ^2)*X
    C = (ωₒ^2 - ωᵣ^2) / (ωₒ*ωᵣ^2*X)
    return L,C
end