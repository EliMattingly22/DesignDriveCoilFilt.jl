using Optim

function PipeFlow(D,P,L;Fluid="Water",
    ΔT = 40,
     FlowType = "Turbulent",
     AbsRoughness = 1.5e-6)
    if Fluid=="Air"
        DynVisc = 20e-6; #Pa*seconds ***Air***
        Dens = 1.225; #***Air***
        C = 1000;#J/(kg*k) AIR
    elseif Fluid=="Oil"
        DynVisc = 5e-3; #Pa*seconds ***Oil***
        Dens = 900; #***Oil***
        C = 1800;#J/(kg*k) AIR
    else()
        C = 4187;#J/(kg*k) WATER
        Dens = 1000;#kg/m^3
        DynVisc = 7.97e-4; #Pa*seconds ***Water***
    end
    
    
    P = P*6894.76 #converting from PSI to Pa
    Re = 1000 #Initial guess
    ϵ = AbsRoughness/D;
    if FlowType == "Laminar"

        f, Re, LinVel = Frict(Re,P,L,D,Dens,DynVisc)

    elseif FlowType == "Turbulent"
        f, Re, LinVel = Frict(ϵ, Re,P,L,D,Dens,DynVisc)
    end
    


    println("Friction factor:$(f)")
    
    println("Re Num:$(Re)")
    
    Area = π/4*D.^2
    MassFlow = Area*LinVel*Dens
    PowerDiss = MassFlow*C*ΔT
    println("Mass Flow:$MassFlow)")
    println("Max power:$(PowerDiss)")
    
    return MassFlow,LinVel,PowerDiss
    end
    
    function Frict(Re,P,L,D,Dens,DynVisc)
        
        f = 64 / Re
        LinVel = sqrt(2*P/(L/D*f*Dens))
        
        for i = 1:10
            Re = Dens*D*LinVel/DynVisc
            f = 64 / Re
            LinVel = sqrt(2*P/(L/D*f*Dens))
        
        end

    return f, Re, LinVel

    end

    function Frict(ϵ,Re,P,L,D,Dens,DynVisc)
        # A = -2*log((ϵ/3.7)-5.02/Re*log(ϵ/3.7+13/Re)); %"A Review of Explicit Friction Factor Equations"  Zigrang and Sylvester eq. 13
        # FFac = A.^(-2);
        
        # FFacEq(F)  = 1.14-2*log10(ϵ+9.35/(Re.*sqrt(F)))-1/sqrt(F)
        FFacEq(F)  = -2*log10( (ϵ/3.7) - 5.02/Re * log10(ϵ/3.7 + 13/Re) ) - 1 /√(F)
        OptimCost(F) = abs.(FFacEq(F))
        f = optimize(OptimCost,0.,1.).minimizer
        LinVel = sqrt(2*P/(L/D*f*Dens))
        
        for i = 1:10
            Re = Dens*D*LinVel/DynVisc
            FFacEq(F)  = -2*log10( (ϵ/3.7) - 5.02/Re * log10(ϵ/3.7 + 13/Re) ) - 1 /√(F)
            OptimCost(F) = abs.(FFacEq(F))
            f = optimize(OptimCost,0.,1.).minimizer
            LinVel = sqrt(2*P/(L/D*f*Dens))
        
        end

        return f, Re, LinVel
        
        
    end