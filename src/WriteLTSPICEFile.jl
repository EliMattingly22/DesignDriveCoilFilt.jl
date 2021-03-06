function writeFilterSPICE(FileName,
    RDrive,
    NumDriveElements,
    LDrive,
    CDrive,
    SerCap,
    CParAct,
    LFilt,
    LFilt2_ESR,
    LFilt1_ESR,
    LFilt2,
    CFilt,
    CFilt2,
    LNotch,
    LNotch_ESR,
    CNotch,
    LNotch2,
    LNotch2_ESR,
    CNotch2,
    RDamp,
    LDrivePair,
    CDrivePair)


    open(FileName,"w") do file
        write(file,
"Version 4
SHEET 1 1972 680
WIRE 352 -80 256 -80
WIRE 512 -80 432 -80
WIRE 1184 -64 1072 -64
WIRE 1376 -64 1264 -64
WIRE -32 64 -128 64
WIRE 112 64 48 64
WIRE 224 64 176 64
WIRE 256 64 256 -80
WIRE 256 64 224 64
WIRE 288 64 256 64
WIRE 432 64 368 64
WIRE 512 64 512 -80
WIRE 512 64 496 64
WIRE 768 64 512 64
WIRE 912 64 768 64
WIRE 1072 64 1072 -64
WIRE 1072 64 912 64
WIRE 1104 64 1072 64
WIRE 1264 64 1184 64
WIRE 1376 64 1376 -64
WIRE 1376 64 1328 64
WIRE 1424 64 1376 64
WIRE 1600 64 1424 64
WIRE 1712 64 1664 64
WIRE 1872 64 1776 64
WIRE 1872 96 1872 64
WIRE -128 128 -128 64
WIRE 224 160 224 64
WIRE 224 160 144 160
WIRE 272 160 224 160
WIRE 768 160 768 64
WIRE 912 160 912 64
WIRE 144 192 144 160
WIRE 272 192 272 160
WIRE 1424 208 1424 64
WIRE -128 240 -128 208
WIRE 1872 240 1872 176
WIRE 768 272 768 240
WIRE 912 272 912 240
WIRE -128 384 -128 320
WIRE 144 384 144 272
WIRE 144 384 -128 384
WIRE 272 384 272 256
WIRE 272 384 144 384
WIRE 768 384 768 336
WIRE 768 384 272 384
WIRE 912 384 912 336
WIRE 912 384 768 384
WIRE 1424 384 1424 272
WIRE 1424 384 912 384
WIRE 1872 384 1872 320
WIRE 1872 384 1424 384
WIRE -128 432 -128 384
FLAG -128 432 0
SYMBOL voltage -128 224 R0
WINDOW 123 0 0 Left 0
WINDOW 39 0 0 Left 0
SYMATTR InstName V1
SYMATTR Value AC 1
SYMBOL ind 64 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LFilt1
SYMATTR Value $(num2SciString(LFilt))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt1_ESR))
SYMBOL cap 176 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CFilt1
SYMATTR Value $(num2SciString(CFilt))
SYMBOL res -144 112 R0
SYMATTR InstName RSrc
SYMATTR Value 1m
SYMBOL ind 384 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LFilt2
SYMATTR Value $(num2SciString(LFilt))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt1_ESR))
SYMBOL cap 496 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CFilt2
SYMATTR Value $(num2SciString(CFilt))
SYMBOL ind 128 176 R0
SYMATTR InstName LFilt3
SYMATTR Value $(num2SciString(LFilt2))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt2_ESR))
SYMBOL cap 256 192 R0
SYMATTR InstName CFilt3
SYMATTR Value $(num2SciString(CFilt2))
SYMBOL ind 752 144 R0
SYMATTR InstName LNotch
SYMATTR Value $(num2SciString(LNotch))
SYMATTR SpiceLine Rser=$(num2SciString(LNotch_ESR))
SYMBOL cap 752 272 R0
SYMATTR InstName CNotch
SYMATTR Value $(num2SciString(CNotch))
SYMBOL cap 1408 208 R0
SYMATTR InstName CPar
SYMATTR Value $(num2SciString(CParAct))
SYMBOL cap 1664 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CSer
SYMATTR Value $(num2SciString(SerCap))
SYMBOL cap 1776 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CDrive
SYMATTR Value $(num2SciString(CDrive/NumDriveElements))
SYMBOL ind 1856 80 R0
SYMATTR InstName LDrive
SYMATTR Value $(num2SciString(LDrive*NumDriveElements))
SYMBOL res 1856 224 R0
SYMATTR InstName RDrive
SYMATTR Value $(num2SciString(RDrive*NumDriveElements))
SYMBOL res 448 -96 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName RDamp1
SYMATTR Value $(num2SciString(RDamp))
SYMBOL ind 1200 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LDrivePair
SYMATTR Value $(num2SciString(LDrivePair))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt1_ESR))
SYMBOL cap 1328 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CDrivePair
SYMATTR Value $(num2SciString(CDrivePair))
SYMBOL res 1280 -80 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName RDamp2
SYMATTR Value $(num2SciString(RDamp))
SYMBOL ind 896 144 R0
SYMATTR InstName LNotch1
SYMATTR Value $(num2SciString(LNotch2))
SYMATTR SpiceLine Rser=$(num2SciString(LNotch2_ESR))
SYMBOL cap 896 272 R0
SYMATTR InstName CNotch1
SYMATTR Value $(num2SciString(CNotch2))
TEXT 264 432 Left 2 !.ac oct 1k 1k 500k")
    end
end


function num2SciString(Val::Real,sigdigits=4)
    Val = round(Val;sigdigits=sigdigits)
    
    if abs(Val)>1e3
        OutVal =  ((Val *1e-3))
         ReturnStr =  "$(OutVal)k"
    elseif abs(Val)>1
        OutVal =  ((Val ))
         ReturnStr =  "$OutVal"
    elseif abs(Val)>1e-3
        OutVal =  ((Val *1e3))
         ReturnStr =  "$(OutVal)m"
    elseif abs(Val)>1e-6
        OutVal =  ((Val *1e6))
         ReturnStr =  "$(OutVal)u"
    elseif abs(Val)>1e-9
        OutVal =  ((Val *1e9))
         ReturnStr =  "$(OutVal)n"
    elseif abs(Val)>1e-12
        OutVal =  ((Val *1e12))
         ReturnStr =  "$(OutVal)p"
    else
        OutVal =  ((Val *1e15))
         ReturnStr =  "$(OutVal)f"
    end

    if length(ReturnStr)>(sigdigits+2)
        ReturnStr = ReturnStr[1:sigdigits+1]*ReturnStr[end]
    end
    return ReturnStr
end

