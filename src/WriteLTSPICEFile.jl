function writeFilterSPICE(FileName,
    RDrive,
    NumDriveElements,
    LDrive,
    CDrive,
    SerCap,
    LTee_2,
    LTee_2_ESR,
    LTee_1,
    LTee_1_ESR,
    CParAct,
    LFilt,
    LFilt2_ESR,
    LFilt1_ESR,
    LFiltMatch_C,
    CFilt,
    CFiltMatch_L,
    LNotch,
    LNotch_ESR,
    CNotch,
    LNotch2,
    LNotch2_ESR,
    CNotch2,
    LNotch_Tune,
    LNotch_Tune_ESR,
    RDamp)


    open(FileName,"w") do file
        write(file,
"Version 4
SHEET 1 1972 680
WIRE 352 -80 256 -80
WIRE 512 -80 432 -80
WIRE 976 -64 864 -64
WIRE 1216 -64 1056 -64
WIRE -32 64 -128 64
WIRE 112 64 48 64
WIRE 224 64 176 64
WIRE 256 64 256 -80
WIRE 256 64 224 64
WIRE 288 64 256 64
WIRE 432 64 368 64
WIRE 512 64 512 -80
WIRE 512 64 496 64
WIRE 576 64 512 64
WIRE 704 64 656 64
WIRE 768 64 704 64
WIRE 864 64 864 -64
WIRE 864 64 768 64
WIRE 896 64 864 64
WIRE 1056 64 976 64
WIRE 1216 64 1216 -64
WIRE 1216 64 1120 64
WIRE 1296 64 1216 64
WIRE 1424 64 1376 64
WIRE 1472 64 1424 64
WIRE 1600 64 1552 64
WIRE 1712 64 1664 64
WIRE 1872 64 1776 64
WIRE 1872 96 1872 64
WIRE -128 128 -128 64
WIRE 224 160 224 64
WIRE 224 160 144 160
WIRE 272 160 224 160
WIRE 704 160 704 64
WIRE 768 160 768 64
WIRE 144 192 144 160
WIRE 272 192 272 160
WIRE 1424 208 1424 64
WIRE -128 240 -128 208
WIRE 1872 240 1872 176
WIRE 704 272 704 240
WIRE 768 272 768 240
WIRE -128 384 -128 320
WIRE 144 384 144 272
WIRE 144 384 -128 384
WIRE 272 384 272 256
WIRE 272 384 144 384
WIRE 704 384 704 336
WIRE 704 384 272 384
WIRE 768 384 768 336
WIRE 768 384 704 384
WIRE 1424 384 1424 272
WIRE 1424 384 768 384
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
SYMATTR InstName CLFilt_C1
SYMATTR Value $(num2SciString(LFiltMatch_C))
SYMBOL res -144 112 R0
SYMATTR InstName RSrc
SYMATTR Value 0.001
SYMBOL ind 384 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LFilt2
SYMATTR Value $(num2SciString(LFilt))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt2_ESR))
SYMBOL cap 496 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CLFilt_C2
SYMATTR Value $(num2SciString(LFiltMatch_C))
SYMBOL ind 128 176 R0
SYMATTR InstName LCFilt_L
SYMATTR Value $(num2SciString(CFiltMatch_L))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt2_ESR))
SYMBOL cap 256 192 R0
SYMATTR InstName CFilt
SYMATTR Value $(num2SciString(CFilt))
SYMBOL ind 752 144 R0
SYMATTR InstName LNotch
SYMATTR Value $(num2SciString(LNotch))
SYMBOL cap 752 272 R0
SYMATTR InstName CNotch
SYMATTR Value $(num2SciString(CNotch))
SYMBOL ind 672 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LNotch_Tune
SYMATTR Value $(num2SciString(LNotch_Tune))
SYMATTR SpiceLine Rser=$(num2SciString(LNotch_Tune_ESR))
SYMBOL ind 1392 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LTee1
SYMATTR Value $(num2SciString(LTee_1))
SYMATTR SpiceLine Rser=$(num2SciString(LTee_1_ESR))
SYMBOL ind 1568 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LTee2
SYMATTR Value $(num2SciString(LTee_2))
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
SYMBOL ind 992 48 R90
WINDOW 0 5 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName LDamp
SYMATTR Value $(num2SciString(LFilt))
SYMATTR SpiceLine Rser=$(num2SciString(LFilt1_ESR))
SYMBOL cap 1120 48 R90
WINDOW 0 0 32 VBottom 2
WINDOW 3 32 32 VTop 2
SYMATTR InstName CDamp
SYMATTR Value $(num2SciString(LFiltMatch_C))
SYMBOL res 1072 -80 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName RDamp2
SYMATTR Value $(num2SciString(RDamp))
SYMBOL ind 688 144 R0
WINDOW 0 -86 31 Left 2
SYMATTR InstName LNotch2
SYMATTR Value $(num2SciString(LNotch2))
SYMATTR SpiceLine Rser=$(num2SciString(LNotch2_ESR))
SYMBOL cap 688 272 R0
WINDOW 0 -94 12 Left 2
SYMATTR InstName CNotch2
SYMATTR Value $(num2SciString(CNotch2))
TEXT 264 432 Left 2 !.ac oct 1k 1k 500k
")
    end
end


function num2SciString(Val::Real,sigdigits=3)
    Val = round(Val;sigdigits=sigdigits)
    
    if Val>1e3
        OutVal =  ((Val *1e-3))
         ReturnStr =  "$(OutVal)k"
    elseif Val>1
        OutVal =  ((Val ))
         ReturnStr =  "$OutVal"
    elseif Val>1e-3
        OutVal =  ((Val *1e3))
         ReturnStr =  "$(OutVal)m"
    elseif Val>1e-6
        OutVal =  ((Val *1e6))
         ReturnStr =  "$(OutVal)u"
    elseif Val>1e-9
        OutVal =  ((Val *1e9))
         ReturnStr =  "$(OutVal)n"
    elseif Val>1e-12
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
