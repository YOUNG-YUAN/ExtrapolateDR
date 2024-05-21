module ExtrapolateDR

    using JAC, Optim

    export  Int64, QuanmtumDefects, resonanceEnergy, LevelEnergy, AutoionizationRate, RadType, 
            TransitionRate, DCstrength, branchRate, DRstrength, DRdata

    export  readLevel, SetDR, EmptyOriginalDR, OriginalDR, EmptyExtrapolation, SolveLmax, SolveQuantumDefects, 
            FitLevelEnergy, FitAIRate, FitRadRate, FitOriginDR, Extrapolate, ExtrapolateAll


    """
    `mutable struct RadType`

        ... defines de-excitation radiative transition type of DR.
        + type          ::String    ...which type of this radiative transition
        + KeySubshell   ::Subshell  ...characteristic subshell of this radiative transition type   
        + Rate          ::Float64   ...gross transition rate of this radiative transition type
    """
    mutable struct RadType
        type            ::String
        keySubshell     ::Union{Subshell, Vector{Subshell}}
        rate            ::Float64
    end

    
    """
    `mutable struct RadRate`
    
        ... defines radiative transition rate of DR which including all de-excitation radiative transition types.

        + totalRate ::Float64               ...gross radiative transtion rate of all radiative transition types    
        + component ::Dict{String, RadType} ...which type of this radiative transition
    """
    mutable struct RadRate
        totalRate   ::Float64
        component   ::Dict{String, RadType}
    end

    """
    `mutable struct DRdata`

        ... defines the result of extrapolation of DR.

        + N                     ::Int64     ...principal quantum number of outest electron
        + Nloss                 ::Int64     ...quantum defects
        + i2J                   ::Int64     ...2J of initial level
        + Energy                ::Float64   ...resonance energy of intermediate state of DR
        + iLevelEnergy          ::Float64   ...resonance energy of initial state of DR
        + mLevelEnergy          ::Float64   ...resonance energy of intermediate(double-excited) state of DR
        + AutoionizationRate    ::Float64   ...autoionization rate of DR.
        + radRate               ::RadRate   ...radiative rate of DR which including all de-excitation transition types.
        + DCStrength            ::Float64   ...dielectronic capture strength of DR
        + branchRate            ::Float64   ...branchrate of DR
        + DRStrength            ::Float64   ...dielectronic recombination strength
    """
    mutable struct DRdata
        N                   ::Int64
        Nloss               ::Int64
        i2J                 ::Int64    
        iLevelEnergy        ::Float64
        mLevelEnergy        ::Float64
        resonanceEnergy     ::Float64
        autoionizationRate  ::Float64
        radRate             ::RadRate
        DCStrength          ::Float64
        branchRate          ::Float64
        DRStrength          ::Float64
    end

    
    """
    `ExtrapolateDR.readLevel(level::Level)`

        ... get the max principal quantum number N of the outest electron.
            
            Return N::Int64
    """
    function readLevel(level::Level)
        CsfsIndex   = argmax(abs.(level.mc))
        Occupation  = level.basis.csfs[CsfsIndex].occupation
        N_index     = maximum(findall(!iszero,Occupation))
        N           = level.basis.subshells[N_index].n
        return N
    end


    """
    `ExtrapolateDR.readLevel(level::Level,TargetSubshell::Subshell)`

        ... determine whether target subshell exists in a specified level.
            
            Return isExist::Bool
    """
    function readLevel(TargetSubshell::Subshell, level::Level)
        CsfsIndex    = argmax(abs.(level.mc))
        Occupation   = level.basis.csfs[CsfsIndex].occupation
        allSubshells = level.basis.subshells[findall(!iszero,Occupation)]
        isExist      = false
        for i = eachindex(allSubshells)
            if allSubshells[i].n != TargetSubshell.n
                continue
            elseif allSubshells[i].kappa != TargetSubshell.kappa
                continue
            else
                isExist = true
                break
            end
        end
        return isExist
    end


    """
    `ExtrapolateDR.readLevel(subshellVec::Vector{Subshell}, level::Level, Option::Symbol=:any)`

        ... Option::Symbol=:any (Default) determine whether at least one subshell of subshellVec belongs to a specified level.
            Option::Symbol=:all determine whether all subshells of subshellVec belong to a specified level.

            Return isExist::Bool
    """
    function readLevel(subshellVec::Vector{Subshell}, level::Level, Option::Symbol=:any)
        if Option == :any   # if any subshell of subshellVec belongs to this level, then return a true, otherwise return a false.
            isExist = false
            for i = eachindex(subshellVec)
                if readLevel(subshellVec[i], level)
                    isExist = true
                    break
                else
                    continue
                end
            end
            return isExist
        elseif Option == :all  # only if all subshells of subshellVec belong to this level, then return a true, otherwise return a false.
            isExist = true
            for i = eachindex(subshellVec)
                if readLevel(subshellVec[i], level)
                    continue
                else
                    isExist = false
                    break
                end
            end
            return isExist
        end
    end


    """
    `ExtrapolateDR.SetDR(SetCaptureType::Any, SetRadType::Any)`

        ... Set all key characteristics of the DR process. For example:
        
            SetCaptureType = [[Subshell(2,1), 10, 20, 115], [Subshell(2,-2), 12, 22, 120)]]

            Subshell(2,1) is 2p_1/2, which means inner-excited subshell of double-excited state is 2p_1/2,
            and 10 means double-excited state which includes 2p_1/2 will allow to autoionize from N = 10,
            and 20 means N of double-excited state which includes 2p_1/2 is up to 20,
            and 115 means your want extrapolate the DR process to N = 115 for those double-excited state which includes 2p_1/2;
            In short, from your calculated DR data, your N = 10~20 for 2p_1/2, and you want to extrapolate 2p_1/2 to N = 115.
            [Subshell(2,-2), 10, 20, 115] is for 2p_3/2, similar to [Subshell(2,-2), 10, 20, 115], which means you calculated N = 12~22 for 2p_3/2, and you want to extrapolate 2p_3/2 to N = 120.

            SetRadType = [["typeI", Subshell(2,-1)], ["typeII"], [Subshell(2,1), Subshell(2,-2)]]

            If a final state of DR includes Subshell(2,-1), then the de-excitation transition from double-excited state to this final state is "typeI".
            Similarly, if a final state of DR includes one of Subshell from [Subshell(2,1), Subshell(2,-2)], then it's "typeII".

            Return DRSetting::Dict
    """
    function SetDR(SetCaptureType::Vector, SetRadType::Vector)
        CaptureSetting      = Dict{Subshell, Dict}()
        for i = eachindex(SetCaptureType)
            # [[Subshell,10,11,115]]
            keySubshell     = SetCaptureType[i][1]
            startN          = SetCaptureType[i][2]
            endN            = SetCaptureType[i][3]
            extrapolateN    = SetCaptureType[i][4]
            CaptureSetting[keySubshell] = Dict("KeySubshell" => keySubshell, "startN" => startN, "endN" => endN, "extrapolateN" => extrapolateN)
        end
        RadSetting   = Dict{String, Dict}()
        for i = eachindex(SetRadType)
            # [["type1",Subshell]]
            type            = SetRadType[i][1]
            keySubshell     = SetRadType[i][2]
            RadSetting[type] = Dict("type" => type, "KeySubshell" => keySubshell)
        end
        DRSetting = Dict("Capture" => CaptureSetting, "Transition" => RadSetting)
        return DRSetting
    end


    """
    `ExtrapolateDR.EmptyOriginalDR(DRSetting::Dict)`

        ... create a empty Dict for storing neccessary data, based on DRSetting, from calculated DR data.
            So you need to use function SetDR(SetCaptureType::Vector, SetRadType::Vector) before using this.
            But you do not need to use this function separately, it will be used in function OriginalDR(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict)
                
            Return Origin::Dict
    """
    function EmptyOriginalDR(DRSetting::Dict)
        captureSetting  = DRSetting["Capture"]
        radSetting      = DRSetting["Transition"]
        Origin = Dict{Subshell, Vector{DRdata}}()  # create an empty Dict for saving original DR data.
        for (key, value) in captureSetting  # CaptureSetting includes all capture types
            emptyDRdataVector   = Vector{DRdata}()
            startN              = value["startN"]
            endN                = value["endN"]
            for i = startN : endN
                N                   = i
                Nloss               = 0.
                i2J                 = 0
                iLevelEnergy        = 0.
                mLevelEnergy        = 0.
                resonanceEnergy     = 0.
                autoionizationRate  = 0.
                    component           = Dict{String,RadType}()
                    for (key, value)    = RadSetting        # Each key corresponds a type.
                        type            = value["type"]
                        keySubshell     = value["KeySubshell"]
                        rate            = 0.
                        radType         = RadType(type, keySubshell, Rate)
                        component[type] = RadType
                    end
                    totRate           = 0.
                radRate             = RadRate(totRate, component)
                DCStrength          = 0.
                branchRate          = 0.
                DRStrength          = 0.
                push!(emptyDRdataVector, DRdata(N, Nloss, i2J, iLevelEnergy, mLevelEnergy, resonanceEnergy, autoionizationRate, radRate, DCStrength, branchRate, DRStrength))
            end
            Origin[key]         = emptyDRdataVector
        end
        return Origin
    end
    

    """
    `ExtrapolateDR.OriginalDR(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict)`

        ... get all neccessary data for extrapolation, based on DRSetting, from already calculated DR data Vector{JAC.Dielectronic.Pathway}.
            So you need to calculate a DR process and to use function SetDR(SetCaptureType::Vector, SetRadType::Vector) before using this.

            Return Origin::Dict
    """
    function OriginalDR(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict)      
        Origin = EmptyOriginalDR(DRSetting) # create an empty Dict for saving original DR data.
        for i = eachindex(DR)
            if DR[i].photonRate.Coulomb == 0.0  # if photonRate = 0.0, that means this JAC.Dielectronic.Pathway is invalid.
                continue
            end
            # Transfer Units
            iLevel          = DR[i].initialLevel
            mLevel          = DR[i].intermediateLevel
            fLevel          = DR[i].finalLevel
            iLevelEnergy    = Defaults.convertUnits("energy: from atomic to eV", iLevel.resonanceEnergy)
            mLevelEnergy    = Defaults.convertUnits("energy: from atomic to eV", mLevel.resonanceEnergy)
            Ed              = Defaults.convertUnits("energy: from atomic to eV", DR[i].electronEnergy)
            captureRate     = Defaults.convertUnits("rate: from atomic", DR[i].captureRate) 
            photonRate      = Defaults.convertUnits("rate: from atomic", DR[i].photonRate.Coulomb)
            
            Ni              = readLevel(mLevel)       
            for (key,_) in Origin
                if readLevel(key, mLevel)  # determine whether KeySubshell exists in this double-excited state
                    index = findfirst(x -> x.N == Ni, Origin[key])  # find the index of Origin[key] corresponding Ni
                    Origin[key][index].i2J                      = iLevel.J.den == 1 ? 2 * iLevel.J.num : iLevel.J.num # 2j of initial state           
                    Origin[key][index].iLevelEnergy             = iLevelEnergy  # constant value
                    if Origin[key][index].mLevelEnergy    == 0.
                        Origin[key][index].mLevelEnergy         = mLevelEnergy
                    else
                        Origin[key][index].mLevelEnergy         = max(mLevelEnergy, Origin[key][index].mLevelEnergy)  # the max mLevelEnergy would be saved, because which is low-quanmtum-defects.
                    end
                    Origin[key][index].resonanceEnergy          = Ed
                    Origin[key][index].autoionizationRate  = Origin[key][index].autoionizationRate + captureRate
                        for (key2,value2) in Origin[key][index].radRate.component
                            if readLevel(value2.KeySubshell, fLevel)     # determine which transition type the Pathway corresponds
                                Origin[key][index].radRate.component[key2].Rate  = Origin[key][index].radRate.component[key2].Rate + photonRate
                                Origin[key][index].radRate.totalRate             = Origin[key][index].radRate.totalRate + photonRate
                            else
                                continue
                            end
                        end
                    Origin[key][index].DCStrength = Origin[key][index].DCStrength + 4.95e-30 * (mLevel.J.num + 1) / (2 * iLevel.J.num + 2) / Ed * captureRate
                    # Branch and DR would not change, keep them zeros, because we don't need them.
                else
                    continue
                end
            end
        end
        return Origin
    end


    """
    `ExtrapolateDR.EmptyOriginalDR(DRSetting::Dict)`

        ... create a empty Dict for storing neccessary data, based on DRSetting, from calculated DR data.
            So you need to use function SetDR(SetCaptureType::Vector, SetRadType::Vector) before using this.
            But you do not need to use this function separately, it will be used in function OriginalDR(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict)
                
            Return Origin::Dict
    """
    function EmptyExtrapolation(DRSetting::Dict)
        captureSetting  = DRSetting["Capture"]
        radSetting      = DRSetting["Transition"]
        Extrapolation = Dict{Subshell, Vector{DRdata}}()  # create an empty Dict for saving original DR data.
        for (key, value) in captureSetting  # CaptureSetting includes all capture types
            emptyDRdataVector   = Vector{DRdata}()
            startN              = value["endN"] + 1     # calculated from "startN" to "endN", so extrapolation is from "endN"+1 to "extrapolateN".
            endN                = value["extrapolateN"]
            for i = startN : endN
                N                   = i
                Nloss               = 0.
                i2J                 = 0               
                iLevelEnergy        = 0.
                mLevelEnergy        = 0.
                resonanceEnergy     = 0.
                autoionizationRate  = 0.
                    component           = Dict{String,RadType}()
                    for (key, value)    = radSetting        # Each key corresponds a type.
                        type            = value["type"]
                        keySubshell     = value["KeySubshell"]
                        rate            = 0.
                        radType  = RadType(type, keySubshell, rate)
                        component[type] = radType
                    end
                    totalRate           = 0.
                radRate             = RadRate(totalRate, component)
                DCStrength          = 0.
                branchRate          = 0.
                DRStrength          = 0.
                push!(emptyDRdataVector, DRdata(N, Nloss, i2J, iLevelEnergy, mLevelEnergy, resonanceEnergy, autoionizationRate, transitionRate, DCStrength, branchRate, DRStrength))
            end
            Extrapolation[key]         = emptyDRdataVector
        end
        return Extrapolation
    end


    """
    `ExtrapolateDR.SolveQuantumLoss()`
        ... solve quantum loss of all levels
            But haven't finished!!!
    """
    function SolveQuantumDefects()
        return 0.0
    end

    
    """
    `ExtrapolateDR.FitLevelEnergy(NVec::Vector{Int64}, mLevelEnergyVec::Vector{Float64}, Z::Int64, guess_EI::Float64=0.)`
    
        ... fit level energy with function: 
            " En = EI - Ryd * Z^2 / n^2 ". 
            Least square method is used.
                       NVec : Vector of N. 
            mLevelEnergyVec : Vector of level energy of double-excited state,
                              N is corresponding to NVec.
                   guess_EI : Vector of fitting initial value, which includes one 
                              element [EI], EI derives from 
                              "En = EI - Ryd * Z^2 / n^2", but you don't need to 
                              customize it.

            return (EI::Float64, fitFuncString::String)
    """
    function FitLevelEnergy(NVec::Vector{Int64}, mLevelEnergyVec::Vector{Float64}, Z::Int64, guess_EI::Vector{Float64}=[0.])
        x = NVec
        y = mLevelEnergyVec
        function RSS(p)   # functon of sum of squares of residuals
            EI = p[1]   # a is fitting parameter
            return sum((y .- (EI .- 13.6057 * Z^2 ./ x.^2)).^2)
        end
        result = optimize(RSS, guess_EI) # least square method
        EI = Optim.minimizer(result)
        fitFuncString = "En = EI - Ryd * Z^2 / n^2"
        return (EI, fitFuncString)
    end


    """
    `ExtrapolateDR.FitAIRate(NVec::Vector{Int64}, AIRateVec::Vector{Float64}, guess_a::Float64=0.)`

        ... fit autoionization rate with function:
            " Aa = a / n^3 ".
            Least square method is used.
                 NVec : Vector of N. 
            AIRateVec : Vector of a sum autoionization rate of some one type,
                        with same N, N is corresponding to NVec.
              guess_a : Vector of fitting initial value, which includes one 
                        element [a], a derives from "Aa = a / n^3" , but you 
                        don't need to customize it.

            return (a::Float64 , fitFuncString::String)
    """
    function FitAIRate(NVec::Vector{Int64}, AIRateVec::Vector{Float64}, guess_a::Vector{Float64}=[0.])
        x = NVec
        y = AIRateVec
        ymax = maximum(y)
        y_  = y / ymax * 1e6  # This transformance is for making this fitting suitable for those large order of magnitude rates.
        function RSS(p)   # functon of sum of squares of residuals
            a = p[1]   # a is fitting parameter
            return sum((y_ .- (a ./ x.^3)).^2)
        end
        result = optimize(RSS, guess_a / ymax * 1e6) # least square method, guessed parameter should be transformed too.
        a = result.minimizer[1]
        a = a * ymax * 1e-6  # transform back
        fitFuncString = "Aa = a / n^3"
        return a, fitFuncString
    end
    

    """
    `SolveLmax(NVec::Vector{Int64}, DCstrengthVec::Vector{Float64}, fitAIRateVec::Vector{Float64}, EdVec::Vector{Float64}, ji:: Int64, jc::Float64)`

        ... solve averaged max orbital angular momentum of double-excited state, which is solved by this function:
                " S_DC = gd / (2gi) * (π^2 * ħ^3 / (me * Eij)) * Aa ",
                "   gd = 2 * (2jc + 1)(Lmax + 1)^2                  ",
            where jc is the inner excited electron's j, and 2jc + 1 is the statistic weight.
            According to π^2 * ħ^3 / me = 4.95e-30 eV^2·cm^-2, So this function can be rewrited by :
                " (Lmax + 1)^2 = gi * Eij * S_DC / (4.95e-30 * (2jc)) ".
            where S_DC derives from sum of all original data's DC strength with same N, Aa is autoionization rate, derives from extrapolation, with corresponding N.
                     NVec : Vector of N. 
            DCstrengthVec : Vector of S_DC, each S_DC is a sum DCstrength 
                            with same N, N is corresponding to NVec.
             fitAIRateVec : Vector of extrapolated autoionization rate, 
                            N is corresponding to NVec.
                    EdVec : Vector of Eij.
                       ji : 2j of initial state of DR
                       jc : 2j of inner excited electron

            return Lmax::Float64
    """
    function SolveLmax(NVec::Vector{Int64}, DCstrengthVec::Vector{Float64}, fitAIRateVec::Vector{Float64}, EdVec::Vector{Float64}, ji:: Int64, jc::Int64)
        gi      = ji + 1  # ji is 2j of initial state of DR, and gi is statistic weight 2j + 1 of initial state.    
        LmaxVec = sqrt.(gi * EdVec .* DCstrengthVec .* NVec.^3 / 4.95e-30 / (jc + 1) / fitAIRateVec) .- 1
        Lmax    = sum(LmaxVec) / length(LmaxVec)
        return Lmax
    end
    

    """
    `ExtrapolateDR.FitRadRate(NVec::Vector{Int64}, RadRateVec::Vector{Float64}, radType::Int64, guess_a::Vector{Float64}=[0., 0.])`

        ... fit radiative transition rate. There are two excited electrons of Double-excited state, 
            so the radiative transition is of different types. Some types can be fitted by function
                " Ar = a / n^3 ",
            and the others can be fitted by function
                " Ar = b ".
            This function will pick a suitable function to fit transition rate of some one type. 
            Least square method is used.
                    NVec : Vector of N. 
            radRateVec : Vector of a sum transition rate of some one type,
                           with same N, N is corresponding to NVec.
                 guess_a : Vector of fitting initial value, which includes two 
                           element [a, b], a derives from "Ar = a / n^3" and 
                           b derives from "Ar = b", but you don't need to 
                           customize it.
            
            return (a::Float64 , fitFuncString::String)
                   OR (b::Float64 , fitFuncString::String)               
    """
    function FitRadRate(NVec::Vector{Int64}, radRateVec::Vector{Float64}, radType::Int64, guess_a::Vector{Float64}=[1., 1.])
        x = NVec
        y = radRateVec
        ymax = maximum(y)
        y_  = y / ymax * 1e4  # This transformance is for making this fitting suitable for those large order of magnitude rates.
        if radType == 1  # typeII, then fit by "Ar = b"
            function RSS(p)
                b = p[1]
                return sum((y_ .- (b)).^2)
            end
            result = optimize(RSS, [guess_a[2] * ymax * 1e4])
            b = result.minimizer[1]
            b = b * ymax * 1e-4  # transform back
            fitFuncString = "Ar = b"
            return b, fitFuncString     # return of function FitRadRate()
        elseif radType == 2  # typeII, then fit by "Ar = a / n^3"
            function RSS2(p)
                a = p[1]
                return sum((y_ .- (a ./ x.^3)).^2)
            end
            result = optimize(RSS2, [guess_a[1] * ymax * 1e4])
            a = result.minimizer[1]
            a = a * ymax * 1e-4  # transform back
            fitFuncString = "Ar = a / n^3"
            return a, fitFuncString     # return of function FitRadRate()
        else
            error("Wrong radType! Available for 1 and 2!")
        end
    end


    """
    `ExtrapolateDR.FitRadRate(Origin::Dict, Z::Int64)`

        ... fit original DR date for exptrapolation, including level energy of double-excited state, 
            autoionization rate, averaged max orbital angular momentum (Lmax), transition rate. Additionaly,
            it saved Z and 2j of initial state of DR.
            Origin : Dict derives from OriginalDR(...)
                 Z : Charges of ion before recombination.

            return FitResult::Dict
    """
    function FitOriginDR(Origin::Dict, Z::Int64)
        FitResult = Dict{Subshell, Dict}()
        for (key1, value1) in Origin  # read one capture type, then fit all one time, read one another capture type, then fit all one another time.
            NVec            = Vector{Int64}()
            mLevelEnergyVec = Vector{Float64}()
            AIRateVec       = Vector{Float64}()
            TransDict       = Dict{String,Vector{Float64}}()
            DCstrengthVec   = Vector{Float64}()
            LmaxVec         = Vector{Float64}()
            for i = eachindex(value1)                
                push!(NVec            , value1[i].N.num)
                push!(mLevelEnergyVec , value1[i].mLevelEnergy)
                push!(AIRateVec       , value1[i].autoionizationRate)
                push!(DCstrengthVec   , value1[i].DCStrength)
                push!(LmaxVec         , 0.)
                for (key2, value2) in value1[i].radRate.component
                    if haskey(TransDict, key2)
                        push!(TransDict[key2], value2.Rate)
                    else
                        TransDict[key2] = Vector{Float64}()
                        push!(TransDict[key2], value2.Rate)
                    end
                end
            end
            # Now we got all data for fitting, of one capture type.                    
            FitResult[key1]                     = Dict{String, Any}()
            FitResult[key1]["Z"]                = Z
            FitResult[key1]["i2J"]              = Origin[key1][1].i2J    # 2j + 1 of initial state
            # Fit energy level
            FitResult[key1]["iLevelEnergy"]     = Origin[key1][1].iLevelEnergy
            FitResult[key1]["EI"]               = Dict{String, Any}()
            (fitCoeff, fitFunc)                 = FitLevelEnergy(NVec, mLevelEnergyVec, Z)
            FitResult[key1]["EI"]["fitCoef"]    = fitCoeff[1]
            FitResult[key1]["EI"]["fitFunc"]    = fitFunc
            # Fit autoionization rate
            FitResult[key1]["AI"]               = Dict{String, Any}()
            (fitCoeff, fitFunc)                 = FitAIRate(NVec, AIRateVec)
            FitResult[key1]["AI"]["fitCoef"]    = fitCoeff
            FitResult[key1]["AI"]["fitFunc"]    = fitFunc
            # Solve Lmax
            #=
            angular momentum quantum number can be obtained by equation:
                kappa = -(l + 1)  =>  l = -kappa - 1   (if j = l + 1/2),
                kappa = l         =>  l = kappa        (if j = l - 1/2).
            jc is 2j of inner excited electron of double-excited state, which can be solve by:
                jc = 2j = 2l + 1 = -2 * kappa - 1 (if kappa < 0)
                jc = 2j = 2l - 1 = kappa - 1 (if kappa > 0)
            =#
            jc                                          = key1.kappa > 0 ? 2 * key1.kappa - 1 : -2 * key1.kappa - 1
            ji                                          = FitResult[key1]["i2J"]      # 2j of initial state
            EdVec                                       = mLevelEnergyVec .- value1[1].iLevelEnergy
            fitAIRateVec                                = FitResult[key1]["AI"]["fitCoef"] ./ NVec
            Lmax                                        = SolveLmax(NVec, DCstrengthVec, fitAIRateVec, EdVec, ji, jc)
            FitResult[key1]["Lmax"]                     = Lmax
            # Fit radiative transition rate
            FitResult[key1]["TR"]                       = Dict{String, Dict}()
            for (key2, _) in TransDict
                FitResult[key1]["TR"][key2]             = Dict{String, Any}()
                if key2 == "typeI"
                    (fitCoeff, fitFunc)                 = FitRadRate(NVec, TransDict[key2], Int64(1))
                elseif key2 == "typeII"
                    (fitCoeff, fitFunc)                 = FitRadRate(NVec, TransDict[key2], Int64(2))
                end
                FitResult[key1]["TR"][key2]["fitCoef"]  = fitCoeff
                FitResult[key1]["TR"][key2]["fitFunc"]  = fitFunc
            end
        end
        return FitResult
    end


    """
    `ExtrapolateDR.Extrapolate(DRSetting::Dict, FitResult::Dict)`

        ... extrapolate the DR process, including level energy of double-excited state, resonance energy of DR,
            autoionization rate, transition rate, DC strength, veraged max orbital angular momentum (Lmax), 
            branch rate of DR, DR strength.
            DRSetting : Dict derives from SetDR(...).
            FitResult : Dict derives from FitOriginDR(...).

            return FitResult::Dict
    """
    function Extrapolate(DRSetting::Dict, FitResult::Dict)
        Extrapolation = EmptyExtrapolation(DRSetting)
        for (key1, _) in Extrapolation
            for i = eachindex(Extrapolation[key1])
                Z                                           = FitResult[key1]["Z"]
                Ni                                          = Extrapolation[key1][i].N.num
                Extrapolation[key1][i].i2J                  = FitResult[key1]["i2J"]    # 2j of initial state
                EI                                          = FitResult[key1]["EI"]["fitCoef"]
                Extrapolation[key1][i].iLevelEnergy         = FitResult[key1]["iLevelEnergy"]
                Extrapolation[key1][i].mLevelEnergy         = EI - 13.6057 * Z^2 ./ Ni.^2
                Extrapolation[key1][i].resonanceEnergy      = Extrapolation[key1][i].mLevelEnergy - Extrapolation[key1][i].iLevelEnergy
                Extrapolation[key1][i].autoionizationRate   = FitResult[key1]["AI"]["fitCoef"] / Ni^3
                for (key2, _) in Extrapolation[key1][i].radRate.component                       
                    if FitResult[key1]["TR"][key2]["fitFunc"] == "Ar = a / n^3"  # fitting function is "Ar = a / n^3"
                        Rate = FitResult[key1]["TR"][key2]["fitCoef"] / Ni^3    
                    else  # fitting function is "Ar = b"
                        Rate = FitResult[key1]["TR"][key2]["fitCoef"]
                    end
                    Extrapolation[key1][i].radRate.component[key2].Rate  = Rate
                    Extrapolation[key1][i].radRate.totalRate             = Extrapolation[key1][i].radRate.totalRate + Rate
                end
                #=
                angular momentum quantum number can be obtained by equation:
                    kappa = -(l + 1)  =>  l = -kappa - 1   (if j = l + 1/2),
                    kappa = l         =>  l = kappa        (if j = l - 1/2).
                jc is 2j of inner excited electron of double-excited state, which can be solve by:
                    jc = 2j = 2l + 1 = -2 * kappa - 1 (if kappa < 0)
                    jc = 2j = 2l - 1 = kappa - 1 (if kappa > 0)
                =#
                jc      = key1.kappa > 0 ? 2 * key1.kappa - 1 : -2 * key1.kappa - 1
                Lmax    = FitResult[key1]["Lmax"]
                gi      = Extrapolation[key1][i].i2J + 1
                Ed      = Extrapolation[key1][i].resonanceEnergy
                Aa      = Extrapolation[key1][i].autoionizationRate
                Ar      = Extrapolation[key1][i].radRate.totalRate
                Extrapolation[key1][i].DCStrength   = 4.95e-30 * (2 * jc + 1) * (Lmax + 1)^2 * Aa / gi / Ed
                Extrapolation[key1][i].branchRate   = Ar / (Aa + Ar)
                Extrapolation[key1][i].DRStrength   = Extrapolation[key1][i].DCStrength * Extrapolation[key1][i].branchRate
            end
        end
        return Extrapolation    
    end

    
    """
    `ExtrapolateDR.Extrapolate(Origin::Dict, DRSetting::Dict, Z::Int64)`

        ...    Origin : Dict derives from SetDR(...).
            FitResult : Dict derives from FitOriginDR(...).
                    Z : Charges of ion before recombination.

            return FitResult::Dict
    """
    function Extrapolate(Origin::Dict, DRSetting::Dict, Z::Int64)
        FitResult       = FitOriginDR(Origin, Z)      
        Extrapolation   = Extrapolate(DRSetting, FitResult, Z)       
        return  Extrapolation
    end


    """
    `ExtrapolateAll(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict, Z::Int64)`

        ... do everything including getting original date from Vector{JAC.Dielectronic.Pathway}, and exptrapolate DR process.
                   DR : Vector{JAC.Dielectronic.Pathway} derives from JAC.Dielectronic.
            DRSetting : Dict derives from SetDR(...).
                    Z : Charges of ion before recombination.

            return Everything::Dict
    """
    function ExtrapolateAll(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict, Z::Int64)       
        Origin          = OriginalDR(DR::Vector{JAC.Dielectronic.Pathway}, DRSetting::Dict)
        FitResult       = FitOriginDR(Origin, Z)
        Extrapolation   = Extrapolate(DRSetting, FitResult)
        ALL             = Origin
        for (key, _) in ALL
            ALL[key] = vcat(ALL[key], Extrapolation[key])
        end
        Everything      = Dict("Origin" => Origin, "Extrapolation" => Extrapolation, "ALL" => ALL, "FitResult" => FitResult)
        return Everything
    end
    # End of module
end