using JAC
using JLD2, Printf
using ExtrapolateDR

data = load("C:\\Users\\Young.Y\\Documents\\Julia_old\\2p_n10-20 DR.jld2")
DR = data["wb"]["dielectronic recombination pathways:"]

DRSetting   = SetDR([[Subshell(2,1),10,20,115],[Subshell(2,-2),10,20,115]], [["typeI",Subshell(2,-1)],["typeII",[Subshell(2,1),Subshell(2,-2)]]])
Origin      = OriginalDR(DR, DRSetting)
FitResult   = FitOriginDR(Origin, 15)
Extrapolation = Extrapolate(DRSetting, FitResult)

resultE = Vector{Float64}()
resultDR = Vector{Float64}()
resultWidth = Vector{Float64}()

for (key, value) in Extrapolation
    for i = eachindex(value)
        push!(resultE, Extrapolation[key][i].Energy.E)
        push!(resultDR, Extrapolation[key][i].DR.strength)
        push!(resultWidth, 0.0)
    end
end

sortindex = sortperm(resultE)
resultE = resultE[sortindex]
resultDR = resultDR[sortindex]
resultWidth = resultWidth[sortindex]

outputPath = "C:/Users/Young.Y/Documents/Julia/ExtrapolateDR/DR20.txt"
open(outputPath, "w") do DRresult
    for i = eachindex(resultE)
        str = @sprintf("%4.4f    %1.20e    %1.5e", resultE[i], resultDR[i], resultWidth[i])
        write(DRresult, str, "\n")
    end
end


