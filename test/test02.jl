using JAC
using JLD2, Plots, Printf  
using ExtrapolateDR

data = load("C:\\Users\\Young.Y\\Documents\\Julia\\ExtrapolateDR\\test\\2p_n10-11_all.jld2")
DR = data["wb"]["dielectronic recombination pathways:"]

DRSetting   = SetDR([[Subshell(2,1),10,11,115],[Subshell(2,-2),10,11,115]], [["typeI",Subshell(2,-1)],["typeII",[Subshell(2,1),Subshell(2,-2)]]])
Everything  = ExtrapolateAll(DR, DRSetting, 15)
EVec = Vector{Float64}()
DRVec = Vector{Float64}()
for (key, value) in Everything["ALL"]   
    for i = eachindex(value)
        push!(EVec, value[i].Energy.E)
        push!(DRVec, value[i].DR.strength)
    end
end
f = open("DR.txt", "w")
for (elem1, elem2) in zip(EVec, DRVec)
    line = @sprintf("%+20.10e     %+20.10e     %2d\n",elem1, elem2, 0)
    write(f, line)
end
close(f)