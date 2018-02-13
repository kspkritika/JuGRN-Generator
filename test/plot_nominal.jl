# script to plot the nominal simulation values -
using PyPlot
include("DataDictionary.jl")

# load the data -
data_archive = readdlm("./tmp/simulation_nominal.dat")

# grab the time and data archive -
time_archive = data_archive[:,1]
state_archive = data_archive[:,2:end]
TD = 1.0 # hr
TI = 4.0 # hr

# load the data_dictionary -
data_dictionary = DataDictionary(0.0,0.0,0.0)
V = data_dictionary["volume_of_single_cell"]
M = data_dictionary["mass_of_single_cell"]
SF = (M/V) # convert to concentration (mumol/L)

# make the plot -
plot((time_archive - TD)*(60),state_archive[:,7]*(SF),linewidth=2.0,color="#000000")
plot((time_archive - TD)*(60),state_archive[:,8]*(SF),linewidth=2.0,color="#707070")
plot((time_archive - TD)*(60),state_archive[:,9]*(SF),linewidth=2.0,color="#A9A9A9")

# add some dots -
length_timepoints = length(time_archive)
sample_index = collect(1:50:length_timepoints)
time_translated = (time_archive[sample_index] - TD)*(60)
plot(time_translated,state_archive[sample_index,7]*(SF),"o",mfc="#FFFFFF",mec="#000000")
plot(time_translated,state_archive[sample_index,8]*(SF),"o",mfc="#FFFFFF",mec="#707070")
plot(time_translated,state_archive[sample_index,9]*(SF),"o",mfc="#FFFFFF",mec="#A9A9A9")

# add the inducer -
inducer_level = 10.0
time_translated = (time_archive - TD)*(60)
inducer_array = Float64[]
for time_value in time_translated

    if (time_value < 0.0)
        push!(inducer_array,0.0)
    elseif (time_value>= 0.0 && time_value <= (TI - TD)*60)
        push!(inducer_array,inducer_level)
    else
        push!(inducer_array,0.0)
    end
end
plot(time_translated,inducer_array,linewidth=2.0,color="#FF0000")

# write the axis -
xlabel("Time (min)",fontsize=18)
ylabel(L"Concentration ($\mu$M)",fontsize=18)

# save to disk -
savefig("./figs/Memory-Nominal.pdf")
