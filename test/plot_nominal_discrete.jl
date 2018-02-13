# script to plot the nominal simulation values -
using PyPlot
include("DataDictionary.jl")

# load the data -
data_archive = readdlm("./tmp/simulation_nominal_discrete.dat")

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
plot((time_archive - TD)*(60),state_archive[:,4]*(SF),"--",linewidth=2.0,color="#000000")   # P1
plot((time_archive - TD)*(60),state_archive[:,5]*(SF),"--",linewidth=2.0,color="#707070")   # P2
plot((time_archive - TD)*(60),state_archive[:,6]*(SF),"--",linewidth=2.0,color="#A9A9A9")   # P3

# write the axis -
xlabel("Time (min)",fontsize=18)
ylabel(L"Concentration ($\mu$M)",fontsize=18)

# save to disk -
savefig("./figs/Memory-Nominal-Discrete.pdf")
run(`pdfcrop ./figs/Memory-Nominal-Discrete.pdf ./figs/Memory-Nominal-Discrete.pdf`)
