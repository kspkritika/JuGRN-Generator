# includes -
include("Include.jl")

# load the data_dictionary -
data_dictionary = DataDictionary(0.0,0.0,0.0)

# where do we find files?
path_to_senstivity_files = "./sensitivity"
file_pattern = "AdjSimulation-P"
time_skip = 40

# calc the sensitivity array -
# (T,SA) = calculate_sensitivity_array(path_to_senstivity_files,file_pattern,time_skip,data_dictionary)
ASSA = calculate_average_scaled_sensitivity_array(path_to_senstivity_files,file_pattern,data_dictionary)
