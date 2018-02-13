# include -
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# setup some constants -
const epsilon = 1e-6

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# ======================================================================== #
# Phase 0: estimate the steady-state -
XSS = estimate_steady_state(epsilon,data_dictionary)

# update the initial conditons to the steady state -
data_dictionary["initial_condition_array"] = XSS
# ======================================================================== #
#
# ======================================================================== #
# Phase 1: Run the model for a few time units
time_start = 0.0
time_stop = 1.0
time_step_size = 0.01

# Solve the model equations in phase 2
(T1,X1) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
# ======================================================================== #
#
# ======================================================================== #
# Phase 2: Run the model for a few time units
time_start = T1[end]
time_stop = T1[end] + 3.0
time_step_size = 0.01

# update the initial condition to X1 -
data_dictionary["initial_condition_array"] = X1[end,:]

# update the W1_RNAP parameter -
data_dictionary["control_parameter_dictionary"]["W_gene_1_RNAP"] = 10.0
data_dictionary["control_parameter_dictionary"]["W_gene_3_gene_2"] = 0.0

# Solve the model equations in phase 2
(T2,X2) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
# ======================================================================== #
#
# ======================================================================== #
# Phase 3: Run the model without inducer -
time_start = T2[end]
time_stop = T2[end] + 8.0
time_step_size = 0.01

# update the initial condition to X1 -
data_dictionary["initial_condition_array"] = X2[end,:]

# update the W1_RNAP parameter -
data_dictionary["control_parameter_dictionary"]["W_gene_1_RNAP"] = 0.0

# Solve the model equations in phase 2
(T3,X3) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
# ======================================================================== #

# Package -
time_archive = [T1 ; T2 ; T3]
state_archive = [X1 ; X2 ; X3]

# write to disk -
data_archive = [time_archive state_archive]
writedlm("./tmp/simulation_broken.dat",data_archive)
