# ----------------------------------------------------------------------------------- #
# Copyright (c) 2018 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2018-02-13T05:26:59.088
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")
	dilution_matrix = readdlm("./Dilution.dat")
	degradation_matrix = readdlm("./Degradation.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of gene lengths -
	gene_coding_length_array = [
		1000	;	# 1	gene_1
		1200	;	# 2	gene_2
		900		;	# 3	gene_3
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 4	1	mRNA_gene_1
		gene_coding_length_array[2]	;	# 5	2	mRNA_gene_2
		gene_coding_length_array[3]	;	# 6	3	mRNA_gene_3
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 7	1	protein_gene_1
		round((0.33)*mRNA_coding_length_array[2])	;	# 8	2	protein_gene_2
		round((0.33)*mRNA_coding_length_array[3])	;	# 9	3	protein_gene_3
	]

	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                             # mum
	mass_of_single_cell = 2.8e-13                   # g
	number_of_rnapII = 4600            	            # copies/cells
	number_of_ribosome = 50000         	            # copies/cells
	mRNA_half_life_TF = 0.083                       # hrs
	protein_half_life = 20                          # hrs
	infrastructure_half_life = 300					# hrs
	doubling_time_cell = 0.33                       # hrs
	max_translation_rate = 16.5                     # aa/sec
	max_transcription_rate = 60.0                   # nt/sec
	transcription_initiation_time_contstant = 400  # sec
	average_transcript_length = 1200   	            # nt
	average_protein_length = 400       	            # aa
	fraction_nucleus = 0.0             	            # dimensionless
	av_number = 6.02e23                             # number/mol
	avg_gene_number = 4                             # number of copies of a gene
	polysome_number = 4					        # number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # mumol/gdw
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # mumol/ddw

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(e,0.5)                           # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(e,0.5)                        # hr^-1
	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(e,0.5)			# hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(e,2)                          # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/mass_of_single_cell)*1e9      # nmol/gdw

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 10*(1/av_number)*(1/mass_of_single_cell)*1e9               	# nmol/gdw
	saturation_translation = 750000*(1/av_number)*(1/mass_of_single_cell)*1e6               	# mumol/gdw
	# -------------------------------------------------------------------------------------------#

	# initial condition array -
	initial_condition_array = [

		avg_gene_concentration	;	# 1	gene_1
		avg_gene_concentration	;	# 2	gene_2
		avg_gene_concentration	;	# 3	gene_3
		0.0	;	# 4	mRNA_gene_1
		0.0	;	# 5	mRNA_gene_2
		0.0	;	# 6	mRNA_gene_3
		0.0	;	# 7	protein_gene_1
		0.0	;	# 8	protein_gene_2
		0.0	;	# 9	protein_gene_3
	]

	binding_parameter_dictionary = Dict{AbstractString,Float64}()
	binding_parameter_dictionary["n_gene_2_gene_1"] = 1.0
	binding_parameter_dictionary["K_gene_2_gene_1"] = 0.05
	binding_parameter_dictionary["n_gene_2_gene_3"] = 1.0
	binding_parameter_dictionary["K_gene_2_gene_3"] = 0.05
	binding_parameter_dictionary["n_gene_3_gene_1"] = 1.0
	binding_parameter_dictionary["K_gene_3_gene_1"] = 0.05
	binding_parameter_dictionary["n_gene_3_gene_2"] = 1.0
	binding_parameter_dictionary["K_gene_3_gene_2"] = 0.05

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{AbstractString,Float64}()
	control_parameter_dictionary["W_gene_1_RNAP"] = 0.0
	control_parameter_dictionary["W_gene_2_RNAP"] = 0.0
	control_parameter_dictionary["W_gene_2_gene_1"] = 10.0
	control_parameter_dictionary["W_gene_2_gene_3"] = 10.0
	control_parameter_dictionary["W_gene_3_RNAP"] = 0.0
	control_parameter_dictionary["W_gene_3_gene_1"] = 10.0
	control_parameter_dictionary["W_gene_3_gene_2"] = 10.0

	# Alias the txtl parameters -
	txtl_parameter_dictionary = Dict{AbstractString,Float64}()
	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation


	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_gene_2_gene_1"	;	# 1
		"K_gene_2_gene_1"	;	# 2
		"n_gene_2_gene_3"	;	# 3
		"K_gene_2_gene_3"	;	# 4
		"n_gene_3_gene_1"	;	# 5
		"K_gene_3_gene_1"	;	# 6
		"n_gene_3_gene_2"	;	# 7
		"K_gene_3_gene_2"	;	# 8
		"W_gene_1_RNAP"	;	# 9
		"W_gene_2_RNAP"	;	# 10
		"W_gene_2_gene_1"	;	# 11
		"W_gene_2_gene_3"	;	# 12
		"W_gene_3_RNAP"	;	# 13
		"W_gene_3_gene_1"	;	# 14
		"W_gene_3_gene_2"	;	# 15
		"rnapII_concentration"	;	# 16
		"ribosome_concentration"	;	# 17
		"degradation_constant_mRNA"	;	# 18
		"degradation_constant_protein"	;	# 19
		"kcat_transcription"	;	# 20
		"kcat_translation"	;	# 21
		"maximum_specific_growth_rate"	;	# 22
		"saturation_constant_transcription"	;	# 23
		"saturation_constant_translation"	;	# 24
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["volume_of_single_cell"] = V
	data_dictionary["mass_of_single_cell"] = mass_of_single_cell
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["average_transcript_length"] = average_transcript_length
	data_dictionary["average_protein_length"] = average_protein_length
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	data_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	data_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	data_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	data_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	data_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	data_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation  # hr^-1
	data_dictionary["kcat_translation_initiation"] = kcat_translation_initiation  # hr^-1
	data_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	data_dictionary["death_rate_constant"] = death_rate_constant
	data_dictionary["avg_gene_concentration"] = avg_gene_concentration
	data_dictionary["saturation_constant_transcription"] = saturation_transcription
	data_dictionary["saturation_constant_translation"] = saturation_translation

	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_matrix"] = dilution_matrix
	data_dictionary["degradation_matrix"] = degradation_matrix

	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["txtl_parameter_dictionary"] = txtl_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
