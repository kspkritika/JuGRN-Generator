# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
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
function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)

	# initial condition array - 
	initial_condition_array = [
		2.0	;	#1	gene_1
		2.0	;	#2	gene_2
		0.0	;	#3	mRNA_gene_1
		0.0	;	#4	mRNA_gene_2
		0.0	;	#5	protein_gene_1
		0.0	;	#6	protein_gene_2
	]

	# array of gene lengths - 
	gene_coding_length_array = [
		15000	;	#1	gene_1
		15000	;	#2	gene_2
	]

	# array of mRNA coding lengths - 
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	#3	1	mRNA_gene_1
		gene_coding_length_array[2]	;	#4	2	mRNA_gene_2
	]

	# array of mRNA coding lengths - 
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[5])	;	#5	1	protein_gene_1
		round((0.33)*mRNA_coding_length_array[6])	;	#6	2	protein_gene_2
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
