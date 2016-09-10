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
function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# initialize the control -
	control_array = zeros(2)

	# Alias the species -
	gene_1 = x[1]
	gene_2 = x[2]
	mRNA_gene_1 = x[3]
	mRNA_gene_2 = x[4]
	protein_gene_1 = x[5]
	protein_gene_2 = x[6]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_1_gene_1 = binding_parameter_dictionary["n_gene_1_gene_1"]
	K_gene_1_gene_1 = binding_parameter_dictionary["K_gene_1_gene_1"]
	n_gene_1_gene_2 = binding_parameter_dictionary["n_gene_1_gene_2"]
	K_gene_1_gene_2 = binding_parameter_dictionary["K_gene_1_gene_2"]
	n_gene_2_gene_1 = binding_parameter_dictionary["n_gene_2_gene_1"]
	K_gene_2_gene_1 = binding_parameter_dictionary["K_gene_2_gene_1"]

	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_gene_1_gene_1 = control_parameter_dictionary["W_gene_1_gene_1"]
	W_gene_1_gene_2 = control_parameter_dictionary["W_gene_1_gene_2"]
	W_gene_2_gene_1 = control_parameter_dictionary["W_gene_2_gene_1"]

	# Control function for gene_1 -
	b_gene_1_gene_1 = ((protein_gene_1)^(n_gene_1_gene_1))/(K_gene_1_gene_1^(n_gene_1_gene_1)+protein_gene_1^(n_gene_1_gene_1))
	b_gene_1_gene_2 = ((protein_gene_2)^(n_gene_1_gene_2))/(K_gene_1_gene_2^(n_gene_1_gene_2)+protein_gene_2^(n_gene_1_gene_2))
	control_array[1] = (W_gene_1_gene_1*b_gene_1_gene_1)/(1+0.5+W_gene_1_gene_1*b_gene_1_gene_1+W_gene_1_gene_2*b_gene_1_gene_2)

	# Control function for gene_2 -
	b_gene_2_gene_1 = ((protein_gene_1)^(n_gene_2_gene_1))/(K_gene_2_gene_1^(n_gene_2_gene_1)+protein_gene_1^(n_gene_2_gene_1))
	control_array[2] = (W_gene_2_gene_1*b_gene_2_gene_1)/(1+0.5+W_gene_2_gene_1*b_gene_2_gene_1)

	return control_array
end
