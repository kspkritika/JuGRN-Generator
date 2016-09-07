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
function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
end

function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

end


function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

end


function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	mRNA_gene_1 = x[3]
	mRNA_gene_2 = x[4]

	# Initialize the degrdation array - 
	degradation_rate_array = zeros(2)
	mRNA_degrdation_constant = data_dictionary["degradation_constant_mRNA"]

	# Calculate the degradation_rate_array - 
	degradation_rate_array[1] = (mRNA_degrdation_constant)*mRNA_gene_1
	degradation_rate_array[2] = (mRNA_degrdation_constant)*mRNA_gene_2

	# return the degrdation rate array - 
	return degradation_rate_array
end

function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	protein_gene_1 = x[5]
	protein_gene_2 = x[6]

	# Initialize the degrdation array - 
	degradation_rate_array = zeros(2)
	protein_degrdation_constant = data_dictionary["degradation_constant_protein"]

	# Calculate the degradation_rate_array - 
	degradation_rate_array[1] = (protein_degrdation_constant)*protein_gene_1
	degradation_rate_array[2] = (protein_degrdation_constant)*protein_gene_2

	# return the degrdation rate array - 
	return degradation_rate_array
end


