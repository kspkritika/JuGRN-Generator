# includes -
include("Kinetics.jl")
include("Control.jl")

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
#
# ----------------------------------------------------------------------------------- #
# Balances: Evaluates model equations given time, state and the data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# t  - current time
# x  - state array
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# dxdt - derivative array at current time step
# ----------------------------------------------------------------------------------- #
function Balances(t,x,data_dictionary)

  # correct for negatives -
  idx_small = find(x.<0)
  x[idx_small] = 0.0

  # Get model matricies and other required data from the data_dictionary -
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
	dilution_matrix = data_dictionary["dilution_matrix"]
	degradation_matrix = data_dictionary["degradation_matrix"]
  mugmax = data_dictionary["maximum_specific_growth_rate"]

  # Calculate the kinetics -
  transcription_rate_array = calculate_transcription_rates(t,x,data_dictionary)
  translation_rate_array = calculate_translation_rates(t,x,data_dictionary)
  mRNA_degradation_rate_array = calculate_mRNA_degradation_rates(t,x,data_dictionary)
  protein_degradation_rate_array = calculate_protein_degradation_rates(t,x,data_dictionary)
  background_transcription_rate_array = calculate_background_transcription_rates(t,x,data_dictionary)

  # Call my control function -
  control_array = Control(t,x,data_dictionary)

  # Compute the modified rate -
  transcription_rate_array = transcription_rate_array.*control_array;

  # rate array -
  txtl_rate_array = [transcription_rate_array ; translation_rate_array]
  degradation_rate_array = [mRNA_degradation_rate_array ; protein_degradation_rate_array]

  # Evaluate the balance equations -
  dxdt_array = stoichiometric_matrix*txtl_rate_array+degradation_matrix*degradation_rate_array+mugmax*dilution_matrix*x+background_transcription_rate_array

  # return -
  return dxdt_array
end
