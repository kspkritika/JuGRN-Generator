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

  # Calculate the kinetics -
  transcription_rate_array = calculate_transcription_rates(t,x,data_dictionary)
  translation_rate_array = calculate_translation_rates(t,x,data_dictionary)
  mRNA_degradation_rate_array = calculate_mRNA_degradation_rates(t,x,data_dictionary)
  protein_degradation_rate_array = calculate_protein_degradation_rates(t,x,data_dictionary)

  # Call my control function -
  control_array = Control(t,x,data_dictionary)

  # Compute the modified rate -
  transcription_rate_array = transcription_rate_array.*control_array;
  dxdt_array = zeros(18)

  # define the balance equations -
  dxdt_array[1]  = -rate_array[2]-rate_array[1]                    # 1 Prothrombin FII
  dxdt_array[2]  = rate_array[2]+rate_array[1]-rate_array[4]       # 2 thrombin FIIa
  dxdt_array[3]  = -rate_array[3]                                  # 3 PC
  dxdt_array[4]  = rate_array[3]-rate_array[14]                    # 4 APC
  dxdt_array[5]  = -rate_array[4]                                  # 5 ATIII
  dxdt_array[6]  = 0.0                                             # 6 TM acts as an enzyme
  dxdt_array[7]  = 0.0                                             # 7 Trigger
  dxdt_array[8] = rate_array[9]-rate_array[12]                     # 8 Fibrin
  dxdt_array[9] = rate_array[10]+rate_array[11]-rate_array[13]     # 9 Plasmin
  dxdt_array[10] = -rate_array[5] - rate_array[18]                 # 10 Fibrinogen
  dxdt_array[11] = -rate_array[10] - rate_array[11]                # 11 Plasminogen
  dxdt_array[12] = -rate_array[15]                                 # 12 tPA
  dxdt_array[13] = -rate_array[16]                                 # 13 uPA
  dxdt_array[14] = 1.0*(rate_array[5]-rate_array[6])               # 14 Fibrin monomer
  dxdt_array[15] = rate_array[6]+rate_array[8]-rate_array[7]       # 15 Protofibril
  dxdt_array[16] = -rate_array[13]                                 # 16 Antiplasmin
  dxdt_array[17] = -rate_array[14]-rate_array[15]-rate_array[16]   # 17 PAI_1
  dxdt_array[18] = 1.0*(rate_array[7]-rate_array[9]-rate_array[17])# 18 Fibers

  # return -
  return dxdt_array
end
