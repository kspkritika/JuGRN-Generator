include("Balances.jl")
using ODE

# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA

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
function SolveBalances(TSTART,TSTOP,Ts,data_dictionary)
# ----------------------------------------------------------------------------------- #
# SolveBalances.jl was generated using the Kwatee code generation system.
# SolveBalances: Solves model equations from TSTART to TSTOP given parameters in data_dictionary.
# Username: jeffreyvarner
# Type: GRN-JULIA
# Version: 1.0
# Generation timestamp: 11-24-2015 11:18:08
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# TSIM - Simulation time vector
# X - Simulation state array (NTIME x NSPECIES)
# ----------------------------------------------------------------------------------- #

# Get required stuff from DataFile struct -
TSIM = collect(TSTART:Ts:TSTOP);
initial_condition_vector = data_dictionary["initial_condition_array"];

# Call the ODE solver -
fbalances(t,y) = Balances(t,y,data_dictionary);
(t,y) = ode23s(fbalances,initial_condition_vector,TSIM;points=:specified);
#(t,y) = ode4s(fbalances,initial_condition_vector,TSIM);

# Map -
# number_of_timesteps = length(t)
# number_of_states = length(initial_condition_vector)
# X = zeros(number_of_timesteps,number_of_states)
# for state_index = 1:number_of_states
#   tmp = map(y->y[state_index],y)
#   for time_index = 1:number_of_timesteps
#     X[time_index,state_index] = tmp[time_index]
#   end
# end

X1 = map(y->y[1],y)
X2 = map(y->y[2],y)
X3 = map(y->y[3],y)
X4 = map(y->y[4],y)
X5 = map(y->y[5],y)
X6 = map(y->y[6],y)
X7 = map(y->y[7],y)
X8 = map(y->y[8],y)
X9 = map(y->y[9],y)
X10 = map(y->y[10],y)
X11 = map(y->y[11],y)
X12 = map(y->y[12],y)
X13 = map(y->y[13],y)
X14 = map(y->y[14],y)
X15 = map(y->y[15],y)
X16 = map(y->y[16],y)
X17 = map(y->y[17],y)
X18 = map(y->y[18],y)

X = [X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18]

idx_n = find(X.<0)
X[idx_n] = 0.0

return (t,X);
end
