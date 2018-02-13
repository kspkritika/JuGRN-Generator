function adj_induction_simulation(time_start,time_stop,time_step,parameter_index,data_dictionary)

  # run the model to steady-state -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS
  number_of_states = length(XSS)
  initial_condition_array = [initial_condition_array ; zeros(number_of_states)]
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # update the W1_RNAP parameter of gene 1 to force induction -
  data_dictionary["control_parameter_dictionary"]["W_gene_1_RNAP"] = 10.0

  # run the model and the sensitivity coefficients -
  (T,X) = SolveAdjBalances(time_start,time_stop,time_step,parameter_index,data_dictionary)

  # return the time and state -
  return (T,X)
end
