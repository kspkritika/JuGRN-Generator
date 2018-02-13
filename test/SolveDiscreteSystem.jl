# include -
include("Include.jl")
include("Discrete.jl")

function main()

    # Script to solve the discrete system of equations -
    # Load the data dictionary -
    data_dictionary = DataDictionary(0.0,0.0,0.0)

    # How many genes do we have?
    gene_coding_length_array = data_dictionary["gene_coding_length_array"]
    number_of_genes = length(gene_coding_length_array)

    # Calculate the system arrays -
    time_step_size = 0.01
    (AM,BM,DM) = generate_discrete_system(data_dictionary,time_step_size)

    # initilize the state/time archives -
    (number_of_species,number_of_species) = size(AM)
    state_archive = zeros(1,number_of_species)
    time_archive = Float64[]

    # what is my global IC -
    IC = data_dictionary["initial_condition_array"][(number_of_genes+1):end]
    for state_index = 1:number_of_species
        state_archive[1,state_index] = IC[state_index]
    end

    # push the initial time -
    push!(time_archive,0.0)

    # ======================================================================== #
    # Phase 1: Run the model for a few time units
    time_start = 0.0
    time_stop = 1.0
    time_step_size = 0.01
    time_array = collect(time_start:time_step_size:time_stop)

    xold = state_archive[end,:]
    for time_value in time_array

        # calculate the new state -
        xnew = evaluate_discrete_system(time_value,xold,AM,BM,data_dictionary)

        # cache the new state and time -
        state_archive = archive_solution_array(state_archive,xnew)
        push!(time_archive,time_value)

        # update the old state -
        xold = xnew
    end
    # ======================================================================== #
    #
    # ======================================================================== #
    # Phase 2: Run the model for a few time units w/P1 expression
    time_start = time_stop
    time_stop = time_start + 3.0
    time_step_size = 0.01
    time_array = collect(time_start:time_step_size:time_stop)

    # update the W1_RNAP parameter -
    data_dictionary["control_parameter_dictionary"]["W_gene_1_RNAP"] = 10.0

    # The end state of Phase 1 is the start to Phase 2 -
    xold = state_archive[end,:]

    # Solve the model equations in phase 2
    for time_value in time_array

        # calculate the new state -
        xnew = evaluate_discrete_system(time_value,xold,AM,BM,data_dictionary)

        # cache the new state and time -
        state_archive = archive_solution_array(state_archive,xnew)
        push!(time_archive,time_value)

        # update the old state -
        xold = xnew
    end
    # ======================================================================== #
    #
    # ======================================================================== #
    # Phase 3: Run the model for a few time units w/o P1 expression
    time_start = time_stop
    time_stop = time_start + 8.0
    time_step_size = 0.01
    time_array = collect(time_start:time_step_size:time_stop)

    # update the W1_RNAP parameter -
    data_dictionary["control_parameter_dictionary"]["W_gene_1_RNAP"] = 0.0

    # The end state of Phase 2 is the start to Phase 3 -
    xold = state_archive[end,:]

    # Solve the model equations in phase 3
    for time_value in time_array

        # calculate the new state -
        xnew = evaluate_discrete_system(time_value,xold,AM,BM,data_dictionary)

        # cache the new state and time -
        state_archive = archive_solution_array(state_archive,xnew)
        push!(time_archive,time_value)

        # update the old state -
        xold = xnew
    end
    # ======================================================================== #

    # return the archives -
    return (time_archive,state_archive)
end

# call main -
(time_archive,state_archive) = main()

# dump results to the tmp folder -
data_archive = [time_archive state_archive]
writedlm("./tmp/simulation_nominal_discrete.dat",data_archive)
