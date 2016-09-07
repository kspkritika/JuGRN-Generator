using ArgParse
using JSON

include("./src/Types.jl")
include("./src/Macros.jl")
include("./src/Parser.jl")
include("./src/Problem.jl")
include("./src/Common.jl")
include("./src/strategy/JuliaStrategy.jl")

# Grab the required functions for code generation -
# const parser_function::Function = parse_vff_file

function parse_commandline()
    settings_object = ArgParseSettings()
    @add_arg_table settings_object begin
      "-o"
        help = "Directory where the Julia model files will be written."
        arg_type = AbstractString
        default = "."

      "-m"
        help = "Path to the grn file written in the vgrn format."
        arg_type = AbstractString
        required = true
    end

    # return a dictionary w/args -
    return parse_args(settings_object)
end


function main()

  # Build the arguement dictionary -
  parsed_args = parse_commandline()

  # Load the statement_vector -
  path_to_model_file = parsed_args["m"]
  statement_vector::Array{VGRNSentence} = parse_grn_file(path_to_model_file)

  # Generate the problem object -
  problem_object = generate_problem_object(statement_vector)

  # initialize the program set -
  component_set = Set{ProgramComponent}()

  # Write the Inputs -
  program_component_inputs = build_inputs_buffer(problem_object)
  push!(component_set,program_component_inputs)

  # Write the data_dictionary -
  program_component_data_dictionary = build_data_dictionary_buffer(problem_object)
  push!(component_set,program_component_data_dictionary)

  # Write the Kinetics -
  program_component_kinetics = build_kinetics_buffer(problem_object)
  push!(component_set,program_component_kinetics)

  # Write the Control -
  program_component_control = build_control_buffer(problem_object)
  push!(component_set,program_component_control)

  # Dump the component_set to disk -
  path_to_output_file = parsed_args["o"]
  write_program_components_to_disk(path_to_output_file,component_set)

  # Transfer distrubtion files to the output -
  transfer_distribution_files("./distribution",path_to_output_file,".jl")
end

main()
