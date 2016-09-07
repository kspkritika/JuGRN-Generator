function build_copyright_header_buffer(problem_object::ProblemObject)

  # What is the current year?
  current_year = string(Dates.year(now()))

  # Get comment data from
  buffer = ""
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  buffer*= "# Copyright (c) $(current_year) Varnerlab\n"
  buffer*= "# Robert Frederick School of Chemical and Biomolecular Engineering\n"
  buffer*= "# Cornell University, Ithaca NY 14850\n"
  buffer*= "#\n"
  buffer*= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
  buffer*= "# of this software and associated documentation files (the \"Software\"), to deal\n"
  buffer*= "# in the Software without restriction, including without limitation the rights\n"
  buffer*= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
  buffer*= "# copies of the Software, and to permit persons to whom the Software is\n"
  buffer*= "# furnished to do so, subject to the following conditions:\n"
  buffer*= "#\n"
  buffer*= "# The above copyright notice and this permission notice shall be included in\n"
  buffer*= "# all copies or substantial portions of the Software.\n"
  buffer*= "#\n"
  buffer*= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
  buffer*= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  buffer*= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
  buffer*= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  buffer*= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
  buffer*= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
  buffer*= "# THE SOFTWARE.\n"
  buffer*= "# ----------------------------------------------------------------------------------- #\n"

  # return -
  return buffer

end

function build_control_buffer(problem_object::ProblemObject)

  filename = "Control.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_data_dictionary_buffer(problem_object::ProblemObject)

  filename = "DataDictionary.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function DataDictionary()\n"
  buffer *= "\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_kinetics_buffer(problem_object::ProblemObject)

  filename = "Kinetics.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"
  buffer *= "\n"

  buffer *= "function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  buffer *= "function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  buffer *= "function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  buffer *= "function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_inputs_buffer(problem_object::ProblemObject)

  filename = "Inputs.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function calculate_input_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "\t# return - \n"
  buffer *= "\treturn zeros(length(x))\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end
