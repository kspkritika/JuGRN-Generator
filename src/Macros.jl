macro include_function(src_file_name,pad_string)

  # create src_buffer -
  src_buffer::Array{AbstractString} = AbstractString[]

  @show pwd()

  # path to distrubtion -
  path_to_src_file = "./include/"*src_file_name*".jl"
  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)

        new_line_with_line_ending = line*"\n"
        push!(src_buffer,new_line_with_line_ending)
    end
  end

  string_value = ""
  for line in src_buffer
    string_value *= pad_string*line
  end

  return string_value
end
