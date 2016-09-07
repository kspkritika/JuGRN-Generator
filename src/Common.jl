function write_program_components_to_disk(file_path::AbstractString,set_of_program_components::Set{ProgramComponent})

  # go through each component, and dump the buffer to disk -
  for program_component in set_of_program_components

    # Get the data -
    filename = program_component.filename
    program_buffer = program_component.buffer

    # build the path -
    path_to_program_file = file_path*"/"*filename

    # Write the file -
    outfile = open(path_to_program_file, "w")
    write(outfile,program_buffer);
    close(outfile);
  end
end

# function build_reaction_comment_string(reaction_object::ReactionObject)
#
#   # Ok, let's build a comment from the reactioj object {Reactants} -> {Products}
#   list_of_reactants = reaction_object.list_of_reactants
#   list_of_products = reaction_object.list_of_products
#
#   # Build reactant string -
#   reactant_buffer = ""
#   for species_object::SpeciesObject in list_of_reactants
#
#     # symbol etc -
#     species_symbol = species_object.species_symbol
#     stoichiometric_coefficient = species_object.stoichiometric_coefficient
#
#     if (stoichiometric_coefficient == 1.0)
#       # fill reactant buffer -
#       reactant_buffer *= "$(species_symbol)+"
#     else
#       # fill reactant buffer -
#       reactant_buffer *= "$(stoichiometric_coefficient)\*$(species_symbol)+"
#     end
#   end
#
#   # Cutoff the trailing +
#   reactant_buffer = reactant_buffer[1:end-1]
#
#   # Check for empty ...
#   if (length(list_of_reactants) == 0)
#     reactant_buffer = "[]"
#   end
#
#   # Build product string -
#   product_buffer = ""
#   for species_object::SpeciesObject in list_of_products
#
#     # symbol etc -
#     species_symbol = species_object.species_symbol
#     stoichiometric_coefficient = species_object.stoichiometric_coefficient
#
#     if (stoichiometric_coefficient == 1.0)
#       # fill product buffer -
#       product_buffer *= "$(species_symbol)+"
#     else
#       # fill product buffer -
#       product_buffer *= "$(stoichiometric_coefficient)\*$(species_symbol)+"
#     end
#   end
#
#   # Cutoff the trailing +
#   product_buffer = product_buffer[1:end-1]
#
#   # Check for empty ...
#   if (length(list_of_products) == 0)
#     product_buffer = "[]"
#   end
#
#   # comment string -
#   comment_buffer = ""
#   comment_buffer = reactant_buffer*" --> "*product_buffer
#   return comment_buffer
# end

function transfer_distribution_files(path_to_distribution_files::AbstractString,
                                      path_to_output_files::AbstractString,
                                      file_extension::AbstractString)


  # Search the directory for src files -
  # load the files -
  searchdir(path,key) = filter(x->contains(x,key),readdir(path))

  # build src file list -
  list_of_src_files = searchdir(path_to_distribution_files,file_extension)

  # go thru the src file list, and copy the files to the output path -
  for src_file in list_of_src_files

    # create src_buffer -
    src_buffer::Array{AbstractString} = AbstractString[]

    # path to distrubtion -
    path_to_src_file = path_to_distribution_files*"/"*src_file
    open(path_to_src_file,"r") do src_file
      for line in eachline(src_file)
        push!(src_buffer,line)
      end
    end

    # Write the file to the output -
    path_to_program_file = path_to_output_files*"/"*src_file
    outfile = open(path_to_program_file, "w")
    write(outfile,src_buffer);
    close(outfile);
  end
end

function partition!(list_of_species::Array{SpeciesObject})

  # ok, frist, we need to split into balanced and unbalanced lists -
  list_of_gene_indexes::Array{Int} = Int[]
  list_of_mRNA_indexes::Array{Int} = Int[]
  list_of_protein_indexes::Array{Int} = Int[]

  for (index,species_object) in enumerate(list_of_species)

    species_type::Symbol = species_object.species_type
    if (species_type == :gene)
      push!(list_of_gene_indexes,index)
    elseif (species_type == :mrna)
      push!(list_of_mRNA_indexes,index)
    elseif (species_type == :protein)
      push!(list_of_protein_indexes,index)
    end
  end

  # combine -
  permutation_index_array = vcat(list_of_gene_indexes,list_of_mRNA_indexes,list_of_protein_indexes)

  # permute the array -
  permute!(list_of_species,permutation_index_array)
end
