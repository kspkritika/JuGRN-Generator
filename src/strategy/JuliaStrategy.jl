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

@debug function build_control_buffer(problem_object::ProblemObject)

  filename = "Control.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # extract list of genes -
  list_of_species = problem_object.list_of_species
  list_of_genes = extract_species_of_type(list_of_species,:gene)

  # get list of connections -
  list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "\t# initialize the control - \n"
  buffer *= "\tcontrol_array = zeros($(length(list_of_genes)))\n"
  buffer *= "\n"

  buffer *= "\t# Alias the species - \n"
  for (index,species_object) in enumerate(list_of_species)

    # Grab the symbol -
    species_symbol = species_object.species_symbol

    # write the record -
    buffer *= "\t$(species_symbol) = x[$(index)]\n"
  end

  buffer *= "\n"
  buffer *= "\t# Alias the binding parameters - \n"
  buffer *= "\tbinding_parameter_dictionary = data_dictionary[\"binding_parameter_dictionary\"]\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tn_$(gene_symbol)_$(actor_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(actor_symbol)\"]\n"
        buffer *= "\tK_$(gene_symbol)_$(actor_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(actor_symbol)\"]\n"
      end
    end

    for connection_object in inhibiting_connections

      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tn_$(gene_symbol)_$(actor_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(actor_symbol)\"]\n"
        buffer *= "\tK_$(gene_symbol)_$(actor_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(actor_symbol)\"]\n"
      end
    end
  end

  buffer *= "\n"
  buffer *= "\t# Alias the control function parameters - \n"
  buffer *= "\tcontrol_parameter_dictionary = data_dictionary[\"control_parameter_dictionary\"]\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tW_$(gene_symbol)_$(actor_symbol) = control_parameter_dictionary[\"W_$(gene_symbol)_$(actor_symbol)\"]\n"
      end
    end

    for connection_object in inhibiting_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tW_$(gene_symbol)_$(actor_symbol) = control_parameter_dictionary[\"W_$(gene_symbol)_$(actor_symbol)\"]\n"
      end
    end
  end

  buffer *= "\n"
  # get list of genes -
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # generate the binding functions -
    buffer *= "\t# Control function for $(gene_symbol) - \n"
    for connection_object in activating_connections

      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list

        actor_symbol = actor_object.species_symbol
        buffer *= "\tb_$(gene_symbol)_$(actor_symbol) = ((protein_$(actor_symbol))^(n_$(gene_symbol)_$(actor_symbol)))/"
        buffer *= "(K_$(gene_symbol)_$(actor_symbol)^(n_$(gene_symbol)_$(actor_symbol))+protein_$(actor_symbol)^(n_$(gene_symbol)_$(actor_symbol)))\n"
      end
    end

    for connection_object in inhibiting_connections

      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list

        actor_symbol = actor_object.species_symbol
        buffer *= "\tb_$(gene_symbol)_$(actor_symbol) = ((protein_$(actor_symbol))^(n_$(gene_symbol)_$(actor_symbol)))/"
        buffer *= "(K_$(gene_symbol)_$(actor_symbol)^(n_$(gene_symbol)_$(actor_symbol))+protein_$(actor_symbol)^(n_$(gene_symbol)_$(actor_symbol)))\n"
      end
    end

    buffer *= "\tcontrol_array[$(index)] = ("
    numerator = ""
    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        numerator *= "W_$(gene_symbol)_$(actor_symbol)*b_$(gene_symbol)_$(actor_symbol)+"
      end
    end
    buffer *= numerator[1:end-1]
    buffer *= ")/(1+"

    demoninator = ""
    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        demoninator *= "W_$(gene_symbol)_$(actor_symbol)*b_$(gene_symbol)_$(actor_symbol)+"
      end
    end

    # ok - do we have inhibitory statements?
    if (isempty(inhibiting_connections) == true)
      buffer *= demoninator[1:end-1]
      buffer *= ")\n"
    else
      buffer *= demoninator
      demoninator = ""
      for connection_object in inhibiting_connections
        # actor -
        actor_list = connection_object.connection_actor_set
        for actor_object in actor_list
          actor_symbol = actor_object.species_symbol
          demoninator *= "W_$(gene_symbol)_$(actor_symbol)*b_$(gene_symbol)_$(actor_symbol)+"
        end
      end
      buffer *= demoninator[1:end-1]
      buffer *= ")\n"
      buffer *= "\n"
    end
  end

  buffer *= "\n"
  buffer *= "\treturn control_array\n"
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

  # get list of species from the po -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)\n"
  buffer *= "\n"
  buffer *= "\t# stoichiometric_matrix and dilution_matrix - \n"
  buffer *= "\tstoichiometric_matrix = readdlm(\"./Network.dat\")\n"
  buffer *= "\tdilution_matrix = readdlm(\"./Dilution.dat\")\n"
  buffer *= "\tdegradation_matrix = readdlm(\"./Degradation.dat\")\n"
  buffer *= "\n"
  buffer *= "\t# array of gene lengths - \n"
  buffer *= "\tgene_coding_length_array = [\n"

  # write out the length of genes -
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :gene)
      buffer *= "\t\t15000\t;\t# $(index)\t$(species_symbol)\n"
    end
  end
  buffer *= "\t]\n"
  buffer *= "\n"
  buffer *= "\t# array of mRNA coding lengths - \n"
  buffer *= "\tmRNA_coding_length_array = [\n"

  # write out the length of genes -
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :mrna)
      buffer *= "\t\tgene_coding_length_array[$(counter)]\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
      counter = counter + 1
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"

  buffer *= "\t# array of mRNA coding lengths - \n"
  buffer *= "\tprotein_coding_length_array = [\n"

  # write out the length of genes -
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :protein)
      buffer *= "\t\tround((0.33)*mRNA_coding_length_array[$(counter)])\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
      counter = counter + 1
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"
  buffer *= @include_function("txtl_constants","\t")
  buffer *= "\n"

  buffer *= "\t# initial condition array - \n"
  buffer *= "\tinitial_condition_array = [\n"

  # write out species -
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :gene)
      buffer *= "\t\tavg_gene_concentration\t;\t# $(index)\t$(species_symbol)\n"
    elseif (species_type == :mrna || species_type == :protein)
      buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"

  # Setup parameters dictionaries -
  list_of_genes = extract_species_of_type(list_of_species,:gene)

  # get list of connections -
  list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections
  buffer *= "\tbinding_parameter_dictionary = Dict{AbstractString,Float64}()\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tbinding_parameter_dictionary[\"n_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
        buffer *= "\tbinding_parameter_dictionary[\"K_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
      end
    end

    for connection_object in inhibiting_connections

      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tbinding_parameter_dictionary[\"n_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
        buffer *= "\tbinding_parameter_dictionary[\"K_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
      end
    end
  end

  buffer *= "\n"
  buffer *= "\t# Alias the control function parameters - \n"
  buffer *= "\tcontrol_parameter_dictionary = Dict{AbstractString,Float64}()\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    for connection_object in activating_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
      end
    end

    for connection_object in inhibiting_connections
      # actor -
      actor_list = connection_object.connection_actor_set
      for actor_object in actor_list
        actor_symbol = actor_object.species_symbol
        buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_$(actor_symbol)\"] = 1.0\n"
      end
    end
  end

  buffer *= "\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = Dict{AbstractString,Any}()\n"
  buffer *= "\tdata_dictionary[\"initial_condition_array\"] = initial_condition_array\n"
  buffer *= "\tdata_dictionary[\"gene_coding_length_array\"] = gene_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"mRNA_coding_length_array\"] = mRNA_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"protein_coding_length_array\"] = protein_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"rnapII_concentration\"] = rnapII_concentration  # \muM \n"
  buffer *= "\tdata_dictionary[\"ribosome_concentration\"] = ribosome_concentration # \muM \n"
  buffer *= "\tdata_dictionary[\"degradation_constant_mRNA\"] = degradation_constant_mRNA  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"degradation_constant_protein\"] = degradation_constant_protein  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"kcat_transcription\"] = kcat_transcription  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"kcat_translation\"] = kcat_translation  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"maximum_specific_growth_rate\"] = maximum_specific_growth_rate  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"death_rate_constant\"] = death_rate_constant \n"
  buffer *= "\tdata_dictionary[\"avg_gene_concentration\"] = avg_gene_concentration \n"
  buffer *= "\tdata_dictionary[\"saturation_constant_transcription\"] = saturation_transcription \n"
  buffer *= "\tdata_dictionary[\"saturation_constant_translation\"] = saturation_translation \n"
  buffer *= "\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"dilution_matrix\"] = dilution_matrix\n"
  buffer *= "\tdata_dictionary[\"degradation_matrix\"] = degradation_matrix\n"
  buffer *= "\n"
  buffer *= "\tdata_dictionary[\"binding_parameter_dictionary\"] = binding_parameter_dictionary\n"
  buffer *= "\tdata_dictionary[\"control_parameter_dictionary\"] = control_parameter_dictionary\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"
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

  # get list of species from the po -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_genes = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :gene)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_genes = number_of_genes + 1
    end
  end
  buffer *="\n"
  buffer *="\t# Initialize the transcription rate - \n"
  buffer *="\ttranscription_rate_array = zeros($(number_of_genes))\n"
  buffer *="\tKSAT = data_dictionary[\"saturation_constant_transcription\"]\n"
  buffer *="\tkcat_transcription = data_dictionary[\"kcat_transcription\"]\n"
  buffer *="\trnapII_concentration = data_dictionary[\"rnapII_concentration\"]\n"
  buffer *="\n"
  buffer *="\t# Populate the transcription rate array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :gene)
      buffer *= "\ttranscription_rate_array[$(counter)] = kcat_transcription*(rnapII_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
      counter = counter + 1
    end
  end

  buffer *="\n"
  buffer *= "\t# return transcription_rate_array - \n"
  buffer *= "\treturn transcription_rate_array\n"

  buffer *= "end\n"
  buffer *= "\n"


  # calculate_background_transcription_rates -
  buffer *= "function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},transcription_rate_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\treturn zeros(length(x))\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  buffer *= "function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_mRNA = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_mRNA = number_of_mRNA + 1
    end
  end
  buffer *="\n"
  buffer *="\t# Initialize the translation rate - \n"
  buffer *="\ttranslation_rate_array = zeros($(number_of_mRNA))\n"
  buffer *="\tKSAT = data_dictionary[\"saturation_constant_translation\"]\n"
  buffer *="\tkcat_translation = data_dictionary[\"kcat_translation\"]\n"
  buffer *="\tribosome_concentration = data_dictionary[\"ribosome_concentration\"]\n"
  buffer *="\n"
  buffer *="\t# Populate the translation rate array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\ttranslation_rate_array[$(counter)] = kcat_translation*(ribosome_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
      counter = counter + 1
    end
  end

  buffer *="\n"
  buffer *= "\t# return translation array - \n"
  buffer *= "\treturn translation_rate_array\n"
  buffer *= "end\n"
  buffer *= "\n"


  # calculate_mRNA_degradation_rates -
  buffer *= "function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_mRNA = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_mRNA = number_of_mRNA + 1
    end
  end
  buffer *= "\n"
  buffer *="\t# Initialize the degrdation array - \n"
  buffer *="\tdegradation_rate_array = zeros($(number_of_mRNA))\n"
  buffer *="\tmRNA_degrdation_constant = data_dictionary[\"degradation_constant_mRNA\"]\n"
  buffer *= "\n"
  buffer *="\t# Calculate the degradation_rate_array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\tdegradation_rate_array[$(counter)] = (mRNA_degrdation_constant)*$(species_symbol)\n"
      counter = counter + 1
    end
  end
  buffer *= "\n"
  buffer *= "\t# return the degrdation rate array - \n"
  buffer *= "\treturn degradation_rate_array\n"
  buffer *= "end\n"
  buffer *= "\n"

  # calculate_protein_degradation_rates
  buffer *= "function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *="\t# Alias the species - \n"
  number_of_proteins = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :protein)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_proteins = number_of_proteins + 1
    end
  end
  buffer *= "\n"
  buffer *="\t# Initialize the degrdation array - \n"
  buffer *="\tdegradation_rate_array = zeros($(number_of_proteins))\n"
  buffer *="\tprotein_degrdation_constant = data_dictionary[\"degradation_constant_protein\"]\n"
  buffer *= "\n"
  buffer *="\t# Calculate the degradation_rate_array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :protein)
      buffer *= "\tdegradation_rate_array[$(counter)] = (protein_degrdation_constant)*$(species_symbol)\n"
      counter = counter + 1
    end
  end
  buffer *= "\n"
  buffer *= "\t# return the degrdation rate array - \n"
  buffer *= "\treturn degradation_rate_array\n"
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
