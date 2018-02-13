include("Include.jl")
using PyCall
@pyimport matplotlib.patches as patches

# what are my colors?
# color_1 = (1/255)*[64,64,64] # dark gray -
# color_1 = (1/255)*[255,99,71] # tomato
color_1 = (1/255)*[0,191,255] # skyblue
color_2 = (1/255)*[224,224,224]
K = 0.45

# load the data_dictionary -
data_dictionary = DataDictionary(0.0,0.0,0.0)

# where do we find files?
path_to_senstivity_files = "./sensitivity"
file_pattern = "AdjSimulation-P"
time_skip = 40

# calc the sensitivity array -
# (T,SA) = calculate_sensitivity_array(path_to_senstivity_files,file_pattern,time_skip,data_dictionary)
ASSA = calculate_average_scaled_sensitivity_array(path_to_senstivity_files,file_pattern,data_dictionary)
(U,S,V) = svd(ASSA,thin=false)

# remove the "outer" data -
#raw_data_array = raw_data_array[2:end-1,2:end-1]
#(U,S,V) = svd(raw_data_array)
#
# # scale -
# scale_factor = maximum(raw_data_array)
# raw_data_array = (1.0/scale_factor)*raw_data_array

# plot the U -
scaled_data_array = abs.(V)

# what is the size of the data array?
(number_of_rows,number_of_cols) = size(scaled_data_array)

# add an extra col for colorbar -
colorbar_col = vec(transpose(linspace(0,1,number_of_rows)))

# main drawing loop -
for col_index = 1:number_of_cols

  data_scaled = scaled_data_array[:,col_index]

  @show data_scaled

  # how many patches per col?
  number_of_patches = length(data_scaled)
  epsilon = 1.2;
  ax = gca()
  for row_index = 1:number_of_patches

    origin_point = [(col_index - 1)+(col_index - 1)*epsilon + 1,(row_index - 1)+(row_index - 1)*epsilon+1];

    # compute color -

    # what color?
    fraction = (data_scaled[row_index]^2)/(K^2+data_scaled[row_index]^2)
    color_value = fraction*color_1+(1-fraction)*color_2
    circle_radius = fraction*0.9+0.1

    # draw the square -
    ax[:add_patch](

               patches.Circle(origin_point,   # (x,y)
                   radius = circle_radius,          # radius
                   facecolor=color_value,
                   edgecolor=color_value,
                   linewidth=0.5,
               )
           )
  end
end

# add the colorbar -
for col_index = 1:1

  Z = colorbar_col

  # scale Z -
  Z = abs.(Z)
  min_value = minimum(Z)
  max_value = maximum(Z)
  data_scaled = (Z-min_value)./(max_value - min_value)

  # how many patches per col?
  number_of_patches = length(data_scaled)
  epsilon = 1.2;
  ax = gca()
  for row_index = 1:number_of_patches

    origin_point = [(col_index - 1)+(col_index - 1)*epsilon+round(2.35*number_of_cols)+1,(row_index - 1)+(row_index - 1)*epsilon+1];

    # what color?
    fraction = (data_scaled[row_index]^2)/(K^2+data_scaled[row_index]^2)
    color_value = fraction*color_1+(1-fraction)*color_2
    circle_radius = fraction*0.9+0.1


    # draw the square -
    ax[:add_patch](

               patches.Circle(origin_point,   # (x,y)
                   radius = circle_radius,       # radius
                   facecolor=color_value,
                   edgecolor=color_value,
                   linewidth=0.5,
               )
           )


  end
end

axis("square")
axis("off")
savefig("./figs/ASSA-SVD-Full-Blue.pdf")
run(`pdfcrop ./figs/ASSA-SVD-Full-Blue.pdf ./figs/ASSA-SVD-Full-Blue.pdf`)
