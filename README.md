# EndOSim

EndOSim is a tool for simulating virtual species occurrence maps for use in species distribution modelling experiments. Often, these types of experiments use simulated species distributions that are obtained by defining a habitat suitabilty map with values ranging between zero and one, interpreted as a probability of occurrence. Bernoulli trials are applied to each cell in the map to determine if the cell is occupied or not.

The limitation of this is that it does not incorporate spatially explicit processes that can influence species distribution characteristics in the real world. This tool gives users a simple way of incorporating movement, neighbourhood effects and species interactions to create more realistic virtual species distributions.

## Setting up a simulation

To begin setting up a simulation, we need to load in a habitat suitability raster in ascii format. It is also important that when loading the raster that we also use the getHeader() function to save the file header. This is important for exporting the virtual species distributions to file after running the simulations.

The header for the file used here is 6 lines in length. This tells the getHeader() function how many lines to read and also how many lines to skip when reading the raster.


```julia
include("D:/git/SpatialVirtualSpecies/src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DataFrames
using PyPlot

header = SpatialVirtualSpecies.getHeader("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability789.asc",6)
suitability = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability789.asc",skipstart=6)

imshow(suitability)
colorbar()
```


![png](img/readme/output_1_0.png)



### Initialise the state layer (species occurrences)

Now we have our habitat suitability layer, we also need to initialise a state layer. The state layer is separate 2D matrix that can be either a zero or a one. If the state of a cell is set to one, then it is occupied by the virtual species, if it is zero then it is unoccupied.

In this example I have randomly set 1% of cells that have suitability values between 0.4-1.0 as occupied. 


```julia
pa = SpatialVirtualSpecies.generateStateLayer(suitability,0.01,(0.4,1.0))
imshow(pa,cmap="gray")
```


![png](img/readme/output_3_0.png)





### Define simulation parameters

Next we define a set of parameters for our simulation. To being with, we'll just use the simplest simulation and focus only on dispersal parameters. To do this, we define three parameters that control the mean number of dispersers from a cell, the probability that a cell will be selected for dispersal and the mean dispersal distance (measured in cells units).

To control the behaviour of how cells are selected for colonisation, a position selector (POS) constructor is created. In this example we will use an exponential function but a linear function is available or you can create your own (See XXX). 


```julia
# Create cartesian index 
paIdx = CartesianIndices(suitability)
# Define dispersal parameters
num_dispersers = 1
prob_dispersal = 0.2
dispersal_dist = 3.0
# Define POS parameter constructors
pos_params = SpatialVirtualSpecies.ExponentialPosSelector(dispersal_dist)
#Create simulation constructor
ca = SpatialVirtualSpecies.SpeciesCellularAutomataSuitabilityWeighted(pa,paIdx,suitability,suitability,pos_params,prob_dispersal,num_dispersers)
```





Both the habitat suitabilty and state layers can be accessed through the cellular automata object (ca)


```julia
imshow(ca.suitabilityActive)
imshow(ca.pa,cmap="gray")
```


![png](img/readme/output_7_0.png)




### Create a simulation function

With our parameters set up, we now need to define the simulation behaviour. For this example we will keep things simple and ensure that for each iteration we apply the colonise function and then the extinction function to our state layers.

You can modify the order of these functions or define your own function that modifies the suitability layer (e.g. temporal trends) or the state layer (e.g. disturbances)


```julia
function simulate(ca,iterations)
    for i in 1:iterations
        SpatialVirtualSpecies.colonise(ca)
        SpatialVirtualSpecies.extinction(ca)
    end
end  
```




### Run the simulation

Now let's run the simulation for 200 timesteps and see what our virtual species distribution looks like


```julia
simulate(ca,200)
imshow(ca.pa,cmap="gray")
```


![png](img/readme/output_11_0.png)





### Write the outputs to file

Using the header we pulled out of the raster file at the start, we can then write this back into an ascii file.


```julia
open("F:/test_sim.asc","w") do io
    write(io,header)
    writedlm(io,ca.pa)
end
```
