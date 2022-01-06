include("../src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DelimitedFiles
function simulate(ca)
    for i in 1:100
        SpatialVirtualSpecies.coloniseNeighbourWeight(ca)
        SpatialVirtualSpecies.extinctionNeighbourWeight(ca)
    end
end
function testDispersalParams(suit)

    nb_getter = SpatialVirtualSpecies.MooreNeighbours(5)
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(5,2,true)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    expoMean = [1.0,3.0,5.0,10.0]
    uniMax = [4.0,5.0,10.0,15.0]
    num_dispersers = [1,3,5]
    prob_disp = [0.1,0.3,0.5]
    for i in 1:2
        if i == 1
            for mn in expoMean
                for prob in prob_disp
                    for dispersers in num_dispersers
                        for i in 1:10
                            pa = ones(400,400)
                            paIdx = CartesianIndices(pa)
                            disp =SpatialVirtualSpecies.ExponentialPosSelector(mn)
                            ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(pa,paIdx,suit,suit,disp,prob,dispersers,nb_getter,0.1,0.1,weightMatrix)
                            simulate(ca)
                            open("F:/PhD/ca_tests/dispersal/sim_expo_mean" * string(mn) * "_prob"*string(prob)*"_disp"*string(dispersers)* "_rep"* string(i) *".asc","w") do io
                                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                                writedlm(io,ca.pa)
                            end
                        end
                    end
                end
            end
        elseif i == 2
            for mx in uniMax
                for prob in prob_disp
                    for dispersers in num_dispersers
                        for i in 1:10
                            pa = ones(400,400)
                            paIdx = CartesianIndices(pa)
                            disp =SpatialVirtualSpecies.UniformPosSelector(mx)
                            ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(pa,paIdx,suit,suit,disp,prob,dispersers,nb_getter,0.1,0.1,weightMatrix)
                            simulate(ca)
                            open("F:/PhD/ca_tests/dispersal/sim_uniform_max" * string(mx) * "_prob"*string(prob)*"_disp"*string(dispersers)* "_rep" * string(i) *".asc","w") do io
                                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                                writedlm(io,ca.pa)
                            end
                        end
                    end
                end
            end
        end
    end
end
function testNeighbourhoodParams(suit)
    nb_size = [1,3,5,7,9]
    nb_weight = [0.0,0.3,0.6,1.0]

    weight_scale = [1,2,3]
    for neighbourhood_size in nb_size
        for neighbourhood_weight in nb_weight
            for scale in weight_scale
                for i in 1:10
                    pa = ones(400,400)
                    paIdx = CartesianIndices(pa)
                    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,scale,true)
                    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
                    nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
                    disp =SpatialVirtualSpecies.ExponentialPosSelector(neighbourhood_size)
                    ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(pa,paIdx,suit,suit,disp,0.2,1,nb_getter,neighbourhood_weight,neighbourhood_weight,weightMatrix)
                    simulate(ca)
                    open("F:/PhD/ca_tests/neighbour/sim_" * string(neighbourhood_size) * "_nbweight"*string(neighbourhood_weight)*"_scale"*string(scale)* "_rep"* string(i) *".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,ca.pa)
                    end
                end
            end
        end
    end
end
# function testInteractionParams()
#
# end
suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability545.asc",skipstart=6)

#testDispersalParams(suit)
testNeighbourhoodParams(suit)
