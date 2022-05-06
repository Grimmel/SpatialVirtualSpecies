include("../src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DelimitedFiles, DataFrames,CSV


function species_1_stabilityTest(suit::Matrix{Float64})
    results = DataFrame(rep = Int64[], timestep = Int64[], nOcc = Float64[])

    neighbourhood_weight = 0.05
    neighbourhood_size = 2
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,1,true)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
    disp =SpatialVirtualSpecies.ExponentialPosSelector(5.0)
    for i in 1:1
        pa = SpatialVirtualSpecies.generateStateLayer(suit)
        paIdx = CartesianIndices(pa)
        ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa,paIdx,suit,suit,disp,0.25,1,nb_getter,neighbourhood_weight,neighbourhood_weight,weightMatrix)
        for j in 1:5000
            SpatialVirtualSpecies.colonise(ca)
            SpatialVirtualSpecies.extinction(ca)
            if j > 10
                nOccupied = sum(ca.pa)
                push!(results,[i j nOccupied])
            end
        end
        open("F:/PhD/ca_tests/frag_05_5_25.asc","w") do io
            write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
            writedlm(io,ca.pa)
        end
    end
    #CSV.write("F:/PhD/transfer_species/species1_v7_frag.csv", results)
end
function species_1_simulation_set(suit::Matrix{Float64},ls_name)
        neighbourhood_weight = 0.05
        neighbourhood_size = 2
        nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,1,true)
        weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
        nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
        disp =SpatialVirtualSpecies.ExponentialPosSelector(5.0)
        for i in 201:500
            pa = SpatialVirtualSpecies.generateStateLayer(suit)
            paIdx = CartesianIndices(pa)
            ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa,paIdx,suit,suit,disp,0.25,1,nb_getter,neighbourhood_weight,neighbourhood_weight,weightMatrix)
            for j in 1:3000
                SpatialVirtualSpecies.colonise(ca)
                SpatialVirtualSpecies.extinction(ca)
            end
            open("F:/PhD/transfer_species/"*ls_name*"/no_interaction/occurrences"*string(i)*".asc","w") do io
                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                writedlm(io,ca.pa)
            end
        end
end
function interaction_simulation_set(suit1::Matrix{Float64},suit2::Matrix{Float64})
        neighbourhood_weight = 0.05
        neighbourhood_size = 2
        nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,1,true)
        weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
        nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
        disp =SpatialVirtualSpecies.ExponentialPosSelector(5.0)
        scale = SpatialVirtualSpecies.LogisticScalingParameters()

        species_interactions = SpatialVirtualSpecies.Interaction(0.7,0.2)

        for i in 201:500
            pa1 = SpatialVirtualSpecies.generateStateLayer(suit1)
            pa1Idx = CartesianIndices(pa1)
            ca1 = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa1,pa1Idx,suit1,suit1,disp,0.25,1,nb_getter,neighbourhood_weight,neighbourhood_weight,scale,scale,weightMatrix)

            pa2 = SpatialVirtualSpecies.generateStateLayer(suit2)
            pa2Idx = CartesianIndices(pa2)
            ca2 = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa2,pa2Idx,suit2,suit2,disp,0.1,1,nb_getter,0.05,0.05,scale,scale,weightMatrix)
            for j in 1:4000
                if j>20
                    SpatialVirtualSpecies.interact(species_interactions,ca1,ca2)
                end
                SpatialVirtualSpecies.colonise(ca1)
                SpatialVirtualSpecies.extinction(ca1)
                SpatialVirtualSpecies.colonise(ca2)
                SpatialVirtualSpecies.extinction(ca2)
            end
            open("F:/PhD/transfer_species/landscape_1/interaction/occurrences"*string(i)*".asc","w") do io
                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                writedlm(io,ca1.pa)
            end
        end
end
function species_int_stabilityTest(suit1::Matrix{Float64},suit2::Matrix{Float64})
    results = DataFrame(rep = Int64[], timestep = Int64[], nOcc = Float64[])

    neighbourhood_weight = 0.05
    neighbourhood_size = 2
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,1,true)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
    disp =SpatialVirtualSpecies.ExponentialPosSelector(5.0)

    species_interactions = SpatialVirtualSpecies.Interaction(0.7,0.2)

    for i in 1:1
        pa1 = SpatialVirtualSpecies.generateStateLayer(suit1)
        pa1Idx = CartesianIndices(pa1)
        ca1 = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa1,pa1Idx,suit1,suit1,disp,0.25,1,nb_getter,neighbourhood_weight,neighbourhood_weight,weightMatrix)

        pa2 = SpatialVirtualSpecies.generateStateLayer(suit2)
        pa2Idx = CartesianIndices(pa2)
        ca2 = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa2,pa2Idx,suit2,suit2,disp,0.1,1,nb_getter,0.05,0.05,weightMatrix)

        for j in 1:5000
            if j>20
                SpatialVirtualSpecies.interact(species_interactions,ca1,ca2)
            end
            SpatialVirtualSpecies.colonise(ca1)
            SpatialVirtualSpecies.extinction(ca1)
            SpatialVirtualSpecies.colonise(ca2)
            SpatialVirtualSpecies.extinction(ca2)
            if j > 10
                nOccupied = sum(ca1.pa)
                push!(results,[i j nOccupied])
            end
        end
        open("F:/PhD/ca_tests/interaction_sp1v5.asc","w") do io
            write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
            writedlm(io,ca1.pa)
        end
        open("F:/PhD/ca_tests/interaction_sp2v5.asc","w") do io
            write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
            writedlm(io,ca2.pa)
        end
    end
    #CSV.write("F:/PhD/transfer_species/species1_int_v1.csv", results)

end
suit_frag = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_223.asc",skipstart=6)
suit_cont = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_789.asc",skipstart=6)
suit_int =  readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_32.asc",skipstart=6)
#species_1_stabilityTest(suit_frag)
#species_1_stabilityTest(suit_cont)
#species_int_stabilityTest(suit_cont,suit_int)
#species_1_simulation_set(suit_cont,"landscape_1")
interaction_simulation_set(suit_cont,suit_int)
