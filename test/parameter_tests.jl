include("../src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DelimitedFiles, DataFrames,CSV
function simulate(ca)
    for i in 1:200
        SpatialVirtualSpecies.colonise(ca)
        SpatialVirtualSpecies.extinction(ca)
    end
end
function testDispersalParams(suit,ls)
    results = DataFrame(disp_mean = Float64[], n_disp = Int64[],prob_disp=Float64[], iteration = Int64[],time=Float64[])

    nb_getter = SpatialVirtualSpecies.MooreNeighbours(5)
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(5,2,true)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    expoMean = [1.0,3.0,5.0]
    num_dispersers = [1,3,5]
    prob_disp = [0.1,0.2,0.4]
    for mn in expoMean
        for prob in prob_disp
            for dispersers in num_dispersers
                for i in 1:10
                    pa = ones(400,400)
                    paIdx = CartesianIndices(pa)
                    disp =SpatialVirtualSpecies.ExponentialPosSelector(mn)
                    ca = SpatialVirtualSpecies.SpeciesCellularAutomataSuitabilityWeighted(pa,paIdx,suit,suit,disp,prob,dispersers)
                    time = @elapsed simulate(ca)
                    push!(results,[mn dispersers prob i time])
                    open("F:/PhD/ca_tests/ls"*string(ls)*"/dispersal/sim_expo_mean" * string(mn) * "_prob"*string(prob)*"_disp"*string(dispersers)* "_rep"* string(i) *".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,ca.pa)
                    end
                end
            end
        end
    end
    CSV.write("F:/PhD/ca_tests/dispersal_ls"*string(ls)*"df.csv", results)
end
function testNeighbourhoodParams(suit,ls)
    results = DataFrame(nb_size = Int64[], sweight = Float64[],cweight=Float64[],iteration = Int64[],time=Float64[])
    nb_size = [1,2,3]
    surv_nb_weight = [0.05,0.1,0.3]
    col_nb_weight = [0.05,0.1,0.3]
    for neighbourhood_size in nb_size
        for surv_weight in surv_nb_weight
            for col_weight in col_nb_weight
                for i in 1:10
                    pa = ones(400,400)
                    paIdx = CartesianIndices(pa)
                    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(neighbourhood_size,1,true)
                    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
                    nb_getter = SpatialVirtualSpecies.MooreNeighbours(neighbourhood_size)
                    disp =SpatialVirtualSpecies.ExponentialPosSelector(1.0)
                    ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourhoodWeighted(pa,paIdx,suit,suit,disp,0.1,1,nb_getter,surv_weight,col_weight,weightMatrix)
                    time = @elapsed simulate(ca)
                    push!(results,[neighbourhood_size surv_weight col_weight i time])
                    open("F:/PhD/ca_tests/ls"*string(ls)*"/neighbour/sim_" * string(neighbourhood_size) * "_sweight"*string(surv_weight)*"cweight"*string(col_weight)* "_rep"* string(i) *".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,ca.pa)
                    end
                end
            end
        end
    end
    CSV.write("F:/PhD/ca_tests/neighbourhood_ls"*string(ls)*"df.csv", results)
end
function simulateInteractions(ca1,ca2,inter)
    for i in 1:200
        if i>20
            SpatialVirtualSpecies.interact(inter,ca1,ca2)
        end
        SpatialVirtualSpecies.colonise(ca1)
        SpatialVirtualSpecies.extinction(ca1)
        SpatialVirtualSpecies.colonise(ca2)
        SpatialVirtualSpecies.extinction(ca2)
    end
end
function testInteractionParams(suit1,suit2)
    results = DataFrame(strength = Float64[], edge = Float64[],iteration = Int64[],time=Float64[])
    for strength in [0.2,0.5,1.0]
        for edge in [0.0,0.3,0.6,1.0]
            for d_prob in [0.1,0.2,0.4]
                species_interactions = SpatialVirtualSpecies.Interaction(strength,edge)
                disp =SpatialVirtualSpecies.ExponentialPosSelector(3)

                for i in 1:10
                    pa_species1 = SpatialVirtualSpecies.generateStateLayer(suit1)
                    pa_species2 = SpatialVirtualSpecies.generateStateLayer(suit2)
                    paIdx = CartesianIndices(pa_species1)
                    ca_species1 = SpatialVirtualSpecies.SpeciesCellularAutomataSuitabilityWeighted(pa_species1,paIdx,suit1,suit1,disp,0.3,5)
                    ca_species2 = SpatialVirtualSpecies.SpeciesCellularAutomataSuitabilityWeighted(pa_species2,paIdx,suit2,suit2,disp,d_prob,1)

                    time = @elapsed simulateInteractions(ca_species1,ca_species2,species_interactions)
                    push!(results,[strength edge i time])
                    open("F:/PhD/ca_tests/interaction/sim3_" * string(strength) * "_edge"*string(edge)*"_disp"*string(d_prob)*"_species1_"*string(i)*".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,ca_species1.pa)
                    end
                    open("F:/PhD/ca_tests/interaction/sim3_" * string(strength) * "_edge"*string(edge)*"_disp"*string(d_prob)*"_species2_"*string(i)*".asc","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,ca_species2.pa)
                    end
                end
            end
        end
    end
    CSV.write("F:/PhD/ca_tests/interaction_df.csv", results)
end
suit_frag = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_223.asc",skipstart=6)
suit_frag_scaled = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_223scaled.asc",skipstart=6)

suit_contig= readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability_789.asc",skipstart=6)
# testDispersalParams(suit_frag,223)
testNeighbourhoodParams(suit_frag,223)
# testDispersalParams(suit_contig,789)
testNeighbourhoodParams(suit_contig,789)
# testInteractionParams(suit_frag,suit_contig)
#testDispersalParams(suit_frag_scaled,"223scaled")
testNeighbourhoodParams(suit_frag_scaled,"223scaled")
