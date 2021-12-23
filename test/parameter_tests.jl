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
    pa = ones(400,400)
    paIdx = CartesianIndices(pa)
    nb_getter = SpatialVirtualSpecies.MooreNeighbours(5)
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(5,2)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    expoMean = [1,3,5,10]
    uniMax = [4,5,10,15]
    for i in 1:2
        if i == 1
            for mn in expoM
                disp =SpatialVirtualSpecies.ExponentialPosSelector(mn)
                ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(pa,paIdx,suit,suit,positionSelector,0.3,1,nb_getter,0.1,0.1,weightMatrix)
                simulate(ca)
                open("D:/testsim30.asc","w") do io
                    write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                    writedlm(io,ca.pa)
                end

            end
        elseif i == 2
            for mx in uniMax
                disp =SpatialVirtualSpecies.UniformPosSelector(mx)
            end
        end
    end
end
function testNeighbourhoodParams()

end
function testInteractionParams()

end





function simulateWeighted()
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability545.asc",skipstart=6)
    pa = ones(400,400)
    paIdx = CartesianIndices(pa)
    positionSelector = SpatialVirtualSpecies.ExponentialPosSelector(5)
    nb_getter = SpatialVirtualSpecies.MooreNeighbours(5)
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(5,2)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    ca = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(pa,paIdx,suit,suit,positionSelector,0.3,1,nb_getter,0.1,0.1,weightMatrix)
    for i in 1:100
        SpatialVirtualSpecies.coloniseNeighbourWeight(ca)
        SpatialVirtualSpecies.extinctionNeighbourWeight(ca)
    end
    open("D:/testsim30.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,ca.pa)
    end
end
