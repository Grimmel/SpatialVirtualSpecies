include("../src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DelimitedFiles
using BenchmarkTools
function simulate()
    print('hello')
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

function simulateInteractions()
    suitA = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability545.asc",skipstart=6)
    paA = ones(400,400)
    paIdxA = CartesianIndices(paA)
    suitB = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability771.asc",skipstart=6)
    paB = ones(400,400)
    paIdxB = CartesianIndices(paB)

    positionSelector = SpatialVirtualSpecies.ExponentialPosSelector(5)
    nb_getter = SpatialVirtualSpecies.MooreNeighbours(5)
    nb_weights_type = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(5,2)
    weightMatrix = SpatialVirtualSpecies.generateWeightMatrix(nb_weights_type)
    caA= SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(paA,paIdxA,suitA,suitA,positionSelector,0.3,1,nb_getter,0.1,0.1,weightMatrix)
    caB = SpatialVirtualSpecies.SpeciesCellularAutomataNeighbourEffect(paB,paIdxB,suitB,suitB,positionSelector,0.3,1,nb_getter,0.1,0.1,weightMatrix)
    inter = SpatialVirtualSpecies.Interaction(0.25,0.5)
    for i in 1:100
        if i>50
            SpatialVirtualSpecies.interact(inter,caA,caB)
        end
        SpatialVirtualSpecies.coloniseNeighbourWeight(caA)
        SpatialVirtualSpecies.extinctionNeighbourWeight(caA)
        SpatialVirtualSpecies.coloniseNeighbourWeight(caB)
        SpatialVirtualSpecies.extinctionNeighbourWeight(caB)
    end
    open("D:/testsim30A.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,caA.pa)
    end
    open("D:/testsim30B.asc","w") do io
        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
        writedlm(io,caB.pa)
    end
end
@profiler simulateInteractions()
simulate()
