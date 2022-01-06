module SpatialVirtualSpecies
include("ColoniseBehaviour.jl")
include("Neighbourhoods.jl")

using Distributions
using Random
using StatsBase

abstract type AbstractSpeciesCellularAutomata end

mutable struct SpeciesCellularAutomata <: AbstractSpeciesCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    suitabilityInitial::Matrix{Float64}
    suitabilityActive::Matrix{Float64}
    # Dispersal parameters
    posSelector::AbstractPositionSelector
    dispersalProbability::Float64
end

mutable struct SpeciesCellularAutomataSuitWeighted <: AbstractSpeciesCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    suitabilityInitial::Matrix{Float64}
    suitabilityActive::Matrix{Float64}
    # Dispersal parameters
    posSelector::AbstractPositionSelector
    dispersalProbability::Float64
    meanNumberDispersers::Int64
end

mutable struct SpeciesCellularAutomataNeighbourEffect <: AbstractSpeciesCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    suitabilityInitial::Matrix{Float64}
    suitabilityActive::Matrix{Float64}
    # Dispersal parameters
    posSelector::AbstractPositionSelector
    dispersalProbability::Float64
    meanNumberDispersers::Int64
    # Weight parameters
    neighbourhoodParams::AbstractNeighbourhood
    neighbourSurvivalWeight::Float64
    neighbourDispersalWeight::Float64
    ### TODO: Check neighbourhoodParams Radius = weightMatrix radius
    weightMatrix::Matrix{Float64}
end

"""
    selectProportion(pa,caIndex,dispersalProbability)

Given a 2D array of occurence, randomly select cells that will disperse at
the next time step. The number of dispersers is determined stochastically
by sampling from a binomial distribution.

"""
function selectProportion(pa::Matrix{Float64},caIndex::Array{CartesianIndex{2},2},probability::Float64)
    nPresences = Int(sum(pa))
    total = Array{CartesianIndex}(undef, nPresences)
    counter = 1
    for idx in caIndex
        if pa[idx] === 1.0
            total[counter] = idx
            counter+=1
        end

    end
    dispersers = sample(total,rand(Binomial(nPresences,probability)))
end
"""
    neighbourHoodWeight(ca)

Select a single cell to colonise

"""
function neighbourHoodWeight(paNeighbourhood::Matrix{Float64},suitNeighbourhood::Matrix{Float64})
    oc = vec(paNeighbourhood)
    su = vec(suitNeighbourhood)
    weight = oc.*su
    #weight = oc.*vec(weightMatrix[idxYMin:idxYMax,idxXMin:idxXMax])
    return(sum(weight))
end

#TODO
# function checkNeighbourhoodWeightDistribution(suitNeighbourhood::Matrix{Float64})

"""
    colonise(ca)

Select a single cell to colonise

"""
function colonise(ca::SpeciesCellularAutomata)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        newXY = selectColonisedCoordinate(ca.posSelector,i,rng)
        if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
            ca.pa[newXY[1],newXY[2]] = 1.0
        end
    end
end

"""
    coloniseSuitWeight(ca)

Randomly selects a number cells that attempt to colonise nearby cells.
A maximum number of dispersers is specified and each cell selected to disperse
has


"""
function coloniseSuitWeight(ca::SpeciesCellularAutomata)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        # Determine maximum number of dispersers scaled by suitabiility
        meanNumDispersers = ca.meanNumberDispersers*ca.suitabilityActive[i]
        if meanNumDispersers > 0.0
            numberDispersers = rand(rng,Poisson(meanNumDispersers))
            for j in 1:numberDispersers
                newXY = selectColonisedCoordinate(ca.posSelector,i,rng)
                if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
                    ca.pa[newXY[1],newXY[2]] = 1.0
                end
            end
        end
    end
end
"""
    coloniseNeighbourWeight(ca)

Randomly selects a number cells that attempt to colonise nearby cells.
Each cell selected to disperse has a random number of attempts at colonising
selected from a Poisson distribution. For each attempt, the new position is
selected randomly and converted to occupied (1.0).

Note: As this is neighbour weighted, the order of cells are randomised to prevent
      any bias due to the iterating over cells in the same order each time the
      function is called. This effectively makes the neighbourhood effects
      asynchronus across the landscape.
      An alternative would be to copy the occurence matrix and refer to this for
      neighbourhood weights. This would effectively make the neighbourhood effects
      synchronus (i.e. all processes happen at the same time, instantly).
"""
function coloniseNeighbourWeight(ca::SpeciesCellularAutomataNeighbourEffect)
    rng = MersenneTwister()
    shp = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    # Random shuffle to avoid grid index bias due to order of applying this function
    idxShuffle = sample(collect(1:1:length(dCells)),length(dCells),replace=false)
    for i in idxShuffle
        cellIndex = dCells[i]
        paNeighbourhood = getNeighbourhood(ca.neighbourhoodParams,ca.pa,cellIndex,shp)
        suitNeighbourhood = getWeightedNeighbourhood(ca.neighbourhoodParams,ca.suitabilityActive,ca.weightMatrix,cellIndex,shp)
        neighbourWeight = neighbourHoodWeight(paNeighbourhood,suitNeighbourhood)
        dispWeight = neighbourWeight * ca.neighbourDispersalWeight
        # Determine maximum number of dispersers scaled by suitabiility
        dispersalMultiplier = 1/(1+MathConstants.e^(-10*(ca.suitabilityActive[cellIndex]-(1-dispWeight))))
        meanNumDispersers = (1 + (ca.suitabilityActive[cellIndex]*dispersalMultiplier))*ca.meanNumberDispersers

        numberDispersers = rand(rng,Poisson(meanNumDispersers))
        for j in 1:numberDispersers
            newXY = selectColonisedCoordinate(ca.posSelector,cellIndex,rng)
            if newXY[2]>=1 && newXY[2] <= shp[2] && newXY[1] >=1 && newXY[1]<=shp[1]
                ca.pa[newXY[1],newXY[2]] = 1.0
            end
        end
    end
end


function extinction(ca::SpeciesCellularAutomata)
    rng = MersenneTwister()
    for idx in ca.pa_cart_index
        if ca.pa[idx] === 1.0
            survived = rand(rng,Bernoulli(ca.suitabilityActive[idx]),1)
            if survived[1] == false
                ca.pa[idx] = 0
            end
        end
    end
end
"""
    extinctionNeighbourWeight(ca)

Each occupied cell is

Note: As this is neighbour weighted, the order of cells are randomised to prevent
      any bias due to the iterating over cells in the same order each time the
      function is called. This effectively makes the neighbourhood effects
      asynchronus across the landscape.
      An alternative would be to copy the occurence matrix and refer to this for
      neighbourhood weights. This would effectively make the neighbourhood effects
      synchronus (i.e. all processes happen at the same time, instantly).
"""
function extinctionNeighbourWeight(ca::SpeciesCellularAutomataNeighbourEffect)
    rng = MersenneTwister()
    # Randomise cell index to prevent bias from the neighbourhood weighting
    # Alternative could be to copy the the pa array.
    shp = size(ca.pa)
    nCells = shp[1]*shp[2]
    idxShuffle = sample(collect(1:1:nCells),nCells,replace=false)
    for i in idxShuffle
        cellIndex = ca.pa_cart_index[i]
        if ca.pa[i] === 1.0
            # Determine the neighbourhood weighted survival probability
            paNeighbourhood = getNeighbourhood(ca.neighbourhoodParams,ca.pa,cellIndex,shp)
            suitNeighbourhood = getWeightedNeighbourhood(ca.neighbourhoodParams,ca.suitabilityActive,ca.weightMatrix,cellIndex,shp)
            neighbourWeight = neighbourHoodWeight(paNeighbourhood,suitNeighbourhood,ca.weightMatrix)
            survWeight = neighbourWeight * ca.neighbourSurvivalWeight
            # Logistic function to scale survival proba
            sf = 1/(1+MathConstants.e^(-10*(ca.suitabilityActive[cellIndex]-(1-survWeight))))
            survivalProbability = ca.suitabilityActive[cellIndex] + (1-ca.suitabilityActive[cellIndex])*sf
            # Determine survival
            survived = sample(rng,[0,1],ProbabilityWeights([1.0-survivalProbability,survivalProbability]))
            if survived === 0
                ca.pa[cellIndex] = 0
            end
        end
    end
end

struct Interaction
    strength::Float64
    competetiveEdge::Float64
    effectSpeciesAonB::Float64
    effectSpeciesBonA::Float64

    function Interaction(strength,competetiveEdge)
        effectSpeciesAonB = strength*competetiveEdge
        effectSpeciesBonA = strength*(1-competetiveEdge)
        return(new(strength,competetiveEdge,effectSpeciesAonB,effectSpeciesBonA))
    end
end
"""
    interact(interaction)

Interact is a function that defines interaction behaviour by modiying suitability
of a cell, depending on if it is colonised by the interacting species. It uses the
initial suitability to ensure the values revert back to the original value if the
interacting species no long occupies the cell in subsequent iterations.

"""
function interact(interaction::Interaction,speciesA::AbstractSpeciesCellularAutomata,speciesB::AbstractSpeciesCellularAutomata)
    speciesA.suitabilityActive = speciesA.suitabilityInitial .* (1 .- (interaction.effectSpeciesBonA .* speciesB.pa))
    speciesB.suitabilityActive = speciesB.suitabilityInitial .* (1 .- (interaction.effectSpeciesAonB .* speciesA.pa))
end
end # module
