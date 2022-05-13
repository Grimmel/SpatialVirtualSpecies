module SpatialVirtualSpecies
include("ColoniseBehaviour.jl")
include("Neighbourhoods.jl")
include("Utils.jl")
using Distributions
using Random
using StatsBase

abstract type AbstractSpeciesCellularAutomata end
abstract type AbstractColoniseParameters end
abstract type AbstractNeighbourhoodParameters end


abstract type AbstractScalingParameters end

mutable struct SpeciesCellularAutomata <: AbstractSpeciesCellularAutomata
    # State layer & control layer
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    survivalControl_Passive::Matrix{Float64}
    coloniseControl_Passive::Matrix{Float64}
    survivalControl_Active::Matrix{Float64}
    coloniseControl_Active::Matrix{Float64}
    SpeciesCellularAutomata(pa,pa_cart_index,survivalControl_Passive,coloniseControl_Passive) = new(pa,pa_cart_index,survivalControl_Passive,coloniseControl_Passive,survivalControl_Passive,coloniseControl_Passive)
end
mutable struct ColoniseParameters <: AbstractColoniseParameters
    posSelector::AbstractPositionSelector
    dispersalProbability::Float64
end
mutable struct ColoniseWeightedParameters <: AbstractColoniseParameters
    posSelector::AbstractPositionSelector
    dispersalProbability::Float64
    meanNumberDispersers::Int64
end
mutable struct NeighbourhoodParameters <: AbstractNeighbourhoodParameters
    neighbourhoodType::AbstractNeighbourhood
    neighbourSurvivalWeight::Float64
    neighbourColoniseWeight::Float64
    coloniseScaleParams::AbstractScalingParameters
    survivalScaleParams::AbstractScalingParameters
    weightMatrix::Matrix{Float64}
end
mutable struct LogisticScalingParameters <: AbstractScalingParameters
    maxVal::Float64
    k::Float64
    LogisticScalingParameters(;maxVal=0.9,k=-10) = new(maxVal,k)
end
"""
    scaleValue(params,suitVal,weightVal)

A function to scale the suitability values by a given weight

"""
function scaleValue(params::LogisticScalingParameters,suitVal::Float64,weightVal::Float64)
    return params.maxVal/(1+MathConstants.e^(params.k*(suitVal-(1-weightVal))))
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
    colonise(ca)

Select a single cell to colonise

"""
function colonise(ca::SpeciesCellularAutomata,cp::ColoniseParameters)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,cp.dispersalProbability)
    for i in dCells
        newXY = selectColonisedCoordinate(cp.posSelector,i,rng)
        if coordInBoundaries(newXY,shape) && (ca.coloniseControl_Active[newXY[1],newXY[2]]>0)
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
function colonise(ca::SpeciesCellularAutomata,cp::ColoniseWeightedParameters)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,cp.dispersalProbability)
    for i in dCells
        # Determine maximum number of dispersers scaled by suitabiility
        meanNumDispersers = cp.meanNumberDispersers*ca.coloniseControl_Active[i]
        if meanNumDispersers > 0.0
            numberDispersers = rand(rng,Poisson(meanNumDispersers))
            for j in 1:numberDispersers
                newXY = selectColonisedCoordinate(cp.posSelector,i,rng)
                if coordInBoundaries(newXY,shape) && ca.coloniseControl_Active[newXY[1],newXY[2]]>0
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
function colonise(ca::SpeciesCellularAutomata,cp::ColoniseWeightedParameters,np::NeighbourhoodParameters)
    rng = MersenneTwister()
    shape = size(ca.pa)
    # Select cells to be sources for colonisation
    dCells = selectProportion(ca.pa,ca.pa_cart_index,cp.dispersalProbability)
    # Random shuffle to avoid grid index bias due to order of applying this function
    idxShuffle = sample(collect(1:1:length(dCells)),length(dCells),replace=false)
    for i in idxShuffle
        cellIndex = dCells[i]
        # Get the occurrences in the neighbourhood
        paNeighbourhood = getNeighbourhood(np.neighbourhoodType,ca.pa,cellIndex,shape)
        # Get the suitabilty neighbourhood weighted by the weight matrix
        suitNeighbourhood = getWeightedNeighbourhood(np.neighbourhoodType,ca.coloniseControl_Active,np.weightMatrix,cellIndex,shape)
        # Calculate the neighbourhood weight
        neighbourWeight = neighbourHoodWeight(paNeighbourhood,suitNeighbourhood)
        colWeight = neighbourWeight * np.neighbourColoniseWeight
        # Determine maximum number of dispersers scaled by suitabiility
        coloniseMultiplier = scaleValue(np.coloniseScaleParams,ca.coloniseControl_Active[cellIndex],colWeight)
        #dispersalMultiplier = 1/(1+MathConstants.e^(-10*(ca.suitabilityActive[cellIndex]-(1-dispWeight))))
        meanNumDispersers = (1 + (ca.coloniseControl_Active[cellIndex]*coloniseMultiplier))*cp.meanNumberDispersers
        numberDispersers = rand(rng,Poisson(meanNumDispersers))
        for j in 1:numberDispersers
            newXY = selectColonisedCoordinate(cp.posSelector,cellIndex,rng)
            if coordInBoundaries(newXY,shape) && ca.coloniseControl_Active[newXY[1],newXY[2]]>0
                ca.pa[newXY[1],newXY[2]] = 1.0
            end
        end
    end
end

function extinction(ca::SpeciesCellularAutomata)
    rng = MersenneTwister()
    for idx in ca.pa_cart_index
        if ca.pa[idx] === 1.0
            survived = rand(rng,Bernoulli(ca.survivalControl_Active[idx]),1)
            if survived[1] == false
                ca.pa[idx] = 0
            end
        end
    end
end
"""
    extinction(ca)

Each occupied cell is

Note: As this is neighbour weighted, the order of cells are randomised to prevent
      any bias due to the iterating over cells in the same order each time the
      function is called. This effectively makes the neighbourhood effects
      asynchronus across the landscape.
      An alternative would be to copy the occurence matrix and refer to this for
      neighbourhood weights. This would effectively make the neighbourhood effects
      synchronus (i.e. all processes happen at the same time, instantly).
"""
function extinction(ca::SpeciesCellularAutomata,np::NeighbourhoodParameters)
    rng = MersenneTwister()
    # Randomise cell index to prevent bias from the neighbourhood weighting
    # Alternative could be to copy the the pa array.
    shape = size(ca.pa)
    nCells = shape[1]*shape[2]
    idxShuffle = sample(collect(1:1:nCells),nCells,replace=false)
    for i in idxShuffle
        cellIndex = ca.pa_cart_index[i]
        if ca.pa[i] === 1.0
            # Determine the neighbourhood weighted survival probability
            paNeighbourhood = getNeighbourhood(np.neighbourhoodType,ca.pa,cellIndex,shape)
            suitNeighbourhood = getWeightedNeighbourhood(np.neighbourhoodType,ca.survivalControl_Active,np.weightMatrix,cellIndex,shape)
            neighbourWeight = neighbourHoodWeight(paNeighbourhood,suitNeighbourhood)
            survWeight = neighbourWeight * np.neighbourSurvivalWeight
            # Scale survival probability
            sf = scaleValue(np.survivalScaleParams,ca.survivalControl_Active[cellIndex],survWeight)
            #sf = 0.9/(1+MathConstants.e^(-10*(ca.suitabilityActive[cellIndex]-(1-survWeight))))
            survivalProbability = ca.survivalControl_Active[cellIndex] + (1-ca.survivalControl_Active[cellIndex])*sf
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
    useActive::Bool
    effectSpeciesAonB::Float64
    effectSpeciesBonA::Float64


    function Interaction(strength,competetiveEdge,useActive)
        effectSpeciesAonB = strength*competetiveEdge
        effectSpeciesBonA = strength*(1-competetiveEdge)
        return(new(strength,competetiveEdge,useActive,effectSpeciesAonB,effectSpeciesBonA))
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
    if interaction.useActive == false
        speciesA.survivalControl_Active = speciesA.survivalControl_Passive .- (speciesA.survivalControl_Passive.* (interaction.effectSpeciesBonA .* speciesB.pa))
        speciesB.survivalControl_Active = speciesB.survivalControl_Passive .- (speciesB.survivalControl_Passive .* (interaction.effectSpeciesAonB .* speciesA.pa))
    else
        speciesA.survivalControl_Active = speciesA.survivalControl_Active .- (speciesA.survivalControl_Active .* (interaction.effectSpeciesBonA .* speciesB.pa))
        speciesB.survivalControl_Active = speciesB.survivalControl_Active .- (speciesB.survivalControl_Active .* (interaction.effectSpeciesAonB .* speciesA.pa))
    end
end
function reset_active_layer(species::AbstractSpeciesCellularAutomata)
    species.survivalControl_Active = species.survivalControl_Passive
end
end # module
