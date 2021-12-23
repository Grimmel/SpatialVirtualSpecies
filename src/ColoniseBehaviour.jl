using Random

abstract type AbstractPositionSelector end

struct ExponentialPosSelector <: AbstractPositionSelector
    meanDistance::Float64
end
struct UniformPosSelector<: AbstractPositionSelector
    maxDistance::Float64
end

"""
    selectColonisedCoordinate(dispParam, meanDistance,cartIdx)

Calculates a new position that does not include a starting position (cartIdx).
This is achived stochastically by sampling distance from an exponential
distribution and direction of movement from a uniform distribution between 0-360.

Note: The selected distance as a constant added to it in order to ensure a
      distance less than 1 does not select a position the same as the starting one.
      This value is +0.7071 because this is the ratio between the distance to
      a direct:diagonal neighbour.
      This ensures that rounding the X and Y values the direct neighbours are
      selected more frequently than the diagonals. Higher values will bias to
      diagonals and lower values will result in selection of the starting position.

"""
function selectColonisedCoordinate(dispParams::ExponentialPosSelector,startPos::CartesianIndex,rng::MersenneTwister)
    # Exponential version
    distance = rand(rng,Exponential(dispParams.meanDistance),1)
    # Add constant to ensure selection of pos outside current cell.
    distance = distance[1]+0.7071
    angle = 360.0*rand(rng)
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + startPos[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + startPos[1]
    coord = (y,x)
    return(coord)
end
function selectColonisedCoordinate(dispParams::UniformPosSelector,startPos::CartesianIndex,rng::MersenneTwister)
    # Uniform version
    distance = rand(rng,1)*dispParams.maxDistance
    # Add constant to ensure selection of pos outside current cell.
    distance = distance[1]+0.7071
    angle = 360.0*rand(rng)
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + startPos[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + startPos[1]
    coord = (y,x)
    return(coord)
end
