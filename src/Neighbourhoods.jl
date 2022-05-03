abstract type AbstractNeighbourhoodWeightGenerator end

struct IDW_MooreNeighbourhoodWeight <:AbstractNeighbourhoodWeightGenerator
    radius::Int64
    scale::Int64
    excludeCentre::Bool
end
"""
    generateWeightMatrix()

Obtains the neighbourhood of cells from a raster (ras) around the coordinate
of the focal point (coord).

"""
function generateWeightMatrix(params::IDW_MooreNeighbourhoodWeight)
    w = 2*params.radius+1
    m = ones(w,w)
    mIdx = CartesianIndices(m)
    for idx in mIdx
        m[idx] = sqrt(((idx[1]-(params.radius+1))^2)+((idx[2]-(params.radius+1))^2))
    end
    m = m.^-params.scale
    if params.excludeCentre== true
        m[params.radius+1,params.radius+1] = 0.0
    end
    return(m)
end

abstract type AbstractNeighbourhood end

struct MooreNeighbours <:AbstractNeighbourhood
    r::Int64
end

# struct MooreToroidalNeighbours<:AbstractNeighbourhood
#
# end

"""
    getNeighbourhood()

Obtains the neighbourhood of cells from a raster (ras) around the coordinate
of the focal point (coord).

"""
function getNeighbourhood(nb_params::MooreNeighbours,ras::Matrix{Float64},coord::CartesianIndex,shape::Tuple)
    # Assumes square neighbourhood
    ymin = coord[1]-nb_params.r
    ymax = coord[1]+nb_params.r
    xmin = coord[2]-nb_params.r
    xmax = coord[2]+nb_params.r
    # Check boundary conditions
    if ymin < 1
        ymin = 1
    end
    if ymax > shape[1]
        ymax = shape[1]
    end
    if xmin < 1
         xmin = 1
      end
    if xmax > shape[2]
        xmax = shape[2]
    end
    return ras[ymin:ymax,xmin:xmax]
end
function getWeightedNeighbourhood(nb_params::MooreNeighbours,ras::Matrix{Float64},weights::Matrix{Float64},coord::CartesianIndex,shape::Tuple)
    # Assumes square neighbourhood
    ymin = coord[1]-nb_params.r
    ymax = coord[1]+nb_params.r
    xmin = coord[2]-nb_params.r
    xmax = coord[2]+nb_params.r
    # N
    idxYMin = 1
    idxYMax = 1+2*nb_params.r
    idxXMin = 1
    idxXMax = 1+2*nb_params.r
    # Check boundary conditions
    if ymin < 1
        ymin = 1
        idxYMin = nb_params.r+(2-coord[1])
    end
    if ymax > shape[1]
        ymax = shape[1]
        idxYMax = (nb_params.r+1)+shape[1]-coord[1]
    end
    if xmin < 1
         xmin = 1
         idxXMin = nb_params.r+(2-coord[2])
      end
    if xmax > shape[2]
        xmax = shape[2]
        idxXMax = (nb_params.r+1)+shape[2]-coord[2]
    end
    return(ras[ymin:ymax,xmin:xmax].*weights[idxYMin:idxYMax,idxXMin:idxXMax])
end
"""
    neighbourHoodWeight(ca)

Calculate the weighting of a neighbourhood

"""
function neighbourHoodWeight(paNeighbourhood::Matrix{Float64},suitNeighbourhood::Matrix{Float64})
    oc = vec(paNeighbourhood)
    su = vec(suitNeighbourhood)
    weight = oc.*su
    return(sum(weight))
end
