
function generateStateLayer(suitability::Matrix{Float64})
    shp = size(suitability)
    return(ones(shp[1],shp[2]))
end

function generateStateLayer(suitability::Matrix{Float64},proportion::Float64,suitabilityRange::Tuple{Float64,Float64} )
    shp = size(suitability)
    idx = CartesianIndices(suitability)
    validIdx = CartesianIndex[]
    state = zeros(shp[1],shp[2])
    for i in idx
        if suitability[i]>=suitabilityRange[1] && suitability[i]<=suitabilityRange[2]
            push!(validIdx,i)
        end
    end
    selectedIdx = sample(1:length(validIdx),Int(round(proportion*length(validIdx),digits=0)))
    selected = validIdx[selectedIdx]
    for i in selected
        state[i] = 1
    end
    return state
end
