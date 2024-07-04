#module BaseExtensions



# the process envisinoed for this file turned out to be unnecessary at present, to solve VecOfVec broadcasting efficiently
# but the repeated forloops have been fit in 1 spot
# however i hope to return here for better shenanigans
export 
  dumloop_add!,
  dumloop_subtract!,
  dumloop_multiply!,
  dumloop_divide!,
  dumloop_product

function dumloop_add!(d::Vec3D, e::Vec3D)
  for i in eachindex(d)
      d[i] .+= e[i]
  end
  return d
end
function dumloop_add!(d::StatVec3D, e::StatVec3D)
  for i in eachindex(d)
      d[i] += e[i]
  end
  return d
end

function dumloop_subtract!(d::Vec3D, e::Vec3D)
  for i in eachindex(d)
      d[i] .-= e[i]
  end
  return d
end
function dumloop_subtract!(d::StatVec3D, e::StatVec3D)
  for i in eachindex(d)
      d[i] -= e[i]
  end
  return d
end

function dumloop_multiply!(d::Vec3D, e::Vec3D)
  for i in eachindex(d)
      d[i] .*= e[i]
  end
  return d
end
function dumloop_multiply!(d::StatVec3D, e::StatVec3D)
  for i in eachindex(d)
      d[i] .*= e[i]
  end
  return d
end
function dumloop_multiply!(d::Vec3D, e::Vector{T}) where T
  for i in eachindex(d)
      d[i] .*= e[i]
  end
  return d
end
function dumloop_multiply!(d::StatVec3D, e::Vector{T}) where T
  for i in eachindex(d)
      d[i] *= e[i]
  end
  return d
end

function dumloop_product(result::StatVec3D, d::StatVec3D, e::StatVec3D )
  for i in eachindex(d)
      result[i] = d[i] .* e[i]
  end
  return result
end
function dumloop_product(result::Vec3D, d::Vec3D, e::Vec3D )
  result = d
  for i in eachindex(d)
    result[i] .*=  e[i]
  end

  return result
end
function dumloop_multiply!(d::Vec3D, e::Number)
  for i in eachindex(d)
      d[i] .*= e
  end
  return d
end
function dumloop_multiply!(d::StatVec3D, e::Number)
  for i in eachindex(d)
      d[i] *= e
  end
  return d
end


function dumloop_divide!(d::Vec3D, e::Vec3D)
  for i in eachindex(d)
      d[i] ./= e[i]
  end
  return d
end

function dumloop_divide!(d::Vec3D, e::Vector{T}) where T
  for i in eachindex(d)
      d[i] ./= e[i]
  end
end
function dumloop_divide!(d::StatVec3D, e::StatVec3D)
  for i in eachindex(d)
      d[i] ./= e[i]
  end
  return d
end
function dumloop_divide!(d::StatVec3D, e::Number)
  for i in eachindex(d)
      d[i] /= e[i]
  end
  return d
end
function dumloop_divide!(d::StatVec3D, e::Vector{T}) where T
  for i in eachindex(d)
      d[i] /= e[i]
  end
  return d
end





#end #module