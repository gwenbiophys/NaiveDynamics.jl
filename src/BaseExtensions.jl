import Base: *, /, -, +

# the process envisinoed for this file turned out to be unnecessary at present, to solve VecOfVec broadcasting efficiently
# but the repeated forloops have been fit in 1 spot
# however i hope to return here for better shenanigans

mutable struct VecOfMVec
    vec::AbstractArray{MVector{3, AbstractFloat}, 1}
end


function dumloop_multiply!(d, e)
    @inbounds for i in eachindex(d)
        d[i] .*= e[i]
    end
end

function dumloop_subtract!(d, e)
  @inbounds for i in eachindex(d)
      d[i] .-= e[i]
  end
end
function dumloop_divide!(d, e)
  @inbounds for i in eachindex(d)
      d[i] ./= e[i]
  end
end
function dumloop_add!(d, e)
  @inbounds for i in eachindex(d)
      d[i] .+= e[i]
  end
end