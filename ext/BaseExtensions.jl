module BaseExtensions



# the process envisinoed for this file turned out to be unnecessary at present, to solve VecOfVec broadcasting efficiently
# but the repeated forloops have been fit in 1 spot
# however i hope to return here for better shenanigans
"""
⊕(d, e) = dumloop_add!(d, e)
⊖(d, e) = dumloop_subtract!(d, e)
⊗(a, b) = dumloop_multiply!(d, e)
⊘(a, b) = dumloop_divide!(d, e)
  \oplus
  \ominus
  \otimes
  \oslash
cat = [[0, 0, 0], [0,0,0]]
dog = [[5, 12, 9], 1, 22, 3]]
cat ⊕ dog
cat

For use with Vectors of Mutable Vectors, 
overwrite the values of the Vector preceding the symbol with the implied arithmetic operation.


"""

⊕(d, e) = dumloop_add!(d, e)
⊖(d, e) = dumloop_subtract!(d, e)
⊗(a, b) = dumloop_multiply!(d, e)
⊘(a, b) = dumloop_divide!(d, e)


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



end #module