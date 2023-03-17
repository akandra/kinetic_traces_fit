module Structs

using Parameters

# structure to deal with fitting parameters
@with_kw mutable struct fitpar

    value::Float64 = 0.0
    min::Float64   = -Inf
    max::Float64   = Inf
    var::Bool      = false
    glbl::Bool     = false

end

end