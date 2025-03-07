@doc raw"""
    abstract type Semifield end

Abstract type for a semifield ``(S, \oplus, \otimes, inv, zero, one, \top)``.
"""
abstract type Semifield{T} end

"""
    val(x::Semifield) → T

Return the "real" value / object wrapped in the semifield type
"""
val(x::Semifield) = x.val

"""
    Base.valtype(::Type{<:Semifield{T}}) where T

Return the type of the value wrapped by the semifield.
"""
Base.valtype(::Type{<:Semifield{T}}) where T = T


"""
    x ⊕ y

Semifield addition.
"""
⊕

"""
    x ⊗ y

Semifield multiplication.
"""
⊗

"""

Semifield inversion.
"""
inv

"""
Top element for a completed semifield ''\top'' (get with \top)
"""
⊤

## FVA: The following seems to be incorrect except on certain semirings/semifields. Mmp is not one of them. 

"""
        Overload the ''Base.zero'' unary constructor to work as well as on types. 
"""
Base.zero(x::Semifield) = zero(typeof(x))

"""
        Overload the ''Base.one'' unary constructor to work on elements as well as on types. 
"""
Base.one(x::Semifield) = one(typeof(x))

"""
        Overload the ''top'' unary constructor to work on elements as well as on types. 
"""
top(x::Semifield) = top(typeof(x))

"""
Scalar product by an integer as repeated addition.

This acts as left and side actions for the semifield as a semimodule.

This operation transforms the semifield into a semimodule, so it should be generalized a bit.  
"""
function Base.:*(i::Integer, s::Semifield)
    i < 0 && throw(ArgumentError("integer has to be positive or zero"))
    iszero(i) && return zero(s)

    res = zero(s)
    for n in 1:i
        res = res ⊕ s
    end
    res
end
Base.:*(s::Semifield, i::Integer) = i * s
#FVA: but this limits the algebra to be that of integer multiplication.
#Also, this is only true in the commutative core of the Semiring, cfr. Golan on semimodules.

"""
Basic convert function for Semifields elements.
"""
Base.convert(T::Type{<:Semifield}, x::Number) = T(x)
Base.convert(T::Type{<:Semifield}, x::Semifield) = T(val(x))

"""
Basic inline equality predicate for semifield elements.
"""
Base.:(==)(x::Semifield, y::Semifield) = val(x) == val(y)

"""
Basic equality predicate for semifield elements.
"""
Base.isequal(x::Semifield, y::Semifield) = isequal(val(x), val(y))

"""
Basic approximate equality predicate for semifield elements.
"""
Base.isapprox(x::Semifield, y::Semifield) = val(x) ≈ val(y)

"""
Basic output primitive for semifield elements.
"""
Base.show(io::IO, x::Semifield) = print(io, val(x))

