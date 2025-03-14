# SPDX-License-Identifier: CECILL-2.1
#FVA: Bring in all of the machinery from Semirings.jl
#FVA: However we define the matrix semiring on dense matrices. 
import LinearAlgebra: diagm
#using Reexport
#@reexport 
using Semirings

"""
    MatrixSemiring{D,T} <: Semiring{T}

Semiring over `D x D` square matrices. `+` is the elementwise addition and `*`
is the standard matrix multiplication.
"""
struct MatrixSemiring{D,T} #<: Semiring{AbstractMatrix{T}}
    #val::T#FVA: this looks buggy. It should be a matrix.
    val::AbstractMatrix{T}#No indication that they have to be a semiring here!

    function MatrixSemiring{D,T}(x::AbstractMatrix{T}) where {D,T}
        size(x, 1) == size(x, 2) || throw(ArgumentError("`x` should be a square matrix"))
        new{D,T}(x)
    end

    MatrixSemiring{D}(x::T) where {D,T} = MatrixSemiring{D,T}(x)
end
#FVA: the following constructor subtypes the Semiring for a specific dimension.
#MatrixSemiring(x::AbstractMatrix) = MatrixSemiring{size(x, 1),typeof(x)}(x)
MatrixSemiring(x::AbstractMatrix) = MatrixSemiring{size(x, 1),eltype(x)}(x)

Semirings.val(x::MatrixSemiring) = x.val
Base.:+(x::MatrixSemiring{D}, y::MatrixSemiring{D}) where D = MatrixSemiring{D}(val(x) + val(y))
Base.:*(x::MatrixSemiring{D}, y::MatrixSemiring{D}) where D = MatrixSemiring{D}(val(x) * val(y))

Base.zero(::Type{MatrixSemiring{D,T}}) where {D,T} = MatrixSemiring{D}(zeros(eltype(T), D, D))
Base.one(::Type{MatrixSemiring{D,T}}) where {D,T} = MatrixSemiring{D,T}(diagm(ones(eltype(T), D)))
#Base.zero(::Type{MatrixSemiring{D,T}}) where {D,T<:AbstractMatrix} = MatrixSemiring{D,T}(spzeros(eltype(T), D, D))
#Base.one(::Type{MatrixSemiring{D,T}}) where {D,T<:AbstractMatrix} = MatrixSemiring{D,T}(sparse([(d, d) for d in 1:D], one(eltype(T)), D, D))

#Base.iszero(x::MatrixSemiring{D,T}) where {D,T<:AbstractMatrix} = nnz(x) == 0 || all(iszero, x.values)
Base.iszero(x::MatrixSemiring) = all(iszero, x)
Base.isone(x::MatrixSemiring{D}) where D = all(d -> isone(x[d,d]), 1:D)

# Direct sum/product
#@inline ⊕(x::MatrixSemiring{P}, y::MatrixSemiring{Q}) where {P,Q} = MatrixSemiring{P+Q}(cat(val(x), val(y); dims = (1, 2)))
#@inline ⊗(x::MatrixSemiring{P}, y::MatrixSemiring{Q}) where {P,Q} = MatrixSemiring{P*Q}(kron(val(x), val(y)))
#=
function Base.Broadcast.broadcasted(::typeof(⊕),
                                    x::SparseArray{Tx,I,N},
                                    y::SparseArray{Ty,I,N}) where {Tx<:MatrixSemiring,Ty<:MatrixSemiring,I,N}
    size(x) == size(y) || throw(DimensionMismatch("`x` and `y` should have the same size"))
    nzcoo = vcat(x.nzcoo, y.nzcoo)
    nzval = vcat(x.val .⊕ Ref(zero(Ty)), Ref(zero(Tx)) .⊕ y.nzval)
    SparseArray{T,I,N}(size(x), nzcoo, nzval)
end
=#
Base.transpose(x::MatrixSemiring{D,T}) where {D,T} = MatrixSemiring{D,T}(copy(transpose(val(x))))

#Base.Broadcast.broadcasted(::typeof(Base.transpose), x::AbstractSparseArray{T,I,N}) where {T<:MatrixSemiring,I,N} =
#    SparseArray{T,I,N}(size(x), x.nzcoo, transpose.(x.nzval))
Base.Broadcast.broadcasted(::typeof(Base.transpose), x::MatrixSemiring{D,T}) where {D,T} =
    MatrixSemiring{D,T}(transpose(val(x)))#FVA: use full-typed constructor.
