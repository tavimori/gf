module gf

using LinearAlgebra

include("DivTable8.jl")
include("MulTable8.jl")

export GF, rank, rref_with_pivots!

struct GF <: Number
    x::Int64
    function GF(x::Integer)
        if x > 255 || x < 0
            error("the field only supports 0-255")
        end
        return new(x)
    end
end

# unary plus
Base.:(+)(a::GF) = a
# unary minus 
Base.:-(a::GF) = a
# binary plus
Base.:(+)(a::GF, b::GF)::GF = GF(a.x ⊻ b.x)
# binary minus
Base.:-(a::GF, b::GF)::GF = GF(a.x ⊻ b.x)

# abs
Base.:abs(a::GF)::Integer = a == GF(0) ? 0 : 1

# # FIXME: this is a trivial norm, perhaps should implement a general p-norm? 
LinearAlgebra.:norm(a::GF) = a == GF(0) ? 0 : 1

function LinearAlgebra.:rank(A::AbstractArray{GF,2}) 
    isempty(A) && return 0 # 0-dimensional case
    # s = svdvals(A)
    # tol = max(atol, rtol*s[1])
    # count(x -> x > tol, s)
    B, p = rref_with_pivots!(A)
    return length(p)
end

# # binary multiplication
Base.:*(a::GF, b::GF)::GF = GF(mulTable8[a.x << 8 | b.x + 1])
Base.:/(a::GF, b::GF)::GF = GF(divTable8[a.x << 8 | b.x + 1])

# # modified from https://github.com/blegat/RowEchelon.jl/blob/master/src/RowEchelon_with_pivots.jl
function rref_with_pivots!(A::Matrix{T}, ɛ=T <: Union{Rational,Integer} ? 0 : eps(norm(A,Inf))) where T
    nr, nc = size(A)
    pivots = Vector{Int64}()
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i:nr,j]))
        mi = mi+i - 1
        if m <= ɛ
            if ɛ > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            for k=j:nc
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j]
            for k = j:nc
                A[i,k] /= d
            end
            for k = 1:nr
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            append!(pivots,j)
            i += 1
            j += 1
        end
    end
    return A, pivots
end


end # module
