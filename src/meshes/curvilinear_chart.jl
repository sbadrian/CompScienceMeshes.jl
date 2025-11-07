#TODO: Make docstring
struct CurvilinearSimplex{U,D,C,N,T} <: AbstractSimplex{U,D}
    vertices::SVector{N,SVector{U,T}}
    ζnodes::NTuple{N,T}               # ← was Float64
end

# endpoints of a 1D curvilinear element (by ζ = 0 and ζ = 1)
@inline function nodes(s::CurvilinearSimplex{U,1}) where {U}
    ζ = s.ζnodes
    Tζ = eltype(ζ)
    i0 = findfirst(==(zero(Tζ)), ζ)
    i1 = findfirst(==(one(Tζ)),  ζ)
    return (s.vertices[i0], s.vertices[i1])
end

# convenience
vertices(s::CurvilinearSimplex) = s.vertices

# value-form
coordtype(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = T
dimension(::CurvilinearSimplex{U,D,C,N,T})  where {U,D,C,N,T} = D
universedimension(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = U
vertextype(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = SVector{U,T}

# type-form
coordtype(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = T
dimension(::Type{CurvilinearSimplex{U,D,C,N,T}})  where {U,D,C,N,T} = D
universedimension(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = U
vertextype(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = SVector{U,T}

function neighborhood(p::C, bary) where {C<:CurvilinearSimplex}
    D = dimension(p)
    T = coordtype(p)
    P = SVector{D,T}
    cart = barytocart(p, T.(bary))
    Q = typeof(cart)
    U = length(cart)
    MeshPointNM{T,C,D,U}(p, P(bary), cart)
end

barytocart(ch::CurvilinearSimplex{U,1}, u)   where {U} = map(ch, u[1])

cartesian(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real})  where {U} = barytocart(ch, u)
parametric(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real}) where {U} = u

#TODO: @Piotr write docstring describing what this function is doing
#TODO: can't we simply make the barytochart function out of this?
#Discuss: map in any case is not a great name. Maybe "pushforward" is better
#Because it maps from the bary coordinates to the physical space
@inline function map(ch::CurvilinearSimplex{U,1,C,N,T}, ζ::Real) where {U,C,N,T}
    z = T(ζ)
    s = zero(SVector{U,T})
    @inbounds for r in 1:N
        s += ch.vertices[r] * _ℓ(ch.ζnodes, r, z)
    end
    s
end

#TODO: @Piotr write docstring describing what this function is doing
@inline function dmap(ch::CurvilinearSimplex{U,1,C,N,T}, ζ::Real) where {U,C,N,T}
    z = T(ζ)
    s = zero(SVector{U,T})
    @inbounds for r in 1:N
        s += ch.vertices[r] * _dℓ(ch.ζnodes, r, z)
    end
    s
end

#TODO: @Piotr write docstring describing what this function is doing
# Lagrange basis on [0,1] through ζnodes
@inline function _ℓ(ζnodes::NTuple{D,T}, r::Int, ζ::T) where {D,T<:Real}
    ζr = ζnodes[r]
    num = one(T); den = one(T)
    @inbounds for j in 1:D
        j == r && continue
        num *= (ζ - ζnodes[j])
        den *= (ζr - ζnodes[j])
    end
    num/den
end

#TODO: @Piotr write docstring describing what this function is doing
@inline function _dℓ(ζnodes::NTuple{D,T}, r::Int, ζ::T) where {D,T<:Real}
    ζr = ζnodes[r]
    s = zero(T)
    @inbounds for k in 1:D
        k == r && continue
        num = one(T); den = (ζr - ζnodes[k])
        @inbounds for j in 1:D
            (j == r || j == k) && continue
            num *= (ζ - ζnodes[j])
            den *= (ζr - ζnodes[j])
        end
        s += num/den
    end
    s
end


# Gmsh-consistent ζ-nodes for a 1D line with N local nodes
gmsh_line_ζnodes(::Val{N}) where {N} =
    N == 2 ? (1.0, 0.0) :
             (1.0, 0.0, ((k)/(N-1) for k in reverse(1:N-2))...)


function chart(
    m::CompScienceMeshes.CurvilinearMesh{U,N,T,O},
    faceid::Int
) where {U,N,T,O}

    vindices = m.faces[faceid]
    X = m.vertices[vindices]
    ζnodes = gmsh_line_ζnodes(Val(N))                   # match that order

    CompScienceMeshes.CurvilinearSimplex{U,1,U-1,N,T}(X, ζnodes)
end

jacobian(ch::CurvilinearSimplex{U,1}, ζ::Real) where {U} = norm(dmap(ch, ζ))
center(ch::CurvilinearSimplex) = map(ch, 0.5)

domain(::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T} = ReferenceSimplex{1,T,2}()

function measure(ch::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T}
    ξ = (-sqrt(T(3)/T(5)), zero(T),  sqrt(T(3)/T(5)))
    w = (T(5)/T(9),       T(8)/T(9), T(5)/T(9))
    s = zero(T)
    @inbounds for k in 1:3
        ζ = (ξ[k] + one(T)) / T(2)
        s += w[k] * jacobian(ch, ζ) / T(2)
    end
    s
end