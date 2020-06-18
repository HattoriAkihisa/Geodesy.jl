abstract type CovMatrix end
abstract type CovMatrixRotation end

function Base.show(io::IO, cm::CovMatrix)
    tn = cm |> typeof |> string
    fn = cm |> typeof |> fieldnames .|> string
    print("$tn($(fn[1])=$(getfield(cm,1)), $(fn[2])=$(getfield(cm,2)), $(fn[3])=$(getfield(cm,3)),
                $(fn[4])=$(getfield(cm,4)), $(fn[5])=$(getfield(cm,5)), $(fn[6])=$(getfield(cm,6)))")
end

function Base.show(io::IO, cmr::CovMatrixRotation)
    print(io, typeof(cmr))
    Base.show(io, getfield(cmr, :R))
end

struct CovECEF{T <: Number} <: CovMatrix
    σ2x::T
    σ2y::T
    σ2z::T
    σxy::T
    σyz::T
    σzx::T

    function CovECEF(args...)
        T = promote_type(typeof.(args)...)
        new{T}(args...)
    end
end

struct CovENU{T <: Number} <: CovMatrix
    σ2e::T
    σ2n::T
    σ2u::T
    σen::T
    σnu::T
    σue::T

    function CovENU(args...)
        T = promote_type(typeof.(args)...)
        new{T}(args...)
    end
end

"""
    matcov(T, Σ)

Returen CovMatrix struct
"""
@inline function matcov(T, Σ)
    return T(Σ[1,1], Σ[2,2], Σ[3,3], Σ[2,1], Σ[3,2], Σ[3,1])
end

"""
    covmat(cm)

Return variance-covariance matrix
"""
@inline function covmat(cm::T) where T <: CovMatrix
    [
        getfield(cm, 1) getfield(cm, 4) getfield(cm, 6)
        getfield(cm, 4) getfield(cm, 2) getfield(cm, 5)
        getfield(cm, 6) getfield(cm, 5) getfield(cm, 3)
    ]
end

struct CovENUfromECEF{T} <: CovMatrixRotation
    R::SMatrix{3,3,T}
    Rt::SMatrix{3,3,T}
end

CovENUfromECEF(R::M) where M <: AbstractMatrix = CovENUfromECEF(SMatrix{3,3}(R), SMatrix{3,3}(transpose(R)))
CovENUfromECEF(trans::ENUfromECEF) = CovENUfromECEF(trans.lat, trans.lon)

CovENUfromECEF(origin::ECEF{T}, lat::T, lon::T) where T <: Number = CovENUfromECEF(lat, lon)
CovENUfromECEF(origin::LLA, datum) = CovENUfromECEF(origin.lat, origin.lon)
CovENUfromECEF(origin::ECEF, datum) = ENUfromECEF(origin, datum) |> CovENUfromECEF
CovENUfromECEF(origin::UTMZ, datum) = ENUfromECEF(origin, datum) |> CovENUfromECEF
CovENUfromECEF(origin::UTM, zone::Integer, isnorth::Bool, datum) = ENUfromECEF(origin, zone, isnorth, datum) |> CovENUfromECEF

function CovENUfromECEF(lat::T, lon::T) where T <: Number
    sinλ, cosλ = sind(lat), cosd(lat)
    sinϕ, cosϕ = sind(lon), cosd(lon)
        
    R = [
               -sinλ          cosλ  0.0
        -cosλ * sinϕ -sinλ * sinϕ cosϕ
         cosλ * cosϕ  sinλ * cosϕ sinϕ
    ] |> SMatrix{3,3}

    return CovENUfromECEF(R, transpose(R))
end

struct CovECEFfromENU{T} <: CovMatrixRotation
    R::SMatrix{3,3,T}
    Rt::SMatrix{3,3,T}
end

CovECEFfromENU(lat::T, lon::T) where T <: Number = inv(CovENUfromECEF(lat, lon))

CovECEFfromENU(origin::ECEF{T}, lat::T, lon::T) where T <: Number = CovECEFfromENU(lat, lon)
CovECEFfromENU(origin::LLA, datum) = CovECEFfromENU(origin.lat, origin.lon)
CovECEFfromENU(origin::ECEF, datum) = ECEFfromENU(origin, datum) |> CovECEFfromENU
CovECEFfromENU(origin::UTMZ, datum) = ECEFfromENU(origin, datum) |> CovECEFfromENU
CovECEFfromENU(origin::UTM, zone::Integer, isnorth::Bool, datum) = ECEFfromENU(origin, zone, isnorth, datum) |> CovECEFfromENU

Base.inv(trans::CovECEFfromENU) = CovENUfromECEF(trans.Rt, trans.R)
Base.inv(trans::CovENUfromECEF) = CovECEFfromENU(trans.Rt, trans.R)

for (before, after) in [("ECEF", "ENU"), ("ENU", "ECEF")]
    CovM = Symbol("Cov" * after * "from" * before)
    arg = Symbol("Cov" * before)
    ret = Symbol("Cov" * after)

    @eval begin
        function (trans::$CovM)(x::$arg)
            Σx = covmat(x)
            Σ  = trans.R * Σx * trans.Rt

            return matcov($ret, Σ)
        end
    end
end