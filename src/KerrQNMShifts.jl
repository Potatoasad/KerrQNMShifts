module KerrQNMShifts

using KerrQuasinormalModes
using PyCall
using CSV
using DataFrames
using Plots
using Cubature
# Write your package code here.

function GetExprMetric(file)
    df = DataFrame(CSV.File(file,types=Dict(1=>String)))
    exprs = []
    for (d,f) ∈ eachrow(df[df.Coefficient .!= "0",[:Index,:Coefficient]])
        #println((d,f))
        d = d == "-" ? "" : d
        d = replace(d,"th"=>"θ")
        ps = "ψ"*d
        push!(exprs,"($f)*$(ps)(r,x)*ψL(r,x)")
    end
    func_def = "$(join(exprs," + "))"
    println(func_def)
    Meta.parse(func_def)
end

function GetExprDerivative(file)
    df = DataFrame(CSV.File(file,types=Dict(1=>String)))
    exprs = []
    for (d,f) ∈ eachrow(df[df.Coefficient .!= "0",:])
        #println((d,f))
        d = d == "-" ? "" : d
        d = replace(d,"th"=>"θ")
        ps = "ψ"*d
        push!(exprs,"($f)*$(ps)(r,x)*ψL(r,x)")
    end
    func_def = "$(join(exprs," + "))"
    println(func_def)
    Meta.parse(func_def)
end


function MakeδT(file)
    thisexpr = GetExprMetric("TuekolskyShifts.csv")
    @eval δT(r,x,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,ψL; M=1) = $thisexpr
    δT
end

#δT = MakeδT("TuekolskyShifts.csv")

function Make∂ωT(file)
    thisexpr2 = GetExprDerivative("TuekolskyFrequencyDerivative.csv")
    @eval ∂ωT(r,x,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,ψL; M=1) = $thisexpr2
    ∂ωT
end

#∂ωT = Make∂ωT("TuekolskyFrequencyDerivative.csv")


#NEW SPEEDUP TRY
function DoubleIntegrand1(ψL,δTs::Function,ψ; β=0.5)
    ψ1 = ∂r(ψ)
    ψ2 = ∂θ(ψ)
    ψ11 = ∂r(ψ1);
    ψ12 = ∂r(ψ2);
    ψ22 = ∂θ(ψ2);
    FF = let r₊ = ψ.R.r₊, r₋ = ψ.R.r₋, a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, δTs=δTs, ψ11 = ψ11,ψ12 = ψ12,ψ22 = ψ22,ψ1 = ψ1,ψ2 = ψ2
    function F(r,z)
        δTs(r,z,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,ψL)*((r-r₊)*(r-r₋))^s
    end
    end
    Δr = (ψ.R.r₊ - ψ.R.r₋)*0.4
    let r₊ = ψ.R.r₊, Δr = Δr
        (t,z) -> β*im*((FF(r₊ - Δr*(im+1.0) + im*β*t/(1-t),z))/((1-t)^2))
    end
end

function DoubleIntegrand2(ψL,δTs::Function,ψ; β=0.5)
    ψ1 = ∂r(ψ)
    ψ2 = ∂θ(ψ)
    ψ11 = ∂r(ψ1);
    ψ12 = ∂r(ψ2);
    ψ22 = ∂θ(ψ2);
    FF = let r₊ = ψ.R.r₊, r₋ = ψ.R.r₋, a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, δTs=δTs, ψ11 = ψ11,ψ12 = ψ12,ψ22 = ψ22,ψ1 = ψ1,ψ2 = ψ2
    function F(r,z)
        δTs(r,z,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,ψL)*((r-r₊)*(r-r₋))^s
    end
    end
    Δr = (ψ.R.r₊ - ψ.R.r₋)*0.4
    let r₊ = ψ.R.r₊, Δr = Δr
        (t,z) -> β*im*((-FF(r₊ - Δr*(im-1.0) + im*β*t/(1-t),z))/((1-t)^2))
    end
end


function CircleIntegrand(ψL,δTs::Function,ψ)
    ψ1 = ∂r(ψ)
    ψ2 = ∂θ(ψ)
    ψ11 = ∂r(ψ1);
    ψ12 = ∂r(ψ2);
    ψ22 = ∂θ(ψ2);
    FF = let r₊ = ψ.R.r₊, r₋ = ψ.R.r₋, a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, δTs=δTs, ψ11 = ψ11,ψ12 = ψ12,ψ22 = ψ22,ψ1 = ψ1,ψ2 = ψ2
    function F(r,z)
        δTs(r,z,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,ψL)*((r-r₊)*(r-r₋))^s
    end
    end
    #(t,z) -> (-2*Δr*FF(r₊ + Δr*(2*t-(im+1.0)),z))
    Δr = (ψ.R.r₊ - ψ.R.r₋)*0.4
    let r₊ = ψ.R.r₊, Δr = Δr
        (t,z) -> (-2*Δr*FF(r₊ + Δr*(1.0-2.0*t-im),z))
    end
end

function importqnm()
    global qnm = pyimport("qnm")
    qnm.download_data()
end


qnm = pyimport("qnm")
qnm.download_data()
function qnmfunctionnew(s,l,m,n,a; qnm=qnm)
    grav_freq = qnm.modes_cache(s=s,l=l,m=m,n=n)
    ω, Alm, Cllʼ = grav_freq(a=a)
    qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω=ω,Alm=Alm,Cllʼ=Cllʼ)
end

function ComputeShiftold(s,l,m,n,a; ϵ = 10^(-10))
    #println("computing qnm")
    ψ = qnmfunctionnew(s,l,m,n,a)
    #println("Onto Defining Functions")
    F₊₁ = DoubleIntegrand1(ψ,δT,ψ)
    F₊₂ = DoubleIntegrand2(ψ,δT,ψ)
    Fₒ = CircleIntegrand(ψ,δT,ψ)
    G₊₁ = DoubleIntegrand1(ψ,∂ωT,ψ)
    G₊₂ = DoubleIntegrand2(ψ,∂ωT,ψ)
    Gₒ = CircleIntegrand(ψ,∂ωT,ψ)
    #println("THe value of F₊ is")
    #println(F₊(2.1,0.2))
    numerator1 = hcubature(t-> F₊₁(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    numerator2 = hcubature(t-> F₊₂(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    numerator3 = hcubature(t-> Fₒ(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    denomenator1 = hcubature(t-> G₊₁(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    denomenator2 = hcubature(t-> G₊₂(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    denomenator3 = hcubature(t-> Gₒ(t[1],t[2]), [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)])
    δω =  -(numerator1[1]+numerator2[1]+numerator3[1])/(denomenator1[1]+denomenator2[1]+denomenator3[1])
    ψ.ω,δω
end

function Converter(f1s::Function)
function f2(x,v)
    kk = size(x)
    thesize = length(kk) == 1 ? 1 : kk[2]
    for i ∈ 1:(thesize)
        a = f1s(x[1,i],x[2,i])
        v[1,i] = real(a)
        v[2,i] = imag(a)
    end
end
f2
end

function ComputeShift(s,l,m,n,a,δT::Function,∂ωT::Function; ϵ = 10^(-10), abstol = 10^(-10))
    #println("computing qnm")
    println("newfunc")
    ψ = qnmfunctionnew(s,l,m,n,a)
    #println("Onto Defining Functions")
    #IntegralCHanger
    Δr = (ψ.R.r₊-ψ.R.r₋)
    wheretosendit = 0.8
    β = Δr*(1/wheretosendit - 1)
    F₊₁ = DoubleIntegrand1(ψ,δT,ψ; β=β)
    F₊₂ = DoubleIntegrand2(ψ,δT,ψ; β=β)
    Fₒ = CircleIntegrand(ψ,δT,ψ)
    G₊₁ = DoubleIntegrand1(ψ,∂ωT,ψ; β=β)
    G₊₂ = DoubleIntegrand2(ψ,∂ωT,ψ; β=β)
    Gₒ = CircleIntegrand(ψ,∂ωT,ψ)
    F₊₁(0.4,0.3) |> println
    F₊₂(0.4,0.3) |> println
    Fₒ(0.4,0.3) |> println
    G₊₁(0.4,0.3) |> println
    G₊₂(0.4,0.3) |> println
    Gₒ(0.4,0.3) |> println
    FF₊₁ = F₊₁ |> Converter;
    FF₊₂ = F₊₂ |> Converter;
    FFₒ = Fₒ  |> Converter;
    GG₊₁ = G₊₁ |> Converter;
    GG₊₂ = G₊₂ |> Converter;
    GGₒ = Gₒ |> Converter;
    #println("THe value of F₊ is
    #println(F₊(2.1,0.2))
    println("starting")
    numerator1 = pcubature_v(2, FF₊₁, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("First Done")
    numerator2 = pcubature_v(2, FF₊₂, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("Second Done")
    tovec = (x -> ([real(x[1]),imag(x[1])],x[2]))
    numerator3 = pcubature_v(2, FFₒ, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("3 Done")
    denomenator1 = pcubature_v(2, GG₊₁, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("4 Done")
    denomenator2 = pcubature_v(2, GG₊₂, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("5 Done")
    denomenator3 = pcubature_v(2, GGₒ, [(0.0+ϵ),-(1.0-ϵ)], [(1.0-ϵ),(1.0-ϵ)]; error_norm = Cubature.PAIRED, abstol=abstol)
    print("6 Done")
    tonum = (x -> x[1]+im*x[2])
    numer = tonum(numerator1[1]+numerator2[1]+numerator3[1]);
    δnumer = sum((numerator1[2]+numerator2[2]+numerator3[2]).^2) |> sqrt;
    denom = tonum(denomenator1[1]+denomenator2[1]+denomenator3[1]);
    δdenom = sum((denomenator1[2]+denomenator2[2]+denomenator3[2]).^2) |> sqrt;
    δω =  -tonum(numerator1[1]+numerator2[1]+numerator3[1])/tonum(denomenator1[1]+denomenator2[1]+denomenator3[1])
    δω_error = sqrt((abs(1/denom)*δnumer)^2 + (abs(numer/(denom^2))*δdenom)^2)
    ψ.ω,δω,δω_error
end


export qnm, MakeδT, Make∂ωT, ComputeShift, importqnm

end
