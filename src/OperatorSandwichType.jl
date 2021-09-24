struct OperatorSandwich{T1,T2}
    Op::OperatorShift
    ψ::T1
    ψL::T2
    δT::Function
    ∂ωT::Function
    weight::Function
end

function OperatorSandwich(ψL,Op::OperatorShift,weight,ψ; params=())
    δTs = Op.δT
    ∂ωTs = Op.∂ωT

    ψ1 = ∂r(ψ)
    ψ2 = ∂θ(ψ)
    ψ11 = ∂r(ψ1);
    ψ12 = ∂r(ψ2);
    ψ22 = ∂θ(ψ2);

    FF1 = let r₊ = ψ.R.r₊, r₋ = ψ.R.r₋, a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, δTs=δTs, ψ11 = ψ11,ψ12 = ψ12,ψ22 = ψ22,ψ1 = ψ1,ψ2 = ψ2
    function F1(r,z)
        ψL(r,z)*δTs(r,z,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,params...)*weight(r,z)
    end
    end

    FF2 = let r₊ = ψ.R.r₊, r₋ = ψ.R.r₋, a = ψ.a, m = ψ.m ,ω = ψ.ω ,s = ψ.s, ψ=ψ, ψL=ψL, δTs=δTs, ψ11 = ψ11,ψ12 = ψ12,ψ22 = ψ22,ψ1 = ψ1,ψ2 = ψ2
    function F2(r,z)
        ψL(r,z)*∂ωTs(r,z,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,params...)*weight(r,z)
    end
    end
    OperatorSandwich(Op,ψ,ψL,FF1,FF2,weight)
end
