function GetExprMetric(file)
    df = DataFrame(CSV.File(file,types=Dict(1=>String)))
    exprs = []
    for (d,f) ∈ eachrow(df[df.Coefficient .!= "0",[:Index,:Coefficient]])
        #println((d,f))
        d = d == "-" ? "" : d
        d = replace(d,"th"=>"θ")
        ps = "ψ"*d
        push!(exprs,"($f)*$(ps)(r,x)")
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
        push!(exprs,"($f)*$(ps)(r,x)")
    end
    func_def = "$(join(exprs," + "))"
    println(func_def)
    Meta.parse(func_def)
end

function MakeδT(file; additional_params = ())
    thisexpr = GetExprMetric(file)
    @eval δT(r,x,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,$(additional_params...); M=1) = $thisexpr
    δT
end

#δT = MakeδT("TuekolskyShifts.csv")

function Make∂ωT(file; additional_params = ())
    thisexpr2 = GetExprDerivative(file)
    @eval ∂ωT(r,x,a,m,ω,s,ψ11,ψ12,ψ22,ψ1,ψ2,ψ,$(additional_params...); M=1) = $thisexpr2
    ∂ωT
end

#∂ωT = Make∂ωT("TuekolskyFrequencyDerivative.csv")
