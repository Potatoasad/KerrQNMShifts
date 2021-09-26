using ContourIntegrals

### Use expressions for the operators, noted in their respective csvs
file1 = "OperatorShifts/KerrNewman/TuekolskyShifts.csv"
file2 = "OperatorShifts/KerrNewman/TuekolskyFrequencyDerivative.csv"
T = OperatorShift(file1,file2)

### Get the integrand in the contour integral for the inner product computation
ψ = qnmfunctionnew(2,2,2,0,0.5)
weight = let r₊ = ψ.R.r₊ , r₋ = ψ.R.r₋ , s = ψ.s
    (r,z) -> ((r-r₊)*(r-r₋))^s
end
Ts = OperatorSandwich(ψ,T,weight,ψ)
δT = Ts.δT;
∂ωT = Ts.∂ωT;

### Compute the expressions so that they compile the first time
δT(2.2,0.1) |> println
∂ωT(2.2,0.1) |> println

### Define the useful contours
r₊ = ψ.R.r₊ ; r₋ = ψ.R.r₋ ; s = ψ.s ; Δr = 0.1*(r₊-r₋); ϵ = eps(0.1);

point1 = r₊ + Δr - Δr*im
point2 = r₊ - Δr - Δr*im

radial1 = SemiInfiniteLine(point1 , point1 + Δr*im , false)
angular = LineSegment(-1.0+100*ϵ , 1.0-100*ϵ , true) #to avoid the NaNs at the edges
C1 = radial1 ⊗ angular

radial2 = LineSegment(point1,point2,true)
C2 = radial2 ⊗ angular

radial3 = SemiInfiniteLine(point2 , point2 + Δr*im , true)
C3 = radial3 ⊗ angular

TheContour = C1⊕C2⊕C3

### Define Transformed Functions (that live on the compactified domain [0,1]⊗[0,1])
δTₒ = TransformIntegrand(δT , thedomain)
∂ωTₒ = TransformIntegrand(∂ωT , thedomain)

### Compute the expressions so that they compile the first time
δTₒ(0.1,0.1) |> println
∂ωTₒ(0.1,0.1) |> println


#=
### Use expressions for the operators, noted in their respective csvs
file1b = "OperatorShifts/KerrNewman0/TuekolskyShifts.csv"
file2b = "OperatorShifts/KerrNewman0/TuekolskyFrequencyDerivative.csv"
Tb = OperatorShift(file1b,file2b)

### Get the integrand in the contour integral for the inner product computation
Tsb = OperatorSandwich(ψ,Tb,weight,ψ)
δTb = Tsb.δT;
∂ωTb = Tsb.∂ωT;

### Compute the expressions so that they compile the first time
δTb(2.2,0.1) |> println
∂ωTb(2.2,0.1) |> println

### Define Transformed Functions (that live on the compactified domain [0,1]⊗[0,1])
δTbₒ = TransformIntegrand(δTb , thedomain)
∂ωTbₒ = TransformIntegrand(∂ωTb , thedomain)

### Compute the expressions so that they compile the first time
δTbₒ(0.1,0.1) |> println
∂ωTbₒ(0.1,0.1) |> println
=#
