struct OperatorShift
    filename::AbstractString
    δT::Function
    ∂ωT::Function
end

function OperatorShift(file::AbstractString; additional_params = ())
    δT = MakeδT(file; additional_params = additional_params)
    ∂ωT = Make∂ωT(file; additional_params = additional_params)
    OperatorShift(file,δT,∂ωT)
end
