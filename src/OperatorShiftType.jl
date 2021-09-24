struct OperatorShift
    filenameδT::AbstractString
    filename∂ωT::AbstractString
    δT::Function
    ∂ωT::Function
end

function OperatorShift(file1::AbstractString,file2::AbstractString; additional_params = ())
    δT = MakeδT(file1; additional_params = additional_params)
    ∂ωT = Make∂ωT(file2; additional_params = additional_params)
    OperatorShift(file1,file2,δT,∂ωT)
end
