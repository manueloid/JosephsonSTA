##{{{using  Polynomials utilities and control ω
"""
polynomial(t::Float64, cp::ControlParameter) return the value of the control polynomial at the time `t` given some parameters
"""
function polynomial(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
    return 6(γ - 1) * (t / cp.final_time)^5 -
           15(γ - 1) * (t / cp.final_time)^4 +
           10(γ - 1) * (t / cp.final_time)^3 +
           1
end

"""
polynomial_1d(t::Float64, cp::ControlParameter) return the value of the second derivative of the polynomial at the time `t` for a given ControlParameter.
"""
function polynomial_1d(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
    return 1 / cp.final_time * (
        30(γ - 1) * (t / cp.final_time)^4 -
        60(γ - 1) * (t / cp.final_time)^3 +
        30(γ - 1) * (t / cp.final_time)^2
    )
end


"""
polynomial_2d(t::Float64, cp::ControlParameter) return the value of the second derivative of the polynomial at the time `t` for a given ControlParameter.
"""
function polynomial_2d(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
    return 1 / cp.final_time^2 * (
        120(γ - 1) * (t / cp.final_time)^3 -
        180(γ - 1) * (t / cp.final_time)^2 +
        60(γ - 1) * (t / cp.final_time)
    )
end

"""
control_ω(t::Float64, cp::ControlParameter) return the control `ω` given the control polynomial 
"""
function control_ω(t::Float64, cp::ControlParameter)
    -polynomial_2d(t, cp) / polynomial(t, cp) + cp.ω0^2 / polynomial(t, cp)^4
end

"""
correction_poly(cp::ControlParameter, correction_vector::Array{Float64,1}) takes the correction values and returns a lagrange polynomial that will then be added to the initial control parameter 
"""
function correction_poly(t::Float64, cp::ControlParameter, correction_vector::Array{Float64,1})
    ys = [0.0; correction_vector; 0.0] # faster than vcat([]...)
    xs = range(0.0, cp.final_time, length=length(ys)) |> collect
    return Lagrange(xs, ys)(t)
end

##}}}
