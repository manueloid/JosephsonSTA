"""
piecewise(t::Float64, cp::ControlParameter, f::Function) take a function `f` and generates a piecewise function of the form 
`f(t) = f(0) if t < 0
		f(t) if 0 <= t <= cp.final_time
		f(cp.final_time) if t > cp.final_time`
"""
function piecewise(t::Float64, cp::ControlParameter, f::Function)
    if t < 0
        return f(0.0)
    elseif t > cp.final_time
        return f(cp.final_time)
    else
        return f(t)
    end
end
"""
polynomial(t::Float64, cp::ControlParameter) return the value of the control polynomial at the time `t` given some parameters as a piecewise function in the whole real line.
"""
function polynomial(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
	f(t) = 6(γ - 1) * (t / cp.final_time)^5 -
           15(γ - 1) * (t / cp.final_time)^4 +
           10(γ - 1) * (t / cp.final_time)^3 +
           1
	return piecewise(t, cp, t -> f(t))
end
"""
polynomial_1d(t::Float64, cp::ControlParameter) return the value of the first derivative of the polynomial at the time `t` for a given ControlParameter as a piecewise function on the whole line.
"""
function polynomial_1d(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
	f(t)= 1 / cp.final_time * (
        30(γ - 1) * (t / cp.final_time)^4 -
        60(γ - 1) * (t / cp.final_time)^3 +
        30(γ - 1) * (t / cp.final_time)^2
    )
	return piecewise(t, cp, t -> f(t))
end
"""
polynomial_2d(t::Float64, cp::ControlParameter) return the value of the second derivative of the polynomial at the time `t` for a given ControlParameter on the whole line
"""
function polynomial_2d(t::Float64, cp::ControlParameter)
    γ = √(cp.ω0 / cp.ωf)
	f(t)=1 / cp.final_time^2 * (
        120(γ - 1) * (t / cp.final_time)^3 -
        180(γ - 1) * (t / cp.final_time)^2 +
        60(γ - 1) * (t / cp.final_time)
    )
	return piecewise(t, cp, t -> f(t))
end
"""
control_ω(t::Float64, cp::ControlParameter) return the control `ω` given the control polynomial 
"""
function control_ω(t::Float64, cp::ControlParameter)
    -polynomial_2d(t, cp) / polynomial(t, cp) + cp.ω0^2 / polynomial(t, cp)^4
end
"""
correction_polyin(cp::ControlParameter, correction_vector::Array{Float64,1}) takes the correction values and returns a lagrange polynomial that will then be added to the initial control parameter.
The function is defined in the interval [0, cp.final_time]
"""
function correction_polyin(t::Float64, cp::ControlParameter, correction_vector::Array{Float64,1})
    ys = [0.0; correction_vector; 0.0] # faster than vcat([]...)
    xs = range(0.0, cp.final_time, length=length(ys)) |> collect
    return Lagrange(xs, ys)(t)
end
"""
correction_poly(cp::ControlParameter, correction_vector::Array{Float64,1}) takes the correction values and returns a lagrange polynomial that will then be added to the initial control parameter.
The function is defined on the whole real line.
"""
correction_poly(t::Float64, cp::ControlParameter, correction_vector::Array{Float64,1}) = piecewise(t, cp, t -> correction_polyin(t, cp, correction_vector))
"""
control_ad(t::Float64, cp::ControlParameter) return the control `ω` for the adiabatic protocol for the given the control parameter. This is just a polynomial that interpolates the the initial and final value with the extra condition that the derivatives needs to be zero at initial and final time.
"""
function control_ad(t::Float64, cp::ControlParameter)
	ω0, ωf, tf = cp.ω0, cp.ωf, cp.final_time
	f(t) = ω0^2 - 3(ω0^2 - ωf^2) * (t / tf)^2 + 2(ω0^2 - ωf^2) * (t / tf)^3
	return piecewise(t, cp, t -> f(t))
end
