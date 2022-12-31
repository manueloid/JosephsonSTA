include("base.jl")
include("polynomials.jl")

using QuadGK
using SpecialPolynomials
##{{{ Fuctions to create the wavefunctions
"""
phase_integral(t::Float64, cp::ControlParameter) return the integral that goes up in the phase part of the wavefunction at a time t return the integral that goes up in the phase part of the wavefunction at a time t.
"""
function phase_integral(t::Float64, cp::ControlParameter)
    return quadgk(τ -> cp.ω0 / polynomial(τ, cp)^2, 0.0, t, rtol=1e-3, atol=1e-3)[1]
end

"""
hermite(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter) returns the hermite polynomial part of the wavefunction
"""
function hermite(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter)
    return SpecialPolynomials.basis(Hermite, n_index)(√(cp.ω0 / (2 * 2.0 / cp.NParticles)) * z / polynomial(t, cp))
end

"""
time_dep_Gaussian(cp::ControlParameter) return the time dependent function in the space dependent part of the wavefunction, the f in  exp(f(t) z²). 
"""
function time_dep_Gaussian(t::Float64, cp::ControlParameter)
    return im * (polynomial_1d(t, cp) / polynomial(t, cp) + im * cp.ω0 / polynomial(t, cp)^2) / (4.0 * 2.0 / cp.NParticles)
end

"""
normalization_factor(cp::ControlParameter) return the normalization factor for the wavefunction
"""
function normalization_factor(cp::ControlParameter)
    return (cp.ω0 / (2.0 * 2.0 / cp.NParticles * pi))^(0.25)
end

"""
create_wavefunction(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter) return the value of the mth state of the wavefunction for a time `t` and in position `z`
"""
function create_wavefunction(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter)
    return normalization_factor(cp) /
           √(2.0^n_index * factorial(n_index) * polynomial(t, cp)) *
           exp(-im * (n_index + 0.5) * phase_integral(t, cp)) *
           exp(time_dep_Gaussian(t, cp) * z^2) *
           hermite(z, t, n_index, cp)
end

"""
ground_state(z::Float64, t::Float64, cp::ControlParameter) return the value of the ground state at a time `t` and in position `z`
"""
function ground_state(z::Float64, t::Float64, cp::ControlParameter)
    return create_wavefunction(z, t, 0, cp)
end

"""
ground_state_1d(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter) return the first derivative of the ground state of the wavefunction for a time `t` and in position `z` 
"""
function ground_state_1d(z::Float64, t::Float64, cp::ControlParameter)
    return ground_state(z, t, cp) * 2.0 * z * time_dep_Gaussian(t, cp)
end

"""
ground_state_1d(z::Float64, t::Float64, n_index::Int64, cp::ControlParameter) return the first derivative of the ground state of the wavefunction for a time `t` and in position `z` 
"""
function ground_state_2d(z::Float64, t::Float64, cp::ControlParameter)
    return ground_state(z, t, cp) * (4z^2 * time_dep_Gaussian(t, cp)^2 + 2time_dep_Gaussian(t, cp))
end
##}}}
## testing
cparam = standard_cp()
using BenchmarkTools
@btime create_wavefunction(0.1, 0.2, 2, cparam)