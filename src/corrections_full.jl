##{{{ Calculation of the Kns
"""
gradient_element(xs::Vector{Float64}, ys::Vector{Float64}, index::Int64, t::Float64) return the nth Lagrange polynomials given the arrays of the points that need to be interpolated
"""
function gradient_element(xs::Vector{Float64}, ys::Vector{Float64}, index::Int64, t::Float64)
    ys[index] = 1.0
    value = Lagrange(xs, ys)(t)
    ys[index] = 0.0
    return value
end
"""
gradient_λ(number_points::Int64, final_time::Float64, t::Float64) return the gradient with respect of lambda at time `t`
"""
function gradient_λ(number_points::Int64, final_time::Float64, t::Float64)
    xs = range(0.0, final_time, number_points + 2) |> collect
    ys = zeros(Float64, number_points + 2)
    return [gradient_element(xs, ys, index, t) for index in 2:number_points+1]
end

"""
This function returns the Kns given the fuction ⟨2|z²|0⟩:= z_square_int
"""
function Kns(cp::ControlParameter, nlambda::Int64=5)
    return 0.25 * quadgk(t -> gradient_λ(nlambda, cp.final_time, t) * z_square_analytic(t, cp), 0.0, cp.final_time)[1]
end
##}}}
##{{{ Common functions for both types of Hamiltonian
bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0 # One-line function that returns the piecewise function bₕ(z)
"""
This function returns the sum of the two following integrals ⟨2|z²|0⟩ and ⟨2|∂²|0⟩ with respect of time
"""
function Gn_analytic(cp::ControlParameter)
    return quadgk(
        t ->
            (2.0 / cp.NParticles)^2 * second_deriv_analytic(t, cp) -
            z_square_analytic(t, cp),
        0.0,
        cp.final_time,
        rtol=1e-5,
    )[1]
end
precompile(Gn_analytic, (ControlParameter,))
##}}}
##{{{Calculation of the Gns for the intermediate Hamiltonian

dbh(z::Float64, h::Float64) = abs(z) <= 1 ? -z * √((1 + z + h) * (1 - z)) : 0
function create_wavefunction_deriv(z::Float64, t::Float64, cp::ControlParameter)
    return create_wavefunction(z, t, 0, cp) * 2.0 * z * time_dep_Gaussian(t, cp)
end
function create_wavefunction_deriv_second(z::Float64, t::Float64, cp::ControlParameter)
    return create_wavefunction(z, t, 0, cp) * (4z^2 * time_dep_Gaussian(t, cp)^2 + 2time_dep_Gaussian(t, cp))
end
function first_element(mth_bra::Int64, cp::ControlParameter)
    return hcubature(
        x -> -(2.0 / cp.NParticles)^2 * dbh(x[1], 0.0) *
             conj(create_wavefunction(x[1], x[2], mth_bra, cp)) *
             create_wavefunction_deriv(x[1], x[2], cp),
        [-1.0, 0.0], [1.0, cp.final_time],
        rtol=1e-5
    )[1]
end
precompile(first_element, (Int64, ControlParameter))
function second_element(mth_bra::Int64, cp::ControlParameter)
    return hcubature(
        x -> -(2.0 / cp.NParticles)^2 * bh(x[1], 0.0) *
             conj(create_wavefunction(x[1], x[2], mth_bra, cp)) *
             create_wavefunction_deriv_second(x[1], x[2], cp),
        [-1.0, 0.0], [1.0, cp.final_time],
        rtol=1e-5
    )[1]
end
precompile(second_element, (Int64, ControlParameter))
function third_element(mth_bra::Int64, cp::ControlParameter)
    return hcubature(
        x -> -2bh(x[1], 0.0) * conj(create_wavefunction(x[1], x[2], mth_bra, cp)) *
             create_wavefunction(x[1], x[2], 0, cp),
        [-1.0, 0.0], [1.0, cp.final_time],
        rtol=1e-5
    )[1]
end
precompile(third_element, (Int64, ControlParameter))

"""
This function returns the integral with respect of time of the terms connected to the bₕ function. Check the notes for more information
"""
function Gn_bh(mth_bra::Int64, cp::ControlParameterInt)
    return first_element(mth_bra, cp) + second_element(mth_bra, cp) + third_element(mth_bra, cp)
end
precompile(Gn_bh, (Int64, ControlParameterInt,))

"""
Gn(m::Int64, cp::ControlParameter) return the value of the singular Gn for a given m.
For m ≢ 2, we know that the analytic part of the integral is 0 and thus we are only interested in Gn_bh
"""
function Gn(m::Int64, cp::ControlParameterInt)
    if m == 2
        result = Gn_bh(m, cp) + Gn_analytic(cp)
    else
        result = Gn_bh(m, cp)
    end
end

precompile(Gn, (Int64, ControlParameterInt))
# finally, as we can check from the notes, we are only interested in the K₂ term as the other ones are all 0. Let us focus only on m = 2 
##}}} ##src
##{{{Calculation of the Gns for the full Hamiltonian ##src
"""
Gn_bh(mth_bra::Int64, cp::ControlParameterFull) integration over the two variables of the term depending on bh(z)
"""
function Gn_bh(mth_bra::Int64, cp::ControlParameterFull)
    h = 2.0 / cp.NParticles
    return -2.0 * hcubature(
        x -> conj(create_wavefunction(x[1], x[2], mth_bra, cp)) *
             bh(x[1], h) *
             ground_state(x[1] + h, x[2], cp),
        [-20.0, 0.0],
        [20.0, cp.final_time],
        rtol=1e-5,
    )[1]
end
precompile(Gn_bh, (Int64, ControlParameterFull))

"""
Gn(m::Int64, cp::ControlParameter) return the value of the singular Gn for a given m.
For m ≢ 2, we know that the analytic part of the integral is 0 and thus we are only interested in Gn_bh
"""
function Gn(m::Int64, cp::ControlParameterFull)
    if m == 2
        result = Gn_bh(m, cp) + Gn_analytic(cp)
    else
        result = Gn_bh(m, cp)
    end
end
precompile(Gn, (Int64, ControlParameter))
##}}} ##src
##{{{Calculation of corrections for the original eSTA##src
"""
This function will evaluate the corrections given the maximum numbers of STA wavefunctions we want to use and the number of points we want to interpolate
"""
function corrections(cp::ControlParameter, nlambda::Int64=5, max_bra::Int64=4)
    gns = 0.0
    g2 = Gn(2, cp)
    for mbra = 4:2:max_bra
        gns += Gn(mbra, cp) |> abs2
    end
    k2 = Kns(cp, nlambda)
    return (gns + abs2(g2)) * real(conj(g2) * k2) /
           ((real(conj(g2) * k2))' * (real(conj(g2) * k2)))
end
precompile(corrections, (ControlParameter, Int64, Int64,))
##}}} ##src
