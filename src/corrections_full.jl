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
function Kns(nlambda::Int64, cp::ControlParameter)
    return 0.25 * quadgk(t -> gradient_λ(nlambda, cp.final_time, t) * z_square_analytic(t, cp), 0.0, cp.final_time)[1]
end
##}}}
##{{{Calculation of the Gns ##src

bh(z::Float64, h::Float64) = abs(z) <= 1 ? √((1 + z + h) * (1 - z)) : 0 # One-line function that returns the piecewise function bₕ(z)

"""
Gn_bh(mth_bra::Int64, cp::ControlParameter) integration over the two variables of the term depending on bh(z)
"""
function Gn_bh(mth_bra::Int64, cp::ControlParameter)
    h = 2.0 / cp.NParticles
    return -2.0 * hcubature(
        x -> conj(create_wavefunction(x[1], x[2], mth_bra, cp)) *
             bh(x[1], 0.0) *
             ground_state(x[1] + h, x[2], cp),
        [-20.0, 0.0],
        [20.0, cp.final_time],
        rtol=1e-5,
    )[1]
end

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

Gn(m::Int64, cparam::ControlParameter) = m == 2 ? Gn_bh(m, cparam) + Gn_analytic(cparam) : Gn_bh(m, cparam)
# finally, as we can check from the notes, we are only interested in the K₂ term as the other ones are all 0. Let us focus only on m = 2
##}}} ##src
##{{{Calculation of corrections for the original eSTA##src
"""
This function will evaluate the corrections given the maximum numbers of STA wavefunctions we want to use and the number of points we want to interpolate
"""
function corrections(nlambda::Int64, max_bra::Int64, cparam::ControlParameter)
    gns = 0.0
    g2 = Gn(2, cparam)
    for mbra = 4:2:max_bra
        gns += Gn(mbra, cparam) |> abs2
    end
    k2 = Kns(nlambda, cparam)
    return (gns + abs2(g2)) * real(conj(g2) * k2) /
           ((real(conj(g2) * k2))' * (real(conj(g2) * k2)))
end
##}}} ##src
