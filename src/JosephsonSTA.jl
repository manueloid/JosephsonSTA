module JosephsonSTA

using QuadGK
using SpecialPolynomials
using QuadGK
using HCubature
using QuantumOptics
using DataFrames
using CSV

##{{{ Including the basic functionalities 
export ControlParameter
export ControlParameterFull
export ControlParameterInt
export offset
export cp_time
export cp_nparticles
include("base.jl")
##}}}

##{{{ Including polynomials and wavefunction utilities
export control_ω
export correction_poly
include("polynomials.jl")
##}}}

##{{{Including the functions to produce the STA wavefunctions
export create_wavefunction
export ground_state
export ground_state_1d
export ground_state_2d
export z_square_analytic
export second_deriv_analytic
include("wavefunctions.jl")
##}}}

export corrections
include("corrections_full.jl")

##{{{ Including utilities to calculate the fidelity of the protocol
export fidelity_single
export fidelity_multiple
include("simulation.jl")
##}}}

end


