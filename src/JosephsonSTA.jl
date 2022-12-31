module JosephsonSTA

using QuadGK
using SpecialPolynomials

##{{{ Including the basic functionalities 
export ControlParameter
export offset
export standard_cp
include("base.jl")
##}}}

##{{{ Including polynomials and wavefunction utilities
export control_ω
include("polynomials.jl")
##}}}

##{{{Including the functions to produce the STA wavefunctions
export create_wavefunction
export ground_state
export ground_state_1d
export ground_state_2d
include("wavefunctions.jl")
##}}}


end


