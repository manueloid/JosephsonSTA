"""
ControlParameter(ω0::Float64, ωf::Float64, final_time::Float64, NParticles::Int64) store all the relevant quantities to the simulation
"""
struct ControlParameter
    ω0::Float64
    ωf::Float64
    final_time::Float64
    NParticles::Int64
end

"""
offset(cp::ControlParameter, offset::Float64) given an initial `ControlParameter` variable, shift the initial and final ω by a factor `offset`
"""
function offset(cp::ControlParameter, offset::Float64)
    offset_cp = ControlParameter(cp.ω0 + offset, cp.ωf + offset, cp.final_time, cp.NParticles)
    return offset_cp
end

"""
standard_cp() utility to return a standard control paramter to play around
"""
function standard_cp()
    return ControlParameter(2.0, 14.2, 0.3, 50)
end
