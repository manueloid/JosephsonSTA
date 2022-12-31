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
standard_cp() utility to return a standard control paramter to play around
"""
function standard_cp()
    return ControlParameter(2.0, 14.2, 0.3, 50)
end

# One liners to create new ControlParameter variables from existing ones
cp_time(cp::ControlParameter, time::Float64) = ControlParameter(cp.ω0, cp.ωf, time, cp.NParticles)
cp_particles(cp::ControlParameter, particles::Int64) = ControlParameter(cp.ω0, cp.ωf, cp.final_time, particles)
offset(cp::ControlParameter, offset::Float64) = ControlParameter(cp.ω0 + offset, cp.ωf + offset, cp.final_time, cp.NParticles)
