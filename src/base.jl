abstract type ControlParameter end
"""
ControlParameterFull(ω0::Float64, ωf::Float64, final_time::Float64, NParticles::Int64) store all the relevant quantities to the simulation
"""
mutable struct ControlParameterFull <: ControlParameter
    ω0::Float64
    ωf::Float64
    final_time::Float64
    NParticles::Int64
end

"""
ControlParameterInt(ω0::Float64, ωf::Float64, final_time::Float64, NParticles::Int64) store all the relevant quantities to the simulation
"""
mutable struct ControlParameterInt <: ControlParameter
    ω0::Float64
    ωf::Float64
    final_time::Float64
    NParticles::Int64
end

# One liners to create new ControlParameter variables from existing ones
ControlParameterFull() = ControlParameterFull(2.0, 2.0√(1.0 + 50.0), 0.3, 50)
ControlParameterInt() = ControlParameterInt(2.0, 2.0√(1.0 + 50.0), 0.3, 50)
ControlParameterFull(tf::Float64, np::Int64) = ControlParameterFull(2.0, 2.0√(1.0 + 50.0),tf, np)
ControlParameterInt(tf::Float64, np::Int64) = ControlParameterInt(2.0, 2.0√(1.0 + 50.0),tf, np)
ControlParameterFull( np::Int64, tf::Float64) = ControlParameterFull(2.0, 2.0√(1.0 + 50.0),tf, np)
ControlParameterInt( np::Int64, tf::Float64) = ControlParameterInt(2.0, 2.0√(1.0 + 50.0),tf, np)

# Updating time and doing some offset
cp_time(cp::ControlParameter, t::Float64) = typeof(cp)(cp.ω0, cp.ωf, t, cp.NParticles)
cp_nparticles(cp::ControlParameter, np::Int64) = typeof(cp)(cp.ω0, cp.ωf, cp.final_time, np)
offset(cp::ControlParameter, offset::Float64) = typeof(cp)(cp.ω0 + offset, cp.ωf + offset, cp.final_time, cp.NParticles)
