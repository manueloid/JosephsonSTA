module JosephsonSTA

export ControlParameter
# Write your package code here.
# Fundamental structure of the whole code is the ControlParameter one that contains all the values I need for the calculation.
struct ControlParameter
    ω0::Float64
    ωf::Float64
    final_time::Float64
    NParticles::Int64
end 

end
using DataFrames
