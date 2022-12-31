
"""
constant_quantities(cp::ControlParameter) define the operators and the initial and final ground states of the protocol.
The result is a Dict.
The idea is to evaluate the fidelity for different final times, so these values will remain the same and there is no need to define them every time.
"""
function constant_quantities(cp::ControlParameter)
    Jz = sigmaz(SpinBasis(cp.NParticles / 2)) / 2 |> dense
    Jx = sigmax(SpinBasis(cp.NParticles / 2)) / 2
    ψ0 = eigenstates(2.0 / cp.NParticles * (cp.ω0^2 / 4.0 - 1.0) * Jz^2 - 2.0 * Jx, 1)[2][1]
    ψf = eigenstates(2.0 / cp.NParticles * (cp.ωf^2 / 4.0 - 1.0) * Jz^2 - 2.0 * Jx, 1)[2][1]
    return (Jz=Jz, Jx=Jx, ψ0=ψ0, ψf=ψf)
end

"""
fidelity_single(nlambda::Int64, maxbra::Int64, cp::ControlParameter) utility to calculate the fidelity on the go, with just the desired ControlParameter variable.
Not suitable to use to calculate the fidelities for different final times as the there is no need to define every time the operators and ground states as they do not depend from final time.
"""
function fidelity_single(cp::ControlParameter, nlambda::Int64=5, maxbra::Int64=4)
    qt = constant_quantities(cp)
    corr = corrections(cp, nlambda, maxbra)
    ω_esta(t) = control_ω(t, cp) - correction_poly(t, cp, corr)
    fidelity(t, psi) = abs2.(dagger(qt.ψf) * psi)
    function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
        return (2.0 / cp.NParticles * (ω_esta(t) / 4.0 - 1.0) * qt.Jz^2 - 2.0 * qt.Jx)
    end
    return timeevolution.schroedinger_dynamic([0.0, final_time], qt.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
end

"""
fidelity_multiple(cp::ControlParameter, final_times::Vector{Float64}, nlambda::Int64 = 5, maxbra::Int64 = 4) loop through a given array of final times to return the corresponding fidelities as an array.
"""
function fidelity_multiple(cp::ControlParameter, final_times::Vector{Float64}, nlambda::Int64=5, maxbra::Int64=4)
    qt = constant_quantities(cp)
    fidelities = ones(length(final_times))
    Threads.@threads for (index, tf) in enumerate(final_times) |> collect
        local cparam = cp_time(cp, tf)
        corr = corrections(cparam, nlambda, maxbra)
        ω_esta(t) = control_ω(t, cp) - correction_poly(t, cp, corr)
        fidelity(t, psi) = abs2.(dagger(qt.ψf) * psi)
        # Here I need to define the control parameter
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
            return (2.0 / cp.NParticles * (ω_esta(t) / 4.0 - 1.0) * qt.Jz^2 - 2.0 * qt.Jx)
        end
        println("Evaluating the fidelity for final time $tf")
        fidelities[index] = timeevolution.schroedinger_dynamic([0.0, tf], qt.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
        println("Done")
    end
    return fidelities
end

"""
fidelity_multiple(data::DataFrame) multiple dispatch of the fidelity_multiple function to act with DataFrame object in case the corrections are stored in some other file.
It works in the same way as the other methods.
"""
function fidelity_multiple(data::DataFrame)
    seed = corrections_control_parameter(data, data.tf[1])[2]
    qt = constant_quantities(seed)
    fidelities = ones(length(data.tf))
    Threads.@threads for (index, final_time) in enumerate(data.tf) |> collect
        local corrections, cparam = corrections_control_parameter(data, final_time)#Here I will define the cparam via some read functions
        ω_esta(t) = control_ω(t, cparam) - correction_poly(t, cparam, corrections)
        fidelity(t, psi) = abs2.(dagger(qt.ψf) * psi)
        # Here I need to define the control parameter
        function H_eSTA(t, psi) # Function to return the time dependent Hamiltonian
            return (2.0 / cparam.NParticles * (ω_esta(t) / 4.0 - 1.0) * qt.Jz^2 - 2.0 * qt.Jx)
        end
        fidelities[index] = timeevolution.schroedinger_dynamic([0.0, final_time], qt.ψ0, H_eSTA; fout=fidelity)[2][end]  # Time evolution where the output is not the resulting state but the fidelity. It helps improving the speed of the calculation
    end
    return fidelities
end

