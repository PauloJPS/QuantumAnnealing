import numpy as np 
import qutip as qt 

def Grover_Hi_Hf(m, N):
    """
        Receive the number of qubits and return the initial Hamiltonian for the Grover problem
        ---------
        Parameters:
            N : (Number of qubits - int)
            m : (index of the market state - int)
        Return:
            Hi : A qutip object representing the initial Hamiltonian for the Grover problem 
            Hf : A qutip object representing the final Hamiltonian for the Grover problem 

    """
    I = GetIdentity(N)
    eigenstates = I.eigenstates()[1]
    psi = np.sqrt(1/2**N) * np.sum(eigenstates)
    ket_m = eigenstates[m]

    Hi = I - psi * psi.dag()
    Hf = I - ket_m * ket_m.dag() 
    return Hi, Hf

def GetIdentity(N):
    """
        Receive the number of qubits and return the Identity matrix for its Hilbert space 
        This is necessary because we want to keep the tensorial composition form 
        Parameters:
            N : (Number of qubits - int)
        Return:
            A qutip object representing a Identity matrix (2^N X 2^N)
    """
    I = np.zeros(N, dtype=object)
    for i in range(N):
        I[i] = qt.identity(2) 
    I = qt.tensor(I)
    return I 

def Spectrum_alongInterpolation(Hi, Hf, psii, As, Bs, tlist):
    """
        Return the spectrum of the time-dependent Hamiltonian for each time interval
        Parameters: 
            Hi : Initial Hamiltonian (qutip object)
            Hf : Final Hamiltonian (qutip object)
            As : Parametrization function (python function)
            Bs : Parametrization function (python function)
            tlist : time list (np.array)
        Return: 
            list of eigenenergies for each time interval in tlist 

    """
    args = {'t_max': tlist[-1]} 
    return np.array([(As(i, args)*Hi + Bs(i, args)*Hf).eigenenergies() for i in tlist])
    

def QuantumAnnealing(Hi, Hf, psii, As, Bs, tlist):
    """
        Return the system states along the interpolation for a given parametrization and a time list    
        Parameters: 
            Hi : Initial Hamiltonian (qutip object)
            Hf : Final Hamiltonian (qutip object)
            psii : Initial State (qutip object)
            As : Parametrization function (python function)
            Bs : Parametrization function (python function)
            tlist : time list (np.array)
        Return: 
            system state in each interpolation time interval 
    """
    Ht = [[Hi, As], [Hf, Bs]]
      
    args = {'t_max': tlist[-1]} 
    result = qt.sesolve(H=Ht, psi0=psii, tlist=tlist,  args=args)
    return result.states 

def GetFidelity(Hi, Hf, psii, psif, As, Bs, tlist):
    """
        Return the system states along the interpolation for a given parametrization and a time list    
        Parameters: 
            Hi : Initial Hamiltonian (qutip object)
            Hf : Final Hamiltonian (qutip object)
            psii : Initial State (qutip object)
            psif : Expected Final State (qutip object)
            As : Parametrization function (python function)
            Bs : Parametrization function (python function)
            tlist : time list (np.array)

        Return:
            The evolution fidelity related to psif 
    """
    psif_ev = QuantumAnnealing(Hi, Hf, psii, As, Bs, tlist)[-1]

    return np.abs((psif.dag() * psif_ev * psif_ev.dag() * psif)[0][0][0])

def Fidelity_X_tmax(Hi, Hf, psii, psif, As, Bs, tmax_list, tlist_nints):
    """
        Return the Fidelity in function of the evolution total time 
        Parameters: 
            Hi : Initial Hamiltonian (qutip object)
            Hf : Final Hamiltonian (qutip object)
            psii : Initial State (qutip object)
            psif : Expected Final State (qutip object)
            As : Parametrization function (python function)
            Bs : Parametrization function (python function)
            tmax_list : total time list (np.array)
            tlist_nints : number of intervals of the evolution time 
        Return:
            fidelity : np.array (len(fidelity) == len(tmax_list))
    """
    N = len(tmax_list)
    fidelity = np.zeros(N)
    for i in range(N):
        tlist = np.linspace(1, tmax_list[i], tlist_nints)
        fidelity[i] = GetFidelity(Hi, Hf, psii, psif, As, Bs, tlist) 
    return fidelity


