import numpy as np
import orb_mech_constants as omc

class OrbitalElements:
    def __init__(self,a,e,i,Omega,omega,nu):
        self.a = a
        self.e = e
        self.i = i
        self.Omega = Omega
        self.omega = omega
        self.nu = nu

class StateVector:
    def __init__(self,R: np.ndarray,V: np.ndarray):
        self.R = R
        self.V = V

def orbital_elements_from_state_vector(state_vector: StateVector):
    '''
    Returns orbital elements given state vector.
    '''
    R = state_vector.R
    V = state_vector.V

    h = np.cross(R,V)
    n = np.cross([0,0,1],h)
    e = (np.cross(V,h) - omc.MU_EARTH * R / np.linalg.norm(R)) / omc.MU_EARTH
    a = 1 / (2 / np.linalg.norm(R) - np.dot(V,V) / omc.MU_EARTH)
    i = np.arccos(h[2] / np.linalg.norm(h))
    if n[1] >= 0:
        Omega = np.arccos(n[0] / np.linalg.norm(n))
    else:
        Omega = 2*np.pi - np.arccos(n[0] / np.linalg.norm(n))
    if e[2] >= 0:
        omega = np.arccos(np.dot(n,e) / (np.linalg.norm(n) * np.linalg.norm(e)))
    else:
        omega = 2*np.pi - np.arccos(np.dot(n,e) / (np.linalg.norm(n) * np.linalg.norm(e)))
    if np.dot(R,V) >= 0:
        nu = np.arccos(np.dot(e,R) / (np.linalg.norm(e) * np.linalg.norm(R)))
    else:
        nu = 2*np.pi - np.arccos(np.dot(e,R) / (np.linalg.norm(e) * np.linalg.norm(R)))

    return OrbitalElements(a,np.linalg.norm(e),i,Omega,omega,nu)

def state_vector_from_orbital_elements(orbital_elements: OrbitalElements):
    '''
    Returns state vector given orbital elements.
    '''
    a = orbital_elements.a
    e = orbital_elements.e
    i = orbital_elements.i
    Omega = orbital_elements.Omega
    omega = orbital_elements.omega
    nu = orbital_elements.nu

    p = a * (1 - e**2)
    R = p / (1 + e * np.cos(nu)) * np.array([np.cos(nu),np.sin(nu),0])
    V = np.sqrt(omc.MU_EARTH / p) * np.array([-np.sin(nu),e + np.cos(nu),0])

    # TODO: Find perifocal to ECI rotation matrix

    return StateVector(R,V)


def gibbs_solver(R1,R2,R3):
    '''
    Returns orbital elements given 3 position vectors.
    '''
    