
import numpy as np
import numpy.random as random
from tqdm import tqdm

def calc_source_E(lattice: np.ndarray, g: float) -> float:
    """Calculate the source term for the Wolff algorithm

    Args:
        lattice (np.ndarray): The lattice of spins
        J (float): The coupling constant for nearest neighbours
        g (float): The anisotropy constant (external field)

    Returns:
        float: The source term
    """
    return  g*np.sum(lattice[:, 2]**2 - 0.5*(lattice[:, 0]**2 + lattice[:, 1]**2))

class Lattice():
    def __init__(self, N: int, rng, start: str="hot"):            
        self.N = N
        self.rng = rng
        if start == "cold":
            self.spins = np.stack([np.zeros((self.N)), 
                                   np.zeros((self.N)),
                                   np.ones((self.N))]).T
        elif start == "hot":
            self.spins = self.rng.random((self.N, 3))
        else:
            raise ValueError("start must be either 'hot' or 'cold'")
        self.spins = self.spins/np.linalg.norm(self.spins, axis=1)[:, None]


    def __len__(self):
        return self.N


    def __repr__(self):
        return f"Lattice of size {self.N}\n with spins:\n{self.spins}"


    def calc_energy(self, J: float, g: float) -> float:
        """Calculate the action of the lattice

        Args:
            J (float): The coupling constant for nearest neighbours
            g (float): The anisotropy constant (external field)

        Returns:
            float: The action
        """
        E = -J*np.sum(np.sum(np.roll(self.spins, -1, axis=0)*self.spins, axis=1))
        E += g*np.sum(
            self.spins[:, 2]**2 - 0.5*(self.spins[:, 0]**2 + self.spins[:, 1]**2))
        return E/self.N


    def dH(self, idx: int, J: float, g: float, beta: float) -> float:
        """Compute the change in energy when considering flipping the spin at
        index idx w.r.t to a randomly chosen unit vector (as shown in
        https://www.physik.uni-leipzig.de/~janke/Paper/spinmc.pdf)

        Args:
            idx (int): The index of the spin to flip
            J (float): The coupling constant for nearest neighbours
            g (float): The anisotropy constant (external field)
            beta (float): The inverse temperature

        Returns:
            float: The change in action
        """
        r = self.rng.random((3,))
        r = r/np.linalg.norm(r)
        s = self.spins[idx]
        sprev = self.spins[(idx-1)%self.N]
        snext = self.spins[(idx+1)%self.N]
        sparallel = np.dot(s, r)*r
        sperpendicular = s-sparallel
        s_tilda = sperpendicular - sparallel
        dH = np.sum(-J*(sprev+snext)*s_tilda)\
            +g*(s_tilda[2]**2-0.5*(s_tilda[0]**2+s_tilda[1]**2))\
            -np.sum(-J*(sprev+snext)*s\
            +g*(s[2]**2-0.5*(s[0]**2+s[1]**2)))
        
        return dH
    

    def MC_step(self, idx: int, J: float, g: float, beta: float) -> None:
        """Perform a single Monte Carlo step

        Args:
            idx (int): The index of the spin to flip
            J (float): The coupling constant for nearest neighbours
            g (float): The anisotropy constant (external field)
            beta (float): The inverse temperature

        Returns:
            None: The spin at index idx is updated internally
        """
        dH = self.dH(idx, J, g, beta)
        chi = self.rng.random()
        if chi < np.exp(-beta*dH): # accept
            self.spins[idx] = -self.spins[idx] # updated internally
        return None


    def makeBonds(self, blob: np.ndarray, sign_lattice: np.ndarray,
                  memoW: np.ndarray, beta: float, ) -> np.ndarray:
        P_bond = 1-np.exp(-2*beta)
        expand = np.array([])
        for idx in blob:
            idx = idx.astype(int)
            prev, next = (idx-1)%self.N, (idx+1)%self.N
            if sign_lattice[prev] == sign_lattice[idx]:
                if self.rng.random()< P_bond:
                    if not memoW[prev]: 
                        if prev not in blob:
                            expand = np.concatenate([expand, np.array([prev])])
            if sign_lattice[next] == sign_lattice[idx]:
                if self.rng.random()< P_bond:
                    if not memoW[next]: 
                        if next not in blob:
                            expand = np.concatenate([expand, np.array([next])])
        return expand
            
            
            
    
    def Wolff_update(self, idx: int, J: float, g: float, beta: float) -> None:
        """Perform a single Monte Carlo step using the Wolff algorithm as shown in 
        http://latt.if.usp.br/cgi-bin-delyra/cntsnd?technical-pages/twawesab/Text.pdf

        Args:
            idx (int): The index of the spin to begin the Wolff algorithm
            J (float): The coupling constant for nearest neighbours
            g (float): The anisotropy constant (external field)
            beta (float): The inverse temperature

        Returns:
            None: The clusters are updated internally
        """
        memoW = np.zeros(self.N, dtype=bool)
        memoW[idx] = True
        blob = np.array([idx])
        while len(blob):
            r = self.rng.random((3,))
            r = r/np.linalg.norm(r)
            
            parallel_component = np.sum(self.spins*r, axis=1)
            perpendicular_component = self.spins - parallel_component[:, None]*r
            sign_lattice = np.sign(parallel_component)
            expand = self.makeBonds(blob, sign_lattice, memoW, beta)
            for idx in blob:
                idx = idx.astype(int)
                sign_lattice[idx] = -sign_lattice[idx]
            newSpins = perpendicular_component \
                + sign_lattice[:, None]*r[None, :]*parallel_component[:, None]
            dS = beta*(calc_source_E(newSpins, g)-calc_source_E(self.spins, g))
            if self.rng.random() < np.exp(-dS):
                self.spins = newSpins
            blob = expand.copy()
            for idx in blob:
                idx = idx.astype(int)
                memoW[idx] = True
        return None
        
    def MC_evolution(self, steps: int, J: float=1, g: float=1,
                     beta: float=1, animate: bool=False, numFrames: int=10) -> dict:
        """Perform Monte Carlo evolution (also known as `sweeps`) according to
        https://www.physik.uni-leipzig.de/~janke/Paper/spinmc.pdf)

        Args:
            steps (int): The number of steps to perform
            J (float): The coupling constant for nearest neighbours
            g (float): The anisotropy constant (external field)
            beta (float): The inverse temperature

        Returns:
            np.ndarray: The action after each step
        """
        actions = np.zeros(steps)
        indices = self.rng.integers(low=0, high=self.N, size=(steps,))
        
        if animate:
            spin_list = []
        
        for i in range(steps):
            newaction = self.calc_energy(J, g)
            actions[i] = newaction
            if animate and i%(steps//numFrames) == 0:
                spin_list.append(self.spins.copy())
            self.Wolff_update(indices[i], J, g, beta)
                
        results = {"actions": actions}
        
        if animate:
            results["spin_list"] = spin_list
        
        return results


def Z(actions):
    return np.mean(np.exp(-actions))

    
    