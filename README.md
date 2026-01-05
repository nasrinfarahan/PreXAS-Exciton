PreXAS-Exciton: BSE-based code for X-ray absorption spectroscopy of photoexcited materials
===========================================================================================

This repository contains a research code for the ab initio simulation of
optical-pump–x-ray-probe absorption in solids,
explicitly including excitonic effects through the Bethe–Salpeter equation
(BSE).

The implementation is based on an all-electron many-body framework and is
designed to work with GW+BSE calculations performed using the EXCITING code.

The theoretical formalism and representative applications are described in:
N. Farahani and D. Popova-Gorelova
Revealing fingerprints of valence excitons in x-ray absorption spectra with the Bethe–Salpeter equation
Phys. Rev. B 110, 235126 (2024)
https://doi.org/10.1103/PhysRevB.110.235126
------------------------------------------------------------------------------
If you use this code or its theoretical framework, please cite: [our initial
publication](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.235126),

Theoretical Background
------------------------------------------------------------------------------
In an optical-pump–x-ray-probe experiment, an optical pump pulse excites the
system from its ground state into a valence-excited many-body state. A
subsequent x-ray probe pulse induces transitions from deep core levels into the
photoexcited electronic structure. The measured signal corresponds to the x-ray
absorption spectrum of a non-equilibrium excitonic state, rather than that
of the ground state.

**Valence-Excited Initial States**

The optically excited states of the system are described as solutions of the
valence Bethe–Salpeter equation

$$
\hat{H}_{BSE}^{val} |\Psi_{\lambda_{I}} \rangle = E_{\lambda_{I}} |\Psi_{\lambda_{I}} \rangle,               \quad (1)
$$

The corresponding excitonic wavefunctions are expanded in the basis of
electron–hole pairs as

$$
|\Psi_{\lambda_{I}} \rangle = \sum_{cvk} A_{cvk}^{\lambda_{I}} \hat{c}^{\dagger}_{ck} \hat{c}_{vk} |\Psi_{0} \rangle,               \quad (2)
$$

where $`A_{cvk}^{\lambda_{I}} `$ are the valence BSE eigenvectors, 
$`E_{\lambda_{I}} `$  the excitation energies, and  $`|\Psi_{0} \rangle`$ the many-body ground state.


**Core-Excited Final States**

The core-excited states excited by the x-ray probe pulse are obtained from the
core-level Bethe–Salpeter equation

$$
\hat{H}_{BSE}^{core} |\Psi_{\lambda_{F}} \rangle = E_{\lambda_{F}} |\Psi_{\lambda_{F}} \rangle,               \quad (3)
$$

with the excitonic wavefunctions

$$
|\Psi_{\lambda_{F}} \rangle = \sum_{c\mu k} A_{c\mu k}^{\lambda_{F}} \hat{c}^{\dagger}_{ck} \hat{c}_{\mu k} |\Psi_{0} \rangle,               \quad (4)
$$

where $`\mu`$ labels core states and $`c`$ conduction states.

**Pump–Probe X-Ray Absorption Cross Section**


Within second-order time-dependent perturbation theory, the x-ray absorption
cross section of a photoexcited system is given by


$$
\sigma = 
			\frac{1}{c}\sum_{\lambda_F} 
			| \sum_{\lambda_I} w_{\lambda_I} \frac{ \langle \Psi_{\lambda_F}  |\hat \psi^\dagger(r_1) e^{i\mathbf k_x\cdot r_1}(\boldsymbol\epsilon_x\cdot \hat{\mathbf p})\hat \psi (r_1) |\Psi_{\lambda_I}\rangle}{\Gamma_F - i(E_{\lambda_I} - E_{\lambda_F}+\omega_x)} 
			|^2,                \quad (5)
$$

where $`\omega_x`$ is the x-ray probe photon energy, and $`w_{\lambda_I}`$
denotes the occupation (or population weight) of the optically excited state $`\lambda_I`$

$$
w_{\lambda_I} = \sqrt{\frac{2\pi I_{\text{o}}}{c\omega_{\text{o}}^2}}\frac{ \sum_{c,v,k} A^{*\lambda_I}_{cvk}\int d^3 r \phi_{ck}^\dagger(r)( \boldsymbol\epsilon_{\text{o}} \cdot \hat {\mathbf p})  \phi_{vk}(r) }{\Gamma_I - i(E_{0} - E_{\lambda_I}+\omega_{\text{o}})},               \quad (6)
$$





**Transition Matrix Elements**


Evaluating the dipole matrix element between valence- and core-excited states
yields

$$
\langle \Psi_{\lambda_F}  |\hat \psi^\dagger(r_1) e^{i\mathbf k_x\cdot r_1}(\boldsymbol\epsilon_x\cdot \hat{\mathbf p})\hat \psi (r_1) |\Psi_{\lambda_I}\rangle = \sum_{c,v,\mu,k } A^{*\lambda_F}_{c \mu k} A^{\lambda_I}_{c v k} \int d^3 r \phi_{v k}^\dagger(r)e^{i \mathbf k_x\cdot r}( \boldsymbol\epsilon_x \cdot \hat {\mathbf p})  \phi_{\mu k}(r) ,               \quad (7)
$$

where  $`\epsilon_x`$  is the polarization vector of the x-ray probe and the
momentum matrix elements are defined as

$$
P_{\mu v k} = \langle \mu k |\hat{p} | vk \rangle ,              \quad (8)
$$

And subsituating $\quad(7)$ into $\quad(5)$, the cross section becomes

$$
\sigma = 
			\frac{1}{c}\sum_{\lambda_F} 
			| \sum_{\lambda_I} w_{\lambda_I} \frac{\sum_{c,v,\mu,k } A^{*\lambda_F}_{c \mu k} A^{\lambda_I}_{c v k} \int d^3 r \phi_{vk}^\dagger(r)e^{i \mathbf k_x\cdot r}( \boldsymbol\epsilon_x \cdot \hat {\mathbf p})  \phi_{\mu k}(r)}{\Gamma_F - i(E_{\lambda_I} - E_{\lambda_F}+\omega_x)}|^2.              \quad (9)
$$


Code Usage
-----------------------------------------------------------------------------
To use **PreXAS‑Exciton**, you need valence‑ and core‑level BSE data by performing two separate BSE runs in exciting and a simple configuration file that defines the pump–probe setup.

start by performing two separate BSE runs that provide the core‑excited eigenstates $`A_{c \mu \mathbf{k}}^
{\lambda_F}`$  and the valence‑excited eigenstates .
These BSE calculations must share the same kk-point grid, although they may include different numbers of unoccupied (empty) states.

**Ground state and GW calculations**
1. Perform a ground-state DFT calculation with exciting using LAPW+lo basis and GGA-PBE.
2. Run a $G_{0}W_{0}$ calculation on a converged k-grid (converge number of empty bands) to obtain quasiparticle energies for BSE input.


**Valence BSE (optical/pump)**

Run a valence BSE (optical) with exciting to obtain:

  1. Exciton energies $E_{\lambda_{I}}$ and excitonic eigenvectors $`A_{c v \mathbf{k}}^{\lambda_I}`$.
  2. Optical momentum matrix elements $`P_{c v k} = \langle c k |\hat{p} | vk\rangle`$.




**Core BSE (x-ray/probe)**

Run a core-level BSE (x-ray) with exciting to obtain:

  1. Exciton energies $E_{\lambda_{F}}$ and excitonic eigenvectors $`A_{c \mu \mathbf{k}}^{\lambda_F}`$.
  2. Core-conduction momentum matrix elements $`P_{c \mu k} = \langle c k |\hat{p} | \mu k\rangle`$.
 

     Both BSE calculations must use the same k-point grid.

**Modified exciting branch for core-valence matrix elements**


Important: Use a modified version of exciting (your custom build) that computes momentum matrix elements between core and valence states $\langle \mu k |\hat{p} | vk \rangle ,$. 

**Running PreXAS-Exciton**

After performing the BSE calculations, gather the following required files in a single directory along with your *input.cfg* file:
​
