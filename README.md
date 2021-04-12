# Rademaker_etal_BilayerFLEX-ME_2021
Source code and data files for the manuscript "Enhanced superconductivity in FeSe/SrTiO3 from the combination of forward scattering phonons and spin fluctuations". 

Reference: Physical Review B 103, 144504 (2021). (doi:10.1103/PhysRevB.103.144504)
Preprint: arXiv:2101.08307 (2021). (https://arxiv.org/abs/2101.08307)

Title: Enhanced superconductivity in FeSe/SrTiO3 from the combination of forward scattering phonons and spin fluctuations

Authors: Louk Rademaker, Gustavo Alvarez-Suchini, Ken Nakatsukasa, Yan Wang, and Steven Johnston

Abstract: We study the effect of combining spin fluctuations and forward scattering electron-phonon (e- ph) coupling on the superconductivity in the FeSe/SrTiO3 system modeled by a phenomenological two-band Hubbard model with long-range e-ph interactions. We treat the electron and phonon degrees of freedom on an equal footing using a fully self-consistent FLEX plus Migdal-Eliashberg calculation, which includes a self-consistent determination of the spin fluctuation spectrum. Based on FeSe monolayers, we focus on the case where one of the bands lies below the Fermi level (i.e. incipient), and demonstrate that the combined interactions can enhance or suppress 𝑇𝑐, depending on their relative strength. For a suitable choice of parameters, the spin-fluctuation mechanism yields a 𝑇𝑐 ≈ 46.8 K incipient 𝑠± superconductor, consistent with surface-doped FeSe thin films. A forward-focused e-ph interaction further enhances the 𝑇𝑐, as observed in monolayer FeSe on SrTiO3.

Contents:
- The folder "Matlab Source_Code" contains the version of the two-orbital FLEX+Forward Scattering code used to produce the results in the paper. Run the two_orb_FLEX.m script to obtain the data.

- "Figure_3_tz_table.pdf" lists all of the t_\perp values used for each value of U throughout the paper. (The same value of U used for lambda_m = 0 is also used for lambda_m != 0.)

- The folders "Figure X", where X = 1-5 contain all scripts and data files needed to reproduce the figures in the paper. 

- "Notes on the FLEX formalism.pdf" contains general notes about our implementation. 
