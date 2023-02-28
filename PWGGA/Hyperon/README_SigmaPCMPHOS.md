# AliAnalysisTaskSigmaPCMPHOS

---

## Description

The following code is used to analyze 13 TeV data to extract $\Sigma^0$ hyperon signal via its electronmagnetic decay $\Sigma^0\to \Lambda \gamma$. 3 files are required:

-  Header (**AliAnalysisTaskSigmaPCMPHOS.h**) file containing class methods' prototypes required for the analysis, see this file for explanation of used variables;
    
- Methods' implementation (**AliAnalysisTaskSigmaPCMPHOS.cxx**) files corresponding to the header files;
    
- Class instance (**AddAnalysisTaskSigmaPCMPHOS.C**) macro;

## Analysis Layout

After the setup is completed, Monte-Carlo events are processed if the simulation flag is true. Photon verteces are found afterwards. Then $\Lambda$ hyperons and $\gamma$ are reconstruted. Finally, $\gamma \gamma$ and $\gamma \Lambda$ invariant mass vs. transverse momentum histograms are filled using photons from the photon conversion method and the photon spectrometer of ALICE.