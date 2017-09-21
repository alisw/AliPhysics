Neutral meson analysis in $pp$ at 13 TeV  
=========================================
Main code to analyse LHC16 proton-proton data. The code should support both AOD and ESD format of input data.

## Usage
This analysis task shoul alwaysd be used with [phos tender](https://github.com/alisw/AliPhysics/tree/master/PWGGA/PHOSTasksPHOS_PbPb/AddAODPHOSTender.C). This is important as tender handles calibrations, geometry etc. Just add the following lines in your analysis code:

```c++
// In your run.C macro

gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
AliPHOSTenderTask * tenderPHOS = AddAODPHOSTender("PHOSTenderTask", "PHOStender", tenderOption, 1, isMC);
AliPHOSTenderSupply * PHOSSupply = tenderPHOS->GetPHOSTenderSupply();
PHOSSupply->ForceUsingBadMap("/path/to/badmap/BadMap_LHC16.root");

if(isMC)
{
    // Use custom Zero Suppression threshold if needed
    Double_t zs_threshold = 0.020;
    PHOSSupply->ApplyZeroSuppression(zs_threshold); 
}

// If you don't strong cuts on clusters, leave 
// this list of bad cells empty.
std::vector<Int_t> cells;

TString files = AddAnalysisTaskPP(AliVEvent::kINT7, // Physics selection 
    "Changed ZS threshold to 20 MeV",               // Comments, tags to distinguish output *.root files
    "SelectionPrefix",                              // Prefix all your output if you want to use add multiple AliAnalysisTaskPPs
    "",                                             // Path to Bad Map, don't use it if you have set one in Tender
    cells,                                          // List of bad cells
    isMC);

// Add files to the grid plugin
//
alienHandler->SetOutputFiles(files);
```

## Class hierarchy

```c++
AliAnalysisTaskPP            //  Main analysis 
└───   MixingSample
└───   PhotonSelection         
        └───   PhysPhotonSelection   // Fills all histograms needed to reconstruct $\pi^{0}$s
        │           DetectorHistogram
        │
        └───   PhotonTimecutSelection  // Study of timing cut efficiency and purity
        └───   QualityPhotonSelection  // QA plots for PHOS clusters
        └───   TagAndProbeSelection
        └───   PhysPhotonSelectionMC   // Applies nonlinearity to PHOS clusters in MC
        └───   MesonSelectionMC        // Neutral meson efficiency study
        │           ParticlesHistogram
        │
        └───   PythiaInfoSelection     // Collects cross section and ntrials data. Needed for jet-jet MC only.

```
