Neutral meson analysis in $pp$ at 13 TeV  
=========================================
Main code to analyse LHC16 proton-proton data. The code should support both AOD and ESD format of input data.

## Usage
This analysis task shoul alwaysd be used with [phos tender](https://github.com/alisw/AliPhysics/tree/master/PWGGA/PHOSTasksPHOS_PbPb/AddAODPHOSTender.C). This is important as tender handles calibrations, geometry etc. Just add the following lines in your analysis code:

```c++
// In your run.C macro

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
├── AliAnalysisTaskPP13
├── AliPP13::ClusterCuts
├── AliPP13::DetectorHistogram
├── AliPP13::MixingSample
├── AliPP13::ParticlesHistogram
└── AliPP13::PhotonSelection
    ├── AliPP13::PhotonSpectrumSelection
    ├── AliPP13::PhotonTimecutStudySelection
    ├── AliPP13::PythiaInfoSelection
    ├── AliPP13::QualityPhotonSelection
    ├── AliPP13::TagAndProbeSelection       // TODO: Check final parametrization
    ├── AliPP13::EfficiencySelectionMC
    └── AliPP13::PhysPhotonSelection
        ├── AliPP13::NonlinearityScanSelection
        └── AliPP13::PhysPhotonSelectionMC
            └── AliPP13::MesonSelectionMC
```
