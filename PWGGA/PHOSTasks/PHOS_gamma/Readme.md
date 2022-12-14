Codes for analysis of photons in proton-proton and lead-lead collisions

AliAnalysisTaskGammaPHOSPP  photon  analysis in proton-proton collisions, Runs 1 and 2, fills QA, one and two photon and track histograms. It is possible to analyze MC data as well.

AliAnalysisTaskPHOSPbPbQARun2 fills QA histograms for Pb-Pb events in Run2

Example how to load AliAnalysisTaskGammaPHOSPP:

    TString tenderOption = isMC ? "Run2NoNcellCut" : "Run2Tune";
    Int_t   tenderPass   = 1;
    TString nonlinearity = isMC ? "Run2TuneMCNoNcell" : "Run2Tune";
     
    AliAnalysisTaskGammaPHOSPP *task;

    task = reinterpret_cast<AliAnalysisTaskGammaPHOSPP*>((gInterpreter->ExecuteMacro(
        Form("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_gamma/AddTaskGammaPHOSPP.C(%d, \"%s\", %d, \"%s\", %d)", 
               isMC, tenderOption.Data(), tenderPass, nonlinearity.Data(), 1))));
