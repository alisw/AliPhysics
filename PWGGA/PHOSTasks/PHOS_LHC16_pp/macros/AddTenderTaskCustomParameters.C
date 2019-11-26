#ifdef __CLING__
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C>
#endif

AliPHOSTenderTask * AddTenderTaskCustomParameters(Bool_t isMC = kFALSE, TString badmap = "BadMap_LHC16-updated.root")
{
    // Copy necessary map from the private directory
    //
    gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/o/okovalen/maps/%s .",badmap.Data()));

    // Restore the analysis manager
    // 
    AliAnalysisManager * manager = AliAnalysisManager::GetAnalysisManager();

    // Setup PHOS Tender
    //
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");

    TString tenderOption = "Run2Default";
    AliPHOSTenderTask * tenderPHOS = AddAODPHOSTender("PHOSTenderTask", "PHOStender", tenderOption, 1, isMC);
    AliPHOSTenderSupply * PHOSSupply = tenderPHOS->GetPHOSTenderSupply();
    // IMPORTANT: Set the map of bad channels
    PHOSSupply->ForceUsingBadMap(badmap.Data());


    if (isMC)
    {
        // Important: Keep track of this variable
        // ZS threshold in unit of GeV
        Double_t zs_threshold = 0.020;
        PHOSSupply->ApplyZeroSuppression(zs_threshold);
        PHOSSupply->SetNonlinearityVersion("Run2TuneMC");
    }

    return tenderPHOS;
}
