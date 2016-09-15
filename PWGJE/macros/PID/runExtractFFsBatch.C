//___________________________________________________________________
Int_t runExtractFFsBatch(Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                         Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                         Int_t rebinZ /*= 1*/, Int_t rebinXi /*= 1*/, Bool_t onlyUseRelevantMCIDforMatrix /*=kFALSE*/)
{
  gROOT->LoadMacro("../trunk/bhess_PID/AliAnalysisTaskPIDV0base.cxx+g");
  gROOT->LoadMacro("../trunk/bhess_PID/AliAnalysisTaskPID.cxx+g");
  
  gROOT->LoadMacro("extractFFs.C+");
  
  const Int_t numCentralities = 9;
  const Int_t centralities[numCentralities + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 80, 100 };
  
  TString centralityString = "";
  
  const Int_t numJetPtBins = 4;
  const Int_t jetPt[numJetPtBins + 3] = { 5, 10, 15, 20, 30, 10, 40 };
  
  TString jetPtString = "";
  
  
  // pp Jets
  ///*
  const TString path = "/hera/alice/bhess/analysis/10d_e.pass2_merged/nclCut";
  Bool_t useCentralities = kFALSE;
  
  const TString rawOutputPathName = "10d_e.pass2_merged/nclCut/bhess_PID_Jets.root";
  const TString listName = "";
  
  for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
    if (iJetPt == numJetPtBins) {
      // 10-40 bin
      iJetPt++;
    }
    jetPtString = Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]);
    
    {
      // Default
      TString fileName = 
        Form("%s/bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Pre-PID
      TString fileName =
        Form("%s/bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Shape
      TString fileName =
        Form("%s/bhess_PID_Jets_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Sigma
      TString fileName =
        Form("%s/bhess_PID_Jets_SystematicsSigmaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
      
      fileName =
        Form("%s/bhess_PID_Jets_SystematicsSigmaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Eta
      TString fileName =
        Form("%s/bhess_PID_Jets_SystematicsEtaDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
      
      fileName =
        Form("%s/bhess_PID_Jets_SystematicsEtaUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Splines
      TString fileName =
        Form("%s/bhess_PID_Jets_SystematicsSplinesDown_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
      
      fileName =
        Form("%s/bhess_PID_Jets_SystematicsSplinesUp_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
    {
      // Muons
      TString fileName =
        Form("%s/bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonsEqualElectrons_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi);
    }
  }
  //*/
  
  // MC_pp jets, default
  /*
  const TString path = "/hera/alice/bhess/analysis/10f6a/nclCut";
  Bool_t useCentralities = kFALSE;
  
  const TString rawOutputPathName = "10f6a/nclCut/bhess_PID_Jets.root";
  const TString listName = "";
  
  for (Int_t iJetPt = 0; iJetPt <= numJetPtBins; iJetPt++) {
    if (iJetPt == numJetPtBins) {
      // 10-40 bin
      iJetPt++;
    }
    jetPtString = Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]);
    
    {
      // Default
      TString fileName = 
        Form("%s/bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi, onlyUseRelevantMCIDforMatrix);
    }
    {
      // Shape
      TString fileName =
        Form("%s/bhess_PID_Jets_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra%s.root",
             path.Data(), jetPtString.Data());
      
      extractFFs(rawOutputPathName, listName, fileName, chargeMode, lowerCentrality, upperCentrality, jetPt[iJetPt], jetPt[iJetPt + 1],
                 rebinZ, rebinXi, onlyUseRelevantMCIDforMatrix);
    }
  }
  */
}