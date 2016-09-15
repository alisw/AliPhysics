// e.g. a 'runExtractFFs.C("finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/Jets/nclCut/noMCidForGen/bhess_PID_Jets.root", "", "finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/Jets/nclCut/noMCidForGen/bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra_jetPt5.0_10.0.root", 0, -2, -2, 5, 10, 1, 2, kFALSE)' -b -q
//___________________________________________________________________
Int_t runExtractFFs(TString pathNameData, TString listName /*= ""*/, TString pathNameFractionsAndYields,
                    Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                    Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                    Double_t lowerJetPt /* = -1*/, Double_t upperJetPt /*= -1*/,
                    Int_t rebinZ /*= 1*/, Int_t rebinXi /*= 1*/, Bool_t onlyUseRelevantMCIDforMatrix /*=kFALSE*/)
{
  gROOT->LoadMacro("FFs/bhess_PID/AliAnalysisTaskPIDV0base.cxx+g");
  gROOT->LoadMacro("FFs/bhess_PID/AliAnalysisTaskPID.cxx+g");
  
  gROOT->LoadMacro("extractFFs.C+");
  
  return extractFFs(pathNameData, listName, pathNameFractionsAndYields, chargeMode, lowerCentrality, upperCentrality, lowerJetPt,
                    upperJetPt, rebinZ, rebinXi, onlyUseRelevantMCIDforMatrix);
}