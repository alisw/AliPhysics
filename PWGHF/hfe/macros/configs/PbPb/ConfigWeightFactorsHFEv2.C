void ConfigWeightFactorsHFEv2(AliAnalysisTaskFlowTPCTOFEPSP *task, Bool_t syst = kFALSE, Int_t collType = 1, TString filename = "nonHFEcorrect.root"){
  //
  // Set weighting factors for nonHFE backgrounds
  // Option "Type": 0 for default PbPb (old weights for Hijing); 1 for Pb-Pb LHC11a10abis; 2 Pb-Pb LHC11a10b_plus; 3 Pb-Pb LHC11a10b_plus for LHC11a10abis ; 4 Pb-Pb LHC11a10b_plus for LHC11a10bbis; 
  //
  //Get the correction factors for Non-HF electron yields from a root-file
  const int nSpec = 9;
  Double_t elecBackGroundWeight[11][9][44][3];//centrality, species, momentum, background level
  for(Int_t iCent = 0; iCent < 11; iCent++){
    for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){
      for(Int_t iBin = 0; iBin < 44; iBin++){
        for(Int_t iError = 0; iError < 3; iError++){
          elecBackGroundWeight[iCent][iSpecies][iBin][iError] = 0;
        }
      }
    }
  }
  const Char_t *backNameMC[9] = {"pion","eta","omega","phi","etap","rho","kaon","k0s","lambda"};
  printf("Take the weights from %s\n",Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s",filename.Data()));
  //printf("Take the weights from %s\n",Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  printf("collType %d\n",collType);
  TFile *weightFile = TFile::Open(Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  //TFile *weightFile = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s",filename.Data()));
  if(weightFile){
    //weightFile->ls();
    if(syst){
      // Not implemented
    }
    for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){//species of decaying mesons
        TH1F *hRatio = 0x0;
        if(collType == 0){
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = (TH1F*)weightFile->Get(Form("hRatio%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio%s%d",backNameMC[iSpecies],iCent));
          }
        }
	else if(collType == 2){
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10abis_%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10abis_%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10abis_%s%d",backNameMC[iSpecies],iCent));
          }
        }
	else if(collType == 3){
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bplus_%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bplus_%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bplus_%s%d",backNameMC[iSpecies],iCent));
          }
        }
	else if(collType == 4){
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_bplusvasabis_%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_bplusvasabis_%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_bplusvasabis_%s%d",backNameMC[iSpecies],iCent));
          }
        }
	else {
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bbis_%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bbis_%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio_PbPb_2.76TeV_LHC11a10bbis_%s%d",backNameMC[iSpecies],iCent));
          }
        }
	for(Int_t iBin = 1; iBin < 45; iBin++){//momentum bin of mother meson
	  elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);
	}
      }
    }
    weightFile->Close();
  }
  else{
    printf("No reference file for background electron weighting found!\n");   
  }
  
  const Double_t binLimit[45] =  {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
  
  for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
    for(Int_t iBin = 1; iBin < 45; iBin++){//for all centralities, pt bins and all meson decays, set weighting factors for daughter electrons
      task->SetBinLimits(iBin-1,binLimit[iBin-1]);
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){
        for(Int_t iError = 0; iError < 3; iError++)
          task->SetElecBackGroundFactors(iBin-1, iSpecies, iCent, iError, elecBackGroundWeight[iCent][iSpecies][iBin-1][iError]);
      }
    }
  }
}
