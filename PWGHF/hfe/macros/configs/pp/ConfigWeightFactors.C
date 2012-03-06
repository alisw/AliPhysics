void ConfigWeightFactors(AliAnalysisTaskHFE *task, Bool_t syst = kFALSE, Bool_t ispp = 1){
  //
  // Set weighting factors for nonHFE backgrounds
  //
  
  //Get the correction factors for Non-HF electron yields from a root-file
  Double_t elecBackGroundWeight[11][6][44][3];//centrality, species, momentum, background level
  for(Int_t iCent = 0; iCent < 11; iCent++){
    for(Int_t iSpecies = 0; iSpecies < 6; iSpecies++){
      for(Int_t iBin = 0; iBin < 44; iBin++){
        for(Int_t iError = 0; iError < 3; iError++){
          elecBackGroundWeight[iCent][iSpecies][iBin][iError] = 0;
        }
      }
    }
  }
  const Char_t *backNameMC[6] = {"pion","eta","omega","phi","etap","rho"};
  TFile *weightFile = TFile::Open(Form("%s/util/hfe/rootfiles/nonHFEcorrect.root", gSystem->Getenv("TRAIN_ROOT")));
  if(weightFile){
    if(syst){
      TH1F *hRelErr[2];//errors for pion yields, which form the correlated component of the relative error for all other decaying mesons
      hRelErr[0] = (TH1F*)weightFile->Get("hErrorspionLower");
      hRelErr[1] = (TH1F*)weightFile->Get("hErrorspionUpper");
    }
    for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
      for(Int_t iSpecies = 0; iSpecies < 6; iSpecies++){//species of decaying mesons
        TH1F *hRatio = 0x0;
        if(ispp){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio%s",backNameMC[iSpecies]));
        }
        else{
          if((iCent == 1)||(iCent == 4)){ 
            hRatio = (TH1F*)weightFile->Get(Form("hRatio%s%d",backNameMC[iSpecies],iCent-1));
          }
          else if(iCent > 7){
            hRatio = weightFile->Get(Form("hRatio%s7",backNameMC[iSpecies]));
          }
          else{
            hRatio = (TH1F*)weightFile->Get(Form("hRatio%s%d",backNameMC[iSpecies],iCent));
          }
        }
	for(Int_t iBin = 1; iBin < 45; iBin++){//momentum bin of mother meson
	  if(iCent == 0){
	    elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);
            if(syst && ispp){
              for(Int_t iError = 0; iError < 2; iError++){//0: best estimate, 1,2: lower, upper uncertainty level
                elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[iError]->GetBinContent(iBin));//Addition of relative errors from histograms with "+", because lower errors are defined as negative numbers in the reference histograms!
	      }
	    }
	  }
          else{
            if(!ispp){
              elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);             
            }           
          }
	}
      }
    }
    weightFile->Close();
  }
  else{
    printf("No reference file for weighting found!\n");   
  }
  
  const Double_t binLimit[45] =  {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum - possibly we won't need the last bin
  
  for(Int_t iBin = 1; iBin < 45; iBin++){//for all pt bins and all meson decays, set weighting factors for daughter electrons
    task->SetBinLimits(iBin,binLimit[iBin-1]);
    for(Int_t iSpecies = 0; iSpecies < 6; iSpecies++){
      for(Int_t iError = 0; iError < 3; iError++)
	task->SetElecBackGroundFactors(iBin, iSpecies, 0, iError, elecBackGroundWeight[0][iSpecies][iBin-1][iError]);
    }
  }
}
