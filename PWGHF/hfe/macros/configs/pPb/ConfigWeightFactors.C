void ConfigWeightFactors(AliAnalysisTaskHFE *task, Bool_t syst = kFALSE, Int_t collType = 1, TString filename = "nonHFEcorrect.root"){
  //
  // Set weighting factors for nonHFE backgrounds
  // Option "collType": 0 for pp 2.76 TeV; 1 for pp 7 TeV; 2 for PbPb; 3 for DPMJET pPb; 4 for HIJING pPb; 5 for DPMJET/HIJING pPb
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
  printf("Take the weights from %s\n",Form("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
  TFile *weightFile = TFile::Open(Form("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
  if(weightFile){
    if(syst){
      TH1F *hRelErr[2][2];//errors for pion yields, which form the correlated component of the relative error for all other decaying mesons, except for eta, which are parameterized independently
      if(collType == 1){
        hRelErr[0][0] = (TH1F*)weightFile->Get("hErrorspionLower");
        hRelErr[0][1] = (TH1F*)weightFile->Get("hErrorspionUpper");
        hRelErr[1][0] = (TH1F*)weightFile->Get("hErrorsetaLower");
        hRelErr[1][1] = (TH1F*)weightFile->Get("hErrorsetaUpper");
      }
      else if(collType == 0){
        hRelErr[0][0] = (TH1F*)weightFile->Get("hErrors_2.76TeV_pionLower");
        hRelErr[0][1] = (TH1F*)weightFile->Get("hErrors_2.76TeV_pionUpper");
        //hRelErr[1][0] = (TH1F*)weightFile->Get("hErrors_2.76TeV_etaLower");
        //hRelErr[1][1] = (TH1F*)weightFile->Get("hErrors_2.76TeV_etaUpper");
      }
      else if(collType == 3 || collType == 4 || collType == 5){
        hRelErr[0][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_pionLower");
        hRelErr[0][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_pionUpper");
        for(int i=0; i<hRelErr[0][0]->GetNbinsX(); i++){ //assign 7% systematic uncertainties for pPb
          hRelErr[0][0]->SetBinContent(i+1, -0.07);
          hRelErr[0][1]->SetBinContent(i+1, 0.07);
        }
      }
    }
    for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
      for(Int_t iSpecies = 0; iSpecies < 6; iSpecies++){//species of decaying mesons
        TH1F *hRatio = 0x0;
        if(collType == 1){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio%s",backNameMC[iSpecies]));
        }
        else if(collType == 0){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_2.76TeV_%s",backNameMC[iSpecies]));
        }
        else if(collType == 3){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_DPMJET_%s",backNameMC[iSpecies]));
        }
        else if(collType == 4){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_%s",backNameMC[iSpecies]));
        }
        else if(collType == 5){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_DvsH_%s",backNameMC[iSpecies]));
        }
        else{
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
	for(Int_t iBin = 1; iBin < 45; iBin++){//momentum bin of mother meson
	  if(iCent == 0){
 	    elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);
            if(syst && (collType != 2)){
              for(Int_t iError = 0; iError < 2; iError++){//0: best estimate, 1,2: lower, upper uncertainty level
                if((iSpecies == 1) && (collType == 1))
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[1][iError]->GetBinContent(iBin));
                else
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[0][iError]->GetBinContent(iBin));//Addition of relative errors from histograms with "+", because lower errors are defined as negative numbers in the reference histograms!
	      }
	    }
	  }
          else{
            if(collType == 2){
              elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);    
            }           
          }
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
      for(Int_t iSpecies = 0; iSpecies < 6; iSpecies++){
        for(Int_t iError = 0; iError < 3; iError++)
          task->SetElecBackGroundFactors(iBin-1, iSpecies, iCent, iError, elecBackGroundWeight[iCent][iSpecies][iBin-1][iError]);
      }
    }
  }
}
