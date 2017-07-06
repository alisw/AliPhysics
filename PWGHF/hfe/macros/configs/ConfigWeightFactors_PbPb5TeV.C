void ConfigWeightFactors_PbPb5TeV(AliAnalysisTaskHFE *task, Bool_t syst = kFALSE, Int_t collType = 1, TString filename = "nonHFEcorrect_PbPb5TeV.root"){
  //
  // Set weighting factors for nonHFE backgrounds
  // Option "collType": 0 for pp 2.76 TeV; 1 for pp 7 TeV; 2 for PbPb; 3 for DPMJET pPb; 4 for HIJING pPb; 5 for DPMJET/HIJING pPb: 6 Pb-Pb LHC11a10abis; 7 Pb-Pb LHC11a10b_plus; 8 Pb-Pb LHC11a10b_plus for LHC11a10abis ; 9 Pb-Pb LHC11a10b_plus for LHC11a10bbis; 10 for pp LHC106; 11 for pp LHC106a; 12 for pp LHC107a_d
  //
  // 60: PbPb LHC16g1 minimum bias MC (mfaggin, 29/06/2017)
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
  /*
        // GSI version
  printf("Take the weights from %s\n",Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  printf("collType %d\n",collType);
  TFile *weightFile = TFile::Open(Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  */
        // GRID version
    printf("Take the weights from %s\n",Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filename.Data()));
  printf("collType %d\n",collType);
//  TFile *weightFile = TFile::Open(Form("%s/PWGHF/hfe/macros/configs/PbPb/%s", gSystem->Getenv("ALICE_PHYSICS"),filename.Data()));
    TFile *weightFile = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filename.Data()));

  if(weightFile){
    if(syst){
      TH1F *hRelErr[9][2];//errors for pion yields, which form the correlated component of the relative error for all other decaying mesons, except for eta, which are parameterized independently
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
	if(collType == 4){
	  hRelErr[6][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_kaonLower");
	  hRelErr[6][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_kaonUpper");
	  hRelErr[7][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_k0sLower");
	  hRelErr[7][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_k0sUpper");
	  hRelErr[8][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_lambdaLower");
	  hRelErr[8][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_lambdaUpper");
	}
      }
    }

    for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){//species of decaying mesons
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
        else if(collType == 10){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6_%s",backNameMC[iSpecies]));
        }
        else if(collType == 11){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6a_%s",backNameMC[iSpecies]));
        }
        else if(collType == 12){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f7a_d_%s",backNameMC[iSpecies]));
        }
        else if(collType == 13){
	      hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6a_10f7a_d_%s",backNameMC[iSpecies]));
	    }
	else if(collType == 14){
	      hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f7a_d_10f6a_%s",backNameMC[iSpecies]));
	}
	else if(collType == 15){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6_10f6a_%s",backNameMC[iSpecies]));
	}
	else if(collType == 16){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6a_10f6_%s",backNameMC[iSpecies]));
        }
	else if(collType == 17){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6_10f7a_d_%s",backNameMC[iSpecies]));
	}
	else if(collType == 18){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f7a_d_10f6_%s",backNameMC[iSpecies]));
        }
	else if(collType == 19){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_11b10a_%s",backNameMC[iSpecies]));
        }
	else if(collType == 20){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_11b10b_%s",backNameMC[iSpecies]));
        }
	else if(collType == 21){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12a9_%s",backNameMC[iSpecies]));
        }
	else if(collType == 22){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12e6_%s",backNameMC[iSpecies]));
        }
	else if(collType == 23){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12f1a_%s",backNameMC[iSpecies]));
        }
	else if(collType == 24){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12f1b_%s",backNameMC[iSpecies]));
        }
	else if(collType == 25){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12f1a_to_12f1b_%s",backNameMC[iSpecies]));
        }
	else if(collType == 26){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12f1b_to_12f1a_%s",backNameMC[iSpecies]));
        }
	else if(collType == 27){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12e6_tu_%s",backNameMC[iSpecies]));
        }
	else if(collType == 28){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_12e6_td_%s",backNameMC[iSpecies]));
        }
	else if(collType == 29){
	  hRatio = (TH1F*)weightFile->Get(Form("hRatio_11b10b_to_12e6_%s",backNameMC[iSpecies]));
        }

        // 60: PbPb LHC16g1 minimum bias MC (mfaggin, 29/06/2017)
        else if(collType == 60){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_16g1_cp_%s",backNameMC[iSpecies]));
          cout << "\n-------------------------------------------------\n";
          cout << "-------------------------------------------------\n";
          cout << "      PbPb LHC16g1 minimum bias MC weights read      ";
          cout << "\n-------------------------------------------------\n";
          cout << "-------------------------------------------------\n";
        }


        else if(collType == 2){
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
	else if(collType == 6){
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
	else if(collType == 7){
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
	else if(collType == 8){
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
	  if(iCent == 0){
              if(!hRatio){
                  printf("ERROR: Histogram for %s empty!\n",backNameMC[iSpecies]);   
                  printf("set weight to 1\n");
              }
              elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);
            if(syst && (collType != 2)){
              for(Int_t iError = 0; iError < 2; iError++){//0: best estimate, 1,2: lower, upper uncertainty level
                if((iSpecies == 1) && (collType == 1))
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[1][iError]->GetBinContent(iBin));
				else if((iSpecies > 5) && (collType ==4))
				   elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[iSpecies][iError]->GetBinContent(iBin));
				else
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[0][iError]->GetBinContent(iBin));//Addition of relative errors from histograms with "+", because lower errors are defined as negative numbers in the reference histograms!
	      }
	    }
	  }
          else{
            if((collType == 2) || (collType > 5)){
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
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){
        for(Int_t iError = 0; iError < 3; iError++)
          task->SetElecBackGroundFactors(iBin-1, iSpecies, iCent, iError, elecBackGroundWeight[iCent][iSpecies][iBin-1][iError]);
      }
    }
  }
}
void ConfigWeightFactors(AliAnalysisTaskHFEMulti *task, Bool_t syst = kFALSE, Int_t collType = 1, TString filename = "nonHFEcorrect_pPb.root"){
  //
  // Set weighting factors for nonHFE backgrounds
  // Option "collType": 0 for pp 2.76 TeV; 1 for pp 7 TeV; 2 for PbPb; 3 for DPMJET pPb; 4 for HIJING pPb; 5 for DPMJET/HIJING pPb: 6 Pb-Pb LHC11a10abis; 7 Pb-Pb LHC11a10b_plus; 8 Pb-Pb LHC11a10b_plus for LHC11a10abis ; 9 Pb-Pb LHC11a10b_plus for LHC11a10bbis; 10 for pp LHC106; 11 for pp LHC106a; 12 for pp LHC107a_d
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
  printf("Take the weights from %s\n",Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  printf("collType %d\n",collType);
  TFile *weightFile = TFile::Open(Form("%s/util/hfe/%s", gSystem->Getenv("TRAIN_ROOT"),filename.Data()));
  if(weightFile){
    if(syst){
      TH1F *hRelErr[9][2];//errors for pion yields, which form the correlated component of the relative error for all other decaying mesons, except for eta, which are parameterized independently
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
      if(collType == 4){
		  hRelErr[6][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_kaonLower");
		  hRelErr[6][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_kaonUpper");
		  hRelErr[7][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_k0sLower");
		  hRelErr[7][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_k0sUpper");
		  hRelErr[8][0] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_lambdaLower");
		  hRelErr[8][1] = (TH1F*)weightFile->Get("hErrors_pPb_5.023TeV_lambdaUpper");

	  }
	 }
    }

    for(Int_t iCent = 0; iCent < 11; iCent++){//centrality bins
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){//species of decaying mesons
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
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_centMB_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 41){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent0005_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 42){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent0510_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 43){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent1020_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 44){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent2040_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 45){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent4060_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 46){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent6080_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 47){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_HIJING_cent80100_%s",backNameMC[iSpecies]));
          if((iSpecies > 5) && (hRatio==NULL)) hRatio = new TH1F;
        }
        else if(collType == 5){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_pPb_5.023TeV_DvsH_%s",backNameMC[iSpecies]));
        }
        else if(collType == 10){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6_%s",backNameMC[iSpecies]));
        }
        else if(collType == 11){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f6a_%s",backNameMC[iSpecies]));
        }
        else if(collType == 12){
          hRatio = (TH1F*)weightFile->Get(Form("hRatio_10f7a_d_%s",backNameMC[iSpecies]));
        }
        else if(collType == 2){
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
	else if(collType == 6){
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
	else if(collType == 7){
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
	else if(collType == 8){
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
	  if(iCent == 0){
 	    elecBackGroundWeight[iCent][iSpecies][iBin-1][0] = hRatio->GetBinContent(iBin);
            if(syst && (collType != 2)){
              for(Int_t iError = 0; iError < 2; iError++){//0: best estimate, 1,2: lower, upper uncertainty level
                if((iSpecies == 1) && (collType == 1))
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[1][iError]->GetBinContent(iBin));
				else if((iSpecies > 5) && (collType ==4))
				   elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[iSpecies][iError]->GetBinContent(iBin));
				else
                  elecBackGroundWeight[iCent][iSpecies][iBin-1][iError+1]=elecBackGroundWeight[iCent][iSpecies][iBin-1][0]*(1+hRelErr[0][iError]->GetBinContent(iBin));//Addition of relative errors from histograms with "+", because lower errors are defined as negative numbers in the reference histograms!
	      }
	    }
	  }
          else{
            if((collType == 2) || (collType > 5)){
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
      for(Int_t iSpecies = 0; iSpecies < nSpec; iSpecies++){
        for(Int_t iError = 0; iError < 3; iError++)
          task->SetElecBackGroundFactors(iBin-1, iSpecies, iCent, iError, elecBackGroundWeight[iCent][iSpecies][iBin-1][iError]);
      }
    }
  }
}
