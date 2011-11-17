//fbellini@cern.ch, 11/11/2011
//-----------------------------------------------------------
Int_t MakeTrendingTOFQA(char * runlist, Int_t year=2011, char *period="LHC11d", char* pass="cpass1", Int_t trainId=0){
	Int_t filesCounter=0;
  
	if (!runlist) {
		printf("Invalid list of runs given as input: nothing done\n");
		return 1;
	}	
	Int_t runNumber;
	char infile[300]; 
	char postFileName[20];
  
	char trendFileName[100]; 
	//Define trending output
	if (trainId==0){
		sprintf(trendFileName,"treeTOFQA_%s_%s.root",period,pass);  
	}else{
		sprintf(trendFileName,"treeTOFQA_QA%i_%s_%s.root",trainId,period,pass);
	}
	TFile * trendFile=new TFile(trendFileName,"recreate");
  
	Double_t avTime=-9999., peakTime=-9999., spreadTime=-9999., peakTimeErr=-9999., spreadTimeErr=-9999., negTimeRatio=-9999.,
		avRawTime=-9999., peakRawTime=-9999., spreadRawTime=-9999., peakRawTimeErr=-9999., spreadRawTimeErr=-9999., 
		avTot=-9999., peakTot=-9999.,spreadTot=-9999.,  peakTotErr=-9999.,spreadTotErr=-9999.,
		orphansRatio=-9999., avL=-9999., negLratio=-9999.,
		effPt1=-9999., effPt2=-9999., matchEffLinFit1Gev=-9999.,matchEffLinFit1GevErr=-9999.;
  
	Double_t avPiDiffTime=-9999.,peakPiDiffTime=-9999., spreadPiDiffTime=-9999.,peakPiDiffTimeErr=-9999., spreadPiDiffTimeErr=-9999.;
  
	Double_t avT0A=-9999.,peakT0A=-9999., spreadT0A=-9999.,peakT0AErr=-9999., spreadT0AErr=-9999.;
	Double_t avT0C=-9999.,peakT0C=-9999., spreadT0C=-9999.,peakT0CErr=-9999., spreadT0CErr=-9999.;
	Double_t avT0AC=-9999.,peakT0AC=-9999., spreadT0AC=-9999.,peakT0ACErr=-9999., spreadT0ACErr=-9999.;
	Double_t avT0res=-9999.,peakT0res=-9999., spreadT0res=-9999.,peakT0resErr=-9999., spreadT0resErr=-9999.;
	Double_t avT0fillRes=-9999.;

	Int_t avMulti=0;
	Float_t fractionEventsWHits=-9999.;
  
	TTree * ttree=new TTree("trendTree","tree of trending variables");
	ttree->Branch("run",&runNumber,"run/I");
	ttree->Branch("avMulti",&avMulti,"avMulti/I");
	ttree->Branch("fractionEventsWHits",&fractionEventsWHits,"fractionEventsWHits/F");
	ttree->Branch("avTime",&avTime,"avTime/D"); //mean time
	ttree->Branch("peakTime",&peakTime,"peakTime/D"); //main peak time after fit
	ttree->Branch("spreadTime",&spreadTime,"spreadTime/D"); //spread of main peak of time after fit
	ttree->Branch("peakTimeErr",&peakTimeErr,"peakTimeErr/D"); //main peak time after fit error
	ttree->Branch("spreadTimeErr",&spreadTimeErr,"spreadTimeErr/D"); //spread of main peak of time after fit error
	ttree->Branch("negTimeRatio",&negTimeRatio,"negTimeRatio/D"); //negative time ratio
  
	ttree->Branch("avRawTime",&avRawTime,"avRawTime/D"); //mean raw time
	ttree->Branch("peakRawTime",&peakRawTime,"peakRawTime/D"); //mean peak of RAW TIME after fit
	ttree->Branch("spreadRawTime",&spreadRawTime,"spreadRawTime/D"); //spread of main peak of raw time after fit
	ttree->Branch("peakRawTimeErr",&peakRawTimeErr,"peakRawTimeErr/D"); //main peak raw  time after fit error
	ttree->Branch("spreadRawTimeErr",&spreadRawTimeErr,"spreadRawTimeErr/D"); //spread of  raw main peak of time after fit error
  
	ttree->Branch("avTot",&avTot,"avTot/D"); //main peak tot
	ttree->Branch("peakTot",&peakTot,"peakTot/D"); // main peak of tot after fit
	ttree->Branch("spreadTot",&spreadTot,"spreadTot/D"); //spread of main peak of tot after fit
	ttree->Branch("peakTotErr",&peakTotErr,"peakTotErr/D"); // main peak of tot after fit
	ttree->Branch("spreadTotErr",&spreadTotErr,"spreadTotErr/D"); //spread of main peak of tot after fit
  
	ttree->Branch("orphansRatio",&orphansRatio,"orphansRatio/D"); //orphans ratio

	ttree->Branch("avL",&avL,"avL/D"); //mean track length
	ttree->Branch("negLratio",&negLratio,"negLratio/D");//ratio of tracks with track length <350 cm

	ttree->Branch("effPt1",&effPt1,"effPt1/D");//matching eff at 1 GeV/c
	ttree->Branch("effPt2",&effPt2,"effPt2/D"); //matching eff at 2 GeV/c
	ttree->Branch("matchEffLinFit1Gev",&matchEffLinFit1Gev,"matchEffLinFit1Gev/D");//matching eff fit param 
	ttree->Branch("matchEffLinFit1GevErr",&matchEffLinFit1GevErr,"matchEffLinFit1GevErr/D");////matching eff fit param error
  
	ttree->Branch("avPiDiffTime",&avPiDiffTime,"avPiDiffTime/D"); //mean t-texp
	ttree->Branch("peakPiDiffTime",&peakPiDiffTime,"peakPiDiffTime/D"); //main peak t-texp after fit
	ttree->Branch("spreadPiDiffTime",&spreadPiDiffTime,"spreadPiDiffTime/D"); //spread of main peak t-texp after fit
	ttree->Branch("peakPiDiffTimeErr",&peakPiDiffTimeErr,"peakPiDiffTimeErr/D"); //main peak t-texp after fit error
	ttree->Branch("spreadPiDiffTimeErr",&spreadPiDiffTimeErr,"spreadPiDiffTimeErr/D"); //spread of main peak of t-texp after fit error

	ttree->Branch("avT0A",&avT0A,"avT0A/D"); //main peak t0A
	ttree->Branch("peakT0A",&peakT0A,"peakT0A/D"); // main peak of t0A after fit
	ttree->Branch("spreadT0A",&spreadT0A,"spreadTot/D"); //spread of main peak of t0A after fit
	ttree->Branch("peakT0AErr",&peakT0AErr,"peakT0AErr/D"); // main peak of t0A after fit
	ttree->Branch("spreadT0AErr",&spreadT0AErr,"spreadT0AErr/D"); //spread of main peak of t0A after fit

	ttree->Branch("avT0C",&avT0C,"avT0C/D"); //main peak t0C
	ttree->Branch("peakT0C",&peakT0C,"peakT0C/D"); // main peak of t0C after fit
	ttree->Branch("spreadT0C",&spreadT0C,"spreadT0C/D"); //spread of main peak of t0C after fit
	ttree->Branch("peakT0CErr",&peakT0CErr,"peakT0CErr/D"); // main peak of t0C after fit
	ttree->Branch("spreadT0CErr",&spreadT0CErr,"spreadT0CErr/D"); //spread of main peak of t0C after fit
 
	ttree->Branch("avT0AC",&avT0AC,"avT0AC/D"); //main peak t0AC
	ttree->Branch("peakT0AC",&peakT0AC,"peakT0AC/D"); // main peak of t0AC after fit
	ttree->Branch("spreadT0AC",&spreadT0AC,"spreadT0AC/D"); //spread of main peak of t0AC after fit
	ttree->Branch("peakT0ACErr",&peakT0ACErr,"peakT0ACErr/D"); // main peak of t0AC after fit
	ttree->Branch("spreadT0ACErr",&spreadT0ACErr,"spreadT0ACErr/D"); //spread of main peak of t0AC after fit
 
	ttree->Branch("avT0res",&avT0res,"avT0res/D"); //main peak t0AC
	ttree->Branch("peakT0res",&peakT0res,"peakT0res/D"); // main peak of t0AC after fit
	ttree->Branch("spreadT0res",&spreadT0res,"spreadT0res/D"); //spread of main peak of t0AC after fit
	ttree->Branch("peakT0resErr",&peakT0resErr,"peakT0resErr/D"); // main peak of t0AC after fit
	ttree->Branch("spreadT0resErr",&spreadT0resErr,"spreadT0resErr/D"); //spread of main peak of t0AC after fit
	ttree->Branch("avT0fillRes",&avT0fillRes,"avT0fillRes/D"); //t0 fill res
  

	FILE * files = fopen(runlist, "r") ; 
	while (fscanf(files,"%d",&runNumber)==1 ){
    
		//reset all variables
		avTime=-9999.; peakTime=-9999.; spreadTime=-9999.; peakTimeErr=-9999.; spreadTimeErr=-9999.;negTimeRatio=-9999.;
		avRawTime=-9999.; peakRawTime=-9999.; spreadRawTime=-9999.; peakRawTimeErr=-9999.; spreadRawTimeErr=-9999.; 
		avTot=-9999.; peakTot=-9999.;spreadTot=-9999.;  peakTotErr=-9999.;spreadTotErr=-9999.;
		orphansRatio=-9999.; avL=-9999.; negLratio=-9999.;
		effPt1=-9999.; effPt2=-9999.; matchEffLinFit1Gev=-9999.;matchEffLinFit1GevErr=-9999.;
    
		avPiDiffTime=-9999.;peakPiDiffTime=-9999.; spreadPiDiffTime=-9999.;peakPiDiffTimeErr=-9999.; spreadPiDiffTimeErr=-9999.;
    
		avT0A=-9999.;peakT0A=-9999.; spreadT0A=-9999.;peakT0AErr=-9999.; spreadT0AErr=-9999.;
		avT0C=-9999.;peakT0C=-9999.; spreadT0C=-9999.;peakT0CErr=-9999.; spreadT0CErr=-9999.;
		avT0AC=-9999.;peakT0AC=-9999.; spreadT0AC=-9999.;peakT0ACErr=-9999.; spreadT0ACErr=-9999.;
		avT0res=-9999.;peakT0res=-9999.; spreadT0res=-9999.;peakT0resErr=-9999.; spreadT0resErr=-9999.;
		avMulti=0; avT0fillRes=0.;
		fractionEventsWHits=-9999.;
    
		//get QAtrain output
		//SetQAtrainOutputName(runNumber, year,period,pass,trainId);
		if (trainId==0){
			sprintf(infile,"alien:///alice/data/%i/%s/000%d/ESDs/%s/QAresults.root",year,period,runNumber,pass);
		} else{
			sprintf(infile,"alien:///alice/data/%i/%s/000%d/ESDs/%s/QA%i/QAresults.root",year,period,runNumber,pass,trainId);
		}
		printf("============== Opening QA file(s) for run %i =======================\n",runNumber);
    
    
		//run post-analysis
		if (RunESDQApostAnalysis(infile,runNumber,kTRUE)==0){
			filesCounter++;
      
			//get post-analysis output
			sprintf(postFileName,"postQA_%i.root",runNumber);
			TFile *postfile=TFile::Open(postFileName);
			if (!postfile) {
				printf("Post-analysis output not found - cannot perform trending analysis. Exiting.");
				return 2;
			} else {
				printf("==============  Retrieving post-analysis for run %i ===============\n",runNumber);
			}
			//postfile->ls();
			//save quantities for trending
      
			//--------------------------------- Multiplicity ----------------------------------//
			TH1F*hMulti=(TH1F*)postfile->Get("hTOFmatchedPerEvt");
			if ((hMulti)&&(hMulti->GetEntries()>0)) {
				avMulti=hMulti->GetMean();
			}
      
			//--------------------------------- T0F signal ----------------------------------//
      
			TH1F*hTime=(TH1F*)postfile->Get("hTOFmatchedESDtime");
			if ((hTime)&&(hTime->GetEntries()>0)) {
				avTime=hTime->GetMean();
				hTime->Fit("landau","","",0.,50.);
				peakTime=(hTime->GetFunction("landau"))->GetParameter(1);
				spreadTime=(hTime->GetFunction("landau"))->GetParameter(2);
				peakTimeErr=(hTime->GetFunction("landau"))->GetParError(1);
				spreadTimeErr=(hTime->GetFunction("landau"))->GetParError(2);
	
				negTimeRatio=((Double_t)hTime->Integral(1,3)*100.)/((Double_t)hTime->Integral());
				printf("Main peak time (landau): mean = %f +- %f\n",peakTime,peakTimeErr );
				printf("Main peak time (landau): spread = %f +- %f\n",spreadTime,spreadTimeErr );
				printf("Ratio of tracks with time<12.5 ns / total = %f\n",negTimeRatio );
	
				//add integral of main peak over total
			}
			printf("---------------------------------------------------------------- \n");
			TH1F * hRawTime = (TH1F*)postfile->Get("hTOFmatchedESDrawTime");
			if ((hRawTime)&&(hRawTime->GetEntries()>0)){
				avRawTime=hRawTime->GetMean();
				hRawTime->Fit("landau","","",200.,250.);
				peakRawTime=(hRawTime->GetFunction("landau"))->GetParameter(1);
				spreadRawTime=(hRawTime->GetFunction("landau"))->GetParameter(2);
				peakRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(1);
				spreadRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(2);
				printf("Main peak raw time (landau): mean = %f +- %f\n",peakTime,peakTimeErr );
				printf("Main peak raw time (landau): spread = %f +- %f\n",spreadRawTime,spreadRawTimeErr );
			}
      
			printf("---------------------------------------------------------------- \n");
      
			TH1F * hTot = (TH1F*)postfile->Get("hTOFmatchedESDToT");
			if ((hTot)&&(hTot->GetEntries()>0)){
				avTot=hTot->GetMean();
				hTot->Fit("gaus","","",0.,50.);
				peakTot=(hTot->GetFunction("gaus"))->GetParameter(1);
				spreadTot=(hTot->GetFunction("gaus"))->GetParameter(2);
				peakTotErr=(hTot->GetFunction("gaus"))->GetParError(1);
				spreadTotErr=(hTot->GetFunction("gaus"))->GetParError(2);
				printf("Main peak ToT (gaus): mean = %f +- %f\n",peakTot,peakTotErr );
				printf("Main peak ToT (gaus): spread = %f +- %f\n",spreadTot,spreadTotErr );	
			}      
			printf("---------------------------------------------------------------- \n");
			TH1F * hOrphansRatio=(TH1F*)postfile->Get("hOrphansRatio");
			if (hOrphansRatio)
				orphansRatio=hOrphansRatio->GetMean();
      
			TH1F * hL=(TH1F*)postfile->Get("hTOFmatchedESDtrkLength");
			if (hL)
				avL=hL->GetMean();
      
			TH1F *hLnegRatio =(TH1F*)postfile->Get("hLnegRatio");
			if (hLnegRatio)
				negLratio=hLnegRatio->GetMean();
      
			//--------------------------------- matching eff ----------------------------------//
     
			//Double_t linFitEff1Param[3]={0.,0.,0.};
			TH1F *hMatchingVsPt =(TH1F*)postfile->Get("hTOFmatchedESDPt");
			if (hMatchingVsPt) {
				if (hMatchingVsPt->GetEntries()>0){
					hMatchingVsPt->Fit("pol0","","",1.0,5.);
					hMatchingVsPt->Draw();
					if (hMatchingVsPt->GetFunction("pol0")){
						matchEffLinFit1Gev=(hMatchingVsPt->GetFunction("pol0"))->GetParameter(0);
						matchEffLinFit1GevErr=(hMatchingVsPt->GetFunction("pol0"))->GetParError(0);	
						printf("Matching efficiency fit param is %f +- %f\n",matchEffLinFit1Gev,matchEffLinFit1GevErr );
					}
				} else {
					printf("WARNING: matching efficiency plot has 0 entries. Skipped!\n");
				}
			} else {
				printf("WARNING: cannot retrieve matching efficiency plot\n");
			}

			//--------------------------------- t-texp ----------------------------------//
     
			TH1F*hPiDiffTime=(TH1F*)postfile->Get("hTOFmatchedExpTimePi");
			if ((hPiDiffTime)&&(hPiDiffTime->GetEntries()>0)) {
				avPiDiffTime=hPiDiffTime->GetMean();
				hPiDiffTime->Fit("gaus","","",-1000.,500.);
				if (hPiDiffTime->GetFunction("gaus")){
					peakPiDiffTime=(hPiDiffTime->GetFunction("gaus"))->GetParameter(1);
					spreadPiDiffTime=(hPiDiffTime->GetFunction("gaus"))->GetParameter(2);
					peakPiDiffTimeErr=(hPiDiffTime->GetFunction("gaus"))->GetParError(1);
					spreadPiDiffTimeErr=(hPiDiffTime->GetFunction("gaus"))->GetParError(2);
					printf("Main peak t-t_exp (gaus): mean = %f +- %f\n",peakPiDiffTime,peakPiDiffTimeErr );
					printf("Main peak t-t_exp (gaus): spread = %f +- %f\n",spreadPiDiffTime,spreadPiDiffTimeErr );
				}
			}
			//--------------------------------- T0 detector ----------------------------------//
      
			TH1F*hT0A=(TH1F*)postfile->Get("hEventT0DetA");
			if ((hT0A)&&(hT0A->GetEntries()>0)) {
				avhT0A=hT0A->GetMean();
				hT0A->Fit("gaus");
				peakT0A=(hT0A->GetFunction("gaus"))->GetParameter(1);
				spreadT0A=(hT0A->GetFunction("gaus"))->GetParameter(2);
				peakT0AErr=(hT0A->GetFunction("gaus"))->GetParError(1);
				spreadT0AErr=(hT0A->GetFunction("gaus"))->GetParError(2);	
				printf("Main peak T0A(gaus): mean = %f +- %f\n",peakT0A,peakT0AErr );
				printf("Main peak T0A (gaus): spread = %f +- %f\n",spreadT0A,spreadT0AErr );	 
				//add integral of main peak over total
			}

			TH1F*hT0C=(TH1F*)postfile->Get("hEventT0DetC");
			if ((hT0C)&&(hT0C->GetEntries()>0)) {
				avhT0C=hT0C->GetMean();
				hT0C->Fit("gaus");
				peakT0C=(hT0C->GetFunction("gaus"))->GetParameter(1);
				spreadT0C=(hT0C->GetFunction("gaus"))->GetParameter(2);
				peakT0CErr=(hT0C->GetFunction("gaus"))->GetParError(1);
				spreadT0CErr=(hT0C->GetFunction("gaus"))->GetParError(2);	
				printf("Main peak T0C(gaus): mean = %f +- %f\n",peakT0C,peakT0CErr );
				printf("Main peak T0C (gaus): spread = %f +- %f\n",spreadT0C,spreadT0CErr );	 
				//add integral of main peak over total
			}

			TH1F*hT0AC=(TH1F*)postfile->Get("hEventT0DetAND");
			if ((hT0AC)&&(hT0AC->GetEntries()>0)) {
				avhT0AC=hT0AC->GetMean();
				hT0AC->Fit("gaus");
				peakT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(1);
				spreadT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(2);
				peakT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(1);
				spreadT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(2);	
				printf("Main peak T0AC(gaus): mean = %f +- %f\n",peakT0AC,peakT0ACErr );
				printf("Main peak T0AC (gaus): spread = %f +- %f\n",spreadT0AC,spreadT0ACErr );	 
			}
      
			TH1F*hT0res=(TH1F*)postfile->Get("hT0DetRes");
			if ((hT0res)&&(hT0res->GetEntries()>0)) {
				avhT0res=hT0res->GetMean();
				hT0res->Fit("gaus");
				peakT0res=(hT0res->GetFunction("gaus"))->GetParameter(1);
				spreadT0res=(hT0res->GetFunction("gaus"))->GetParameter(2);
				peakT0resErr=(hT0res->GetFunction("gaus"))->GetParError(1);
				spreadT0resErr=(hT0res->GetFunction("gaus"))->GetParError(2);	
				printf("Main peak T0res(gaus): mean = %f +- %f\n",peakT0res,peakT0resErr );
				printf("Main peak T0res (gaus): spread = %f +- %f\n",spreadT0res,spreadT0resErr );	 
				//add integral of main peak over total
			}

			TH1F*hT0fillRes=(TH1F*)postfile->Get("hT0fillRes");
			if ((hT0fillRes)&&(hT0fillRes->GetEntries()>0)) {
				avT0fillRes=hT0fillRes->GetMean();
			}
      
			ttree->Fill();
			printf("==============  Saving trending quantities for run %i ===============\n",runNumber);
			if (postfile->IsOpen()) {
				printf("Trying to close\n");
				postfile->Close();
			}  
		} else {
			printf("==============   QA for run %i not available - skipping ===============\n",runNumber);
		}
    
	} 
	printf("Number of files processed = %i\n",filesCounter);
  
	trendFile->cd();
	ttree->Write();
	trendFile->Close();

	return  MakeTrendingFromTreeWithErrors(trendFileName); 
}

//______________________________________________________________________________

Int_t MakeTrendingFromTreeWithErrors(char* trendFileName=NULL){

	if (!trendFileName) 
		return 3;
  
	TFile *fin=TFile::Open(trendFileName);
	if (!fin) 
		return 4;
	Int_t runNumber;
	Double_t avTime=0., peakTime=0., spreadTime=0., peakTimeErr=0., spreadTimeErr=0.,negTimeRatio=0.,
		avRawTime=0., peakRawTime=0., spreadRawTime=0., peakRawTimeErr=0., spreadRawTimeErr=0., 
		avTot=0., peakTot=0.,spreadTot=0.,  peakTotErr=0.,spreadTotErr=0.,
		orphansRatio=0., avL=0., negLratio=0.,
		effPt1=0., effPt2=0., matchEffLinFit1Gev=0.,matchEffLinFit1GevErr=0.;
	Double_t avDiffTime=0.,peakDiffTime=0., spreadDiffTime=0.,peakDiffTimeErr=0., spreadDiffTimeErr=0.,avT0fillRes=0.;
   
	Double_t avT0A=0.,peakT0A=0., spreadT0A=0.,peakT0AErr=0., spreadT0AErr=0.;
	Double_t avT0C=0.,peakT0C=0., spreadT0C=0.,peakT0CErr=0., spreadT0CErr=0.;
	Double_t avT0AC=0.,peakT0AC=0., spreadT0AC=0.,peakT0ACErr=0., spreadT0ACErr=0.;
	Double_t avT0res=0.,peakT0res=0., spreadT0res=0.,peakT0resErr=0., spreadT0resErr=0.;
	Int_t avMulti=0;
	Float_t fractionEventsWHits=0;
   
	TTree * ttree=(TTree*)fin->Get("trendTree");
	ttree->SetBranchAddress("run",&runNumber);
   
	ttree->SetBranchAddress("avTime",&avTime); //mean time
	ttree->SetBranchAddress("peakTime",&peakTime); //main peak time after fit
	ttree->SetBranchAddress("spreadTime",&spreadTime); //spread of main peak of time after fit
	ttree->SetBranchAddress("peakTimeErr",&peakTimeErr); //main peak time after fit error
	ttree->SetBranchAddress("spreadTimeErr",&spreadTimeErr); //spread of main peak of time after fit error
   
	ttree->SetBranchAddress("negTimeRatio",&negTimeRatio); //negative time ratio
   
	ttree->SetBranchAddress("avRawTime",&avRawTime); //mean raw time
	ttree->SetBranchAddress("peakRawTime",&peakRawTime); //mean peak of raw time after fit
	ttree->SetBranchAddress("spreadRawTime",&spreadRawTime); //spread of main peak of raw time after fit
	ttree->SetBranchAddress("peakRawTimeErr",&peakRawTimeErr); //main peak raw  time after fit error
	ttree->SetBranchAddress("spreadRawTimeErr",&spreadRawTimeErr); //spread of  raw main peak of time after fit error
  
	ttree->SetBranchAddress("avTot",&avTot); //main peak tot
	ttree->SetBranchAddress("peakTot",&peakTot); // main peak of tot after fit
	ttree->SetBranchAddress("spreadTot",&spreadTot); //spread of main peak of tot after fit
	ttree->SetBranchAddress("peakTotErr",&peakTotErr); // main peak of tot after fit
	ttree->SetBranchAddress("spreadTotErr",&spreadTotErr); //spread of main peak of tot after fit
   
	ttree->SetBranchAddress("orphansRatio",&orphansRatio); //orphans ratio
   
	ttree->SetBranchAddress("avL",&avL); //mean track length
	ttree->SetBranchAddress("negLratio",&negLratio);//ratio of tracks with track length <350 cm
   
	ttree->SetBranchAddress("effPt1",&effPt1);//matching eff at 1 GeV/c
	ttree->SetBranchAddress("effPt2",&effPt2); //matching eff at 2 GeV/c
	ttree->SetBranchAddress("matchEffLinFit1Gev",&matchEffLinFit1Gev);//matching eff fit param 
	ttree->SetBranchAddress("matchEffLinFit1GevErr",&matchEffLinFit1GevErr);////matching eff fit param error
   
	ttree->SetBranchAddress("avPiDiffTime",&avDiffTime); //mean t-texp
	ttree->SetBranchAddress("peakPiDiffTime",&peakDiffTime); //main peak t-texp after fit
	ttree->SetBranchAddress("spreadPiDiffTime",&spreadDiffTime); //spread of main peak t-texp after fit
	ttree->SetBranchAddress("peakPiDiffTimeErr",&peakDiffTimeErr); //main peak t-texp after fit error
	ttree->SetBranchAddress("spreadPiDiffTimeErr",&spreadDiffTimeErr); //spread of main peak of t-texp after fit error
   
	ttree->SetBranchAddress("avT0A",&avT0A); //main peak t0A
	ttree->SetBranchAddress("peakT0A",&peakT0A); // main peak of t0A after fit
	ttree->SetBranchAddress("spreadT0A",&spreadT0A); //spread of main peak of t0A after fit
	ttree->SetBranchAddress("peakT0AErr",&peakT0AErr); // main peak of t0A after fit
	ttree->SetBranchAddress("spreadT0AErr",&spreadT0AErr); //spread of main peak of t0A after fit
   
	ttree->SetBranchAddress("avT0C",&avT0C); //main peak t0C
	ttree->SetBranchAddress("peakT0C",&peakT0C); // main peak of t0C after fit
	ttree->SetBranchAddress("spreadT0C",&spreadT0C); //spread of main peak of t0C after fit
	ttree->SetBranchAddress("peakT0CErr",&peakT0CErr); // main peak of t0C after fit
	ttree->SetBranchAddress("spreadT0CErr",&spreadT0CErr); //spread of main peak of t0C after fit
   
	ttree->SetBranchAddress("avT0AC",&avT0AC); //main peak t0AC
	ttree->SetBranchAddress("peakT0AC",&peakT0AC); // main peak of t0AC after fit
	ttree->SetBranchAddress("spreadT0AC",&spreadT0AC); //spread of main peak of t0AC after fit
	ttree->SetBranchAddress("peakT0ACErr",&peakT0ACErr); // main peak of t0AC after fit
	ttree->SetBranchAddress("spreadT0ACErr",&spreadT0ACErr); //spread of main peak of t0AC after fit
   
	ttree->SetBranchAddress("avT0res",&avT0res); //main peak t0AC
	ttree->SetBranchAddress("peakT0res",&peakT0res); // main peak of t0AC after fit
	ttree->SetBranchAddress("spreadT0res",&spreadT0res); //spread of main peak of t0AC after fit
	ttree->SetBranchAddress("peakT0resErr",&peakT0resErr); // main peak of t0AC after fit
	ttree->SetBranchAddress("spreadT0resErr",&spreadT0resErr); //spread of main peak of t0AC after fit
	ttree->SetBranchAddress("avT0fillRes",&avT0fillRes); //t0 fill res

	Int_t nRuns=ttree->GetEntries();
	TList lista;
   
	TH1F * hAvDiffTimeVsRun=new TH1F("hAvDiffTimeVsRun","hAvDiffTimeVsRun;run;<t^{TOF}-t_{exp,#pi}> (ps)",nRuns,0., nRuns);//, 600, 0. , 600.);
	hAvDiffTimeVsRun->SetDrawOption("E1");
	hAvDiffTimeVsRun->SetMarkerStyle(20);
	hAvDiffTimeVsRun->SetMarkerColor(kBlue);
	//   hAvTimeVsRun->GetYaxis()->SetRangeUser(0.0, 50.0);

	TH1F * hPeakDiffTimeVsRun=new TH1F("hPeakDiffTimeVsRun","hPeakDiffTimeVsRun (gaussian fit) ;run; <t^{TOF}-t_{exp,#pi}> (ps)",nRuns,0., nRuns);//,600, 0. , 600. );
	hPeakDiffTimeVsRun->SetDrawOption("E1");
	hPeakDiffTimeVsRun->SetMarkerStyle(20);
	hPeakDiffTimeVsRun->SetMarkerColor(kBlue);
   
	TH1F * hSpreadDiffTimeVsRun=new TH1F("hSpreadDiffTimeVsRun","hSpreadDiffTimeVsRun (gaussian fit);run; #sigma(t^{TOF}-t_{exp,#pi}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
	hSpreadDiffTimeVsRun->SetDrawOption("E1");
	hSpreadDiffTimeVsRun->SetMarkerStyle(20);
	hSpreadDiffTimeVsRun->SetMarkerColor(kBlue);

	TH1F * hAvTimeVsRun=new TH1F("hAvTimeVsRun","hAvTimeVsRun;run;<t^{TOF}> (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
	hAvTimeVsRun->SetDrawOption("E1");
	hAvTimeVsRun->SetMarkerStyle(20);
	hAvTimeVsRun->SetMarkerColor(kBlue);
	//   hAvTimeVsRun->GetYaxis()->SetRangeUser(0.0, 50.0);

	TH1F * hPeakTimeVsRun=new TH1F("hPeakTimeVsRun","hPeakTimeVsRun (gaussian fit);run;t_{peak}^{TOF} (ns)",nRuns,0., nRuns);//,600, 0. , 600. );
	hPeakTimeVsRun->SetDrawOption("E1");
	hPeakTimeVsRun->SetMarkerStyle(20);
	hPeakTimeVsRun->SetMarkerColor(kBlue);
   
	TH1F * hSpreadTimeVsRun=new TH1F("hSpreadTimeVsRun","hSpreadTimeVsRun (gaussian fit);run; #sigma(t^{TOF}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
	hSpreadTimeVsRun->SetDrawOption("E1");
	hSpreadTimeVsRun->SetMarkerStyle(20);
	hSpreadTimeVsRun->SetMarkerColor(kBlue);
  
  
	TH1F * hAvRawTimeVsRun=new TH1F("hAvRawTimeVsRun","hAvRawTimeVsRun;run;<t_{raw}^{TOF}> (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
	hAvRawTimeVsRun->SetDrawOption("E1");
	hAvRawTimeVsRun->SetMarkerStyle(21);
	hAvRawTimeVsRun->SetMarkerColor(kGreen);

	TH1F * hPeakRawTimeVsRun=new TH1F("hPeakRawTimeVsRun","hPeakRawTimeVsRun (gaussian fit);run;t_{peak,raw}^{TOF} (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
	hPeakRawTimeVsRun->SetDrawOption("E1");
	hPeakRawTimeVsRun->SetMarkerStyle(21);
	hPeakRawTimeVsRun->SetMarkerColor(kGreen);

	TH1F * hSpreadRawTimeVsRun=new TH1F("hSpreadRawTimeVsRun","hSpreadRawTimeVsRun (gaussian fit);run;#sigma(t_{raw}^{TOF}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
	hSpreadRawTimeVsRun->SetDrawOption("E1");
	hSpreadRawTimeVsRun->SetMarkerStyle(21);
	hSpreadRawTimeVsRun->SetMarkerColor(kGreen);
   
	TH1F * hAvTotVsRun=new TH1F("hAvTotVsRun","hAvTotVsRun;run;<ToT> (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
	hAvTotVsRun->SetDrawOption("E1");
	hAvTotVsRun->SetMarkerStyle(22);
   
	TH1F * hPeakTotVsRun=new TH1F("hPeakTotVsRun","hPeakTotVsRun (gaussian fit);run;ToT_{peak} (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
	hPeakTotVsRun->SetDrawOption("E1");
	hPeakTotVsRun->SetMarkerStyle(22);
   
	TH1F * hSpreadTotVsRun=new TH1F("hSpreadTotVsRun","hSpreadTotVsRun (gaussian fit);#sigma(ToT) (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
	hSpreadTotVsRun->SetDrawOption("E1");
	hSpreadTotVsRun->SetMarkerStyle(22);
   
	TH1F * hNegTimeRatioVsRun=new TH1F("hNegTimeRatioVsRun","hNegTimeRatioVsRun;run;ratio of tracks with t^{TOF}<12.5 ns (%)",nRuns, 0., nRuns);//, 100, 0. , 100.);
	hNegTimeRatioVsRun->SetDrawOption("E");

	TH1F * hOrphansRatioVsRun=new TH1F("hOrphansRatioVsRun","hOrphansRatioVsRun; run; ratio of orphans (%)",nRuns, 0., nRuns);//, 1000, 0. , 100.);
	hOrphansRatioVsRun->SetDrawOption("E");

	TH1F * hMeanLVsRun=new TH1F("hMeanLVsRun","hMeanLVsRun;run; <L> (cm)",nRuns, 0., nRuns);//, 350, 350. , 700.);
	hMeanLVsRun->SetDrawOption("E");
	TH1F * hNegLRatioVsRun=new TH1F("hNegLRatioVsRun","hNegLRatioVsRun;run; ratio of tracks with L<350 cm (%)",nRuns, 0., nRuns);//, 1000, 0. , 100.);
	hNegLRatioVsRun->SetDrawOption("E");
	TH1F * hMatchEffVsRun=new TH1F("hMatchEffVsRun","hMatchEffVsRun;run;matching efficiency (pT>1.5 GeV/c)",nRuns, 0., nRuns);//, 100, 0. , 1.);
	hMatchEffVsRun->SetDrawOption("E");
	TH1F * hMatchEffVsRun1=new TH1F("hMatchEffVsRun1","hMatchEffVsRun;run;matching efficiency (pT>1.0 GeV/c)",nRuns, 0., nRuns);
	hMatchEffVsRun1->SetDrawOption("E");
	TH1F * hPeakT0AVsRun=new TH1F("hPeakT0AVsRun","hPeakT0AVsRun (gaussian fit);run;t0A (ps)",nRuns,0., nRuns);
	TH1F * hPeakT0CVsRun=new TH1F("hPeakT0CVsRun","hPeakT0CVsRun (gaussian fit);run;t0AC (ps)",nRuns,0., nRuns);
	TH1F * hPeakT0ACVsRun=new TH1F("hPeakT0ACVsRun","hPeakT0ACVsRun (gaussian fit);run;t0AC (ps)",nRuns,0., nRuns);
	TH1F * hT0fillResVsRun=new TH1F("hT0fillResVsRun","hT0fillResVsRun;run;t0_spread (ps)",nRuns,0., nRuns);
  
	lista.Add(hAvDiffTimeVsRun);
	lista.Add(hPeakDiffTimeVsRun);
	lista.Add(hSpreadDiffTimeVsRun);
	lista.Add(hAvTimeVsRun);
	lista.Add(hPeakTimeVsRun);
	lista.Add(hSpreadTimeVsRun);
	lista.Add(  hAvRawTimeVsRun);
	lista.Add(  hPeakRawTimeVsRun);
	lista.Add(  hSpreadRawTimeVsRun); 
	lista.Add(  hAvTotVsRun);
	lista.Add(  hPeakTotVsRun);
	lista.Add(  hSpreadTotVsRun);
	lista.Add(  hNegTimeRatioVsRun);
	lista.Add(  hOrphansRatioVsRun);
	lista.Add( hMeanLVsRun);
	lista.Add(  hNegLRatioVsRun);
	lista.Add(  hMatchEffVsRun);
	lista.Add(hPeakT0AVsRun);
	lista.Add(hPeakT0CVsRun);
	lista.Add(hPeakT0ACVsRun);
	lista.Add(hT0fillResVsRun);
	char runlabel[6];
   
	for (Int_t irun=0;irun<nRuns;irun++){
		ttree->GetEntry(irun);
    
		sprintf(runlabel,"%i",runNumber);
    
		hAvDiffTimeVsRun->SetBinContent(irun+1, avDiffTime);
		hAvDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hPeakDiffTimeVsRun->SetBinContent(irun+1,peakDiffTime);
		hPeakDiffTimeVsRun->SetBinError(irun+1,peakDiffTimeErr);
		hPeakDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hSpreadDiffTimeVsRun->SetBinContent(irun+1,spreadDiffTime);
		hSpreadDiffTimeVsRun->SetBinError(irun+1,spreadDiffTimeErr);
		hSpreadDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hAvTimeVsRun->SetBinContent(irun+1, avTime);
		hAvTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hPeakTimeVsRun->SetBinContent(irun+1,peakTime);
		hPeakTimeVsRun->SetBinError(irun+1,peakTimeErr);
		hPeakTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hSpreadTimeVsRun->SetBinContent(irun+1,spreadTime);
		hSpreadTimeVsRun->SetBinError(irun+1,spreadTimeErr);
		hSpreadTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hAvRawTimeVsRun->SetBinContent(irun+1, avRawTime);
		hAvRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hPeakRawTimeVsRun->SetBinContent(irun+1,peakRawTime);
		hPeakRawTimeVsRun->SetBinError(irun+1,peakRawTimeErr);
		hPeakRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hSpreadRawTimeVsRun->SetBinContent(irun+1,spreadRawTime);
		hSpreadRawTimeVsRun->SetBinError(irun+1,spreadRawTimeErr);
		hSpreadRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hAvTotVsRun->SetBinContent(irun,avTot);
		hAvTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hPeakTotVsRun->SetBinContent(irun,peakTot);
		hPeakTotVsRun->SetBinError(irun+1,peakTotErr);
		hPeakTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hSpreadTotVsRun->SetBinContent(irun,spreadTot);
		hSpreadTotVsRun->SetBinError(irun+1,spreadTotErr);
		hSpreadTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hNegTimeRatioVsRun->SetBinContent(irun,negTimeRatio);
		hNegTimeRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hOrphansRatioVsRun->SetBinContent(irun,orphansRatio);
		hOrphansRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hMeanLVsRun->SetBinContent(irun,avL);
		hMeanLVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hNegLRatioVsRun->SetBinContent(irun,negLratio);
		hNegLRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hMatchEffVsRun->SetBinContent(irun,matchEffLinFit1Gev);
		hMatchEffVsRun->SetBinError(irun+1,matchEffLinFit1GevErr);
		hMatchEffVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hPeakT0AVsRun->SetBinContent(irun+1,peakT0A);
		hPeakT0AVsRun->SetBinError(irun+1,spreadT0A);
		hPeakT0AVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hPeakT0CVsRun->SetBinContent(irun+1,peakT0C);
		hPeakT0CVsRun->SetBinError(irun+1,spreadT0C);
		hPeakT0CVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
		hPeakT0ACVsRun->SetBinContent(irun+1,peakT0AC);
		hPeakT0ACVsRun->SetBinError(irun+1,spreadT0AC);
		hPeakT0ACVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

		hT0fillResVsRun->SetBinContent(irun+1,avT0fillRes);
		hT0fillResVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

	}
  
	char  outfilename[200];
	sprintf(outfilename, "trend_%s",trendFileName);
	TFile * fout=new TFile(outfilename,"recreate");
	fout->cd();
	lista.Write();
	fout->Close();

	TString plotDir(Form("PlotsTrending"));
	gSystem->Exec(Form("mkdir %s",plotDir.Data()));

	TCanvas* cAvDiffTimeVsRun = new TCanvas("cAvDiffTimeVsRun","cAvDiffTimeVsRun",50,50,750,550);
	gPad->SetGridx();
	gPad->SetGridy();
	hAvDiffTimeVsRun->Draw();
	cAvDiffTimeVsRun->Print(Form("%s/cAvDiffTimeVsRun.png",plotDir.Data()));

	TCanvas* cPeakDiffTimeVsRun = new TCanvas("cPeakDiffTimeVsRun","cPeakDiffTimeVsRun", 50,50,750,550);
	hPeakDiffTimeVsRun->Draw();
	cPeakDiffTimeVsRun->Print(Form("%s/cPeakDiffTimeVsRun",plotDir.Data()));

	TCanvas* cSpreadDiffTimeVsRun = new TCanvas("cSpreadDiffTimeVsRun","cSpreadDiffTimeVsRun", 50,50,750,550);
	hSpreadDiffTimeVsRun->Draw();
	cSpreadDiffTimeVsRun->Print(Form("%s/cSpreadDiffTimeVsRun",plotDir.Data()));

	TCanvas* cAvTimeVsRun = new TCanvas("cAvTimeVsRun","cAvTimeVsRun", 50,50,750,550);
	hAvTimeVsRun->Draw();
	cAvTimeVsRun->Print(Form("%s/cAvTimeVsRun",plotDir.Data()));

	TCanvas* cPeakTimeVsRun = new TCanvas("cPeakTimeVsRun","cPeakTimeVsRun", 50,50,750,550);
	hPeakTimeVsRun->Draw();
	cPeakTimeVsRun->Print(Form("%s/cPeakTimeVsRun",plotDir.Data()));

	TCanvas* cSpreadTimeVsRun = new TCanvas("cSpreadTimeVsRun","cSpreadTimeVsRun", 50,50,750,550);
	hSpreadTimeVsRun->Draw();
	cSpreadTimeVsRun->Print(Form("%s/cSpreadTimeVsRun",plotDir.Data()));

	TCanvas* cAvRawTimeVsRun = new TCanvas("cAvRawTimeVsRun","cAvRawTimeVsRun", 50,50,750,550);
	hAvRawTimeVsRun->Draw();
	cAvRawTimeVsRun->Print(Form("%s/cAvRawTimeVsRun",plotDir.Data()));

	TCanvas* cPeakRawTimeVsRun = new TCanvas("cPeakRawTimeVsRun","cPeakRawTimeVsRun", 50,50,750,550);
	hPeakRawTimeVsRun->Draw();
	cPeakRawTimeVsRun->Print(Form("%s/cPeakRawTimeVsRun",plotDir.Data()));

	TCanvas* cSpreadRawTimeVsRun = new TCanvas("cSpreadRawTimeVsRun","cSpreadRawTimeVsRun", 50,50,750,550);
	hSpreadRawTimeVsRun->Draw();
	cSpreadRawTimeVsRun->Print(Form("%s/cSpreadRawTimeVsRun",plotDir.Data()));

	TCanvas* cAvTotVsRun = new TCanvas("cAvTotVsRun","cAvTotVsRun", 50,50,750,550);
	hAvTotVsRun->Draw();
	cAvTotVsRun->Print(Form("%s/cAvTotVsRun",plotDir.Data()));

	TCanvas* cPeakTotVsRun = new TCanvas("cPeakTotVsRun","cPeakTotVsRun", 50,50,750,550);
	hPeakTotVsRun->Draw();
	cPeakTotVsRun->Print(Form("%s/cPeakTotVsRun",plotDir.Data()));

	TCanvas* cSpreadTotVsRun = new TCanvas("cSpreadTotVsRun","cSpreadTotVsRun", 50,50,750,550);
	hSpreadTotVsRun->Draw();
	cSpreadTotVsRun->Print(Form("%s/cSpreadTotVsRun",plotDir.Data()));

	TCanvas* cNegTimeRatioVsRun = new TCanvas("cNegTimeRatioVsRun","cNegTimeRatioVsRun", 50,50,750,550);
	hNegTimeRatioVsRun->Draw();
	cNegTimeRatioVsRun->Print(Form("%s/cNegTimeRatioVsRun",plotDir.Data()));

	TCanvas* cOrphansRatioVsRun = new TCanvas("cOrphansRatioVsRun","cOrphansRatioVsRun", 50,50,750,550);
	hOrphansRatioVsRun->Draw();
	cOrphansRatioVsRun->Print(Form("%s/cOrphansRatioVsRun",plotDir.Data()));

	TCanvas* cMeanLVsRun = new TCanvas("cMeanLVsRun","cMeanLVsRun", 50,50,750,550);
	hMeanLVsRun->Draw();
	cMeanLVsRun->Print(Form("%s/cMeanLVsRun",plotDir.Data()));

	TCanvas* cNegLRatioVsRun = new TCanvas("cNegLRatioVsRun","cNegLRatioVsRun", 50,50,750,550);
	hNegLRatioVsRun->Draw();
	cNegLRatioVsRun->Print(Form("%s/cNegLRatioVsRun",plotDir.Data()));

	TCanvas* cMatchEffVsRun = new TCanvas("cMatchEffVsRun","cMatchEffVsRun", 50,50,750,550);
	hMatchEffVsRun->Draw();
	cMatchEffVsRun->Print(Form("%s/cMatchEffVsRun",plotDir.Data()));

	TCanvas* cPeakT0AVsRun = new TCanvas("cPeakT0AVsRun","cPeakT0AVsRun", 50,50,750,550);
	hPeakT0AVsRun->Draw();
	cPeakT0AVsRun->Print(Form("%s/cPeakT0AVsRun.png",plotDir.Data()));

	TCanvas* cPeakT0CVsRun = new TCanvas("cPeakT0CVsRun","cPeakT0CVsRun", 50,50,750,550);
	hPeakT0CVsRun->Draw();
	cPeakT0CVsRun->Print(Form("%s/cPeakT0CVsRun.png",plotDir.Data()));

	TCanvas* cPeakT0ACVsRun = new TCanvas("cPeakT0ACVsRun","cPeakT0ACVsRun", 50,50,750,550);
	hPeakT0ACVsRun->Draw();
	cPeakT0ACVsRun->Print(Form("%s/cPeakT0ACVsRun.png",plotDir.Data()));

	TCanvas* cT0fillResVsRun = new TCanvas("cT0fillResVsRun","cT0fillResVsRun", 50,50,750,550);
	hT0fillResVsRun->Draw();
	cT0fillResVsRun->Print(Form("%s/cT0fillResVsRun.png",plotDir.Data()));


	return 0;
}

//------------------------------------------------------------------------------------

Int_t RunESDQApostAnalysis(char *qafilename=NULL, Int_t runNumber=-1, Bool_t IsOnGrid=kFALSE, Bool_t canvasE=kFALSE) {
  
	Bool_t debug=kFALSE;
  
	/*access qa PWG1 output files - saved locally or on grid as specified by the second argument */
  
	char defaultQAoutput[30]="QAresults.root";
	if (IsOnGrid) TGrid::Connect("alien://");
	TFile * fin= TFile::Open(qafilename,"r");
	printf("Opening file %s\n",qafilename);
	if (!fin) {
		printf("ERROR: QA output not found. Exiting with status -1\n");
		return -1;
	} else {
		printf("INFO: QA output file %s open. \n",fin->GetName());
	}
  
 
	//access histograms lists
	char tofQAdirName[15]="TOF_Performance";
	char genListName[15]="cGeneralTOFqa";
	char t0ListName[15]="cTimeZeroTOFqa";
	char pidListName[15]="cPIDTOFqa"; 
	char posListName[15]="cPositiveTOFqa";
	char negListName[15]="cNegativeTOFqa";
  
	TDirectoryFile * tofQAdir=(TDirectoryFile*)fin->Get(tofQAdirName);
	if(debug){
		printf("------------------------------------------------------------------\n");
		tofQAdir->ls();
		printf("------------------------------------------------------------------\n");
	}
  
	TList * generalList=(TList*)tofQAdir->Get(genListName);
	TList  *timeZeroList=(TList*)tofQAdir->Get(t0ListName);
	TList  *pidList=(TList*)tofQAdir->Get(pidListName);
	TList  *posList=(TList*)tofQAdir->Get(posListName);
	TList  *negList=(TList*)tofQAdir->Get(negListName);

	if (!generalList) printf("WARNING: general QA histograms absent or not accessible\n");
	if (!timeZeroList) printf("WARNING: timeZero QA histograms absent or not accessible\n");
	if (!pidList) printf("WARNING: PID QA histograms absent or not accessible\n");
	if (!posList) printf("WARNING: general QA histograms for positive tracks absent or not accessible\n");
	if (!negList) printf("WARNING: general QA histograms for negative tracks absent or not accessible\n");

	if ( (!generalList) && (!timeZeroList) && (!pidList) ){
		printf("ERROR: no QA available \n");
		return 1;
	}
  
	if (debug){
		generalList->ls();
		printf("------------------------------------------------------------------\n");
		timeZeroList->ls();
		printf("------------------------------------------------------------------\n");
		pidList->ls();
		printf("------------------------------------------------------------------\n");
	}

	/*set graphic style*/
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTitleBorderSize(0)  ;
	gStyle->SetTitleFont(42);
	gStyle->SetTextFont(42);
	//gStyle->SetPalette(1);
	gStyle->SetStatColor(kWhite); 
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptStat(0);
  
	//DEFINE OUTPUT FILE 
	char outfilename[60];
	sprintf(outfilename,"postQA_%i.root",runNumber);
	TFile * fout=new TFile(outfilename,"recreate");

	TH1I* hRunNumber=new TH1I("hRunNumber","hRunNumber;run",1,runNumber,runNumber+1);
	hRunNumber->Fill(runNumber);

	//-------------------------------------------------------------
	/*GENERAL MONITOR - MULTIPLICITY*/
	TH1F*hMulti= (TH1F*) generalList->At(0);

	TH1F* hFractionEventsWhits=new TH1F("hFractionEventsWhits","hFractionEventsWhits;fraction of events with hits (%)",200,0.,100.);
	Float_t fraction=0.0;
	if (hMulti->GetEntries()>0.0) fraction=((Float_t) hMulti->GetBinContent(1))/((Float_t) hMulti->GetEntries());
	else fraction=0.0;
	hFractionEventsWhits->Fill(fraction*100.);

	//-------------------------------------------------------------
	/*GENERAL MONITOR - TIMING AND GEOMETRY*/
	TH1F * hTime = (TH1F*) generalList->At(1);
	hTime->SetMarkerStyle(8);
	hTime->SetMarkerSize(0.7);
	hTime->SetMarkerColor(kBlue);
	hTime->SetLineColor(kBlue);
	hTime->SetFillColor(kBlue);
	hTime->SetFillStyle(3001);
	hTime->Rebin(2);
	hTime->GetYaxis()->SetTitle("matched tracks"); 
	hTime->GetYaxis()->SetTitleOffset(1.35);
  
	TH1F * hRawTime = (TH1F*) generalList->At(2);
	hRawTime->SetMarkerStyle(21);
	hRawTime->SetMarkerSize(0.7);
	hRawTime->SetMarkerColor(kGreen);
	hRawTime->SetLineColor(kGreen);
	hRawTime->SetFillColor(kGreen);
	hRawTime->SetFillStyle(3001); 
	hRawTime->Rebin(2);
	hRawTime->GetYaxis()->SetTitle("matched tracks");
	hRawTime->GetYaxis()->SetTitleOffset(1.35);
  
	TLegend *lTime = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
	lTime->SetTextSize(0.04281433);
	lTime->AddEntry(hRawTime, "raw","L");
	lTime->AddEntry(hTime, "ESD","L"); 
	lTime->SetFillColor(kWhite);
	lTime->SetShadowColor(0);

	TH1F * hTot = (TH1F*) generalList->At(3);
	hTot->SetMarkerStyle(8);
	hTot->SetMarkerSize(0.7);
	hTot->SetMarkerColor(kViolet-3);
	hTot->SetLineColor(kViolet-3);
	hTot->SetFillColor(kViolet-3);
	hTot->SetFillStyle(3001);
	//hTime->SetDrawOption();
	hTot->GetYaxis()->SetTitle("matched tracks");
	hTot->GetYaxis()->SetTitleOffset(1.35);
  
  
	char orphansTxt[200];
	Float_t orphansRatio=0.0;
	if (hTot->GetEntries()>1){
		orphansRatio=((Float_t) hTot->GetBinContent(1))/((Float_t) hTot->GetEntries()) ;
	}
	sprintf(orphansTxt,"orphans/matched tracks = %.4f ",orphansRatio);
	TH1F * hOrphansRatio=new TH1F("hOrphansRatio","Percentage of signals with only leading edge; percentage (%)",1000,0.,100.);
	hOrphansRatio->Fill(orphansRatio*100.);  


	TPaveText *tOrphans = new TPaveText(0.38,0.63,0.88,0.7, "NDC");
	//NDC sets coords relative to pad
	tOrphans->SetBorderSize(0);
	tOrphans->SetTextSize(0.045);
	tOrphans->SetFillColor(0); //white background
	tOrphans->SetTextAlign(12);
	tOrphans->SetTextColor(kViolet-3);
	tOrphans->AddText(orphansTxt);
      
	TH1F * hL = (TH1F*) generalList->At(4);
	hL->SetMarkerStyle(8);
	hL->SetMarkerSize(0.7);
	hL->SetMarkerColor(kOrange-3);
	hL->SetLineColor(kOrange-3);
	hL->SetFillColor(kOrange-3);
	hL->SetFillStyle(3001);
	//hTime->SetDrawOption();
	hL->GetYaxis()->SetTitle("matched tracks");
	hL->GetYaxis()->SetTitleOffset(1.35);

	char negLengthTxt[200];
	Float_t negLengthRatio=0.0;
	if (hL->GetEntries()>1){
		negLengthRatio=(hL->Integral(1,750))/((Float_t) hL->GetEntries()) ;
	}
	sprintf(negLengthTxt,"tracks with L<350cm /matched tracks = %.5f ",negLengthRatio);
	TH1F * hLnegRatio=new TH1F("hLnegRatio","Ratio of TOF-matched tracks with L<350cm; ratio (%)",10000,0.,100.);
	hLnegRatio->Fill(negLengthRatio*100);
  


	TPaveText *tLength = new TPaveText(0.15,0.83,0.65,0.87, "NDC");
	//NDC sets coords relative to pad
	tLength->SetBorderSize(0);
	tLength->SetTextSize(0.04);
	tLength->SetFillColor(0); //white background
	tLength->SetTextAlign(11);
	tLength->SetTextColor(kOrange-3);
	tLength->AddText(negLengthTxt);
 
	if (canvasE){
		TCanvas *cTrackProperties= new TCanvas("cTrackProperties","summary of matched tracks properties",900,900);
		cTrackProperties->Divide(2,2);
		cTrackProperties->cd(1);
		gPad->SetLogy();
		gPad->SetGridx();
		gPad->SetGridy();
		hTime->Draw("BAR");
		hRawTime ->Draw("BARsame");
		lTime->Draw();  
		cTrackProperties->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		hTot->Draw("BAR");
		tOrphans->Draw(); 
		cTrackProperties->cd(3);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		hL->Draw("BAR");
		tLength->Draw(); 
	}
	fout->cd();
	hRunNumber->Write();
	hMulti->Write();
	hFractionEventsWhits->Write();
	hTime->Write();
	hRawTime->Write();
	hTot->Write();
	hOrphansRatio->Write();
	hL->Write();
	hLnegRatio->Write();
  
	TH2F* hDxPos4profile = (TH2F*) generalList->At(14);
	TH2F* hDxNeg4profile = (TH2F*) generalList->At(13);
    
	char profilename[30];
	const Int_t ybinMin = 0;
	const Int_t ybinMax =hDxPos4profile->GetYaxis()->GetNbins() ;
	sprintf(profilename,"profDxPos");
	TProfile * profDxPos = (TProfile*)hDxPos4profile->ProfileX(profilename, ybinMin,ybinMax); 
	sprintf(profilename,"profDxNeg");
	profDxPos->SetLineWidth(2);
	TProfile * profDxNeg = (TProfile*)hDxNeg4profile->ProfileX(profilename, ybinMin, ybinMax); 
	profDxNeg->SetLineWidth(2);  
 
	TH1 *profRatioPosOverNegDx = (TH1*) profDxPos->Clone();
	profRatioPosOverNegDx->SetName("profRatioPosOverNegDx");
	profRatioPosOverNegDx->Divide((TH1*) profDxNeg);
	profRatioPosOverNegDx->GetYaxis()->SetRangeUser(-5.,5.);
	profRatioPosOverNegDx->GetXaxis()->SetRangeUser(0.,2.);
	if (canvasE){
		TCanvas *residuals= new TCanvas("residuals","residuals",900,450);
		residuals->Divide(2,1);
		residuals->cd(1);
		gPad->SetGridx();
		gPad->SetGridy();
		hDxPos4profile->GetYaxis()->SetRangeUser(-5.,5.);
		hDxPos4profile->Draw("colz");
		profDxPos->SetLineColor(kRed);
		profDxPos ->Draw("same");
		residuals->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		hDxNeg4profile->GetYaxis()->SetRangeUser(-5.,5.); 
		hDxNeg4profile->Draw("colz");
		profDxNeg->SetLineColor(kBlue);
		profDxNeg->Draw("same"); 
	}

	fout->cd();
	hDxPos4profile->Write();
	hDxNeg4profile->Write();
	profDxPos->Write();
	profDxNeg->Write();
	profRatioPosOverNegDx->Write();
	//-------------------------------------------------------------
	/* T0 DETECTOR MONITOR*/

	TH1F * hT0AC = (TH1F*) timeZeroList->At(0);
	hT0AC->SetMarkerStyle(8);
	hT0AC->SetMarkerSize(0.7);
	hT0AC->SetMarkerColor(kRed);
	hT0AC->SetLineColor(kRed);
	hT0AC->SetFillColor(kRed);
	hT0AC->SetFillStyle(1001);
	hT0AC->Rebin(2);
	hT0AC->GetYaxis()->SetTitle("events"); 
	hT0AC->GetYaxis()->SetTitleOffset(1.35);
	hT0AC->GetXaxis()->SetLabelSize(0.03);
 
	TH1F * hT0A = (TH1F*) timeZeroList->At(1);
	hT0A->SetMarkerStyle(8);
	hT0A->SetMarkerSize(0.7);
	hT0A->SetMarkerColor(kBlue);
	hT0A->SetLineColor(kBlue);
	hT0A->SetFillColor(kBlue);
	hT0A->SetFillStyle(1001);
	hT0A->Rebin(2);
	hT0A->GetYaxis()->SetTitle("events"); 
	hT0A->GetYaxis()->SetTitleOffset(1.35);
	hT0A->GetXaxis()->SetLabelSize(0.03);

	TH1F * hT0C = (TH1F*) timeZeroList->At(2);
	hT0C->SetMarkerStyle(8);
	hT0C->SetMarkerSize(0.7);
	hT0C->SetMarkerColor(kGreen);
	hT0C->SetLineColor(kGreen);
	hT0C->SetFillColor(kGreen);
	hT0C->SetFillStyle(1001);
	hT0C->Rebin(2);
	hT0C->GetYaxis()->SetTitle("events"); 
	hT0C->GetYaxis()->SetTitleOffset(1.35);
	hT0C->GetXaxis()->SetLabelSize(0.03);
 
	TLegend *lT0 = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
	lT0->SetTextSize(0.041);
	lT0->AddEntry(hT0AC, "T0 A&C","L");
	lT0->AddEntry(hT0A, "T0 A","L"); 
	lT0->AddEntry(hT0C, "T0 C","L");
	lT0->SetFillColor(kWhite);
	lT0->SetShadowColor(0);
  
	TH1F * hT0res = (TH1F*) timeZeroList->At(3);
	hT0res->GetXaxis()->SetLabelSize(0.03);
	if (canvasE){
		TCanvas *cT0detector= new TCanvas("cT0detector","T0 detector",900,450);
		cT0detector->Divide(2,1);
		cT0detector->cd(1);
		gPad->SetGridx();
		gPad->SetGridy();
		hT0AC->Draw("BAR");
		hT0AC->SetTitle("timeZero measured by T0 detector");
		hT0A ->Draw("BARsame");
		hT0C ->Draw("BARsame");
		lT0->Draw();  
		cT0detector->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		// gPad->SetLogy();
		hT0res->Draw();
		// myText1->Draw(); 
		// cTrackProperties->cd(3);
		// gPad->SetLogy();
		// hL->Draw("BAR");
		// myText2->Draw(); 
	}

	TH1F * hT0fillRes = (TH1F*) timeZeroList->At(8);
	hT0fillRes->GetXaxis()->SetLabelSize(0.03);

	fout->cd();
	hT0AC->Write();
	hT0A->Write();
	hT0C->Write();
	hT0res->Write();
	hT0fillRes->Write();
	lT0->Write();
  

	//-------------------------------------------------------------
	/*MATCHING EFFICIENCY  MONITOR*/

	//TH1F * hMatchingVsP =new TH1F("hMatchingVsP","Matching probability vs. P; P(GeV/c); matching probability", 50, 0., 5. );
  
	TH1F * hMatchingVsPt =new TH1F("hMatchingVsPt","Matching probability vs. Pt; Pt(GeV/c); matching probability", 50, 0., 5. );
  
	TH1F * hMatchingVsEta =new TH1F("hMatchingVsEta","Matching probability vs. #\Eta; #\Eta; matching probability", 20, -1., 1.);
  
	TH1F * hMatchingVsPhi =new TH1F("hMatchingVsPhi","Matching probability vs. Phi; Phi(rad); matching probability", 628, 0., 6.28);
  
	/*/matching as function of p
	  TH1F * hDenom=(TH1F*)generalList->At(9); 
	  if (hDenom) {  
	  hMatchingVsP=(TH1F*) generalList->At(5); 
	  hMatchingVsP->Rebin(5);
	  hDenom->Rebin(5);
	  hMatchingVsP->Divide(hDenom);
	  hMatchingVsP->GetYaxis()->SetTitle("matching efficiency");
	  hMatchingVsP->SetTitle("TOF matching efficiency as function of momentum");
	  }*/
  
	//matching as function of pT
  
	// hDenom->Clear();
	TH1F * hDenom=(TH1F*)generalList->At(10); 
	if (hDenom) {  
		hMatchingVsPt=(TH1F*) generalList->At(6); 
		hMatchingVsPt->Rebin(5);
		hDenom->Rebin(5);
		hMatchingVsPt->Divide(hDenom);
		hMatchingVsPt->GetYaxis()->SetTitle("matching efficiency");
		hMatchingVsPt->SetTitle("TOF matching efficiency as function of transverse momentum");
		hMatchingVsPt->GetYaxis()->SetRangeUser(0,1.2); 
	}
  
	//matching as function of eta
	hDenom->Clear();
	hDenom=(TH1F*)generalList->At(11); 
	if (hDenom) {  
		hMatchingVsEta=(TH1F*) generalList->At(7); 
		hMatchingVsEta->Rebin(5);
		hDenom->Rebin(5);
		hMatchingVsEta->Divide(hDenom);
		hMatchingVsEta->GetXaxis()->SetRangeUser(-1,1);
		hMatchingVsEta->GetYaxis()->SetTitle("matching efficiency");
		hMatchingVsEta->GetYaxis()->SetRangeUser(0,1.2);
		hMatchingVsEta->SetTitle("TOF matching efficiency as function of pseudorapidity");
	}
	//matching as function of phi
	hDenom->Clear();
	hDenom=(TH1F*)generalList->At(12); 
	if (hDenom) {  
		hMatchingVsPhi=(TH1F*) generalList->At(8); 
		//hMatchingVsPhi->Rebin(5);
		//hDenom->Rebin(5);
		hMatchingVsPhi->Divide(hDenom);
		hMatchingVsPhi->GetYaxis()->SetTitle("matching efficiency");
		hMatchingVsPhi->SetTitle("TOF matching efficiency as function of phi");
		hMatchingVsPhi->GetYaxis()->SetRangeUser(0,1.2);
	}
	if (  canvasE){
		TCanvas *cMatchingPerformance= new TCanvas("cMatchingPerformance","summary of matching performance",700,400);
		cMatchingPerformance->Divide(2,2);
		cMatchingPerformance->cd(1);
		gPad->SetGridx();
		gPad->SetGridy();
		hMatchingVsPt->Draw();
		cMatchingPerformance->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		hMatchingVsEta->Draw();
		cMatchingPerformance->cd(3);
		gPad->SetGridx();
		gPad->SetGridy();
		hMatchingVsPhi->Draw();
	}
	fout->cd();
	hMatchingVsPt->Write();
	hMatchingVsEta->Write();
	hMatchingVsPhi->Write();
  
	//----------------------------------------------------
	/* PID PERFORMANCE MONITOR */

	TH2F * hBetaP=(TH2F*)pidList->At(0);
	hBetaP->GetYaxis()->SetRangeUser(0.,1.2);
	TH1F * hMass=(TH1F*)pidList->At(1);
	//hMass->SetMarkerColor(kBlue);
	//hMass->SetLineColor(kBlue);
	hMass->SetFillColor(kAzure+10);
	hMass->SetFillStyle(1001);
	hMass->Rebin(2);
	hMass->GetYaxis()->SetTitle("tracks"); 
	hMass->GetYaxis()->SetTitleOffset(1.35);
	hMass->GetXaxis()->SetLabelSize(0.03);

	TH1F * hPionDiff=(TH1F*)pidList->At(3); 
  

	TH2F * hDiffTimePi=(TH2F*)pidList->At(4); 
	hDiffTimePi->GetYaxis()->SetRangeUser(-2000.,2000.);
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimePi");
 
	TProfile * profDiffTimePi = (TProfile*)hDiffTimePi->ProfileX(profilename, 490, 510); 
	profDiffTimePi->SetLineWidth(2);
	profDiffTimePi->SetLineColor(kRed+2); 

	TH2F * hDiffTimePiTh=(TH2F*)pidList->At(6); 
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimePiTh");
	hDiffTimePiTh->GetYaxis()->SetRangeUser(-2000.,2000.); 
  
	TProfile * profDiffTimePiTh = (TProfile*)hDiffTimePiTh->ProfileX(profilename, 490, 510);
	profDiffTimePiTh->SetLineWidth(2);
	profDiffTimePiTh->SetLineColor(kRed+2);

	TH2F * hDiffTimeKa=(TH2F*)pidList->At(9); 
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimeKa");
	hDiffTimeKa->GetYaxis()->SetRangeUser(-2000.,2000.);
  
	TProfile * profDiffTimeKa = (TProfile*)hDiffTimeKa->ProfileX(profilename, 490, 510); 
	profDiffTimeKa->SetLineWidth(2);
	profDiffTimeKa->SetLineColor(kBlue);  

	TH2F * hDiffTimeKaTh=(TH2F*)pidList->At(11); 
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimeKaTh");
	hDiffTimeKaTh->GetYaxis()->SetRangeUser(-2000.,2000.);
	TProfile * profDiffTimeKaTh = (TProfile*)hDiffTimeKaTh->ProfileX(profilename, 490, 510); 
	profDiffTimeKaTh->SetLineWidth(2);
	profDiffTimeKaTh->SetLineColor(kBlue);
  
	TH2F * hDiffTimePro=(TH2F*)pidList->At(14); 
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimePro");
	hDiffTimePro->GetYaxis()->SetRangeUser(-2000.,2000.);
	TProfile * profDiffTimePro = (TProfile*)hDiffTimePro->ProfileX(profilename, 490, 510); 
	profDiffTimePro->SetLineWidth(2);
	profDiffTimePro->SetLineColor(kGreen+2);  

	TH2F * hDiffTimeProTh=(TH2F*)pidList->At(16); 
	//hDiffTime->GetYaxis()->Rebin(2);//1 bin=10 ps
	sprintf(profilename,"profDiffTimeProTh");
	hDiffTimeProTh->GetYaxis()->SetRangeUser(-2000.,2000.);
	TProfile * profDiffTimeProTh = (TProfile*)hDiffTimeProTh->ProfileX(profilename, 490, 510);
	profDiffTimeProTh->SetLineWidth(2);
	profDiffTimeProTh->SetLineColor(kGreen+2);

	if (canvasE){
		TCanvas *cPidPerformance= new TCanvas("cPidPerformance","summary of pid performance",800,800);
		cPidPerformance->Divide(2,1);
		cPidPerformance->cd(1);
		gPad->SetGridy();
		gPad->SetGridx();
		gPad->SetLogz();
		hBetaP->Draw("colz"); 
  
		cPidPerformance->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		hMass->Draw("HIST BAR");
  
		TLegend * lPid=new TLegend(0.75,0.75,0.95,0.95,"PID");
		lPid->AddEntry(profDiffTimePi,"#pi^{#pm}","l");
		lPid->AddEntry(profDiffTimeKa,"K^{#pm}","l");
		lPid->AddEntry(profDiffTimePro,"p^{#pm}","l");

		gStyle->SetOptStat(10);
		TCanvas *cPidPerformance2= new TCanvas("cPidPerformance2","summary of pid performance - tracking",700,700);
		cPidPerformance2->Divide(2,2);
		cPidPerformance2->cd(1);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hDiffTimePi->Draw("colz");
		// profDiffTimePi->Draw("same");
  
		cPidPerformance2->cd(2);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hDiffTimeKa->Draw("colz");
		//profDiffTimeKa->Draw("same");
		cPidPerformance2->cd(3);
		gPad->SetLogz();
		hDiffTimePro->Draw("colz");
		//  profDiffTimePro->Draw("same");
 
		// cPidPerformance2->cd(4);
		//  profDiffTimePi->Draw();
		//  profDiffTimeKa->Draw("same");
		//  profDiffTimePro->Draw("same");
		//  lPid->Draw("same");
  
		TCanvas *cPidPerformanceTh= new TCanvas("cPidPerformanceTh","summary of pid performance - theoretical times",700,700);
		cPidPerformanceTh->Divide(2,2);
		cPidPerformanceTh->cd(1);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hDiffTimePiTh->Draw("colz");
		profDiffTimePiTh->Draw("same");
		cPidPerformanceTh->cd(2);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hDiffTimeKaTh->Draw("colz");
		profDiffTimeKaTh->Draw("same");
		cPidPerformanceTh->cd(3);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hDiffTimeProTh->Draw("colz");
		profDiffTimeProTh->Draw("same");
		// cPidPerformanceTh->cd(4);
		// profDiffTimePiTh->Draw();
		// profDiffTimeKaTh->Draw("same");
		// profDiffTimeProTh->Draw("same");
		// lPid->Draw("same");
	}
  
	TH1F * hPionDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimePi"); 
	TH1F * hKaonDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimeKa"); 
	TH1F * hProtonDiff=(TH1F*)pidList->FindObject("hTOFmatchedExpTimePro"); 
 

	fout->cd();
	hPionDiff->Write();
	hKaonDiff->Write();
	hProtonDiff->Write();

	hBetaP->Write();
	hMass->Write();
	hDiffTimePi->Write();
	profDiffTimePi->Write();
	hDiffTimeKa->Write();
	profDiffTimeKa->Write();
	hDiffTimePro->Write();
	profDiffTimePro->Write();
	//lPid->Draw();
	hDiffTimePiTh->Write();
	profDiffTimePiTh->Write();
	hDiffTimeKaTh->Write();
	profDiffTimeKaTh->Write();
	hDiffTimeProTh->Write();
	profDiffTimeProTh->Write();
  
	//SIGMAS PID
	TH2F * hSigmaPi=(TH2F*)pidList->At(7); 
	sprintf(profilename,"profSigmaPi");
	hSigmaPi->GetYaxis()->SetRangeUser(-5.,5.);
	TProfile * profSigmaPi = (TProfile*)hSigmaPi->ProfileX(profilename); 
	profSigmaPi->SetLineWidth(2);
	profSigmaPi->SetLineColor(kRed+2); 

	TH2F * hSigmaKa=(TH2F*)pidList->At(12); 
	sprintf(profilename,"profSigmaKa");
	hSigmaKa->GetYaxis()->SetRangeUser(-5.,5.);
	TProfile * profSigmaKa = (TProfile*)hSigmaKa->ProfileX(profilename); 
	profSigmaKa->SetLineWidth(2);
	profSigmaKa->SetLineColor(kBlue);  

	TH2F * hSigmaPro=(TH2F*)pidList->At(17); 
	sprintf(profilename,"profSigmaPro");
	hSigmaPro->GetYaxis()->SetRangeUser(-5.,5.);
	TProfile * profSigmaPro = (TProfile*)hSigmaPro->ProfileX(profilename); 
	profSigmaPro->SetLineWidth(2);
	profSigmaPro->SetLineColor(kGreen+2);  


	if (canvasE){
	  
		TLegend * lSigmaPid=new TLegend(0.75,0.75,0.95,0.95,"#sigma_{PID}");
		lSigmaPid->AddEntry(profSigmaPi,"#pi^{#pm}","l");
		lSigmaPid->AddEntry(profSigmaKa,"K^{#pm}","l");
		lSigmaPid->AddEntry(profSigmaPro,"p^{#pm}","l");
	  
		TCanvas *cPidPerformance3= new TCanvas("cPidPerformance3","summary of pid performance - sigmas",1200,400);
		cPidPerformance3->Divide(3,1);
		cPidPerformance3->cd(1);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
	  
		hSigmaPi->Draw("colz");
		profSigmaPi->Draw("same");
	  
		cPidPerformance3->cd(2);
		gPad->SetLogz();
		gPad->SetGridx();
		gPad->SetGridy();
		hSigmaKa->Draw("colz");
		profSigmaKa->Draw("same");
	  
		cPidPerformance3->cd(3);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogz();
		hSigmaPro->Draw("colz");
		profSigmaPro->Draw("same");
	}

	fout->cd();
	hSigmaPi->Write();
	profSigmaPi->Write();
	hSigmaKa->Write();
	profSigmaKa->Write();
	hSigmaPro->Write();
	profSigmaPro->Write();

	if (canvasE){
		TH2F* hTOFmatchedDzVsStrip = (TH2F*)generalList->FindObject("hTOFmatchedDzVsStrip");
		TCanvas* cProfile = new TCanvas("cProfile","cProfile",50,50,750,550);
		gStyle->SetOptStat(0);
		hTOFmatchedDzVsStrip->Draw("colz");
		Int_t binmin = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(-3);
		Int_t binmax = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(3);
		TProfile* hDzProfile = (TProfile*)hTOFmatchedDzVsStrip->ProfileX("hDzProfile",binmin, binmax);
		hDzProfile->SetLineWidth(3);
		hDzProfile->Draw("same");
		cProfile->SetGridx();
		cProfile->SetGridy();
		TString plotDir(Form("Plots_run%d",runNumber));
		gSystem->Exec(Form("mkdir %s",plotDir.Data()));
		cPidPerformance3->Print(Form("%s/PID_sigmas.png",plotDir.Data()));
		cPidPerformance->Print(Form("%s/PID.png",plotDir.Data()));
		cPidPerformanceTh->Print(Form("%s/PID_thereticalTimes.png",plotDir.Data()));
		cPidPerformance2->Print(Form("%s/PID_ExpTimes.png",plotDir.Data()));
		cMatchingPerformance->Print(Form("%s/Matching.png",plotDir.Data()));
		cT0detector->Print(Form("%s/T0Detector.png",plotDir.Data()));
		cTrackProperties->Print(Form("%s/TrackProperties.png",plotDir.Data()));
		residuals->Print(Form("%s/Residuals.png",plotDir.Data()));
		cProfile->Print(Form("%s/ProfileDZvsStripNumber.png",plotDir.Data()));
	}

	return 0;
}

//----------------------------------------------------------
char * SetQAtrainOutputName(Int_t run=0,Int_t year=2011,char *period="LHC11a", char* pass="cpass1",Int_t trainId=76){
  
	char infile[200];
	sprintf(infile,"alien:///alice/data/%i/%s/000%d/ESDs/%s/QA%i/QAresults.root",year,period,runNumber,pass,trainId);
	return infile;
  
}



