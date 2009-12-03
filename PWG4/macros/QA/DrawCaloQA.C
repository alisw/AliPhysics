//Histograms
//CaloClusters 
TH1F * fhE  ; //! E distribution, Reco
TH1F * fhPt ; //! pT distribution, Reco
TH1F * fhPhi; //! phi distribution, Reco 
TH1F * fhEta; //! eta distribution, Reco 
TH2F * fhEtaPhi  ; //! eta vs phi, Reco 
TH1F * fhECharged  ; //! E distribution, Reco, matched with track
TH1F * fhPtCharged ; //! pT distribution, Reco, matched with track
TH1F * fhPhiCharged; //! phi distribution, Reco, matched with track 
TH1F * fhEtaCharged; //! eta distribution, Reco, matched with track 
TH2F * fhEtaPhiCharged  ; //! eta vs phi, Reco, matched with track
TH1F * fhEChargedNoOut  ; //! E distribution, Reco, matched with track, no outer param
TH1F * fhPtChargedNoOut ; //! pT distribution, Reco, matched with track, no outer param
TH1F * fhPhiChargedNoOut; //! phi distribution, Reco, matched with track, no outer param
TH1F * fhEtaChargedNoOut; //! eta distribution, Reco, matched with track, no outer param
TH2F * fhEtaPhiChargedNoOut  ; //! eta vs phi, Reco, matched with track, no outer param
TH1F * fhDeltaE  ; //! MC-Reco E distribution	
TH1F * fhDeltaPt ; //! MC-Reco pT distribution
TH1F * fhDeltaPhi; //! MC-Reco phi distribution
TH1F * fhDeltaEta; //! MC-Reco eta distribution
TH1F * fhRatioE  ; //! Reco/MC E distribution	
TH1F * fhRatioPt ; //! Reco/MC pT distribution
TH1F * fhRatioPhi; //! Reco/MC phi distribution
TH1F * fhRatioEta; //! Reco/MC eta distribution
TH2F * fh2E  ; //! E distribution, Reco vs MC
TH2F * fh2Pt ; //! pT distribution, Reco vs MC
TH2F * fh2Phi; //! phi distribution, Reco vs MC
TH2F * fh2Eta; //! eta distribution, Reco vs MC
TH2F * fhIM; //! cluster pairs invariant mass
TH2F * fhAsym; //! cluster pairs invariant mass	
TH2F * fhNCellsPerCluster; //! N cells per cluster	
TH1F * fhNClusters; //! Number of clusters

//Calo Cells
TH1F * fhNCells; //! Number of towers/crystals with signal
TH1F * fhAmplitude; //! Amplitude measured in towers/crystals

//Calorimeters Correlation
TH2F * fhCaloCorrNClusters; // EMCAL vs PHOS, number of clusters	
TH2F * fhCaloCorrEClusters; // EMCAL vs PHOS, total measured cluster energy
TH2F * fhCaloCorrNCells; // EMCAL vs PHOS, number of cells
TH2F * fhCaloCorrECells; // EMCAL vs PHOS,  total measured cell energy

//MC  
TH1F *fhGenGamPt  ; // pt of primary gamma
TH1F *fhGenGamEta ; // eta of primart gamma
TH1F *fhGenGamPhi ; // phi of primary gamma	
TH1F *fhGenPi0Pt  ; // pt of primary pi0
TH1F *fhGenPi0Eta ; // eta of primart pi0
TH1F *fhGenPi0Phi ; // phi of primary pi0	
TH1F *fhGenEtaPt  ; // pt of primary eta
TH1F *fhGenEtaEta ; // eta of primart eta
TH1F *fhGenEtaPhi ; // phi of primary eta
TH1F *fhGenOmegaPt  ; // pt of primary omega
TH1F *fhGenOmegaEta ; // eta of primart omega
TH1F *fhGenOmegaPhi ; // phi of primary omega
TH1F *fhGenElePt  ; // pt of primary electron
TH1F *fhGenEleEta ; // eta of primart electron
TH1F *fhGenElePhi ; // phi of primary electron

//TH3F * fhEMVxyz    ; // Electromagnetic particle production vertex
TH2F * fhEMVxyz    ; // Electromagnetic particle production vertex
TH2F * fhEMR       ; // Electromagnetic distance to vertex vs rec energy  
//TH3F * fhHaVxyz    ; // Hadron production vertex
TH2F * fhHaVxyz    ; // Hadron production vertex
TH2F * fhHaR       ; // Hadron distance to vertex vs rec energy  

TH2F * fhGamE  ; //! E distribution of generated photons, Reco
TH2F * fhGamPt ; //! pT distribution of generated photons, Reco
TH2F * fhGamPhi; //! phi distribution of generated photon, Reco 
TH2F * fhGamEta; //! eta distribution of generated photons, Reco 	
TH1F * fhGamDeltaE  ; //! MC-Reco E distribution of generated photons	
TH1F * fhGamDeltaPt ; //! MC-Reco pT distribution of generated photons
TH1F * fhGamDeltaPhi; //! MC-Reco phi distribution of generated photons
TH1F * fhGamDeltaEta; //! MC-Reco eta distribution of generated photons
TH1F * fhGamRatioE  ; //! Reco/MC E distribution of generated photons	
TH1F * fhGamRatioPt ; //! Reco/MC pT distribution of generated photons
TH1F * fhGamRatioPhi; //! Reco/MC phi distribution of generated photons
TH1F * fhGamRatioEta; //! Reco/MC eta distribution of generated photons
TH2F * fhEleE  ; //! E distribution of generated electrons, Reco
TH2F * fhElePt ; //! pT distribution of generated electrons, Reco
TH2F * fhElePhi; //! phi distribution of generated electron, Reco 
TH2F * fhEleEta; //! eta distribution of generated electrons, Reco 		
TH2F * fhPi0E  ; //! E distribution of generated pi0, Reco, gamma decay overlapped
TH2F * fhPi0Pt ; //! pT distribution of generated pi0, Reco, gamma decay overlapped
TH2F * fhPi0Phi; //! phi distribution of generated pi0, Reco, gamma decay overlapped
TH2F * fhPi0Eta; //! eta distribution of generated pi0, Reco, gamma decay overlapped
TH2F * fhNeHadE  ; //! E distribution of generated neutral hadron, Reco
TH2F * fhNeHadPt ; //! pT distribution of generated neutral hadron, Reco
TH2F * fhNeHadPhi; //! phi distribution of generated neutral hadron, Reco 
TH2F * fhNeHadEta; //! eta distribution of generated neutral hadron, Reco 	
TH2F * fhChHadE  ; //! E distribution of generated charged hadron, Reco
TH2F * fhChHadPt ; //! pT distribution of generated charged hadron, Reco
TH2F * fhChHadPhi; //! phi distribution of generated charged hadron, Reco 
TH2F * fhChHadEta; //! eta distribution of generated charged hadron, Reco 

TH2F * fhGamECharged  ; //! E distribution of generated photons, Reco, track matched cluster
TH2F * fhGamPtCharged ; //! pT distribution of generated photons, Reco, track matched cluster
TH2F * fhGamPhiCharged; //! phi distribution of generated photon, Reco, track matched cluster 
TH2F * fhGamEtaCharged; //! eta distribution of generated photons, Reco, track matched cluster 
TH2F * fhEleECharged  ; //! E distribution of generated electrons, Reco, track matched cluster
TH2F * fhElePtCharged ; //! pT distribution of generated electrons, Reco, track matched cluster
TH2F * fhElePhiCharged; //! phi distribution of generated electron, Reco, track matched cluster 
TH2F * fhEleEtaCharged; //! eta distribution of generated electrons, Reco, track matched cluster 		
TH2F * fhPi0ECharged  ; //! E distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
TH2F * fhPi0PtCharged ; //! pT distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
TH2F * fhPi0PhiCharged; //! phi distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
TH2F * fhPi0EtaCharged; //! eta distribution of generated pi0, Reco, gamma decay overlapped, track matched cluster
TH2F * fhNeHadECharged  ; //! E distribution of generated neutral hadron, Reco, track matched cluster
TH2F * fhNeHadPtCharged ; //! pT distribution of generated neutral hadron, Reco, track matched cluster
TH2F * fhNeHadPhiCharged; //! phi distribution of generated neutral hadron, Reco , track matched cluster
TH2F * fhNeHadEtaCharged; //! eta distribution of generated neutral hadron, Reco, track matched cluster 	
TH2F * fhChHadECharged  ; //! E distribution of generated charged hadron, Reco, track matched cluster
TH2F * fhChHadPtCharged ; //! pT distribution of generated charged hadron, Reco, track matched cluster
TH2F * fhChHadPhiCharged; //! phi distribution of generated charged hadron, Reco, track matched cluster 
TH2F * fhChHadEtaCharged; //! eta distribution of generated charged hadron, Reco, track matched cluster 	

TH1F *fhGenGamAccE   ; // E of primary gamma
TH1F *fhGenGamAccPt  ; // pt of primary gamma
TH1F *fhGenGamAccEta ; // eta of primart gamma
TH1F *fhGenGamAccPhi ; // phi of primary gamma	
TH1F *fhGenPi0AccE   ; // E of primary pi0
TH1F *fhGenPi0AccPt  ; // pt of primary pi0
TH1F *fhGenPi0AccEta ; // eta of primart pi0
TH1F *fhGenPi0AccPhi ; // phi of primary pi0		

//Histograms for track-matching
TH2F *fh1pOverE;     //! p/E for track-cluster matches
TH1F *fh1dR;         //! distance between projected track and cluster
TH2F *fh2EledEdx;    //! dE/dx vs. momentum for electron candidates
TH2F *fh2MatchdEdx;  //! dE/dx vs. momentum for all matches
TH2F *fhMCEle1pOverE;     //! p/E for track-cluster matches, MC electrons
TH1F *fhMCEle1dR;         //! distance between projected track and cluster, MC electrons
TH2F *fhMCEle2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC electrons	

TH2F *fhMCChHad1pOverE;     //! p/E for track-cluster matches, MC charged hadrons
TH1F *fhMCChHad1dR;         //! distance between projected track and cluster, MC charged hadrons
TH2F *fhMCChHad2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC charged

TH2F *fhMCNeutral1pOverE;     //! p/E for track-cluster matches, MC neutral
TH1F *fhMCNeutral1dR;         //! distance between projected track and cluster, MC neutral
TH2F *fhMCNeutral2MatchdEdx;  //! dE/dx vs. momentum for all matches, MC neutral	

TH2F *fh1pOverER02;           //! p/E for track-cluster matches, dR > 0.2	
TH2F *fhMCEle1pOverER02;      //! p/E for track-cluster matches, dR > 0.2, MC electrons
TH2F *fhMCChHad1pOverER02;    //! p/E for track-cluster matches, dR > 0.2, MC charged hadrons
TH2F *fhMCNeutral1pOverER02;  //! p/E for track-cluster matches, dR > 0.2, MC neutral

//________________________________________________________________________
void ReadHistograms(TString name, Bool_t isDataMC, Bool_t  fCorrelateCalos)
{
	TFile *f = new TFile("Calo.Performance.root","read");
	TList* outputList = f->Get("Calo.Performance");
	
	// Histograms of this analsys are kept in the same list as other analysis, recover the position of
	// the first one and then add the next 
	Int_t index = outputList->IndexOf(outputList->FindObject(name+"hE"));
	//printf("Calo: %s, index: %d\n",fCalorimeter.Data(),index);
	//Read histograms, must be in the same order as in GetCreateOutputObject.
	fhE      = (TH1F *) outputList->At(index++); 	
	fhPt     = (TH1F *) outputList->At(index++); 
	fhPhi    = (TH1F *) outputList->At(index++); 
	fhEta    = (TH1F *) outputList->At(index++);
	fhEtaPhi = (TH2F *) outputList->At(index++);
	
	fhECharged      = (TH1F *) outputList->At(index++); 	
	fhPtCharged     = (TH1F *) outputList->At(index++); 
	fhPhiCharged    = (TH1F *) outputList->At(index++); 
	fhEtaCharged    = (TH1F *) outputList->At(index++);
	fhEtaPhiCharged = (TH2F *) outputList->At(index++);
	
	fhEChargedNoOut      = (TH1F *) outputList->At(index++); 	
	fhPtChargedNoOut     = (TH1F *) outputList->At(index++); 
	fhPhiChargedNoOut    = (TH1F *) outputList->At(index++); 
	fhEtaChargedNoOut    = (TH1F *) outputList->At(index++);
	fhEtaPhiChargedNoOut = (TH2F *) outputList->At(index++);
	
	fh1pOverE    = (TH2F *) outputList->At(index++);
	fh1dR        = (TH1F *) outputList->At(index++);
	fh2MatchdEdx = (TH2F *) outputList->At(index++);
	fh2EledEdx   = (TH2F *) outputList->At(index++);
	fh1pOverER02 = (TH2F *) outputList->At(index++); 

	fhIM     = (TH2F *) outputList->At(index++);
	fhAsym   = (TH2F *) outputList->At(index++);
	
	fhNCellsPerCluster = (TH2F *) outputList->At(index++);
	fhNClusters  = (TH1F *) outputList->At(index++); 
	fhNCells     = (TH1F *) outputList->At(index++); 
	fhAmplitude  = (TH1F *) outputList->At(index++); 
	
	if(fCorrelateCalos){
		fhCaloCorrNClusters = (TH2F *) outputList->At(index++);
		fhCaloCorrEClusters = (TH2F *) outputList->At(index++); 
		fhCaloCorrNCells    = (TH2F *) outputList->At(index++); 
		fhCaloCorrECells    = (TH2F *) outputList->At(index++); 
	}
	

	if(isDataMC){
		fhDeltaE   = (TH1F *) outputList->At(index++); 
		fhDeltaPt  = (TH1F *) outputList->At(index++); 
		fhDeltaPhi = (TH1F *) outputList->At(index++); 
		fhDeltaEta = (TH1F *) outputList->At(index++); 
		
		fhRatioE   = (TH1F *) outputList->At(index++); 
		fhRatioPt  = (TH1F *) outputList->At(index++); 
		fhRatioPhi = (TH1F *) outputList->At(index++); 
		fhRatioEta = (TH1F *) outputList->At(index++); 
		
		fh2E       = (TH2F *) outputList->At(index++); 
		fh2Pt      = (TH2F *) outputList->At(index++); 
		fh2Phi     = (TH2F *) outputList->At(index++); 
		fh2Eta     = (TH2F *) outputList->At(index++); 
		
		fhGamE     = (TH2F *) outputList->At(index++); 
		fhGamPt    = (TH2F *) outputList->At(index++); 
		fhGamPhi   = (TH2F *) outputList->At(index++); 
		fhGamEta   = (TH2F *) outputList->At(index++); 
		
		fhGamDeltaE   = (TH1F *) outputList->At(index++); 
		fhGamDeltaPt  = (TH1F *) outputList->At(index++); 
		fhGamDeltaPhi = (TH1F *) outputList->At(index++); 
		fhGamDeltaEta = (TH1F *) outputList->At(index++); 
		
		fhGamRatioE   = (TH1F *) outputList->At(index++); 
		fhGamRatioPt  = (TH1F *) outputList->At(index++); 
		fhGamRatioPhi = (TH1F *) outputList->At(index++); 
		fhGamRatioEta = (TH1F *) outputList->At(index++); 
				
		fhPi0E     = (TH2F *) outputList->At(index++); 
		fhPi0Pt    = (TH2F *) outputList->At(index++); 
		fhPi0Phi   = (TH2F *) outputList->At(index++); 
		fhPi0Eta   = (TH2F *) outputList->At(index++); 		
		
		fhEleE     = (TH2F *) outputList->At(index++); 
		fhElePt    = (TH2F *) outputList->At(index++); 
		fhElePhi   = (TH2F *) outputList->At(index++); 
		fhEleEta   = (TH2F *) outputList->At(index++); 		
		
		fhNeHadE     = (TH2F *) outputList->At(index++); 
		fhNeHadPt    = (TH2F *) outputList->At(index++); 
		fhNeHadPhi   = (TH2F *) outputList->At(index++); 
		fhNeHadEta   = (TH2F *) outputList->At(index++); 		
		
		fhChHadE     = (TH2F *) outputList->At(index++); 
		fhChHadPt    = (TH2F *) outputList->At(index++); 
		fhChHadPhi   = (TH2F *) outputList->At(index++); 
		fhChHadEta   = (TH2F *) outputList->At(index++); 					
		fhGamECharged     = (TH2F *) outputList->At(index++); 
		fhGamPtCharged    = (TH2F *) outputList->At(index++); 
		fhGamPhiCharged   = (TH2F *) outputList->At(index++); 
		fhGamEtaCharged   = (TH2F *) outputList->At(index++); 
		
		fhPi0ECharged     = (TH2F *) outputList->At(index++); 
		fhPi0PtCharged    = (TH2F *) outputList->At(index++); 
		fhPi0PhiCharged   = (TH2F *) outputList->At(index++); 
		fhPi0EtaCharged   = (TH2F *) outputList->At(index++); 		
		
		fhEleECharged     = (TH2F *) outputList->At(index++); 
		fhElePtCharged    = (TH2F *) outputList->At(index++); 
		fhElePhiCharged   = (TH2F *) outputList->At(index++); 
		fhEleEtaCharged   = (TH2F *) outputList->At(index++); 		
		
		fhNeHadECharged     = (TH2F *) outputList->At(index++); 
		fhNeHadPtCharged    = (TH2F *) outputList->At(index++); 
		fhNeHadPhiCharged   = (TH2F *) outputList->At(index++); 
		fhNeHadEtaCharged   = (TH2F *) outputList->At(index++); 		
		
		fhChHadECharged     = (TH2F *) outputList->At(index++); 
		fhChHadPtCharged    = (TH2F *) outputList->At(index++); 
		fhChHadPhiCharged   = (TH2F *) outputList->At(index++); 
		fhChHadEtaCharged   = (TH2F *) outputList->At(index++); 				
		
//		fhEMVxyz     = (TH3F *) outputList->At(index++); 
//		fhHaVxyz     = (TH3F *) outputList->At(index++); 
		
		fhEMVxyz     = (TH2F *) outputList->At(index++); 
		fhHaVxyz     = (TH2F *) outputList->At(index++); 
		fhEMR        = (TH2F *) outputList->At(index++); 
		fhHaR        = (TH2F *) outputList->At(index++); 
		
		fhGenGamPt    = (TH1F *) outputList->At(index++); 
		fhGenGamEta   = (TH1F *) outputList->At(index++); 
		fhGenGamPhi   = (TH1F *) outputList->At(index++); 
		
		fhGenPi0Pt    = (TH1F *) outputList->At(index++); 
		fhGenPi0Eta   = (TH1F *) outputList->At(index++); 
		fhGenPi0Phi   = (TH1F *) outputList->At(index++); 
		
		fhGenEtaPt    = (TH1F *) outputList->At(index++); 
		fhGenEtaEta   = (TH1F *) outputList->At(index++); 
		fhGenEtaPhi   = (TH1F *) outputList->At(index++); 
		
		fhGenOmegaPt  = (TH1F *) outputList->At(index++); 
		fhGenOmegaEta = (TH1F *) outputList->At(index++); 
		fhGenOmegaPhi = (TH1F *) outputList->At(index++); 
		
		fhGenElePt    = (TH1F *) outputList->At(index++); 
		fhGenEleEta   = (TH1F *) outputList->At(index++); 
		fhGenElePhi   = (TH1F *) outputList->At(index++); 
		
		fhGenGamAccE   = (TH1F *) outputList->At(index++); 		
		fhGenGamAccPt  = (TH1F *) outputList->At(index++); 
		fhGenGamAccEta = (TH1F *) outputList->At(index++); 
		fhGenGamAccPhi = (TH1F *) outputList->At(index++); 
		
		fhGenPi0AccE   = (TH1F *) outputList->At(index++); 		
		fhGenPi0AccPt  = (TH1F *) outputList->At(index++); 
		fhGenPi0AccEta = (TH1F *) outputList->At(index++); 
		fhGenPi0AccPhi = (TH1F *) outputList->At(index++); 
	     
		//Track matching
		fhMCEle1pOverE =    (TH2F *) outputList->At(index++);
		fhMCEle1dR =        (TH1F *) outputList->At(index++);
		fhMCEle2MatchdEdx = (TH2F *) outputList->At(index++);
	
		fhMCChHad1pOverE =    (TH2F *) outputList->At(index++);
		fhMCChHad1dR =        (TH1F *) outputList->At(index++);
		fhMCChHad2MatchdEdx = (TH2F *) outputList->At(index++);
	
		fhMCNeutral1pOverE    = (TH2F *) outputList->At(index++);
		fhMCNeutral1dR        = (TH1F *) outputList->At(index++);
		fhMCNeutral2MatchdEdx = (TH2F *) outputList->At(index++);
	
		fhMCEle1pOverER02     =    (TH2F *) outputList->At(index++);
		fhMCChHad1pOverER02   =    (TH2F *) outputList->At(index++);
		fhMCNeutral1pOverER02 =    (TH2F *) outputList->At(index++);
	}//Is data MC
}


//__________________________________________________________________
void  DrawCaloQA(TString fCalorimeter = "PHOS", Bool_t isDataMC = kFALSE) 
{
	//Macro for replotting Calorimeter QA histograms
	Bool_t fCorrelateCalos = kFALSE;
	//By default when PHOS QA is called correlation not done
	if(fCalorimeter == "PHOS") fCorrelateCalos = kFALSE;
	else fCorrelateCalos = kTRUE;
	cout<<"Calo? "<<fCalorimeter<<" Correlate plots? "<<fCorrelateCalos<<" MC plots?"<<isDataMC<<endl;
	//Do some plots to end
	gROOT->Macro("./style.C");//Set different root style parameters
	//Recover histograms from output histograms list, needed for distributed analysis.	
	ReadHistograms(fCalorimeter+"_", isDataMC, fCorrelateCalos);
	Float_t minx = 0;
	Float_t maxx = 10;
	
	char name[128];
	char cname[128];
	
	//Reconstructed distributions
	//printf("c1\n");
	sprintf(cname,"QA_%s_rec",fCalorimeter.Data());
	TCanvas  * c = new TCanvas(cname, "Reconstructed distributions", 400, 400) ;
	c->Divide(2, 2);
	
	fhE->SetAxisRange(minx,maxx,"X");
	fhPt->SetAxisRange(minx,maxx,"X");
	
	c->cd(1) ; 
	gPad->SetLogy();
	fhE->SetLineColor(4);
	fhE->Draw();
	
	c->cd(2) ; 
	gPad->SetLogy();
	fhPt->SetLineColor(4);
	fhPt->Draw();
	
	c->cd(3) ; 
	fhPhi->SetLineColor(4);
	fhPhi->Draw();
	
	c->cd(4) ; 
	fhEta->SetLineColor(4);
	fhEta->Draw();
	
	sprintf(name,"QA_%s_ReconstructedDistributions.eps",fCalorimeter.Data());
	c->Print(name);
	
	//Reconstructed distributions, matched with tracks
	//printf("c2\n");
	sprintf(cname,"QA_%s_rectrackmatch",fCalorimeter.Data());
	TCanvas  * c2 = new TCanvas(cname, "Reconstructed distributions, matched with tracks", 400, 400) ;
	c2->Divide(2, 2);
	
	fhECharged->SetAxisRange(minx,maxx,"X");
	fhPtCharged->SetAxisRange(minx,maxx,"X");

	c2->cd(1) ; 
	gPad->SetLogy();
	fhECharged->SetLineColor(4);
	fhECharged->Draw();
	
	c2->cd(2) ; 
	gPad->SetLogy();
	fhPtCharged->SetLineColor(4);
	fhPtCharged->Draw();
	
	c2->cd(3) ; 
	fhPhiCharged->SetLineColor(4);
	fhPhiCharged->Draw();
	
	c2->cd(4) ; 
	fhEtaCharged->SetLineColor(4);
	fhEtaCharged->Draw();
	
	sprintf(name,"QA_%s_ReconstructedDistributions_TrackMatched.eps",fCalorimeter.Data());
	c2->Print(name);
	
	TH1F *	hEChargedClone   = (TH1F*)   fhECharged->Clone("EChargedClone");
	TH1F *	hPtChargedClone  = (TH1F*)   fhPtCharged->Clone("PtChargedClone");
	TH1F *	hEtaChargedClone = (TH1F*)   fhEtaCharged->Clone("EtaChargedClone");
	TH1F *	hPhiChargedClone = (TH1F*)   fhPhiCharged->Clone("PhiChargedClone");
	
	TH1F *	hEChargedClone2   = (TH1F*)   fhECharged->Clone("EChargedClone2");
	TH1F *	hPtChargedClone2  = (TH1F*)   fhPtCharged->Clone("PtChargedClone2");
	TH1F *	hEtaChargedClone2 = (TH1F*)   fhEtaCharged->Clone("EtaChargedClone2");
	TH1F *	hPhiChargedClone2 = (TH1F*)   fhPhiCharged->Clone("PhiChargedClone2");
	
	//Ratio: reconstructed track matched/ all reconstructed
	//printf("c3\n");
	sprintf(cname,"QA_%s_rectrackmatchrat",fCalorimeter.Data());
	TCanvas  * c3 = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed", 400, 400) ;
	c3->Divide(2, 2);
	
	c3->cd(1) ;
	gPad->SetLogy();
	hEChargedClone->SetTitleOffset(1.6,"Y");
	hEChargedClone->SetYTitle("track matched / all   ");
	hEChargedClone->SetXTitle("E (GeV)");
	hEChargedClone->Divide(fhE);
	hEChargedClone->Draw();
	
	c3->cd(2) ; 
	gPad->SetLogy();
	hPtChargedClone->SetTitleOffset(1.6,"Y");
	hPtChargedClone->SetYTitle("track matched / all   ");
	hPtChargedClone->SetXTitle("p_{T} (GeV/c)");
	hPtChargedClone->Divide(fhPt);
	hPtChargedClone->Draw();
	
	c3->cd(3) ;
	gPad->SetLogy();
	hPhiChargedClone->SetTitleOffset(1.6,"Y");
	hPhiChargedClone->SetYTitle("track matched / all   ");
	hPhiChargedClone->SetXTitle("#phi (rad)");
	hPhiChargedClone->Divide(fhPhi);
	hPhiChargedClone->Draw();
	
	c3->cd(4) ; 
	gPad->SetLogy();
	hEtaChargedClone->SetTitleOffset(1.6,"Y");
	hEtaChargedClone->SetYTitle("track matched / all   ");
	hEtaChargedClone->SetXTitle("#eta");
	hEtaChargedClone->Divide(fhEta);
	hEtaChargedClone->Draw();
	
	sprintf(name,"QA_%s_RatioReconstructedMatchedDistributions.eps",fCalorimeter.Data());
	c3->Print(name);
	
	//Ratio: reconstructed track matched (minus no track param) / all
	//printf("c333\n");
	sprintf(cname,"QA_%s_rectrackmatchratout",fCalorimeter.Data());
	TCanvas  * c333 = new TCanvas(cname, "Ratio: reconstructed track matched (with outer track param)/ all", 400, 400) ;
	c333->Divide(2, 2);
	
	c333->cd(1) ;
	hEChargedClone2->Add(fhEChargedNoOut,-1);
	hEChargedClone2->SetYTitle("track matched / all");
	hEChargedClone2->SetXTitle("E (GeV)");
	hEChargedClone2->Divide(fhE);
	hEChargedClone2->Draw();
	
	c333->cd(2) ; 
	hPtChargedClone2->Add(fhPtChargedNoOut,-1);
	hPtChargedClone2->SetYTitle("track matched / all");
	hPtChargedClone2->SetXTitle("p_{T} (GeV/c)");
	hPtChargedClone2->Divide(fhPt);
	hPtChargedClone2->Draw();
	
	c333->cd(3) ;
	hPhiChargedClone2->Add(fhPhiChargedNoOut,-1);
	hPhiChargedClone2->SetYTitle("track matched / all");
	hPhiChargedClone2->SetXTitle("#phi (rad)");
	hPhiChargedClone2->Divide(fhPhi);
	hPhiChargedClone2->Draw();
	
	c333->cd(4) ; 
	hEtaChargedClone2->Add(fhEtaChargedNoOut,-1);
	hEtaChargedClone2->SetYTitle("track matched / all");
	hEtaChargedClone2->SetXTitle("#eta");
	hEtaChargedClone2->Divide(fhEta);
	hEtaChargedClone2->Draw();
	
	sprintf(name,"QA_%s_RatioReconstructedMatchedDistributionsOuter.eps",fCalorimeter.Data());
	c333->Print(name);
	
//	//Reconstructed distributions, matched with tracks but no outer param
//	//printf("c2\n");
//	sprintf(cname,"QA_%s_rectrackmatch_noout",fCalorimeter.Data());
//	TCanvas  * c22 = new TCanvas(cname, "Reconstructed distributions, matched with tracks, no outer track param", 400, 400) ;
//	c22->Divide(2, 2);
//	
//	c22->cd(1) ; 
//	gPad->SetLogy();
//	fhEChargedNoOut->SetLineColor(4);
//	fhEChargedNoOut->Draw();
//	
//	c22->cd(2) ; 
//	gPad->SetLogy();
//	fhPtChargedNoOut->SetLineColor(4);
//	fhPtChargedNoOut->Draw();
//	
//	c22->cd(3) ; 
//	fhPhiChargedNoOut->SetLineColor(4);
//	fhPhiChargedNoOut->Draw();
//	
//	c22->cd(4) ; 
//	fhEtaChargedNoOut->SetLineColor(4);
//	fhEtaChargedNoOut->Draw();
//	
//	sprintf(name,"QA_%s_ReconstructedDistributions_TrackMatched_NoOutParam.eps",fCalorimeter.Data());
//	c22->Print(name);
	
	//Ratio: reconstructed track matched/ all reconstructed
	//printf("c3\n");
	
//	TH1F *	hEChargedNoOutClone   = (TH1F*)   fhEChargedNoOut->Clone("EChargedNoOutClone");
//	TH1F *	hPtChargedNoOutClone  = (TH1F*)   fhPtChargedNoOut->Clone("PtChargedNoOutClone");
//	TH1F *	hEtaChargedNoOutClone = (TH1F*)   fhEtaChargedNoOut->Clone("EtaChargedNoOutClone");
//	TH1F *	hPhiChargedNoOutClone = (TH1F*)   fhPhiChargedNoOut->Clone("PhiChargedNoOutClone");	
//	
//	sprintf(cname,"QA_%s_rectrackmatchratnoout",fCalorimeter.Data());
//	TCanvas  * c33 = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed", 400, 400) ;
//	c33->Divide(2, 2);
//	
//	c33->cd(1) ;
//	hEChargedNoOutClone->SetYTitle("track matched no out/ all matched");
//	hEChargedNoOutClone->SetXTitle("E (GeV)");
//	hEChargedNoOutClone->Divide(fhECharged);
//	hEChargedNoOutClone->Draw();
//	
//	c33->cd(2) ; 
//	hPtChargedNoOutClone->SetYTitle("track matched no out / all matched");
//	hPtChargedNoOutClone->SetXTitle("p_{T} (GeV/c)");
//	hPtChargedNoOutClone->Divide(fhPtCharged);
//	hPtChargedNoOutClone->Draw();
//	
//	c33->cd(3) ;
//	hPhiChargedNoOutClone->SetYTitle("track matched no out/ all matched");
//	hPhiChargedNoOutClone->SetXTitle("#phi (rad)");
//	hPhiChargedNoOutClone->Divide(fhPhiCharged);
//	hPhiChargedNoOutClone->Draw();
//	
//	c33->cd(4) ; 
//	hEtaChargedNoOutClone->SetYTitle("track matched no out/ all matched");
//	hEtaChargedNoOutClone->SetXTitle("#eta");
//	hEtaChargedNoOutClone->Divide(fhEtaCharged);
//	hEtaChargedNoOutClone->Draw();
//	
//	sprintf(name,"QA_%s_RatioMatchedDistributionsAllToNoOut.eps",fCalorimeter.Data());
//	c33->Print(name);
	
	

	//TRACK MATCHING P/E distributions
	//printf("cPoverE\n");
	sprintf(cname,"QA_%s_trkmatch",fCalorimeter.Data());
	TCanvas *cme = new TCanvas(cname,"Track-matching distributions", 400, 400);
	cme->Divide(2,2);
	
	TLegend pLegendpE0(0.6,0.55,0.9,0.8);
	pLegendpE0.SetTextSize(0.04);
	pLegendpE0.AddEntry(fh1pOverE,"all","L");
	pLegendpE0.AddEntry(fh1pOverER02,"dR < 0.02","L");
	pLegendpE0.SetFillColor(10);
	pLegendpE0.SetBorderSize(1);
	//pLegendpE0.Draw();
	
	cme->cd(1);
	gPad->SetLogy();
	fh1pOverE->SetTitle("Track matches p/E");
	fh1pOverE->Draw();
	fh1pOverER02->SetLineColor(4);
	fh1pOverER02->Draw("same");
	pLegendpE0.Draw();
	
	cme->cd(2);
	gPad->SetLogy();
	fh1dR->Draw();
	
	cme->cd(3);
	fh2MatchdEdx->Draw();
	
	cme->cd(4);
	fh2EledEdx->Draw();
	
	sprintf(name,"QA_%s_TrackMatchingEleDist.eps",fCalorimeter.Data());
	cme->Print(name);       

	//eta vs phi
	//printf("c4\n");
	sprintf(cname,"QA_%s_etavsphi",fCalorimeter.Data());
//	TCanvas  * c4 = new TCanvas(cname, "reconstructed #eta vs #phi", 600, 200) ;
//	c4->Divide(3, 1);
	
	TCanvas  * c4 = new TCanvas(cname, "reconstructed #eta vs #phi", 400, 200) ;
	c4->Divide(2, 1);
	
	c4->cd(1) ;
	fhEtaPhi->Draw("cont");
	
	c4->cd(2) ; 
	fhEtaPhiCharged->Draw("cont");
	
//	c4->cd(3) ; 
//	fhEtaPhiChargedNoOut->Draw("cont");

	sprintf(name,"QA_%s_ReconstructedEtaVsPhi.eps",fCalorimeter.Data());
	c4->Print(name);
	
	//Invariant mass
	Int_t binmin = -1;
	Int_t binmax = -1;
	//printf("c5\n");
	if(fhIM->GetEntries() > 1){
		Int_t nebins  = fhIM->GetNbinsX();
		Int_t emax = (Int_t) fhIM->GetXaxis()->GetXmax();
		Int_t emin = (Int_t) fhIM->GetXaxis()->GetXmin();
		if (emin != 0 ) printf("emin != 0 \n");
		//printf("IM: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
		
		sprintf(cname,"QA_%s_IM",fCalorimeter.Data());
		//	printf("c5\n");
		TCanvas  * c5 = new TCanvas(cname, "Invariant mass", 400, 400) ;
		c5->Divide(2, 2);
		
		c5->cd(1) ; 
		fhIM->SetLineColor(4);
		fhIM->Draw();
		
		c5->cd(2) ; 
		binmin = 0;
		binmax =  (Int_t) (5-emin)*nebins/emax;
		TH1D *pyim5 = fhIM->ProjectionY("pyim5",binmin,binmax);
		pyim5->SetTitle("E_{pair} < 5 GeV");
		pyim5->SetLineColor(4);
		pyim5->Draw();
		
		c5->cd(3) ; 
		binmin =  (Int_t) (5-emin)*nebins/emax;
		binmax =  (Int_t) (10-emin)*nebins/emax;
		TH1D *pyim510 = fhIM->ProjectionY("pyim5_10",binmin,binmax);
		pyim510->SetTitle("5 < E_{pair} < 10 GeV");
		pyim510->SetLineColor(4);
		pyim510->Draw();
		
		c5->cd(4) ;
		binmin =  (Int_t) (10-emin)*nebins/emax;
		binmax = -1;
		TH1D *pyim10 = fhIM->ProjectionY("pyim10",binmin,binmax);
		pyim10->SetTitle("E_{pair} > 10 GeV");
		pyim10->SetLineColor(4);
		pyim10->Draw();
		
		sprintf(name,"QA_%s_InvariantMass.eps",fCalorimeter.Data());
		c5->Print(name);
	}
	
	//Asymmetry
	//printf("c5b\n");
	if(fhAsym->GetEntries() > 1){
		Int_t nebins  = fhAsym->GetNbinsX();
		Int_t emax = (Int_t) fhAsym->GetXaxis()->GetXmax();
		Int_t emin = (Int_t) fhAsym->GetXaxis()->GetXmin();
		if (emin != 0 ) printf("emin != 0 \n");
		//printf("Asym: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
		
		sprintf(cname,"QA_%s_Asym",fCalorimeter.Data());
		//	printf("c5\n");
		TCanvas  * c5b = new TCanvas(cname, "Asymmetry", 400, 400) ;
		c5b->Divide(2, 2);
		
		c5b->cd(1) ; 
		fhAsym->SetTitleOffset(1.6,"Y");
		fhAsym->SetLineColor(4);
		fhAsym->Draw();
		
		c5b->cd(2) ; 
		binmin = 0;
		binmax = (Int_t) (5-emin)*nebins/emax;
		TH1D *pyAsym5 = fhAsym->ProjectionY("pyAsym5",binmin,binmax);
		pyAsym5->SetTitle("E_{pair} < 5 GeV");
		pyAsym5->SetLineColor(4);
		pyAsym5->Draw();
		
		c5b->cd(3) ; 
		binmin = (Int_t) (5-emin)*nebins/emax;
		binmax = (Int_t) (10-emin)*nebins/emax;
		TH1D *pyAsym510 = fhAsym->ProjectionY("pyAsym5_10",binmin,binmax);
		pyAsym510->SetTitle("5 < E_{pair} < 10 GeV");
		pyAsym510->SetLineColor(4);
		pyAsym510->Draw();
		
		c5b->cd(4) ;
		binmin = (Int_t) (10-emin)*nebins/emax;
		binmax = -1;
		TH1D *pyAsym10 = fhAsym->ProjectionY("pyAsym10",binmin,binmax);
		pyAsym10->SetTitle("E_{pair} > 10 GeV");
		pyAsym10->SetLineColor(4);
		pyAsym10->Draw();
		
		sprintf(name,"QA_%s_Asymmetry.eps",fCalorimeter.Data());
		c5b->Print(name);
	}
		
	
	// CaloClusters CaloCells
	//printf("c9\n");
	sprintf(cname,"QA_%s_nclustercellsamp",fCalorimeter.Data());
	TCanvas  * c9 = new TCanvas(cname, " CaloClusters and CaloCells", 400, 400) ;
	c9->Divide(2, 2);
	
	c9->cd(1) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNClusters->SetLineColor(4);
	fhNClusters->Draw();
	
	c9->cd(2) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNCells->SetLineColor(4);
	fhNCells->Draw();
	
	c9->cd(3) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNCellsPerCluster->SetLineColor(4);
	fhNCellsPerCluster->Draw();
	
	c9->cd(4) ; 
	gPad->SetLogy();
	//gPad->SetLogx();
	fhAmplitude->SetLineColor(4);
	fhAmplitude->Draw();
	
	sprintf(name,"QA_%s_CaloClustersAndCaloCells.eps",fCalorimeter.Data());
	c9->Print(name);
	//printf("corr\n");
	if(fCorrelateCalos){
		//Calorimeter Correlation, PHOS vs EMCAL
		sprintf(cname,"QA_%s_CaloCorr_EMCALvsPHOS",fCalorimeter.Data());
		TCanvas  * ccorr = new TCanvas(cname, " EMCAL vs PHOS", 400, 400) ;
		ccorr->Divide(2, 2);
		
		ccorr->cd(1) ; 
		//gPad->SetLogy();
		//gPad->SetLogx();
		fhCaloCorrNClusters ->Draw();
		
		ccorr->cd(2) ; 
		//gPad->SetLogy();
		//gPad->SetLogx();
		fhCaloCorrNCells->Draw();
		
		ccorr->cd(3) ; 
		//gPad->SetLogy();
		//gPad->SetLogx();
		fhCaloCorrEClusters->Draw();
		
		ccorr->cd(4) ; 
		//gPad->SetLogy();
		//gPad->SetLogx();
		fhCaloCorrECells->Draw();
		
		sprintf(name,"QA_%s_CaloCorr_EMCALvsPHOS.eps",fCalorimeter.Data());
		ccorr->Print(name); printf("Plot: %s\n",name);
	}
	
	
	if(isDataMC){
	  //Reconstructed vs MC distributions
	  //printf("c6\n");
	  sprintf(cname,"QA_%s_recvsmc",fCalorimeter.Data());
	  TCanvas  * c6 = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
	  c6->Divide(2, 2);
	  
	  fh2E->SetAxisRange(minx,maxx,"X");
	  fh2Pt->SetAxisRange(minx,maxx,"X");
	  fh2E->SetAxisRange(minx,maxx,"Y");
	  fh2Pt->SetAxisRange(minx,maxx,"Y");
	  
	  c6->cd(1) ; 
	  fh2E->SetTitleOffset(1.6,"Y");
	  fh2E->SetLineColor(4);
	  fh2E->Draw();
	  
	  c6->cd(2) ; 
	  fh2Pt->SetTitleOffset(1.6,"Y");
	  fh2Pt->SetLineColor(4);
	  fh2Pt->Draw();
	  
	  c6->cd(3) ; 
	  fh2Phi->SetTitleOffset(1.6,"Y");
	  fh2Phi->SetLineColor(4);
	  fh2Phi->Draw();
	  
	  c6->cd(4) ; 
	  fh2Eta->SetTitleOffset(1.6,"Y");
	  fh2Eta->SetLineColor(4);
	  fh2Eta->Draw();
	  
	  sprintf(name,"QA_%s_ReconstructedVSMCDistributions.eps",fCalorimeter.Data());
	  c6->Print(name);	
	  
	  //Reconstructed vs MC distributions
	  //printf("c6\n");
	  sprintf(cname,"QA_%s_gamrecvsmc",fCalorimeter.Data());
	  TCanvas  * c6Gam = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
	  c6Gam->Divide(2, 2);
	  
	  fhGamE->SetAxisRange(minx,maxx,"X");
	  fhGamPt->SetAxisRange(minx,maxx,"X");
	  
	  c6Gam->cd(1) ; 
	  fhGamE->Draw();
	  
	  c6Gam->cd(2) ; 
	  fhGamPt->Draw();
	  
	  c6Gam->cd(3) ; 
	  fhGamPhi->Draw();
	  
	  c6Gam->cd(4) ; 
	  fhGamEta->Draw();
	  
	  sprintf(name,"QA_%s_GammaReconstructedVSMCDistributions.eps",fCalorimeter.Data());
	  c6->Print(name);	
	  
	  //Generated - reconstructed  
	  //printf("c7\n");
	  sprintf(cname,"QA_%s_diffgenrec",fCalorimeter.Data());
	  TCanvas  * c7 = new TCanvas(cname, "generated - reconstructed", 400, 400) ;
	  c7->Divide(2, 2);
	  
	  c7->cd(1) ; 
	  gPad->SetLogy();
	  fhGamDeltaE->SetLineColor(4);
	  fhDeltaE->Draw();
	  fhGamDeltaE->Draw("same");
	  
	  TLegend pLegendd(0.65,0.55,0.9,0.8);
	  pLegendd.SetTextSize(0.06);
	  pLegendd.AddEntry(fhDeltaE,"all","L");
	  pLegendd.AddEntry(fhGamDeltaE,"from  #gamma","L");
	  pLegendd.SetFillColor(10);
	  pLegendd.SetBorderSize(1);
	  pLegendd.Draw();
	  
	  c7->cd(2) ; 
	  gPad->SetLogy();
	  fhGamDeltaPt->SetLineColor(4);
	  fhDeltaPt->Draw();
	  fhGamDeltaPt->Draw("same");
	  
	  c7->cd(3) ; 
	  gPad->SetLogy();
	  fhGamDeltaPhi->SetLineColor(4);
	  fhDeltaPhi->Draw();
	  fhGamDeltaPhi->Draw("same");
	  
	  c7->cd(4) ; 
	  gPad->SetLogy();
	  fhGamDeltaEta->SetLineColor(4);
	  fhDeltaEta->Draw();
	  fhGamDeltaEta->Draw("same");
	  
	  sprintf(name,"QA_%s_DiffGeneratedReconstructed.eps",fCalorimeter.Data());
	  c7->Print(name);
	  
	  // Reconstructed / Generated 
	  //printf("c8\n");
	  sprintf(cname,"QA_%s_ratiorecgen",fCalorimeter.Data());
	  TCanvas  * c8 = new TCanvas(cname, " reconstructed / generated", 400, 400) ;
	  c8->Divide(2, 2);
	  
	  c8->cd(1) ; 
	  gPad->SetLogy();
	  fhGamRatioE->SetLineColor(4);
	  fhRatioE->Draw();
	  fhGamRatioE->Draw("same");
	  
	  TLegend pLegendr(0.65,0.55,0.9,0.8);
	  pLegendr.SetTextSize(0.06);
	  pLegendr.AddEntry(fhRatioE,"all","L");
	  pLegendr.AddEntry(fhGamRatioE,"from  #gamma","L");
	  pLegendr.SetFillColor(10);
	  pLegendr.SetBorderSize(1);
	  pLegendr.Draw();
	  
	  c8->cd(2) ; 
	  gPad->SetLogy();
	  fhGamRatioPt->SetLineColor(4);
	  fhRatioPt->Draw();
	  fhGamRatioPt->Draw("same");
	  
	  c8->cd(3) ; 
	  fhGamRatioPhi->SetLineColor(4);
	  fhRatioPhi->Draw();
	  fhGamRatioPhi->Draw("same");
	  
	  c8->cd(4) ; 
	  fhGamRatioEta->SetLineColor(4);
	  fhRatioEta->Draw();
	  fhGamRatioEta->Draw("same");
	  
	  sprintf(name,"QA_%s_ReconstructedDivGenerated.eps",fCalorimeter.Data());
	  c8->Print(name);
	  
	  
	  //Generated distributions
	  //printf("c1\n");
	  sprintf(cname,"QA_%s_gen",fCalorimeter.Data());
	  TCanvas  * c10 = new TCanvas(cname, "Generated distributions", 600, 200) ;
	  c10->Divide(3, 1);
	  
	  c10->cd(1) ; 
	  gPad->SetLogy();
	  TH1F * haxispt  = (TH1F*) fhGenPi0Pt->Clone("axispt");  
	  haxispt->SetTitle("Generated Particles p_{T}, |#eta| < 1");
	  fhGenPi0Pt->SetLineColor(1);
	  fhGenGamPt->SetLineColor(4);
	  fhGenEtaPt->SetLineColor(2);
	  fhGenOmegaPt->SetLineColor(7);
	  fhGenElePt->SetLineColor(6);
	  
	  //Select the maximum of the histogram to show all lines.
	  if(fhGenPi0Pt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
	     fhGenPi0Pt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenElePt->GetMaximum())
	    haxispt->SetMaximum(fhGenPi0Pt->GetMaximum());
	  else if(fhGenGamPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
		  fhGenGamPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenElePt->GetMaximum())
	    haxispt->SetMaximum(fhGenGamPt->GetMaximum());
	  else if(fhGenEtaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && 
		  fhGenEtaPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenElePt->GetMaximum())
	    haxispt->SetMaximum(fhGenEtaPt->GetMaximum());	
	  else if(fhGenOmegaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
		  fhGenOmegaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenElePt->GetMaximum())
	    haxispt->SetMaximum(fhGenOmegaPt->GetMaximum());
	  else if(fhGenElePt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
		  fhGenElePt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenGamPt->GetMaximum())
	    haxispt->SetMaximum(fhGenElePt->GetMaximum());
	  haxispt->SetMinimum(1);
	  haxispt->SetAxisRange(minx,maxx,"X");	
	  haxispt->Draw("axis");
	  fhGenPi0Pt->Draw("same");
	  fhGenGamPt->Draw("same");
	  fhGenEtaPt->Draw("same");
	  fhGenOmegaPt->Draw("same");
	  fhGenElePt->Draw("same");
	  
	  TLegend pLegend(0.85,0.65,0.95,0.93);
	  pLegend.SetTextSize(0.06);
	  pLegend.AddEntry(fhGenPi0Pt,"  #pi^{0}","L");
	  pLegend.AddEntry(fhGenGamPt,"  #gamma","L");
	  pLegend.AddEntry(fhGenEtaPt,"  #eta","L");
	  pLegend.AddEntry(fhGenOmegaPt,"  #omega","L");
	  pLegend.AddEntry(fhGenElePt,"  e^{#pm}","L");
	  pLegend.SetFillColor(10);
	  pLegend.SetBorderSize(1);
	  pLegend.Draw();
	  
	  c10->cd(2) ;
	  gPad->SetLogy();
	  TH1F * haxiseta  = (TH1F*) fhGenPi0Eta->Clone("axiseta");  
	  haxiseta->SetTitle("Generated Particles #eta, |#eta| < 1");
	  fhGenPi0Eta->SetLineColor(1);
	  fhGenGamEta->SetLineColor(4);
	  fhGenEtaEta->SetLineColor(2);
	  fhGenOmegaEta->SetLineColor(7);
	  fhGenEleEta->SetLineColor(6);
	  //Select the maximum of the histogram to show all lines.
	  if(fhGenPi0Eta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
	     fhGenPi0Eta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEleEta->GetMaximum())
	    haxiseta->SetMaximum(fhGenPi0Eta->GetMaximum());
	  else if(fhGenGamEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
		  fhGenGamEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEleEta->GetMaximum())
	    haxiseta->SetMaximum(fhGenGamEta->GetMaximum());
	  else if(fhGenEtaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && 
		  fhGenEtaEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
	    haxiseta->SetMaximum(fhGenEtaEta->GetMaximum());	
	  else if(fhGenOmegaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
		  fhGenOmegaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
	    haxiseta->SetMaximum(fhGenOmegaEta->GetMaximum());
	  else if(fhGenEleEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
		  fhGenEleEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenGamEta->GetMaximum())
	    haxiseta->SetMaximum(fhGenEleEta->GetMaximum());
	  haxiseta->SetMinimum(100);
	  haxiseta->Draw("axis");
	  fhGenPi0Eta->Draw("same");
	  fhGenGamEta->Draw("same");
	  fhGenEtaEta->Draw("same");
	  fhGenOmegaEta->Draw("same");
	  fhGenEleEta->Draw("same");
	  
	  
	  c10->cd(3) ; 
	  gPad->SetLogy();
	  TH1F * haxisphi  = (TH1F*) fhGenPi0Phi->Clone("axisphi");  
	  haxisphi->SetTitle("Generated Particles #phi, |#eta| < 1");
	  fhGenPi0Phi->SetLineColor(1);
	  fhGenGamPhi->SetLineColor(4);
	  fhGenEtaPhi->SetLineColor(2);
	  fhGenOmegaPhi->SetLineColor(7);
	  fhGenElePhi->SetLineColor(6);
	  //Select the maximum of the histogram to show all lines.
	  if(fhGenPi0Phi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
	     fhGenPi0Phi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenElePhi->GetMaximum())
	    haxisphi->SetMaximum(fhGenPi0Phi->GetMaximum());
	  else if(fhGenGamPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
		  fhGenGamPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
	    haxisphi->SetMaximum(fhGenGamPhi->GetMaximum());
	  else if(fhGenEtaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && 
		  fhGenEtaPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
	    haxisphi->SetMaximum(fhGenEtaPhi->GetMaximum());	
	  else if(fhGenOmegaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
		  fhGenOmegaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
	    haxisphi->SetMaximum(fhGenOmegaPhi->GetMaximum());
	  else if(fhGenElePhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
		  fhGenElePhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenGamPhi->GetMaximum())
	    haxisphi->SetMaximum(fhGenElePhi->GetMaximum());
	  haxisphi->SetMinimum(100);
	  haxisphi->Draw("axis");
	  fhGenPi0Phi->Draw("same");
	  fhGenGamPhi->Draw("same");
	  fhGenEtaPhi->Draw("same");
	  fhGenOmegaPhi->Draw("same");
	  fhGenElePhi->Draw("same");
	  
	  sprintf(name,"QA_%s_GeneratedDistributions.eps",fCalorimeter.Data());
	  c10->Print(name);
	  
	  
	  //Reconstructed clusters depending on its original particle.
	  //printf("c1\n");
	  sprintf(cname,"QA_%s_recgenid",fCalorimeter.Data());
	  TCanvas  * c11 = new TCanvas(cname, "Reconstructed particles, function of their original particle ID", 400, 400) ;
	  c11->Divide(2, 2);
	  
	  
	  c11->cd(1) ; 
	  gPad->SetLogy();
	  TH1F * hGamE   = (TH1F*) fhGamE->ProjectionX("hGamE",-1,-1);
	  TH1F * hPi0E   = (TH1F*) fhPi0E->ProjectionX("hPi0E",-1,-1);
	  TH1F * hEleE   = (TH1F*) fhEleE->ProjectionX("hEleE",-1,-1);
	  TH1F * hNeHadE = (TH1F*) fhNeHadE->ProjectionX("hNeHadE",-1,-1);
	  TH1F * hChHadE = (TH1F*) fhChHadE->ProjectionX("hChHadE",-1,-1);
	  TH1F * haxisE  = (TH1F*) hPi0E->Clone("axisE");  
	  haxisE->SetTitle("Reconstructed particles E, function of their original particle ID");
	  haxisE->SetAxisRange(minx,maxx,"X");	
	  
	  hPi0E->SetLineColor(1);
	  hGamE->SetLineColor(4);
	  hNeHadE->SetLineColor(2);
	  hChHadE->SetLineColor(7);
	  hEleE->SetLineColor(6);
	  
	  //Select the maximum of the histogram to show all lines.
	  if(hPi0E->GetMaximum() >= hGamE->GetMaximum() && hPi0E->GetMaximum() >= hNeHadE->GetMaximum() && 
	     hPi0E->GetMaximum() >= hChHadE->GetMaximum() && hPi0E->GetMaximum() >= hEleE->GetMaximum())
	    haxisE->SetMaximum(hPi0E->GetMaximum());
	  else if(hGamE->GetMaximum() >= hPi0E->GetMaximum() && hGamE->GetMaximum() >= hNeHadE->GetMaximum() && 
		  hGamE->GetMaximum() >= hChHadE->GetMaximum() && hGamE->GetMaximum() >= hEleE->GetMaximum())
	    haxisE->SetMaximum(hGamE->GetMaximum());
	  else if(hNeHadE->GetMaximum() >= hPi0E->GetMaximum() && hNeHadE->GetMaximum() >= hGamE->GetMaximum() && 
		  hNeHadE->GetMaximum() >= hChHadE->GetMaximum() && hNeHadE->GetMaximum() >= hEleE->GetMaximum())
	    haxisE->SetMaximum(hNeHadE->GetMaximum());	
	  else if(hChHadE->GetMaximum() >= hPi0E->GetMaximum() && hChHadE->GetMaximum() >= hNeHadE->GetMaximum() && 
		  hChHadE->GetMaximum() >= hGamE->GetMaximum() && hChHadE->GetMaximum() >= hEleE->GetMaximum())
	    haxisE->SetMaximum(hChHadE->GetMaximum());
	  else if(hEleE->GetMaximum() >= hPi0E->GetMaximum() && hEleE->GetMaximum() >= hNeHadE->GetMaximum() && 
		  hEleE->GetMaximum() >= hChHadE->GetMaximum() && hEleE->GetMaximum() >= hGamE->GetMaximum())
	    haxisE->SetMaximum(hEleE->GetMaximum());
	  haxisE->SetXTitle("E (GeV)");
	  haxisE->SetMinimum(1);
	  haxisE->Draw("axis");
	  hPi0E->Draw("same");
	  hGamE->Draw("same");
	  hNeHadE->Draw("same");
	  hChHadE->Draw("same");
	  hEleE->Draw("same");
	  
	  TLegend pLegend2(0.8,0.65,0.95,0.93);
	  pLegend2.SetTextSize(0.06);
	  pLegend2.AddEntry(hPi0E,"  #pi^{0}","L");
	  pLegend2.AddEntry(hGamE,"  #gamma","L");
	  pLegend2.AddEntry(hEleE,"  e^{#pm}","L");
	  pLegend2.AddEntry(hChHadE,"  h^{#pm}","L");
	  pLegend2.AddEntry(hNeHadE,"  h^{0}","L");
	  pLegend2.SetFillColor(10);
	  pLegend2.SetBorderSize(1);
	  pLegend2.Draw();
	  
	  
	  c11->cd(2) ; 
	  gPad->SetLogy();
	  //printf("%s, %s, %s, %s, %s\n",fhGamPt->GetName(),fhPi0Pt->GetName(),fhElePt->GetName(),fhNeHadPt->GetName(), fhChHadPt->GetName());
	  TH1F * hGamPt   = (TH1F*) fhGamPt->ProjectionX("hGamPt",-1,-1);
	  TH1F * hPi0Pt   = (TH1F*) fhPi0Pt->ProjectionX("hPi0Pt",-1,-1);
	  TH1F * hElePt   = (TH1F*) fhElePt->ProjectionX("hElePt",-1,-1);
	  TH1F * hNeHadPt = (TH1F*) fhNeHadPt->ProjectionX("hNeHadPt",-1,-1);
	  TH1F * hChHadPt = (TH1F*) fhChHadPt->ProjectionX("hChHadPt",-1,-1);
	  haxispt  = (TH1F*) hPi0Pt->Clone("axispt");  
	  haxispt->SetTitle("Reconstructed particles p_{T}, function of their original particle ID");
	  hPi0Pt->SetLineColor(1);
	  hGamPt->SetLineColor(4);
	  hNeHadPt->SetLineColor(2);
	  hChHadPt->SetLineColor(7);
	  hElePt->SetLineColor(6);
	  
	  //Select the maximum of the histogram to show all lines.
	  if(hPi0Pt->GetMaximum() >= hGamPt->GetMaximum() && hPi0Pt->GetMaximum() >= hNeHadPt->GetMaximum() && 
	     hPi0Pt->GetMaximum() >= hChHadPt->GetMaximum() && hPi0Pt->GetMaximum() >= hElePt->GetMaximum())
	    haxispt->SetMaximum(hPi0Pt->GetMaximum());
	  else if(hGamPt->GetMaximum() >= hPi0Pt->GetMaximum() && hGamPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
		  hGamPt->GetMaximum() >= hChHadPt->GetMaximum() && hGamPt->GetMaximum() >= hElePt->GetMaximum())
	    haxispt->SetMaximum(hGamPt->GetMaximum());
	  else if(hNeHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hNeHadPt->GetMaximum() >= hGamPt->GetMaximum() && 
		  hNeHadPt->GetMaximum() >= hChHadPt->GetMaximum() && hNeHadPt->GetMaximum() >= hElePt->GetMaximum())
	    haxispt->SetMaximum(hNeHadPt->GetMaximum());	
	  else if(hChHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hChHadPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
		  hChHadPt->GetMaximum() >= hGamPt->GetMaximum() && hChHadPt->GetMaximum() >= hElePt->GetMaximum())
	    haxispt->SetMaximum(hChHadPt->GetMaximum());
	  else if(hElePt->GetMaximum() >= hPi0Pt->GetMaximum() && hElePt->GetMaximum() >= hNeHadPt->GetMaximum() && 
		  hElePt->GetMaximum() >= hChHadPt->GetMaximum() && hElePt->GetMaximum() >= hGamPt->GetMaximum())
	    haxispt->SetMaximum(hElePt->GetMaximum());
	  haxispt->SetXTitle("p_{T} (GeV/c)");
	  haxispt->SetMinimum(1);
	  haxispt->SetAxisRange(minx,maxx,"X");
	  haxispt->Draw("axis");
	  hPi0Pt->Draw("same");
	  hGamPt->Draw("same");
	  hNeHadPt->Draw("same");
	  hChHadPt->Draw("same");
	  hElePt->Draw("same");
	  
	  
	  c11->cd(3) ;
	  gPad->SetLogy();
	  
	  TH1F * hGamEta   = (TH1F*) fhGamEta->ProjectionX("hGamEta",-1,-1);
	  TH1F * hPi0Eta   = (TH1F*) fhPi0Eta->ProjectionX("hPi0Eta",-1,-1);
	  TH1F * hEleEta   = (TH1F*) fhEleEta->ProjectionX("hEleEta",-1,-1);
	  TH1F * hNeHadEta = (TH1F*) fhNeHadEta->ProjectionX("hNeHadEta",-1,-1);
	  TH1F * hChHadEta = (TH1F*) fhChHadEta->ProjectionX("hChHadEta",-1,-1);
	  haxiseta  = (TH1F*) hPi0Eta->Clone("axiseta");  
	  haxiseta->SetTitle("Reconstructed particles #eta, function of their original particle ID");
	  hPi0Eta->SetLineColor(1);
	  hGamEta->SetLineColor(4);
	  hNeHadEta->SetLineColor(2);
	  hChHadEta->SetLineColor(7);
	  hEleEta->SetLineColor(6);
	  //Select the maximum of the histogram to show all lines.
	  if(hPi0Eta->GetMaximum() >= hGamEta->GetMaximum() && hPi0Eta->GetMaximum() >= hNeHadEta->GetMaximum() && 
	     hPi0Eta->GetMaximum() >= hChHadEta->GetMaximum() && hPi0Eta->GetMaximum() >= hEleEta->GetMaximum())
	    haxiseta->SetMaximum(hPi0Eta->GetMaximum());
	  else if(hGamEta->GetMaximum() >= hPi0Eta->GetMaximum() && hGamEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
		  hGamEta->GetMaximum() >= hChHadEta->GetMaximum() && hGamEta->GetMaximum() >= hEleEta->GetMaximum())
	    haxiseta->SetMaximum(hGamEta->GetMaximum());
	  else if(hNeHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hNeHadEta->GetMaximum() >= hGamEta->GetMaximum() && 
		  hNeHadEta->GetMaximum() >= hChHadEta->GetMaximum() && hNeHadEta->GetMaximum() >= hEleEta->GetMaximum())
	    haxiseta->SetMaximum(hNeHadEta->GetMaximum());	
	  else if(hChHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hChHadEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
		  hChHadEta->GetMaximum() >= hGamEta->GetMaximum() && hChHadEta->GetMaximum() >= hEleEta->GetMaximum())
	    haxiseta->SetMaximum(hChHadEta->GetMaximum());
	  else if(hEleEta->GetMaximum() >= hPi0Eta->GetMaximum() && hEleEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
		  hEleEta->GetMaximum() >= hChHadEta->GetMaximum() && hEleEta->GetMaximum() >= hGamEta->GetMaximum())
	    haxiseta->SetMaximum(hEleEta->GetMaximum());
	  
	  haxiseta->SetXTitle("#eta");
	  haxiseta->Draw("axis");
	  hPi0Eta->Draw("same");
	  hGamEta->Draw("same");
	  hNeHadEta->Draw("same");
	  hChHadEta->Draw("same");
	  hEleEta->Draw("same");
	  
	  
	  c11->cd(4) ; 
	  gPad->SetLogy();
	  TH1F * hGamPhi   = (TH1F*) fhGamPhi->ProjectionX("hGamPhi",-1,-1);
	  TH1F * hPi0Phi   = (TH1F*) fhPi0Phi->ProjectionX("hPi0Phi",-1,-1);
	  TH1F * hElePhi   = (TH1F*) fhElePhi->ProjectionX("hElePhi",-1,-1);
	  TH1F * hNeHadPhi = (TH1F*) fhNeHadPhi->ProjectionX("hNeHadPhi",-1,-1);
	  TH1F * hChHadPhi = (TH1F*) fhChHadPhi->ProjectionX("hChHadPhi",-1,-1);
	  haxisphi  = (TH1F*) hPi0Phi->Clone("axisphi");  
	  haxisphi->SetTitle("Reconstructed particles #phi, function of their original particle ID");
	  
	  hPi0Phi->SetLineColor(1);
	  hGamPhi->SetLineColor(4);
	  hNeHadPhi->SetLineColor(2);
	  hChHadPhi->SetLineColor(7);
	  hElePhi->SetLineColor(6);
	  //Select the maximum of the histogram to show all lines.
	  if(hPi0Phi->GetMaximum() >= hGamPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
	     hPi0Phi->GetMaximum() >= hChHadPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hElePhi->GetMaximum())
	    haxisphi->SetMaximum(hPi0Phi->GetMaximum());
	  else if(hGamPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hGamPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
		  hGamPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hGamPhi->GetMaximum() >= hElePhi->GetMaximum())
	    haxisphi->SetMaximum(hGamPhi->GetMaximum());
	  else if(hNeHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hNeHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && 
		  hNeHadPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hNeHadPhi->GetMaximum() >= hElePhi->GetMaximum())
	    haxisphi->SetMaximum(hNeHadPhi->GetMaximum());	
	  else if(hChHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hChHadPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
		  hChHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && hChHadPhi->GetMaximum() >= hElePhi->GetMaximum())
	    haxisphi->SetMaximum(hChHadPhi->GetMaximum());
	  else if(hElePhi->GetMaximum() >= hPi0Phi->GetMaximum() && hElePhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
		  hElePhi->GetMaximum() >= hChHadPhi->GetMaximum() && hElePhi->GetMaximum() >= hGamPhi->GetMaximum())
	    haxisphi->SetMaximum(hElePhi->GetMaximum());
	  haxisphi->SetXTitle("#phi (rad)");
	  haxisphi->Draw("axis");
	  hPi0Phi->Draw("same");
	  hGamPhi->Draw("same");
	  hNeHadPhi->Draw("same");
	  hChHadPhi->Draw("same");
	  hElePhi->Draw("same");
	  
	  sprintf(name,"QA_%s_RecDistributionsGenID.eps",fCalorimeter.Data());
	  c11->Print(name);
	  
	  
	  //Ratio reconstructed clusters / generated particles in acceptance, for different particle ID
	  //printf("c1\n");
	  
	  TH1F *	hPi0EClone   = (TH1F*)   hPi0E  ->Clone("hPi0EClone");
	  TH1F *	hGamEClone   = (TH1F*)   hGamE  ->Clone("hGamEClone");
	  TH1F *	hPi0PtClone  = (TH1F*)   hPi0Pt ->Clone("hPi0PtClone");
	  TH1F *	hGamPtClone  = (TH1F*)   hGamPt ->Clone("hGamPtClone");	
	  TH1F *	hPi0EtaClone = (TH1F*)   hPi0Eta->Clone("hPi0EtaClone");
	  TH1F *	hGamEtaClone = (TH1F*)   hGamEta->Clone("hGamEtaClone");	
	  TH1F *	hPi0PhiClone = (TH1F*)   hPi0Phi->Clone("hPi0PhiClone");
	  TH1F *	hGamPhiClone = (TH1F*)   hGamPhi->Clone("hGamPhiClone");	
	  
	  sprintf(cname,"QA_%s_recgenidratio",fCalorimeter.Data());
	  TCanvas  * c12 = new TCanvas(cname, "Ratio reconstructed clusters / generated particles in acceptance, for different particle ID", 400, 400) ;
	  c12->Divide(2, 2);
	  
	  c12->cd(1) ; 
	  gPad->SetLogy();
	  haxisE->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	  hPi0EClone->Divide(fhGenPi0AccE);
	  hGamEClone->Divide(fhGenGamAccE);
	  haxisE->SetMaximum(1000);
	  haxisE->SetMinimum(1e-2);
	  haxisE->SetXTitle("E (GeV)");
	  haxisE->SetYTitle("ratio = rec/gen");
	  haxisE->Draw("axis");
	  hPi0E->Draw("same");
	  hGamE->Draw("same");
	
	  TLegend pLegend3(0.75,0.2,0.9,0.4);
	  pLegend3.SetTextSize(0.06);
	  pLegend3.AddEntry(hPi0EClone,"  #pi^{0}","L");
	  pLegend3.AddEntry(hGamEClone,"  #gamma","L");
	  pLegend3.SetFillColor(10);
	  pLegend3.SetBorderSize(1);
	  pLegend3.Draw();
	  
	  c12->cd(2) ; 
	  gPad->SetLogy();
	  haxispt->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	  hPi0PtClone->Divide(fhGenPi0AccPt);
	  hGamPtClone->Divide(fhGenGamAccPt);
	  haxispt->SetMaximum(5);
	  haxispt->SetMinimum(1e-2);
	  haxispt->SetXTitle("p_{T} (GeV/c)");
	  haxispt->SetYTitle("ratio = rec/gen");
	  haxispt->Draw("axis");
	  hPi0PtClone->Draw("same");
	  hGamPtClone->Draw("same");
	  
	  c12->cd(3) ;
	  gPad->SetLogy();
	  
	  haxiseta->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	  hPi0EtaClone->Divide(fhGenPi0AccEta);
	  hGamEtaClone->Divide(fhGenGamAccEta);
	  haxiseta->SetMaximum(1.2);
	  haxiseta->SetMinimum(1e-2);
	  haxiseta->SetYTitle("ratio = rec/gen");
	  haxiseta->SetXTitle("#eta");
	  haxiseta->Draw("axis");
	  hPi0EtaClone->Draw("same");
	  hGamEtaClone->Draw("same");
	  
	  
	  c12->cd(4) ; 
	  gPad->SetLogy();
	  haxisphi->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	  hPi0PhiClone->Divide(fhGenPi0AccPhi);
	  hGamPhiClone->Divide(fhGenGamAccPhi);
	  haxisphi->SetYTitle("ratio = rec/gen");
	  haxisphi->SetXTitle("#phi (rad)");
	  haxisphi->SetMaximum(1.2);
	  haxisphi->SetMinimum(1e-2);
	  haxisphi->Draw("axis");
	  hPi0PhiClone->Draw("same");
	  hGamPhiClone->Draw("same");
	  
	  sprintf(name,"QA_%s_EfficiencyGenID.eps",fCalorimeter.Data());
	  c12->Print(name);
	  
	  
	  //Reconstructed distributions
	  //printf("c1\n");
	  sprintf(cname,"QA_%s_vertex",fCalorimeter.Data());
	  TCanvas  * c13 = new TCanvas(cname, "Particle vertex", 400, 400) ;
	  c13->Divide(2, 2);
	  
	  c13->cd(1) ; 
	  //gPad->SetLogy();
	  fhEMVxyz->SetTitleOffset(1.6,"Y");
	  fhEMVxyz->Draw();
	  
	  c13->cd(2) ; 
	  //gPad->SetLogy();
	  fhHaVxyz->SetTitleOffset(1.6,"Y");
	  fhHaVxyz->Draw();
	  
	  c13->cd(3) ;
	  gPad->SetLogy();
	  TH1F * hEMR = (TH1F*) fhEMR->ProjectionY("hEM",-1,-1); 
	  hEMR->SetLineColor(4);
	  hEMR->Draw();
	  
	  c13->cd(4) ; 
	  gPad->SetLogy();
	  TH1F * hHaR = (TH1F*) fhHaR->ProjectionY("hHa",-1,-1); 
	  hHaR->SetLineColor(4);
	  hHaR->Draw();
	  
	  
	  sprintf(name,"QA_%s_ParticleVertex.eps",fCalorimeter.Data());
	  c13->Print(name);
	  
	  
	  //Track-matching distributions
	  
	  
	  //Reconstructed distributions, matched with tracks, generated particle dependence
	  //printf("c2\n");
	  sprintf(cname,"QA_%s_rectrackmatchGenID",fCalorimeter.Data());
	  TCanvas  * c22ch = new TCanvas(cname, "Reconstructed distributions, matched with tracks, for different particle ID", 400, 400) ;
	  c22ch->Divide(2, 2);
	  
	  c22ch->cd(1) ; 
	  
	  TH1F * hGamECharged   = (TH1F*) fhGamECharged->ProjectionX("hGamECharged",-1,-1);
	  TH1F * hPi0ECharged   = (TH1F*) fhPi0ECharged->ProjectionX("hPi0ECharged",-1,-1);
	  TH1F * hEleECharged   = (TH1F*) fhEleECharged->ProjectionX("hEleECharged",-1,-1);
	  TH1F * hNeHadECharged = (TH1F*) fhNeHadECharged->ProjectionX("hNeHadECharged",-1,-1);
	  TH1F * hChHadECharged = (TH1F*) fhChHadECharged->ProjectionX("hChHadECharged",-1,-1);
	  hPi0ECharged->SetLineColor(1);
	  hGamECharged->SetLineColor(4);
	  hNeHadECharged->SetLineColor(2);
	  hChHadECharged->SetLineColor(7);
	  hEleECharged->SetLineColor(6);	
	  gPad->SetLogy();
	  fhECharged->SetLineColor(3);
	  fhECharged->SetMinimum(0.5);
	  fhECharged->Draw();
	  hPi0ECharged->Draw("same");
	  hGamECharged->Draw("same");
	  hNeHadECharged->Draw("same");
	  hChHadECharged->Draw("same");
	  hEleECharged->Draw("same");
	  TLegend pLegend22(0.75,0.45,0.9,0.8);
	  pLegend22.SetTextSize(0.06);
	  pLegend22.AddEntry(fhECharged,"all","L");
	  pLegend22.AddEntry(hPi0ECharged,"#pi^{0}","L");
	  pLegend22.AddEntry(hGamECharged,"#gamma","L");
	  pLegend22.AddEntry(hEleECharged,"e^{#pm}","L");
	  pLegend22.AddEntry(hChHadECharged,"h^{#pm}","L");
	  pLegend22.AddEntry(hNeHadECharged,"h^{0}","L");
	  pLegend22.SetFillColor(10);
	  pLegend22.SetBorderSize(1);
	  pLegend22.Draw();
	  
	  c22ch->cd(2) ; 
	  
	  TH1F * hGamPtCharged   = (TH1F*) fhGamPtCharged->ProjectionX("hGamPtCharged",-1,-1);
	  TH1F * hPi0PtCharged   = (TH1F*) fhPi0PtCharged->ProjectionX("hPi0PtCharged",-1,-1);
	  TH1F * hElePtCharged   = (TH1F*) fhElePtCharged->ProjectionX("hElePtCharged",-1,-1);
	  TH1F * hNeHadPtCharged = (TH1F*) fhNeHadPtCharged->ProjectionX("hNeHadPtCharged",-1,-1);
	  TH1F * hChHadPtCharged = (TH1F*) fhChHadPtCharged->ProjectionX("hChHadPtCharged",-1,-1);
	  hPi0PtCharged->SetLineColor(1);
	  hGamPtCharged->SetLineColor(4);
	  hNeHadPtCharged->SetLineColor(2);
	  hChHadPtCharged->SetLineColor(7);
	  hElePtCharged->SetLineColor(6);	
	  gPad->SetLogy();
	  fhPtCharged->SetLineColor(3);
	  fhPtCharged->SetMinimum(0.5);
	  fhPtCharged->Draw();
	  hPi0PtCharged->Draw("same");
	  hGamPtCharged->Draw("same");
	  hNeHadPtCharged->Draw("same");
	  hChHadPtCharged->Draw("same");
	  hElePtCharged->Draw("same");	
	  
	  c22ch->cd(4) ; 
	  
	  TH1F * hGamEtaCharged   = (TH1F*) fhGamEtaCharged->ProjectionX("hGamEtaCharged",-1,-1);
	  TH1F * hPi0EtaCharged   = (TH1F*) fhPi0EtaCharged->ProjectionX("hPi0EtaCharged",-1,-1);
	  TH1F * hEleEtaCharged   = (TH1F*) fhEleEtaCharged->ProjectionX("hEleEtaCharged",-1,-1);
	  TH1F * hNeHadEtaCharged = (TH1F*) fhNeHadEtaCharged->ProjectionX("hNeHadEtaCharged",-1,-1);
	  TH1F * hChHadEtaCharged = (TH1F*) fhChHadEtaCharged->ProjectionX("hChHadEtaCharged",-1,-1);
	  hPi0EtaCharged->SetLineColor(1);
	  hGamEtaCharged->SetLineColor(4);
	  hNeHadEtaCharged->SetLineColor(2);
	  hChHadEtaCharged->SetLineColor(7);
	  hEleEtaCharged->SetLineColor(6);	
	  gPad->SetLogy();
	  fhEtaCharged->SetLineColor(3);
	  fhEtaCharged->SetMinimum(0.5);
	  fhEtaCharged->Draw();
	  hPi0EtaCharged->Draw("same");
	  hGamEtaCharged->Draw("same");
	  hNeHadEtaCharged->Draw("same");
	  hChHadEtaCharged->Draw("same");
	  hEleEtaCharged->Draw("same");
	  
	  c22ch->cd(3) ; 
	  
	  TH1F * hGamPhiCharged   = (TH1F*) fhGamPhiCharged->ProjectionX("hGamPhiCharged",-1,-1);
	  TH1F * hPi0PhiCharged   = (TH1F*) fhPi0PhiCharged->ProjectionX("hPi0PhiCharged",-1,-1);
	  TH1F * hElePhiCharged   = (TH1F*) fhElePhiCharged->ProjectionX("hElePhiCharged",-1,-1);
	  TH1F * hNeHadPhiCharged = (TH1F*) fhNeHadPhiCharged->ProjectionX("hNeHadPhiCharged",-1,-1);
	  TH1F * hChHadPhiCharged = (TH1F*) fhChHadPhiCharged->ProjectionX("hChHadPhiCharged",-1,-1);
	  hPi0PhiCharged->SetLineColor(1);
	  hGamPhiCharged->SetLineColor(4);
	  hNeHadPhiCharged->SetLineColor(2);
	  hChHadPhiCharged->SetLineColor(7);
	  hElePhiCharged->SetLineColor(6);	
	  gPad->SetLogy();
	  fhPhiCharged->SetLineColor(3);
	  fhPhiCharged->SetMinimum(0.5);
	  fhPhiCharged->Draw();
	  hPi0PhiCharged->Draw("same");
	  hGamPhiCharged->Draw("same");
	  hNeHadPhiCharged->Draw("same");
	  hChHadPhiCharged->Draw("same");
	  hElePhiCharged->Draw("same");
	  
	  
	  sprintf(name,"QA_%s_ReconstructedDistributions_TrackMatchedGenID.eps",fCalorimeter.Data());
	  c22ch->Print(name);
	  
	  TH1F *	hGamEChargedClone   = (TH1F*)   hGamECharged->Clone("GamEChargedClone");
	  TH1F *	hGamPtChargedClone  = (TH1F*)   hGamPtCharged->Clone("GamPtChargedClone");
	  TH1F *	hGamEtaChargedClone = (TH1F*)   hGamEtaCharged->Clone("GamEtaChargedClone");
	  TH1F *	hGamPhiChargedClone = (TH1F*)   hGamPhiCharged->Clone("GamPhiChargedClone");
	  
	  TH1F *	hPi0EChargedClone   = (TH1F*)   hPi0ECharged->Clone("Pi0EChargedClone");
	  TH1F *	hPi0PtChargedClone  = (TH1F*)   hPi0PtCharged->Clone("Pi0PtChargedClone");
	  TH1F *	hPi0EtaChargedClone = (TH1F*)   hPi0EtaCharged->Clone("Pi0EtaChargedClone");
	  TH1F *	hPi0PhiChargedClone = (TH1F*)   hPi0PhiCharged->Clone("Pi0PhiChargedClone");
	
	  TH1F *	hEleEChargedClone   = (TH1F*)   hEleECharged->Clone("EleEChargedClone");
	  TH1F *	hElePtChargedClone  = (TH1F*)   hElePtCharged->Clone("ElePtChargedClone");
	  TH1F *	hEleEtaChargedClone = (TH1F*)   hEleEtaCharged->Clone("EleEtaChargedClone");
	  TH1F *	hElePhiChargedClone = (TH1F*)   hElePhiCharged->Clone("ElePhiChargedClone");	
	  
	  TH1F *	hNeHadEChargedClone   = (TH1F*)   hNeHadECharged->Clone("NeHadEChargedClone");
	  TH1F *	hNeHadPtChargedClone  = (TH1F*)   hNeHadPtCharged->Clone("NeHadPtChargedClone");
	  TH1F *	hNeHadEtaChargedClone = (TH1F*)   hNeHadEtaCharged->Clone("NeHadEtaChargedClone");
	  TH1F *	hNeHadPhiChargedClone = (TH1F*)   hNeHadPhiCharged->Clone("NeHadPhiChargedClone");
	  
	  TH1F *	hChHadEChargedClone   = (TH1F*)   hChHadECharged->Clone("ChHadEChargedClone");
	  TH1F *	hChHadPtChargedClone  = (TH1F*)   hChHadPtCharged->Clone("ChHadPtChargedClone");
	  TH1F *	hChHadEtaChargedClone = (TH1F*)   hChHadEtaCharged->Clone("ChHadEtaChargedClone");
	  TH1F *	hChHadPhiChargedClone = (TH1F*)   hChHadPhiCharged->Clone("ChHadPhiChargedClone");	
	  
	  //Ratio: reconstructed track matched/ all reconstructed
	  //printf("c3\n");
	  sprintf(cname,"QA_%s_rectrackmatchratGenID",fCalorimeter.Data());
	  TCanvas  * c3ch = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed, for different particle ID", 400, 400) ;
	  c3ch->Divide(2, 2);
	  
	  c3ch->cd(1) ;
	  hEChargedClone->SetMaximum(1.2);
	  hEChargedClone->SetMinimum(0.001);	
	  hEChargedClone->SetLineColor(3);
	  hEChargedClone->SetYTitle("track matched / all");
	  hPi0EChargedClone->Divide(hPi0E);
	  hGamEChargedClone->Divide(hGamE);
	  hEleEChargedClone->Divide(hEleE);
	  hNeHadEChargedClone->Divide(hNeHadE);
	  hChHadEChargedClone->Divide(hChHadE);
	  hEChargedClone->Draw();
	  hPi0EChargedClone->Draw("same");
	  hGamEChargedClone->Draw("same");
	  hEleEChargedClone->Draw("same");
	  hNeHadEChargedClone->Draw("same");
	  hChHadEChargedClone->Draw("same");
	  
	  TLegend pLegend3ch(0.75,0.45,0.9,0.8);
	  pLegend3ch.SetTextSize(0.06);
	  pLegend3ch.AddEntry(hEChargedClone,"all","L");
	  pLegend3ch.AddEntry(hPi0EChargedClone,"#pi^{0}","L");
	  pLegend3ch.AddEntry(hGamEChargedClone,"#gamma","L");
	  pLegend3ch.AddEntry(hEleEChargedClone,"e^{#pm}","L");
	  pLegend3ch.AddEntry(hChHadEChargedClone,"h^{#pm}","L");
	  pLegend3ch.AddEntry(hNeHadEChargedClone,"h^{0}","L");
	  pLegend3ch.SetFillColor(10);
	  pLegend3ch.SetBorderSize(1);
	  pLegend3ch.Draw();
	  
	  c3ch->cd(2) ;
	  hPtChargedClone->SetMaximum(1.2);
	  hPtChargedClone->SetMinimum(0.001);	
	  hPtChargedClone->SetLineColor(3);
	  hPtChargedClone->SetYTitle("track matched / all");
	  hPi0PtChargedClone->Divide(hPi0Pt);
	  hGamPtChargedClone->Divide(hGamPt);
	  hElePtChargedClone->Divide(hElePt);
	  hNeHadPtChargedClone->Divide(hNeHadPt);
	  hChHadPtChargedClone->Divide(hChHadPt);
	  hPtChargedClone->Draw();
	  hPi0PtChargedClone->Draw("same");
	  hGamPtChargedClone->Draw("same");
	  hElePtChargedClone->Draw("same");
	  hNeHadPtChargedClone->Draw("same");
	  hChHadPtChargedClone->Draw("same");
	  
	  c3ch->cd(4) ;
	  hEtaChargedClone->SetMaximum(1.2);
	  hEtaChargedClone->SetMinimum(0.001);	
	  hEtaChargedClone->SetLineColor(3);
	  hEtaChargedClone->SetYTitle("track matched / all");
	  hPi0EtaChargedClone->Divide(hPi0Eta);
	  hGamEtaChargedClone->Divide(hGamEta);
	  hEleEtaChargedClone->Divide(hEleEta);
	  hNeHadEtaChargedClone->Divide(hNeHadEta);
	  hChHadEtaChargedClone->Divide(hChHadEta);
	  hEtaChargedClone->Draw();
	  hPi0EtaChargedClone->Draw("same");
	  hGamEtaChargedClone->Draw("same");
	  hEleEtaChargedClone->Draw("same");
	  hNeHadEtaChargedClone->Draw("same");
	  hChHadEtaChargedClone->Draw("same");
	  
	  c3ch->cd(3) ;
	  hPhiChargedClone->SetMaximum(1.2);
	  hPhiChargedClone->SetMinimum(0.001);
	  hPhiChargedClone->SetLineColor(3);
	  hPhiChargedClone->SetYTitle("track matched / all");
	  hPi0PhiChargedClone->Divide(hPi0Phi);
	  hGamPhiChargedClone->Divide(hGamPhi);
	  hElePhiChargedClone->Divide(hElePhi);
	  hNeHadPhiChargedClone->Divide(hNeHadPhi);
	  hChHadPhiChargedClone->Divide(hChHadPhi);
	  hPhiChargedClone->Draw();
	  hPi0PhiChargedClone->Draw("same");
	  hGamPhiChargedClone->Draw("same");
	  hElePhiChargedClone->Draw("same");
	  hNeHadPhiChargedClone->Draw("same");
	  hChHadPhiChargedClone->Draw("same");
	  
	  sprintf(name,"QA_%s_RatioReconstructedMatchedDistributionsGenID.eps",fCalorimeter.Data());
	  c3ch->Print(name);	  
	  
	  sprintf(cname,"QA_%s_trkmatchMCEle",fCalorimeter.Data());
	  TCanvas *cmemc = new TCanvas(cname,"Track-matching distributions from MC electrons", 600, 200);
	  cmemc->Divide(3,1);
	  
	  cmemc->cd(1);
	  gPad->SetLogy();
	  fhMCEle1pOverE->Draw();
	  fhMCEle1pOverER02->SetLineColor(4);
	  fhMCEle1pOverE->SetLineColor(1);
	  fhMCEle1pOverER02->Draw("same");
	  pLegendpE0.Draw();
	  
	  cmemc->cd(2);
	  gPad->SetLogy();
	  fhMCEle1dR->Draw();
	  
	  cmemc->cd(3);
	  fhMCEle2MatchdEdx->Draw();
	  
	  sprintf(name,"QA_%s_TrackMatchingDistMCEle.eps",fCalorimeter.Data());
	  cmemc->Print(name);  
	  
	  sprintf(cname,"QA_%s_trkmatchMCChHad",fCalorimeter.Data());
	  TCanvas *cmemchad = new TCanvas(cname,"Track-matching distributions from MC charged hadrons", 600, 200);
	  cmemchad->Divide(3,1);
	  
	  cmemchad->cd(1);
	  gPad->SetLogy();
	  fhMCChHad1pOverE->Draw();
	  fhMCChHad1pOverE->SetLineColor(1);
	  fhMCChHad1pOverER02->SetLineColor(4);
	  fhMCChHad1pOverER02->Draw("same");
	  pLegendpE0.Draw();
	  
	  cmemchad->cd(2);
	  gPad->SetLogy();
	  fhMCChHad1dR->Draw();
	  
	  cmemchad->cd(3);
	  fhMCChHad2MatchdEdx->Draw();
	  
	  sprintf(name,"QA_%s_TrackMatchingDistMCChHad.eps",fCalorimeter.Data());
	  cmemchad->Print(name);       
	  
	  sprintf(cname,"QA_%s_trkmatchMCNeutral",fCalorimeter.Data());
	  TCanvas *cmemcn = new TCanvas(cname,"Track-matching distributions from MC neutrals", 600, 200);
	  cmemcn->Divide(3,1);
	  
	  cmemcn->cd(1);
	  gPad->SetLogy();
	  fhMCNeutral1pOverE->Draw();
	  fhMCNeutral1pOverE->SetLineColor(1);
	  fhMCNeutral1pOverER02->SetLineColor(4);
	  fhMCNeutral1pOverER02->Draw("same");
	  pLegendpE0.Draw();
	  
	  cmemcn->cd(2);
	  gPad->SetLogy();
	  fhMCNeutral1dR->Draw();
	  
	  cmemcn->cd(3);
	  fhMCNeutral2MatchdEdx->Draw();
	  
	  sprintf(name,"QA_%s_TrackMatchingDistMCNeutral.eps",fCalorimeter.Data());
	  cmemcn->Print(name);       
	  
	  sprintf(cname,"QA_%s_trkmatchpE",fCalorimeter.Data());
	  TCanvas *cmpoe = new TCanvas(cname,"Track-matching distributions, p/E", 400, 200);
	  cmpoe->Divide(2,1);
	  
	  cmpoe->cd(1);
	  gPad->SetLogy();
	  fh1pOverE->SetLineColor(1);
	  fhMCEle1pOverE->SetLineColor(4);
	  fhMCChHad1pOverE->SetLineColor(2);
	  fhMCNeutral1pOverE->SetLineColor(7);
	  fh1pOverE->SetMinimum(0.5);
	  fh1pOverE->Draw();
	  fhMCEle1pOverE->Draw("same");
	  fhMCChHad1pOverE->Draw("same");
	  fhMCNeutral1pOverE->Draw("same");
	  TLegend pLegendpE(0.65,0.55,0.9,0.8);
	  pLegendpE.SetTextSize(0.06);
	  pLegendpE.AddEntry(fh1pOverE,"all","L");
	  pLegendpE.AddEntry(fhMCEle1pOverE,"e^{#pm}","L");
	  pLegendpE.AddEntry(fhMCChHad1pOverE,"h^{#pm}","L");
	  pLegendpE.AddEntry(fhMCNeutral1pOverE,"neutrals","L");
	  pLegendpE.SetFillColor(10);
	  pLegendpE.SetBorderSize(1);
	  pLegendpE.Draw();
	  
	  cmpoe->cd(2);
	  gPad->SetLogy();
	  fh1pOverER02->SetTitle("Track matches p/E, dR<0.2");
	  fh1pOverER02->SetLineColor(1);
	  fhMCEle1pOverER02->SetLineColor(4);
	  fhMCChHad1pOverER02->SetLineColor(2);
	  fhMCNeutral1pOverER02->SetLineColor(7);
	  fh1pOverER02->SetMaximum(fh1pOverE->GetMaximum());
	  fh1pOverER02->SetMinimum(0.5);
	  fh1pOverER02->Draw();
	  fhMCEle1pOverER02->Draw("same");
	  fhMCChHad1pOverER02->Draw("same");
	  fhMCNeutral1pOverER02->Draw("same");
	  
	  //		TLegend pLegendpE2(0.65,0.55,0.9,0.8);
	  //		pLegendpE2.SetTextSize(0.06);
	  //		pLegendpE2.SetHeader("dR < 0.02");
	  //		pLegendpE2.SetFillColor(10);
	  //		pLegendpE2.SetBorderSize(1);
	  //		pLegendpE2.Draw();
	//Track matching histograms
	sprintf(name,"QA_%s_TrackMatchingPOverE.eps",fCalorimeter.Data());
	cmpoe->Print(name);       	
	}//MC related histograms

	//	char line[1024] ; 
	//	sprintf(line, ".!tar -zcf QA_%s_%s.tar.gz *%s*.eps", fCalorimeter.Data(), GetName(),fCalorimeter.Data()) ; 
	//	gROOT->ProcessLine(line);
	//	sprintf(line, ".!rm -fR *.eps"); 
	//	gROOT->ProcessLine(line);
	//	
	//	printf("AliAnaCalorimeterQA::Terminate() - !! All the pdf files are in QA_%s_%s.tar.gz !!!\n",  fCalorimeter.Data(), GetName());
	
}

