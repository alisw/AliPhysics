void V0_efficiencies()

{      
// Caution: time windows and cuts on signals are HARD-coded 
//          they need to be adjusted to the current configuration

        gROOT->Reset();
	rl = AliRunLoader::Open("galice.root");
	rl->LoadgAlice();
	gAlice = rl->GetAliRun();
        Int_t nevent = rl->GetNumberOfEvents();
	cout<<" --- Number of events in file = "<< nevent <<" ---"<<endl;
	
//___________________________________________________
// Book HISTOGRAMS 
      
    	TH2F *hMapV0L   = new TH2F("hMapV0L","V0L",161,-0.8,0.8,161,-0.8,0.8);
    	TH2F *hMapV0R   = new TH2F("hMapV0R","V0R",161,-0.8,0.8,161,-0.8,0.8);	 
	TH1F *hpdgV0    = new TH1F("hpdgV0","pdgV0",1001,-500,500);  
	TH1F *hvertexZ  = new TH1F("hvertexZ","vertex Z",972,-40.,40.); 
	TH1F *hTimC     = new TH1F("hTimC","Time of Flight in V0C",500,0,49);
	TH1F *hTimA     = new TH1F("hTimA","Time of Flight in V0A",500,0,49);
	TH1F *hCell     = new TH1F("hCell","Cell Number",100,0,99);
	TH2F *hCellADC  = new TH2F("hCellADC","ADC vs Cell",100,0,99,1000,0,3999);
 
        TH1F *hMul0     = new TH1F("hMul0 ","Multiplicity in V0",80,0,79);
	TH1F *hMulC0    = new TH1F("hMulC0","Multiplicity in V0C",50,0,49);
        TH1F *hMulA0    = new TH1F("hMulA0","Multiplicity in V0A",50,0,49);
        TH1F *hMulAnd0  = new TH1F("hMulAnd0","Trigger and",50,0,49);
	
//___________________________________________________
	
	AliLoader *ld = rl->GetLoader("VZEROLoader");
	ld->GetHitsDataLoader()->Load("READ");
	
        rl->LoadHeader();
  
        // to access the particle  Stack
        rl->LoadKinematics("READ");

	AliVZERO  *v0 = (AliVZERO*)gAlice->GetDetector("VZERO");
	
	bool    flagV0L  = false;
	bool    flagV0R  = false;	
	Float_t timeV0L  = 1e12;
	Float_t timeV0R  = 1e12;	
	
	Int_t   fNdigits = 0;
	Int_t   fMulA    = 0;
	Int_t   fMulC    = 0;
        Int_t   fDigits  = 0;  
        Float_t fPhotoCathodeEfficiency =   0.18;
	
	
	Int_t   nVOL     = 0;
	Int_t   nVOR     = 0;
	Int_t   nVOLetR  = 0;
	Int_t   nVOLouR  = 0;
	
	TDatabasePDG * pdgdb = TDatabasePDG::Instance();


        Float_t fPMVoltage =  768.0;
        Float_t fPMGain1   = TMath::Power((fPMVoltage / 112.5) ,7.04277);
	Float_t fPMGain[64];
	Float_t cPM[64]; 
	
	for(Int_t ii=0; ii<64; ii++){
	fPMGain[ii] = gRandom->Gaus(fPMGain1, fPMGain1/5); // 20% uncertainty on the PM gain
        cPM[ii] = fPhotoCathodeEfficiency * fPMGain[ii];}
	
	Float_t kC     =  1e-11;
        Float_t kthau  =  2.1*1e-9;
        Float_t ktheta =  50.0 * kC;
        Float_t kQe    =  1.6e-19;
	Float_t coef   =  200.0;       //  p-p
//	Float_t coef   =    1.0;       // Pb-Pb
		
        Int_t map[80];     
        for(Int_t i=0; i<80; i++) map[i] = 0;
        Int_t hit[80];     
        for(Int_t i=0; i<80; i++) hit[i] = 0;
		
	for(Int_t event = 0; event < nevent; event++){
	Float_t timeV0L = 1e12;
	Float_t timeV0R = 1e12;	
	
	    rl->GetEvent(event);
            for(Int_t i=0; i<80; i++) map[i] = 0;
            for(Int_t i=0; i<80; i++) hit[i] = 0;
	    ld->LoadHits();
	    v0->SetTreeAddress();
	    TTree *treeH = ld->TreeH();
	    Int_t ntracks = treeH->GetEntries();
	
	    for (Int_t itrack=0; itrack<ntracks;itrack++) {
    
		v0->ResetHits();
		treeH->GetEvent(itrack);
		if(v0){
		  for (AliVZEROhit *vhit=(AliVZEROhit*)v0->FirstHit(itrack);
			  vhit; vhit=(AliVZEROhit*)v0->NextHit())
		  {
//			vhit->Dump();
                        hpdgV0->Fill(vhit->TrackPiD());
			hvertexZ->Fill(vhit->Vz()/100.);
			
			Float_t dt_scintillator = gRandom->Gaus(0,0.7);		
		        Float_t time = dt_scintillator + 1e9*vhit->Tof();
			
                        if(vhit->Z() > 0){
			if (time < 9 || time > 13) {continue;}		//time window V0A : 11 +/- 2 ns
			
			flagV0L = true;	
			hTimA->Fill(time);	
			if(time < timeV0L) timeV0L = time;	
		        hMapV0L->Fill(vhit->X()/100.,vhit->Y()/100.);}
					 
                        if(vhit->Z() < 0){
			if (time < 1 ||time > 5) {continue;}		//time window V0C : 3 +/- 2 ns
			
			flagV0R = true;	
			hTimC->Fill(time);
			if(time < timeV0R) timeV0R = time;
		        hMapV0R->Fill(vhit->X()/100.,vhit->Y()/100.);}
			
			Int_t nPhot = vhit->Nphot();
	                Int_t cell  = vhit->Cell();                                    
	                map[cell] += nPhot;
	                hit[cell] ++;
			
		  }
		} 
	     }
	     	
	       
        Int_t map2[64];      // cell to digits
        Int_t hit2[64];      // cell to digits
	Int_t j;
        Float_t time1[80];  
        Float_t time2[64];  
 	
	for (j=0; j<16; j++){
	map2[j] = map [j];
	hit2[j] = hit [j];
	//time2[j]=time1[j];
	}
	
	for (j=0; j<16; j++){
	map2[j+16] = map [2*j+16]+map [2*j+17];
	hit2[j+16] = hit [2*j+16]+hit [2*j+17];
	//time2[j+16]= TMath::Min(time1 [16 + 2*j], time1 [16 + 2*j + 1]);
	}
	
	for (j=32; j<64; j++){
	map2[j] =map[j+16];
	hit2[j] =hit[j+16];
	//time2[j] =time[j+16];
	}
		
	   fNdigits = 0;
	   fMulC    = 0;
	   fMulA    = 0;       
           for (Int_t i=0; i<64; i++) { 
                 Float_t q1 = Float_t ( map2[i] )* cPM[i] * kQe;
                 Float_t noise = gRandom->Gaus(10.5,3.22);
                 Float_t pmResponse  =  q1/kC*TMath::Power(ktheta/kthau,1/(1-ktheta/kthau)) 
                              + noise*1e-3;
                 map2[i] = Int_t( pmResponse * coef);
		 int test = 0;
                 if(map2[i] > 10) {map2[i] = Int_t( gRandom->Gaus(map2[i], 80));}		// charge smearing of MIP/4 -> sigma = MIP/4 = 120 (MIP = 480)
                 if(map2[i] > 240) {                                        			// cut at MIP/2 = 240
                       hCell->Fill(float(i));
		       hCellADC->Fill(float(i),map2[i]);		       
                       fNdigits++;
		       if(i<32) fMulC++;
		       if(i>31) fMulA++;
		       }
		       
	     } 
	     hMul0->Fill(fNdigits);
	     hMulC0->Fill(fMulC);
	     hMulA0->Fill(fMulA);
	     hMulAnd0->Fill(TMath::Min(fMulA, fMulC));
	     	     	
	     if(fMulA > 0){	
		nVOL++;}
		
	     if(fMulC > 0){	
		nVOR++;}
		
	     if(fMulA > 0 && fMulC > 0){	
		nVOLetR++;}
				
	     if(fMulA > 0 || fMulC > 0){	
		nVOLouR++;}      
 	
	     if(event%100==0) cout <<" event    = " << event <<endl;
//	     cout <<" multi   = " << fNdigits <<endl;   
        }
  	
	cout <<" nVOA     = " << nVOL <<endl;
	cout <<" nVOC     = " << nVOR <<endl;
	cout <<" nVOAandC = " << nVOLetR <<endl;
	cout <<" nVOAorC  = " << nVOLouR <<endl;
  		
//__________________________________________________
//      Fill root file

        TFile *histoFile = new TFile("Efficiencies.root","RECREATE");
   
    	hMapV0L->Write();
    	hMapV0R->Write();
	hpdgV0->Write();
	hvertexZ->Write();
	hTimC->Write();
	hTimA->Write();
	hCell->Write();
	hCellADC->Write();

	hMul0->Write();
	hMulC0->Write();
	hMulA0->Write();
	hMulAnd0->Write();
	histoFile->Close();   
	
}
