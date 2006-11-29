void analHits (const char *filename="galice.root",Int_t evNumber=0, char *opt="Liny"){
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L analHits.C               //this loads the macro in memory
//     Root > analHits();                 //by default process first event   
//     Root > analHits("galice2.root",2); //process third event from 
//                                          galice2.root file.
//Begin_Html
/*
<img src="picts/analHits.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////
  if(gAlice){
    delete gAlice;
    gAlice=0;
  }
  else{
    // Dynamically link some shared libs
    if(gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    } // end if
  }
// Connect the Root Galice file containing Geometry, Kine and Hits
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if(!file) file = new TFile(filename);

// Get AliRun object from file or create it if not on file
    if(!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if(gAlice) printf("AliRun object found on file\n");
	if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    } // end if !gAlice
      
// Set event pointer to this event
    Int_t nparticles = gAlice->GetEvent(evNumber);
    if (nparticles <= 0){
	cout << "No particles found for event " << evNumber;
	cout << " in file " << filename << endl;
	return;
    } // end if nparticles <=0

// Pointer to specific detector hits.
    AliFMDhit    *fmdHit;
    AliITShit    *itsHit;
    AliMUONHit   *muonHit;
    AliPHOSHit   *phosHit;
    AliPMDhit    *pmdHit;
    AliHMPIDHit   *richHit;
    AliSTARThit  *startHit;
    AliTOFhit    *tofHit;
    AliTPChit    *tpcHit;
    AliTPCTrackHits *tpc2Hit;
    AliTRDhit    *trdHit;
    AliCASTORhit *castorHit;
    AliZDCHit    *zdcHit;
    AliEMCALHit  *emcalHit;

// Get pointers to ALL Alice detectors and Hits containers
    AliFMD    *FMD    = (AliFMD*)    gAlice->GetDetector("FMD");
    AliITS    *ITS    = (AliITS*)    gAlice->GetDetector("ITS");
    AliMUON   *MUON   = (AliMUON*)   gAlice->GetDetector("MUON");
    AliPHOS   *PHOS   = (AliPHOS*)   gAlice->GetDetector("PHOS");
    AliPMD    *PMD    = (AliPMD*)    gAlice->GetDetector("PMD");
    AliHMPID   *HMPID   = (AliHMPID*)   gAlice->GetDetector("HMPID");
    AliSTART  *START  = (AliSTART*)  gAlice->GetDetector("START");
    AliTOF    *TOF    = (AliTOF*)    gAlice->GetDetector("TOF");
    AliTPC    *TPC    = (AliTPC*)    gAlice->GetDetector("TPC");
    AliTPC    *TPC    = (AliTPC*)    gAlice->GetDetector("TPC");
    AliTRD    *TRD    = (AliTRD*)    gAlice->GetDetector("TRD");
    AliCASTOR *CASTOR = (AliCASTOR*) gAlice->GetDetector("CASTOR");
    AliZDC    *ZDC    = (AliZDC*)    gAlice->GetDetector("ZDC");
    AliEMCAL  *EMCAL  = (AliEMCAL*)  gAlice->GetDetector("EMCAL");

// Get pointer to the particles
//    TClonesArray *Particles = gAlice->Particles();
//    TParticle    *part;

    // Create histograms
    if(FMD)   TH1F *hFMD   = new TH1F("hFMD"   ,"Hit Radius",100,0.,100.);
    if(ITS)   TH1F *hITS   = new TH1F("hITS"   ,"Ionization",100,0.,3.e-3);
    if(MUON)  TH1F *hMUON  = new TH1F("hMUON"  ,"Hit Radius",100,0.,500.);
    if(PHOS)   TH1F *hPHOS  = new TH1F("hPHOS"  ,"Energy Dep.",100,0.,0.5);
    if(PMD)     TH1F *hPMD   = new TH1F("hPMD"   ,"Energy Dep.",100,0.,1.e+5);
    if(HMPID)  TH1F *hHMPID  = new TH1F("hHMPID"  ,"Energy loss",100,0.,1.e-5);
    if(START) TH1F *hSTART = new TH1F("hSTART" ,"Time of Flight",100,0.,10.);
    if(TOF)   TH1F *hTOF   = new TH1F("hTOF"   ,"Time of Flight",100,0.,1.e-5);
    if(TPC)   TH1F *hTPC   = new TH1F("hTPC"   ,"Charge",100,0.,70.2);
    if(TRD)   TH1F *hTRD   = new TH1F("hTRD"   ,"Charge",100,0.,10.);
    if(CASTOR)TH1F *hCASTOR= new TH1F("hCASTOR","Hit Radius",100,0.,10.);
    if(ZDC)   TH1F *hZDC   = new TH1F("hZDC"   ,"Energy",100,0.,5.);
    if(EMCAL) TH1F *hEMCAL = new TH1F("hEMCAL" ,"Energy",100,0.,2.);
//    TH1F *hTPAR  = new TH1F("hTPAR" ,"?",6,1,7);  
    Int_t track,ntracks = gAlice->TreeH()->GetEntries();
// Start loop on tracks in the hits containers
    for(track=0; track<ntracks;track++){
      //MI change
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(track);
      if(FMD){
	    for(fmdHit=(AliFMDhit*)FMD->FirstHit(-1);fmdHit;
		fmdHit=(AliFMDhit*)FMD->NextHit()){
		hFMD->Fill(TMath::Hypot(fmdHit->X(),fmdHit->Y()));
	    } // end for fmdHit
	} // end if FMD
	if(ITS){
	    for(itsHit=(AliITShit*)ITS->FirstHit(-1);itsHit;
		itsHit=(AliITShit*)ITS->NextHit()){
		if(itsHit->GetIonization()>0.0){//only after a step in the ITS
		    hITS->Fill(itsHit->GetIonization());
		} // end if
	    } // end for itsHit
	} // end if ITS
	if(MUON){
	    for(muonHit=(AliMUONHit*)MUON->FirstHit(-1);muonHit;
		muonHit=(AliMUONHit*)MUON->NextHit()){
		hMUON->Fill(TMath::Hypot(muonHit->X(),muonHit->Y()));
	    } // end for muonHit
	} // end if MUON
	if(PHOS){
	    for(phosHit=(AliPHOSHit*)PHOS->FirstHit(-1);phosHit;
		phosHit=(AliPHOSHit*)PHOS->NextHit()){
		hPHOS->Fill(phosHit->GetEnergy());
	    } // end for phosHit
	} // end if PHOS
	if(PMD){
	    for(pmdHit=(AliPMDhit*)PMD->FirstHit(-1);pmdHit;
		pmdHit=(AliPMDhit*)PMD->NextHit()){
		hPMD->Fill(pmdHit->GetEnergy());
	    } // end for pmdHit
	} // end if PMD
	if(HMPID){
	    for(richHit=(AliHMPIDHit*)HMPID->FirstHit(-1);richHit;
		richHit=(AliHMPIDHit*)HMPID->NextHit()){
		hHMPID->Fill(richHit->fEloss);
	    } // end for richHit
	} // end if HMPID
	if(START){
	    for(startHit=(AliSTARThit*)START->FirstHit(-1);startHit;
		startHit=(AliSTARThit*)START->NextHit()){
		hSTART->Fill(startHit->fTime);
	    } // end for startHit
	} // end if START
	if(TOF){
	    for(tofHit=(AliTOFhit*)TOF->FirstHit(-1);tofHit;
		tofHit=(AliTOFhit*)TOF->NextHit()){
		hTOF->Fill(tofHit->GetTof());
	    } // end for tofHit
	} // end if TOF

	if(TRD) {
	    for(trdHit=(AliTRDhit*)TRD->FirstHit(-1);trdHit;
		trdHit=(AliTRDhit*)TRD->NextHit()) {
		hTRD->Fill((Float_t)(trdHit->GetCharge()));
	    } // end for
	} // end if TRD
	if(CASTOR) {
	    for(castorHit=(AliCASTORhit*)CASTOR->FirstHit(-1);castorHit;
		castorHit=(AliCASTORhit*)CASTOR->NextHit()) {
		hCASTOR->Fill(TMath::Hypot(castorHit->X(),castorHit->Y()));
	    } // end for
	} // end if CASTOR
	if(ZDC){
	    for(zdcHit=(AliZDCHit*)ZDC->FirstHit(-1);zdcHit;
		zdcHit=(AliZDCHit*)ZDC->NextHit()){
		hZDC->Fill(zdcHit->GetEnergy());
	    } // end for zdcdHit
	} // end if ZDC
	if(TPC) {	 
	  for(tpcHit=(AliTPChit*)TPC->FirstHit(-1);tpcHit;
	      tpcHit=(AliTPChit*)TPC->NextHit()) {
	    hTPC->Fill((Float_t)(tpcHit->fQ));
	  } // end for tpcHit
	} // end if TPC
	if(EMCAL) {
	  for(emcalHit=(AliEMCALHit*)EMCAL->FirstHit(-1);emcalHit;
	      emcalHit=(AliEMCALHit*)EMCAL->NextHit()) {
	    hEMCAL->Fill((Float_t)(emcalHit->GetEnergy()));
	  } // end for tpcHit
	} // end if TPC

    } // end for track

//Create a canvas, set the view range, show histograms
    TCanvas *c0 = new TCanvas("c0","Alice Detectors",400,10,600,700);
    if(opt=="Logy")c0->SetLogy();
    if(FMD){
	hFMD->SetFillColor(42);
	hFMD->Draw();
	c0->Print("analHitsFMD.ps");
    } // end if FMD
    if(ITS){
	hITS->SetFillColor(42);
	hITS->Draw();
	c0->SaveAs("analHitsITS.ps");
    } // end if ITS
    if(MUON){
	hMUON->SetFillColor(42);
	hMUON->Draw();
	c0->SaveAs("analHitsMUON.ps");
    } // end if MUON
    if(PHOS){
	hPHOS->SetFillColor(42);
	hPHOS->Draw();
	c0->SaveAs("analHitsPHOS.ps");
    } // end if PHOS
    if(PMD){
	hPMD->SetFillColor(42);
	hPMD->Draw();
	c0->SaveAs("analHitsPMD.ps");
    } // end if PMD
    if(HMPID){
	hHMPID->SetFillColor(42);
	hHMPID->Draw();
	c0->SaveAs("analHitsHMPID.ps");
    } // end if HMPID
    if(START){
	hSTART->SetFillColor(42);
	hSTART->Draw();
	c0->SaveAs("analHitsSTART.ps");
    } // end if START
    if(TOF){
	hTOF->SetFillColor(42);
	hTOF->Draw();
	c0->SaveAs("analHitsTOF.ps");
    } // end if TOF
    if(TPC){
	hTPC->SetFillColor(42);
	hTPC->Draw();
	c0->SaveAs("analHitsTPC.ps");
    } // end if TPC
    if(TRD){
	hTRD->SetFillColor(42);
	hTRD->Draw();
	c0->SaveAs("analHitsTRD.ps");
    } // end if TRD
    if(CASTOR){
	hCASTOR->SetFillColor(42);
	hCASTOR->Draw();
	c0->SaveAs("analHitsCASTOR.ps");
    } // end if TRD
    if(ZDC){
	hZDC->SetFillColor(42);
	hZDC->Draw();
	c0->SaveAs("analHitsZDC.ps");
    } // end if ZDC
    if(EMCAL){
	hEMCAL->SetFillColor(42);
	hEMCAL->Draw();
	c0->SaveAs("analHitsEMCAL.ps");
    } // end if ZDC

// Clean Up
    /*
    if(FMD)    delete hFMD;
    if(ITS)    delete hITS;
    if(MUON)   delete hMUON;
    if(PHOS)   delete hPHOS;
    if(PMD)    delete hPMD;
    if(HMPID)   delete hHMPID;
    if(START)  delete hSTART;
    if(TOF)    delete hTOF;
    if(TPC)    delete hTPC;
    if(TRD)    delete hTRD;
    if(CASTOR) delete hCASTOR;
    if(ZDC)    delete hZDC;
    */
}
