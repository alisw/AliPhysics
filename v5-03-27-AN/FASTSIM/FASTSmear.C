TH1F *hpt;
void FASTSmear(char *finname="galice_upsi.root",char *foutname="FAST_upsi.root",Float_t bkg=0,Float_t nevmax=100000000){
 
  printf ("processing file %s\n",finname); 
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }    
  // create the structure to hold the variables for the branch
  
  struct fast_t {
    Float_t accp;
    Float_t effp;
    Float_t trigeffp;
    Float_t pgenp;
    Float_t psmearp; 
    Float_t thetagenp;
    Float_t thetasmearp;
    Float_t phigenp;
    Float_t phismearp;
    Float_t accm;
    Float_t effm;
    Float_t trigeffm;
    Float_t pgenm;
    Float_t psmearm; 
    Float_t thetagenm;
    Float_t thetasmearm;
    Float_t phigenm;
    Float_t phismearm;
  };

  fast_t fast;

  // open the input file
  TFile *fin=new TFile(finname);
  gAlice=(AliRun*)fin->Get("gAlice");

  // create the output file, its tree and branch 
  TFile *fout = new TFile(foutname,"RECREATE");
  TTree *FASTtrack = new TTree("FASTtrack","mu tracks smeared with fastsim");
  FASTtrack->Branch("fast",&fast.accp,"accp:effp:trigeffp:pgenp:psmearp:thetagenp:thetasmearp:phigenp:phismearp:accm:effm:trigeffm:pgenm:psmearm:thetagenm:thetasmearm:phigenm:phismearm");

  // create the fast tracker

  AliFastMuonTriggerEff *trigeff = new AliFastMuonTriggerEff();
  AliFastMuonTrackingAcc *acc = new AliFastMuonTrackingAcc();
  AliFastMuonTrackingEff *eff = new AliFastMuonTrackingEff();
  AliFastMuonTrackingRes *res = new AliFastMuonTrackingRes();
  acc->SetBackground(bkg);
  eff->SetBackground(bkg);
  res->SetBackground(bkg);  
  acc->Init(); 
  eff->Init(); 
  res->Init(); 
  AliFastDetector* tracker = new AliFastDetector("Tracker", "Muon Tracker");
  tracker->AddResponse(trigeff);
  tracker->AddResponse(acc);
  tracker->AddResponse(eff);
  tracker->AddResponse(res);
  tracker->Init();
  tracker->Dump();

  // loop over the events

  Int_t nev=AliRunLoader::GetNumberOfEvents();
  TParticle *mup, *mum;

  if (nev>nevmax) nev = nevmax;
  
  table = AliMUONFastTracking::Instance(); 
  printf ("background set to %g \n",table->GetBackground()); 
  for (Int_t iev=0; iev<nev; iev++) {
    Int_t npart=gAlice->GetEvent(iev);
    //    npart=500;
    for (Int_t ipart=0; ipart<npart; ipart+=3) {
      if (!(ipart%30)) printf ("Event #%d upsilon #%d \n",iev,ipart/3);
      part = gAlice->Particle(ipart+1);
      if (part->GetPdgCode() == 13) mum = part;
      else if (part->GetPdgCode() == -13) mup = part;
      part = gAlice->Particle(ipart+2);
      if (part->GetPdgCode() == 13) mum = part;
      else if (part->GetPdgCode() == -13) mup = part;

      // the mu+
      printf ("background set to %g \n",table->GetBackground()); 
      res->SetCharge(1);
      eff->SetCharge(1);
      acc->SetCharge(1);
      fast.pgenp     = mup->P();
      Double_t ptp    = fast.pgenp * TMath::Sin(mup->Theta());
      fast.thetagenp = 180./TMath::Pi() * mup->Theta();
      fast.phigenp   = 180./TMath::Pi() * mup->Phi();
      if (fast.phigenp>180) fast.phigenp -=360;
      res->Evaluate(fast.pgenp,fast.thetagenp,fast.phigenp,
		    fast.psmearp,fast.thetasmearp,fast.phismearp);
      fast.effp = eff->Evaluate(ptp,fast.thetagenp,fast.phigenp);
      fast.accp = acc->Evaluate(ptp,fast.thetagenp,fast.phigenp);
      fast.trigeffp = trigeff->Evaluate(1,ptp,fast.thetagenp,fast.phigenp);

      // the mu- 
      res->SetCharge(-1);
      acc->SetCharge(-1);
      eff->SetCharge(-1);
      fast.pgenm     = mum->P();
      Double_t ptm    = fast.pgenm * TMath::Sin(mum->Theta());
      fast.thetagenm = 180./TMath::Pi() * mum->Theta();
      fast.phigenm   = 180./TMath::Pi() * mum->Phi();
      if (fast.phigenm>180) fast.phigenm -=360;
      res->Evaluate(fast.pgenm,fast.thetagenm,fast.phigenm,
		    fast.psmearm,fast.thetasmearm,fast.phismearm);
      fast.effm = eff->Evaluate(ptm,fast.thetagenm,fast.phigenm);
      fast.accm = acc->Evaluate(ptm,fast.thetagenm,fast.phigenm);
      fast.trigeffm = trigeff->Evaluate(-1,ptm,fast.thetagenm,fast.phigenm);

      // fill the tree
      FASTtrack->Fill();
    }
  } 
  fout->Write();
  fout->Close();
  fin->Close();
}



