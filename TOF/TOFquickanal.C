Int_t TOFquickanal(Int_t eventNumber = 0)
{
  /////////////////////////////////////////////////////////////////////////
  //   This macro is a small example of a ROOT macro
  //   illustrating how to read the output of GALICE
  //   and fill some histograms concerning the TOF Hit Tree.
  //
  //     Root > .L TOFquickanal.C   //this loads the macro in memory
  //     Root > TOFquickanal();     //by default process first event
  //     Root > TOFquickanal(2);    //process third event
  //Begin_Html
  /*
    <img src="picts/TOFquickanal.gif">
  */
  //End_Html
  //
  // Author: F. Pierella , Bologna University 12-04-2001
  // Updated to the new I/O by: A. De Caro, C. Zampolli
  /////////////////////////////////////////////////////////////////////////
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  Int_t rc = 0;
  
  AliRunLoader *rl =AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,"update");
  if (!rl) 
    {
      cerr << "Can't load RunLoader from file!\n";
      rc = 1;
      return rc;
    }

  rl->LoadgAlice();
  gAlice=rl->GetAliRun();

  if (!gAlice)
    {
      cerr << "<TOFquickanal> AliRun object not found on file \n";
      rc = 2;
      return rc;
    }

  // Get the pointer to the TOF detector
  AliLoader *tofl = rl->GetLoader("TOFLoader");
  AliTOF * tof = (AliTOF*) gAlice->GetDetector("TOF");
  if (tof == 0x0 || tofl == 0x0) {
    cerr << "<TOFquickanal> Can not find TOF or TOFLoader\n";
    rc = 3;
    return rc;
  }

  //=======> Create histograms
  //---> Time of Flight for Primary Particles (ns)
  TH1F *htofprim = new TH1F("htofprim","Time of Flight for Primary Particles",100,0.,100.);
  //--->Time of Flight for Secondary Particles (ns)
  TH1F *htofsec  = new TH1F("htofsec","Time of Flight for Secondary Particles",100,0.,100.);
  
  //---> r (radius) coordinate of production in the ALICE frame for secondary particles that produce at 
  //     least one TOF-hit (cm) - cylindrical coordinate system assumed, primary plus secondary-
  TH1F *hradius = new TH1F("hradius","r (radius) coordinate at the production vertex for secondary particles with at least one TOF-Hit",50,0.,500.);
  
  //---> Momentum of primary particles that produce (at least) one TOF-hit when the hit
  //     is produced (Gev/c)
  TH1F *htofmom  = new TH1F("htofmom","Momentum of primary particles when the Hit is produced",50,0.,5.);
  
  //---> Momentum of primary particles that produce (at least) one TOF-hit at the production vertex
  //     (Gev/c)
  TH1F *hprodmom  = new TH1F("hprodmom","Momentum of primary particles (with at least one TOF hit) at the production ",50,0.,5.); 
  
  //---> Theta of production for primary particles that produce (at least) one TOF-hit (deg)
  TH1F *hprodthe  = new TH1F("hprodthe","Theta of primary particles (with at least one TOF hit) at the production ",90,0.,180.);
  
  //---> Phi of production for primary particles that produce (at least) one TOF-hit (deg)
  TH1F *hprodphi  = new TH1F("hprodphi","Phi of primary particles (with at least one TOF hit) at the production ",180,-180.,180.);
  
  //---> z Coordinate of the TOF Hit (z beam axis) - primary plus secondary - (cm)
  TH1F *hzcoor = new TH1F("hzcoor","z Coordinate of the TOF Hit",800,-400.,400.);
  
  //---> Incidence Angle of the particle on the pad (or strip) (deg)  - primary plus secondary - 
  TH1F *hincangle = new TH1F("hincangle","Incidence Angle of the particle on the strip",90,0.,180.);
  
  printf ("Processing event %d \n", eventNumber);
  rl->GetEvent(eventNumber);
  
  // Get pointers to Alice detectors and Hits containers
  tofl->LoadHits();
  TTree *TH = tofl->TreeH();
  tof->SetTreeAddress();
  if (!TH) {
    cout << "<TOFquickanal> No hit tree found" << endl;
    rc = 4;
    return rc;
  }
  
  // Import the Kine Tree for the event eventNumber in the file  
  rl->LoadHeader();
  rl->LoadKinematics();
  //AliStack * stack = rl->Stack();
  
  Int_t ntracks = TH->GetEntries();
  cout<<" ntracks = "<<ntracks<<endl;
  
  AliTOFhitT0 *tofHit;
  
  // Start loop on tracks in the hits containers
  for (Int_t track=0; track<ntracks;track++) {
    
    tof->ResetHits();
    TH->GetEvent(track);
    
    for(tofHit=(AliTOFhitT0*)tof->FirstHit(track); tofHit; tofHit=(AliTOFhitT0*)tof->NextHit()) {
      
      Float_t toflight = tofHit->GetTof();
      toflight        *= 1.E+09;  // conversion from s to ns
      Double_t tofmom  = tofHit->GetMom();
      
      Int_t ipart = tofHit->Track();
      TParticle *particle = gAlice->Particle(ipart);
      if (particle->GetFirstMother() < 0) {
	htofprim->Fill(toflight);
	htofmom->Fill(tofmom); 
      } else {
	htofsec->Fill(toflight); 
      }
      
      Double_t zcoor = tofHit->Z();
      hzcoor->Fill(zcoor);
      
      Double_t incangle = tofHit->GetIncA();
      hincangle->Fill(incangle);
      
      Double_t xcoor  = particle->Vx();
      Double_t ycoor  = particle->Vy();
      Double_t radius = TMath::Sqrt(xcoor*xcoor+ycoor*ycoor);
      if (particle->GetFirstMother() >= 0) hradius->Fill(radius);
      
      Double_t prodmom = particle->P();        
      if (prodmom!=0.) {
	Double_t dummy = (particle->Pz())/prodmom;
	Double_t prodthe = TMath::ACos(dummy);
	prodthe *= 57.29578; // conversion from rad to deg
	if (particle->GetFirstMother() < 0) hprodthe->Fill(prodthe);
      } // theta at production vertex
      
      if (particle->GetFirstMother() < 0) {         
	hprodmom->Fill(prodmom);
	Double_t dummypx = particle->Px();
	Double_t dummypy = particle->Py();
	Double_t prodphi = TMath::ATan2(dummypy,dummypx);
	prodphi *= 57.29578; // conversion from rad to deg
	hprodphi->Fill(prodphi);
      } // phi at production vertex
    } // close loop on TOF-hits
  } // close loop on tracks in the hits containers
  
  //Create  canvas, set the view range, show histograms
  TCanvas *c1 = new TCanvas("c1","Alice TOF hits quick analysis",400,10,600,700);
  c1->cd();
  hprodmom->Draw();
  
  TCanvas *c2 = new TCanvas("c2","Alice TOF hits quick analysis",400,10,600,700);
  c2->cd();
  hprodthe->Draw();
  
  TCanvas *c3 = new TCanvas("c3","Alice TOF hits quick analysis",400,10,600,700);
  c3->cd();
  hprodphi->Draw();
  
  TCanvas *c4 = new TCanvas("c4","Alice TOF hits quick analysis",400,10,600,700);
  c4->cd();
  hzcoor->Draw();
  
  TCanvas *c5 = new TCanvas("c5","Alice TOF hits quick analysis",400,10,600,700);
  c5->cd();
  hradius->Draw();
  
  TCanvas *c6 = new TCanvas("c6","Alice TOF hits quick analysis",400,10,600,700);
  c6->cd();
  htofprim->Draw();
  
  TCanvas *c7 = new TCanvas("c7","Alice TOF hits quick analysis",400,10,600,700);
  c7->cd();
  htofsec->Draw();
  
  
  TCanvas *c8 = new TCanvas("c8","Alice TOF hits quick analysis",400,10,600,700);
  c8->cd();
  htofmom->Draw();
  
  TCanvas *c9 = new TCanvas("c9","Alice TOF hits quick analysis",400,10,600,700);
  c9->cd();
  hincangle->Draw();
  
  //tofl->UnloadHits();
  //rl->UnloadHeader();
  //rl->UnloadgAlice();
  //rl->UnloadKinematics();

  return rc;

}
