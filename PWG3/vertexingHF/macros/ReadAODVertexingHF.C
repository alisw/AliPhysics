void ReadAODVertexingHF(const char *aodFileName="AliAOD.root",
			const char *aodHFFileName="AliAOD.VertexingHF.root") 
{
  //
  // Example macro to read D0->Kpi candidates from AOD (having the
  // standard AOD + a friend heavy-flavour AOD) and apply cuts
  // Origin: A.Dainese
  //

  Bool_t useParFiles=kFALSE;
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  // create a test histogram
  TH2F *hCPtaVSd0d0 = new TH2F("hCPtaVSd0d0","D^{0} correlation plot",1000,-50000,50000,1000,-1,1);
  hCPtaVSd0d0->SetXTitle("Product of impact parameters [#mu m^{2}]");
  hCPtaVSd0d0->SetYTitle("Cosine of pointing angle");
  TH1F *hMass = new TH1F("hMass","D^{0} mass plot",100,1.7,2);
  hMass->SetXTitle("Invariant mass [GeV]");
  hMass->SetYTitle("Entries");
  TH1F *hSecVtxZ = new TH1F("hSecVtxZ","D^{0} decay vertex z",1000,-10,10);
  hSecVtxZ->SetXTitle("z of decay vertex [cm]");
  hSecVtxZ->SetYTitle("Entries");
  TH1F *hDeltaMassDstar = new TH1F("hDeltaMassDstar","D* delta mass plot",100,0,0.3);
  hDeltaMassDstar->SetXTitle("M(Kpipi)-M(Kpi) [GeV]");
  hDeltaMassDstar->SetYTitle("Entries");


  // open input file and get the TTree
  TFile inFile(aodFileName,"READ");
  if (!inFile.IsOpen()) return;

  TTree *aodTree = (TTree*)inFile.Get("aodTree");
  aodTree->AddFriend("aodTree",aodHFFileName);

  AliAODEvent *aod = new AliAODEvent();

  aod->ReadFromTree(aodTree);

  // load heavy flavour vertices
  TClonesArray *arrayVerticesHF = 
    (TClonesArray*)aod->GetList()->FindObject("VerticesHF"); 
  
  // load D0->Kpi candidates
  TClonesArray *arrayD0toKpi = 
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi"); 
     
  // load 3prong candidates
  TClonesArray *array3Prong = 
    (TClonesArray*)aod->GetList()->FindObject("Charm3Prong"); 
    
  // load D* candidates
  TClonesArray *arrayDstar = 
    (TClonesArray*)aod->GetList()->FindObject("Dstar"); 
    
  // load cascade (V0+track) candidates
  TClonesArray *arrayCascades = 
    (TClonesArray*)aod->GetList()->FindObject("CascadesHF"); 
    

  Double_t cutsD0[9]=
  // cutsD0[0] = inv. mass half width [GeV]
  // cutsD0[1] = dca [cm]
  // cutsD0[2] = cosThetaStar
  // cutsD0[3] = pTK [GeV/c]
  // cutsD0[4] = pTPi [GeV/c]
  // cutsD0[5] = d0K [cm]   upper limit!
  // cutsD0[6] = d0Pi [cm]  upper limit!
  // cutsD0[7] = d0d0 [cm^2]
  // cutsD0[8] = cosThetaPoint
                     {1000.,
		      100000.,
		      1.1,
		      0.,
		      0.,
		      100000.,
		      100000.,
		      100000000.,
		      -1.1}; 
  Double_t cutsDstar[5]=
    // (to be passed to AliAODRecoCascadeHF::SelectDstar())
    // 0 = inv. mass half width of D* [GeV]
    // 1 = half width of (M_Kpipi-M_Kpi) [GeV]
    // 2 = PtMin of pi_s [GeV/c]
    // 3 = PtMax of pi_s [GeV/c]
    // 4 = theta, angle between the trace of pi_s and D0 decay plane [rad]
                     {999999.,
		      999999.,
		      0.1, 
		      1.0, 
		      0.5};

  Double_t cutsLctoV0[9]=// cuts on Lambdac candidates to V0+bachelor
                        // (to be passed to AliAODRecoDecayHF3Prong::SelectLctoV0())
                        // 0 = inv. mass half width in K0s hypothesis [GeV]   
                        // 1 = inv. mass half width in Lambda hypothesis [GeV]   
                        // 2 = inv. mass V0 in K0s hypothesis half width [GeV]   
                        // 3 = inv. mass V0 in Lambda hypothesis half width [GeV]   
                        // 4 = pT min Bachelor track [GeV/c]
                        // 5 = pT min V0-Positive track [GeV/c]
                        // 6 = pT min V0-Negative track [GeV/c]
                        // 7 = dca cut on the V0 (cm)
                        // 8 = dca cut on the cascade (cm)
    {2.0,2.0,1.0,1.0,0.0,0.0,0.0,1000.,1000.};

  Int_t nTotHF=0,nTotD0toKpi=0,nTotDstar=0,nTot3Prong=0,nTotCasc=0;
  AliAODVertex *vtx1=0;

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout<<"\n------------ Event: "<<nEv<<" ------------------"<<endl;

    // read event
    aodTree->GetEvent(nEv);
    //aodTree->BranchRef();

    //print event info
    aod->GetHeader()->Print();

    // primary vertex
    vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    vtx1->Print();

    // make trkIDtoEntry register (temporary)
    Int_t trkIDtoEntry[100000];
    for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
      AliAODTrack *track = aod->GetTrack(it);
      trkIDtoEntry[track->GetID()]=it;
    }
    

    // loop over D0->Kpi candidates
    Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
    nTotD0toKpi += nD0toKpi;
    cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
    
    for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);

      Bool_t unsetvtx=kFALSE;
      if(!d->GetOwnPrimaryVtx()) {
	d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	unsetvtx=kTRUE;
      }
      Int_t okD0=0,okD0bar=0; 
      if(d->SelectD0(cutsD0,okD0,okD0bar)) {
	//cout<<1e8*d->Prodd0d0()<<endl;
	hMass->Fill(d->InvMassD0(),0.5);
	hMass->Fill(d->InvMassD0bar(),0.5);
	hCPtaVSd0d0->Fill(1e8*d->Prodd0d0(),d->CosPointingAngle());
	hSecVtxZ->Fill(d->GetSecVtxZ());
	//cout<<d->GetSecVtxX() <<endl;

	// get daughter AOD tracks
	AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
	AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
	if(!trk0 || !trk1) {
	  trk0=aod->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
	  trk1=aod->GetTrack(trkIDtoEntry[d->GetProngID(1)]);
	  cout<<"references to standard AOD not available"<<endl;
	}
	cout<<"pt of positive track: "<<trk0->Pt()<<endl;

	// make a AliNeutralTrackParam from the D0 
	// and calculate impact parameters to primary vertex
	AliNeutralTrackParam trackD0(d);
	cout<<"pt of D0: "<<d->Pt()<<" (AliAODRecoDecay); "<<trackD0.Pt()<<" (track)"<<endl;
	//trackD0.Print();
	Double_t dz[2],covdz[3];
	trackD0.PropagateToDCA(vtx1,aod->GetMagneticField(),1000.,dz,covdz);
	cout<<"D0 impact parameter rphi: "<<dz[0]<<" +- "<<TMath::Sqrt(covdz[0])<<endl;
      }
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
    }


    // loop over D* candidates
    Int_t nDstar = arrayDstar->GetEntriesFast();
    nTotDstar += nDstar;
    cout<<"Number of D*->D0pi: "<<nDstar<<endl;
    
    for (Int_t iDstar = 0; iDstar < nDstar; iDstar++) {
      AliAODRecoCascadeHF *c = (AliAODRecoCascadeHF*)arrayDstar->UncheckedAt(iDstar);
      Bool_t unsetvtx=kFALSE;
      if(!c->GetOwnPrimaryVtx()) {
	c->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	c->Get2Prong()->SetOwnPrimaryVtx(vtx1);
	unsetvtx=kTRUE;
      }
      if(c->SelectDstar(cutsDstar,cutsD0)) {
	hDeltaMassDstar->Fill(c->DeltaInvMass());
	// get daughters
	AliAODTrack *trk = (AliAODTrack*)c->GetBachelor();
	AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)c->Get2Prong();
	cout<<"pt of soft pi: "<<trk->Pt()<<endl;
	cout<<"pt of D0: "<<d->Pt()<<endl;
      }

      if(unsetvtx) {
	c->UnsetOwnPrimaryVtx();
	c->Get2Prong()->UnsetOwnPrimaryVtx();
      }
    }

    
    // count 3prong candidates
    Int_t n3Prong = array3Prong->GetEntriesFast();
    nTot3Prong += n3Prong;
    cout<<"Number of Charm->3Prong: "<<n3Prong<<endl;
    
    // loop over HF vertices
    Int_t nVtxsHF = arrayVerticesHF->GetEntriesFast();
    nTotHF += nVtxsHF;
    cout<<"Number of heavy-flavour vertices: "<<nVtxsHF<<endl;
    for (Int_t iVtx = 0; iVtx < nVtxsHF; iVtx++) {
      AliAODVertex *vtxHF = (AliAODVertex*)arrayVerticesHF->UncheckedAt(iVtx);
      // print info
      //cout << iVtx << ": vertex z position: " << vtxHF->GetZ() << endl;
    }


    // count cascade candidates
    if (arrayCascades){
      Int_t nCasc = arrayCascades->GetEntriesFast();
      nTotCasc+=nCasc;
      cout << "Number of Cascades: "<<nCasc<<endl;
    }
    
  }
  
  printf("\n Total HF vertices: %d\n",nTotHF);
  printf("\n Total D0->Kpi: %d\n",nTotD0toKpi);
  printf("\n Total D*->D0pi: %d\n",nTotDstar);
  printf("\n Total Charm->3Prong: %d\n",nTot3Prong);
  if (arrayCascades) printf("\n Total Cascades: %d\n",nTotCasc);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,700);
  c1->Divide(2,2);
  c1->cd(1);
  hCPtaVSd0d0->Draw("colz");
  c1->cd(2);
  hMass->SetFillColor(4);
  hMass->Draw();
  c1->cd(3);
  hSecVtxZ->SetFillColor(2);
  hSecVtxZ->Draw();
  c1->cd(4);
  hDeltaMassDstar->SetFillColor(3);
  hDeltaMassDstar->Draw();

  return;
}
//------------------------------------------------------------------------
void ReadAODVertexingHFsa(const char *aodHFFileName="AliAOD.VertexingHF.sa.root") 
{
  //
  // Example macro to read D0->Kpi candidates from a stand-alone
  // heavy-flavour AOD (i.e. without standard AOD) and apply cuts
  // Origin: A.Dainese
  //
  Bool_t useParFiles=kFALSE;
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  // create a test histogram
  TH2F *hCPtaVSd0d0 = new TH2F("hCPtaVSd0d0","D^{0} correlation plot",1000,-50000,50000,1000,-1,1);
  hCPtaVSd0d0->SetXTitle("Product of impact parameters [#mu m^{2}]");
  hCPtaVSd0d0->SetYTitle("Cosine of pointing angle");
  TH1F *hMass = new TH1F("hMass","D^{0} mass plot",100,1.7,2);
  hMass->SetXTitle("Invariant mass [GeV]");
  hMass->SetYTitle("Entries");
  TH1F *hSecVtxZ = new TH1F("hSecVtxZ","D^{0} decay vertex z",1000,-10,10);
  hSecVtxZ->SetXTitle("z of decay vertex [cm]");
  hSecVtxZ->SetYTitle("Entries");

  // open input file and get the TTree
  TFile inFile(aodHFFileName,"READ");
  if (!inFile.IsOpen()) return;

  TTree *aodTree = (TTree*)inFile.Get("aodTree");

  AliAODEvent *aod = new AliAODEvent();

  aod->ReadFromTree(aodTree);

  // load heavy flavour vertices
  TClonesArray *arrayVerticesHF = 
    (TClonesArray*)aod->GetList()->FindObject("VerticesHF"); 
  
  // load D0->Kpi candidates
  TClonesArray *arrayD0toKpi = 
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi"); 
     

  Double_t cutsD0[9]=
  // cutsD0[0] = inv. mass half width [GeV]
  // cutsD0[1] = dca [cm]
  // cutsD0[2] = cosThetaStar
  // cutsD0[3] = pTK [GeV/c]
  // cutsD0[4] = pTPi [GeV/c]
  // cutsD0[5] = d0K [cm]   upper limit!
  // cutsD0[6] = d0Pi [cm]  upper limit!
  // cutsD0[7] = d0d0 [cm^2]
  // cutsD0[8] = cosThetaPoint
                     {1000.,
		      100000.,
		      1.1,
		      0.,
		      0.,
		      100000.,
		      100000.,
		      100000000.,
		      -1.1}; 

  Int_t nTotHF=0,nTotD0toKpi=0;
  AliAODVertex *vtx1=0;

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout<<"\n------------ Event: "<<nEv<<" ------------------"<<endl;

    // read event
    aodTree->GetEvent(nEv);
    //aodTree->BranchRef();

    
    // loop over D0->Kpi candidates
    Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
    nTotD0toKpi += nD0toKpi;
    cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
    
    for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
      Int_t okD0=0,okD0bar=0; 
      if(d->SelectD0(cutsD0,okD0,okD0bar)) {
	//cout<<1e8*d->Prodd0d0()<<endl;
	hMass->Fill(d->InvMassD0(),0.5);
	hMass->Fill(d->InvMassD0bar(),0.5);
	hCPtaVSd0d0->Fill(1e8*d->Prodd0d0(),d->CosPointingAngle());
	hSecVtxZ->Fill(d->GetSecVtxZ());
	//cout<<d->GetSecVtxZ() <<endl;

      }
    }
    
  }
  
  printf("\n Total D0->Kpi: %d\n",nTotD0toKpi);

  TCanvas *c = new TCanvas("c","c",0,0,1000,1000);
  c->Divide(2,2);
  c->cd(1);
  hCPtaVSd0d0->Draw("colz");
  c->cd(2);
  hMass->SetFillColor(4);
  hMass->Draw();
  c->cd(3);
  hSecVtxZ->SetFillColor(2);
  hSecVtxZ->Draw();

  return;
}
