void ReadAODVertexingHF(const char *aodFileName="AliAOD.root",
			const char *aodHFFileName="AliAOD.VertexingHF.root") 
{
  //
  // Example macro to read D0->Kpi candidates from AOD (having the
  // standard AOD + a friend heavy-flavour AOD) and apply cuts
  // Origin: A.Dainese
  //
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3vertexingHF.so");

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

  Int_t nTotHF=0,nTotD0toKpi=0,nTot3Prong=0;
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
      }
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
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
    
  }
  
  printf("\n Total HF vertices: %d\n",nTotHF);
  printf("\n Total D0->Kpi: %d\n",nTotD0toKpi);
  printf("\n Total Charm->3Prong: %d\n",nTot3Prong);

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
//------------------------------------------------------------------------
void ReadAODVertexingHFsa(const char *aodHFFileName="AliAOD.VertexingHF.sa.root") 
{
  //
  // Example macro to read D0->Kpi candidates from a stand-alone
  // heavy-flavour AOD (i.e. without standard AOD) and apply cuts
  // Origin: A.Dainese
  //
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");

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
