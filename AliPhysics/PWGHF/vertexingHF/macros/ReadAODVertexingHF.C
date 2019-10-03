void ReadAODVertexingHF(const char *aodFileName="AliAOD.root",
			const char *aodHFFileName="AliAOD.VertexingHF.root") 
{
  //
  // Example macro to read D0->Kpi candidates from AOD (having the
  // standard AOD + a friend heavy-flavour AOD) and apply cuts
  // Origin: A.Dainese
  //

  TStopwatch t;
  t.Start();

  Bool_t useParFiles=kFALSE;
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");
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
    


  Int_t nTotHF=0,nTotD0toKpi=0,nTotDstar=0,nTot3Prong=0,nTotCasc=0;
  AliAODVertex *vtx1=0;

  AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi("CutsD0toKpi");
  cutsD0toKpi->SetStandardCutsPP2010();
  cutsD0toKpi->SetRemoveDaughtersFromPrim(kFALSE);

  AliRDHFCutsDStartoKpipi *cutsDStar = new AliRDHFCutsDStartoKpipi("CutsDStartoKpipi");
  cutsDStar->SetStandardCutsPP2010();

  AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = new AliRDHFCutsDplustoKpipi("CutsDplustoKpipi");
  cutsDplustoKpipi->SetStandardCutsPP2010();
  cutsDplustoKpipi->SetRemoveDaughtersFromPrim(kFALSE);
  
  Int_t nTot3ProngSele=0;
  Int_t nTotD0toKpiSele=0;
  Int_t nTotDStarSele=0;

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

    /*
    // make trkIDtoEntry register (temporary)
    Int_t trkIDtoEntry[100000];
    for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
      AliAODTrack *track = aod->GetTrack(it);
      trkIDtoEntry[track->GetID()]=it;
    }
    */

    // Fix references to daughter tracks
    //AliAnalysisVertexingHF *fixer = new AliAnalysisVertexingHF();
    //fixer->FixReferences(aod);
    //delete fixer;
    //
    //AliRDHFCutsD0toKpi *mycuts=new AliRDHFCutsD0toKpi();
    //mycuts->SetFixRefs(kTRUE);
    //mycuts->IsEventSelected(aod);

    // loop over D0->Kpi candidates
    Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
    nTotD0toKpi += nD0toKpi;
    cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
    
    for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);

      d->SetOwnPrimaryVtx(vtx1);

      /*
	// get daughter AOD tracks
	AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
	AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);

	if(trk0->GetStatus()) printf("ok %d\n",iD0toKpi);
      */
      printf("D0 %d\n",iD0toKpi);
      if(cutsD0toKpi->IsSelected(d,AliRDHFCuts::kAll)) {printf("D0 %d passed\n",iD0toKpi); nTotD0toKpiSele++;}

    }


    // loop over D* candidates
    Int_t nDstar = arrayDstar->GetEntriesFast();
    nTotDstar += nDstar;
    cout<<"Number of D*->D0pi: "<<nDstar<<endl;
    

  
    for (Int_t iDstar = 0; iDstar < nDstar; iDstar++) {
      AliAODRecoCascadeHF *c = (AliAODRecoCascadeHF*)arrayDstar->UncheckedAt(iDstar);
      printf("D* %d\n",iDstar);
      if(cutsDStar->IsSelected(c,AliRDHFCuts::kCandidate)) {printf("D* %d passed\n",iDstar); nTotDStarSele++;}

    }

  

    // count 3prong candidates
    Int_t n3Prong = array3Prong->GetEntriesFast();
    nTot3Prong += n3Prong;
    cout<<"Number of Charm->3Prong: "<<n3Prong<<endl;

    for (Int_t i3p = 0; i3p < n3Prong; i3p++) {
      AliAODRecoDecayHF3Prong *ccc = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3p);
      printf("3p %d\n",i3p);
      //if(cutsDplustoKpipi->IsSelected(ccc,AliRDHFCuts::kCandidate)) {printf("3p %d passed\n",i3p); nTot3ProngSele++;}

    }

    /*
    // loop over HF vertices
    Int_t nVtxsHF = arrayVerticesHF->GetEntriesFast();
    nTotHF += nVtxsHF;
    cout<<"Number of heavy-flavour vertices: "<<nVtxsHF<<endl;
    for (Int_t iVtx = 0; iVtx <0; iVtx++) {
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
    */
  }
  
  printf("\n Total HF vertices: %d\n",nTotHF);
  printf("\n Total D0->Kpi: %d; selected %d\n",nTotD0toKpi,nTotD0toKpiSele);
  printf("\n Total D*->D0pi: %d; selected %d\n",nTotDstar,nTotDStarSele);
  printf("\n Total Charm->3Prong: %d; selected %d\n",nTot3Prong,nTot3ProngSele);
  if (arrayCascades) printf("\n Total Cascades: %d\n",nTotCasc);

  /*
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
  */

  t.Stop();
  t.Print();

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
