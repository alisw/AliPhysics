void AliVertexerTracksTest(TString outname="AliVertexerTracksTest.root",
			   Bool_t tpconly=kTRUE,
			   Bool_t useMeanVtx=kFALSE,
			   Double_t nSigma=3.,
			   Int_t itsrefit=1,
			   Double_t maxd0z0=0.5,
			   Int_t minitscls=5,
			   Int_t mintracks=1,
			   Bool_t onlyfit=kFALSE) {
//------------------------------------------------------------------------
// Run vertex reconstruction and store results in histos and ntuple
//------------------------------------------------------------------------
// WITHOUT Kinematics.root
//
  AliGeomManager::LoadGeometry("geometry.root");

  if (gAlice) {
    delete AliRunLoader::GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval) {
    cerr<<"LoadgAlice returned error"<<endl;
    delete rl;
    return;
  }
  retval = rl->LoadHeader();
  if (retval) {
    cerr<<"LoadHeader returned error"<<endl;
    delete rl;
    return;
  }
  gAlice=rl->GetAliRun();

  //rl->LoadKinematics();

  // Get field from galice.root
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);
 
  //**** Switch on the PID class
  AliPID pid;

  Int_t nEvents = (Int_t)rl->GetNumberOfEvents();

  Double_t truevtx[3];
  Double_t vtx[3]; 
  Double_t errvtx[3]; 
  Double_t diffX,diffY,diffZ;
  Double_t diffXerr,diffYerr,diffZerr;
  Float_t multiplicity;
  Float_t chi2red;
  Int_t ntracklets;
  Int_t nc;
  Int_t triggered;

  // Histograms
  TH1F *hdiffX = new TH1F("hdiffX","x_{found} - x_{true} [#mum]",1000,-5000,5000);
  TH1F *hdiffY = new TH1F("hdiffY","y_{found} - y_{true} [#mum]",1000,-5000,5000);
  TH1F *hdiffZ = new TH1F("hdiffZ","z_{found} - z_{true} [#mum]",1000,-5000,5000);

  TH1F *hdiffXerr = new TH1F("hdiffXerr","#Delta x/#sigma_{x}",100,-5,5);
  TH1F *hdiffYerr = new TH1F("hdiffYerr","#Delta y/#sigma_{y}",100,-5,5);
  TH1F *hdiffZerr = new TH1F("hdiffZerr","#Delta z/#sigma_{z}",100,-5,5);

  //TNtple
  TNtuple *ntVtxRes = new TNtuple("ntVtxRes","Vertex Resolution's Ntupla with multiplicity","triggered:ntracks:nitstracks:nitstracks5or6:diffX:diffY:diffZ:diffXerr:diffYerr:diffZerr:multiplicity:nc:zMC:chi2red");



  // -----------   Create vertexer --------------------------
  AliESDVertex *initVertex = 0;
  if(useMeanVtx) {
    // open the mean vertex
    TFile *invtx = new TFile("AliESDVertexMean.root");
    initVertex = (AliESDVertex*)invtx->Get("vtxmean");
    invtx->Close();
    delete invtx;
  } else {
    Double_t pos[3]={0.5,0.5,0.}; 
    Double_t err[3]={3.,3.,5.};
    initVertex = new AliESDVertex(pos,err);
  }
  AliVertexerTracks *vertexer = new AliVertexerTracks(AliTracker::GetBz());
  vertexer->SetDebug(1); // set to 1 to see what it does
  vertexer->SetVtxStart(initVertex);
  if(!useMeanVtx) vertexer->SetConstraintOff();
  if(onlyfit) vertexer->SetOnlyFitter();
  vertexer->SetNSigmad0(nSigma);
  if(!itsrefit || tpconly) vertexer->SetITSrefitNotRequired();
  vertexer->SetMinTracks(mintracks);
  vertexer->SetMinITSClusters(minitscls);
  vertexer->SetMaxd0z0(maxd0z0);
  //cout<<" Nsigma:  "<<vertexer->GetNSigmad0()<<endl;
  //----------------------------------------------------------
 

  if(gSystem->AccessPathName("AliESDs_itstpc.root",kFileExists)) {
    printf("AliESDs_itstpc.root not found!\n"); 
    return;
  }
  TFile *fin = TFile::Open("AliESDs_itstpc.root");
  AliESDEvent *esd = new AliESDEvent();
  TTree *esdTree = (TTree*)fin->Get("esdTree");
  if(!esdTree) return;
  esd->ReadFromTree(esdTree);
  Int_t events = esdTree->GetEntries();
  printf("Number of events in ESD tree: %d\n",events);

  TArrayF o;

  for(Int_t iev=0; iev<events; iev++) { //LOOP ON EVENTS
    printf("-------------- EVENT   %d --------------------\n",iev);

    diffX=999999;
    diffY=999999;
    diffZ=999999;
    diffXerr=999999;
    diffYerr=999999;
    diffZerr=999999;
    nc=0;

    //=================================================
    //         LOOK IN THE SIMULATION
    // true vertex position
    Int_t npart = (Int_t)gAlice->GetEvent(iev);
    printf("particles  %d\n",npart);
    AliHeader *header = (AliHeader*)rl->GetHeader();
    AliGenPythiaEventHeader *pyheader = (AliGenPythiaEventHeader*)header->GenEventHeader();
    pyheader->PrimaryVertex(o);
    truevtx[0] = o[0];
    truevtx[1] = o[1];
    truevtx[2] = o[2];
    printf("True Vertex: (%f, %f, %f)\n",o[0],o[1],o[2]);

    Double_t dNchdy = 0.;
    multiplicity=0.;
    Int_t nitstracks = 0;
    Int_t nitstracks1 = 0;
    Int_t nitstracks2 = 0;
    Int_t nitstracks3 = 0;
    Int_t nitstracks4 = 0;
    Int_t nitstracks5 = 0;
    Int_t nitstracks6 = 0;
    Int_t nitstracks5or6 = 0;


    esdTree->GetEvent(iev);
    triggered=0;
    chi2red=0.;
    ULong64_t triggerMask = esd->GetTriggerMask();
    if(triggerMask&32 && ((triggerMask&1) || (triggerMask&2))) triggered=1;

    // get # of tracklets
    AliMultiplicity *alimult = esd->GetMultiplicity();
    if(alimult) {
      ntracklets = alimult->GetNumberOfTracklets();
      for(Int_t l=0;l<alimult->GetNumberOfTracklets();l++) 
	if(alimult->GetDeltaPhi(l)<-9998.) ntracklets--;
    } else {
      ntracklets = 0;
    }
    printf("Number of SPD tracklets in ESD %d  :  %d\n",iev,ntracklets);
    multiplicity = (Float_t)ntracklets;

    Int_t ntracks = esd->GetNumberOfTracks();
    printf("Number of tracks in ESD %d   :   %d\n",iev,ntracks);
    if(ntracks==0) { 
      ntVtxRes->Fill(triggered,ntracks,nitstracks,nitstracks5or6,diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr,multiplicity,nc,truevtx[2],chi2red); 
      continue; 
    }
    
    // VERTEX    
    AliESDVertex *vertex = 0;
    if(!tpconly) {
      vertex = (AliESDVertex*)vertexer->FindPrimaryVertex(esd);
    } else {
      TObjArray trkArrayOrig(ntracks);
      UShort_t *idOrig = new UShort_t[ntracks];
      const Double_t kRadius  = 2.8; //something less than the beam pipe radius
      const Double_t kMaxStep = 5;   //max step over the material
      Bool_t ok;
      Int_t nTrksOrig=0;
      for(Int_t i=0; i<ntracks; i++) {
	AliESDtrack *esdt = esd->GetTrack(i);
	AliExternalTrackParam *tpcTrack =
	  (AliExternalTrackParam *)esdt->GetTPCInnerParam();
	ok = kFALSE;
	if (tpcTrack)
	  ok = AliTracker::
	    PropagateTrackTo(tpcTrack,kRadius,esdt->GetMass(),kRadius,kTRUE);
	if (ok) {
	  trkArrayOrig.AddLast(tpcTrack);
	  idOrig[nTrksOrig]=(UShort_t)esdt->GetID();
	  nTrksOrig++;
	}
      }
      
      vertex = (AliESDVertex*)vertexer->FindPrimaryVertex(&trkArrayOrig,idOrig);
      delete [] idOrig;
    }


    nc = (Int_t)vertex->GetNContributors();
    if(nc>=2) chi2red = vertex->GetChi2toNDF(); 
    printf("nc = %d\n",nc);


    // count ITS tracks
    for(Int_t itrk=0; itrk<ntracks; itrk++) {
      AliESDtrack *esdtrack = (AliESDtrack*)esd->GetTrack(itrk);
      UInt_t status = esdtrack->GetStatus();
      // only tracks found also in ITS
      if(! (status&AliESDtrack::kITSin) ) continue;      
      nitstracks++;
      Int_t npts = (Int_t)esdtrack->GetNcls(0);
      if(npts==6) {nitstracks6++;nitstracks5or6++;}
      if(npts==5) {nitstracks5++;nitstracks5or6++;}
      if(npts==4) nitstracks4++;
      if(npts==3) nitstracks3++;
      if(npts==2) nitstracks2++;
      if(npts==1) nitstracks1++;
    }
    printf("Number of kITSin tracks in ESD %d   :   %d\n",iev,nitstracks);
    printf("           6 pts: %d\n",nitstracks6);
    printf("           5 pts: %d\n",nitstracks5);
    printf("           4 pts: %d\n",nitstracks4);
    printf("           3 pts: %d\n",nitstracks3);
    printf("           2 pts: %d\n",nitstracks2);
    printf("           1 pts: %d\n",nitstracks1);


    if(nc>=1) {
      vertex->GetXYZ(vtx);
      vertex->GetSigmaXYZ(errvtx); 
      diffX = 10000.* (vtx[0] - truevtx[0]);
      diffY = 10000.* (vtx[1] - truevtx[1]);
      diffZ = 10000.* (vtx[2] - truevtx[2]);
      hdiffX->Fill(diffX);
      hdiffY->Fill(diffY);
      hdiffZ->Fill(diffZ);
      diffXerr = diffX/errvtx[0]/10000.;
      diffYerr = diffY/errvtx[1]/10000.;
      diffZerr = diffZ/errvtx[2]/10000.;
      hdiffXerr->Fill(diffXerr);
      hdiffYerr->Fill(diffYerr);
      hdiffZerr->Fill(diffZerr);
    } 

    ntVtxRes->Fill(triggered,ntracks,nitstracks,nitstracks5or6,diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr,multiplicity,(Float_t)nc,truevtx[2],chi2red);   

  } // END LOOP ON EVENTS


  delete esdTree;
  delete fin;
  delete vertexer;
  delete initVertex;

  // Output file
  TFile *fout  = new TFile(outname.Data(),"recreate");
  ntVtxRes->Write();
  hdiffX->Write();
  hdiffY->Write();
  hdiffZ->Write();
  hdiffXerr->Write();
  hdiffYerr->Write();
  hdiffZerr->Write();
  fout->Close();

  return;
} 
//------------------------------------------------------------------------
void VertexForOneEvent(Int_t iev=0,
		       Double_t nSigma=3.,
		       Bool_t useMeanVtx=kFALSE) {
//------------------------------------------------------------------------
// Run vertex reconstruction for event iev
//------------------------------------------------------------------------

  if (gAlice) {
    delete AliRunLoader::GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval) {
    cerr<<"LoadgAlice returned error"<<endl;
    delete rl;
    return;
  }
  retval = rl->LoadHeader();
  if (retval) {
    cerr<<"LoadHeader returned error"<<endl;
    delete rl;
    return;
  }
  gAlice=rl->GetAliRun();
  //rl->LoadKinematics();
  // Get field from galice.root
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);
 
  Double_t truevtx[3];
  Double_t vtx[3]; 
  Double_t errvtx[3]; 
  Double_t diffX,diffY,diffZ;
  Double_t diffXerr,diffYerr,diffZerr;
  Double_t multiplicity;
  Int_t nc;

  // -----------   Create vertexer --------------------------
  AliESDVertex *initVertex = 0;
  if(useMeanVtx) {
    // open the mean vertex
    TFile *invtx = new TFile("AliESDVertexMean.root");
    initVertex = (AliESDVertex*)invtx->Get("vtxmean");
    invtx->Close();
    delete invtx;
  } else {
    Double_t pos[3]={0.,-0.,0.}; 
    Double_t err[3]={3.,3.,5.};
    initVertex = new AliESDVertex(pos,err);
  }
  AliVertexerTracks *vertexer = new AliVertexerTracks(AliTracker::GetBz());
  //vertexer->SetITSrefitNotRequired();
  vertexer->SetVtxStart(initVertex);
  if(!useMeanVtx) vertexer->SetConstraintOff();
  //vertexer->SetMinTracks(2);
  //vertexer->SetNSigmad0(nSigma);
  //cout<<" Nsigma:  "<<vertexer->GetNSigmad0()<<endl;
  //vertexer->SetFinderAlgorithm(2);
  //vertexer->SetDCAcut(0.1); // 1 mm
  vertexer->SetDebug(1); // set to 1 to see what it does
  //----------------------------------------------------------
 
  if(gSystem->AccessPathName("AliESDs.root",kFileExists)) {
    printf("AliESDs.root not found!\n"); 
    return;
  }
  TFile *fin = TFile::Open("AliESDs.root");
  AliESDEvent *esd = new AliESDEvent();
  TTree *esdTree = (TTree*)fin->Get("esdTree");
  if(!esdTree) return;
  esd->ReadFromTree(esdTree);

  TArrayF o;

  nc=0;

  //=================================================
  //         LOOK IN THE SIMULATION
  // true vertex position
  Int_t npart = (Int_t)gAlice->GetEvent(iev);
  printf("particles  %d\n",npart);
  AliHeader *header = (AliHeader*)rl->GetHeader();
  AliGenPythiaEventHeader *pyheader = (AliGenPythiaEventHeader*)header->GenEventHeader();
  pyheader->PrimaryVertex(o);
  truevtx[0] = o[0];
  truevtx[1] = o[1];
  truevtx[2] = o[2];
  printf("True Vertex: (%f, %f, %f)\n",o[0],o[1],o[2]);

  esdTree->GetEvent(iev);
  Int_t ntracks = esd->GetNumberOfTracks();
  printf("Number of tracks in ESD %d   :   %d\n",iev,ntracks);
    
  // VERTEX    
  AliESDVertex *vertex = (AliESDVertex*)vertexer->FindPrimaryVertex(esd);
  nc = (Int_t)vertex->GetNContributors();


  Int_t nitstracks = 0;
  Int_t nitstracks1 = 0;
  Int_t nitstracks2 = 0;
  Int_t nitstracks3 = 0;
  Int_t nitstracks4 = 0;
  Int_t nitstracks5 = 0;
  Int_t nitstracks6 = 0;
  Int_t nitstracks5or6 = 0;

  // count ITS tracks
  for(Int_t itrk=0; itrk<ntracks; itrk++) {
    AliESDtrack *esdtrack = (AliESDtrack*)esd->GetTrack(itrk);
    UInt_t status = esdtrack->GetStatus();
    // only tracks found also in ITS
    if(! (status&AliESDtrack::kITSrefit) ) continue;      
    nitstracks++;
    Int_t npts = (Int_t)esdtrack->GetNcls(0);
    if(npts==6) {nitstracks6++;nitstracks5or6++;}
    if(npts==5) {nitstracks5++;nitstracks5or6++;}
    if(npts==4) nitstracks4++;
    if(npts==3) nitstracks3++;
    if(npts==2) nitstracks2++;
    if(npts==1) nitstracks1++;
  }
  printf("Number of kITSrefit tracks in ESD %d   :   %d\n",iev,nitstracks);
  printf("           6 pts: %d\n",nitstracks6);
  printf("           5 pts: %d\n",nitstracks5);
  printf("           4 pts: %d\n",nitstracks4);
  printf("           3 pts: %d\n",nitstracks3);
  printf("           2 pts: %d\n",nitstracks2);
  printf("           1 pts: %d\n",nitstracks1);
  

  delete esdTree;
  delete fin;
  delete vertexer;
  delete initVertex;

  return;
} 
//----------------------------------------------------------------------------
Double_t FitFunc(Double_t *x,Double_t *par) {
  return par[0]+par[1]/TMath::Sqrt(x[0]);
}
Int_t GetBin(Float_t mult) {
  if(mult<5.)  return 0;
  if(mult<7.)  return 1;
  if(mult<12.) return 2;
  if(mult<15.) return 3;
  if(mult<22.) return 4;
  return 5;
}
Int_t GetBinTracklets(Float_t tracklets) {
  if(tracklets<2.*5.)  return 0;
  if(tracklets<2.*7.)  return 1;
  if(tracklets<2.*12.) return 2;
  if(tracklets<2.*15.) return 3;
  if(tracklets<2.*22.) return 4;
  return 5;
}
Int_t GetBinZ(Float_t z) {
  if(z<-13.) return 0;
  if(z<-11.) return 1;
  if(z<-9.) return 2;
  if(z<-7.) return 3;
  if(z<-5.) return 4;
  if(z<-3.) return 5;
  if(z<-1.) return 6;
  if(z<1.) return 7;
  if(z<3.) return 8;
  if(z<5.) return 9;
  if(z<7.) return 10;
  if(z<9.) return 11;
  if(z<11.) return 12;
  if(z<13.) return 13;
  return 14;
}
//----------------------------------------------------------------------------
void PlotVtxRes(TString inname="AliVertexerTracksTest.root",
		Bool_t tpconly=kFALSE) {
  //---------------------------------------------------------------------
  // Plot efficiency, resolutions, pulls 
  // [at least 10k pp events are needed]
  //---------------------------------------------------------------------

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(10001);

  Float_t zMCmin=0.;
  Float_t zMCmax=15.;
  Float_t ncminforzshift=1.;
  Float_t maxdiffx = 1e3;
  Float_t rangefactor=1.;
  if(tpconly) {
    maxdiffx *= 1e2;
    rangefactor = 15.;
  }

  Int_t   nEvVtx=0;
  Int_t   nEvNoVtx=0;
  Int_t   ev,nUsedTrks;
  Float_t nESDTrks,nTrks,nTrks5or6,ntracklets;
  Float_t diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr,zMC,nc;
  Float_t triggered;


  //
  // HISTOGRAMS
  //
  TH1F *hdx = new TH1F("hdx","",50,-600*rangefactor,600*rangefactor);
  hdx->SetXTitle("#Delta x [#mu m]");
  hdx->SetYTitle("events");
  //
  TH1F *hdy = new TH1F("hdy","",50,-600*rangefactor,600*rangefactor);
  hdy->SetXTitle("#Delta y [#mu m]");
  hdy->SetYTitle("events");
  //
  TH1F *hdz = new TH1F("hdz","",200,-600*rangefactor,600*rangefactor);
  hdz->SetXTitle("#Delta z [#mu m]");
  hdz->SetYTitle("events");
  //
  TH1F *hzMCVtx = new TH1F("hzMCVtx","events w/ vertex",30,-15,15);
  hzMCVtx->SetXTitle("z_{true} [cm]");
  hzMCVtx->SetYTitle("events");
  hzMCVtx->Sumw2();
  //
  TH1F *hzMCNoVtx = new TH1F("hzMCNoVtx","events w/o vertex",30,-15,15);
  hzMCNoVtx->SetXTitle("z_{true} [cm]");
  hzMCNoVtx->SetYTitle("events");
  hzMCNoVtx->Sumw2();
  //
  TH1F *hTrackletsVtx = new TH1F("hTrackletsVtx","events w/ vertex",51,-0.5,50.5);
  hTrackletsVtx->SetXTitle("SPD tracklets");
  hTrackletsVtx->SetYTitle("events");
  //
  TH1F *hTrackletsNoVtx = new TH1F("hTrackletsNoVtx","events w/o vertex",51,-0.5,50.5);
  hTrackletsNoVtx->SetXTitle("SPD tracklets");
  hTrackletsNoVtx->SetYTitle("events");
  //
  TH1F *hTracksVtx = new TH1F("hTracksVtx","events w/ vertex",51,-0.5,50.5);
  hTracksVtx->SetXTitle("Number of reco tracks in TPC&ITS");
  hTracksVtx->SetYTitle("events");
  //
  TH1F *hTracksNoVtx = new TH1F("hTracksNoVtx","events w/o vertex",51,-0.5,50.5);
  hTracksNoVtx->SetXTitle("Number of reco tracks in TPC&ITS");
  hTracksNoVtx->SetYTitle("events");
  //
  TH1F *hTracksVtx5or6 = new TH1F("hTracksVtx5or6","events w/ vertex",51,-0.5,50.5);
  hTracksVtx5or6->SetXTitle("Number of reco tracks in TPC&ITS(cls>=5)");
  hTracksVtx5or6->SetYTitle("events");
  //
  TH1F *hTracksNoVtx5or6 = new TH1F("hTracksNoVtx5or6","events w/o vertex",51,-0.5,50.5);
  hTracksNoVtx5or6->SetXTitle("Number of reco tracks in TPC&ITS(cls>=5)");
  hTracksNoVtx5or6->SetYTitle("events");
  //
  TH1F *hPullsx = new TH1F("hPullsx","",50,-7.,7.);
  hPullsx->SetXTitle("#Delta x/#sqrt{<x,x>}");
  hPullsx->SetYTitle("events");
  //
  TH1F *hPullsy = new TH1F("hPullsy","",50,-7.,7.);
  hPullsy->SetXTitle("#Delta y/#sqrt{<y,y>}");
  hPullsy->SetYTitle("events");
  //
  TH1F *hPullsz = new TH1F("hPullsz","",50,-7.,7.);
  hPullsz->SetXTitle("#Delta z/#sqrt{<z,z>}");
  hPullsz->SetYTitle("events");
  //
  TProfile *hntrackletsVSzMC = new TProfile("hntrackletsVSzMC","",30,-15,15,0,30);
  hntrackletsVSzMC->SetXTitle("z_{true} [cm]");
  hntrackletsVSzMC->SetYTitle("<SPD tracklets>");
  //
  TProfile *hnitstracksVSzMC = new TProfile("hnitstracksVSzMC","",30,-15,15,0,30);
  hnitstracksVSzMC->SetXTitle("z_{true} [cm]");
  hnitstracksVSzMC->SetYTitle("<tracks TPC&ITS>");
  //
  TProfile *hnitstracks5or6VSzMC = new TProfile("hnitstracks5or6VSzMC","",30,-15,15,0,30);
  hnitstracks5or6VSzMC->SetXTitle("z_{true} [cm]");
  hnitstracks5or6VSzMC->SetYTitle("<tracks TPC&ITS(cls>=5)>");

  Double_t mult[6]={0,0,0,0,0,0};
  Double_t mult2[6]={0,0,0,0,0,0};
  Double_t sigmamult[6]={0,0,0,0,0,0};
  Double_t sigmamult2[6]={0,0,0,0,0,0};
  Double_t tottracks[6]={0,0,0,0,0,0};
  Double_t fittrks[6]={0,0,0,0,0,0};
  Double_t tracks[6]={0,0,0,0,0,0};
  Double_t etracks[6]={0,0,0,0,0,0};
  Double_t xres[6]={0,0,0,0,0,0};
  Double_t exres[6]={0,0,0,0,0,0};
  Double_t yres[6]={0,0,0,0,0,0};
  Double_t eyres[6]={0,0,0,0,0,0};
  Double_t zres[6]={0,0,0,0,0,0};
  Double_t ezres[6]={0,0,0,0,0,0};
  Double_t xpull[6]={0,0,0,0,0,0};
  Double_t expull[6]={0,0,0,0,0,0};
  Double_t ypull[6]={0,0,0,0,0,0};
  Double_t eypull[6]={0,0,0,0,0,0};
  Double_t zpull[6]={0,0,0,0,0,0};
  Double_t ezpull[6]={0,0,0,0,0,0};
  Int_t    evts[6]={0,0,0,0,0,0};
  Int_t    totevts[6]={0,0,0,0,0,0};
  Int_t    vtx[6]={0,0,0,0,0,0};
  Int_t    all[6]={0,0,0,0,0,0};
  Double_t probVtx[6]={0,0,0,0,0,0};
  Double_t eprobVtx[6]={0,0,0,0,0,0};
  Int_t    bin,zbin;
  //
  // OTHER HISTOGRAMS
  //
  TH1F *hdx0 = new TH1F("hdx0","x resolution - bin 0",50,-500*rangefactor,500*rangefactor);
  TH1F *hdx1 = new TH1F("hdx1","x resolution - bin 1",50,-500*rangefactor,500*rangefactor);
  TH1F *hdx2 = new TH1F("hdx2","x resolution - bin 2",50,-500*rangefactor,500*rangefactor);
  TH1F *hdx3 = new TH1F("hdx3","x resolution - bin 3",50,-400*rangefactor,400*rangefactor);
  TH1F *hdx4 = new TH1F("hdx4","x resolution - bin 4",50,-300*rangefactor,300*rangefactor);
  TH1F *hdx5 = new TH1F("hdx5","x resolution - bin 5",50,-300*rangefactor,300*rangefactor);
  TH1F *hdy0 = new TH1F("hdy0","y resolution - bin 0",50,-500*rangefactor,500*rangefactor);
  TH1F *hdy1 = new TH1F("hdy1","y resolution - bin 1",50,-500*rangefactor,500*rangefactor);
  TH1F *hdy2 = new TH1F("hdy2","y resolution - bin 2",50,-500*rangefactor,500*rangefactor);
  TH1F *hdy3 = new TH1F("hdy3","y resolution - bin 3",50,-400*rangefactor,400*rangefactor);
  TH1F *hdy4 = new TH1F("hdy4","y resolution - bin 4",50,-300*rangefactor,300*rangefactor);
  TH1F *hdy5 = new TH1F("hdy5","y resolution - bin 5",50,-300*rangefactor,300*rangefactor);
  TH1F *hdz0 = new TH1F("hdz0","z resolution - bin 0",50,-1000*rangefactor,1000*rangefactor);
  TH1F *hdz1 = new TH1F("hdz1","z resolution - bin 1",50,-800*rangefactor,800*rangefactor);
  TH1F *hdz2 = new TH1F("hdz2","z resolution - bin 2",50,-800*rangefactor,800*rangefactor);
  TH1F *hdz3 = new TH1F("hdz3","z resolution - bin 3",50,-600*rangefactor,600*rangefactor);
  TH1F *hdz4 = new TH1F("hdz4","z resolution - bin 4",50,-500*rangefactor,500*rangefactor);
  TH1F *hdz5 = new TH1F("hdz5","z resolution - bin 5",50,-500*rangefactor,500*rangefactor);
  //
  TH1F *hdz_z = NULL; hdz_z = new TH1F[15];
  for(Int_t l=0;l<15;l++) hdz_z[l] = *hdz;


  TH1F *hpx0 = new TH1F("hpx0","x pull - bin 0",50,-7,7);
  TH1F *hpx1 = new TH1F("hpx1","x pull - bin 1",50,-7,7);
  TH1F *hpx2 = new TH1F("hpx2","x pull - bin 2",50,-7,7);
  TH1F *hpx3 = new TH1F("hpx3","x pull - bin 3",50,-7,7);
  TH1F *hpx4 = new TH1F("hpx4","x pull - bin 4",50,-7,7);
  TH1F *hpx5 = new TH1F("hpx5","x pull - bin 5",50,-7,7);
  TH1F *hpy0 = new TH1F("hpy0","y pull - bin 0",50,-7,7);
  TH1F *hpy1 = new TH1F("hpy1","y pull - bin 1",50,-7,7);
  TH1F *hpy2 = new TH1F("hpy2","y pull - bin 2",50,-7,7);
  TH1F *hpy3 = new TH1F("hpy3","y pull - bin 3",50,-7,7);
  TH1F *hpy4 = new TH1F("hpy4","y pull - bin 4",50,-7,7);
  TH1F *hpy5 = new TH1F("hpy5","y pull - bin 5",50,-7,7);
  TH1F *hpz0 = new TH1F("hpz0","z pull - bin 0",50,-7,7);
  TH1F *hpz1 = new TH1F("hpz1","z pull - bin 1",50,-7,7);
  TH1F *hpz2 = new TH1F("hpz2","z pull - bin 2",50,-7,7);
  TH1F *hpz3 = new TH1F("hpz3","z pull - bin 3",50,-7,7);
  TH1F *hpz4 = new TH1F("hpz4","z pull - bin 4",50,-7,7);
  TH1F *hpz5 = new TH1F("hpz5","z pull - bin 5",50,-7,7);


  TH1F *hmult = new TH1F("hmult","hmult",50,-0.5,49.5);

  TFile *in = new TFile(inname.Data());
  TNtuple *nt = (TNtuple*)in->Get("ntVtxRes");
  nt->SetBranchAddress("triggered",&triggered);
  nt->SetBranchAddress("ntracks",&nESDTrks);
  nt->SetBranchAddress("nitstracks",&nTrks);
  nt->SetBranchAddress("nitstracks5or6",&nTrks5or6);
  nt->SetBranchAddress("diffX",&diffX);
  nt->SetBranchAddress("diffY",&diffY);
  nt->SetBranchAddress("diffZ",&diffZ);
  nt->SetBranchAddress("diffXerr",&diffXerr);
  nt->SetBranchAddress("diffYerr",&diffYerr);
  nt->SetBranchAddress("diffZerr",&diffZerr);
  nt->SetBranchAddress("multiplicity",&ntracklets);
  nt->SetBranchAddress("zMC",&zMC);
  nt->SetBranchAddress("nc",&nc);
  Int_t entries = (Int_t)nt->GetEntries();
  Int_t nbytes=0;

  for(Int_t i=0; i<entries; i++) {
    nbytes += nt->GetEvent(i);
    // only triggered events
    //    if(triggered<0.5) continue;
    // consider only events with zMCmin<|zMC|<zMCmax
    if(TMath::Abs(zMC)<zMCmin || TMath::Abs(zMC)>zMCmax) continue;
    //
    bin = GetBinTracklets(ntracklets);
    zbin = GetBinZ(zMC);
    mult[bin] += ntracklets;
    hmult->Fill(ntracklets);
    totevts[bin]++;
    hntrackletsVSzMC->Fill(zMC,ntracklets);
    hnitstracksVSzMC->Fill(zMC,nTrks);
    hnitstracks5or6VSzMC->Fill(zMC,nTrks5or6);
    if(TMath::Abs(diffX) < maxdiffx) { // vtx found - closer than maxdiffx micron
      nEvVtx++;
      mult2[bin] += ntracklets;
      vtx[bin]++;
      tottracks[bin] += nTrks;
      evts[bin]++;
      if(bin==0) hdx0->Fill(diffX);
      if(bin==1) hdx1->Fill(diffX);
      if(bin==2) hdx2->Fill(diffX);
      if(bin==3) hdx3->Fill(diffX);
      if(bin==4) hdx4->Fill(diffX);
      if(bin==5) hdx5->Fill(diffX);
      if(bin==0) hdy0->Fill(diffY);
      if(bin==1) hdy1->Fill(diffY);
      if(bin==2) hdy2->Fill(diffY);
      if(bin==3) hdy3->Fill(diffY);
      if(bin==4) hdy4->Fill(diffY);
      if(bin==5) hdy5->Fill(diffY);
      if(bin==0) hdz0->Fill(diffZ);
      if(bin==1) hdz1->Fill(diffZ);
      if(bin==2) hdz2->Fill(diffZ);
      if(bin==3) hdz3->Fill(diffZ);
      if(bin==4) hdz4->Fill(diffZ);
      if(bin==5) hdz5->Fill(diffZ);
      hdx->Fill(diffX);
      hdy->Fill(diffY);
      hdz->Fill(diffZ);
      if(nc>=ncminforzshift) hdz_z[zbin].Fill(diffZ);
      hzMCVtx->Fill(zMC);
      hTrackletsVtx->Fill(ntracklets);
      hTracksVtx->Fill(nTrks);
      hTracksVtx5or6->Fill(nTrks5or6);
      hPullsx->Fill(diffXerr);
      hPullsy->Fill(diffYerr);
      hPullsz->Fill(diffZerr);
      if(bin==0) hpx0->Fill(diffXerr);
      if(bin==1) hpx1->Fill(diffXerr);
      if(bin==2) hpx2->Fill(diffXerr);
      if(bin==3) hpx3->Fill(diffXerr);
      if(bin==4) hpx4->Fill(diffXerr);
      if(bin==5) hpx5->Fill(diffXerr);
      if(bin==0) hpy0->Fill(diffYerr);
      if(bin==1) hpy1->Fill(diffYerr);
      if(bin==2) hpy2->Fill(diffYerr);
      if(bin==3) hpy3->Fill(diffYerr);
      if(bin==4) hpy4->Fill(diffYerr);
      if(bin==5) hpy5->Fill(diffYerr);
      if(bin==0) hpz0->Fill(diffZerr);
      if(bin==1) hpz1->Fill(diffZerr);
      if(bin==2) hpz2->Fill(diffZerr);
      if(bin==3) hpz3->Fill(diffZerr);
      if(bin==4) hpz4->Fill(diffZerr);
      if(bin==5) hpz5->Fill(diffZerr);
    } else {
      nEvNoVtx++;
      hzMCNoVtx->Fill(zMC);
      hTrackletsNoVtx->Fill(ntracklets);
      hTracksNoVtx->Fill(nTrks);
      hTracksNoVtx5or6->Fill(nTrks5or6);
    }
    all[bin]++;
  }


  for(Int_t k=0;k<6;k++) {
    mult2[k]      /= evts[k];
    mult[k]      /= totevts[k];
    tottracks[k] /= evts[k];
    probVtx[k]    = (Double_t)vtx[k]/all[k];
    eprobVtx[k]   = (Double_t)vtx[k]/all[k]/all[k]+(Double_t)vtx[k]*vtx[k]/all[k]/all[k]/all[k];
    eprobVtx[k] = TMath::Sqrt(eprobVtx[k]); 
  }

  // calculate spread in multiplicity
  for(Int_t i=0; i<entries; i++) {
    nbytes += nt->GetEvent(i);
    bin = GetBinTracklets(ntracklets);
    sigmamult[bin] += (ntracklets-mult[bin])*(ntracklets-mult[bin])/totevts[bin];
    if(diffX < 90000.) sigmamult2[bin] += (ntracklets-mult2[bin])*(ntracklets-mult2[bin])/evts[bin];
  }

  for(Int_t k=0;k<6;k++) {
    sigmamult[k]  = TMath::Sqrt(sigmamult[k]);
    sigmamult2[k] = TMath::Sqrt(sigmamult2[k]);
  }

  // fraction of found vertices
  Float_t frac = (Float_t)nEvVtx/(nEvVtx+nEvNoVtx);
  printf(" Fraction of events with vertex: %f\n",frac);


  // FIT FUNCTIONS
  TF1 *g = new TF1("g","gaus");
  g->SetLineWidth(3);
  TF1 *line = new TF1("line","pol1");
  line->SetLineWidth(3);
  TF1 *func = new TF1("func",FitFunc,2,80,2);
  func->SetLineWidth(1);

  Double_t x1,y1,sigmafit;
  Char_t message[100];

  //
  //            DRAW PLOTS
  //
  gStyle->SetOptFit(1111);
  
  // tracks VS zMC
  TCanvas *c0 = new TCanvas("c0","c0",0,0,1000,400);
  c0->Divide(3,1);
  c0->cd(1); 
  hntrackletsVSzMC->SetMinimum(0);
  hntrackletsVSzMC->SetMaximum(14);
  hntrackletsVSzMC->Draw("profs");
  c0->cd(2); 
  hnitstracksVSzMC->SetMinimum(0);
  hnitstracksVSzMC->SetMaximum(14);
  hnitstracksVSzMC->Draw("profs");
  c0->cd(3); 
  hnitstracks5or6VSzMC->SetMinimum(0);
  hnitstracks5or6VSzMC->SetMaximum(14);
  hnitstracks5or6VSzMC->Draw("profs");


  // Tracks in ITS for events w/ and w/o vertex
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,500);
  c1->Divide(2,1);
  c1->cd(1); 
  hTracksVtx->Draw();
  c1->cd(2); 
  hTracksNoVtx->Draw();

  // probability vs ntracklets
  TCanvas *c1a = new TCanvas("c1a","c1a",0,0,500,500);
  TH1F *hTotTracklets = (TH1F*)hTrackletsVtx->Clone("hTotTracklets");
  hTotTracklets->Add(hTrackletsNoVtx);
  TH1F *hProbTracklets = (TH1F*)hTrackletsVtx->Clone("hProbTracklets");
  hProbTracklets->Divide(hTotTracklets);
  hProbTracklets->SetYTitle("Probability");
  hProbTracklets->SetTitle("Probability to find the vertex");
  hProbTracklets->Draw();

  // probability vs ntracks
  TCanvas *c1b = new TCanvas("c1b","c1b",0,0,500,500);
  TH1F *hTotTracks = (TH1F*)hTracksVtx->Clone("hTotTracks");
  hTotTracks->Add(hTracksNoVtx);
  TH1F *hProbTracks = (TH1F*)hTracksVtx->Clone("hProbTracks");
  hProbTracks->Divide(hTotTracks);
  hProbTracks->SetYTitle("Probability");
  hProbTracks->SetTitle("Probability to find the vertex");
  hProbTracks->Draw();

  TCanvas *c1c = new TCanvas("c1c","c1c",0,0,500,500);
  TH1F *hTotTracks5or6 = (TH1F*)hTracksVtx5or6->Clone("hTotTracks5or6");
  hTotTracks5or6->Add(hTracksNoVtx5or6);
  TH1F *hProbTracks5or6 = (TH1F*)hTracksVtx5or6->Clone("hProbTracks5or6");
  hProbTracks5or6->Divide(hTotTracks5or6);
  hProbTracks5or6->SetYTitle("Probability");
  hProbTracks5or6->SetTitle("Probability to find the vertex");
  hProbTracks5or6->Draw();


  // probability vs zMC
  TCanvas *c1e = new TCanvas("c1e","c1e",0,0,500,500);
  TH1F *hTotzMC = (TH1F*)hzMCVtx->Clone("hTotzMC");
  hTotzMC->Add(hzMCNoVtx);
  TH1F *hProbzMC = (TH1F*)hzMCVtx->Clone("hProbzMC");
  hProbzMC->Divide(hTotzMC);
  hProbzMC->SetYTitle("Probability");
  hProbzMC->SetTitle("Probability to find the vertex");
  hProbzMC->Draw();

  // Global resolutions
  TCanvas *c2 = new TCanvas("c2","c2",0,0,1000,400);
  c2->Divide(3,1);
  c2->cd(1); 
  hdx->Draw();
  g->SetRange(-1.*hdx->GetRMS(),+1.*hdx->GetRMS());
  hdx->Fit("g","R");
  c2->cd(2); 
  hdy->Draw();
  g->SetRange(-1.*hdy->GetRMS(),+1.*hdy->GetRMS());
  hdy->Fit("g","R");
  c2->cd(3); 
  hdz->Draw();
  g->SetRange(-1.*hdz->GetRMS(),+1.*hdz->GetRMS());
  hdz->Fit("g","R");

  // Pulls
  TCanvas *c4 = new TCanvas("c4","c4",0,0,1000,400);
  c4->Divide(3,1);
  c4->cd(1); 
  hPullsx->Draw();
  g->SetRange(-3.*hPullsx->GetRMS(),+3.*hPullsx->GetRMS());
  hPullsx->Fit("g","R");
  c4->cd(2); 
  hPullsy->Draw();
  g->SetRange(-3.*hPullsy->GetRMS(),+3.*hPullsy->GetRMS());
  hPullsy->Fit("g","R");
  c4->cd(3); 
  hPullsz->Draw();
  g->SetRange(-3.*hPullsz->GetRMS(),+3.*hPullsz->GetRMS());
  hPullsz->Fit("g","R");
  
  // probability VS multiplicity
  TCanvas *c5 = new TCanvas("c5","c5",0,0,600,500);
  TH2F *fr5 = new TH2F("fr5","",2,0,80,2,0,1.1); 
  fr5->SetXTitle("SPD tracklets");
  fr5->SetYTitle("Probability");
  fr5->Draw();
  TGraphErrors *gr5 = new TGraphErrors(6,mult,probVtx,sigmamult,eprobVtx);
  gr5->Draw("p");
  gr5->SetMarkerStyle(21);

  //
  // resolutions VS multiplicity
  //
  // X
  TCanvas *ax = new TCanvas("ax","ax",0,0,900,700);
  ax->Divide(3,2);
  ax->cd(1); 
  hdx0->Draw(); 
  g->SetRange(-3.*hdx0->GetRMS(),+3.*hdx0->GetRMS());
  hdx0->Fit("g","Q"); 
  xres[0]=g->GetParameter(2);
  exres[0]=g->GetParError(2);
  ax->cd(2); 
  hdx1->Draw(); 
  g->SetRange(-3.*hdx1->GetRMS(),+3.*hdx1->GetRMS());
  hdx1->Fit("g","Q"); 
  xres[1]=g->GetParameter(2);
  exres[1]=g->GetParError(2);
  ax->cd(3); 
  hdx2->Draw(); 
  g->SetRange(-3.*hdx2->GetRMS(),+3.*hdx2->GetRMS());
  hdx2->Fit("g","Q"); 
  xres[2]=g->GetParameter(2);
  exres[2]=g->GetParError(2);
  ax->cd(4); 
  hdx3->Draw();
  g->SetRange(-3.*hdx3->GetRMS(),+3.*hdx3->GetRMS());
  hdx3->Fit("g","Q"); 
  xres[3]=g->GetParameter(2);
  exres[3]=g->GetParError(2);
  ax->cd(5); 
  hdx4->Draw();
  g->SetRange(-3.*hdx4->GetRMS(),+3.*hdx4->GetRMS());
  hdx4->Fit("g","Q");
  xres[4]=g->GetParameter(2);
  exres[4]=g->GetParError(2);
  ax->cd(6); 
  hdx5->Draw();
  g->SetRange(-3.*hdx5->GetRMS(),+3.*hdx5->GetRMS());
  hdx5->Fit("g","Q");
  xres[5]=g->GetParameter(2);
  exres[5]=g->GetParError(2);
  // Y
  TCanvas *ay = new TCanvas("ay","ay",0,0,900,700);
  ay->Divide(3,2);
  ay->cd(1); 
  hdy0->Draw(); 
  g->SetRange(-3.*hdy0->GetRMS(),+3.*hdy0->GetRMS());
  hdy0->Fit("g","Q"); 
  yres[0]=g->GetParameter(2);
  eyres[0]=g->GetParError(2);
  ay->cd(2); 
  hdy1->Draw(); 
  g->SetRange(-3.*hdy1->GetRMS(),+3.*hdy1->GetRMS());
  hdy1->Fit("g","Q"); 
  yres[1]=g->GetParameter(2);
  eyres[1]=g->GetParError(2);
  ay->cd(3); 
  hdy2->Draw(); 
  g->SetRange(-3.*hdy2->GetRMS(),+3.*hdy2->GetRMS());
  hdy2->Fit("g","Q"); 
  yres[2]=g->GetParameter(2);
  eyres[2]=g->GetParError(2);
  ay->cd(4); 
  hdy3->Draw();
  g->SetRange(-3.*hdy3->GetRMS(),+3.*hdy3->GetRMS());
  hdy3->Fit("g","Q"); 
  yres[3]=g->GetParameter(2);
  eyres[3]=g->GetParError(2);
  ay->cd(5); 
  hdy4->Draw();
  g->SetRange(-3.*hdy4->GetRMS(),+3.*hdy4->GetRMS());
  hdy4->Fit("g","Q");
  yres[4]=g->GetParameter(2);
  eyres[4]=g->GetParError(2);
  ay->cd(6); 
  hdy5->Draw();
  g->SetRange(-3.*hdy5->GetRMS(),+3.*hdy5->GetRMS());
  hdy5->Fit("g","Q");
  yres[5]=g->GetParameter(2);
  eyres[5]=g->GetParError(2);
  // Z
  TCanvas *az = new TCanvas("az","az",0,0,900,700);
  az->Divide(3,2);
  az->cd(1); 
  hdz0->Draw(); 
  g->SetRange(-3.*hdz0->GetRMS(),+3.*hdz0->GetRMS());
  hdz0->Fit("g","Q"); 
  zres[0]=g->GetParameter(2);
  ezres[0]=g->GetParError(2);
  az->cd(2); 
  hdz1->Draw(); 
  g->SetRange(-3.*hdz1->GetRMS(),+3.*hdz1->GetRMS());
  hdz1->Fit("g","Q"); 
  zres[1]=g->GetParameter(2);
  ezres[1]=g->GetParError(2);
  az->cd(3); 
  hdz2->Draw(); 
  g->SetRange(-3.*hdz2->GetRMS(),+3.*hdz2->GetRMS());
  hdz2->Fit("g","Q"); 
  zres[2]=g->GetParameter(2);
  ezres[2]=g->GetParError(2);
  az->cd(4); 
  hdz3->Draw();
  g->SetRange(-3.*hdz3->GetRMS(),+3.*hdz3->GetRMS());
  hdz3->Fit("g","Q"); 
  zres[3]=g->GetParameter(2);
  ezres[3]=g->GetParError(2);
  az->cd(5); 
  hdz4->Draw();
  g->SetRange(-3.*hdz4->GetRMS(),+3.*hdz4->GetRMS());
  hdz4->Fit("g","Q");
  zres[4]=g->GetParameter(2);
  ezres[4]=g->GetParError(2);
  az->cd(6); 
  hdz5->Draw();
  g->SetRange(-3.*hdz5->GetRMS(),+3.*hdz5->GetRMS());
  hdz5->Fit("g","Q");
  zres[5]=g->GetParameter(2);
  ezres[5]=g->GetParError(2);


  //
  // pulls VS multiplicity
  //
  // X
  TCanvas *bx = new TCanvas("bx","bx",0,0,900,700);
  bx->Divide(3,2);
  bx->cd(1); 
  hpx0->Draw(); 
  g->SetRange(-3.*hpx0->GetRMS(),+3.*hpx0->GetRMS());
  hpx0->Fit("g","Q"); 
  xpull[0]=g->GetParameter(2);
  expull[0]=g->GetParError(2);
  bx->cd(2); 
  hpx1->Draw(); 
  g->SetRange(-3.*hpx1->GetRMS(),+3.*hpx1->GetRMS());
  hpx1->Fit("g","Q"); 
  xpull[1]=g->GetParameter(2);
  expull[1]=g->GetParError(2);
  bx->cd(3); 
  hpx2->Draw(); 
  g->SetRange(-3.*hpx2->GetRMS(),+3.*hpx2->GetRMS());
  hpx2->Fit("g","Q"); 
  xpull[2]=g->GetParameter(2);
  expull[2]=g->GetParError(2);
  bx->cd(4); 
  hpx3->Draw();
  g->SetRange(-3.*hpx3->GetRMS(),+3.*hpx3->GetRMS());
  hpx3->Fit("g","Q"); 
  xpull[3]=g->GetParameter(2);
  expull[3]=g->GetParError(2);
  bx->cd(5); 
  hpx4->Draw();
  g->SetRange(-3.*hpx4->GetRMS(),+3.*hpx4->GetRMS());
  hpx4->Fit("g","Q");
  xpull[4]=g->GetParameter(2);
  expull[4]=g->GetParError(2);
  bx->cd(6); 
  hpx5->Draw();
  g->SetRange(-3.*hpx5->GetRMS(),+3.*hpx5->GetRMS());
  hpx5->Fit("g","Q");
  xpull[5]=g->GetParameter(2);
  expull[5]=g->GetParError(2);
  // Y
  TCanvas *by = new TCanvas("by","by",0,0,900,700);
  by->Divide(3,2);
  by->cd(1); 
  hpy0->Draw(); 
  g->SetRange(-3.*hpy0->GetRMS(),+3.*hpy0->GetRMS());
  hpy0->Fit("g","Q"); 
  ypull[0]=g->GetParameter(2);
  eypull[0]=g->GetParError(2);
  by->cd(2); 
  hpy1->Draw(); 
  g->SetRange(-3.*hpy1->GetRMS(),+3.*hpy1->GetRMS());
  hpy1->Fit("g","Q"); 
  ypull[1]=g->GetParameter(2);
  eypull[1]=g->GetParError(2);
  by->cd(3); 
  hpy2->Draw(); 
  g->SetRange(-3.*hpy2->GetRMS(),+3.*hpy2->GetRMS());
  hpy2->Fit("g","Q"); 
  ypull[2]=g->GetParameter(2);
  eypull[2]=g->GetParError(2);
  by->cd(4); 
  hpy3->Draw();
  g->SetRange(-3.*hpy3->GetRMS(),+3.*hpy3->GetRMS());
  hpy3->Fit("g","Q"); 
  ypull[3]=g->GetParameter(2);
  eypull[3]=g->GetParError(2);
  by->cd(5); 
  hpy4->Draw();
  g->SetRange(-3.*hpy4->GetRMS(),+3.*hpy4->GetRMS());
  hpy4->Fit("g","Q");
  ypull[4]=g->GetParameter(2);
  eypull[4]=g->GetParError(2);
  by->cd(6); 
  hpy5->Draw();
  g->SetRange(-3.*hpy5->GetRMS(),+3.*hpy5->GetRMS());
  hpy5->Fit("g","Q");
  ypull[5]=g->GetParameter(2);
  eypull[5]=g->GetParError(2);
  // Z
  TCanvas *bz = new TCanvas("bz","bz",0,0,900,700);
  bz->Divide(3,2);
  bz->cd(1); 
  hpz0->Draw(); 
  g->SetRange(-3.*hpz0->GetRMS(),+3.*hpz0->GetRMS());
  hpz0->Fit("g","Q"); 
  zpull[0]=g->GetParameter(2);
  ezpull[0]=g->GetParError(2);
  bz->cd(2); 
  hpz1->Draw(); 
  g->SetRange(-3.*hpz1->GetRMS(),+3.*hpz1->GetRMS());
  hpz1->Fit("g","Q"); 
  zpull[1]=g->GetParameter(2);
  ezpull[1]=g->GetParError(2);
  bz->cd(3); 
  hpz2->Draw(); 
  g->SetRange(-3.*hpz2->GetRMS(),+3.*hpz2->GetRMS());
  hpz2->Fit("g","Q"); 
  zpull[2]=g->GetParameter(2);
  ezpull[2]=g->GetParError(2);
  bz->cd(4); 
  hpz3->Draw();
  g->SetRange(-3.*hpz3->GetRMS(),+3.*hpz3->GetRMS());
  hpz3->Fit("g","Q"); 
  zpull[3]=g->GetParameter(2);
  ezpull[3]=g->GetParError(2);
  bz->cd(5); 
  hpz4->Draw();
  g->SetRange(-3.*hpz4->GetRMS(),+3.*hpz4->GetRMS());
  hpz4->Fit("g","Q");
  zpull[4]=g->GetParameter(2);
  ezpull[4]=g->GetParError(2);
  bz->cd(6); 
  hpz5->Draw();
  g->SetRange(-3.*hpz5->GetRMS(),+3.*hpz5->GetRMS());
  hpz5->Fit("g","Q");
  zpull[5]=g->GetParameter(2);
  ezpull[5]=g->GetParError(2);

  gStyle->SetOptFit(0);

  // res VS ntracklets
  TCanvas *c6 = new TCanvas("c6","c6",0,0,600,500);
  TH2F *fr6 = new TH2F("fr6","",2,0,80,2,0,240*rangefactor); 
  fr6->SetXTitle("SPD tracklets");
  fr6->SetYTitle("#sigma [#mu m]");
  fr6->Draw();
  sigmamult2[0]=sigmamult2[1]=sigmamult2[2]=sigmamult2[3]=sigmamult2[4]=sigmamult2[5]=0.;
  TGraphErrors *gr6x = new TGraphErrors(6,mult2,xres,sigmamult2,exres);
  gr6x->Draw("p");
  gr6x->SetMarkerStyle(22);
  //gr6x->Fit("func","E,R");
  TGraphErrors *gr6z = new TGraphErrors(6,mult2,zres,sigmamult2,ezres);
  gr6z->Draw("p");
  gr6z->SetMarkerStyle(26);
  //gr6z->Fit("func","E,R");
  TLegend *leg6 = new TLegend(0.6,0.75,0.9,0.9);
  leg6->AddEntry(gr6x,"x/y coordinate","p");
  leg6->AddEntry(gr6z,"z coordinate","p");
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->Draw();

  // pull VS ntracklets
  TCanvas *c8 = new TCanvas("c8","c8",0,0,600,500);
  TH2F *fr8 = new TH2F("fr8","",2,0,80,2,0,2.); 
  fr8->SetXTitle("SPD tracklets");
  fr8->SetYTitle("pull");
  fr8->Draw();
  TGraphErrors *gr8x = new TGraphErrors(6,mult2,xpull,sigmamult2,expull);
  gr8x->Draw("p");
  gr8x->SetMarkerStyle(22);
  TGraphErrors *gr8z = new TGraphErrors(6,mult2,zpull,sigmamult2,ezpull);
  gr8z->Draw("p");
  gr8z->SetMarkerStyle(26);
  TLegend *leg8 = new TLegend(0.6,0.75,0.9,0.9);
  leg8->AddEntry(gr6x,"x/y coordinate","p");
  leg8->AddEntry(gr6z,"z coordinate","p");
  leg8->SetFillStyle(0);
  leg8->SetBorderSize(0);
  leg8->Draw();
  TLine *l8 = new TLine(0,1,80,1);
  l8->SetLineStyle(2);
  l8->Draw();


  // mean z res VS zMC
  Float_t zmc[15]={-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14};
  Float_t ezmc[15]; 
  Float_t dzmean[15],edzmean[15];
  Int_t dummy;
  TCanvas *zz = new TCanvas("zz","zz",0,0,1000,1000);
  zz->Divide(5,3);
  for(l=0;l<15;l++) {
    zz->cd(l+1);
    hdz_z[l].Draw(); 
    g->SetRange(-3.*hdz_z[l].GetRMS(),+3.*hdz_z[l].GetRMS());
    hdz_z[l].Fit("g","Q"); 
    dzmean[l]=g->GetParameter(1);
    edzmean[l]=g->GetParError(1);
    ezmc[l]=1.;
  }
  TCanvas *zzz = new TCanvas("zzz","zzz",0,0,600,500);
  TH2F *fr9 = new TH2F("fr9","",2,-15,15,2,-50*rangefactor,50*rangefactor); 
  fr9->SetXTitle("z_{true} [cm]");
  fr9->SetYTitle("<z_{rec} - z_{true}> [#mu m]");
  fr9->Draw();
  TGraphErrors *grz = new TGraphErrors(15,zmc,dzmean,ezmc,edzmean);
  grz->Draw("pz");
  grz->SetMarkerStyle(20);


  return;
}
//----------------------------------------------------------------------------
void SaveFigures(TString name="dummy") {

  TString namefull;

  namefull = name.Data();
  namefull.Prepend("plots/residualsx_");
  namefull.Append(".gif");
  ax->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsx_");
  namefull.Append(".ps");
  ax->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsx_");
  namefull.Append(".eps");
  ax->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsx_");
  namefull.Append(".root");
  ax->Print(namefull.Data());

  namefull = name.Data();
  namefull.Prepend("plots/pullsx_");
  namefull.Append(".gif");
  bx->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsx_");
  namefull.Append(".ps");
  bx->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsx_");
  namefull.Append(".eps");
  bx->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsx_");
  namefull.Append(".root");
  bx->Print(namefull.Data());

  namefull = name.Data();
  namefull.Prepend("plots/residualsz_");
  namefull.Append(".gif");
  az->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsz_");
  namefull.Append(".ps");
  az->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsz_");
  namefull.Append(".eps");
  az->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/residualsz_");
  namefull.Append(".root");
  az->Print(namefull.Data());

  namefull = name.Data();
  namefull.Prepend("plots/pullsz_");
  namefull.Append(".gif");
  bz->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsz_");
  namefull.Append(".ps");
  bz->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsz_");
  namefull.Append(".eps");
  bz->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsz_");
  namefull.Append(".root");
  bz->Print(namefull.Data());

  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets_");
  namefull.Append(".gif");
  c1a->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets_");
  namefull.Append(".ps");
  c1a->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets_");
  namefull.Append(".eps");
  c1a->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets_");
  namefull.Append(".root");
  c1a->Print(namefull.Data());


  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks_");
  namefull.Append(".gif");
  c1b->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks_");
  namefull.Append(".ps");
  c1b->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks_");
  namefull.Append(".eps");
  c1b->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks_");
  namefull.Append(".root");
  c1b->Print(namefull.Data());


  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks5or6_");
  namefull.Append(".gif");
  c1c->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks5or6_");
  namefull.Append(".ps");
  c1c->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks5or6_");
  namefull.Append(".eps");
  c1c->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSnitstracks5or6_");
  namefull.Append(".root");
  c1c->Print(namefull.Data());


  namefull = name.Data();
  namefull.Prepend("plots/effVSzMC_");
  namefull.Append(".gif");
  c1e->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSzMC_");
  namefull.Append(".ps");
  c1e->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSzMC_");
  namefull.Append(".eps");
  c1e->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSzMC_");
  namefull.Append(".root");
  c1e->Print(namefull.Data());



  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets2_");
  namefull.Append(".gif");
  c5->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets2_");
  namefull.Append(".ps");
  c5->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets2_");
  namefull.Append(".eps");
  c5->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/effVSntracklets2_");
  namefull.Append(".root");
  c5->Print(namefull.Data());


  namefull = name.Data();
  namefull.Prepend("plots/sigmaVSntracklets_");
  namefull.Append(".gif");
  c6->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/sigmaVSntracklets_");
  namefull.Append(".ps");
  c6->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/sigmaVSntracklets_");
  namefull.Append(".eps");
  c6->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/sigmaVSntracklets_");
  namefull.Append(".root");
  c6->Print(namefull.Data());


  namefull = name.Data();
  namefull.Prepend("plots/pullsVSntracklets_");
  namefull.Append(".gif");
  c8->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsVSntracklets_");
  namefull.Append(".ps");
  c8->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsVSntracklets_");
  namefull.Append(".eps");
  c8->Print(namefull.Data());
  namefull = name.Data();
  namefull.Prepend("plots/pullsVSntracklets_");
  namefull.Append(".root");
  c8->Print(namefull.Data());


  return;
}
//----------------------------------------------------------------------------
void TestRmTracks(Int_t iev=0) {

  if (gAlice) {
    delete AliRunLoader::GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval) {
    cerr<<"LoadgAlice returned error"<<endl;
    delete rl;
    return;
  }
  retval = rl->LoadHeader();
  if (retval) {
    cerr<<"LoadHeader returned error"<<endl;
    delete rl;
    return;
  }
  gAlice=rl->GetAliRun();
  // Get field from galice.root
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);
 
  Double_t truevtx[3];
  Double_t vtx[3]; 
  Double_t errvtx[3]; 
  Double_t diffX,diffY,diffZ;
  Double_t diffXerr,diffYerr,diffZerr;
  Double_t multiplicity;
  Int_t nc;

  // -----------   Create vertexer --------------------------
  AliESDVertex *initVertex = 0;
  TFile *invtx = new TFile("AliESDVertexMean.root");
  initVertex = (AliESDVertex*)invtx->Get("vtxmean");
  invtx->Close();
  delete invtx;

  AliVertexerTracks *vertexer = new AliVertexerTracks(AliTracker::GetBz());
  vertexer->SetVtxStart(initVertex);
  vertexer->SetOnlyFitter();
  vertexer->SetDebug(0); // set to 1 to see what it does

  Float_t diamondxy[2];
  diamondxy[0]=initVertex->GetXv();
  diamondxy[1]=initVertex->GetYv();
  //----------------------------------------------------------
 
  if(gSystem->AccessPathName("AliESDs.root",kFileExists)) {
    printf("AliESDs.root not found!\n"); 
    return;
  }
  TFile *fin = TFile::Open("AliESDs.root");
  AliESDEvent *esd = new AliESDEvent();
  TTree *esdTree = (TTree*)fin->Get("esdTree");
  if(!esdTree) return;
  esd->ReadFromTree(esdTree);

  TArrayF o;

  nc=0;


  esdTree->GetEvent(iev);
  Int_t ntracks = esd->GetNumberOfTracks();
  printf("Number of tracks in ESD %d   :   %d\n",iev,ntracks);
  
  // example: do vertex without tracks 1 and 2 (AliESDtrack::GetID())
  //
  // 1: tell vertexer to skip those two tracks
  //
  Int_t nskip=2;
  Int_t skipped[2];
  skipped[0]=1;
  skipped[1]=2;
  vertexer->SetSkipTracks(nskip,skipped);  
  AliESDVertex *vertex1 = (AliESDVertex*)vertexer->FindPrimaryVertex(esd);
  vertex1->PrintStatus();
  vertex1->PrintIndices();
  delete vertex1;
  vertex1 = 0;
  //
  // 2: do vertex with all tracks
  //
  skipped[0]=-99;
  skipped[1]=-99;
  vertexer->SetSkipTracks(nskip,skipped);  
  AliESDVertex *vertex2 = (AliESDVertex*)vertexer->FindPrimaryVertex(esd);
  vertex2->PrintStatus();
  vertex2->PrintIndices();
  //
  // 3: now remove those two tracks from vertex2
  //
  TTree *trksTree = new TTree("trksTree","tree with tracks to be removed");
  AliESDtrack *esdTrack = 0;
  trksTree->Branch("tracks","AliESDtrack",&esdTrack);
  AliESDtrack *et = esd->GetTrack(1);
  esdTrack = new AliESDtrack(*et);
  trksTree->Fill();
  delete esdTrack;
  AliESDtrack *et = esd->GetTrack(2);
  esdTrack = new AliESDtrack(*et);
  trksTree->Fill();
  delete esdTrack;
  AliVertexerTracks vt(AliTracker::GetBz());
  AliESDVertex *vertex3 = vt.RemoveTracksFromVertex(vertex2,trksTree,diamondxy);
  vertex3->PrintStatus();
  vertex3->PrintIndices();
  delete vertex2;
  vertex2 = 0;
  delete vertex3;
  vertex3 = 0;


  delete esdTree;
  delete fin;
  vertexer =0;
  delete vertexer;
  delete initVertex;

  return;
} 
//----------------------------------------------------------------------------
void AliComputeVtxMeanFromESD(TString file="AliESDs.root",  
			      Int_t mincontr=20) {
  //-----------------------------------------------------------------------
  // Compute weighted mean and cov. matrix from a set of AliESD
  //-----------------------------------------------------------------------

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);


  Double_t vtx[3],covvtx[6];
  TTree *esdTree = 0;
  AliESDEvent *esd = new AliESDEvent();
  AliESDVertex *vertex = 0;
  TString vtitle;
  Int_t nc,events,total=0;

  Double_t avx=0.;
  Double_t avy=0.;
  Double_t avz=0.;
  Double_t wgtavx=0.;
  Double_t wgtavy=0.;
  Double_t wgtavz=0.;
  Double_t sum=0.;
  Double_t sumwgtx=0.;
  Double_t sumwgty=0.;
  Double_t sumwgtz=0.;
  Double_t rmsx=0;
  Double_t rmsy=0;
  Double_t rmsz=0;
  Double_t varx=0;
  Double_t vary=0;
  Double_t varz=0;
  Double_t covxy=0;
  Double_t covxz=0;
  Double_t covyz=0;
  Double_t eavx,eavy,eavz,ewgtavx,ewgtavy,ewgtavz;

  TH1F* hx = new TH1F("hx","",200,-1000,1000);
  hx->SetXTitle("vertex x [#mu m]");
  hx->SetYTitle("events");
  TH1F* hy = new TH1F("hy","",200,-1000,1000);
  hy->SetXTitle("vertex y [#mu m]");
  hy->SetYTitle("events");
  TH1F* hz = new TH1F("hz","",200,-10,10);
  hz->SetXTitle("vertex z [cm]");
  hz->SetYTitle("events");


  for(Int_t i=1; i<100; i++) {
    TString inname="event.";
    inname+=i;
    inname.Append("/");
    inname.Append(file.Data());
    if(gSystem->AccessPathName(inname.Data())) continue;
    TFile *fin = TFile::Open(inname.Data());
    esdTree = (TTree*)fin->Get("esdTree");
    if(!esdTree) continue;
    esd->ReadFromTree(esdTree);
    events = esdTree->GetEntries();
    printf("Number of events in ESD tree: %d\n",events);
    total += events;

    for(Int_t iev=0; iev<events; iev++) { //LOOP ON EVENTS
      esdTree->GetEvent(iev);
      vertex = (AliESDVertex*)esd->GetPrimaryVertex();
      if(!vertex) continue;
      nc = vertex->GetNContributors();
      if(nc<mincontr) continue;
      vertex->GetXYZ(vtx);
      vertex->GetSigmaXYZ(errvtx); 
      //printf("%f %f %f %f %f %f\n",vtx[0],vtx[1],vtx[2],errvtx[0],errvtx[1],errvtx[2]);

      avx += vtx[0];
      avy += vtx[1];
      avz += vtx[2];
      sum += 1.;
      wgtavx += vtx[0]/errvtx[0]/errvtx[0];
      wgtavy += vtx[1]/errvtx[1]/errvtx[1];
      wgtavz += vtx[2]/errvtx[2]/errvtx[2];
      sumwgtx += 1./errvtx[0]/errvtx[0];
      sumwgty += 1./errvtx[1]/errvtx[1];
      sumwgtz += 1./errvtx[2]/errvtx[2];
     
      hx->Fill(1e4*vtx[0]);
      hy->Fill(1e4*vtx[1]);
      hz->Fill(vtx[2]);

      vertex = 0;
    }
    esdTree=0;
    fin->Close();
    delete fin;
  }
  
  avx /= sum;
  avy /= sum;
  avz /= sum;
  wgtavx /= sumwgtx;
  wgtavy /= sumwgty;
  wgtavz /= sumwgtz;
  ewgtavx = 1./TMath::Sqrt(sumwgtx);
  ewgtavy = 1./TMath::Sqrt(sumwgty);
  ewgtavz = 1./TMath::Sqrt(sumwgtz);
  
  
  for(Int_t i=1; i<100; i++) {
    TString inname="event.";
    inname+=i;
    inname.Append("/");
    inname.Append(file.Data());
    if(gSystem->AccessPathName(inname.Data())) continue;
    TFile *fin = TFile::Open(inname.Data());
    esdTree = (TTree*)fin->Get("esdTree");
    if(!esdTree) continue;
    esd->ReadFromTree(esdTree);
    events = esdTree->GetEntries();
    for(Int_t iev=0; iev<events; iev++) { //LOOP ON EVENTS
      esdTree->GetEvent(iev);
      vertex = (AliESDVertex*)esd->GetPrimaryVertex();
      if(!vertex) continue;
      nc = vertex->GetNContributors();
      if(nc<mincontr) continue;
      vertex->GetXYZ(vtx);
      vertex->GetSigmaXYZ(errvtx); 

      varx += (vtx[0]-avx)*(vtx[0]-avx);
      vary += (vtx[1]-avy)*(vtx[1]-avy);
      varz += (vtx[2]-avz)*(vtx[2]-avz);
      covxy += (vtx[0]-avx)*(vtx[1]-avy);
      covxz += (vtx[0]-avx)*(vtx[2]-avz);
      covyz += (vtx[1]-avy)*(vtx[2]-avz);

      vertex = 0;
    }
  }
  
  varx /= sum;
  vary /= sum;
  varz /= sum;
  covxy /= sum;
  covxz /= sum;
  covyz /= sum;
  rmsx = TMath::Sqrt(varx);
  rmsy = TMath::Sqrt(vary);
  rmsz = TMath::Sqrt(varz);
  eavx = rmsx/TMath::Sqrt(sum);
  eavy = rmsy/TMath::Sqrt(sum);
  eavz = rmsz/TMath::Sqrt(sum);
  

  printf("\n\nNumber of events: Total %d, Used %d\n",total,sum);
  printf("Average:\n x = (%f +- %f) mum\n y = (%f +- %f) mum\n z = (%f +- %f) mu\n",1e4*avx,1e4*eavx,1e4*avy,1e4*eavy,1e4*avz,1e4*eavz);
  printf("Weighted Average:\n x = (%f +- %f) mum\n y = (%f +- %f) mum\n z = (%f +- %f) mum\n",1e4*wgtavx,1e4*ewgtavx,1e4*wgtavy,1e4*ewgtavy,1e4*wgtavz,1e4*ewgtavz);
  printf("RMS:\n x = %f mum\n y = %f mum\n z = %f mum\n",1e4*rmsx,1e4*rmsy,1e4*rmsz);

  TCanvas *c = new TCanvas("c","c",0,0,1000,500);
  c->Divide(3,1);
  c->cd(1);
  hx->Draw();
  TF1 *gx = new TF1("gx","gaus",-1000,1000);
  gx->SetLineColor(2);
  hx->Fit(gx);
  TF1 *gxx = (TF1*)gx->Clone("gxx");
  gxx->FixParameter(2,50.);
  gxx->SetLineStyle(2);
  gxx->Draw("same");
  c->cd(2);
  hy->Draw();
  TF1 *gy = new TF1("gy","gaus",-1000,1000);
  gy->SetLineColor(2);
  hy->Fit(gy);
  TF1 *gyy = (TF1*)gy->Clone("gyy");
  gyy->FixParameter(2,50.);
  gyy->SetLineStyle(2);
  gyy->Draw("same");
  c->cd(3);
  hz->Draw();
  TF1 *gz = new TF1("gz","gaus",-10,10);
  gz->SetLineColor(2);
  hz->Fit(gz);
  TF1 *gzz = (TF1*)gz->Clone("gzz");
  gzz->FixParameter(2,5.3);
  gzz->SetLineStyle(2);
  gzz->Draw("same");

  vtx[0] = wgtavx;
  vtx[1] = wgtavy;
  vtx[2] = wgtavz;
  covvtx[0] = varx;
  covvtx[1] = covxy;
  covvtx[2] = vary;
  covvtx[3] = covxz;
  covvtx[4] = covyz;
  covvtx[5] = 5.3;
  AliESDVertex *vtxmean = new AliESDVertex(vtx,covvtx,0.,(Int_t)sum);
  vtxmean->SetTitle("Vertex from weighted average");
  TFile *out = new TFile("AliESDVertexMean.root","recreate");
  vtxmean->Write("vtxmean");
  out->Close();

  return;
}
//----------------------------------------------------------------------------
void AliComputeVtxMeanFromTree(TString file="AliESDVertices.root",  
			       Int_t mincontr=20) {
  //-----------------------------------------------------------------------
  // Compute weighted mean and cov. matrix from a set of AliESDVertices
  //-----------------------------------------------------------------------
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);


  Double_t vtx[3],covvtx[6];
  TString vtitle;
  Int_t nc,total=0;

  Double_t avx=0.;
  Double_t avy=0.;
  Double_t avz=0.;
  Double_t wgtavx=0.;
  Double_t wgtavy=0.;
  Double_t wgtavz=0.;
  Double_t sum=0.;
  Double_t sumwgtx=0.;
  Double_t sumwgty=0.;
  Double_t sumwgtz=0.;
  Double_t rmsx=0;
  Double_t rmsy=0;
  Double_t rmsz=0;
  Double_t eavx,eavy,eavz,ewgtavx,ewgtavy,ewgtavz;

  TH1F* hx = new TH1F("hx","",200,-1000,1000);
  hx->SetXTitle("vertex x [#mu m]");
  hx->SetYTitle("events");
  TH1F* hy = new TH1F("hy","",200,-1000,1000);
  hy->SetXTitle("vertex y [#mu m]");
  hy->SetYTitle("events");
  TH1F* hz = new TH1F("hz","",200,-10,10);
  hz->SetXTitle("vertex z [cm]");
  hz->SetYTitle("events");


  if(gSystem->AccessPathName(file.Data())) continue;
  TFile *fin = TFile::Open(file.Data());
  TTree *tree = (TTree*)fin->Get("TreeV");
  AliESDVertex *vertex = 0;
  tree->SetBranchAddress("AliESDVertex",&vertex);
  Int_t entries = (Int_t)tree->GetEntries();

  for(Int_t i=0; i<entries; i++) {
    total += 1;
    tree->GetEvent(i);
    nc = vertex->GetNContributors();
    if(nc<mincontr) continue;
    vertex->GetXYZ(vtx);
    vertex->GetSigmaXYZ(errvtx); 
    //printf("%f %f %f %f %f %f\n",vtx[0],vtx[1],vtx[2],errvtx[0],errvtx[1],errvtx[2]);

    avx += vtx[0];
    avy += vtx[1];
    avz += vtx[2];
    sum += 1.;
    wgtavx += vtx[0]/errvtx[0]/errvtx[0];
    wgtavy += vtx[1]/errvtx[1]/errvtx[1];
    wgtavz += vtx[2]/errvtx[2]/errvtx[2];
    sumwgtx += 1./errvtx[0]/errvtx[0];
    sumwgty += 1./errvtx[1]/errvtx[1];
    sumwgtz += 1./errvtx[2]/errvtx[2];
     
    hx->Fill(1e4*vtx[0]);
    hy->Fill(1e4*vtx[1]);
    hz->Fill(vtx[2]);
    
    vertex = 0;
  }
  
  avx /= sum;
  avy /= sum;
  avz /= sum;
  wgtavx /= sumwgtx;
  wgtavy /= sumwgty;
  wgtavz /= sumwgtz;
  ewgtavx = 1./TMath::Sqrt(sumwgtx);
  ewgtavy = 1./TMath::Sqrt(sumwgty);
  ewgtavz = 1./TMath::Sqrt(sumwgtz);
  
  
  for(Int_t i=0; i<entries; i++) {
    tree->GetEvent(i);
    nc = vertex->GetNContributors();
    if(nc<mincontr) continue;
    vertex->GetXYZ(vtx);
    vertex->GetSigmaXYZ(errvtx); 

    varx += (vtx[0]-avx)*(vtx[0]-avx);
    vary += (vtx[1]-avy)*(vtx[1]-avy);
    varz += (vtx[2]-avz)*(vtx[2]-avz);
    covxy += (vtx[0]-avx)*(vtx[1]-avy);
    covxz += (vtx[0]-avx)*(vtx[2]-avz);
    covyz += (vtx[1]-avy)*(vtx[2]-avz);

    vertex = 0;
  }
  
  varx /= sum;
  vary /= sum;
  varz /= sum;
  covxy /= sum;
  covxz /= sum;
  covyz /= sum;
  rmsx = TMath::Sqrt(varx);
  rmsy = TMath::Sqrt(vary);
  rmsz = TMath::Sqrt(varz);
  eavx = rmsx/TMath::Sqrt(sum);
  eavy = rmsy/TMath::Sqrt(sum);
  eavz = rmsz/TMath::Sqrt(sum);
  

  printf("\n\nNumber of events: Total %d, Used %d\n",total,sum);
  printf("Average:\n x = (%f +- %f) mum\n y = (%f +- %f) mum\n z = (%f +- %f) mu\n",1e4*avx,1e4*eavx,1e4*avy,1e4*eavy,1e4*avz,1e4*eavz);
  printf("Weighted Average:\n x = (%f +- %f) mum\n y = (%f +- %f) mum\n z = (%f +- %f) mum\n",1e4*wgtavx,1e4*ewgtavx,1e4*wgtavy,1e4*ewgtavy,1e4*wgtavz,1e4*ewgtavz);
  printf("RMS:\n x = %f mum\n y = %f mum\n z = %f mum\n",1e4*rmsx,1e4*rmsy,1e4*rmsz);

  TCanvas *c = new TCanvas("c","c",0,0,1000,500);
  c->Divide(3,1);
  c->cd(1);
  hx->Draw();
  TF1 *gx = new TF1("gx","gaus",-1000,1000);
  gx->SetLineColor(2);
  hx->Fit(gx);
  TF1 *gxx = (TF1*)gx->Clone("gxx");
  gxx->FixParameter(2,50.);
  gxx->SetLineStyle(2);
  gxx->Draw("same");
  c->cd(2);
  hy->Draw();
  TF1 *gy = new TF1("gy","gaus",-1000,1000);
  gy->SetLineColor(2);
  hy->Fit(gy);
  TF1 *gyy = (TF1*)gy->Clone("gyy");
  gyy->FixParameter(2,50.);
  gyy->SetLineStyle(2);
  gyy->Draw("same");
  c->cd(3);
  hz->Draw();
  TF1 *gz = new TF1("gz","gaus",-10,10);
  gz->SetLineColor(2);
  hz->Fit(gz);
  TF1 *gzz = (TF1*)gz->Clone("gzz");
  gzz->FixParameter(2,5.3);
  gzz->SetLineStyle(2);
  gzz->Draw("same");

  vtx[0] = wgtavx;
  vtx[1] = wgtavy;
  vtx[2] = wgtavz;
  covvtx[0] = varx;
  covvtx[1] = covxy;
  covvtx[2] = vary;
  covvtx[3] = covxz;
  covvtx[4] = covyz;
  covvtx[5] = 5.3;
  AliESDVertex *vtxmean = new AliESDVertex(vtx,covvtx,0.,(Int_t)sum);
  vtxmean->SetTitle("Vertex from weighted average");
  TFile *out = new TFile("AliESDVertexMean.root","recreate");
  vtxmean->Write("vtxmean");
  out->Close();

  return;
}
//---------------------------------------------------------------------------

