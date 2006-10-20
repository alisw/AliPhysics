//------------------------------------------------------------------------
// Test macro for AliVertexerTracks.
// Contains several functions. 
//
// A. Dainese - INFN Legnaro
//------------------------------------------------------------------------
void AliVertexerTracksTest(Double_t nSigma=3.,
			   Bool_t useMeanVtx=kTRUE,
			   TString outname="AliVertexerTracksTest.root") {
//------------------------------------------------------------------------
// Run vertex reconstruction and store results in histos and ntuple
//------------------------------------------------------------------------

  if (gAlice) {
    delete gAlice->GetRunLoader();
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
  rl->LoadKinematics();
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
  Double_t multiplicity;
  Int_t nc;

  // Histograms
  TH1F *hdiffX = new TH1F("hdiffX","x_{found} - x_{true} [#mum]",1000,-5000,5000);
  TH1F *hdiffY = new TH1F("hdiffY","y_{found} - y_{true} [#mum]",1000,-5000,5000);
  TH1F *hdiffZ = new TH1F("hdiffZ","z_{found} - z_{true} [#mum]",1000,-5000,5000);

  TH1F *hdiffXerr = new TH1F("hdiffXerr","#Delta x/#sigma_{x}",100,-5,5);
  TH1F *hdiffYerr = new TH1F("hdiffYerr","#Delta y/#sigma_{y}",100,-5,5);
  TH1F *hdiffZerr = new TH1F("hdiffZerr","#Delta z/#sigma_{z}",100,-5,5);

  //TNtple
  TNtuple *ntVtxRes = new TNtuple("ntVtxRes","Vertex Resolution's Ntupla with multiplicity","ntracks:nitstracks5or6:nitstracksFromStrange:nitstracksFromStrange5or6:diffX:diffY:diffZ:diffXerr:diffYerr:diffZerr:multiplicity:nc");



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
  AliVertexerTracks *vertexer = new AliVertexerTracks();
  vertexer->SetITSrefitNotRequired();
  vertexer->SetVtxStart(initVertex);
  vertexer->SetMinTracks(2);
  vertexer->SetNSigmad0(nSigma);
  //cout<<" Nsigma:  "<<vertexer->GetNSigmad0()<<endl;
  vertexer->SetDebug(1); // set to 1 to see what it does
  //----------------------------------------------------------
 
  if(gSystem->AccessPathName("AliESDs.root",kFileExists)) {
    printf("AliESDs.root not found!\n"); 
    return;
  }
  TFile *fin = TFile::Open("AliESDs.root");
  TTree *esdTree = (TTree*)fin->Get("esdTree");
  if(!esdTree) return;
  AliESD *esd = 0;
  esdTree->SetBranchAddress("ESD",&esd);
  Int_t events = esdTree->GetEntries();
  printf("Number of events in ESD tree: %d\n",events);

  TArrayF o;

  for(Int_t iev=0; iev<events; iev++) { //LOOP ON EVENTS
    printf("-------------- EVENT   %d --------------------\n",iev);

    diffX=99999;
    diffY=99999;
    diffZ=99999;
    diffXerr=99999;
    diffYerr=99999;
    diffZerr=99999;
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
    Int_t nitstracksFromStrange = 0;
    Int_t nitstracksFromStrange5or6 = 0;
    Int_t nitstracks1 = 0;
    Int_t nitstracks2 = 0;
    Int_t nitstracks3 = 0;
    Int_t nitstracks4 = 0;
    Int_t nitstracks5 = 0;
    Int_t nitstracks6 = 0;
    Int_t nitstracks5or6 = 0;
  
    // loop on particles
    for(Int_t pa=0; pa<npart; pa++) {   
      TParticle *part = (TParticle*)gAlice->GetMCApp()->Particle(pa); 
      Int_t pdg = part->GetPdgCode();      
      Int_t apdg = TMath::Abs(pdg);
      Double_t energy  = part->Energy();
      if(energy>6900.) continue; // reject incoming protons
      Double_t pz = part->Pz();
      Double_t y = 0.5*TMath::Log((energy+pz+1.e-13)/(energy-pz+1.e-13)); 
                                                 
      // consider charged pi,K,p,el,mu       
      if(apdg!=11 && apdg!=13 && apdg!=211 && apdg!=321 && apdg!=2212) continue;
      // reject secondaries
      if(TMath::Sqrt((part->Vx()-o[0])*(part->Vx()-o[0])+(part->Vy()-o[1])*(part->Vy()-o[1]))>0.0010) continue;
      if(TMath::Abs(y)<1.0) dNchdy += 0.5; // count 1/2 of particles in |y|<1
    } // end loop on particles
    multiplicity=(Double_t)dNchdy;

    printf(" dNch/dy = %f\n",dNchdy);
    //===============================================================


    esdTree->GetEvent(iev);
    Int_t ntracks = esd->GetNumberOfTracks();
    printf("Number of tracks in ESD %d   :   %d\n",iev,ntracks);
    if(ntracks==0) { 
      ntVtxRes->Fill(ntracks,nitstracks5or6,nitstracksFromStrange,nitstracksFromStrange5or6,diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr,multiplicity,nc); 
      continue; 
    }
    
    // VERTEX    
    AliESDVertex *vertex = (AliESDVertex*)vertexer->FindPrimaryVertex(esd);
    //AliESDVertex *vertex = (AliESDVertex*)vertexer->FindPrimaryVertexOld(esd);

    nc = (Int_t)vertex->GetNContributors();
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
      TParticle *part = (TParticle*)gAlice->GetMCApp()->Particle(TMath::Abs(esdtrack->GetLabel()));
      if(part->GetFirstMother()>=0) {
	Int_t mpdg = TMath::Abs(gAlice->GetMCApp()->Particle(part->GetFirstMother())->GetPdgCode());
	//cout<<TMath::Abs(esdtrack->GetLabel())<<"   "<<mpdg<<endl;
	if(mpdg==321||mpdg==311||mpdg==310||mpdg==3122||mpdg==3312||mpdg==3334) {
	  nitstracksFromStrange++;
	  if(npts>=5) nitstracksFromStrange5or6++;
	}
      }
    }
    printf("Number of kITSin tracks in ESD %d   :   %d\n",iev,nitstracks);
    printf("           6 pts: %d\n",nitstracks6);
    printf("           5 pts: %d\n",nitstracks5);
    printf("           4 pts: %d\n",nitstracks4);
    printf("           3 pts: %d\n",nitstracks3);
    printf("           2 pts: %d\n",nitstracks2);
    printf("           1 pts: %d\n",nitstracks1);
    printf("Number of kITSin tracks from s:   %d\n",nitstracksFromStrange);


    if(nc>=2) {
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

    ntVtxRes->Fill(nitstracks,nitstracks5or6,nitstracksFromStrange,nitstracksFromStrange5or6,diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr,multiplicity,(Float_t)nc);   

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
void PlotVtxRes(const Char_t *inname="AliVertexerTracksTest.root") {
  //---------------------------------------------------------------------
  // Plot efficiency, resolutions, pulls 
  // [at least 10k pp events are needed]
  //---------------------------------------------------------------------

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(10001);

  Int_t   nEvVtx=0;
  Int_t   nEvNoVtx=0;
  Int_t   ev,nUsedTrks;
  Float_t nTrks,nTrksFromStrange,nTrks5or6,nTrksFromStrange5or6,dNchdy;
  Float_t diffX,diffY,diffZ,diffXerr,diffYerr,diffZerr;

  //
  // HISTOGRAMS
  //
  TH1F *hdx = new TH1F("hdx","",50,-600,600);
  hdx->SetXTitle("#Delta x [#mu m]");
  hdx->SetYTitle("events");
  //
  TH1F *hdy = new TH1F("hdy","",50,-600,600);
  hdy->SetXTitle("#Delta y [#mu m]");
  hdy->SetYTitle("events");
  //
  TH1F *hdz = new TH1F("hdz","",50,-600,600);
  hdz->SetXTitle("#Delta z [#mu m]");
  hdz->SetYTitle("events");
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
  hTracksVtx5or6->SetXTitle("Number of reco tracks in TPC&ITS(pts>=5)");
  hTracksVtx5or6->SetYTitle("events");
  //
  TH1F *hTracksNoVtx5or6 = new TH1F("hTracksNoVtx5or6","events w/o vertex",51,-0.5,50.5);
  hTracksNoVtx5or6->SetXTitle("Number of reco tracks in TPC&ITS(pts>=5)");
  hTracksNoVtx5or6->SetYTitle("events");
  //
  TH1F *hTracksVtx5or6nonS = new TH1F("hTracksVtx5or6nonS","events w/ vertex",51,-0.5,50.5);
  hTracksVtx5or6nonS->SetXTitle("Number of reco tracks in TPC&ITS(pts>=5) NOT from s");
  hTracksVtx5or6nonS->SetYTitle("events");
  //
  TH1F *hTracksNoVtx5or6nonS = new TH1F("hTracksNoVtx5or6nonS","events w/o vertex",51,-0.5,50.5);
  hTracksNoVtx5or6nonS->SetXTitle("Number of reco tracks in TPC&ITS(pts>=5) NOT from s");
  hTracksNoVtx5or6nonS->SetYTitle("events");
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
  Int_t    bin;
  //
  // OTHER HISTOGRAMS
  //
  TH1F *hdx0 = new TH1F("hdx0","x resolution - bin 0",50,-500,500);
  TH1F *hdx1 = new TH1F("hdx1","x resolution - bin 1",50,-500,500);
  TH1F *hdx2 = new TH1F("hdx2","x resolution - bin 2",50,-500,500);
  TH1F *hdx3 = new TH1F("hdx3","x resolution - bin 3",50,-400,400);
  TH1F *hdx4 = new TH1F("hdx4","x resolution - bin 4",50,-300,300);
  TH1F *hdx5 = new TH1F("hdx5","x resolution - bin 5",50,-300,300);
  TH1F *hdy0 = new TH1F("hdy0","y resolution - bin 0",50,-500,500);
  TH1F *hdy1 = new TH1F("hdy1","y resolution - bin 1",50,-500,500);
  TH1F *hdy2 = new TH1F("hdy2","y resolution - bin 2",50,-500,500);
  TH1F *hdy3 = new TH1F("hdy3","y resolution - bin 3",50,-400,400);
  TH1F *hdy4 = new TH1F("hdy4","y resolution - bin 4",50,-300,300);
  TH1F *hdy5 = new TH1F("hdy5","y resolution - bin 5",50,-300,300);
  TH1F *hdz0 = new TH1F("hdz0","z resolution - bin 0",50,-1000,1000);
  TH1F *hdz1 = new TH1F("hdz1","z resolution - bin 1",50,-800,800);
  TH1F *hdz2 = new TH1F("hdz2","z resolution - bin 2",50,-800,800);
  TH1F *hdz3 = new TH1F("hdz3","z resolution - bin 3",50,-600,600);
  TH1F *hdz4 = new TH1F("hdz4","z resolution - bin 4",50,-500,500);
  TH1F *hdz5 = new TH1F("hdz5","z resolution - bin 5",50,-500,500);

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

  TFile *in = new TFile(inname);
  TNtuple *nt = (TNtuple*)in->Get("ntVtxRes");
  nt->SetBranchAddress("ntracks",&nTrks);
  nt->SetBranchAddress("nitstracks5or6",&nTrks5or6);
  nt->SetBranchAddress("nitstracksFromStrange",&nTrksFromStrange);
  nt->SetBranchAddress("nitstracksFromStrange5or6",&nTrksFromStrange5or6);
  nt->SetBranchAddress("diffX",&diffX);
  nt->SetBranchAddress("diffY",&diffY);
  nt->SetBranchAddress("diffZ",&diffZ);
  nt->SetBranchAddress("diffXerr",&diffXerr);
  nt->SetBranchAddress("diffYerr",&diffYerr);
  nt->SetBranchAddress("diffZerr",&diffZerr);
  nt->SetBranchAddress("multiplicity",&dNchdy);
  Int_t entries = (Int_t)nt->GetEntries();
  Int_t nbytes=0;

  for(Int_t i=0; i<entries; i++) {
    nbytes += nt->GetEvent(i);
    bin = GetBin(dNchdy);
    mult[bin] += dNchdy;
    hmult->Fill(dNchdy);
    totevts[bin]++;
    if(diffX < 90000.) { // vtx found
      nEvVtx++;
      mult2[bin] += dNchdy;
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
      hTracksVtx->Fill(nTrks);
      hTracksVtx5or6->Fill(nTrks5or6);
      hTracksVtx5or6nonS->Fill(nTrks5or6-nTrksFromStrange5or6);
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
      hTracksNoVtx->Fill(nTrks);
      hTracksNoVtx5or6->Fill(nTrks5or6);
      hTracksNoVtx5or6nonS->Fill(nTrks5or6-nTrksFromStrange5or6);
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
    bin = GetBin(dNchdy);
    sigmamult[bin] += (dNchdy-mult[bin])*(dNchdy-mult[bin])/totevts[bin];
    if(diffX < 90000.) sigmamult2[bin] += (dNchdy-mult2[bin])*(dNchdy-mult2[bin])/evts[bin];
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
  TF1 *func = new TF1("func",FitFunc,2,35,2);
  func->SetLineWidth(1);

  Double_t x1,y1,sigmafit;
  Char_t message[100];

  //
  //            DRAW PLOTS
  //
  gStyle->SetOptFit(1111);
  
  // Tracks in ITS for events w/ and w/o vertex
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,500);
  c1->Divide(2,1);
  c1->cd(1); 
  hTracksVtx->Draw();
  c1->cd(2); 
  hTracksNoVtx->Draw();

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

  TCanvas *c1d = new TCanvas("c1d","c1d",0,0,500,500);
  TH1F *hTotTracks5or6nonS = (TH1F*)hTracksVtx5or6nonS->Clone("hTotTracks5or6nonS");
  hTotTracks5or6nonS->Add(hTracksNoVtx5or6nonS);
  TH1F *hProbTracks5or6nonS = (TH1F*)hTracksVtx5or6nonS->Clone("hProbTracks5or6nonS");
  hProbTracks5or6nonS->Divide(hTotTracks5or6nonS);
  hProbTracks5or6nonS->SetYTitle("Probability");
  hProbTracks5or6nonS->SetTitle("Probability to find the vertex");
  hProbTracks5or6nonS->Draw();

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
  TH2F *fr5 = new TH2F("fr5","",2,0,35,2,0,1.1); 
  fr5->SetXTitle("dN_{ch}/dy");
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

  // res VS dNchdy
  TCanvas *c6 = new TCanvas("c6","c6",0,0,600,500);
  TH2F *fr6 = new TH2F("fr6","",2,0,35,2,0,200); 
  fr6->SetXTitle("dN_{ch}/dy");
  fr6->SetYTitle("#sigma [#mu m]");
  fr6->Draw();
  sigmamult2[0]=sigmamult2[1]=sigmamult2[2]=sigmamult2[3]=sigmamult2[4]=sigmamult2[5]=0.;
  TGraphErrors *gr6x = new TGraphErrors(6,mult2,xres,sigmamult2,exres);
  gr6x->Draw("p");
  gr6x->SetMarkerStyle(22);
  gr6x->Fit("func","E,R");
  TGraphErrors *gr6z = new TGraphErrors(6,mult2,zres,sigmamult2,ezres);
  gr6z->Draw("p");
  gr6z->SetMarkerStyle(26);
  gr6z->Fit("func","E,R");
  TLegend *leg6 = new TLegend(0.6,0.75,0.9,0.9);
  leg6->AddEntry(gr6x,"x/y coordinate","p");
  leg6->AddEntry(gr6z,"z coordinate","p");
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->Draw();

  // pull VS dNchdy
  TCanvas *c8 = new TCanvas("c8","c8",0,0,600,500);
  TH2F *fr8 = new TH2F("fr8","",2,0,35,2,0,2.); 
  fr8->SetXTitle("dN_{ch}/dy");
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
  TLine *l8 = new TLine(0,1,35,1);
  l8->SetLineStyle(2);
  l8->Draw();

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
  AliESD *esd = 0;
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
    esdTree->SetBranchAddress("ESD",&esd);
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
    esdTree->SetBranchAddress("ESD",&esd);
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
