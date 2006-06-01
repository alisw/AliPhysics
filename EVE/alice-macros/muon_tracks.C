void muon_tracks()
{

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();

  rl->LoadTracks("MUON");

  TTree* tt = rl->GetTreeT("MUON", false);

  Alieve::MUONDigitsInfo* di = new Alieve::MUONDigitsInfo();
  di->SetTTree(tt);

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  AliMUON *pMUON = (AliMUON*)gAlice->GetModule("MUON");
  const AliMUONGeometryTransformer* kGeomTransformer = pMUON->GetGeometryTransformer();

  gStyle->SetPalette(1, 0);

  char name[128];
  char title[128];

  /* TRACKS */

  /* Get the z positions of the trigger chambers */

  Float_t xg, yg, zg[4];

  // MT11 = ST6CH1, detElemId = 1100 ... 1117
  // MT12 = ST6CH2, detElemId = 1200 ... 1217
  // MT21 = ST7CH1, detElemId = 1300 ... 1317
  // MT22 = ST7CH2, detElemId = 1400 ... 1417

  kGeomTransformer->Local2Global(1100, 0, 0, 0, xg, yg, zg[0]);
  kGeomTransformer->Local2Global(1200, 0, 0, 0, xg, yg, zg[1]);
  kGeomTransformer->Local2Global(1300, 0, 0, 0, xg, yg, zg[2]);
  kGeomTransformer->Local2Global(1400, 0, 0, 0, xg, yg, zg[3]);

  /* from tracking chambers */

  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);

  TClonesArray *tracks = 0;
  tt->SetBranchAddress("MUONTrack",&tracks);
  tt->GetEntry(0);  // load event 0

  Int_t ntracks = tracks->GetEntriesFast();
  printf("Found %d tracks. \n",ntracks);

  Reve::TrackList* cont = new Reve::TrackList("M-Tracks"); 
  cont->SetMainColor(Color_t(6));
  
  TGListTreeItem *holder = gReve->AddRenderElement(cont);

  Float_t xRec, xRec0;
  Float_t yRec, yRec0;
  Float_t zRec, zRec0;
  
  AliMUONTrack *mt;  
  Reve::RecTrack  rt;
  Int_t count = 0;
  for (Int_t n = 0; n < ntracks; n++) {
    
    mt = (AliMUONTrack*) tracks->At(n);

    printf("Match trigger %d \n",mt->GetMatchTrigger());

    rt.label = n;

    Reve::Track* track = new Reve::Track(&rt, cont->GetRnrStyle());

    if (mt->GetMatchTrigger()) {
      track->SetName(Form("MUONTrack %2d (MT)", rt.label));
      track->SetLineColor(7);
    } else {
      track->SetName(Form("MUONTrack %2d     ", rt.label));
      track->SetLineColor(6);
    }

    AliMUONTrackParam *trackParam = mt->GetTrackParamAtVertex(); 
    xRec0  = trackParam->GetNonBendingCoor();
    yRec0  = trackParam->GetBendingCoor();
    zRec0  = trackParam->GetZ();

    track->SetPoint(count,zRec0,yRec0,xRec0);
    count++;
    
    Float_t xr[20], yr[20], zr[20];
    for (Int_t i = 0; i < 10; i++) xr[i]=yr[i]=zr[i]=0.0;

    Int_t nTrackHits = mt->GetNTrackHits();
    printf("Nhits = %d \n",nTrackHits);
    for (Int_t iHit = 0; iHit < nTrackHits; iHit++){
      trackParamAtHit = mt->GetTrackParamAtHit();
      trackParam = (AliMUONTrackParam*) trackParamAtHit->At(iHit); 
      xRec  = trackParam->GetNonBendingCoor();
      yRec  = trackParam->GetBendingCoor();
      zRec  = trackParam->GetZ();

      //printf("Hit %d x %f y %f z %f \n",iHit,xRec,yRec,zRec);

      xr[iHit] = xRec;
      yr[iHit] = yRec;
      zr[iHit] = zRec;

      track->SetPoint(count,zRec,yRec,xRec);
      count++;
    
    }

    Float_t xrc[20], yrc[20], zrc[20];
    Int_t nrc = 0;
    if (mt->GetMatchTrigger()) {

      for (Int_t i = 0; i < nTrackHits; i++) {
	if (TMath::Abs(zr[i]) > 1000.0) {
	  //printf("Hit %d x %f y %f z %f \n",iHit,xr[i],yr[i],zr[i]);
	  xrc[nrc] = xr[i];
	  yrc[nrc] = yr[i];
	  zrc[nrc] = zr[i];
	  nrc++;
	}
      }

      Double_t xv, yv;
      Float_t ax, bx, ay, by;
      
      // fit x-z
      smatrix.Zero();
      sums.Zero();
      for (Int_t i = 0; i < nrc; i++) {
	xv = (Double_t)zrc[i];
	yv = (Double_t)xrc[i];
	//printf("x-z: xv %f yv %f \n",xv,yv);
	smatrix(0,0) += 1.0;
	smatrix(1,1) += xv*xv;
	smatrix(0,1) += xv;
	smatrix(1,0) += xv;
	sums(0,0)    += yv;
	sums(1,0)    += xv*yv;
      }
      res = smatrix.Invert() * sums;
      ax = res(0,0);
      bx = res(1,0);

      // fit y-z
      smatrix.Zero();
      sums.Zero();
      for (Int_t i = 0; i < nrc; i++) {
	xv = (Double_t)zrc[i];
	yv = (Double_t)yrc[i];
	//printf("y-z: xv %f yv %f \n",xv,yv);
	smatrix(0,0) += 1.0;
	smatrix(1,1) += xv*xv;
	smatrix(0,1) += xv;
	smatrix(1,0) += xv;
	sums(0,0)    += yv;
	sums(1,0)    += xv*yv;
      }
      res = smatrix.Invert() * sums;
      ay = res(0,0);
      by = res(1,0);

      Float_t xtc, ytc, ztc;
      for (Int_t itc = 0; itc < 4; itc++) {

	ztc = zg[itc];
	ytc = ay+by*zg[itc];
	xtc = ax+bx*zg[itc];

	//printf("tc: x %f y %f z %f \n",xtc,ytc,ztc);

	track->SetPoint(count,ztc,ytc,xtc);
	count++;

      }

    }

    cont->AddElement(track);

    gReve->AddRenderElement(holder, track);

  }

  gReve->DrawRenderElement(cont);
  
}
