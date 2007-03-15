#include "TGLViewer.h"

namespace Alieve {
class MUONData;
class Event;
}

Alieve::MUONData* g_muon_data = 0;
Alieve::Event*    g_muon_last_event = 0;

Int_t g_currentEvent = -1;
Bool_t g_fromRaw = kFALSE;

void MUON_display(Bool_t fromRaw = kFALSE, Bool_t showTracks = kTRUE)
{
 
  TTree* dt = 0;
  TTree* ct = 0;

  if (Alieve::gEvent == 0) {
    printf("No alieve event: use alieve_init(...) \n");
    return;
  }

  if (g_currentEvent == Alieve::gEvent->GetEventId()) {
    if (g_fromRaw == fromRaw) {
      printf("Same event... \n");
      return;
    } else {
      if (g_fromRaw) {
	printf("Same event with digits.\n");
	Alieve::gEvent->GotoEvent(g_currentEvent);
      } else {
	printf("Same event with raw.\n");
	Alieve::gEvent->GotoEvent(g_currentEvent);
      }
    }
  }

  g_fromRaw = fromRaw;

  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/raw.root");

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  g_muon_data = new Alieve::MUONData;

  if (!fromRaw) {
    rl->LoadDigits("MUON");
    dt = rl->GetTreeD("MUON", false);
    if (dt == 0) {
      cout << "No digits produced!" << endl;
    } else {
      cout << "Display aliroot digits!" << endl;
      g_muon_data->LoadDigits(dt);
    }
  } else {
    if (gSystem->AccessPathName(dataPath.Data(),kFileExists)) {
      cout << "No raw data produced!" << endl;
    } else {
      cout << "Display raw digits!" << endl;
      g_muon_data->LoadRaw(dataPath.Data());
    }
  }

  rl->LoadRecPoints("MUON");
  ct = rl->GetTreeR("MUON", false);
  g_muon_data->LoadRecPoints(ct);

  g_muon_last_event = Alieve::gEvent;
  
  g_currentEvent = g_muon_last_event->GetEventId();
  
  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();

  Reve::RenderElementList* l = new Reve::RenderElementList("MUONChambers");
  l->SetTitle("MUON chambers");
  l->SetMainColor(Color_t(2));
  gReve->AddRenderElement(l);

  for (Int_t ic = 0; ic < 14; ic++) {

    Alieve::MUONChamber* mucha = new Alieve::MUONChamber();
    
    mucha->SetFrameColor(2);
    mucha->SetChamberID(ic);
    mucha->SetDataSource(g_muon_data);

    gReve->AddRenderElement(l,mucha);

  }

  if (showTracks) MUON_tracks();

  gReve->EnableRedraw();
  gReve->Redraw3D(kTRUE);
  /*
  TGLViewer* view = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
  view->ResetCamerasAfterNextUpdate();
  gReve->GetGLCanvas()->Modified();
  gReve->GetGLCanvas()->Update();
  */
}

//_____________________________________________________________________________
void MUON_tracks() {

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadTracks("MUON");
  TTree* tt = rl->GetTreeT("MUON", false);

  AliMagF* mag_field = Alieve::Event::AssertMagField();

  TClonesArray *tracks = 0;
  tt->SetBranchAddress("MUONTrack",&tracks);
  tt->GetEntry(0);

  Int_t ntracks = tracks->GetEntriesFast();
  //printf("Found %d tracks. \n",ntracks);

  Reve::TrackList* lt = new Reve::TrackList("M-Tracks"); 
  lt->SetMainColor(Color_t(6));
  
  gReve->AddRenderElement(lt);

  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);

  Float_t xRec, xRec0;
  Float_t yRec, yRec0;
  Float_t zRec, zRec0;
  
  Float_t zg[4] = { -1603.5, -1620.5, -1703.5, -1720.5 };

  AliMUONTrack *mt;  
  Reve::RecTrack  rt;
  Int_t count;
  for (Int_t n = 0; n < ntracks; n++) {
    
    count = 0;

    mt = (AliMUONTrack*) tracks->At(n);

    //printf("Match trigger %d \n",mt->GetMatchTrigger());

    rt.label = n;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeTrack(mt, mag_field);

    gReve->AddRenderElement(lt, track);

    continue;

    /* old tracks display */

    Reve::Track* rtrack = new Reve::Track(&rt, lt->GetRnrStyle());

    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    char form[1000];
    if (mt->GetMatchTrigger()) {
      //PH      track->SetName(Form("MUONTrack %2d (MT)", rt.label));
      sprintf(form,"MUONTrack %2d (MT)", rt.label);
      rtrack->SetName(form);
      rtrack->SetLineStyle(1);
      rtrack->SetLineColor(4);
    } else {
      //PH      track->SetName(Form("MUONTrack %2d     ", rt.label));
      sprintf(form,"MUONTrack %2d     ", rt.label);
      rtrack->SetName(form);
      rtrack->SetLineStyle(2);
      rtrack->SetLineColor(4);
    }

    AliMUONTrackParam *trackParam = mt->GetTrackParamAtVertex(); 
    xRec0  = trackParam->GetNonBendingCoor();
    yRec0  = trackParam->GetBendingCoor();
    zRec0  = trackParam->GetZ();

    rtrack->SetPoint(count,xRec0,yRec0,zRec0);
    count++;
    
    Float_t xr[20], yr[20], zr[20];
    for (Int_t i = 0; i < 10; i++) xr[i]=yr[i]=zr[i]=0.0;

    Int_t nTrackHits = mt->GetNTrackHits();
    //printf("Nhits = %d \n",nTrackHits);
    TClonesArray* trackParamAtHit;
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

      rtrack->SetPoint(count,xRec,yRec,zRec);
      count++;
    
    }

    Float_t xrc[20], yrc[20], zrc[20];
    Int_t nrc = 0;
    if (mt->GetMatchTrigger() && 1) {

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
      
      if (nrc < 2) continue;

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
      for (Int_t ii = 0; ii < 4; ii++) {

	ztc = zg[ii];
	ytc = ay+by*zg[ii];
	xtc = ax+bx*zg[ii];

	//printf("tc: x %f y %f z %f \n",xtc,ytc,ztc);

	rtrack->SetPoint(count,xtc,ytc,ztc);
	count++;

      }

    }  // end match trigger

    gReve->AddRenderElement(lt, rtrack);

  }  // end track loop

}

