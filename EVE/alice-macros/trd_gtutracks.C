#ifndef __CINT__
#include "TEveManager.h"
#include "TEveLine.h"
#include "TEveStraightLineSet.h"
#include "TClonesArray.h"
#include "EveBase/AliEveEventManager.h"

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliDataLoader.h"
#include "AliTreeLoader.h"
#include "TRD/AliTRDarrayADC.h"
#include "EveDet/AliEveTRDData.h"
#include "TRD/AliTRDtrackletWord.h"
#include "TRD/AliTRDtrackletMCM.h"
#include "TRD/AliTRDtrackGTU.h"
#include "TRD/AliTRDtrackletGTU.h"
#include "TParticlePDG.h"
#endif

TEveElementList *
trd_gtutracks()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  AliLoader *loader = rl->GetLoader("TRDLoader");
  AliDataLoader *dl = loader->GetDataLoader("gtutracks");
  if (!dl) {
    dl = new AliDataLoader("TRD.GtuTracks.root","gtutracks", "gtutracks");
    rl->GetLoader("TRDLoader")->AddDataLoader(dl);
  }

  dl->Load();
  TTree *trktree = dl->Tree();
  if (!trktree) {
    printf("No GTU track tree");
    return 0x0;
  }
  AliTRDtrackGTU *trk = 0x0;
  trktree->SetBranchAddress("TRDtrackGTU", &trk);
  
  gEve->DisableRedraw();
  TEveElementList* listOfTracks = new TEveElementList("GTU Tracks");
  gEve->AddElement(listOfTracks);

  AliTRDgeometry *geo = new AliTRDgeometry;

  //  printf("found %i tracks\n", trktree->GetEntriesFast());
  for (Int_t i = 0; i < trktree->GetEntriesFast(); i++) {
    trktree->GetEntry(i);
    if (!trk)
      continue;
    Int_t sector = trk->GetSector();
    Double_t *x = new Double_t[3];
    Double_t *p = new Double_t[3];
    Double_t *p2 = new Double_t[3];
    x[0] = x[1] = x[2] = 0;
    x[1] = trk->GetA() / 2 * 160e-4;
    geo->RotateBack(sector*30, x, p);
    x[0] = 400;
    x[1] = x[0] * trk->GetB() + trk->GetA() / 2 * 160e-4;
    x[2] = -1 * x[0] * trk->GetC() / TMath::Tan( -2.0 / 180.0 * TMath::Pi() );
    geo->RotateBack(sector*30, x, p2);
    TEveLine *track = new TEveLine(Form("GTU track pt: %4.2f (%i)", trk->GetPt(), trk->GetLabel()));
    track->SetMainColor((Color_t) 4);
    track->SetNextPoint(p[0], p[1], p[2]);
    track->SetNextPoint(p2[0], p2[1], p2[2]);
    gEve->AddElement(track, listOfTracks);
    delete[] x;
    delete[] p;
    delete[] p2;
//    TEveStraightLineSet* trkl = new TEveStraightLineSet("TRD Tracklets");
//    trkl->SetMainColor((Color_t) 2);
//    gEve->AddElement(trkl, track);
//    Int_t det = trk->GetSector() * 30 + trk->GetStack() * 6;
//    printf("Track in SM: %i, Stack: %i => Det: %i\n", trk->GetSector(), trk->GetStack(), det);
//      for (Int_t layer = 0; layer < 6; layer++) {
//        printf("checking for tracklet in layer %i\n", layer);
//	if (trk->IsTrackletInLayer(layer)) {
//          printf("tracklet in layer %i\n", layer);
//	  AliTRDtrackletGTU *trklet = trk->GetTracklet(layer);
//	  Double_t *x = new Double_t[3];
//	  x[0] = geo->GetTime0(layer);
//	  AliTRDpadPlane *pp = geo->GetPadPlane(layer, trk->GetStack());
//	  x[1] = trklet->GetYbin() * 0.0160; 
//	  x[2] = pp->GetRowPos(trklet->GetZbin()) - pp->GetRowSize(trklet->GetZbin())/2;
//	  Double_t *p = new Double_t[3];
//	  Double_t *p2 = new Double_t[3];
//	  geo->RotateBack(det, x, p);
//	  x[0] -= 10;
//	  x[1] += 10 * trklet->GetdYdX();
//	  x[2] *= x[0]/(x[0]+10);
//	  geo->RotateBack(det, x, p2);
//	  trkl->AddLine(p[0], p[1], p[2], p2[0], p2[1], p2[2]);
//	  delete[] x;
//	  delete[] p;
//	  delete[] p2;
//
//	}
//      }
  }

  gEve->EnableRedraw();
  gEve->Redraw3D();

  return listOfTracks;
}

