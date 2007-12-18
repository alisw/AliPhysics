// $Header$

// Functions to read data from HOMER.
//
// Setup: edit location of HLT configuration in first line of
// homer_display(). This is a temporary solution.
//
// Run as: alieve command_queue.C+ hlt_structs.C+ homer_display.C
//
// nextEvent() will get next event from HOMER.


#include "TTimer.h"
#include "TGLViewer.h"
#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"
//#include "AliEVEHOMERManager.h"

class AliRawReaderMemory;

class AliEVEHOMERManager;
class AliHLTHOMERBlockDesc;

namespace Reve {
class PointSet;
class TrackList;
class Track;
}

namespace Alieve {
class TPCLoader;
class TPCData;
class TPCSector2D;
class TPCSector3D;
}

using namespace Reve;
using namespace Alieve;

TPCLoader*  loader  = 0;
TPCData*    tpcdata = 0;
PointSet*   tpc_cls = 0;
TrackList*  tpc_trk = 0;

AliRawReaderMemory* memreader = 0;
AliEVEHOMERManager* homerM = 0;

Int_t    event  = -1;

TTimer   timer;
TTimer   event_timer;

TThread* ldthread = 0;

TRandom  rnd(0);

TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;


//****************************************************************************
void nextEvent();

//****************************************************************************
void process_tpc_clusters(AliHLTHOMERBlockDesc* b);

//****************************************************************************
void homer_display()
{
  homerM = new AliEVEHOMERManager("/local/home/hlt/TPC-SCC1-Generate.xml");
  //  homerM = new AliEVEHOMERManager("/local/home/hlt/sampleConfig2.xml");

  gReve->AddToListTree(homerM, kTRUE);

  homerM->CreateHOMERSourcesList();
  //  homerM->SelectRawTPC();
  homerM->SelectClusterTPC();
  homerM->SelectESDTPC();
  homerM->ConnectHOMER();

  memreader = new AliRawReaderMemory(0, 0); 
  gStyle->SetPalette(1, 0);

  loader = new TPCLoader;
  loader->SetDoubleSR(kTRUE);
  loader->SetInitParams(40, 900, 10, 100); // Sector params (mint, maxt, thr, maxval)
  //loader->SetInitParams(40, 1023, 10); // Sector params (mint, maxt, thr)
  tpcdata = loader->GetData();
  tpcdata->SetLoadPedestal(0);
  tpcdata->SetLoadThreshold(0);
  // tpcdata->SetAutoPedestal(kTRUE); // For non-zero suppressed data.
  tpcdata->SetAutoPedestal(kFALSE);
  gReve->AddRenderElement(loader);

  tpc_cls = new Reve::PointSet("TPC Clusters");
  tpc_cls->SetMainColor((Color_t)kRed);
  gReve->AddRenderElement(tpc_cls);

  tpc_trk = new TrackList("TPC Tracks");
  gReve->AddRenderElement(tpc_trk);
  tpc_trk->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = tpc_trk->GetRnrStyle();
  rnrStyle->SetMagField( 5 );

  nextEvent();
}


//****************************************************************************
void nextEvent()
{
  tpcdata->DropAllSectors();
  tpc_cls->Reset();
  tpc_trk->DestroyElements();

  homerM->NextEvent();
  TIter next(homerM->GetBlockList());
  AliHLTHOMERBlockDesc* b = 0;
  while ((b = (AliHLTHOMERBlockDesc*)next())) {

    //    printf("Q - %s\n", b->GetDataType().Data());

    if (b->GetDataType().CompareTo("CLUSTERS") == 0) {
      process_tpc_clusters(b);
      tpc_cls->ElementChanged();
    }

    else if (b->GetDataType().CompareTo("DDL_RAW") == 0) {
      Int_t slice = b->GetSubDetector().Atoi();
      Int_t patch = b->GetSubSubDetector().Atoi();
      Int_t eqid = 768 + patch;
      if (patch >= 2)
	eqid += 4*slice + 70;
      else
	eqid += 2*slice;

      //printf("%d %d %d -- %p %d\n", slice, patch, eqid, b->GetData(), b->GetSize());

      memreader->SetMemory(b->GetData(), b->GetSize());
      memreader->SetEquipmentID(eqid);
      memreader->Reset();

      AliTPCRawStream input(memreader);
      input.SetOldRCUFormat(kTRUE);
      //input.SetOldRCUFormat(kFALSE);
      memreader->Select("TPC"); // ("TPC", firstRCU, lastRCU);
      tpcdata->LoadRaw(input, kTRUE, kTRUE);
    }
  }

  loader->UpdateSectors(kTRUE); // true -> delete non present
  tpc_cls->ResetBBox();
  tpc_trk->MakeTracks();

  gReve->Redraw3D(1, 1);
}

//****************************************************************************
void process_tpc_clusters(AliHLTHOMERBlockDesc* b)
{
  AliHLTTPCClusterData    *cd = (AliHLTTPCClusterData*) b->GetData();
  UChar_t *data = (UChar_t*) cd->fSpacePoints;

  //  printf("XXX %p %d; sizeof=%d, calcsize=%d\n", b->GetData(), cd->fSpacePointCnt,
  //	 sizeof(AliHLTTPCSpacePointData), (b->GetSize() - 4)/cd->fSpacePointCnt);

  Int_t   slice = b->GetSubDetector().Atoi();
  Float_t phi = (slice+0.5)*TMath::Pi()/9.0;
  Float_t cos = TMath::Cos(phi);
  Float_t sin = TMath::Sin(phi);

  for (Int_t i = 0; i < cd->fSpacePointCnt; ++i, data += sizeof(AliHLTTPCSpacePointData)) {
    AliHLTTPCSpacePointData *sp = (AliHLTTPCSpacePointData *) data;
    //if (i % 100 == 0)
    //printf("  %4d  %6.2f, %6.2f, %6.2f\n", i, sp->fX, sp->fY, sp->fZ);
    tpc_cls->SetNextPoint(cos*sp->fX - sin*sp->fY,
			  sin*sp->fX + cos*sp->fY,
			  sp->fZ);
  }
}

//****************************************************************************
Reve::Track* esd_make_track(Reve::TrackRnrStyle*   rnrStyle,
			    Int_t                  index,
			    AliESDtrack*           at,
			    AliExternalTrackParam* tp=0)
{
  // Helper function
  Double_t        pbuf[3], vbuf[3];
  Reve::RecTrack  rt;

  if(tp == 0) tp = at;

  rt.label  = at->GetLabel();
  rt.index  = index;
  rt.status = (Int_t) at->GetStatus();
  rt.sign   = tp->GetSign();
  tp->GetXYZ(vbuf);
  rt.V.Set(vbuf);
  tp->GetPxPyPz(pbuf);
  rt.P.Set(pbuf);
  Double_t ep = at->GetP(), mc = at->GetMass();
  rt.beta = ep/TMath::Sqrt(ep*ep + mc*mc);
 
  Reve::Track* track = new Reve::Track(&rt, rnrStyle);
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  track->SetName(Form("ESDTrack %d", rt.label));
  //PH  track->SetTitle(Form("pT=%.3f, pZ=%.3f; V=(%.3f, %.3f, %.3f)",
  //PH		       rt.sign*TMath::Hypot(rt.P.x, rt.P.y), rt.P.z,
  //PH		       rt.V.x, rt.V.y, rt.V.z));
  char form[1000];
  sprintf(form,"Track %d", rt.index);
  track->SetName(form);
  track->SetStdTitle();
  return track;
}

//****************************************************************************
void process_tpc_tracks(AliHLTHOMERBlockDesc* b)
{
  AliESDEvent* esd = (AliESDEvent*) b->GetTObject();

  Reve::TrackRnrStyle* rnrStyle = tpc_trk->GetRnrStyle();

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++)
  {
    AliESDtrack           *at = esd->GetTrack(n);
    AliExternalTrackParam *tp = at;

    Reve::Track* track = esd_make_track(rnrStyle, n, at, tp);
    track->SetAttLineAttMarker(tpc_trk);
    gReve->AddRenderElement(track, tpc_trk);
  }

}

//****************************************************************************
void process_tpc_xxxx(AliESDEvent* esd)
{
  Reve::TrackRnrStyle* rnrStyle = tpc_trk->GetRnrStyle();

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++)
  {
    AliESDtrack           *at = esd->GetTrack(n);
    AliExternalTrackParam *tp = at;

    Reve::Track* track = esd_make_track(rnrStyle, n, at, tp);
    track->SetAttLineAttMarker(tpc_trk);
    gReve->AddRenderElement(track, tpc_trk);
  }

}

//****************************************************************************
/*
  // Getting esd

f=TFile::Open("AliESDs.root")
AliESDEvent * esd = new AliESDEvent
esd->ReadFromTree(esdTree)
esdTree->GetEntry(0)
esd->GetNumberOfTracks()
///////////////////

root [3] AliESDEvent * esd = new AliESDEvent
root [4] esd->ReadFromTree(esdTree)
W-TStreamerInfo::BuildOld: Cannot convert AliESDCaloCluster::fTracksMatched from type:TArrayS* to type:TArrayI, skip element
W-TStreamerInfo::BuildOld: Cannot convert AliESDCaloCluster::fLabels from type:TArrayS* to type:TArrayI, skip element
STEER/AliESDEvent.cxx 1040 AliESDEvent::ReadFromTree() TList contains less than the standard contents 21 < 22
root [5] esdTree->GetEntry(0)
(Int_t)(764314)
root [6] esd->GetNumberOfTracks()
(const Int_t)(275)
root [7] process_tpc_xxxx(esd)
root [8] gReve->Redraw3D(
void Redraw3D(Bool_t resetCameras = kFALSE, Bool_t dropLogicals = kFALSE)
root [8] gReve->Redraw3D(1,1)
root [9] trk_cnt->Elem
variable "trk_cnt" not defined.

variable "trk_cnt->Elem" not defined.
root [9] trk_cnt->Eleme
variable "trk_cnt" not defined.

variable "trk_cnt->Eleme" not defined.
root [9] tpc_trk->ElementChanged(
void ElementChanged(Bool_t update_scenes = kTRUE, Bool_t redraw = kFALSE)
root [9] tpc_trk->ElementChanged(1,1)
root [10] tpc_trk->SelectByP
SelectByPt
SelectByPt
SelectByP
SelectByP
root [10] tpc_trk->SelectByPt(
void SelectByPt(Float_t min_pt, Float_t max_pt)
void SelectByPt(Float_t min_pt, Float_t max_pt, Reve::RenderElement* el)
root [10] tpc_trk->SelectByPt(0,1000000000000)
root [11] tpc_trk->MakeTracks(
void MakeTracks(Bool_t recurse = kTRUE)
root [11] tpc_trk->MakeTracks()
root [12]                



*/




//****************************************************************************
void loopEvent()
{
  event_timer.SetCommand("nextEvent()");
  event_timer.Start(60);
}

//****************************************************************************
void stopLoopEvent()
{
  event_timer.Stop();
}
