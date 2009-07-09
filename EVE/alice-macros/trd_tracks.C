#ifndef __CINT__
#include <TGLViewer.h>
#include <TEveManager.h>
#include <EveBase/AliEveEventManager.h>
#include "TRD/AliTRDdataArrayI.h"
#include "TRD/AliTRDarrayADC.h"
#include <EveDet/AliEveTRDTrackList.h>

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "TRD/AliTRDReconstructor.h"
#include "TRD/AliTRDtrackV1.h"
#endif

void trd_tracks(TEveElement *cont = 0)
{
  // Link data containers
  AliESDfriend *eventESDfriend = 0x0;
  if(!(eventESDfriend = AliEveEventManager::AssertESDfriend())){
    Warning("trd_tracks", "AliESDfriend not found");
    return;
  }

  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliEveEventManager::AssertGeometry();

  AliTRDReconstructor *reco = new AliTRDReconstructor();
  reco->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  reco->SetOption("!nn");

  AliEveTRDTrackList *tracks = new AliEveTRDTrackList("TRD Tracks");
  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++){
    AliESDtrack* esdTrack = esd->GetTrack(n);
    AliESDfriendTrack *friendTrack = eventESDfriend->GetTrack(n);

    TObject *cal = 0x0;
    Int_t ical = 0;
    while((cal = friendTrack->GetCalibObject(ical++))){
      if(strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      AliTRDtrackV1 *trackObj = dynamic_cast<AliTRDtrackV1 *>(cal);
      trackObj->SetReconstructor(reco);
      AliEveTRDTrack *trackEve = new AliEveTRDTrack(trackObj);
      tracks->AddElement(trackEve);
      trackEve->SetESDstatus(esdTrack->GetStatus());
      trackEve->SetName(Form("[%4d] %s", n, trackEve->GetName()));
    }
  }

  delete reco;

  tracks->SetTitle(Form("Tracks %d", tracks->NumChildren()));
  tracks->StampObjProps();
  gEve->AddElement(tracks, cont);

  gEve->Redraw3D();

  TGLViewer *v = gEve->GetDefaultGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  ((TGLOrthoCamera&)v->CurrentCamera()).SetEnableRotate(kTRUE);
  v->UpdateScene();

  return;
}
