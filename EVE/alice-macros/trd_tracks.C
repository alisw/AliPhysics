#ifndef __CINT__
#include <TGLViewer.h>
#include <TEveManager.h>
#include <EveBase/AliEveEventManager.h>
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
  AliESDEvent* esd(AliEveEventManager::AssertESD());
  AliESDfriend *esdFriend(AliEveEventManager::AssertESDfriend());
  if(!esd || !esdFriend){
    Warning("trd_tracks", "Full ESD data missing.");
    return;
  }

/*  AliEveEventManager::AssertGeometry();
  AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
*/
  AliTRDrecoParam *trdRecoParam = AliTRDrecoParam::GetLowFluxParam();
  trdRecoParam->SetPIDNeuralNetwork();
  AliTRDReconstructor *reco = new AliTRDReconstructor();
  reco->SetRecoParam(trdRecoParam);

  AliEveTRDTrackList *tracks = new AliEveTRDTrackList("TRD Tracks");
  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++){
    AliESDtrack* esdTrack(esd->GetTrack(n));
    if(!esdTrack) continue;
    AliESDfriendTrack *friendTrack(esdTrack->GetFriendTrack());
    if(!friendTrack) continue;
    //Info("trd_tracks", Form("Track[%3d] esd[%p] friend[%p]", n, (void*)esdTrack, (void*)friendTrack));

    TObject *cal(NULL);
    Int_t ical(0);
    while((cal = friendTrack->GetCalibObject(ical++))){
      //Info("trd_tracks", Form("  Obj[%d] %s", ical-1, cal->IsA()->GetName()));
      if(strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      AliTRDtrackV1 *trackObj = dynamic_cast<AliTRDtrackV1 *>(cal);
      trackObj->SetReconstructor(reco);
      AliEveTRDTrack *trackEve = new AliEveTRDTrack(new AliTRDtrackV1(*trackObj));
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

  //  TGLViewer *v = gEve->GetDefaultGLViewer();
  //  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  //  ((TGLOrthoCamera&)v->CurrentCamera()).SetEnableRotate(kTRUE);
  //  v->UpdateScene();

  return;
}
