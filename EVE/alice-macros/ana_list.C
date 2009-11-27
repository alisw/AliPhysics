#ifndef __CINT__
#include <TGLViewer.h>
#include <TEveManager.h>
#include <EveBase/AliEveEventManager.h>
#include "TRD/AliTRDarrayADC.h"
#include <EveDet/AliEveListAnalyser.h>

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "TRD/AliTRDReconstructor.h"
#include "TRD/AliTRDtrackV1.h"
#endif

void ana_list(TEveElement *cont = 0)
{
  // Link data containers
  AliESDfriend *eventESDfriend = 0x0;
  if(!(eventESDfriend = AliEveEventManager::AssertESDfriend())){
    Warning("trd_loadObjectList", "AliESDfriend not found");
    return;
  }

  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliEveEventManager::AssertGeometry();

  AliTRDrecoParam *trdRecoParam = AliTRDrecoParam::GetLowFluxParam();
  if (!trdRecoParam)
  {
    printf("Could not load AliTRDrecoParam\n");
    return;
  }
  trdRecoParam->SetPIDNeuralNetwork();
  AliTRDReconstructor *reco = new AliTRDReconstructor();
  if (!reco)
  {
    printf("Could not load AliTRDReconstructor\n");
    return;
  }
  reco->SetRecoParam(trdRecoParam);


  AliEveListAnalyser *objects = new AliEveListAnalyser("TRD Analysis Object");

  for (Int_t n = 0; n < esd->GetNumberOfTracks(); n++)
  {
    AliESDtrack* esdTrack = esd->GetTrack(n);
    AliESDfriendTrack *friendTrack = eventESDfriend->GetTrack(n);

    if (!esdTrack || !friendTrack)
    {
      printf("Problem with track %d\n", n);
      continue;
    }

    TObject *cal = 0x0;
    Int_t ical = 0;

    while((cal = friendTrack->GetCalibObject(ical++))){
      if(strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      AliTRDtrackV1 *trackObj = dynamic_cast<AliTRDtrackV1 *>(cal);
      if (!trackObj)
      {
        printf("Cast to AliTRDtrackV1 failed!\n");
        continue;
      }
      trackObj->SetReconstructor(reco);
      AliEveTRDTrack *trackEve = new AliEveTRDTrack(trackObj);
      if (!trackEve)
      {
        printf("Cast to AliEveTRDTrack failed!\n");
        continue;
      }
      objects->AddElement(trackEve);
      trackEve->SetESDstatus(esdTrack->GetStatus());
      trackEve->SetName(Form("[%4d] %s", n, trackEve->GetName()));
    }
  }

  delete reco;

  objects->SetTitle(Form("Objects %d", objects->NumChildren()));
  objects->StampObjProps();
 
  gEve->AddElement(objects, cont);


  gEve->Redraw3D();

  //  TGLViewer *v = gEve->GetDefaultGLViewer();
  //  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  //  ((TGLOrthoCamera&)v->CurrentCamera()).SetEnableRotate(kTRUE);
  //  v->UpdateScene();

  return;
}
