// How to fill the list analyser with objects via a macro: Please see example below!
// How to use the list analyser (and especially this macro):

// -- Primary selection (to add e.g. tracklets, tracks, etc.):
// - Load the objects you want to analyse with a sufficient macro (e.g. for tracks and tracklets you can use ana_list_load_tracks.C)
// - Run this macro (ana_list.C)
// - In the tab "eve" in the browser select the list analyser (called "Analysis objects" in the standard case)
// - Select the "list" tab in the editor of this object.
// - Click the button "start"
// - Select the objects you want to add by left-clicking on them in the viewer (or as well in the browser (left panel))
// - When you have finished adding the desired objects, click the button "stop"
// Use the list analyser "as usual" (see class documentation)

// -- Secondary selection (to add e.g. single clusters (of a TEvePointSet) or digits (of a TEveQuadSet)):
// To add e.g. AliTrackPoints or tracks you can do the following:
// - Run e.g. the macro "esd_tracks.C" 
// - Select some track you want to analyse as follows: Hold "shift" und right-click on the track (sometimes you have to hold the right mouse button). The menu pops up
// -> Select "ImportClustersFromIndex"
// -> Do this for all tracks you want to analyse.
// - Run this macro (ana_list.C)
// - In the tab "eve" in the browser select the list analyser (called "Analysis objects" in the standard case)
// - Select the "list" tab in the editor of this object.
// - Click the button "start"
// - Select e.g. clusters by holding "ctrl"+"alt" (depending on your system, holding the "alt" only can also be fine) and left-clicking on the desired cluster
// - When you have finished adding the desired objects, click the button "stop"
// Use the list analyser "as usual" (see class documentation)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>

#include <AliTRDarrayADC.h>
#include <AliEveListAnalyser.h>
#endif

void ana_list(TEveElement *cont = 0)
{
  AliEveListAnalyser * objList = new AliEveListAnalyser("Analysis objects");
  
  objList->SetTitle("Analysis objects (0)");

  gEve->AddElement(objList, cont);

  gEve->Redraw3D();
}


// Example about filling the list analyser with objects via a macro. You can use the following example macro to load the list with trd tracks:
/*
// Launches the list analyser and loads tracks into the list.
// If you already have a list (or rather a list analyser) with objects and you want to add the tracks to this list, do the following:
// Right-click the list analyser in the eve-browser and select "ExportToCint". Choose e.g. "list" for the name in the window that pops up.
// In the console type ".x ana_list_load_tracks.C(list)"
// For more information please see "ana_list.C" or have a look at the class documentation.

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

void ana_list_load_trd_tracks(AliEveListAnalyser* objects = 0, TEveElement *cont = 0)
{
  // Link data containers
  AliESDfriend *eventESDfriend = 0x0;
  if(!(eventESDfriend = AliEveEventManager::AssertESDfriend())){
    Warning("ana_list_load_tracks", "AliESDfriend not found");
    return;
  }

  AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

  AliEveEventManager::Instance()->AssertGeometry();

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


  if (!objects)
  {
    objects = new AliEveListAnalyser("Analysis Objects");
    gEve->AddElement(objects, cont);
  }

  // Number of elements already in the list
  Int_t nOld = objects->NumChildren();

  Int_t count = 0;

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
    AliTRDtrackV1 *trackObj = 0x0;
    AliEveTRDTrack *trackEve = 0x0;
    Int_t ical = 0;

    while((cal = friendTrack->GetCalibObject(ical++))){
      if(strcmp(cal->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      trackObj = new AliTRDtrackV1(*((AliTRDtrackV1*)cal));
      if (!trackObj)
      {
        printf("Cast to AliTRDtrackV1 failed!\n");
        continue;
      }
      if (trackObj->IsOwner()) trackObj->SetOwner();
      trackObj->SetReconstructor(reco);   
      
      trackEve = new AliEveTRDTrack(trackObj);
      if (!trackEve)
      {
        prWintf("Cast to AliEveTRDTrack failed!\n");
        continue;
      }
      objects->AddElement(trackEve);
      trackEve->SetESDstatus(esdTrack->GetStatus());
      trackEve->SetName(Form("[%4d] %s", count + nOld, trackEve->GetName()));
      ++count;
    }
  }

  objects->SetTitle(Form("Objects %d", objects->NumChildren()));
  objects->StampObjProps();

  gEve->Redraw3D();

  return;
}
*/
