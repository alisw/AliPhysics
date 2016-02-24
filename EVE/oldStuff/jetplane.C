// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TGLViewer.h>
#include <TStyle.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveViewer.h>
#include <TEveWindow.h>
#include <TEveScene.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliEveEventManager.h>
#include <AliEveJetPlane.h>
#endif

TEveViewer *gJPView   = 0;
TEveScene  *gJPScene  = 0;

AliEveJetPlane* jetplane()
{
  if (gJPView == 0)
  {
    TEveWindowSlot *slot    = 0;
    TEveBrowser    *browser = gEve->GetBrowser();

    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->MakeCurrent();
    gJPView  = gEve->SpawnNewViewer("JetPlane", "");
    gJPScene = gEve->SpawnNewScene("JetPlane", "Scene holding elements of the jet-plane view.");
    gJPView->AddScene(gJPScene);

    gJPView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
  }

  AliAODEvent* aod = AliEveEventManager::AssertAOD();
  
  // We have event id everywhere now.
  Int_t iev = AliEveEventManager::Instance()->GetEventId();

  gStyle->SetPalette(1, 0);

  AliEveJetPlane* jp = new AliEveJetPlane(iev);
  jp->SetPickable(kTRUE);

  // Read Jets in current event

  TClonesArray* jets = aod->GetJets();
  Int_t njets = jets->GetEntries();
  Info("jetplane", "Event: %5d Number of jets: %5d \n", iev, njets);

  for (Int_t ij = 0; ij < njets; ij++)
  {
    AliAODJet *jet = (AliAODJet*) jets->At(ij);
    jp->AddJet(jet);
  }

  // Read tracks in current event

  TClonesArray* tracks = aod->GetTracks();
  Int_t ntracks = tracks->GetEntries();
  Info("jetplane", "Event: %5d Number of tracks: %5d \n", iev, ntracks);

  for (Int_t ij = 0; ij < ntracks; ij++)
  {
    AliAODTrack* track = (AliAODTrack*) tracks->At(ij);
    jp->AddTrack(track);
  }

  jp->CreateArrows();

  // Render Jet Plane
  gJPScene->AddElement(jp);
  AliEveEventManager::Instance()->RegisterTransient(jp);



    
  gEve->Redraw3D();

  return jp;
}
