// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file MUON_displayData.C
///
/// \author B. Vulpescu, LPC

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <Rtypes.h>
#include <TTree.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveTrack.h>

#include <AliLog.h>
#include <AliMpSegmentation.h>
#include <AliMpDDLStore.h>
#include <AliMpCDB.h>
#include <AliMUONTrackExtrap.h>
#include <AliRunLoader.h>
#include <AliESDEvent.h>
#include <AliESDMuonTrack.h>
#include <AliEveEventManager.h>
#include <AliEveMUONData.h>
#include <AliEveMUONChamber.h>
#include <AliEveMUONTrack.h>
#endif

class AliEveMUONData;
class AliEveEventManager;

AliEveMUONData     *g_muon_data       = 0;

Int_t  g_currentEvent = -1;
Bool_t g_fromRaw      = kFALSE;

void MUON_ESD_tracks();

void MUON_displayData(Bool_t fromRaw = kFALSE, Bool_t showTracks = kTRUE, Bool_t clustersFromESD = kTRUE)
{
  //
  // display from real data, eventually with recreated digits
  // tracks: ESD
  // 

  if (!AliMpSegmentation::Instance()) AliMpCDB::LoadMpSegmentation();
  if (!AliMpDDLStore::Instance())     AliMpCDB::LoadDDLStore();

  // set the magnetic field for track extrapolations
  AliEveEventManager::AssertMagField();
  AliMUONTrackExtrap::SetField();

  TTree* dt = 0;
  TTree* ct = 0;
  TTree* ht = 0;

  if (AliEveEventManager::Instance() == 0) {
    printf("No alieve event: use alieve_init(...) \n");
    return;
  }

  if (g_currentEvent == AliEveEventManager::Instance()->GetEventId()) {
    if (g_fromRaw == fromRaw) {
      printf("Same event... \n");
      return;
    } else {
      if (g_fromRaw) {
	printf("Same event with digits.\n");
	AliEveEventManager::Instance()->GotoEvent(g_currentEvent);
      } else {
	printf("Same event with raw.\n");
	AliEveEventManager::Instance()->GotoEvent(g_currentEvent);
      }
    }
  }

  g_fromRaw = fromRaw;

  TString dataPath = TString(AliEveEventManager::Instance()->GetTitle());
  dataPath.Append("/raw.root");

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  g_muon_data = new AliEveMUONData;

  if (!fromRaw) {
    rl->LoadDigits("MUON");
    dt = rl->GetTreeD("MUON", false);
    if (dt == 0) {
      AliInfoGeneral("MUON_displayData.C", "No digits produced!");
    } else {
      AliInfoGeneral("MUON_displayData.C", "With aliroot digits!");
      g_muon_data->LoadDigits(dt);
    }
  } else {
    if (gSystem->AccessPathName(dataPath.Data(),kFileExists)) {
      AliInfoGeneral("MUON_displayData.C", "No raw data produced!");
    } else {
      AliInfoGeneral("MUON_displayData.C", "With raw digits!");
      g_muon_data->LoadRaw(dataPath.Data());
    }
  }

  TString esdDataPath = TString(AliEveEventManager::Instance()->GetTitle());
  esdDataPath.Append("/AliESDs.root");
  if (clustersFromESD) {
    g_muon_data->LoadRecPointsFromESD(esdDataPath.Data());
  } else {
    rl->LoadRecPoints("MUON");
    ct = rl->GetTreeR("MUON", false);
    g_muon_data->LoadRecPoints(ct);
  }
  
  g_currentEvent = AliEveEventManager::Instance()->GetEventId();

  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("MUONChambers");
  l->SetTitle("MUON chambers");
  l->SetMainColor(2);
  gEve->AddElement(l);

  for (Int_t ic = 0; ic < 14; ic++)
  {
    AliEveMUONChamber* mucha = new AliEveMUONChamber(ic);

    mucha->SetFrameColor(2);
    mucha->SetChamberID(ic);

    mucha->SetDataSource(g_muon_data);

    l->AddElement(mucha);
  }

  if (showTracks) {
    MUON_ESD_tracks();
  }

  gEve->Redraw3D(kTRUE);
  gEve->EnableRedraw();
}

//______________________________________________________________________________
void MUON_ESD_tracks()
{
  AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

  // TEveTrackList* lt = new TEveTrackList("ESD-Tracks");
  AliEveMUONTrackList* lt = new AliEveMUONTrackList("ESD-Tracks");
  lt->SetMainColor(6);
  //lt->SetMUON();

  gEve->AddElement(lt);

  AliESDMuonTrack *mt;
  TEveRecTrack rt;
  Int_t nMuonTracks = esd->GetNumberOfMuonTracks();
  Int_t nTrack = 0;
  for (Int_t n = 0; n < nMuonTracks; n++)
  {
    mt = esd->GetMuonTrack(n);

    if (mt->GetNHit() == 0) continue;
    nTrack++;

    rt.fLabel = n;

    AliEveMUONTrack* track = new AliEveMUONTrack(&rt, lt->GetPropagator());

    track->MakeESDTrack(mt);

    lt->AddElement(track);
  }
  lt->HackMomentumLimits();
}

