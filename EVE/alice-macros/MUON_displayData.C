// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class AliEveMUONData;
class AliEveEventManager;

AliEveMUONData     *g_muon_data       = 0;

Int_t  g_currentEvent = -1;
Bool_t g_fromRaw      = kFALSE;

AliMagFMaps *g_field = 0;

void MUON_displayData(Bool_t fromRaw = kFALSE, Bool_t showTracks = kTRUE, Bool_t clustersFromESD = kTRUE)
{
  //
  // display from real data, eventually with recreated digits
  // tracks: ESD
  // 

  if (!AliMpSegmentation::Instance()) AliMpCDB::LoadMpSegmentation();
  if (!AliMpDDLStore::Instance())     AliMpCDB::LoadDDLStore();

  if (g_field == 0) {
    printf("Loading field map...\n");
    g_field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
    AliTracker::SetFieldMap(g_field, kFALSE);
    AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());
  }

  TTree* dt = 0;
  TTree* ct = 0;
  TTree* ht = 0;

  if (AliEveEventManager::GetMaster() == 0) {
    printf("No alieve event: use alieve_init(...) \n");
    return;
  }

  if (g_currentEvent == AliEveEventManager::GetMaster()->GetEventId()) {
    if (g_fromRaw == fromRaw) {
      printf("Same event... \n");
      return;
    } else {
      if (g_fromRaw) {
	printf("Same event with digits.\n");
	AliEveEventManager::GetMaster()->GotoEvent(g_currentEvent);
      } else {
	printf("Same event with raw.\n");
	AliEveEventManager::GetMaster()->GotoEvent(g_currentEvent);
      }
    }
  }

  g_fromRaw = fromRaw;

  TString dataPath = TString(AliEveEventManager::GetMaster()->GetTitle());
  dataPath.Append("/raw.root");

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  g_muon_data = new AliEveMUONData;

  if (!fromRaw) {
    rl->LoadDigits("MUON");
    dt = rl->GetTreeD("MUON", false);
    if (dt == 0) {
      cout << "No digits produced!" << endl;
    } else {
      cout << "With aliroot digits!" << endl;
      g_muon_data->LoadDigits(dt);
    }
  } else {
    if (gSystem->AccessPathName(dataPath.Data(),kFileExists)) {
      cout << "No raw data produced!" << endl;
    } else {
      cout << "With raw digits!" << endl;
      g_muon_data->LoadRaw(dataPath.Data());
    }
  }

  TString esdDataPath = TString(AliEveEventManager::GetMaster()->GetTitle());
  esdDataPath.Append("/AliESDs.root");
  if (clustersFromESD) {
    g_muon_data->LoadRecPointsFromESD(esdDataPath.Data());
  } else {
    rl->LoadRecPoints("MUON");
    ct = rl->GetTreeR("MUON", false);
    g_muon_data->LoadRecPoints(ct);
  }
  
  g_currentEvent = AliEveEventManager::GetMaster()->GetEventId();

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
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  // TEveTrackList* lt = new TEveTrackList("ESD-Tracks");
  AliEveMUONTrackList* lt = new TEveTrackList("ESD-Tracks");
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

