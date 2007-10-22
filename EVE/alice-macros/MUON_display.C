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
 
  if (!AliMpSegmentation::Instance()) AliMpCDB::LoadMpSegmentation();  
  if (!AliMpDDLStore::Instance())     AliMpCDB::LoadDDLStore();

  TTree* dt = 0;
  TTree* ct = 0;
  TTree* ht = 0;

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
  dataPath.Append("/rawmuon.root");

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  g_muon_data = new Alieve::MUONData;
  
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
  
  rl->LoadRecPoints("MUON");
  ct = rl->GetTreeR("MUON", false);
  g_muon_data->LoadRecPoints(ct);
  
  rl->LoadHits("MUON");
  ht = rl->GetTreeH("MUON", false);
  g_muon_data->LoadHits(ht);
  
  g_muon_last_event = Alieve::gEvent;
  
  g_currentEvent = g_muon_last_event->GetEventId();
  
  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();
  
  Reve::RenderElementList* l = new Reve::RenderElementList("MUONChambers");
  l->SetTitle("MUON chambers");
  l->SetMainColor(Color_t(2));
  gReve->AddRenderElement(l);
  
  for (Int_t ic = 0; ic < 14; ic++) {

    Alieve::MUONChamber* mucha = new Alieve::MUONChamber(ic);
    
    mucha->SetFrameColor(2);
    mucha->SetChamberID(ic);
    
    mucha->SetDataSource(g_muon_data);

    gReve->AddRenderElement(mucha, l);

  }

  if (showTracks) {
    MUON_tracks();
    MUON_trigger_tracks();
    MUON_ESD_tracks();
    MUON_Ref_tracks();
    MUON_MC_tracks();
  }

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

  TClonesArray *tracks = 0;
  tt->SetBranchAddress("MUONTrack",&tracks);
  tt->GetEntry(0);

  Int_t ntracks = tracks->GetEntriesFast();
  //printf("Found %d tracks. \n",ntracks);

  Reve::TrackList* lt = new Reve::TrackList("M-Tracks"); 
  lt->SetMainColor(Color_t(6));
  //lt->SetMUON();  

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

    rt.label = n;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeMUONTrack(mt);

    gReve->AddRenderElement(track, lt);

  }

  rl->UnloadTracks("MUON");

}

//_____________________________________________________________________________
void MUON_trigger_tracks() {

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadTracks("MUON");
  TTree* tt = rl->GetTreeT("MUON", false);

  TClonesArray *tracks = 0;
  tt->SetBranchAddress("MUONTriggerTrack",&tracks);
  tt->GetEntry(0);

  Int_t ntracks = tracks->GetEntriesFast();
  //printf("Found %d tracks. \n",ntracks);

  Reve::TrackList* lt = new Reve::TrackList("MT-Tracks"); 
  lt->SetMainColor(Color_t(4));
  //lt->SetMUON();  

  gReve->AddRenderElement(lt);

  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);

  Float_t xRec, xRec0;
  Float_t yRec, yRec0;
  Float_t zRec, zRec0;
  
  Float_t zg[4] = { -1603.5, -1620.5, -1703.5, -1720.5 };

  AliMUONTriggerTrack *mt;  
  Reve::RecTrack  rt;
  Int_t count;
  for (Int_t n = 0; n < ntracks; n++) {
    
    count = 0;

    mt = (AliMUONTriggerTrack*) tracks->At(n);

    rt.label = n;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeMUONTriggerTrack(mt);

    gReve->AddRenderElement(track, lt);

  }

  rl->UnloadTracks("MUON");

}

//_____________________________________________________________________________
void MUON_ESD_tracks() {

  AliESDEvent* esd = Alieve::Event::AssertESD();

  Reve::TrackList* lt = new Reve::TrackList("ESD-Tracks"); 
  lt->SetMainColor(Color_t(6));
  //lt->SetMUON();

  gReve->AddRenderElement(lt);

  AliESDMuonTrack *mt;
  Reve::RecTrack rt;
  Int_t nMuonTracks = esd->GetNumberOfMuonTracks();
  for (Int_t n = 0; n < nMuonTracks; n++) {

    mt = esd->GetMuonTrack(n);

    rt.label = n;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeESDTrack(mt);

    gReve->AddRenderElement(track, lt);

  }

}

//_____________________________________________________________________________
void MUON_Ref_tracks() {

  TString dataPath = TString(Alieve::gEvent->GetTitle());
  dataPath.Append("/galice.root");

  AliMUONRecoCheck recoCheck(dataPath.Data(),dataPath.Data());
  AliMUONVTrackStore* trackRefStore = recoCheck.ReconstructibleTracks(Alieve::gEvent->GetEventId());
  TIter next(trackRefStore->CreateIterator());
  AliMUONTrack* trackRef;
  
  Reve::TrackList* lt = new Reve::TrackList("Ref-Tracks"); 
  lt->SetMainColor(Color_t(6));

  gReve->AddRenderElement(lt);

  Reve::RecTrack rt;
  Int_t i = 0;  
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {

    rt.label = i++;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeRefTrack(trackRef);

    gReve->AddRenderElement(track, lt);

  }

}

//_____________________________________________________________________________
void MUON_MC_tracks() {

  Double_t RADDEG = 180.0/TMath::Pi();

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  Int_t nPrimary = stack->GetNprimary();
  Int_t nTracks  = stack->GetNtrack();

  Reve::TrackList* lt = new Reve::TrackList("MC-Tracks"); 
  lt->SetMainColor(Color_t(6));
  //lt->SetMUON();

  gReve->AddRenderElement(lt);

  Int_t pdgCode;
  TParticle *part;
  Reve::RecTrack rt;

  Int_t nHitTracks = g_muon_data->GetNTrackList();
  Int_t index;
  for (Int_t i = 0; i < nHitTracks; i++) {

    index = g_muon_data->GetTrack(i);
    if (index >= nTracks) {
      cout << "Hit track index larger than number in stack!" << endl;
      continue;
    }

    part = stack->Particle(index);
    if (part->P() < 0.001) continue;  // skip momenta < 1.0 MeV/c
    rt.label = i;

    Alieve::MUONTrack* track = new Alieve::MUONTrack(&rt, lt->GetRnrStyle());

    track->MakeMCTrack(part);

    gReve->AddRenderElement(track, lt);

  }

  rl->UnloadKinematics();

}

