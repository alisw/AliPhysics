// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
AliEveJetPlane* jetplane(Int_t iev)
{	
  TFile* f       = new TFile("aod.root");
  TTree* treeAOD = (TTree*) f->Get("AOD");
  AliAODEvent* aod = new AliAODEvent();
  aod->ReadFromTree(treeAOD);
  treeAOD->GetEntry(iev);

  AliEveJetPlane* jp = new AliEveJetPlane(iev);

  // Read Jets in current event

  TClonesArray* jets = aod->GetJets();
  Int_t njets = jets->GetEntries();
  printf("AliEveEventManager: %5d Number of jets: %5d \n", iev, njets);

  for (Int_t ij = 0; ij < njets; ij++)
  {
    AliAODJet jet = (AliAODJet) jets->At(ij);
    jp->AddJet(jet);
  }

  // Read tracks in current event

  TClonesArray* tracks = aod->GetTracks();
  Int_t ntracks = tracks->GetEntries();
  printf("AliEveEventManager: %5d Number of tracks: %5d \n", iev, ntracks);

  for (Int_t ij = 0; ij < ntracks; ij++)
  {
    AliAODTrack track = (AliAODTrack) tracks->At(ij);
    jp->AddTrack(track);
  }

  // Render Jet Plane
  gStyle->SetPalette(1, 0);
  gEve->AddElement(jp);
  gEve->Redraw3D();

  return jp;
}
