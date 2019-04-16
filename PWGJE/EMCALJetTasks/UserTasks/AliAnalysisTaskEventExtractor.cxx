/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TGrid.h>
#include <TFile.h>
#include <TSystem.h>

#include "AliMCEvent.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"

#include "AliVTrack.h"
#include "AliVHeader.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskEventExtractor.h"

/// \cond CLASSIMP
ClassImp(AliEventTree)
/// \endcond

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEventExtractor)
/// \endcond


//________________________________________________________________________
AliEventTree::AliEventTree() : TNamed("EventTree", "EventTree"), fTree(0), fInitialized(kFALSE)
{
  // Create tree buffer: Up to 50000 particles + tracks per event
  fBuffer_Particles_Pt      = new Float_t[50000];
  fBuffer_Particles_Eta     = new Float_t[50000];
  fBuffer_Particles_Phi     = new Float_t[50000];
  fBuffer_Particles_Charge  = new Int_t[50000];
  fBuffer_Particles_Label   = new Int_t[50000];
  fBuffer_Particles_PdgCode = new Int_t[50000];
  fBuffer_Particles_E       = new Float_t[50000];

  fBuffer_Tracks_Pt         = new Float_t[50000];
  fBuffer_Tracks_Eta        = new Float_t[50000];
  fBuffer_Tracks_Phi        = new Float_t[50000];
  fBuffer_Tracks_Charge     = new Int_t[50000];
  fBuffer_Tracks_Label      = new Int_t[50000];
}

//________________________________________________________________________
AliEventTree::AliEventTree(const char* name) : TNamed(name, name), fTree(0), fInitialized(kFALSE)
{
  // Create tree buffer: Up to 50000 particles + tracks per event
  fBuffer_Particles_Pt      = new Float_t[50000];
  fBuffer_Particles_Eta     = new Float_t[50000];
  fBuffer_Particles_Phi     = new Float_t[50000];
  fBuffer_Particles_Charge  = new Int_t[50000];
  fBuffer_Particles_Label   = new Int_t[50000];
  fBuffer_Particles_PdgCode = new Int_t[50000];
  fBuffer_Particles_E       = new Float_t[50000];

  fBuffer_Tracks_Pt         = new Float_t[50000];
  fBuffer_Tracks_Eta        = new Float_t[50000];
  fBuffer_Tracks_Phi        = new Float_t[50000];
  fBuffer_Tracks_Charge     = new Int_t[50000];
  fBuffer_Tracks_Label      = new Int_t[50000];
}

//________________________________________________________________________
void AliEventTree::InitializeTree(Bool_t saveTracks)
{
  // Create the tree with active branches

  if(fInitialized)
    AliFatal("Tree is already initialized.");

  // ### Prepare the event tree
  fTree = new TTree(Form("EventTree_%s", GetName()), "");
  fTree->Branch("Event_Number_Particles",&fBuffer_Event_Number_Particles,"Event_Number_Particles/I");
  fTree->Branch("Event_Number_Tracks",&fBuffer_Event_Number_Tracks,"Event_Number_Tracks/I");
  fTree->Branch("Event_ID",&fBuffer_Event_ID,"Event_ID/L");

  fTree->Branch("Particles_Pt",fBuffer_Particles_Pt,"Particles_Pt[Event_Number_Particles]/F");
  fTree->Branch("Particles_Eta",fBuffer_Particles_Eta,"Particles_Eta[Event_Number_Particles]/F");
  fTree->Branch("Particles_Phi",fBuffer_Particles_Phi,"Particles_Phi[Event_Number_Particles]/F");
  fTree->Branch("Particles_Charge",fBuffer_Particles_Charge,"Particles_Charge[Event_Number_Particles]/S");
  fTree->Branch("Particles_Label",fBuffer_Particles_Label,"Particles_Label[Event_Number_Particles]/I");
  fTree->Branch("Particles_PdgCode",fBuffer_Particles_PdgCode,"Particles_PdgCode[Event_Number_Particles]/I");
  fTree->Branch("Particles_E",fBuffer_Particles_E,"Particles_E[Event_Number_Particles]/F");
  if(saveTracks)
  {
    fTree->Branch("Tracks_Pt",fBuffer_Tracks_Pt,"Tracks_Pt[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Eta",fBuffer_Tracks_Eta,"Tracks_Eta[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Phi",fBuffer_Tracks_Phi,"Tracks_Phi[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Charge",fBuffer_Tracks_Charge,"Tracks_Charge[Event_Number_Tracks]/S");
    fTree->Branch("Tracks_Label",fBuffer_Tracks_Label,"Tracks_Label[Event_Number_Tracks]/I");
  }

  fInitialized = kTRUE;
}

//________________________________________________________________________
Bool_t AliEventTree::AddEventToTree(Long64_t eventID, std::vector<AliVParticle*> particles, std::vector<AliVTrack*> tracks)
{
  if(!fInitialized)
    AliFatal("Tree is not initialized.");


  fBuffer_Event_Number_Particles  = 0;
  fBuffer_Event_Number_Tracks     = 0;
  fBuffer_Event_ID                = eventID;

  // Fill buffer tracks & particles
  for(Int_t i = 0; i < tracks.size(); i++)
  {
    if(!tracks[i]) continue;
    fBuffer_Tracks_Pt[fBuffer_Event_Number_Tracks]     = tracks[i]->Pt();
    fBuffer_Tracks_Eta[fBuffer_Event_Number_Tracks]    = tracks[i]->Eta();
    fBuffer_Tracks_Phi[fBuffer_Event_Number_Tracks]    = tracks[i]->Phi();
    fBuffer_Tracks_Charge[fBuffer_Event_Number_Tracks] = tracks[i]->Charge();
    fBuffer_Tracks_Label[fBuffer_Event_Number_Tracks]  = tracks[i]->GetLabel();
    fBuffer_Event_Number_Tracks++;
  }
  for(Int_t i = 0; i < particles.size(); i++)
  {
    if(!particles[i]) continue;
    fBuffer_Particles_Pt[fBuffer_Event_Number_Particles]       = particles[i]->Pt();
    fBuffer_Particles_Eta[fBuffer_Event_Number_Particles]      = particles[i]->Eta();
    fBuffer_Particles_Phi[fBuffer_Event_Number_Particles]      = particles[i]->Phi();
    fBuffer_Particles_Charge[fBuffer_Event_Number_Particles]   = particles[i]->Charge();
    fBuffer_Particles_Label[fBuffer_Event_Number_Particles]    = particles[i]->GetLabel();
    fBuffer_Particles_PdgCode[fBuffer_Event_Number_Particles]  = particles[i]->PdgCode();
    fBuffer_Particles_E[fBuffer_Event_Number_Particles]        = particles[i]->E();
    fBuffer_Event_Number_Particles++;
  }

  // Add all buffers to tree
  fTree->Fill();

  return kTRUE;
}

//________________________________________________________________________
AliAnalysisTaskEventExtractor::AliAnalysisTaskEventExtractor() :
  AliAnalysisTaskEmcal("AliAnalysisTaskEventExtractor", kTRUE),
  fSaveTracks(kTRUE),
  fEventPercentage(1.0),
  fMCParticleArrayName("mcparticles"),
  fTrackArrayName("tracks"),
  fRandomSeed(0),
  fCustomStartupScript(),
  fTree(0),
  fMultiplicity(0),
  fMCParticleArray(0),
  fTrackArray(0),
  fRandomGenerator(0)
{
  fRandomGenerator = new TRandom3();

  SetMakeGeneralHistograms(kTRUE);
  fTree = new AliEventTree(GetName());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEventExtractor::AliAnalysisTaskEventExtractor(const char *name) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fSaveTracks(kTRUE),
  fEventPercentage(1.0),
  fMCParticleArrayName("mcparticles"),
  fTrackArrayName("tracks"),
  fRandomSeed(0),
  fCustomStartupScript(),
  fTree(0),
  fMultiplicity(0),
  fMCParticleArray(0),
  fTrackArray(0),
  fRandomGenerator(0)
{
  fRandomGenerator = new TRandom3();

  SetMakeGeneralHistograms(kTRUE);
  fTree = new AliEventTree(GetName());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEventExtractor::~AliAnalysisTaskEventExtractor()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEventExtractor::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  fRandomGenerator->SetSeed(fRandomSeed);

  // ### Initialize the event tree (settings must all be given at this stage)
  fTree->InitializeTree(fSaveTracks);

  OpenFile(2);
  PostData(2, fTree->GetTreePointer());
}

//________________________________________________________________________
void AliAnalysisTaskEventExtractor::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  if (fMCParticleArrayName != "")
    fMCParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticleArrayName.Data()));
  if (fTrackArrayName != "")
    fTrackArray      = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrayName.Data()));

  if(!fMCParticleArray)
    AliFatal(Form("MC particle array not valid"));
  if((fTrackArrayName != "") && !fTrackArray)
    AliFatal(Form("Track array not valid, but demanded"));

  // ### Execute shell script at startup
  if(!fCustomStartupScript.IsNull())
  {
    TGrid::Connect("alien://");
    TFile::Cp(fCustomStartupScript.Data(), "./myScript.sh");
    gSystem->Exec("bash ./myScript.sh");
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEventExtractor::Run()
{
  // ################################### EVENT SELECTION
  if(fRandomGenerator->Rndm() >= fEventPercentage)
    return kFALSE;

  // ################################### EVENT PROPERTIES

  // Load vertex if possible
  Long64_t eventID = 0;
  const AliVVertex* myVertex = InputEvent()->GetPrimaryVertex();
  if(!myVertex && MCEvent())
    myVertex = MCEvent()->GetPrimaryVertex();
  Double_t vtx[3] = {0, 0, 0};
  if(myVertex)
  {
    vtx[0] = myVertex->GetX(); vtx[1] = myVertex->GetY(); vtx[2] = myVertex->GetZ();
  }

  // Get event ID from header
  AliVHeader* eventIDHeader = InputEvent()->GetHeader();
  if(eventIDHeader)
    eventID = eventIDHeader->GetEventIdAsLong();

  // Now save particles + tracks to std vectors and add them to the tree
  std::vector<AliVParticle*> particles;
  std::vector<AliVTrack*>    tracks;
  for(Int_t iPart=0; iPart<fMCParticleArray->GetEntries();iPart++)
  {
    AliVParticle* part = (AliVParticle*)fMCParticleArray->At(iPart);
    if(!part) continue;
    if(!part->IsPhysicalPrimary()) continue;
    particles.push_back(part);
  }
  for(Int_t iTrack=0; iTrack<fTrackArray->GetEntries();iTrack++)
  {
    AliVTrack* track = (AliVTrack*)fTrackArray->At(iTrack);
    if(!track) continue;
    tracks.push_back(track);
  }
  fTree->AddEventToTree(eventID, particles, tracks);

  return kTRUE;
}


//########################################################################
// HELPERS
//########################################################################

//________________________________________________________________________
inline void AliAnalysisTaskEventExtractor::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskEventExtractor::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskEventExtractor::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}

//________________________________________________________________________
inline void AliAnalysisTaskEventExtractor::FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add)
{
  TH3* tmpHist = static_cast<TH3*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  if(add)
    tmpHist->Fill(x,y,z,add);
  else
    tmpHist->Fill(x,y,z);
}


//________________________________________________________________________
template <class T> T* AliAnalysisTaskEventExtractor::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEventExtractor::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEventExtractor::AddHistogram3D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, Int_t zBins, Double_t zMin, Double_t zMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax, zBins, zMin, zMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
void AliAnalysisTaskEventExtractor::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

// ### ADDTASK MACRO
//________________________________________________________________________
AliAnalysisTaskEventExtractor* AliAnalysisTaskEventExtractor::AddTaskEventExtractor(TString trackArray, TString particleArray, const char* taskNameSuffix)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     suffix             = "";
  if(taskNameSuffix != 0)
    suffix = taskNameSuffix;

  // ###### Task name
  TString name("AliAnalysisTaskEventExtractor");
  if (suffix != "") {
    name += "_";
    name += suffix;
  }

  // ###### Setup task with default settings
  AliAnalysisTaskEventExtractor* myTask = new AliAnalysisTaskEventExtractor(name);
  myTask->SetNeedEmcalGeom(kFALSE);
  myTask->SetVzRange(-10.,10.);
  mgr->AddTask(myTask);

  // ###### Connect inputs/outputs
  mgr->ConnectInput  (myTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (myTask, 1, mgr->CreateContainer(Form("%s_histos", name.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EventExtractor", mgr->GetCommonFileName())) );
  mgr->ConnectOutput (myTask, 2, mgr->CreateContainer(Form("%s_tree", name.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return myTask;
}
