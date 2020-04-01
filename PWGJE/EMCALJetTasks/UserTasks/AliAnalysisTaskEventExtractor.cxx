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

#include "AliAODTrack.h"
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
  fBuffer_Particles_IsPrimary = new Bool_t[50000];

  fBuffer_Tracks_Pt         = new Float_t[50000];
  fBuffer_Tracks_Eta        = new Float_t[50000];
  fBuffer_Tracks_Phi        = new Float_t[50000];
  fBuffer_Tracks_Charge     = new Int_t[50000];
  fBuffer_Tracks_Label      = new Int_t[50000];
  fBuffer_Tracks_IsGlobal   = new Bool_t[50000];
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
  fBuffer_Particles_IsPrimary = new Bool_t[50000];

  fBuffer_Tracks_Pt         = new Float_t[50000];
  fBuffer_Tracks_Eta        = new Float_t[50000];
  fBuffer_Tracks_Phi        = new Float_t[50000];
  fBuffer_Tracks_Charge     = new Int_t[50000];
  fBuffer_Tracks_Label      = new Int_t[50000];
  fBuffer_Tracks_IsGlobal   = new Bool_t[50000];

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
  fTree->Branch("Event_Vertex_X",&fBuffer_Event_Vertex_X,"Event_Vertex_X/F");
  fTree->Branch("Event_Vertex_Y",&fBuffer_Event_Vertex_Y,"Event_Vertex_Y/F");
  fTree->Branch("Event_Vertex_Z",&fBuffer_Event_Vertex_Z,"Event_Vertex_Z/F");
  fTree->Branch("Event_ID",&fBuffer_Event_ID,"Event_ID/L");

  fTree->Branch("Particles_Pt",fBuffer_Particles_Pt,"Particles_Pt[Event_Number_Particles]/F");
  fTree->Branch("Particles_Eta",fBuffer_Particles_Eta,"Particles_Eta[Event_Number_Particles]/F");
  fTree->Branch("Particles_Phi",fBuffer_Particles_Phi,"Particles_Phi[Event_Number_Particles]/F");
  fTree->Branch("Particles_Charge",fBuffer_Particles_Charge,"Particles_Charge[Event_Number_Particles]/S");
  fTree->Branch("Particles_Label",fBuffer_Particles_Label,"Particles_Label[Event_Number_Particles]/I");
  fTree->Branch("Particles_PdgCode",fBuffer_Particles_PdgCode,"Particles_PdgCode[Event_Number_Particles]/I");
  fTree->Branch("Particles_E",fBuffer_Particles_E,"Particles_E[Event_Number_Particles]/F");
  fTree->Branch("Particles_IsPrimary",fBuffer_Particles_IsPrimary,"Particles_IsPrimary[Event_Number_Particles]/O");

  if(saveTracks)
  {
    fTree->Branch("Tracks_Pt",fBuffer_Tracks_Pt,"Tracks_Pt[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Eta",fBuffer_Tracks_Eta,"Tracks_Eta[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Phi",fBuffer_Tracks_Phi,"Tracks_Phi[Event_Number_Tracks]/F");
    fTree->Branch("Tracks_Charge",fBuffer_Tracks_Charge,"Tracks_Charge[Event_Number_Tracks]/S");
    fTree->Branch("Tracks_Label",fBuffer_Tracks_Label,"Tracks_Label[Event_Number_Tracks]/I");
    fTree->Branch("Tracks_IsGlobal",fBuffer_Tracks_IsGlobal,"Tracks_IsGlobal[Event_Number_Tracks]/O");
  }

  fInitialized = kTRUE;
}

//________________________________________________________________________
Bool_t AliEventTree::AddEventToTree(Long64_t eventID, std::vector<Double_t> vertex, std::vector<AliVParticle*> particles, std::vector<AliVTrack*> tracks)
{
  if(!fInitialized)
    AliFatal("Tree is not initialized.");


  fBuffer_Event_Number_Particles  = 0;
  fBuffer_Event_Number_Tracks     = 0;
  fBuffer_Event_ID                = eventID;
  fBuffer_Event_Vertex_X          = vertex[0];
  fBuffer_Event_Vertex_Y          = vertex[1];
  fBuffer_Event_Vertex_Z          = vertex[2];

  // Fill buffer tracks & particles
  for(Size_t i = 0; i < tracks.size(); i++)
  {
    if(!tracks[i]) continue;
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(tracks[i]);
    fBuffer_Tracks_Pt[fBuffer_Event_Number_Tracks]     = tracks[i]->Pt();
    fBuffer_Tracks_Eta[fBuffer_Event_Number_Tracks]    = tracks[i]->Eta();
    fBuffer_Tracks_Phi[fBuffer_Event_Number_Tracks]    = tracks[i]->Phi();
    fBuffer_Tracks_Charge[fBuffer_Event_Number_Tracks] = tracks[i]->Charge();
    fBuffer_Tracks_Label[fBuffer_Event_Number_Tracks]  = tracks[i]->GetLabel();
    fBuffer_Tracks_IsGlobal[fBuffer_Event_Number_Tracks] = (aodTrack) ? aodTrack->IsGlobalConstrained() : kFALSE;
    fBuffer_Event_Number_Tracks++;
  }
  for(Size_t i = 0; i < particles.size(); i++)
  {
    if(!particles[i]) continue;
    fBuffer_Particles_Pt[fBuffer_Event_Number_Particles]       = particles[i]->Pt();
    fBuffer_Particles_Eta[fBuffer_Event_Number_Particles]      = particles[i]->Eta();
    fBuffer_Particles_Phi[fBuffer_Event_Number_Particles]      = particles[i]->Phi();
    fBuffer_Particles_Charge[fBuffer_Event_Number_Particles]   = particles[i]->Charge();
    fBuffer_Particles_Label[fBuffer_Event_Number_Particles]    = particles[i]->GetLabel();
    fBuffer_Particles_PdgCode[fBuffer_Event_Number_Particles]  = particles[i]->PdgCode();
    fBuffer_Particles_E[fBuffer_Event_Number_Particles]        = particles[i]->E();
    fBuffer_Particles_IsPrimary[fBuffer_Event_Number_Particles]= particles[i]->IsPhysicalPrimary();
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
  fRandomSeed(0),
  fCustomStartupScript(),
  fTree(0),
  fMultiplicity(0),
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
  fRandomSeed(0),
  fCustomStartupScript(),
  fTree(0),
  fMultiplicity(0),
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

  if(!GetParticleContainer(0))
    AliFatal(Form("MC particle container not valid"));
  if(!GetTrackContainer(1) && fSaveTracks)
    AliFatal(Form("Track container not valid, but demanded"));

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
  std::vector<Double_t> vtx;
  if(myVertex)
  {
    vtx.push_back(myVertex->GetX()); vtx.push_back(myVertex->GetY()); vtx.push_back(myVertex->GetZ());
  }

  // Get event ID from header
  AliVHeader* eventIDHeader = InputEvent()->GetHeader();
  if(eventIDHeader)
    eventID = eventIDHeader->GetEventIdAsLong();

  // Now save particles + tracks to std vectors and add them to the tree
  std::vector<AliVParticle*> particles;
  std::vector<AliVTrack*>    tracks;

  for (auto part : GetParticleContainer(0)->accepted()) {
    particles.push_back(part);
  }

  if(GetParticleContainer(1) && fSaveTracks)
  {
    for (auto track : GetTrackContainer(1)->accepted()) {
      tracks.push_back(track);
    }
  }
  fTree->AddEventToTree(eventID, vtx, particles, tracks);

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

  myTask->AddMCParticleContainer(particleArray);
  myTask->AddTrackContainer(trackArray);

  // Configure the eta limits for the particle and track containers.
  // They are configured here so we could modify them later if so desired.
  // The default values are selected because the eta distribution at particle level
  // is quite similar to the detector level eta distribution. However, we add some
  // additional margin to be certain that we aren't artificially distorting the particle
  // level distribution.
  double etaMin = -2.0;
  double etaMax = 2.0;
  myTask->GetParticleContainer(0)->SetEtaLimits(1.5*etaMin, 1.5*etaMax);
  if (myTask->GetTrackContainer(1))
  {
    myTask->GetTrackContainer(1)->SetEtaLimits(etaMin, etaMax);
  }

  mgr->AddTask(myTask);

  // ###### Connect inputs/outputs
  mgr->ConnectInput  (myTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (myTask, 1, mgr->CreateContainer(Form("%s_histos", name.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EventExtractor", mgr->GetCommonFileName())) );
  mgr->ConnectOutput (myTask, 2, mgr->CreateContainer(Form("%s_tree", name.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return myTask;
}
