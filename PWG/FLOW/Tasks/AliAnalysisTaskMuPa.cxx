/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
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

/************************************** 
*   TBI add description eventually    * 
**************************************/ 
  
#include <Riostream.h>
#include <AliAnalysisTaskMuPa.h>
#include <AliLog.h>
#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include "AliMultSelection.h"
#include <TFile.h>
#include <unistd.h>

#include <unistd.h>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMuPa)

//================================================================================================================

AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(const char *name): 
 AliAnalysisTaskSE(name), 
 fBaseList(NULL),
 fBasePro(NULL),
 fFillQAHistograms(kFALSE),
 fTerminateAfterQA(kFALSE),

 // QA:
 fQAList(NULL),
 fQAPro(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 fQAFilterBits(NULL),

 // Control event histograms:
 fControlEventHistogramsList(NULL),
 fControlEventHistogramsPro(NULL),
 fMultiplicity(0.),
 fMultiplicityHist(NULL),
 fSelectedTracks(0),
 fSelectedTracksHist(NULL),
 fCentrality(0.),

 // Control particle histograms:
 fControlParticleHistogramsList(NULL),
 fControlParticleHistogramsPro(NULL),
 fGlobalTracksAOD(NULL),
 fFilterBit(128), // defaulted to TPC-only

 // Q-vectors:
 fQvectorList(NULL),
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kTRUE),
 fMaxHarmonic(6),
 fMaxCorrelator(8),
 // Correlations:
 fCorrelationsList(NULL),        
 fCorrelationsFlagsPro(NULL), 
 fCalculateCorrelations(kTRUE),

 // Nested loops:
 fNestedLoopsList(NULL),        
 fNestedLoopsFlagsPro(NULL), 
 fCalculateNestedLoops(kFALSE),

 // Final results:
 fFinalResultsList(NULL),

 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateFile(NULL),
 fMaxNumberOfEvents(-44),
 fBailOutFile(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(const char *name)");

  // Base list:
  fBaseList = new TList();
  fBaseList->SetName("outputMuPaAnalysis");
  fBaseList->SetOwner(kTRUE);

  // Initialize all non-built in types:
  this->InitializeNonBuiltInTypes();

  // Initialize all arrays:
  this->InitializeArrays();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reser  TList *fQAList;   // base list to hold all QA output object
  // Output slot #1 writes into a TList container

  DefineOutput(1, TList::Class());  

} // AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(const char *name): 

//================================================================================================================

AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(): 
 AliAnalysisTaskSE(),
 fBaseList(NULL),
 fBasePro(NULL),
 fFillQAHistograms(kFALSE),
 fTerminateAfterQA(kFALSE),

 // QA:
 fQAList(NULL),
 fQAPro(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 fQAFilterBits(NULL),

 // Control event histograms:
 fControlEventHistogramsList(NULL),
 fControlEventHistogramsPro(NULL),
 fMultiplicity(0.),
 fMultiplicityHist(NULL),
 fSelectedTracks(0),
 fSelectedTracksHist(NULL),
 fCentrality(0.),

 // Control particle histograms:
 fControlParticleHistogramsList(NULL),
 fControlParticleHistogramsPro(NULL),
 fGlobalTracksAOD(NULL),
 fFilterBit(-44),

 // Q-vectors:
 fQvectorList(NULL),
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kTRUE),
 fMaxHarmonic(6),
 fMaxCorrelator(8),
 // Correlations:
 fCorrelationsList(NULL),        
 fCorrelationsFlagsPro(NULL), 
 fCalculateCorrelations(kTRUE),

 // Nested loops:
 fNestedLoopsList(NULL),        
 fNestedLoopsFlagsPro(NULL), 
 fCalculateNestedLoops(kFALSE),

 // Final results:
 fFinalResultsList(NULL),

 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateFile(NULL),
 fMaxNumberOfEvents(-44),
 fBailOutFile(NULL)
{
  // Dummy constructor.
 
  // Important: arrays have to be initialized also in the dummy constructor. 

  // Initialize all non-built in types:
  this->InitializeNonBuiltInTypes();

  // Initialize all arrays:
  this->InitializeArrays();

  AliDebug(2,"AliAnalysisTaskMuPa::AliAnalysisTaskMuPa()");

} // AliAnalysisTaskMuPa::AliAnalysisTaskMuPa():

//================================================================================================================

AliAnalysisTaskMuPa::~AliAnalysisTaskMuPa()
{
 // Destructor.

 if(fBaseList) delete fBaseList;

 // if(fGlobalTracksAOD) delete fGlobalTracksAOD;
  
} // AliAnalysisTaskMuPa::~AliAnalysisTaskMuPa()

//================================================================================================================

void AliAnalysisTaskMuPa::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Check before bookings if all the values user has provided via setters make sense;
 // b) Trick to avoid name clashes, part 1;
 // c) Book base profile;
 // d) Book and nest all lists;
 // e) Book all objects;
 // *) Trick to avoid name clashes, part 2.
 
 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Check before bookings if all the values user has provided via setters make sense:
 this->InsanityChecks();
 
 // b) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // c) Book base profile:€5
 this->BookBaseProfile();

 // d) Book and nest all lists:
 this->BookAndNestAllLists();

 // e) Book all objects:
 this->BookQAHistograms();
 this->BookControlEventHistograms();
 this->BookControlParticleHistograms();
 this->BookQvectorHistograms();
 this->BookCorrelationsHistograms();
 this->BookNestedLoopsHistograms();
 this->BookFinalResultsHistograms();

 // TBI 20210513 unclassified:
 fGlobalTracksAOD = new TExMap();

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fBaseList);

} // void AliAnalysisTaskMuPa::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMuPa::UserExec(Option_t *) 
{
 // Main loop (called for each event).
 // a) Get pointer to AOD event;
 // b) QA;
 // c) Filter from this event only what is needed, and store that info in local data members;
 // d) Fill control event histograms before cuts;
 // e) Event cuts;
 // f) Fill control event histograms after cuts;
 // g) Start analysis over AODs;
 // h) Fill e-b-e quantities;
 // i) Calculate correlations;
 // j) Calculate nested loops;
 // k) Reset event-by-event objects;
 // *) PostData.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}

 // c) Filter out "normal global" tracks for default analysis and cut on their number:
 //    'TPC-only' tracks and 'global constrained to vertex' come with negative ID, and are therefore not stored in fGlobalTracksAOD
 this->GlobalTracksAOD(aAOD); // [0] stands for default analysis
 if(0 == fGlobalTracksAOD->GetSize()) return; // yes, go to next event TBI 20210513 re-think this line, perhaps add some further check

 // b) QA:
 if(fFillQAHistograms){QA(aAOD);}
 if(fTerminateAfterQA){return;} // TBI 20210513 this can be optimized further, i.e. to force jump immediately to Terminate(...)

 // c) Filter from this event only what is needed, and store that info in local data members:
 this->FilterEvent(aAOD);

 // d) Fill control event histograms before cuts:
 this->FillControlEventHistograms(aAOD,BEFORE,RECO);
 
 // e) Event cuts:
 if(!SurvivesEventCuts(aAOD)){return;}

 // f) Fill control event histograms after cuts:
 this->FillControlEventHistograms(aAOD,AFTER,RECO);

 // g) Start analysis over AODs:
 Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
 Double_t dPt = 0., wPt = 1.; // transverse momentum and corresponding pT weight
 Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
 //Double_t dRapidity = 0., wRapidty = 1.; // rapidity and corresponding rapidity weight TBI_20200612 enable eventually
 Double_t wToPowerP = 1.; // weight raised to power p
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 fSelectedTracks = 0; // counter for tracks which survived all cuts, and which are added to Q-vectors
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to "a track" (i.e. any track)
  if(!aTrack){continue;}
  
  if(!aTrack->TestFilterBit(fFilterBit)) continue; // TBI 20210514 this one has to be here, to avoid double or triple counting. Is it better to filter out all tracks in TArrayD?

  // Particle histograms and track cuts:
  this->FillControlParticleHistograms(aTrack,0,RECO);
  if(!SurvivesParticleCuts(aTrack)){continue;}
  if(fSelectedTracks >= fSelectedTracksCuts[1]){break;}
  this->FillControlParticleHistograms(aTrack,1,RECO);

  // Finally, calculate Q-vectors:
  dPhi = aTrack->Phi();
  dPt = aTrack->Pt(); dPt += 0.; // TBI 20210515 shut down the compiler warnings temporarily
  dEta = aTrack->Eta(); dEta += 0.; // TBI 20210515 shut down the compiler warnings temporarily

  for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  {
   for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
   {
    //if(fUseWeights[0]||fUseWeights[1]||fUseWeights[2]||fUseWeights[3]){wToPowerP = pow(wPhi*wPt*wEta,wp);} 
    //if(fUseWeights[0]||fUseWeights[1]||fUseWeights[2]||fUseWeights[3]){wToPowerP = pow(wPhi*wPt*wEta*wRapidity,wp);} // TBI_20200612 same as above, just taking also wRapidity
    fQvector[h][wp] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));    
   } // for(Int_t wp=0;wp<fMaxCorrelator+1;wp++)
  } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)   

  // Nested loops containers:
  if(fCalculateNestedLoops)
  {
   if(ftaNestedLoops[0]){ftaNestedLoops[0]->AddAt(dPhi,fSelectedTracks);} 
   if(ftaNestedLoops[1]){ftaNestedLoops[1]->AddAt(wPhi*wPt*wEta,fSelectedTracks);} 
  }

  fSelectedTracks++; // number of selected tracks, i.e. number of tracks added to Q-vector
  fMultiplicity += wPhi*wPt*wEta; // this is really the multiplicity. Only if weights are unit, fMultiplicity = fSelectedTracks

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks

 // TBI_20210515 I am a bit incosistent here, because control histogram have been filled nevertheles.
 /*
 if(fSelectedTracks < fSelectedTracksCuts[0])
 {
  this->ResetEventByEventQuantities();
  return;
 }
 */  

 // h) Fill e-b-e quantities:
 fMultiplicityHist->Fill(fMultiplicity);
 fSelectedTracksHist->Fill(fSelectedTracks);

 // i) Calculate correlations:
 if(fCalculateCorrelations){this->CalculateCorrelations();}

 // j) Calculate nested loops:
 if(fCalculateNestedLoops){this->CalculateNestedLoops();}

 // k) Reset event-by-event objects:
 this->ResetEventByEventQuantities();

 // *) PostData:
 PostData(1,fBaseList);

 // *) Online monitoring:
 cout<<"\n\033[1;32m"<<Form("INFO: Done with event #%d",(Int_t)fSelectedTracksHist->GetEntries())<<"\033[0m\n"<<endl; 
 if(fOnlineMonitoring){this->OnlineMonitoring();}

} // void AliAnalysisTaskMuPa::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskMuPa::Terminate(Option_t *)
{
 // Accessing the merged output list.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 fBaseList = (TList*)GetOutputData(1);
 if(!fBaseList){exit(1);}

 // Do some calculation in offline mode here:

 if(fCalculateNestedLoops){this->ComparisonNestedLoopsVsCorrelations();}

} // end of void AliAnalysisTaskMuPa::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskMuPa::ResetEventByEventQuantities()
{
 // Reset all global event-by-event quantities here:

 // a) Multiplicity;
 // b) Centrality;
 // c) Q-vector;
 // d) Reset ebe containers for nested loops.
 
 // a) Multiplicity:
 fMultiplicity = 0.;
 fSelectedTracks = 0;

 // b) Centrality:
 fCentrality = 0.;

 // c) Q-vector:
 if(fCalculateQvector)
 {
  for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  {
   for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
   {
    fQvector[h][wp] = TComplex(0.,0.); 
   }
  }
 } // if(fCalculateQvector)

 // d) Reset ebe containers for nested loops:
 if(fCalculateNestedLoops)
 {
  if(ftaNestedLoops[0]){ftaNestedLoops[0]->Reset();} 
  if(ftaNestedLoops[1]){ftaNestedLoops[1]->Reset();}  
 }

} // void AliAnalysisTaskMuPa::ResetEventByEventQuantities()

//================================================================================================================

void AliAnalysisTaskMuPa::OnlineMonitoring()
{
 // Per request, do some online monitoring.

 // a) Update regularly the output file. Only the events which survive cuts are counted. # of events are entries in fSelectedTracksHist;
 // b) Bail out after specified number of events.

 // a) Update regularly the output file:
 if(fUpdateOutputFile)
 {
  if(!fSelectedTracksHist){cout<<__LINE__<<endl;exit(1);} // TBI 20210515 is this histogram available also if I am doing only QA?
  Int_t currentEventNo = (Int_t)fSelectedTracksHist->GetEntries();
  if(0 == currentEventNo % fUpdateFrequency)
  {
   //cout<<Form("nEvts: %d",currentEventNo)<<endl;
   cout<<Form("\nPer request, updating after %d events the file %s .\n",currentEventNo,fUpdateFile->Data())<<endl;
   sleep(2);
   TFile *f = new TFile(fUpdateFile->Data(),"recreate");
   fBaseList->Write(fBaseList->GetName(),TObject::kSingleKey);
   f->Close();
  }
 } // if(fUpdateOutputFile)

 // b) Bail out after specified number of events:
 if(fMaxNumberOfEvents > 0)
 {
  if(!fSelectedTracksHist){cout<<__LINE__<<endl;exit(1);} // TBI 20210515 is this histogram available also if I am doing only QA?
  Int_t currentEventNo = (Int_t)fSelectedTracksHist->GetEntries();
  if(fMaxNumberOfEvents == currentEventNo)
  {
   cout<<Form("\nPer request, bailing out after %d events in the file %s .\n",fMaxNumberOfEvents,fBailOutFile->Data())<<endl;
   sleep(2);
   TFile *f = new TFile(fBailOutFile->Data(),"recreate");
   fBaseList->Write(fBaseList->GetName(),TObject::kSingleKey);
   f->Close();
   exit(1);
  }
 } // if(fMaxNumberOfEvents > 0)

} // void AliAnalysisTaskMuPa::OnlineMonitoring()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeNonBuiltInTypes()
{
 // Initialize all data members which are not built-in types in this method.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;
 
 fCentralityEstimator = TString("V0M"); // by default, we use V0M as centrality estimator. Can be customized with task->SetCentralityEstimator("V0M") 
 fPeriod = TString("not set"); // can be customized with e.g. task->SetPeriod("LHC10h");
 fAODNumber = TString("not set"); // can be customized with e.g. task->SetAODNumber("AOD160");

} // void AliAnalysisTaskMuPa::InitializeNonBuiltInTypes()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 this->InitializeArraysForQAHistograms();
 this->InitializeArraysForControlEventHistograms();
 this->InitializeArraysForControlParticleHistograms();
 this->InitializeArraysForQvectors();
 this->InitializeArraysForCorrelationsHistograms();
 this->InitializeArraysForNestedLoopsHistograms();
 this->InitializeArraysForCommonLabels();

} // void AliAnalysisTaskMuPa::InitializeArrays()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForQvectors()
{
 // TBI 20210514 finalize documentation

 // Initialize all arrays for Q-vectors:
 for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++) 
 {
  for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
  {
   fQvector[h][wp] = TComplex(0.,0.);
  }
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForQvectors()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForNestedLoopsHistograms()
{
 // TBI 20210515 finalize documentation

 // Initialize all arrays for nested loops:
 for(Int_t k=0;k<4;k++) // order [2p=0,4p=1,6p=2,8p=3]
  { 
   for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
   {
    for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
    {
     fNestedLoopsPro[k][n][v] = NULL;
    } // for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality] 
   } // for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
  } // for(Int_t n=0;n<6;n++) // harmonics [n=1,n=2,...,n=6]

 ftaNestedLoops[0] = NULL;
 ftaNestedLoops[1] = NULL;

} // void AliAnalysisTaskMuPa::InitializeArraysForNestedLoopsHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForCorrelationsHistograms()
{
 // TBI 20210514 finalize documentation

 // Initialize all arrays for correlations (and nested loops):
 for(Int_t k=0;k<4;k++) // order [2p=0,4p=1,6p=2,8p=3]
  { 
   for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
   {
    for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
    {
     fCorrelationsPro[k][n][v] = NULL;
     //fNestedLoopsPro[k][n][v] = NULL;
    } // for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality] 
   } // for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
  } // for(Int_t n=0;n<6;n++) // harmonics [n=1,n=2,...,n=6]

} // void AliAnalysisTaskMuPa::InitializeArraysForCorrelationsHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForQAHistograms()
{
 // Initialize all arrays for QA histograms.

 // a) Centrality;
 // b) Kinematics for specified filter bits.

 // a) Centrality:
 for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
 {
  fQACentralityHist[ce1] = NULL;
  for(Int_t ce2=0;ce2<gCentralityEstimators;ce2++)
  {
   fQACentralityCorrHist[ce1][ce2] = NULL;
  }
 }
 
 // b) Kinematics for specified filter bits (use in combination with SetQAFilterBits(...))
 for(Int_t fb=0;fb<gFilterBits;fb++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4
   {
    fQAKinematicsFilterBits[fb][rs][kv] = NULL;
   }   
  } 
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForQAHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForControlEventHistograms()
{
 // Initialize all arrays for control event histograms. Cuts hardwired here are default event cuts.

 // a) Multiplicity.
 // b) Centrality;
 // c) Vertex;
 // d) Remaining event histograms.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Multiplicity:
 fMultiplicityBins[0] = 3000.;
 fMultiplicityBins[1] = 0.;
 fMultiplicityBins[2] = 3000.;
 fMultiplicityCuts[0] = 0.;
 fMultiplicityCuts[1] = 3000.; 
 fSelectedTracksCuts[0] = -1;
 fSelectedTracksCuts[1] = 1e6; // TBI 20210515 is this safe?

 // b) Centrality:
 fCentralityBins[0] = 100.;
 fCentralityBins[1] = 0.;
 fCentralityBins[2] = 100.;
 for(Int_t ba=0;ba<2;ba++) // before/after cuts
 { 
  fCentralityHist[ba] = NULL;
 } // for(Int_t ba=0;ba<2;ba++) // before/after cuts
 fCentralityCuts[0] = 0.;
 fCentralityCuts[1] = 100.;
 
 // c) Vertex:
 fVertexBins[0] = 20000.;
 fVertexBins[1] = -50.;
 fVertexBins[2] = 50.;
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    fVertexHist[ba][rs][xyz] = NULL; 
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t xyz=0;xyz<3;xyz++)
 fVertexCuts[X][0] = -10.;
 fVertexCuts[X][1] = 10.;
 fVertexCuts[Y][0] = -10.;
 fVertexCuts[Y][1] = 10.;
 fVertexCuts[Z][0] = -10.;
 fVertexCuts[Z][1] = 10.;

 fNContributorsBins[0] = 3000.;
 fNContributorsBins[1] = 0.;
 fNContributorsBins[2] = 3000.;
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   fNContributorsHist[ba][rs] = NULL; 
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t xyz=0;xyz<3;xyz++)
 fNContributorsCuts[0] = 2; // min
 fNContributorsCuts[1] = 1e6; // max, this one is typically left open

 // d) Remaining event histograms:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t oeh=0;oeh<gEventHistograms;oeh++)
  {
   fEventHistograms[ba][oeh] = NULL;
  }
 } // for(Int_t ba=0;ba<2;ba++)

 // defaults bins:
 fEventBins[MagneticField][0] = 4;
 fEventBins[MagneticField][1] = -2.;
 fEventBins[MagneticField][2] = 2.;
 fEventBins[PrimaryVertex][0] = 1;
 fEventBins[PrimaryVertex][1] = 0.;
 fEventBins[PrimaryVertex][2] = 1.;

 // default cuts:
 fEventCuts[MagneticField][0] = -1.;
 fEventCuts[MagneticField][1] = 1.;
 fEventCuts[PrimaryVertex][0] = 0.;
 fEventCuts[PrimaryVertex][1] = 1.;

} // void AliAnalysisTaskMuPa::InitializeArraysForControlEventHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForControlParticleHistograms()
{
 // Initialize all arrays for control particle histograms. Cuts hardwired here are default particle cuts. 

 // a) Kinematics; 
 // b) DCA;
 // c) Remaining particle histograms; 
 // d) The default particle cuts. 

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Kinematics:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4
   {
    fKinematicsHist[ba][rs][kv] = NULL; 
   } // for(Int_t kv=0;kv<gKinematicVariables;kv++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++)

 // default phi binning:
 fKinematicsBins[PHI][0] = 360;
 fKinematicsBins[PHI][1] = 0.;
 fKinematicsBins[PHI][2] = TMath::TwoPi();

 // default pt binning:
 fKinematicsBins[PT][0] = 1000;
 fKinematicsBins[PT][1] = 0.;
 fKinematicsBins[PT][2] = 100.;

 // default eta binning:
 fKinematicsBins[ETA][0] = 200;
 fKinematicsBins[ETA][1] = -2.;
 fKinematicsBins[ETA][2] = 2.;

 // default e binning:
 fKinematicsBins[E][0] = 1000;
 fKinematicsBins[E][1] = 0.;
 fKinematicsBins[E][2] = 100.;

 // default charge binning:
 fKinematicsBins[CHARGE][0] = 5;
 fKinematicsBins[CHARGE][1] = -2.;
 fKinematicsBins[CHARGE][2] = 3.;

 // b) DCA:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t xyTz=0;xyTz<2;xyTz++)
   {
    fDCAHist[ba][rs][xyTz] = NULL; 
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++)

 // default DCA xy binning:
 fDCABins[0][0] = 1000;
 fDCABins[0][1] = -20.;
 fDCABins[0][2] = 20.;

 // default DCA z binning:
 fDCABins[1][0] = 1000;
 fDCABins[1][1] = -20.;
 fDCABins[1][2] = 20.;

 // c) Remaining particle histograms:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t t=0;t<gParticleHistograms;t++) // type, see enum 'eParticle'
   {
    fParticleHist[ba][rs][t] = NULL; 
   } // for(Int_t t=0;t<gParticleHistograms;t++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++)

 // default binning, see enum 'eParticle':
 fParticleBins[TPCNcls][0] = 200;
 fParticleBins[TPCNcls][1] = 0.;
 fParticleBins[TPCNcls][2] = 200.;

 fParticleBins[TPCnclsS][0] = 200;
 fParticleBins[TPCnclsS][1] = 0.;
 fParticleBins[TPCnclsS][2] = 200.;

 fParticleBins[TPCnclsFractionShared][0] = 200;
 fParticleBins[TPCnclsFractionShared][1] = 0.;
 fParticleBins[TPCnclsFractionShared][2] = 200.;

 fParticleBins[TPCNCrossedRows][0] = 200;
 fParticleBins[TPCNCrossedRows][1] = 0.;
 fParticleBins[TPCNCrossedRows][2] = 200.;

 fParticleBins[TPCChi2perNDF][0] = 100;
 fParticleBins[TPCChi2perNDF][1] = 0.;
 fParticleBins[TPCChi2perNDF][2] = 10.;

 fParticleBins[TPCFoundFraction][0] = 100;
 fParticleBins[TPCFoundFraction][1] = 0.;
 fParticleBins[TPCFoundFraction][2] = 2.;

 fParticleBins[Chi2TPCConstrainedVsGlobal][0] = 400;
 fParticleBins[Chi2TPCConstrainedVsGlobal][1] = -2.;
 fParticleBins[Chi2TPCConstrainedVsGlobal][2] = 2.;

 fParticleBins[ITSNcls][0] = 100;
 fParticleBins[ITSNcls][1] = 0.;
 fParticleBins[ITSNcls][2] = 10.;

 fParticleBins[ITSChi2perNDF][0] = 100;
 fParticleBins[ITSChi2perNDF][1] = 0.;
 fParticleBins[ITSChi2perNDF][2] = 1.;

 // d) The default particle cuts: 
 //  These cuts are used by default, if you want other cuts, indicate that via dedicated setters.
 //  d1) Kinematics:
 //    default phi cuts:
 fKinematicsCuts[PHI][0] = fKinematicsBins[PHI][1];
 fKinematicsCuts[PHI][1] = fKinematicsBins[PHI][2];

 //    default pt cuts:
 fKinematicsCuts[PT][0] = fKinematicsBins[PT][1];
 fKinematicsCuts[PT][1] = fKinematicsBins[PT][2];

 //    default eta cuts:
 fKinematicsCuts[ETA][0] = fKinematicsBins[ETA][1];
 fKinematicsCuts[ETA][1] = fKinematicsBins[ETA][2];

 //    default e cuts:
 fKinematicsCuts[E][0] = fKinematicsBins[E][1];
 fKinematicsCuts[E][1] = fKinematicsBins[E][2];

 //    default charge cuts:
 fKinematicsCuts[CHARGE][0] = fKinematicsBins[CHARGE][1];
 fKinematicsCuts[CHARGE][1] = fKinematicsBins[CHARGE][2];

 //  d2) DCA:
 //    default DCA xy cuts:
 fDCACuts[0][0] = fDCABins[0][1];
 fDCACuts[0][1] = fDCABins[0][2];

 //    default DCA z cuts:
 fDCACuts[1][0] = fDCABins[1][1];
 fDCACuts[1][1] = fDCABins[1][2];

 //  d3) Remaining:
 for(Int_t t=0;t<gParticleHistograms;t++)
 {
  fParticleCuts[t][0] = fParticleBins[t][1];
  fParticleCuts[t][1] = fParticleBins[t][2];
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForControlParticleHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForCommonLabels()
{
 // Initialize all arrays holding common labels.

 fBeforeAfterColor[0] = kRed; 
 fBeforeAfterColor[1] = kBlue; 

} // void AliAnalysisTaskMuPa::InitializeArraysForCommonLabels()

//=======================================================================================================================

void AliAnalysisTaskMuPa::InsanityChecks()
{
 // Check before bookings if all the values user has provided via setters make sense.

 return; // 20210517 disabled temporarily 

 // a) Multiplicity;
 // b) Centrality.

  cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) TBI Multiplicity:
 // ...

 // b) Centrality:
 if(((Int_t)fCentralityBins[0])<=0){cout<<__LINE__<<endl;exit(1);}
 if(fCentralityBins[1]<0. || fCentralityBins[1] > 100.){cout<<__LINE__<<endl;exit(1);}
 if(fCentralityBins[2]<=0. || fCentralityBins[2] > 100.){cout<<__LINE__<<endl;exit(1);}
 if(fCentralityBins[2]<=fCentralityBins[1]){cout<<__LINE__<<endl;exit(1);}

 // To do:
 // 1/ cout<<"FATAL: Only \"x\", \"y\" or \"z\" are supported for the 1st argument in this setter."<<endl; exit(0); TBI 20210512 move this to the insanity checks
 // 2/ the same thing for SetKinematicsBins(..) and SetKinematicsCuts(..)

} // void AliAnalysisTaskMuPa::InsanityChecks()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookBaseProfile()
{
 // Book base profile which keeps flags relevant for the whole analysis.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 fBasePro = new TProfile("fBasePro","flags for the whole analysis",4,0.,4.);
 fBasePro->SetStats(kFALSE);
 fBasePro->GetXaxis()->SetBinLabel(1,Form("fPeriod = %s",fPeriod.Data()));
 fBasePro->GetXaxis()->SetBinLabel(2,Form("fAODNumber = %s",fAODNumber.Data())); // TBI 20210513 
 //fBasePro->GetXaxis()->SetBinLabel(3,"fFillQAhistograms"); fBasePro->Fill(2.5,fFillQAHistograms);
 //fBasePro->GetXaxis()->SetBinLabel(4,"fTerminateAfterQA"); fBasePro->Fill(3.5,fTerminateAfterQA);
 fBaseList->Add(fBasePro);

} // void AliAnalysisTaskMuPa::BookBaseProfile()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fBaseList.

 // a) Book and nest lists for QA histograms;
 // b) Book and nest lists for control event histograms;
 // c) Book and nest lists for control particle histograms;
 // d) Book and nest lists for Q-vectorsl;
 // e) Book and nest all lists for correlations;
 // f) Book and nest all lists for nested loops;

 // *) Book and nest lists for final results.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 if(!fBaseList){cout<<__LINE__<<endl;exit(1);}

 // a) Book and nest lists for QA histograms:
 fQAList = new TList();
 fQAList->SetName("QA");
 fQAList->SetOwner(kTRUE);
 fBaseList->Add(fQAList);

 // b) Book and nest lists for control histograms:
 fControlEventHistogramsList = new TList();
 fControlEventHistogramsList->SetName("ControlEventHistograms");
 fControlEventHistogramsList->SetOwner(kTRUE);
 fBaseList->Add(fControlEventHistogramsList);

 // c) Book and nest lists for control particle histograms;
 fControlParticleHistogramsList = new TList();
 fControlParticleHistogramsList->SetName("ControlParticleHistograms");
 fControlParticleHistogramsList->SetOwner(kTRUE);
 fBaseList->Add(fControlParticleHistogramsList);

 // d) Book and nest lists for Q-vectors:
 fQvectorList = new TList();
 fQvectorList->SetName("Q-vectors");
 fQvectorList->SetOwner(kTRUE);
 fBaseList->Add(fQvectorList);

 // e) Book and nest all lists for correlations:
 fCorrelationsList = new TList();
 fCorrelationsList->SetName("Correlations");
 fCorrelationsList->SetOwner(kTRUE);
 fBaseList->Add(fCorrelationsList);

 // f) Book and nest all lists for nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("NestedLoops");
 fNestedLoopsList->SetOwner(kTRUE);
 fBaseList->Add(fNestedLoopsList);

 // ...

 return;

 // * ) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fBaseList->Add(fFinalResultsList);

} // void AliAnalysisTaskMuPa::BookAndNestAllLists()

//=======================================================================================================================


void AliAnalysisTaskMuPa::BookQAHistograms()
{
 // Book all QA histograms.

 // a) Book the profile holding flags;
 // b) Common local style and labels;
 // c) Centrality;
 // d) Particles;

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Book the profile holding flags:
 fQAPro = new TProfile("fQAPro","flags for QA histograms",1,0.,1.);
 fQAPro->SetStats(kFALSE);
 fQAList->Add(fQAPro);

 if(!fFillQAHistograms){return;}

 // b) Common local style and labels:
 Int_t bLineColor = kGreen+2;
 Int_t bFillColor = kGreen-10;
 TString sce[gCentralityEstimators] = {"V0M", "CL0", "CL1", "TRK", "V0A", "V0B", "TKL"}; // keep this in sync with enum eCentralityEstimator in .h
 TString sba[2] = {"before particle cuts","after particle cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString skv[gKinematicVariables] = {"#varphi","p_{T}","#eta","energy","charge"};

 for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
 {
  fQACentralityHist[ce1] = new TH1D(Form("fQACentralityHist[%d]",ce1),Form("%s",sce[ce1].Data()),(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]);
  //fQACentralityHist[ce1]->SetStats(kFALSE);
  fQACentralityHist[ce1]->SetLineColor(bLineColor);
  fQACentralityHist[ce1]->SetFillColor(bFillColor);
  fQACentralityHist[ce1]->GetXaxis()->SetTitle("centrality");
  fQAList->Add(fQACentralityHist[ce1]);

  // 2D correlation plot, only upper diagonal is booked:
  for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
  {
   fQACentralityCorrHist[ce1][ce2] = new TH2D(Form("fQACentralityCorrHist[%d][%d]",ce1,ce2),Form("%s vs. %s",sce[ce1].Data(),sce[ce2].Data()),
                                   (Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2],(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]);
   //fQACentralityCorrHist[ce1][ce2]->SetLineColor(bLineColor);
   //fQACentralityCorrHist[ce1][ce2]->SetFillColor(bFillColor);
   fQACentralityCorrHist[ce1][ce2]->SetOption("col");
   fQACentralityCorrHist[ce1][ce2]->GetXaxis()->SetTitle(Form("centrality %s",sce[ce1].Data()));
   fQACentralityCorrHist[ce1][ce2]->GetYaxis()->SetTitle(Form("centrality %s",sce[ce2].Data()));
   fQAList->Add(fQACentralityCorrHist[ce1][ce2]);
  } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
 } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)

 // Particles:
 fQAFilterBitScan = new TH1I("fQAFilterBitScan","fQAFilterBitScan",gFilterBits,0,gFilterBits);
 fQAFilterBitScan->SetStats(kFALSE);
 fQAFilterBitScan->SetLineColor(bLineColor);
 fQAFilterBitScan->SetFillColor(bFillColor);
 fQAFilterBitScan->SetXTitle("FilterBit");
 fQAFilterBitScan->SetYTitle("# particles");
 for(Int_t fb=0;fb<gFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
 {
  fQAFilterBitScan->GetXaxis()->SetBinLabel(fb+1,Form("%d",1<<fb));
 }
 fQAList->Add(fQAFilterBitScan);
 
 Int_t nIDsMax = 20000;
 fQAIDvsFilterBit = new TH2I("fQAIDvsFilterBit","fQAIDvsFilterBit",gFilterBits,0,gFilterBits,2*nIDsMax,-nIDsMax,nIDsMax);
 fQAIDvsFilterBit->SetOption("col");
 fQAIDvsFilterBit->SetXTitle("FilterBit");
 fQAIDvsFilterBit->SetYTitle("atrack->GetID()");
 for(Int_t fb=0;fb<gFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
 {
  fQAIDvsFilterBit->GetXaxis()->SetBinLabel(fb+1,Form("%d",1<<fb));
 }
 fQAList->Add(fQAIDvsFilterBit);

 // b) Kinematics for specified filter bits (use in combination with SetQAFilterBits(...))
 for(Int_t fb=0;fb<fQAFilterBits->GetSize();fb++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4
   {
    fQAKinematicsFilterBits[fb][rs][kv] = new TH1D(Form("fQAKinematicsFilterBits[%d][%d][%d]",fb,rs,kv),Form("Filter bit: %d, %s, %s",(Int_t)fQAFilterBits->GetAt(fb),srs[rs].Data(),skv[kv].Data()),(Int_t)fKinematicsBins[kv][0],fKinematicsBins[kv][1],fKinematicsBins[kv][2]);
    fQAKinematicsFilterBits[fb][rs][kv]->SetXTitle(skv[kv].Data());
    fQAKinematicsFilterBits[fb][rs][kv]->SetLineColor(kBlack);
    fQAKinematicsFilterBits[fb][rs][kv]->SetFillColor(kGray);
    fQAKinematicsFilterBits[fb][rs][kv]->SetMinimum(0.);
    fQAList->Add(fQAKinematicsFilterBits[fb][rs][kv]);
   }   
  } 
 }

} // void AliAnalysisTaskMuPa::BookQAHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookControlEventHistograms()
{
 // Book all control event histograms.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Multiplicity;
 // d) Centrality;
 // e) Vertex;
 // f) Remaining event histograms.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Book the profile holding flags:
 fControlEventHistogramsPro = new TProfile("fControlEventHistogramsPro","flags for control event histograms",1,0.,1.);
 fControlEventHistogramsPro->SetStats(kFALSE);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(1,Form("fCentralityEstimator = %s",fCentralityEstimator.Data()));
 fControlEventHistogramsList->Add(fControlEventHistogramsPro);

 // b) Common local labels:
 //    Remark: Keep them in sync with enums in .h
 TString sxyz[3] = {"x","y","z"};
 TString sba[2] = {"before event cuts","after event cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString stype[gEventHistograms] = {"MagneticField","PrimaryVertex"};

 // c) Multiplicity:
 fMultiplicityHist = new TH1D("fMultiplicityHist","multiplicity = sum of particle weights of tracks in Q-vector",(Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2]);
 //fMultiplicityHist->SetStats(kFALSE);
 fMultiplicityHist->SetLineColor(COLOR);
 fMultiplicityHist->SetFillColor(FILLCOLOR);
 fMultiplicityHist->GetXaxis()->SetTitle("Multiplicity");
 fControlEventHistogramsList->Add(fMultiplicityHist);

 fSelectedTracksHist = new TH1I("fSelectedTracksHist","integer counter of selected tracks added in Q-vector",(Int_t)fMultiplicityBins[0],(Int_t)fMultiplicityBins[1],(Int_t)fMultiplicityBins[2]);
 //fSelectedTracksHist->SetStats(kFALSE);
 fSelectedTracksHist->SetLineColor(COLOR);
 fSelectedTracksHist->SetFillColor(FILLCOLOR);
 fSelectedTracksHist->GetXaxis()->SetTitle("integer counter of selected tracks ");
 fControlEventHistogramsList->Add(fSelectedTracksHist);

 // d) Centrality:
 for(Int_t ba=0;ba<2;ba++) // before/after cuts
 { 
  fCentralityHist[ba] = new TH1D(Form("fCentralityHist[%d]",ba),Form("%s, %s",sba[ba].Data(),fCentralityEstimator.Data()),(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]);
  //fCentralityHist[ba]->SetStats(kFALSE);
  fCentralityHist[ba]->SetLineColor(fBeforeAfterColor[ba]);
  fCentralityHist[ba]->SetFillColor(fBeforeAfterColor[ba]-10);
  fCentralityHist[ba]->GetXaxis()->SetTitle("centrality");
  fControlEventHistogramsList->Add(fCentralityHist[ba]);
 } // for(Int_t ba=0;ba<2;ba++) // before/after cuts

 // e) Vertex:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    fVertexHist[ba][rs][xyz] = new TH1D(Form("fVertexHist[%d][%d][%d]",ba,rs,xyz),Form("%s, %s",sba[ba].Data(),srs[rs].Data()),(Int_t)fVertexBins[0],fVertexBins[1],fVertexBins[2]); 
    //fVertexHist[ba][rs][xyz]->SetStats(kFALSE);
    fVertexHist[ba][rs][xyz]->GetXaxis()->SetTitle(Form("V_{%s}",sxyz[xyz].Data()));
    fVertexHist[ba][rs][xyz]->SetLineColor(fBeforeAfterColor[ba]);
    fVertexHist[ba][rs][xyz]->SetFillColor(fBeforeAfterColor[ba]-10);
    fControlEventHistogramsList->Add(fVertexHist[ba][rs][xyz]); 
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t xyz=0;xyz<3;xyz++)

 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   fNContributorsHist[ba][rs] = new TH1I(Form("fNContributorsHist[%d][%d]",ba,rs),Form("%s, %s",sba[ba].Data(),srs[rs].Data()),(Int_t)fNContributorsBins[0],fNContributorsBins[1],fNContributorsBins[2]); 
   //fNContributorsHist[ba][rs]->SetStats(kFALSE);
   fNContributorsHist[ba][rs]->GetXaxis()->SetTitle("avtx->GetNContributors()");
   fNContributorsHist[ba][rs]->SetLineColor(fBeforeAfterColor[ba]);
   fNContributorsHist[ba][rs]->SetFillColor(fBeforeAfterColor[ba]-10);
   fControlEventHistogramsList->Add(fNContributorsHist[ba][rs]); 
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t xyz=0;xyz<3;xyz++)

 // f) Remaining event histograms:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t t=0;t<gEventHistograms;t++) // type, see enum 'eEvent'
  {
   fEventHistograms[ba][t] = new TH1D(Form("fEventHistograms[%d][%d]",ba,t),Form("%s, %s",sba[ba].Data(),stype[t].Data()),(Int_t)fEventBins[t][0],fEventBins[t][1],fEventBins[t][2]); 
   fEventHistograms[ba][t]->SetLineColor(fBeforeAfterColor[ba]);  
   fEventHistograms[ba][t]->SetFillColor(fBeforeAfterColor[ba]-10);
   fControlEventHistogramsList->Add(fEventHistograms[ba][t]); 
  } // for(Int_t t=0;t<gEventHistograms;t++) // type, see enum 'eEvent'
 } // for(Int_t ba=0;ba<2;ba++)

} // void AliAnalysisTaskMuPa::BookControlEventHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookControlParticleHistograms()
{
 // Book all control particle histograms.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Kinematics;
 // d) DCA;
 // e) Remaining.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Book the profile holding flags:
 fControlParticleHistogramsPro = new TProfile("fControlParticleHistogramsPro","flags for control particle histograms",1,0.,1.);
 fControlParticleHistogramsPro->SetStats(kFALSE);
 fControlParticleHistogramsPro->GetXaxis()->SetBinLabel(1,"fFilterBit"); fControlParticleHistogramsPro->Fill(0.5,fFilterBit);
 fControlParticleHistogramsList->Add(fControlParticleHistogramsPro);

 // b) Common local labels:
 //    Remark: Keep them in sync with enums in .h
 TString sxyTz[2] = {"xy","z"};
 TString sba[2] = {"before particle cuts","after particle cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString skv[gKinematicVariables] = {"#varphi","p_{T}","#eta","energy","charge"};
 TString stype[gParticleHistograms] = {"TPCNcls","TPCnclsS","TPCnclsFractionShared","TPCNCrossedRows","TPCChi2perNDF","TPCFoundFraction","Chi2TPCConstrainedVsGlobal","ITSNcls","ITSChi2perNDF"};

 // c) Kinematics:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4 TBI 20210512 this is not enforced to be in sync with the definition of enums
   {
    fKinematicsHist[ba][rs][kv] = new TH1D(Form("fKinematicsHist[%d][%d][%d]",ba,rs,kv),Form("%s, %s",sba[ba].Data(),srs[rs].Data()),(Int_t)fKinematicsBins[kv][0],fKinematicsBins[kv][1],fKinematicsBins[kv][2]); 
    //fKinematicsHist[ba][rs][kv]->SetStats(kFALSE);
    fKinematicsHist[ba][rs][kv]->GetXaxis()->SetTitle(skv[kv].Data());
    fKinematicsHist[ba][rs][kv]->SetMinimum(0.);
    fKinematicsHist[ba][rs][kv]->SetLineColor(fBeforeAfterColor[ba]);
    fKinematicsHist[ba][rs][kv]->SetFillColor(fBeforeAfterColor[ba]-10);
    fControlParticleHistogramsList->Add(fKinematicsHist[ba][rs][kv]); 
   } // for(Int_t kv=0;kv<gKinematicVariables;kv++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++) 

 // d) DCA:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t xyTz=0;xyTz<2;xyTz++) 
   {
    fDCAHist[ba][rs][xyTz] = new TH1D(Form("fDCAHist[%d][%d][%d]",ba,rs,xyTz),Form("%s, %s",sba[ba].Data(),srs[rs].Data()),(Int_t)fDCABins[xyTz][0],fDCABins[xyTz][1],fDCABins[xyTz][2]); 
    //fDCAHist[ba][rs][xyTz]->SetStats(kFALSE);
    fDCAHist[ba][rs][xyTz]->GetXaxis()->SetTitle(sxyTz[xyTz].Data());
    fDCAHist[ba][rs][xyTz]->SetMinimum(0.);
    fDCAHist[ba][rs][xyTz]->SetLineColor(fBeforeAfterColor[ba]);
    fDCAHist[ba][rs][xyTz]->SetFillColor(fBeforeAfterColor[ba]-10);
    fControlParticleHistogramsList->Add(fDCAHist[ba][rs][xyTz]); 
   } // for(Int_t xyTz=0;xyTz<gKinematicVariables;xyTz++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++) 

 // e) Remaining:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   for(Int_t t=0;t<gParticleHistograms;t++) // type, see enum 'eParticle'
   {
    fParticleHist[ba][rs][t] = new TH1D(Form("fParticleHist[%d][%d][%d]",ba,rs,t),Form("%s, %s, %s",stype[t].Data(),sba[ba].Data(),srs[rs].Data()),(Int_t)fParticleBins[t][0],fParticleBins[t][1],fParticleBins[t][2]);  
    fParticleHist[ba][rs][t]->GetXaxis()->SetTitle(stype[t].Data());
    fParticleHist[ba][rs][t]->SetMinimum(0.);  
    fParticleHist[ba][rs][t]->SetLineColor(fBeforeAfterColor[ba]);
    fParticleHist[ba][rs][t]->SetFillColor(fBeforeAfterColor[ba]-10);
    fControlParticleHistogramsList->Add(fParticleHist[ba][rs][t]); 
   } // for(Int_t t=0;t<gParticleHistograms;t++)
  } // for(Int_t rs=0;rs<2;rs++)
 } // for(Int_t ba=0;ba<2;ba++)

} // void AliAnalysisTaskMuPa::BookControlParticleHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookCorrelationsHistograms()
{
 // Book all correlations histograms.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Histograms.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Book the profile holding flags:
 fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro","flags for correlations",1,0.,1.);
 fCorrelationsFlagsPro->SetStats(kFALSE);
 fCorrelationsFlagsPro->GetXaxis()->SetLabelSize(0.05); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateCorrelations");  
 fCorrelationsFlagsPro->Fill(0.5,fCalculateCorrelations);
 fCorrelationsList->Add(fCorrelationsFlagsPro);

 if(!fCalculateCorrelations){return;}

 // b) Common local labels:
 TString oVariable[4] = {"#varphi_{1}-#varphi_{2}","#varphBookQvectorHistogramsi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
                         "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
                         "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-#varphi_{7}-#varphi_{8}"};
 Int_t vvvariableNBins[3] = {1,(Int_t)fMultiplicityBins[0],(Int_t)fCentralityBins[0]};
 Double_t vvvariableMinMax[3][2] = { {0.,1.}, // integrated 
                                    {fMultiplicityBins[1],fMultiplicityBins[2]}, // multiplicity
                                    {fCentralityBins[1],fCentralityBins[2]} // centrality
                                   };
 TString vvVariable[3] = {"integrated","multiplicity","centrality"};

 for(Int_t k=0;k<4;k++) // order [2p=0,4p=1,6p=2,8p=3]
 {  
  for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
  {
   for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
   {
    fCorrelationsPro[k][n][v] = new TProfile(Form("fCorrelationsPro[%d][%d][%d]",k,n,v),Form("#LT#LTcos[%s(%s)]#GT#GT",1==n+1?"":Form("%d",n+1),oVariable[k].Data()),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
    fCorrelationsPro[k][n][v]->SetStats(kFALSE);
    fCorrelationsPro[k][n][v]->Sumw2();
    fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(vvVariable[v].Data());
    //fCorrelationsPro[k][n][v]->SetFillColor(colorsW[v]-10);
    //fCorrelationsPro[k][n][v]->SetLineColor(colorsW[v]);
    fCorrelationsList->Add(fCorrelationsPro[k][n][v]);
   } // for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
  } // for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
 } // for(Int_t n=0;n<6;n++) // harmonics [n=1,n=2,...,n=6]

} // void AliAnalysisTaskMuPa::BookCorrelationsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookNestedLoopsHistograms()
{
 // Book all nested loops histograms.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Histograms.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Book the profile holding flags:
 fNestedLoopsFlagsPro = new TProfile("fNestedLoopsFlagsPro","flags for nested loops",1,0.,1.);
 fNestedLoopsFlagsPro->SetStats(kFALSE);
 fNestedLoopsFlagsPro->GetXaxis()->SetLabelSize(0.05); 
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateNestedLoops");  
 fNestedLoopsFlagsPro->Fill(0.5,fCalculateNestedLoops);
 fNestedLoopsList->Add(fNestedLoopsFlagsPro);

 if(!fCalculateNestedLoops){return;}

 const Int_t maxSize = 20000;
 ftaNestedLoops[0] = new TArrayD(maxSize); // ebe container for azimuthal angles 
 ftaNestedLoops[1] = new TArrayD(maxSize); // ebe container for particle weights (product of all)  

 // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms())
 TString oVariable[4] = {"#varphi_{1}-#varphi_{2}","#varphBookQvectorHistogramsi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
                         "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
                         "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-#varphi_{7}-#varphi_{8}"};
 Int_t vvvariableNBins[3] = {1,(Int_t)fMultiplicityBins[0],(Int_t)fCentralityBins[0]};
 Double_t vvvariableMinMax[3][2] = { {0.,1.}, // integrated 
                                    {fMultiplicityBins[1],fMultiplicityBins[2]}, // multiplicity
                                    {fCentralityBins[1],fCentralityBins[2]} // centrality
                                   };
 TString vvVariable[3] = {"integrated","multiplicity","centrality"};

 for(Int_t k=0;k<4;k++) // order [2p=0,4p=1,6p=2,8p=3]
 {  
  for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
  {
   for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
   {
    fNestedLoopsPro[k][n][v] = new TProfile(Form("fNestedLoopsPro[%d][%d][%d]",k,n,v),Form("#LT#LTcos[%s(%s)]#GT#GT",1==n+1?"":Form("%d",n+1),oVariable[k].Data()),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
    fNestedLoopsPro[k][n][v]->SetStats(kFALSE);
    fNestedLoopsPro[k][n][v]->Sumw2();
    fNestedLoopsPro[k][n][v]->GetXaxis()->SetTitle(vvVariable[v].Data());
    //fNestedLoopsPro[k][n][v]->SetFillColor(colorsW[v]-10);
    //fNestedLoopsPro[k][n][v]->SetLineColor(colorsW[v]);
    fNestedLoopsList->Add(fNestedLoopsPro[k][n][v]);
   } // for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
  } // for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
 } // for(Int_t n=0;n<6;n++) // harmonics [n=1,n=2,...,n=6]

} // void AliAnalysisTaskMuPa::BookNestedLoopsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl; 

} // void AliAnalysisTaskMuPa::BookFinalResultsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookQvectorHistograms()
{
 // TBI 20210514 documentation

 // ... Flags for Q-vectors: TBI 20210513 temporarily here
 fQvectorFlagsPro = new TProfile("fQvectorFlagsPro","flags for Q-vector objects",3,0.,3.);
 fQvectorFlagsPro->SetStats(kFALSE);
 fQvectorFlagsPro->SetLineColor(COLOR);
 fQvectorFlagsPro->SetFillColor(FILLCOLOR);
 fQvectorFlagsPro->GetXaxis()->SetLabelSize(0.05); 
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateQvector");  
 fQvectorFlagsPro->Fill(0.5,fCalculateQvector);
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(2,"fMaxHarmonic");  
 fQvectorFlagsPro->Fill(1.5,fMaxHarmonic);
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(3,"fMaxCorrelator");  
 fQvectorFlagsPro->Fill(2.5,fMaxCorrelator);
 fQvectorList->Add(fQvectorFlagsPro);

} // void AliAnalysisTaskMuPa::BookQvectorHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::GetPointers(TList *baseList)
{
 // Get all pointers. This method is  essential only for merging and boostrap.

 // a) Check the pointer for base list fBaseList;
 // b) Get pointer for profile holding internal flags and set again all flags;
 // c) Get pointers for control event histograms;
 // d) Get pointers control particle histograms;

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Check the pointer for base list fBaseList:
 fBaseList = baseList; 
 if(!fBaseList){cout<<__LINE__<<endl;exit(1);}

 // b) Get pointer for profile holding internal flags and set again all flags:
 fBasePro = dynamic_cast<TProfile*>(fBaseList->FindObject("fBasePro"));
 if(!fBasePro){cout<<__LINE__<<endl;exit(1);}
 // old example: fCalculateCorrelations = (Bool_t)fBasePro->GetBinContent(1);

 // c) Get pointers for control event histograms:
 this->GetPointersForControlEventHistograms();
 
 // c) Get pointers for control particle histograms:
 this->GetPointersForControlParticleHistograms();

} // void AliAnalysisTaskMuPa::GetPointers(TList *baseList)

//=======================================================================================================================

void AliAnalysisTaskMuPa::GetPointersForControlEventHistograms(void)
{
 // Get pointers for control event histograms:

 // a) Get pointer for fControlEventHistogramsList:
 // b) Get pointer for fControlEventHistogramsPro;
 // c) Re-initiate all relevant flags; 
 // d) Get all specific pointers.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Get pointer for fControlEventHistogramsList:
 fControlEventHistogramsList = dynamic_cast<TList*>(fBaseList->FindObject("ControlEventHistograms"));
 if(!fControlEventHistogramsList){cout<<__LINE__<<endl;exit(1);}

 // b) Get pointer for fControlEventHistogramsPro:
 fControlEventHistogramsPro = dynamic_cast<TProfile*>(fControlEventHistogramsList->FindObject("fControlEventHistogramsPro"));
 if(!fControlEventHistogramsPro){cout<<__LINE__<<endl;exit(1);}

 // c) Re-initiate all relevant flags:
 // old example: fCalculateCorrelations = (Bool_t)fControlEventHistogramsPro->GetBinContent(1);

 // d) Get all specific pointers:
 // ... TBI multiplicity

 // Centrality:
 for(Int_t ba=0;ba<2;ba++) // before/after cuts
 { 
  fCentralityHist[ba] = dynamic_cast<TH1D*>(fControlEventHistogramsList->FindObject(Form("fCentralityHist[%d]",ba)));
  if(!fCentralityHist[ba])
  {
   cout<<Form("fCentralityHist[%d]",ba)<<endl;
   cout<<__LINE__<<endl;exit(1); 
  }  
 } // for(Int_t ba=0;ba<2;ba++) // before/after cuts
  
} // void AliAnalysisTaskMuPa::GetPointersForControlEventHistograms(void)

//=======================================================================================================================

void AliAnalysisTaskMuPa::GetPointersForControlParticleHistograms(void)
{
 // Get pointers for control event histograms:

 // a) Get pointer for fControlParticleHistogramsList:
 // b) Get pointer for fCorrelationsFlagsPro;
 // c) Re-initiate all relevant flags; 
 // d) Get all specific pointers.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Get pointer for fControlParticleHistogramsList:
 fControlParticleHistogramsList = dynamic_cast<TList*>(fBaseList->FindObject("ControlParticleHistograms"));
 if(!fControlParticleHistogramsList){cout<<__LINE__<<endl;exit(1);}

 // b) Get pointer for fCorrelationsFlagsPro:
 fControlParticleHistogramsPro = dynamic_cast<TProfile*>(fControlParticleHistogramsList->FindObject("fControlParticleHistogramsPro"));
 if(!fControlParticleHistogramsPro){cout<<__LINE__<<endl;exit(1);}

 // c) Re-initiate all relevant flags:
 // old example: fCalculateCorrelations = (Bool_t)fControlParticleHistogramsPro->GetBinContent(1);
 
 // d) Get all specific pointers:
 // ...

} // void AliAnalysisTaskMuPa::GetPointersForControlParticleHistograms(void)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FilterEvent(AliVEvent *ave)
{
 // Filter from this event the ebe values for local data members.

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Centrality;

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  // Centrality:
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  fCentrality = ams->GetMultiplicityPercentile(fCentralityEstimator.Data());

  // ... 
 } // if(aAOD)

 return;

} // void AliAnalysisTaskMuPa::FilterEvent(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)
{
 // Fill control histograms before event cuts (ba = 0), or after (ba = 1). For reconstructed data (rs = 0), or simulated (rs = 1)
 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Centrality;
 // c) Vertex;
 // d) Remaining event distributions.

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  // b) Centrality:
  if(fCentralityHist[ba]){fCentralityHist[ba]->Fill(fCentrality);}

  // c) Vertex: TH1D *fVertexHist[2][2][3]; //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) return; 
  if(fVertexHist[ba][rs][X]){fVertexHist[ba][rs][X]->Fill(avtx->GetX());}
  if(fVertexHist[ba][rs][Y]){fVertexHist[ba][rs][Y]->Fill(avtx->GetY());}
  if(fVertexHist[ba][rs][Z]){fVertexHist[ba][rs][Z]->Fill(avtx->GetZ());}
  if(fNContributorsHist[ba][rs]){fNContributorsHist[ba][rs]->Fill(avtx->GetNContributors());}

  // d) Remaining event distributions:
  if(fEventHistograms[ba][MagneticField]){fEventHistograms[ba][MagneticField]->Fill(aAOD->GetMagneticField());}
  if(fEventHistograms[ba][PrimaryVertex]){fEventHistograms[ba][PrimaryVertex]->Fill(0.44);} // here we only count # of events with valid pointer aAOD->GetPrimaryVertex()
 } // if(aAOD)

} // AliAnalysisTaskMuPa::FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FillControlParticleHistograms(AliAODTrack *aTrack, const Int_t ba, const Int_t rs)
{
 // Fill control histograms before particles cuts (ba = 0), or after (ba = 1). For reconstructed data (rs = 0), or simulated (rs = 1)
 // More detailed treatment of AOD track cuts can be found in Tasks/AliFlowTrackCuts.cxx
 
 //cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // Kinematics:
 if(fKinematicsHist[ba][rs][PHI]){fKinematicsHist[ba][rs][PHI]->Fill(aTrack->Phi());}
 if(fKinematicsHist[ba][rs][PT]){fKinematicsHist[ba][rs][PT]->Fill(aTrack->Pt());}
 if(fKinematicsHist[ba][rs][ETA]){fKinematicsHist[ba][rs][ETA]->Fill(aTrack->Eta());}
 if(fKinematicsHist[ba][rs][E]){fKinematicsHist[ba][rs][E]->Fill(aTrack->E());}
 if(fKinematicsHist[ba][rs][CHARGE]){fKinematicsHist[ba][rs][CHARGE]->Fill(aTrack->Charge());}

 // DCA:
 if(fDCAHist[ba][rs][0]){fDCAHist[ba][rs][0]->Fill(aTrack->DCA());} // "xy"
 if(fDCAHist[ba][rs][1]){fDCAHist[ba][rs][1]->Fill(aTrack->ZAtDCA());} // "z"

 // Remaining:
 if(fParticleHist[ba][rs][TPCNcls]){fParticleHist[ba][rs][TPCNcls]->Fill(aTrack->GetTPCNcls());}
 if(fParticleHist[ba][rs][TPCnclsS]){fParticleHist[ba][rs][TPCnclsS]->Fill(aTrack->GetTPCnclsS());}
 if(fParticleHist[ba][rs][TPCnclsFractionShared]){if(TMath::Abs(aTrack->GetTPCnclsS())>0.){fParticleHist[ba][rs][TPCnclsFractionShared]->Fill(aTrack->GetTPCNcls()/aTrack->GetTPCnclsS());}}
 if(fParticleHist[ba][rs][TPCNCrossedRows]){fParticleHist[ba][rs][TPCNCrossedRows]->Fill(aTrack->GetTPCNCrossedRows());}
 if(fParticleHist[ba][rs][TPCChi2perNDF]){fParticleHist[ba][rs][TPCChi2perNDF]->Fill(aTrack->Chi2perNDF());}
 if(fParticleHist[ba][rs][TPCFoundFraction]){fParticleHist[ba][rs][TPCFoundFraction]->Fill( aTrack->GetTPCFoundFraction());}
 if(fParticleHist[ba][rs][Chi2TPCConstrainedVsGlobal]){fParticleHist[ba][rs][Chi2TPCConstrainedVsGlobal]->Fill(aTrack->GetChi2TPCConstrainedVsGlobal());}
 if(fParticleHist[ba][rs][ITSNcls]){fParticleHist[ba][rs][ITSNcls]->Fill(aTrack->GetITSNcls());}
 if(fParticleHist[ba][rs][ITSChi2perNDF]){if(TMath::Abs(aTrack->GetITSNcls())>0.){fParticleHist[ba][rs][ITSChi2perNDF]->Fill(aTrack->GetITSchi2()/aTrack->GetITSNcls());}}

} // AliAnalysisTaskMuPa::FillControlParticleHistograms(AliAODTrack *aTrack, const Int_t ba, const Int_t rs)

//=======================================================================================================================

Bool_t AliAnalysisTaskMuPa::SurvivesEventCuts(AliVEvent *ave)
{
 // Check if the current event survives event cuts. TBI 20210512 works only for 'reco' at the moment

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) AOD.

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  // Centrality:
  if(fCentrality < fCentralityCuts[0]) return kFALSE;
  if(fCentrality > fCentralityCuts[1]) return kFALSE;
 
  // Vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) return kFALSE; 
  if((Int_t)avtx->GetNContributors()<fNContributorsCuts[0]) return kFALSE;
  if((Int_t)avtx->GetNContributors()>fNContributorsCuts[1]) return kFALSE;
  if(avtx->GetX() < fVertexCuts[X][0]) return kFALSE;
  if(avtx->GetX() > fVertexCuts[X][1]) return kFALSE; 
  if(avtx->GetY() < fVertexCuts[Y][0]) return kFALSE;
  if(avtx->GetY() > fVertexCuts[Y][1]) return kFALSE;
  if(avtx->GetZ() < fVertexCuts[Z][0]) return kFALSE;
  if(avtx->GetZ() > fVertexCuts[Z][1]) return kFALSE;

  // Remaining event cuts:
  // TBI 20210518 magnetic field value

 } // if(aAOD)

 return kTRUE;

} // AliAnalysisTaskMuPa::SurvivesEventCuts(AliVEvent *ave)

//=======================================================================================================================

Bool_t AliAnalysisTaskMuPa::SurvivesParticleCuts(AliAODTrack *aTrack)
{
 // Check if the current partice survives the specific track cuts (e.g. applied only on TPC-only tracks).

 //cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // TBI 20210517 most likely fUseKinematicsCuts is obsolete, since now I use the default cuts by default.

 //if(fUseKinematicsCuts[PHI])
 {
  if(aTrack->Phi() < fKinematicsCuts[PHI][0]) return kFALSE;
  if(aTrack->Phi() >= fKinematicsCuts[PHI][1]) return kFALSE;
 }

 //if(fUseKinematicsCuts[PT])
 {
  if(aTrack->Pt() < fKinematicsCuts[PT][0]) return kFALSE;
  if(aTrack->Pt() >= fKinematicsCuts[PT][1]) return kFALSE;
 }

 //if(fUseKinematicsCuts[ETA])
 {
  if(aTrack->Eta() < fKinematicsCuts[ETA][0]) return kFALSE;
  if(aTrack->Eta() >= fKinematicsCuts[ETA][1]) return kFALSE;
 }

 //if(fUseKinematicsCuts[E])
 {
  if(aTrack->E() < fKinematicsCuts[E][0]) return kFALSE;
  if(aTrack->E() >= fKinematicsCuts[E][1]) return kFALSE;
 }

 //if(fUseKinematicsCuts[CHARGE])
 {
  if(aTrack->Charge() < fKinematicsCuts[CHARGE][0]) return kFALSE;
  if(aTrack->Charge() >= fKinematicsCuts[CHARGE][1]) return kFALSE;
 }

 // DCA xy:
 if(aTrack->DCA() < fDCACuts[0][0]) return kFALSE;
 if(aTrack->DCA() >= fDCACuts[0][1]) return kFALSE;
 
 // DCA z:
 if(aTrack->ZAtDCA() < fDCACuts[1][0]) return kFALSE;
 if(aTrack->ZAtDCA() >= fDCACuts[1][1]) return kFALSE;

 // Remaining:
 if(aTrack->GetTPCNcls() < fParticleCuts[TPCNcls][0]) return kFALSE;
 if(aTrack->GetTPCNcls() >= fParticleCuts[TPCNcls][1]) return kFALSE;
 if(aTrack->GetTPCnclsS() < fParticleCuts[TPCnclsS][0]) return kFALSE;
 if(aTrack->GetTPCnclsS() >= fParticleCuts[TPCnclsS][1]) return kFALSE;
 if(TMath::Abs(aTrack->GetTPCnclsS())>0.)
 {
  if(aTrack->GetTPCNcls()/aTrack->GetTPCnclsS() < fParticleCuts[TPCnclsFractionShared][0]) return kFALSE;
  if(aTrack->GetTPCNcls()/aTrack->GetTPCnclsS() >= fParticleCuts[TPCnclsFractionShared][1]) return kFALSE;
 }
 if(aTrack->GetTPCNCrossedRows() < fParticleCuts[TPCNCrossedRows][0]) return kFALSE;
 if(aTrack->GetTPCNCrossedRows() >= fParticleCuts[TPCNCrossedRows][1]) return kFALSE;
 if(aTrack->Chi2perNDF() < fParticleCuts[TPCChi2perNDF][0]) return kFALSE;
 if(aTrack->Chi2perNDF() >= fParticleCuts[TPCChi2perNDF][1]) return kFALSE;
 if(aTrack->GetTPCFoundFraction() < fParticleCuts[TPCFoundFraction][0]) return kFALSE;
 if(aTrack->GetTPCFoundFraction() >= fParticleCuts[TPCFoundFraction][1]) return kFALSE;
 if(aTrack->GetChi2TPCConstrainedVsGlobal() < fParticleCuts[Chi2TPCConstrainedVsGlobal][0]) return kFALSE;
 if(aTrack->GetChi2TPCConstrainedVsGlobal() >= fParticleCuts[Chi2TPCConstrainedVsGlobal][1]) return kFALSE;
 if(aTrack->GetITSNcls() < fParticleCuts[ITSNcls][0]) return kFALSE;
 if(aTrack->GetITSNcls() >= fParticleCuts[ITSNcls][1]) return kFALSE;
 if(TMath::Abs(aTrack->GetITSNcls())>0.)
 {
  if(aTrack->GetITSchi2()/aTrack->GetITSNcls() < fParticleCuts[ITSChi2perNDF][0]) return kFALSE;
  if(aTrack->GetITSchi2()/aTrack->GetITSNcls() >= fParticleCuts[ITSChi2perNDF][1]) return kFALSE;
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMuPa::SurvivesParticleCuts(AliAODTrack *aTrack)

//=======================================================================================================================

void AliAnalysisTaskMuPa::QA(AliVEvent *ave)
{
 // Fill all additional QA stuff here, that is already not filled in Control Event Histograms or Control Particle Histograms. 

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Centrality;

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  // a) Centrality: 
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  TString sce[gCentralityEstimators] = {"V0M", "CL0", "CL1", "TRK", "V0A", "V0B", "TKL"}; // keep this in sync with enum eCentralityEstimator in .h
  
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   if(fQACentralityHist[ce1]){fQACentralityHist[ce1]->Fill(ams->GetMultiplicityPercentile(sce[ce1].Data()));}
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    if(fQACentralityCorrHist[ce1][ce2]){fQACentralityCorrHist[ce1][ce2]->Fill(ams->GetMultiplicityPercentile(sce[ce1].Data()),ams->GetMultiplicityPercentile(sce[ce2].Data()));}
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
   
  // *) Start loop over AOD tracks:
  Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  {
   AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to "a track" (i.e. any track)
   if(!aTrack){continue;}
   
   // Filter bit distribution:
   for(Int_t fb=0;fb<gFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
   {
    if(aTrack->TestFilterBit(1<<fb))
    {
     fQAFilterBitScan->Fill(fb);
     fQAIDvsFilterBit->Fill(fb,aTrack->GetID());
    }
   }

   // Filter bit kinematics:
   for(Int_t fb=0;fb<fQAFilterBits->GetSize();fb++)
   {
    if(aTrack->TestFilterBit((Int_t)fQAFilterBits->GetAt(fb)))
    {
     fQAKinematicsFilterBits[fb][0][PHI]->Fill(aTrack->Phi());
     fQAKinematicsFilterBits[fb][0][PT]->Fill(aTrack->Pt());
     fQAKinematicsFilterBits[fb][0][ETA]->Fill(aTrack->Eta());
     fQAKinematicsFilterBits[fb][0][E]->Fill(aTrack->E());
     fQAKinematicsFilterBits[fb][0][CHARGE]->Fill(aTrack->Charge());
    }
   }
   
  } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks

 } // if(aAOD)

} // void AliAnalysisTaskMuPa::QA(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMuPa::GlobalTracksAOD(AliAODEvent *aAOD)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAOD.

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: There is a performance penalty when fGlobalTracksAOD[1] and fGlobalTracksAOD[2] needed for mixed events are calculated.
 //           Yes, I can get them directly from fGlobalTracksAOD[0], without calling this method for them again. TBI today

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 if(0 != fGlobalTracksAOD->GetSize()){fGlobalTracksAOD->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAOD->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMuPa::GlobalTracksAOD(AliAODEvent *aAOD)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Q(Int_t n, Int_t wp)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][wp];} 
 return TComplex::Conjugate(fQvector[-n][wp]);
 
} // TComplex AliAnalysisTaskMuPa::Q(Int_t n, Int_t wp)


//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::One(Int_t n1)
{
 // Generic expression <exp[i(n1*phi1)]>. TBI comment

 TComplex one = Q(n1,1);

 return one;

} // TComplex AliAnalysisTaskMuPa::One(Int_t n1)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Two(Int_t n1, Int_t n2)
{
 // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.

 TComplex two = Q(n1,1)*Q(n2,1)-Q(n1+n2,2);

 return two;

} // TComplex AliAnalysisTaskMuPa::Two(Int_t n1, Int_t n2)

//=======================================================================================================================

void AliAnalysisTaskMuPa::CalculateCorrelations() 
{
 // Calculate analytically multiparticle correlations from Q-vectors. 
 // By default, only correlations for which all harmonics are the same are evaluated. 

 for(Int_t h=1;h<=fMaxHarmonic;h++) // harmonic
 {
  // 2p:
  //cout<<"   => CalculateCorrelations(void), 2p .... "<<endl;
  if(fMultiplicity<2){return;}
  TComplex twoC = Two(h,-h)/Two(0,0).Re(); // cos
  //TComplex twoS = Two(h,-h)/Two(0,0).Im(); // sin
  Double_t wTwo = Two(0,0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
  // integrated:
  if(fCorrelationsPro[0][h-1][0]){fCorrelationsPro[0][h-1][0]->Fill(0.5,twoC,wTwo);}
  // vs. multiplicity:
  if(fCorrelationsPro[0][h-1][1]){fCorrelationsPro[0][h-1][1]->Fill(fMultiplicity+0.5,twoC,wTwo);}
  // vs. centrality:
  if(fCorrelationsPro[0][h-1][2]){fCorrelationsPro[0][h-1][2]->Fill(fCentrality,twoC,wTwo);}
 } 

} // void AliAnalysisTaskMuPa::CalculateCorrelations() 

//=======================================================================================================================

void AliAnalysisTaskMuPa::CalculateNestedLoops() 
{
 // Calculate correlations with nested loops. 

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 cout<<"fSelectedTracks = "<<fSelectedTracks<<endl;

 // 2p:
 if(fSelectedTracks<2){return;}
 cout<<"      CalculateNestedLoops(void), 2-p correlations .... "<<endl;
 for(int i1=0; i1<fSelectedTracks; i1++)
 {
  Double_t dPhi1 = ftaNestedLoops[0]->GetAt(i1);
  Double_t dW1 = ftaNestedLoops[1]->GetAt(i1);
  for(int i2=0; i2<fSelectedTracks; i2++)
  {
   if(i2==i1){continue;}
   Double_t dPhi2 = ftaNestedLoops[0]->GetAt(i2);
   Double_t dW2 = ftaNestedLoops[1]->GetAt(i2);
   for(int h=0; h<6; h++)
   {
    // fill cos, 2p, integreated: 
    fNestedLoopsPro[0][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);
    // fill cos, 2p, vs. M: 
    fNestedLoopsPro[0][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);
    // fill cos, 2p, vs. centrality: 
    fNestedLoopsPro[0][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);
   } // for(int h=1; h<=6; h++)
  } // for(int i2=0; i2<nTracks; ++i2)
 } // for(int i1=0; i1<nTracks; ++i1)

 return; // TBI 20210515 remove when you add support for below

 // 4p:
 if(fSelectedTracks<4){return;}
 cout<<"      CalculateNestedLoops(void), 4-p correlations .... "<<endl;
 for(int i1=0; i1<fSelectedTracks; i1++)
 {
  Double_t dPhi1 = ftaNestedLoops[0]->GetAt(i1);
  Double_t dW1 = ftaNestedLoops[1]->GetAt(i1);
  for(int i2=0; i2<fSelectedTracks; i2++)
  {
   if(i2==i1){continue;}
   Double_t dPhi2 = ftaNestedLoops[0]->GetAt(i2);
   Double_t dW2 = ftaNestedLoops[1]->GetAt(i2);
   for(int i3=0; i3<fSelectedTracks; i3++)
   {
    if(i3==i1||i3==i2){continue;}
    Double_t dPhi3 = ftaNestedLoops[0]->GetAt(i3);
    Double_t dW3 = ftaNestedLoops[1]->GetAt(i3);
    for(int i4=0; i4<fSelectedTracks; i4++)
    {
     if(i4==i1||i4==i2||i4==i3){continue;}
     Double_t dPhi4 = ftaNestedLoops[0]->GetAt(i4);
     Double_t dW4 = ftaNestedLoops[1]->GetAt(i4);
     for(int h=0; h<6; h++)
     {
      // fill cos, 4p, integreated: 
      fNestedLoopsPro[1][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);
      // fill cos, 4p, all harmonics, vs. M: 
      fNestedLoopsPro[1][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);
      // fill cos, 4p, all harmonics, vs. centrality: 
      fNestedLoopsPro[1][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);
     } // for(int h=0; h<6; h++)
    } // for(int i4=0; i4<fSelectedTracks; i4++)   
   } // for(int i3=0; i3<fSelectedTracks; i3++)
  } // for(int i2=0; i2<nTracks; ++i2)
 } // for(int i1=0; i1<nTracks; ++i1)

 // 6p:
 if(fSelectedTracks<6){return;}
 cout<<"      CalculateNestedLoops(void), 6-p correlations .... "<<endl;
 for(int i1=0; i1<fSelectedTracks; i1++)
 {
  Double_t dPhi1 = ftaNestedLoops[0]->GetAt(i1);
  Double_t dW1 = ftaNestedLoops[1]->GetAt(i1);
  for(int i2=0; i2<fSelectedTracks; i2++)
  {
   if(i2==i1){continue;}
   Double_t dPhi2 = ftaNestedLoops[0]->GetAt(i2);
   Double_t dW2 = ftaNestedLoops[1]->GetAt(i2);
   for(int i3=0; i3<fSelectedTracks; i3++)
   {
    if(i3==i1||i3==i2){continue;}
    Double_t dPhi3 = ftaNestedLoops[0]->GetAt(i3);
    Double_t dW3 = ftaNestedLoops[1]->GetAt(i3);
    for(int i4=0; i4<fSelectedTracks; i4++)
    {
     if(i4==i1||i4==i2||i4==i3){continue;}
     Double_t dPhi4 = ftaNestedLoops[0]->GetAt(i4);
     Double_t dW4 = ftaNestedLoops[1]->GetAt(i4);
     for(int i5=0; i5<fSelectedTracks; i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
      Double_t dPhi5 = ftaNestedLoops[0]->GetAt(i5);
      Double_t dW5 = ftaNestedLoops[1]->GetAt(i5);
      for(int i6=0; i6<fSelectedTracks; i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
       Double_t dPhi6 = ftaNestedLoops[0]->GetAt(i6);
       Double_t dW6 = ftaNestedLoops[1]->GetAt(i6);
       for(int h=0; h<6; h++)
       {
        // fill cos, 6p, integreated: 
        fNestedLoopsPro[2][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);
        // fill cos, 6p, all harmonics, vs. M: 
        fNestedLoopsPro[2][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);
        // fill cos, 6p, all harmonics, vs. M: 
        fNestedLoopsPro[2][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);
       } // for(int h=0; h<6; h++)
      } // if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
     } // if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
    } // for(int i4=0; i4<fSelectedTracks; i4++)   
   } // for(int i3=0; i3<fSelectedTracks; i3++)
  } // for(int i2=0; i2<nTracks; ++i2)
 } // for(int i1=0; i1<nTracks; ++i1)

 // 8p:
 if(fSelectedTracks<8){return;}
 cout<<"      CalculateNestedLoops(void), 8-p correlations .... "<<endl;
 for(int i1=0; i1<fSelectedTracks; i1++)
 {
  Double_t dPhi1 = ftaNestedLoops[0]->GetAt(i1);
  Double_t dW1 = ftaNestedLoops[1]->GetAt(i1);
  for(int i2=0; i2<fSelectedTracks; i2++)
  {
   if(i2==i1){continue;}
   Double_t dPhi2 = ftaNestedLoops[0]->GetAt(i2);
   Double_t dW2 = ftaNestedLoops[1]->GetAt(i2);
   for(int i3=0; i3<fSelectedTracks; i3++)
   {
    if(i3==i1||i3==i2){continue;}
    Double_t dPhi3 = ftaNestedLoops[0]->GetAt(i3);
    Double_t dW3 = ftaNestedLoops[1]->GetAt(i3);
    for(int i4=0; i4<fSelectedTracks; i4++)
    {
     if(i4==i1||i4==i2||i4==i3){continue;}
     Double_t dPhi4 = ftaNestedLoops[0]->GetAt(i4);
     Double_t dW4 = ftaNestedLoops[1]->GetAt(i4);
     for(int i5=0; i5<fSelectedTracks; i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
      Double_t dPhi5 = ftaNestedLoops[0]->GetAt(i5);
      Double_t dW5 = ftaNestedLoops[1]->GetAt(i5);
      for(int i6=0; i6<fSelectedTracks; i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
       Double_t dPhi6 = ftaNestedLoops[0]->GetAt(i6);
       Double_t dW6 = ftaNestedLoops[1]->GetAt(i6);
       for(int i7=0; i7<fSelectedTracks; i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
        Double_t dPhi7 = ftaNestedLoops[0]->GetAt(i7);
        Double_t dW7 = ftaNestedLoops[1]->GetAt(i7);
        for(int i8=0; i8<fSelectedTracks; i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7){continue;}
         Double_t dPhi8 = ftaNestedLoops[0]->GetAt(i8);
         Double_t dW8 = ftaNestedLoops[1]->GetAt(i8);
         for(int h=0; h<6; h++)
         {
          // fill cos, 8p, integreated: 
          fNestedLoopsPro[3][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);
          // fill cos, 8p, all harmonics, vs. M: 
          fNestedLoopsPro[3][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);
          // fill cos, 8p, all harmonics, vs. M: 
          fNestedLoopsPro[3][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);
         } // for(int h=0; h<6; h++)
        } // for(int i8=0; i8<fSelectedTracks; i8++)
       } // for(int i7=0; i7<fSelectedTracks; i7++)
      } // if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
     } // if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
    } // for(int i4=0; i4<fSelectedTracks; i4++)   
   } // for(int i3=0; i3<fSelectedTracks; i3++)
  } // for(int i2=0; i2<nTracks; ++i2)
 } // for(int i1=0; i1<nTracks; ++i1)

} // void AliAnalysisTaskMuPa::CalculateNestedLoops(void) 

//=======================================================================================================================

void AliAnalysisTaskMuPa::ComparisonNestedLoopsVsCorrelations()
{
 // Make a ratio fNestedLoopsPro[....]/fCorrelationsPro[....]. If results are the same, these ratios must be 1.

 // a) Integrated comparison;
 // b) Comparison vs. multiplicity;
 // c) Comparison vs. centrality;

 cout<<"\n\033[1;32m"<<__PRETTY_FUNCTION__<<"\033[0m\n"<<endl;

 Int_t nBinsQV = -44;
 Int_t nBinsNL = -44;
 Double_t valueQV = 0.;
 Double_t valueNL = 0.;

 // a) Integrated comparison:
 nBinsQV = fCorrelationsPro[0][0][0]->GetNbinsX();
 nBinsNL = fNestedLoopsPro[0][0][0]->GetNbinsX();
 if(nBinsQV != nBinsNL){cout<<__LINE__<<endl; exit(1);}
 cout<<endl;
 cout<<"   [0] : integrated"<<endl;
 for(Int_t o=0;o<4;o++)
 {
  cout<<Form("   ==== %d-particle correlations ====",2*(o+1))<<endl;
  for(Int_t h=0;h<6;h++)
  {
   for(Int_t b=1;b<=nBinsQV;b++)
   {
    if(fCorrelationsPro[o][h][0]){valueQV = fCorrelationsPro[o][h][0]->GetBinContent(b);}
    if(fNestedLoopsPro[o][h][0]){valueNL = fNestedLoopsPro[o][h][0]->GetBinContent(b);}
    if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
    {
     cout<<Form("   h=%d, Q-vectors:    ",h+1)<<valueQV<<endl; 
     cout<<Form("   h=%d, Nested loops: ",h+1)<<valueNL<<endl; 
     if(TMath::Abs(valueQV-valueNL)>1.e-5)
     {          
      cout<<Form("[%d][%d][%d]",o,h,0)<<endl; cout<<__LINE__<<endl; exit(1);
     }
    } // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
   } // for(Int_t b=1;b<=nBinsQV;b++) 
  } // for(Int_t h=0;h<6;h++)
  cout<<endl;
 } // for(Int_t o=0;o<4;o++) 

 cout<<endl;

 // b) Comparison vs. multiplicity:
 nBinsQV = fCorrelationsPro[0][0][1]->GetNbinsX();
 nBinsNL = fNestedLoopsPro[0][0][1]->GetNbinsX();
 if(nBinsQV != nBinsNL){cout<<__LINE__<<endl; exit(1);}
 cout<<endl;
 cout<<"   [1] : vs. multiplicity"<<endl;
 for(Int_t o=0;o<4;o++)
 {
  cout<<Form("   ==== %d-particle correlations ====",2*(o+1))<<endl;
  for(Int_t h=0;h<6;h++)
  {
   for(Int_t b=1;b<=nBinsQV;b++)
   {
    if(fCorrelationsPro[o][h][1]){valueQV = fCorrelationsPro[o][h][1]->GetBinContent(b);}
    if(fNestedLoopsPro[o][h][1]){valueNL = fNestedLoopsPro[o][h][1]->GetBinContent(b);}
    if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
    {
     cout<<Form("   h=%d, b=%d, Q-vectors:    ",h+1,b)<<valueQV<<endl; 
     cout<<Form("   h=%d, b=%d, Nested loops: ",h+1,b)<<valueNL<<endl; 
     if(TMath::Abs(valueQV-valueNL)>1.e-5)
     {          
      cout<<Form("[%d][%d][%d]",o,h,1)<<endl; cout<<__LINE__<<endl; exit(1);
     }
    } // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
   } // for(Int_t b=1;b<=nBinsQV;b++) 
  } // for(Int_t h=0;h<6;h++)
  cout<<endl;
 } // for(Int_t o=0;o<4;o++) 

 cout<<endl;

 // c) Comparison vs. centrality:
 nBinsQV = fCorrelationsPro[0][0][2]->GetNbinsX();
 nBinsNL = fNestedLoopsPro[0][0][2]->GetNbinsX();
 if(nBinsQV != nBinsNL){cout<<__LINE__<<endl; exit(1);}
 cout<<endl;
 cout<<"   [2] : vs. centrality"<<endl;
 for(Int_t o=0;o<4;o++)
 {
  cout<<Form("   ==== %d-particle correlations ====",2*(o+1))<<endl;
  for(Int_t h=0;h<6;h++)
  {
   for(Int_t b=1;b<=nBinsQV;b++)
   {
    if(fCorrelationsPro[o][h][2]){valueQV = fCorrelationsPro[o][h][2]->GetBinContent(b);}
    if(fNestedLoopsPro[o][h][2]){valueNL = fNestedLoopsPro[o][h][2]->GetBinContent(b);}
    if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
    {
     cout<<Form("   h=%d, b=%d, Q-vectors:    ",h+1,b)<<valueQV<<endl; 
     cout<<Form("   h=%d, b=%d, Nested loops: ",h+1,b)<<valueNL<<endl; 
     if(TMath::Abs(valueQV-valueNL)>1.e-5)
     {          
      cout<<Form("[%d][%d][%d]",o,h,2)<<endl; cout<<__LINE__<<endl; exit(1);
     }
    } // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
   } // for(Int_t b=1;b<=nBinsQV;b++) 
  } // for(Int_t h=0;h<6;h++)
  cout<<endl;
 } // for(Int_t o=0;o<4;o++) 


} // void AliAnalysisTaskMuPa::ComparisonNestedLoopsVsCorrelations(void)




