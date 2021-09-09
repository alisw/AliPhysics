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
#include <AliMultSelection.h>
#include <TFile.h>
#include <unistd.h>
#include <TRandom3.h>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

ClassImp(AliAnalysisTaskMuPa)

//================================================================================================================

AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(const char *name): 
 AliAnalysisTaskSE(name), 
 fBaseList(NULL),
 fBasePro(NULL),
 fRealData(kTRUE),
 fUseFisherYates(kFALSE),
 fRandomIndices(NULL),
 fFillQAHistograms(kFALSE),
 fFillQAHistogramsAll(kFALSE),
 fTerminateAfterQA(kFALSE),
 fVerbose(kFALSE),
 fEventCounter(0),
 fRandomSeed(0),
 fUseTrigger(kFALSE),

 // QA:
 fQAList(NULL),
 fQAPro(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 fQAFilterBits(NULL),
 fQAAnomalousEvents(NULL),
 fQACheckSelfCorrelations(kFALSE),
 fQAEventCutCounter(NULL),
 fQASequentialEventCutCounter(NULL),

 // Control event histograms:
 fControlEventHistogramsList(NULL),
 fControlEventHistogramsPro(NULL),
 fMultiplicity(0.),
 fMultiplicityHist(NULL),
 fSelectedTracks(0),
 fSelectedTracksHist(NULL),
 fUseSelectedTracksCuts(kFALSE),
 fSelContrTreshold(0.),
 fCentrality(0.),
 fUseCentralityCuts(kFALSE),
 fCentralityCorrelationCutVersion(0),
 fUseNContributorsCuts(kFALSE),
 fMinVertexDistance(1.e-6),
 fUseMinVertexDistanceCut(kFALSE),

 // Control particle histograms:
 fControlParticleHistogramsList(NULL),
 fControlParticleHistogramsPro(NULL),
 fSimReco(NULL),
 fUseFakeTracks(kFALSE),
 fGlobalTracksAOD(NULL),
 fFilterGlobalTracksAOD(kFALSE),
 fFilterBit(128), // defaulted to TPC-only
 fUseOnlyPrimaries(kTRUE), 
 fAtLeastOnePointInTheSPD(kFALSE),
 fIgnoreGlobalConstrained(kFALSE),

 // Q-vectors:
 fQvectorList(NULL),
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kTRUE),
 //fMaxHarmonic(6),
 fMaxHarmonic(gMaxHarmonic),
 //fMaxCorrelator(8),
 fMaxCorrelator(gMaxCorrelator),

 // Particle weights:
 fWeightsList(NULL),
 fWeightsFlagsPro(NULL),

 // Centrality weights:
 fCentralityWeightsList(NULL),
 fCentralityWeightsFlagsPro(NULL),

 // Correlations:
 fCorrelationsList(NULL),        
 fCorrelationsFlagsPro(NULL), 
 fCalculateCorrelations(kTRUE),

 // Nested loops:
 fNestedLoopsList(NULL),        
 fNestedLoopsFlagsPro(NULL), 
 fCalculateNestedLoops(kFALSE), 
 fCalculateCustomNestedLoop(kFALSE),

 // Toy NUA:
 fToyNUAList(NULL), 
 fToyNUAFlagsPro(NULL),   

 // Internal validation:
 fInternalValidationList(NULL),
 fInternalValidationFlagsPro(NULL),
 fUseInternalValidation(kFALSE),
 fRescaleWithTheoreticalInput(kFALSE),
 fnEventsInternalValidation(1e4),
 fInternalValidationHarmonics(NULL),

 // Test0:
 fTest0List(NULL), 
 fTest0FlagsPro(NULL),   
 fCalculateTest0(kFALSE),
 fFileWithLabels(NULL),

 // Final results:
 fFinalResultsList(NULL),

 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateFile(NULL),
 fMaxNumberOfEvents(-44),
 fBailOutFile(NULL),

 // *.) Debugging:
 fProcessOnlySpecifiedEvent(kFALSE),
 fRun(0),
 fBunchCross(0),
 fOrbit(0),
 fPeriod(0),
 fPrintEventInfo(kFALSE)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskMuPa::AliAnalysisTaskMuPa(const char *name)");

  // Base list:
  fBaseList = new TList();
  fBaseList->SetName("outputMuPaAnalysis");
  fBaseList->SetOwner(kTRUE);

  // Initialize all non-built in types:
  fTaskName = name;
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
 fRealData(kTRUE),
 fUseFisherYates(kFALSE),
 fRandomIndices(NULL),
 fFillQAHistograms(kFALSE),
 fFillQAHistogramsAll(kFALSE),
 fTerminateAfterQA(kFALSE),
 fVerbose(kFALSE),
 fEventCounter(0),
 fRandomSeed(0),
 fUseTrigger(kFALSE),

 // QA:
 fQAList(NULL),
 fQAPro(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 fQAFilterBits(NULL),
 fQAAnomalousEvents(NULL),
 fQACheckSelfCorrelations(kFALSE),
 fQAEventCutCounter(NULL),
 fQASequentialEventCutCounter(NULL),
 
 // Control event histograms:
 fControlEventHistogramsList(NULL),
 fControlEventHistogramsPro(NULL),
 fMultiplicity(0.),
 fMultiplicityHist(NULL),
 fSelectedTracks(0),
 fSelectedTracksHist(NULL),
 fUseSelectedTracksCuts(kFALSE),
 fSelContrTreshold(0.),
 fCentrality(0.), 
 fUseCentralityCuts(kFALSE),
 fCentralityCorrelationCutVersion(0),
 fUseNContributorsCuts(kFALSE),
 fMinVertexDistance(0.),
 fUseMinVertexDistanceCut(kFALSE),

 // Control particle histograms:
 fControlParticleHistogramsList(NULL),
 fControlParticleHistogramsPro(NULL),
 fSimReco(NULL),
 fUseFakeTracks(kFALSE),
 fGlobalTracksAOD(NULL),
 fFilterGlobalTracksAOD(kFALSE),
 fFilterBit(-44),
 fUseOnlyPrimaries(kFALSE), 
 fAtLeastOnePointInTheSPD(kFALSE),
 fIgnoreGlobalConstrained(kFALSE),

 // Q-vectors:
 fQvectorList(NULL),
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kTRUE),
 //fMaxHarmonic(6),
 fMaxHarmonic(gMaxHarmonic),
 //fMaxCorrelator(8),
 fMaxCorrelator(gMaxCorrelator),

 // Particle weights:
 fWeightsList(NULL),
 fWeightsFlagsPro(NULL),

 // Centrality weights:
 fCentralityWeightsList(NULL),
 fCentralityWeightsFlagsPro(NULL),

 // Correlations:
 fCorrelationsList(NULL),        
 fCorrelationsFlagsPro(NULL), 
 fCalculateCorrelations(kTRUE),

 // Nested loops:
 fNestedLoopsList(NULL),        
 fNestedLoopsFlagsPro(NULL), 
 fCalculateNestedLoops(kFALSE),
 fCalculateCustomNestedLoop(kFALSE),

 // Toy NUA:
 fToyNUAList(NULL), 
 fToyNUAFlagsPro(NULL),   

 // Internal validation:
 fInternalValidationList(NULL),
 fInternalValidationFlagsPro(NULL),
 fUseInternalValidation(kFALSE),
 fRescaleWithTheoreticalInput(kFALSE),
 fnEventsInternalValidation(0),
 fInternalValidationHarmonics(NULL),

 // Test0:
 fTest0List(NULL), 
 fTest0FlagsPro(NULL),   
 fCalculateTest0(kFALSE),
 fFileWithLabels(NULL),

 // Final results:
 fFinalResultsList(NULL),

 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateFile(NULL),
 fMaxNumberOfEvents(-44),
 fBailOutFile(NULL),

 // *.) Debugging:
 fProcessOnlySpecifiedEvent(kFALSE),
 fRun(0),
 fBunchCross(0),
 fOrbit(0),
 fPeriod(0),
 fPrintEventInfo(kFALSE)
{
  // Dummy constructor.
 
  // Important: arrays have to be initialized also in the dummy constructor. 

  // Initialize all non-built in types:
  fTaskName = "";
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

 if(fFilterGlobalTracksAOD)
 { 
  if(fGlobalTracksAOD) delete fGlobalTracksAOD;
 } 
 
 if(fSimReco) delete fSimReco;

 if(fUseFisherYates) delete fRandomIndices;

} // AliAnalysisTaskMuPa::~AliAnalysisTaskMuPa()

//================================================================================================================

void AliAnalysisTaskMuPa::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Check before bookings if all the values user has provided via setters make sense;
 // b) Trick to avoid name clashes, part 1;
 // c) Book random generator;
 // d) Book base profile;
 // e) Book and nest all lists;
 // f) Book all objects;
 // g) Book all look-up tables;
 // *) Trick to avoid name clashes, part 2.
 
 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Check before bookings if all the values user has provided via setters make sense:
 this->InsanityChecks();
 
 // b) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // c) Book random generator:
 delete gRandom;
 gRandom = new TRandom3(fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

 // d) Book base profiles:
 this->BookBaseProfile();

 // e) Book and nest all lists:
 this->BookAndNestAllLists();

 // f) Book all objects:
 this->BookQAHistograms();
 this->BookControlEventHistograms();
 this->BookControlParticleHistograms();
 this->BookQvectorHistograms();
 this->BookWeightsHistograms();
 this->BookCentralityWeightsHistograms();
 this->BookCorrelationsHistograms();
 this->BookNestedLoopsHistograms(); 
 this->BookToyNUAHistograms(); 
 this->BookInternalValidationHistograms(); 
 this->BookTest0Histograms(); 
 this->BookFinalResultsHistograms();

 // g) Book all look-up tables:
 fSimReco = new TExMap(); // TBI 20210802 this doesn't have to be booked by default, only when running over MC
 if(fFilterGlobalTracksAOD)
 {
  fGlobalTracksAOD = new TExMap();
 }

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fBaseList);

} // void AliAnalysisTaskMuPa::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMuPa::UserExec(Option_t *) 
{
 // Main loop (called for each event).
 // a) Get pointer to AOD event;
 // b) Fisher-Yates algorithm;
 // c) Filter out "normal global" tracks for default analysis and cut on their number (needed only for PID studies);
 // d) Filter from this event only what is needed, and store that info in local data members;
 // e) QA (if enabled);
 // f) Fill control event histograms before cuts;
 // g) Event cuts;
 // h) Fill control event histograms after cuts;
 // i) Look up table;
 // j) Start analysis over AODs;
 // k) Fill e-b-e quantities;
 // l) Calculate correlations;
 // m) Calculate nested loops;
 // n) Calculate Test0;
 // o) Reset event-by-event objects;
 // *) PostData.

 Green(__PRETTY_FUNCTION__); 
 fEventCounter++;

 if(fUseInternalValidation){this->InternalValidation(); return;} // only do internal validation for all implemented correlations against the theoretical values

 // a) Get pointer to AOD event:
 AliMCEvent *aMC = MCEvent();                                  // from TaskSE
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE

 if(!aAOD){return;} // This means that at the moment I cannot process only kine, i.e. I process only reco, or both kine + reco
 if(fPrintEventInfo){this->PrintEventInfo(aAOD);}
 if(fProcessOnlySpecifiedEvent){if(!this->SpecifiedEvent(aAOD)){return;}} 

 // b) Fisher-Yates algorithm:
 if(fUseFisherYates){this->RandomIndices(aAOD);}

 // c) Filter out "normal global" tracks for default analysis and cut on their number (needed only for PID studies):
 //    'TPC-only' tracks and 'global constrained to vertex' come with negative ID, and are therefore not stored in fGlobalTracksAOD
 if(fFilterGlobalTracksAOD)
 {
  this->GlobalTracksAOD(aAOD); 
  if(0 == fGlobalTracksAOD->GetSize()) return; // yes, go to next event TBI 20210513 re-think this line, perhaps add some further check
 }

 // d) Filter from this event only what is needed, and store that info in local data members:
 this->FilterEvent(aAOD);

 // e) QA (if enabled):
 if(fFillQAHistograms)
 {
  FillQAHistograms(aAOD,BEFORE,RECO); 
  if(SurvivesEventCuts(aAOD)){FillQAHistograms(aAOD,AFTER,RECO);}
  // Nothing is need at the moment for SIM, otherwise use in the same spirit:
  //  if(aMC){FillQAHistograms(aMC,BEFORE,SIM);} 
  //  if(aMC){FillQAHistograms(aMC,AFTER,SIM);}  
  // Just make it in sync. with call to SurvivesEventCuts(aAOD), so that ControlEvent histos below are filled consistently

  // Special QA cases when both AliMCEvent and AliAODEvent are needed:
  if(SurvivesEventCuts(aAOD))
  {
   if(aAOD && aMC){FillQAHistograms(aAOD,aMC);}
  }

  if(fTerminateAfterQA){return;}
 } // if(fFillQAHistograms)

 // f) Fill control event histograms before cuts:
 if(aMC){this->FillControlEventHistograms(aAOD,BEFORE,SIM);}
 this->FillControlEventHistograms(aAOD,BEFORE,RECO);
 
 // g) Event cuts:
 if(!SurvivesEventCuts(aAOD)){return;} // TBI 20210531 add possibility to run only on sim, when this needs to be generalized. For kine+reco, I do not need to cut only on kine

 // h) Fill control event histograms after cuts:
 if(aMC){this->FillControlEventHistograms(aAOD,AFTER,SIM);}
 this->FillControlEventHistograms(aAOD,AFTER,RECO);

 // i) Look up table:
 if(aMC){this->MakeLookUpTable(aAOD,aMC);}

 // j) Start analysis over AODs:
 Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
 Double_t dPt = 0., wPt = 1.; // transverse momentum and corresponding pT weight
 Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
 Double_t wToPowerP = 1.; // weight raised to power p
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 Int_t nSelectedTracksCounter = 0; // needed only to fill nested loops containers
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aodTrack = NULL;
  if(!fUseFisherYates)
  {
   aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to "a track" (i.e. any track)
  }
  else
  {
   aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack((Int_t)fRandomIndices->GetAt(iTrack)));
  }
  if(!aodTrack){continue;}
 
  // Particle histograms and track cuts:
  if(!aMC)
  { 
   // fill directly only reco:
   this->FillControlParticleHistograms(aodTrack,BEFORE,RECO);
   if(!SurvivesParticleCuts(aodTrack)){continue;}
   this->FillControlParticleHistograms(aodTrack,AFTER,RECO);
  }

  // Finally, calculate Q-vectors:
  dPhi = aodTrack->Phi();
  dPt = aodTrack->Pt(); dPt += 0.; // TBI 20210515 shut down the compiler warnings temporarily
  dEta = aodTrack->Eta(); dEta += 0.; // TBI 20210515 shut down the compiler warnings temporarily

  // Particle weights:
  if(fUseWeights[0])
  {
   wPhi = Weight(dPhi,"phi"); // corresponding phi weight
   if(!(wPhi > 0.))
   {
    cout<<"wPhi is not positive, skipping this particle for the time being..."<<endl;
    cout<<Form("iTrack= %d\ndPhi = %f\nwPhi = %f",iTrack,dPhi,wPhi)<<endl;
    sleep(2);
    continue;
   } 
  } // if(fUseWeights[0])
  if(fUseWeights[1])
  {
   wPt = Weight(dPt,"pt"); // corresponding pt weight
   if(!(wPt > 0.))
   {
    cout<<"wPt is not positive, skipping this particle for the time being..."<<endl;
    cout<<Form("iTrack= %d\ndPt = %f\nwPt = %f",iTrack,dPt,wPt)<<endl;
    sleep(2);
    continue;
   } 
  } // if(fUseWeights[1])
  if(fUseWeights[2])
  {
   wEta = Weight(dEta,"eta"); // corresponding eta weight
   if(!(wEta > 0.))
   {
    cout<<"wEta is not positive, skipping this particle for the time being..."<<endl;
    cout<<Form("iTrack= %d\ndEta = %f\nwEta = %f",iTrack,dEta,wEta)<<endl;
    sleep(2);
    continue;
   } 
  } // if(fUseWeights[2])

  for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  {
   for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
   {
    if(fUseWeights[0]||fUseWeights[1]||fUseWeights[2]){wToPowerP = pow(wPhi*wPt*wEta,wp);} 
    fQvector[h][wp] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));    
   } // for(Int_t wp=0;wp<fMaxCorrelator+1;wp++)
  } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)   

  // Nested loops containers:
  if(fCalculateNestedLoops||fCalculateCustomNestedLoop)
  {
   if(ftaNestedLoops[0]){ftaNestedLoops[0]->AddAt(dPhi,nSelectedTracksCounter);} 
   if(ftaNestedLoops[1]){ftaNestedLoops[1]->AddAt(wPhi*wPt*wEta,nSelectedTracksCounter);} 
   nSelectedTracksCounter++;
  }

  // Sum of particle weights:
  fMultiplicity += wPhi*wPt*wEta; // only if weights are unit, fMultiplicity = fSelectedTracks

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks

 // Efficiency corrections are obtained in the loop below.
 // If Monte Carlo information is available, fill control particle histograms by using the look up table, to get all efficiencies.
 if(aMC)
 {
  Int_t nLabels = aMC->GetNumberOfTracks(); // total number of Monte Carlo tracks
  for(Int_t iLabel=0;iLabel<nLabels;iLabel++) // starting a loop over all Monte Carlo tracks (i.e. labels in this context)
  {
   AliAODMCParticle *aodmcParticle = (AliAODMCParticle*)aMC->GetTrack(iLabel);
   if(!aodmcParticle){continue;}  
   // Get the corresponding reco particle:
   // fSimReco->Add(label,iTrack); // "key" = label, "value" = iTrack
   AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(fSimReco->GetValue(iLabel)));
   
   if(fUseFisherYates){cout<<__LINE__<<endl;exit(1);} // TBI 20210810 check and validate if also here Fisher-Yates needs to be applied

   // Track cuts are applied at MC particle.
   // a) Sim and reco distributions before cuts:
   this->FillControlParticleHistograms(aodmcParticle,BEFORE,SIM);
   if(aodTrack)
   {
    this->FillControlParticleHistograms(aodTrack,BEFORE,RECO); // this one is with positive label, therefore this is not a fake track
   }   
   else if(!aodTrack && fUseFakeTracks)
   {
    aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(fSimReco->GetValue(TMath::Abs(iLabel))));
    this->FillControlParticleHistograms(aodTrack,BEFORE,RECO); // this one is with negative label, therefore this is a fake track
   }   
   if(!aodTrack){continue;}

   // b) Cut:
   if(!SurvivesParticleCuts(aodmcParticle)){continue;} // for Monte Carlo particles, apply only the basic kinematic cuts + cut on charge

   // c) Sim and reco distributions after cuts:
   this->FillControlParticleHistograms(aodmcParticle,AFTER,SIM); // these are then generated particles in the desired phase-space window
   if(!SurvivesParticleCuts(aodTrack)){continue;} // if Monte Carlo particle after reconstruction went out of phase-space window, we do not reconstruct it
   this->FillControlParticleHistograms(aodTrack,AFTER,RECO);  

  } // for(Int_t iLabel=0;iLabel<nLabels;iLabel++) // starting a loop over all Monte Carlo tracks
 } // if(aMC)

 // k) Fill e-b-e quantities:
 if(fMultiplicityHist){fMultiplicityHist->Fill(fMultiplicity);}
 if(fSelectedTracksHist){fSelectedTracksHist->Fill(fSelectedTracks);}

 // l) Calculate correlations:
 if(fCalculateCorrelations){this->CalculateCorrelations();}

 // m) Calculate nested loops:
 if(fCalculateNestedLoops){this->CalculateNestedLoops();}

 // n) Calculate Test0:
 if(fCalculateTest0){this->CalculateTest0();}

 // o) Reset event-by-event objects:
 this->ResetEventByEventQuantities();

 // *) PostData:
 PostData(1,fBaseList);

 // *) Online monitoring:
 Green(Form("INFO: Total event counter: %d\n      Selected event counter: %d (%.2f%%)",fEventCounter,(Int_t)fSelectedTracksHist->GetEntries(),100.*(Int_t)fSelectedTracksHist->GetEntries()/fEventCounter)); 
 if(fOnlineMonitoring){this->OnlineMonitoring();}

} // void AliAnalysisTaskMuPa::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskMuPa::Terminate(Option_t *)
{
 // Accessing the merged output list.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) "online" mode:
 if(!fBaseList)
 {
  fBaseList = (TList*)GetOutputData(1);
  if(!fBaseList){cout<<__LINE__<<endl;exit(1);}
  this->GetPointers(fBaseList);
 }

 // b) "offline" mode:
 // 1/ get 'list' from external ROOT file 'mergedFile'
 // 2/ load flow library and initialize the task: 
 // 	AliAnalysisTaskMuPa *mupa = new AliAnalysisTaskMuPa("MuPa");
 //     mupa->GetPointers(list);
 //     mupa->Terminate(NULL);
 // 3/ close the curtains:
 //     mergedFile->Close();  
 
 // c) Do some calculation in terminate here:
 if(fCalculateNestedLoops){this->ComparisonNestedLoopsVsCorrelations();}

} // end of void AliAnalysisTaskMuPa::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskMuPa::ResetEventByEventQuantities()
{
 // Reset all global event-by-event quantities here:

 // a) Multiplicities;
 // b) Centrality;
 // c) Q-vector;
 // d) Reset ebe containers for nested loops;
 // e) Fisher-Yates algorithm.
 
 // a) Multiplicities:
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
 if(fCalculateNestedLoops||fCalculateCustomNestedLoop)
 {
  if(ftaNestedLoops[0]){ftaNestedLoops[0]->Reset();} 
  if(ftaNestedLoops[1]){ftaNestedLoops[1]->Reset();}  
 }

 // e) Fisher-Yates algorithm:
 if(fUseFisherYates)
 {
  delete fRandomIndices;
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
  if(!fSelectedTracksHist){cout<<__LINE__<<endl;exit(1);}  
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
  if(!fSelectedTracksHist){cout<<__LINE__<<endl;exit(1);}  
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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 fDataTakingPeriod = TString("not set"); // can be customized with e.g. task->SetDataTakingPeriod("LHC10h");
 fAODNumber = TString("not set"); // can be customized with e.g. task->SetAODNumber("AOD160"); 
 fRunNumber = TString("not set"); // can be customized with e.g. task->SetRunNumber("000123456"); 
 fCentralityEstimator = TString("V0M"); // by default, we use V0M as centrality estimator. Can be customized with task->SetCentralityEstimator("V0M") 
 fTrigger = TString("not set"); 

} // void AliAnalysisTaskMuPa::InitializeNonBuiltInTypes()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 this->InitializeArraysForQAHistograms();
 this->InitializeArraysForControlEventHistograms();
 this->InitializeArraysForControlParticleHistograms();
 this->InitializeArraysForQvectors();
 this->InitializeArraysForWeights();
 this->InitializeArraysForCentralityWeights();
 this->InitializeArraysForCorrelationsHistograms();
 this->InitializeArraysForNestedLoopsHistograms();
 this->InitializeArraysForToyNUA();
 this->InitializeArraysForInternalValidation();
 this->InitializeArraysForTest0();
 this->InitializeArraysForCommonLabels();

} // void AliAnalysisTaskMuPa::InitializeArrays()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForQvectors()
{
 // Initialize all arrays for Q-vectors.

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

void AliAnalysisTaskMuPa::InitializeArraysForWeights()
{
 // Initialize all arrays for particle weights.

 for(Int_t w=0;w<gWeights;w++) 
 {
  fUseWeights[w] = kFALSE;
  fWeightsHist[w] = NULL;
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForWeights()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForCentralityWeights()
{
 // Initialize all arrays for centrality weights.

 for(Int_t w=0;w<gCentralityEstimators;w++) 
 {
  fUseCentralityWeights[w] = kFALSE;
  fCentralityWeightsHist[w] = NULL;
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForCentralityWeights()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForNestedLoopsHistograms()
{
 // Initialize all arrays of histograms holding results for nested loops.

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

void AliAnalysisTaskMuPa::InitializeArraysForToyNUA()
{
 // Initialize all arrays for Toy NUA.

 for(Int_t kv=0;kv<gKinematicVariables;kv++) 
 { 
  fUseToyNUA[kv] = kFALSE;
  fToyNUACuts[kv][0] = 0.;
  fToyNUACuts[kv][1] = 0.;
  fToyNUACuts[kv][2] = 0.;
 } 

} // void AliAnalysisTaskMuPa::InitializeArraysForToyNUA()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForInternalValidation()
{
 // Initialize all arrays for internal validation.

 // Multiplicity range: min <= M < max
 fMultRangeInternalValidation[0] = 1000; // min
 fMultRangeInternalValidation[1] = 1001; // max

} // void AliAnalysisTaskMuPa::InitializeArraysForInternalValidation()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForTest0()
{
 // Initialize all arrays for Test0

 for(Int_t mo=0;mo<gMaxCorrelator;mo++) 
 { 
  for(Int_t mi=0;mi<gMaxIndex;mi++) 
  { 
   fTest0Labels[mo][mi] = NULL;
   for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
   { 
    fTest0Pro[mo][mi][v] = NULL;
   } 
  }
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForTest0()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForCorrelationsHistograms()
{
 // Initialize all arrays for histograms holding correlations.

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

 // a) Centrality and multiplicity;
 // b) Kinematics for specified filter bits;
 // c) Check for self-correlations;
 // d) Particle cut counters;
 // e) Trigger counter.

 // a) Centrality and multiplicity:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   fQACentralityHist[ce1][ba] = NULL;
   for(Int_t ce2=0;ce2<gCentralityEstimators;ce2++)
   {
    fQACentralityCorrHist[ce1][ce2][ba] = NULL;
   }
  }
  fQAMultiplicityCorrHist[ba] = NULL;

  for(Int_t gc=0;gc<gGenericCorrelations;gc++)
  {
   fQAGenericCorrHist[gc][ba] = NULL;
  }
 } //for(Int_t ba=0;ba<2;ba++)

 // b) Kinematics for specified filter bits (use in combination with SetQAFilterBits(...))
 for(Int_t fb=0;fb<gFilterBits;fb++)
 {
  for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4
  {
   fQAKinematicsFilterBits[fb][kv] = NULL;   
  } 
 }

 // c) Check for self-correlations:
 for(Int_t sc=0;sc<gQASelfCorrelations;sc++)
 {
  fQASelfCorrelations[sc] = NULL;
  fQASimRecoSelfCorrelations[sc] = NULL;
 }

 // d) Particle cut counters:
 for(Int_t rs=0;rs<2;rs++)
 {
  fQAParticleCutCounter[rs] = NULL;
 }

 // e) Trigger counter:
 for(Int_t ba=0;ba<2;ba++)
 {
  fQATrigger[ba] = NULL;
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForQAHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForControlEventHistograms()
{
 // Initialize all arrays for control event histograms. Cuts hardwired here are default event cuts.

 // a) Multiplicities.
 // b) Centrality;
 // c) Vertex;
 // d) Remaining event histograms;
 // e) Generic correlations.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Multiplicities:
 fMultiplicityBins[0] = 3000.;
 fMultiplicityBins[1] = 0.;
 fMultiplicityBins[2] = 3000.;
 fSelectedTracksCuts[0] = 0;
 fSelectedTracksCuts[1] = 1e6;
 for(Int_t m=0;m<gCentralMultiplicity;m++)
 { 
  for(Int_t ba=0;ba<2;ba++) // before/after cuts
  { 
   fCentralMultiplicityHist[m][ba] = NULL;
  }
 } // for(Int_t m=0;m<gCentralMultiplicity;m++)

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
 for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++) // first centrality estimator
 { 
  for(Int_t ce2=0;ce2<gCentralityEstimators;ce2++) // second centrality estimator
  {
   fUseCentralityCorrelationsCuts[ce1][ce2] = kFALSE; 
   fCentralityCorrelationsCuts[ce1][ce2] = 100.; 
  } // for(Int_t ce2=0;ce2<gCentralityEstimators;ce2++) // second centrality estimator
 } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++) // first centrality estimator
 
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
 fVertexCuts[X][0] = -1.;
 fVertexCuts[X][1] = 1.;
 fVertexCuts[Y][0] = -1.;
 fVertexCuts[Y][1] = 1.;
 fVertexCuts[Z][0] = -10.;
 fVertexCuts[Z][1] = 10.;
 fUseVertexCuts[X] = kFALSE;
 fUseVertexCuts[Y] = kFALSE;
 fUseVertexCuts[Z] = kFALSE;

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
 fEventBins[MagneticField][0] = 20;
 fEventBins[MagneticField][1] = -20.;
 fEventBins[MagneticField][2] = 10.;
 fEventBins[PrimaryVertex][0] = 1;
 fEventBins[PrimaryVertex][1] = 0.;
 fEventBins[PrimaryVertex][2] = 1.;

 // default cuts:
 fEventCuts[MagneticField][0] = -5.5;
 fEventCuts[MagneticField][1] = 5.5;
 fUseEventCuts[MagneticField] = kFALSE;
 fEventCuts[PrimaryVertex][0] = -44;  
 fEventCuts[PrimaryVertex][1] = -44;
 fUseEventCuts[PrimaryVertex] = kFALSE;

 // e) Generic correlations:
 for(Int_t gc=0;gc<gGenericCorrelations;gc++)
 {
  fUseGenericCorrelationsCuts[gc] = kFALSE;
 }

} // void AliAnalysisTaskMuPa::InitializeArraysForControlEventHistograms()

//================================================================================================================

void AliAnalysisTaskMuPa::InitializeArraysForControlParticleHistograms()
{
 // Initialize all arrays for control particle histograms. Cuts hardwired here are default particle cuts. 

 // a) Kinematics; 
 // b) DCA;
 // c) Remaining particle histograms; 
 // d) The default particle cuts. 

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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

 // non-equal bins:
 for(Int_t kv=0;kv<gKinematicVariables;kv++)
 {
  fNonEqualKinematicsnBins[kv] = 0;
  fNonEqualKinematicsRanges[kv] = TArrayD(0,NULL); // yes!
  fUseNonEqualKinematicsBins[kv] = kFALSE;
 }

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

 fParticleBins[TPCnclsFractionShared][0] = 2000;
 fParticleBins[TPCnclsFractionShared][1] = 0.;
 fParticleBins[TPCnclsFractionShared][2] = 2.;

 fParticleBins[TPCNCrossedRows][0] = 200;
 fParticleBins[TPCNCrossedRows][1] = 0.;
 fParticleBins[TPCNCrossedRows][2] = 200.;

 fParticleBins[TPCChi2perNDF][0] = 1000;
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

 fParticleBins[TPCNclsF][0] = 100;
 fParticleBins[TPCNclsF][1] = 0.;
 fParticleBins[TPCNclsF][2] = 1.;

 fParticleBins[HasPointOnITSLayer][0] = 6;
 fParticleBins[HasPointOnITSLayer][1] = 0.;
 fParticleBins[HasPointOnITSLayer][2] = 6.;

 fParticleBins[IsGlobalConstrained][0] = 2;
 fParticleBins[IsGlobalConstrained][1] = 0.;
 fParticleBins[IsGlobalConstrained][2] = 2.;

 // d) The default particle cuts: 
 //  These cuts are used by default, if you want other cuts, indicate that via dedicated setters.
 //  d1) Kinematics:
 //    default phi cuts:
 fKinematicsCuts[PHI][0] = fKinematicsBins[PHI][1];
 fKinematicsCuts[PHI][1] = fKinematicsBins[PHI][2];
 fUseKinematicsCuts[PHI] = kFALSE;

 //    default pt cuts:
 fKinematicsCuts[PT][0] = fKinematicsBins[PT][1];
 fKinematicsCuts[PT][1] = fKinematicsBins[PT][2];
 fUseKinematicsCuts[PT] = kFALSE;

 //    default eta cuts:
 fKinematicsCuts[ETA][0] = fKinematicsBins[ETA][1];
 fKinematicsCuts[ETA][1] = fKinematicsBins[ETA][2];
 fUseKinematicsCuts[ETA] = kFALSE;

 //    default e cuts:
 fKinematicsCuts[E][0] = fKinematicsBins[E][1];
 fKinematicsCuts[E][1] = fKinematicsBins[E][2];
 fUseKinematicsCuts[E] = kFALSE;

 //    default charge cuts:
 fKinematicsCuts[CHARGE][0] = fKinematicsBins[CHARGE][1];
 fKinematicsCuts[CHARGE][1] = fKinematicsBins[CHARGE][2];
 fUseKinematicsCuts[CHARGE] = kFALSE;

 //  d2) DCA:
 //    default DCA xy cuts:
 fDCACuts[0][0] = fDCABins[0][1];
 fDCACuts[0][1] = fDCABins[0][2];
 //    default DCA z cuts:
 fDCACuts[1][0] = fDCABins[1][1];
 fDCACuts[1][1] = fDCABins[1][2];
 //    corresponding booleans:
 fUseDCACuts[0] = kFALSE;
 fUseDCACuts[1] = kFALSE;

 //  d3) Remaining particle cuts:
 for(Int_t t=0;t<gParticleHistograms;t++)
 {
  fParticleCuts[t][0] = fParticleBins[t][1];
  fParticleCuts[t][1] = fParticleBins[t][2];
  fUseParticleCuts[t] = kFALSE;
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
 // Check before bookings if all the values the user has provided via setters make sense.

 // a) Multiplicity;
 // b) Centrality.
 // c) Centrality weights;
 // d) Toy NUA;
 // e) Supported triggers;
 // f) Internal validation.

 Green(__PRETTY_FUNCTION__);

 // a) Multiplicity:
 if(fSelectedTracksCuts[0]<0){cout<<__LINE__<<endl;exit(1);} 
 if(fSelectedTracksCuts[1]<=fSelectedTracksCuts[0]){cout<<__LINE__<<endl;exit(1);}
 if(fSelectedTracksCuts[0]<fMaxCorrelator){cout<<__LINE__<<endl;exit(1);} 

 // b) Centrality:
 if(((Int_t)fCentralityBins[0])<=0){cout<<__LINE__<<endl;exit(1);}
 if(fCentralityBins[1]<0. || fCentralityBins[1] > 100.){cout<<__LINE__<<endl;exit(1);}
 if(fCentralityBins[2]<0. || fCentralityBins[2] > 100.){cout<<__LINE__<<endl;exit(1);} 
 if(fCentralityBins[2]<=fCentralityBins[1]){cout<<__LINE__<<endl;exit(1);}
 if(!(fCentralityEstimator.EqualTo("V0M")||fCentralityEstimator.EqualTo("SPDTracklets")||
      fCentralityEstimator.EqualTo("CL0")||fCentralityEstimator.EqualTo("CL1"))){cout<<__LINE__<<endl;exit(1);}

 // c) Centrality weights:
 Int_t sum = 0;
 for(Int_t ce=0;ce<gCentralityEstimators;ce++)
 {
  sum += (Int_t) fUseCentralityWeights[ce];
 }
 if(sum>1){Red("\n Usage of only one centrality weight is supported at the moment.\n");cout<<__LINE__<<endl;exit(1);}

 // d) Toy NUA:
 if(fUseToyNUA[PHI])
 {
  if(fToyNUACuts[PHI][1]<0.||fToyNUACuts[PHI][2]>TMath::TwoPi()){Red("\n In ToyNUA phi range is out of boundaries\n");cout<<__LINE__<<endl;exit(1);}
 }
 if(fUseToyNUA[PT])
 {
  if(fToyNUACuts[PT][1]<0.||fToyNUACuts[PT][2]>100.){Red("\n In ToyNUA pt range is out of boundaries\n");cout<<__LINE__<<endl;exit(1);}
 }
 if(fUseToyNUA[ETA])
 {
  if(fToyNUACuts[ETA][1]<-1.||fToyNUACuts[ETA][2]>1.){Red("\n In ToyNUA eta range is out of boundaries\n");cout<<__LINE__<<endl;exit(1);}
 }
 // TBI 20210810 add similar protection also for E and CHARGE

 // e) Supported triggers:
 if(fUseTrigger && !(fTrigger.EqualTo("kAny") || fTrigger.EqualTo("kMB") || fTrigger.EqualTo("kINT7") || fTrigger.EqualTo("kCentral") || fTrigger.EqualTo("kSemiCentral") || 
      fTrigger.EqualTo("kCentral_kMB") || fTrigger.EqualTo("kSemiCentral_kMB")))
 {
  Red(Form("The trigger \"%s\" is not supported (yet).",fTrigger.Data()));
  cout<<__LINE__<<endl;
  exit(1);
 } 

 // f) Internal validation:
 if(fUseInternalValidation)
 {
  if(fnEventsInternalValidation < 1){cout<<__LINE__<<endl;exit(1);}
  if(fMultRangeInternalValidation[0] >= fMultRangeInternalValidation[1]){cout<<__LINE__<<endl;exit(1);}
  if(!fInternalValidationHarmonics){cout<<__LINE__<<endl;exit(1);}
  for(Int_t h=0;h<fInternalValidationHarmonics->GetSize();h++)
  {
   if(TMath::Abs(fInternalValidationHarmonics->GetAt(h)) > 0.5)
   {
    cout<<Form("v_%d = %f is too large to be taken seriously.",h+1,fInternalValidationHarmonics->GetAt(h))<<endl;
    cout<<__LINE__<<endl;exit(1);
   }
  }
  if(fRescaleWithTheoreticalInput && fCalculateNestedLoops){cout<<__LINE__<<endl;exit(1);}
 } 

 //Green("=> Done with InsanityChecks()!");

} // void AliAnalysisTaskMuPa::InsanityChecks()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookBaseProfile()
{
 // Book base profile which keeps flags relevant for the whole analysis.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 fBasePro = new TProfile("fBasePro","flags for the whole analysis",13,0.,13.);
 fBasePro->SetStats(kFALSE);
 fBasePro->SetLineColor(COLOR);
 fBasePro->SetFillColor(FILLCOLOR);
 fBasePro->GetXaxis()->SetBinLabel(1,Form("fTaskName = %s",fTaskName.Data()));
 fBasePro->GetXaxis()->SetBinLabel(2,Form("fDataTakingPeriod = %s",fDataTakingPeriod.Data()));
 fBasePro->GetXaxis()->SetBinLabel(3,Form("fAODNumber = %s",fAODNumber.Data()));
 fBasePro->GetXaxis()->SetBinLabel(4,Form("fRunNumber = %s",fRunNumber.Data()));
 fBasePro->GetXaxis()->SetBinLabel(5,"fFillQAhistograms"); fBasePro->Fill(3.5,fFillQAHistograms);
 fBasePro->GetXaxis()->SetBinLabel(6,"fFillQAhistogramsAll"); fBasePro->Fill(4.5,fFillQAHistogramsAll);
 fBasePro->GetXaxis()->SetBinLabel(7,"fTerminateAfterQA"); fBasePro->Fill(5.5,fTerminateAfterQA);
 fBasePro->GetXaxis()->SetBinLabel(8,"fVerbose"); fBasePro->Fill(6.5,fVerbose);
 fBasePro->GetXaxis()->SetBinLabel(9,"fRealData"); fBasePro->Fill(7.5,fRealData);
 fBasePro->GetXaxis()->SetBinLabel(10,"fUseFisherYates"); fBasePro->Fill(8.5,fUseFisherYates);
 fBasePro->GetXaxis()->SetBinLabel(11,"fRandomSeed"); fBasePro->Fill(9.5,fRandomSeed);
 fBasePro->GetXaxis()->SetBinLabel(12,Form("fTrigger = %s",fTrigger.Data()));
 fBasePro->GetXaxis()->SetBinLabel(13,"fUseTrigger"); fBasePro->Fill(11.5,fUseTrigger);
 fBaseList->Add(fBasePro);

} // void AliAnalysisTaskMuPa::BookBaseProfile()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fBaseList.

 // a) QA histograms;
 // b) Control event histograms;
 // c) Control particle histograms;
 // d) Q-vectors;
 // e) Weights;
 // f) Centrality weights;
 // g) Correlations;
 // h) Nested loops;
 // i) Toy NUA;
 // j) Internal validation;
 // k) Test0.

 // *) Book and nest lists for final results.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 if(!fBaseList){cout<<__LINE__<<endl;exit(1);}

 // a) QA histograms:
 fQAList = new TList();
 fQAList->SetName("QA");
 fQAList->SetOwner(kTRUE);
 fBaseList->Add(fQAList);

 // b) Control histograms:
 fControlEventHistogramsList = new TList();
 fControlEventHistogramsList->SetName("ControlEventHistograms");
 fControlEventHistogramsList->SetOwner(kTRUE);
 fBaseList->Add(fControlEventHistogramsList);

 // c) Control particle histograms;
 fControlParticleHistogramsList = new TList();
 fControlParticleHistogramsList->SetName("ControlParticleHistograms");
 fControlParticleHistogramsList->SetOwner(kTRUE);
 fBaseList->Add(fControlParticleHistogramsList);

 // d) Q-vectors:
 fQvectorList = new TList();
 fQvectorList->SetName("Q-vectors");
 fQvectorList->SetOwner(kTRUE);
 fBaseList->Add(fQvectorList);

 // e) Weights:
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);
 fBaseList->Add(fWeightsList);

 // f) Centrality weights:
 fCentralityWeightsList = new TList();
 fCentralityWeightsList->SetName("CentralityWeights");
 fCentralityWeightsList->SetOwner(kTRUE);
 fBaseList->Add(fCentralityWeightsList);

 // g) Correlations:
 fCorrelationsList = new TList();
 fCorrelationsList->SetName("Correlations");
 fCorrelationsList->SetOwner(kTRUE);
 fBaseList->Add(fCorrelationsList);

 // h) Nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("NestedLoops");
 fNestedLoopsList->SetOwner(kTRUE);
 fBaseList->Add(fNestedLoopsList);

 // i) Toy NUA:
 fToyNUAList = new TList();
 fToyNUAList->SetName("ToyNUA");
 fToyNUAList->SetOwner(kTRUE);
 fBaseList->Add(fToyNUAList);

 // j) Internal validation:
 fInternalValidationList = new TList();
 fInternalValidationList->SetName("InternalValidation");
 fInternalValidationList->SetOwner(kTRUE);
 fBaseList->Add(fInternalValidationList);

 // k) Test0:
 fTest0List = new TList();
 fTest0List->SetName("Test0");
 fTest0List->SetOwner(kTRUE);
 fBaseList->Add(fTest0List);

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
 // e) Anomalous events;
 // f) Check for self-correlations;
 // g) Event cut counter;
 // h) Particle cut counter;
 // i) Trigger counter.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fQAPro = new TProfile("fQAPro","flags for QA histograms",1,0.,1.);
 fQAPro->SetStats(kFALSE);
 fQAList->Add(fQAPro);

 if(!fFillQAHistograms){return;}

 // b) Common local style and labels:
 Int_t bLineColor = kGreen+2;
 Int_t bFillColor = kGreen-10;
 TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
 TString sba[2] = {"before event cuts","after event cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString skv[gKinematicVariables] = {"#varphi","p_{T}","#eta","energy","charge"};
 TString sae[gQAAnomalousEvents] = {"|vertex| = 0."};
 TString ssc[gQASelfCorrelations] = {"#varphi","p_{T}","#eta"};
 TString secc[gQAEventCutCounter] = {"nEvts","minSelected","maxSelected","minCent","maxCent","!avtx","minNContr","maxNContr","fMinVertexDistance",
                                     "minVtxX","maxVtxX","minVtxY","maxVtxY","minVtxZ","maxVtxZ",
                                     "centCorrCut[0][1]","centCorrCut[0][2]","centCorrCut[0][3]","centCorrCut[1][2]","centCorrCut[1][3]","centCorrCut[2][3]","centFlattening","nSel.vsNcontr."};
 TString spcc[gQAParticleCutCounter] = {"fFilterBit","kPrimary","phi_min","phi_max","pt_min","pt_max","eta_min","eta_max","e_min","e_max","charge","charge_min","charge_max",
                                        "dca_xy_min","dca_xy_max","dca_z_min","dca_z_max","TPCNcls_min","TPCNcls_max","TPCnclsS_min","TPCnclsS_max",
                                        "TPCnclsFractionShared_min","TPCnclsFractionShared_max","TPCNCrossedRows_min","TPCNCrossedRows_max","Chi2perNDF_min","Chi2perNDF_max",
                                        "TPCFoundFraction_min","TPCFoundFraction_max","Chi2TPCConstrainedVsGlobal_min","Chi2TPCConstrainedVsGlobal_max",
                                        "ITSNcls_min","ITSNcls_max","ITSChi2perNDF_min","ITSChi2perNDF_max","TPCNclsF_min","TPCNclsF_max","fAtLeastOnePointInTheSPD","fIgnoreGlobalConstrained"};

 for(Int_t ba=0;ba<2;ba++)
 {
  if(!fFillQAHistogramsAll){break;} 
  // Centrality correlations:
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   fQACentralityHist[ce1][ba] = new TH1D(Form("fQACentralityHist[%d][%d]",ce1,ba),Form("%s,%s",sce[ce1].Data(),sba[ba].Data()),(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]);
   //fQACentralityHist[ce1][ba]->SetStats(kFALSE);
   fQACentralityHist[ce1][ba]->SetLineColor(fBeforeAfterColor[ba]);
   fQACentralityHist[ce1][ba]->SetFillColor(fBeforeAfterColor[ba]-10);
   fQACentralityHist[ce1][ba]->GetXaxis()->SetTitle("centrality");
   fQAList->Add(fQACentralityHist[ce1][ba]);

   // 2D correlation plot, only upper diagonal is booked:
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    fQACentralityCorrHist[ce1][ce2][ba] = new TH2D(Form("fQACentralityCorrHist[%d][%d][%d]",ce1,ce2,ba),Form("%s vs. %s, %s",sce[ce1].Data(),sce[ce2].Data(),sba[ba].Data()),
                                          100,0.,100.,100,0.,100.);
    fQACentralityCorrHist[ce1][ce2][ba]->SetOption("col");
    fQACentralityCorrHist[ce1][ce2][ba]->GetXaxis()->SetTitle(Form("centrality %s",sce[ce1].Data()));
    fQACentralityCorrHist[ce1][ce2][ba]->GetYaxis()->SetTitle(Form("centrality %s",sce[ce2].Data()));
    fQAList->Add(fQACentralityCorrHist[ce1][ce2][ba]);
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)

  // Multiplicity correlations:
  fQAMultiplicityCorrHist[ba] = new TH2D(Form("fQAMultiplicityCorrHist[%d]",ba),Form("reference mult vs. nSelectedTracks, %s",sba[ba].Data()),
                                        (Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2],(Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2]);
  fQAMultiplicityCorrHist[ba]->SetOption("col");
  fQAMultiplicityCorrHist[ba]->GetXaxis()->SetTitle("RefMultComb08");
  fQAMultiplicityCorrHist[ba]->GetYaxis()->SetTitle("nSelectedTracks");
  fQAList->Add(fQAMultiplicityCorrHist[ba]);

  // Generic correlations:
  // 0: fSelectedTracks vs. avtx->GetNContributors()
  fQAGenericCorrHist[0][ba] = new TH2D(Form("fQAGenericCorrHist[0][%d]",ba),"fSelectedTracks vs. avtx->GetNContributors()",
                                        (Int_t)fMultiplicityBins[0],(Int_t)fMultiplicityBins[1],(Int_t)fMultiplicityBins[2],
                                        (Int_t)fNContributorsBins[0],fNContributorsBins[1],fNContributorsBins[2]);
  fQAGenericCorrHist[0][ba]->GetXaxis()->SetTitle("fSelectedTracks");
  fQAGenericCorrHist[0][ba]->GetYaxis()->SetTitle("avtx->GetNContributors()");

  // 1: ...

  // Common booking for generic correlations:
  for(Int_t gc=0;gc<gGenericCorrelations;gc++)
  {
   fQAGenericCorrHist[gc][ba]->SetOption("col");
   fQAGenericCorrHist[gc][ba]->SetLineColor(fBeforeAfterColor[ba]);
   fQAGenericCorrHist[gc][ba]->SetFillColor(fBeforeAfterColor[ba]-10);
   fQAList->Add(fQAGenericCorrHist[gc][ba]);
  }
 } // for(Int_t ba=0;ba<2;ba++)

 // Particles:
 if(fFillQAHistogramsAll)
 {
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
 } 

 if(fFillQAHistogramsAll)
 {
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
 }

 // b) Kinematics for specified filter bits (use in combination with SetQAFilterBits(...))
 for(Int_t fb=0;fb<fQAFilterBits->GetSize();fb++)
 {
  if(!fFillQAHistogramsAll){break;} 
  for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4
  {
   fQAKinematicsFilterBits[fb][kv] = new TH1D(Form("fQAKinematicsFilterBits[%d][%d]",fb,kv),Form("Filter bit: %d, %s",(Int_t)fQAFilterBits->GetAt(fb),skv[kv].Data()),(Int_t)fKinematicsBins[kv][0],fKinematicsBins[kv][1],fKinematicsBins[kv][2]);
   fQAKinematicsFilterBits[fb][kv]->SetXTitle(skv[kv].Data());
   fQAKinematicsFilterBits[fb][kv]->SetLineColor(COLOR);
   fQAKinematicsFilterBits[fb][kv]->SetFillColor(FILLCOLOR);
   fQAKinematicsFilterBits[fb][kv]->SetMinimum(0.);
   fQAList->Add(fQAKinematicsFilterBits[fb][kv]);   
  } 
 }

 // e) Anomalous events:
 // 0 : count events with |vertex| = 0
 // 1 : ...
 if(fFillQAHistogramsAll)
 {
  fQAAnomalousEvents = new TH1I("fQAAnomalousEvents","counter",gQAAnomalousEvents,0,gQAAnomalousEvents);
  fQAAnomalousEvents->SetLineColor(COLOR);
  fQAAnomalousEvents->SetFillColor(FILLCOLOR);
  for(Int_t ae=1;ae<=gQAAnomalousEvents;ae++)
  {
   fQAAnomalousEvents->GetXaxis()->SetBinLabel(ae,sae[ae-1].Data());
  }
  fQAAnomalousEvents->SetMinimum(0.);
  fQAList->Add(fQAAnomalousEvents);
 } 

 // f) Check for self-correlations:
 for(Int_t sc=0;sc<gQASelfCorrelations;sc++)
 {
  if(!fQACheckSelfCorrelations){break;}

  // Self-correlations in reconstructed sample:
  fQASelfCorrelations[sc] = new TH1D(Form("fQASelfCorrelations[%d]",sc),Form("Check for self-correlations in: %s_{1} - %s_{2}",ssc[sc].Data(),ssc[sc].Data()),200,-0.1,0.1); // TBI 20210526 hw limits
  fQASelfCorrelations[sc]->SetXTitle(ssc[sc].Data());
  fQASelfCorrelations[sc]->SetLineColor(COLOR);
  fQASelfCorrelations[sc]->SetFillColor(FILLCOLOR);
  fQASelfCorrelations[sc]->SetMinimum(0.);
  fQAList->Add(fQASelfCorrelations[sc]);

  // Self-correlations between simulated and reconstructed sample:
  fQASimRecoSelfCorrelations[sc] = new TH1D(Form("fQASimRecoSelfCorrelations[%d]",sc),Form("Check for sim-reco self-correlations in: %s_{1} - %s_{2}",ssc[sc].Data(),ssc[sc].Data()),200,-1.,1.); // TBI 20210526 hw limits
  fQASimRecoSelfCorrelations[sc]->SetXTitle(ssc[sc].Data());
  fQASimRecoSelfCorrelations[sc]->SetLineColor(COLOR);
  fQASimRecoSelfCorrelations[sc]->SetFillColor(FILLCOLOR);
  fQASimRecoSelfCorrelations[sc]->SetMinimum(0.);
  fQAList->Add(fQASimRecoSelfCorrelations[sc]);

 } // for(Int_t sc=0;sc<gQASelfCorrelations;sc++)

 // g) Event cut counter (calculate for the each event cut how many events are cut away):
 fQAEventCutCounter = new TH1I("fQAEventCutCounter","counter",gQAEventCutCounter,0,gQAEventCutCounter);
 fQAEventCutCounter->SetLineColor(COLOR);
 fQAEventCutCounter->SetFillColor(FILLCOLOR);
 for(Int_t ae=1;ae<=gQAEventCutCounter;ae++)
 {
  fQAEventCutCounter->GetXaxis()->SetBinLabel(ae,secc[ae-1].Data());
 }
 fQAEventCutCounter->SetMinimum(0.);
 fQAList->Add(fQAEventCutCounter);

 // sequential version:
 fQASequentialEventCutCounter = new TH1I("fQASequentialEventCutCounter","sequential counter",gQAEventCutCounter,0,gQAEventCutCounter);
 fQASequentialEventCutCounter->SetLineColor(COLOR);
 fQASequentialEventCutCounter->SetFillColor(FILLCOLOR);
 for(Int_t ae=1;ae<=gQAEventCutCounter;ae++)
 {
  fQASequentialEventCutCounter->GetXaxis()->SetBinLabel(ae,secc[ae-1].Data());
 }
 fQASequentialEventCutCounter->SetMinimum(0.);
 fQAList->Add(fQASequentialEventCutCounter);

 // h) Particle cut counter:
 for(Int_t rs=0;rs<2;rs++)
 {
  if(fRealData && 1==rs){continue;}
  fQAParticleCutCounter[rs] = new TH1I(Form("fQAParticleCutCounter[%d]",rs),Form("particle cut counter, %s",0==rs?"reconstructed":"simulated"),gQAParticleCutCounter,0,gQAParticleCutCounter);
  fQAParticleCutCounter[rs]->SetLineColor(COLOR);
  fQAParticleCutCounter[rs]->SetFillColor(FILLCOLOR);
  for(Int_t ae=1;ae<=gQAParticleCutCounter;ae++)
  {
   fQAParticleCutCounter[rs]->GetXaxis()->SetBinLabel(ae,spcc[ae-1].Data());
  }
  fQAParticleCutCounter[rs]->SetMinimum(0.);
  fQAList->Add(fQAParticleCutCounter[rs]);
 }

 // i) Trigger counter:
 for(Int_t ba=0;ba<2;ba++)
 {
  fQATrigger[ba] = new TH1I(Form("fQATrigger[%d]",ba),Form("trigger counter, %s",0==ba?"before":"after"),7,0,7);
  fQATrigger[ba]->SetLineColor(fBeforeAfterColor[ba]);
  fQATrigger[ba]->SetFillColor(fBeforeAfterColor[ba]-10);
  fQATrigger[ba]->GetXaxis()->SetBinLabel(1,"kAny");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(2,"kMB");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(3,"kINT7");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(4,"kCentral");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(5,"kSemiCentral");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(6,"kCentral_kMB");
  fQATrigger[ba]->GetXaxis()->SetBinLabel(7,"kSemiCentral_kMB");
  fQATrigger[ba]->SetMinimum(0.);
  fQAList->Add(fQATrigger[ba]);
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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fControlEventHistogramsPro = new TProfile("fControlEventHistogramsPro","flags for control event histograms",22,0.,22.);
 fControlEventHistogramsPro->SetStats(kFALSE);
 fControlEventHistogramsPro->SetLineColor(COLOR);
 fControlEventHistogramsPro->SetFillColor(FILLCOLOR);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(1,Form("fCentralityEstimator = %s",fCentralityEstimator.Data()));
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(2,"minSelected"); fControlEventHistogramsPro->Fill(1.5,fSelectedTracksCuts[0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(3,"maxSelected"); fControlEventHistogramsPro->Fill(2.5,fSelectedTracksCuts[1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(4,"minCent"); fControlEventHistogramsPro->Fill(3.5,fCentralityCuts[0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(5,"maxCent"); fControlEventHistogramsPro->Fill(4.5,fCentralityCuts[1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(6,"minNContr"); fControlEventHistogramsPro->Fill(5.5,fNContributorsCuts[0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(7,"maxNContr"); fControlEventHistogramsPro->Fill(6.5,fNContributorsCuts[1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(8,"fMinVertexDistance"); fControlEventHistogramsPro->Fill(7.5,fMinVertexDistance);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(9,"minVtxX"); fControlEventHistogramsPro->Fill(8.5,fVertexCuts[X][0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(10,"maxVtxX"); fControlEventHistogramsPro->Fill(9.5,fVertexCuts[X][1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(11,"minVtxY"); fControlEventHistogramsPro->Fill(10.5,fVertexCuts[Y][0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(12,"maxVtxY"); fControlEventHistogramsPro->Fill(11.5,fVertexCuts[Y][1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(13,"minVtxZ"); fControlEventHistogramsPro->Fill(12.5,fVertexCuts[Z][0]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(14,"maxVtxZ"); fControlEventHistogramsPro->Fill(13.5,fVertexCuts[Z][1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(15,"centCorrCut[0][1]"); fControlEventHistogramsPro->Fill(14.5,fCentralityCorrelationsCuts[0][1]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(16,"centCorrCut[0][2]"); fControlEventHistogramsPro->Fill(15.5,fCentralityCorrelationsCuts[0][2]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(17,"centCorrCut[0][3]"); fControlEventHistogramsPro->Fill(16.5,fCentralityCorrelationsCuts[0][3]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(18,"centCorrCut[1][2]"); fControlEventHistogramsPro->Fill(17.5,fCentralityCorrelationsCuts[1][2]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(19,"centCorrCut[1][3]"); fControlEventHistogramsPro->Fill(18.5,fCentralityCorrelationsCuts[1][3]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(20,"centCorrCut[2][3]"); fControlEventHistogramsPro->Fill(19.5,fCentralityCorrelationsCuts[2][3]);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(21,"fCentralityCorrelationCutVersion"); fControlEventHistogramsPro->Fill(20.5,fCentralityCorrelationCutVersion);
 fControlEventHistogramsPro->GetXaxis()->SetBinLabel(22,"fUseGenericCorrelationsCuts[0]"); fControlEventHistogramsPro->Fill(21.5,fUseGenericCorrelationsCuts[0]);
 fControlEventHistogramsList->Add(fControlEventHistogramsPro);

 // b) Common local labels:
 //    Remark: Keep them in sync with enums in .h
 TString sxyz[3] = {"x","y","z"};
 TString sba[2] = {"before event cuts","after event cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString stype[gEventHistograms] = {"MagneticField","PrimaryVertex"};
 TString smult[gCentralMultiplicity] = {"RefMultComb08"};

 // c) Multiplicity:
 fMultiplicityHist = new TH1D("fMultiplicityHist",Form("%s, multiplicity = sum of particle weights in Q-vector",fRunNumber.Data()),(Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2]);
 //fMultiplicityHist->SetStats(kFALSE);
 fMultiplicityHist->SetLineColor(COLOR);
 fMultiplicityHist->SetFillColor(FILLCOLOR);
 fMultiplicityHist->GetXaxis()->SetTitle("Multiplicity");
 fControlEventHistogramsList->Add(fMultiplicityHist);

 fSelectedTracksHist = new TH1I("fSelectedTracksHist",Form("%s, selected particles in Q-vector, after all event and particle cuts",fRunNumber.Data()),(Int_t)fMultiplicityBins[0],(Int_t)fMultiplicityBins[1],(Int_t)fMultiplicityBins[2]);
 //fSelectedTracksHist->SetStats(kFALSE);
 fSelectedTracksHist->SetLineColor(COLOR);
 fSelectedTracksHist->SetFillColor(FILLCOLOR);
 fSelectedTracksHist->GetXaxis()->SetTitle("integer counter of selected tracks ");
 fControlEventHistogramsList->Add(fSelectedTracksHist);

 for(Int_t m=0;m<gCentralMultiplicity;m++)
 {
  for(Int_t ba=0;ba<2;ba++) // before/after cuts
  { 
   fCentralMultiplicityHist[m][ba] = new TH1D(Form("fCentralMultiplicityHist[%d][%d]",m,ba),Form("%s, %s",fRunNumber.Data(),sba[ba].Data()),(Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2]);
   fCentralMultiplicityHist[m][ba]->SetStats(kFALSE);
   fCentralMultiplicityHist[m][ba]->SetLineColor(fBeforeAfterColor[ba]);
   fCentralMultiplicityHist[m][ba]->SetFillColor(fBeforeAfterColor[ba]-10);
   fCentralMultiplicityHist[m][ba]->GetXaxis()->SetTitle(smult[m].Data());
   fControlEventHistogramsList->Add(fCentralMultiplicityHist[m][ba]);
  }
 }

 // d) Centrality:
 for(Int_t ba=0;ba<2;ba++) // before/after cuts
 { 
  fCentralityHist[ba] = new TH1D(Form("fCentralityHist[%d]",ba),Form("%s, %s, %s",fRunNumber.Data(),sba[ba].Data(),fCentralityEstimator.Data()),(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]); // keep changes title in sync. with MakeCentralityWeights.C macro
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
   if(fRealData && 1==rs){continue;}
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    fVertexHist[ba][rs][xyz] = new TH1D(Form("fVertexHist[%d][%d][%d]",ba,rs,xyz),Form("%s, %s, %s",fRunNumber.Data(),sba[ba].Data(),srs[rs].Data()),(Int_t)fVertexBins[0],
                                        xyz==0||xyz==1 ? fVertexBins[1]/10. : fVertexBins[1], xyz==0||xyz==1 ? fVertexBins[2]/10 : fVertexBins[2]);
                                        // Above line: since spread of x and y components is order of magnitude smaller than z, the bin range is adjusted manually    
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
   if(fRealData && 1==rs){continue;}
   fNContributorsHist[ba][rs] = new TH1I(Form("fNContributorsHist[%d][%d]",ba,rs),Form("%s, %s, %s",fRunNumber.Data(),sba[ba].Data(),srs[rs].Data()),(Int_t)fNContributorsBins[0],fNContributorsBins[1],fNContributorsBins[2]); 
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
   fEventHistograms[ba][t] = new TH1D(Form("fEventHistograms[%d][%d]",ba,t),Form("%s, %s, %s",fRunNumber.Data(),sba[ba].Data(),stype[t].Data()),(Int_t)fEventBins[t][0],fEventBins[t][1],fEventBins[t][2]); 
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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fControlParticleHistogramsPro = new TProfile("fControlParticleHistogramsPro","flags for control particle histograms",4,0.,4.);
 fControlParticleHistogramsPro->SetLineColor(COLOR);
 fControlParticleHistogramsPro->SetFillColor(FILLCOLOR);
 fControlParticleHistogramsPro->SetStats(kFALSE);
 fControlParticleHistogramsPro->GetXaxis()->SetBinLabel(1,"fFilterBit"); fControlParticleHistogramsPro->Fill(0.5,fFilterBit);
 fControlParticleHistogramsPro->GetXaxis()->SetBinLabel(2,"fUseOnlyPrimaries"); fControlParticleHistogramsPro->Fill(1.5,fUseOnlyPrimaries);
 fControlParticleHistogramsPro->GetXaxis()->SetBinLabel(3,"fAtLeastOnePointInTheSPD"); fControlParticleHistogramsPro->Fill(2.5,fAtLeastOnePointInTheSPD);
 fControlParticleHistogramsPro->GetXaxis()->SetBinLabel(4,"fIgnoreGlobalConstrained"); fControlParticleHistogramsPro->Fill(3.5,fIgnoreGlobalConstrained);
 fControlParticleHistogramsList->Add(fControlParticleHistogramsPro);

 // b) Common local labels:
 //    Remark: Keep them in sync with enums in .h
 TString sxyTz[2] = {"xy","z"};
 TString sba[2] = {"before particle cuts","after particle cuts"};
 TString srs[2] = {"reconstructed","simulated"};
 TString skv[gKinematicVariables] = {"#varphi","p_{T}","#eta","energy","charge"};
 TString stype[gParticleHistograms] = {"TPCNcls","TPCnclsS","TPCnclsFractionShared","TPCNCrossedRows","TPCChi2perNDF","TPCFoundFraction","Chi2TPCConstrainedVsGlobal","ITSNcls","ITSChi2perNDF",
                                       "TPCNclsF","HasPointOnITSLayer","IsGlobalConstrained"};

 // c) Kinematics:
 for(Int_t ba=0;ba<2;ba++)
 {
  for(Int_t rs=0;rs<2;rs++)
  {
   if(fRealData && 1==rs){continue;}
   for(Int_t kv=0;kv<gKinematicVariables;kv++) // PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4 TBI 20210512 this is not enforced to be in sync with the definition of enums
   {
    fKinematicsHist[ba][rs][kv] = new TH1D(Form("fKinematicsHist[%d][%d][%d]",ba,rs,kv),Form("%s, %d, %s, %s",fRunNumber.Data(),fFilterBit,sba[ba].Data(),srs[rs].Data()),(Int_t)fKinematicsBins[kv][0],fKinematicsBins[kv][1],fKinematicsBins[kv][2]); 
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
   if(fRealData && 1==rs){continue;}
   for(Int_t xyTz=0;xyTz<2;xyTz++) 
   {
    fDCAHist[ba][rs][xyTz] = new TH1D(Form("fDCAHist[%d][%d][%d]",ba,rs,xyTz),Form("%s, %d, %s, %s",fRunNumber.Data(),fFilterBit,sba[ba].Data(),srs[rs].Data()),(Int_t)fDCABins[xyTz][0],fDCABins[xyTz][1],fDCABins[xyTz][2]); 
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
   if(fRealData && 1==rs){continue;}
   for(Int_t t=0;t<gParticleHistograms;t++) // type, see enum 'eParticle'
   {
    fParticleHist[ba][rs][t] = new TH1D(Form("fParticleHist[%d][%d][%d]",ba,rs,t),Form("%s, %d, %s, %s",fRunNumber.Data(),fFilterBit,sba[ba].Data(),srs[rs].Data()),(Int_t)fParticleBins[t][0],fParticleBins[t][1],fParticleBins[t][2]);  
    fParticleHist[ba][rs][t]->GetXaxis()->SetTitle(stype[t].Data());
    fParticleHist[ba][rs][t]->SetMinimum(0.);  
    fParticleHist[ba][rs][t]->SetLineColor(fBeforeAfterColor[ba]);
    fParticleHist[ba][rs][t]->SetFillColor(fBeforeAfterColor[ba]-10);
    fControlParticleHistogramsList->Add(fParticleHist[ba][rs][t]); 
   } // for(Int_t t=0;t<gParticleHistograms;t++)

   // fine tuning for some histograms:
   fParticleHist[ba][rs][IsGlobalConstrained]->GetXaxis()->SetBinLabel(1,"IsGlobalConstrained = kTRUE"); 
   fParticleHist[ba][rs][IsGlobalConstrained]->GetXaxis()->SetBinLabel(2,"IsGlobalConstrained = kFALSE"); 

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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro","flags for correlations",1,0.,1.);
 fCorrelationsFlagsPro->SetStats(kFALSE);
 fCorrelationsFlagsPro->GetXaxis()->SetLabelSize(0.05); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateCorrelations");  
 fCorrelationsFlagsPro->Fill(0.5,fCalculateCorrelations);
 fCorrelationsList->Add(fCorrelationsFlagsPro);

 if(!fCalculateCorrelations){return;}

 // b) Common local labels:
 TString oVariable[4] = {"#varphi_{1}-#varphi_{2}","#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fNestedLoopsFlagsPro = new TProfile("fNestedLoopsFlagsPro","flags for nested loops",2,0.,2.);
 fNestedLoopsFlagsPro->SetStats(kFALSE);
 fNestedLoopsFlagsPro->GetXaxis()->SetLabelSize(0.05); 
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateNestedLoops");  
 fNestedLoopsFlagsPro->Fill(0.5,fCalculateNestedLoops);
 fNestedLoopsFlagsPro->Fill(1.5,fCalculateCustomNestedLoop);
 fNestedLoopsList->Add(fNestedLoopsFlagsPro);

 if(!(fCalculateNestedLoops||fCalculateCustomNestedLoop)){return;}

 const Int_t maxSize = 20000;
 ftaNestedLoops[0] = new TArrayD(maxSize); // ebe container for azimuthal angles 
 ftaNestedLoops[1] = new TArrayD(maxSize); // ebe container for particle weights (product of all)  

 if(!fCalculateNestedLoops){return;} // TBI this is safe enough, I think?

 // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms())
 TString oVariable[4] = {"#varphi_{1}-#varphi_{2}","#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
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

void AliAnalysisTaskMuPa::BookToyNUAHistograms()
{
 // Book all toy NUA histograms.

 // a) Book the profile holding flags;

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fToyNUAFlagsPro = new TProfile("fToyNUAFlagsPro","flags for toy NUA",5,0.,5.);
 fToyNUAFlagsPro->SetStats(kFALSE);
 fToyNUAFlagsPro->SetLineColor(COLOR);
 fToyNUAFlagsPro->SetFillColor(FILLCOLOR);
 fToyNUAFlagsPro->GetXaxis()->SetLabelSize(0.04); 
 fToyNUAFlagsPro->GetXaxis()->SetBinLabel(1,"fUseToyNUA[PHI]"); fToyNUAFlagsPro->Fill(0.5,fUseToyNUA[PHI]);
 fToyNUAFlagsPro->GetXaxis()->SetBinLabel(2,"fUseToyNUA[PT]"); fToyNUAFlagsPro->Fill(1.5,fUseToyNUA[PT]);
 fToyNUAFlagsPro->GetXaxis()->SetBinLabel(3,"fUseToyNUA[ETA]"); fToyNUAFlagsPro->Fill(2.5,fUseToyNUA[ETA]);
 fToyNUAFlagsPro->GetXaxis()->SetBinLabel(4,"fUseToyNUA[E]"); fToyNUAFlagsPro->Fill(3.5,fUseToyNUA[E]);
 fToyNUAFlagsPro->GetXaxis()->SetBinLabel(5,"fUseToyNUA[CHARGE]"); fToyNUAFlagsPro->Fill(4.5,fUseToyNUA[CHARGE]);
 fToyNUAList->Add(fToyNUAFlagsPro);

} // void AliAnalysisTaskMuPa::BookToyNUAHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookInternalValidationHistograms()
{
 // Book all internal validation histograms.

 // a) Book the profile holding flags;

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fInternalValidationFlagsPro = new TProfile("fInternalValidationFlagsPro","flags for internal validation",7,0.,7.);
 fInternalValidationFlagsPro->SetStats(kFALSE);
 fInternalValidationFlagsPro->SetLineColor(COLOR);
 fInternalValidationFlagsPro->SetFillColor(FILLCOLOR);
 fInternalValidationFlagsPro->GetXaxis()->SetLabelSize(0.04); 
 fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(1,"fUseInternalValidation"); fInternalValidationFlagsPro->Fill(0.5,fUseInternalValidation);
 for(Int_t e=0;e<7;e++) // hardwired support for harmonics v1-v6
 {
  fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(e+2,Form("v_{%d}",e+1)); 
  if(fInternalValidationHarmonics && (e+1 <= fInternalValidationHarmonics->GetSize()))
  {
   fInternalValidationFlagsPro->Fill(e+1+0.5,fInternalValidationHarmonics->GetAt(e));
  }
  else
  {
   fInternalValidationFlagsPro->Fill(e+1+0.5,0.0); // default values of all harmonics are 0
  }
 } // for(Int_t e=0;e<fInternalValidationHarmonics->GetSize()) 
 fInternalValidationList->Add(fInternalValidationFlagsPro);

} // void AliAnalysisTaskMuPa::BookInternalValidationHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookTest0Histograms()
{
 // Book all Test0 histograms.

 // a) Book the profile holding flags;
 // b) Generate all labels;
 // c) ...

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Book the profile holding flags:
 fTest0FlagsPro = new TProfile("fTest0FlagsPro","flags for Test0",5,0.,5.);
 fTest0FlagsPro->SetStats(kFALSE);
 fTest0FlagsPro->GetXaxis()->SetLabelSize(0.04); 
 fTest0FlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateTest0"); fTest0FlagsPro->Fill(0.5,fCalculateTest0);
 fTest0List->Add(fTest0FlagsPro);

 if(!fCalculateTest0){return;}
 
 // b) Generate all labels:
 this->GenerateCorrelationsLabels();  

 // c) ...
 Int_t vvvariableNBins[3] = {1,(Int_t)fMultiplicityBins[0],(Int_t)fCentralityBins[0]};
 Double_t vvvariableMinMax[3][2] = { {0.,1.}, // integrated 
                                    {fMultiplicityBins[1],fMultiplicityBins[2]}, // multiplicity
                                    {fCentralityBins[1],fCentralityBins[2]} // centrality
                                   };
 TString vvVariable[3] = {"integrated","multiplicity","centrality"};

 for(Int_t mo=0;mo<gMaxCorrelator;mo++) 
 { 
  for(Int_t mi=0;mi<gMaxIndex;mi++) 
  { 
   if(fTest0Labels[mo][mi]) // book only explicitly requested correlators, other via external file, or with hardcoded loops 
   {
    for(Int_t v=0;v<3;v++) 
    { 
     fTest0Pro[mo][mi][v] = new TProfile(Form("fTest0Pro[%d][%d][%d]",mo,mi,v),fTest0Labels[mo][mi]->Data(),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
     fTest0Pro[mo][mi][v]->SetStats(kFALSE);
     fTest0Pro[mo][mi][v]->Sumw2();
     fTest0Pro[mo][mi][v]->GetXaxis()->SetTitle(vvVariable[v].Data());
     fTest0List->Add(fTest0Pro[mo][mi][v]);
    }
   } // if(fTest0Labels[mo][mi]) 
  }
 }

} // void AliAnalysisTaskMuPa::BookTest0Histograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 if(fVerbose){Green(__PRETTY_FUNCTION__);} 

} // void AliAnalysisTaskMuPa::BookFinalResultsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookQvectorHistograms()
{
 // Book all Q-vector histograms.

 // a) Book the profile holding flags;
 // b) ...

 if(fVerbose){Green(__PRETTY_FUNCTION__);} 

 // a) Book the profile holding flags:
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

 // b) ...

} // void AliAnalysisTaskMuPa::BookQvectorHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookWeightsHistograms()
{
 // Book all objects for particle weights.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Histograms.

 if(fVerbose){Green(__PRETTY_FUNCTION__);} 

 // a) Book the profile holding flags:
 fWeightsFlagsPro = new TProfile("fWeightsFlagsPro","flags for particle weights",3,0.,3.);
 fWeightsFlagsPro->SetStats(kFALSE);
 fWeightsFlagsPro->SetLineColor(COLOR);
 fWeightsFlagsPro->SetFillColor(FILLCOLOR);
 fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);  
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(1,"w_{#varphi}");  
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(2,"w_{p_{t}}");  
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(3,"w_{#eta}"); 
 for(Int_t w=0;w<gWeights;w++) // use weights [phi,pt,eta]
 { 
  if(fUseWeights[w])fWeightsFlagsPro->Fill(w+0.5,1.);
 }
 fWeightsList->Add(fWeightsFlagsPro);

 // b) Common local labels:
 TString sVariable[gWeights] = {"#varphi","p_{t}","#eta"}; // [phi,pt,eta,rapidity]
 TString sWeights[gWeights] = {"w_{#varphi}","w_{p_{t}}","w_{#eta}"};

 // c) Histograms:
 for(Int_t w=0;w<gWeights;w++) // use weights [phi,pt,eta]
 {
  if(!fUseWeights[w]){continue;}
  if(!fWeightsHist[w]) // yes, because these histos are cloned from the external ones, see SetWeightsHist(TH1D* const hist, const char *variable)
  {
   fWeightsHist[w] = new TH1D(Form("fWeightsHist[%d]",w),"",(Int_t)fKinematicsBins[w][0],fKinematicsBins[w][1],fKinematicsBins[w][2]);
   fWeightsHist[w]->SetTitle(Form("Particle weights for %s",sWeights[w].Data()));
   fWeightsHist[w]->SetStats(kFALSE);
   fWeightsHist[w]->GetXaxis()->SetTitle(sVariable[w].Data());
   fWeightsHist[w]->SetFillColor(FILLCOLOR);
   fWeightsHist[w]->SetLineColor(COLOR);
  }
  fWeightsList->Add(fWeightsHist[w]);
 } // for(Int_t w=0;w<gWeights;w++) // use weights [phi,pt,eta]

} // void AliAnalysisTaskMuPa::BookWeightsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::BookCentralityWeightsHistograms()
{
 // Book all objects for centrality weights.

 // a) Book the profile holding flags;
 // b) Common local labels;
 // c) Histograms.

 if(fVerbose){Green(__PRETTY_FUNCTION__);} 

 // a) Book the profile holding flags:
 fCentralityWeightsFlagsPro = new TProfile("fCentralityWeightsFlagsPro","flags for centrality weights",4,0.,4.);
 fCentralityWeightsFlagsPro->SetStats(kFALSE);
 fCentralityWeightsFlagsPro->SetLineColor(COLOR);
 fCentralityWeightsFlagsPro->SetFillColor(FILLCOLOR);
 fCentralityWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);  
 fCentralityWeightsFlagsPro->GetXaxis()->SetBinLabel(1,"w_{V0M}");   \
 fCentralityWeightsFlagsPro->GetXaxis()->SetBinLabel(2,"w_{SPDTracklets}");  
 fCentralityWeightsFlagsPro->GetXaxis()->SetBinLabel(3,"w_{CL0}"); 
 fCentralityWeightsFlagsPro->GetXaxis()->SetBinLabel(4,"w_{CL1}"); 
 for(Int_t w=0;w<gCentralityEstimators;w++) // use weights [phi,pt,eta]
 { 
  if(fUseCentralityWeights[w])fCentralityWeightsFlagsPro->Fill(w+0.5,1.);
 }
 fCentralityWeightsList->Add(fCentralityWeightsFlagsPro);

 // b) Common local labels:
 TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
 TString sCentralityWeights[gCentralityEstimators] = {"w_{V0M}","w_{SPDTracklets}","w_{CL0}","w_{CL1}"};

 // c) Histograms:
 for(Int_t w=0;w<gCentralityEstimators;w++) // use centrality weights ["V0M", "SPDTracklets", "CL0", "CL1"]
 {
  if(!fUseCentralityWeights[w]){continue;}
  if(!fCentralityWeightsHist[w]) // yes, because these histos are cloned from the exteral ones, see SetCentralityWeightsHist(TH1D* const hist, const char *estimator)
  {
   fCentralityWeightsHist[w] = new TH1D(Form("fCentralityWeightsHist[%d]",w),"",(Int_t)fCentralityBins[0],fCentralityBins[1],fCentralityBins[2]);
   fCentralityWeightsHist[w]->SetTitle(Form("Centrality weights for %s",sCentralityWeights[w].Data()));
   fCentralityWeightsHist[w]->SetStats(kFALSE);
   fCentralityWeightsHist[w]->GetXaxis()->SetTitle(sce[w].Data());
   fCentralityWeightsHist[w]->SetFillColor(FILLCOLOR);
   fCentralityWeightsHist[w]->SetLineColor(COLOR);
  }
  fCentralityWeightsList->Add(fCentralityWeightsHist[w]);
 } // for(Int_t w=0;w<gCentralityEstimators;w++) // use centrality weights ["V0M", "SPDTracklets", "CL0", "CL1"]

 // Quick sanity check for consistent binning:
 for(Int_t w=0;w<gCentralityEstimators;w++) // use centrality weights ["V0M", "SPDTracklets", "CL0", "CL1"]
 {
  if(!fUseCentralityWeights[w]){continue;}
  if((Int_t)(fCentralityWeightsHist[w]->GetNbinsX()) != (Int_t)(fCentralityHist[0]->GetNbinsX()))
  {
   cout<<Form("w = %d",w)<<endl;
   cout<<__LINE__<<endl;exit(1);
  }
 } // for(Int_t w=0;w<gCentralityEstimators;w++) // use centrality weights ["V0M", "SPDTracklets", "CL0", "CL1"]

} // void AliAnalysisTaskMuPa::BookCentralityWeightsHistograms()

//=======================================================================================================================

void AliAnalysisTaskMuPa::GetPointers(TList *baseList)
{
 // Get all pointers. This method is essential only for merging and boostrap.

 // a) Check the pointer for base list fBaseList;
 // b) Get pointer for profile holding internal flags and set again all flags;
 // c) Get pointers for control event histograms;
 // d) Get pointers control particle histograms;

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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

  // Count number of particles which survive all particle cuts, and are selected for the main analysis and added to Q-vectors:
  fSelectedTracks = 0;
  Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  {
   AliAODTrack *aodTrack = NULL;
   if(!fUseFisherYates)
   {
    aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to "a track" (i.e. any track)
   }
   else
   {
    aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack((Int_t)fRandomIndices->GetAt(iTrack)));
   }
   if(!aodTrack){continue;}
   if(!SurvivesParticleCuts(aodTrack)){continue;}
   fSelectedTracks++;
  } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 } // if(aAOD)

 return;

} // void AliAnalysisTaskMuPa::FilterEvent(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)
{
 // Fill control histograms before event cuts (ba = 0), or after (ba = 1). For reconstructed data (rs = 0), or simulated (rs = 1)
 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Vertex;
 // c) Centrality;
 // d) Reference multiplicity (maintained centrally);
 // e) Remaining event distributions.

 //if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  // b) Vertex: TH1D *fVertexHist[2][2][3]; //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) return; 
  if(fVertexHist[ba][rs][X]){fVertexHist[ba][rs][X]->Fill(avtx->GetX());}
  if(fVertexHist[ba][rs][Y]){fVertexHist[ba][rs][Y]->Fill(avtx->GetY());}
  if(fVertexHist[ba][rs][Z]){fVertexHist[ba][rs][Z]->Fill(avtx->GetZ());}
  if(fNContributorsHist[ba][rs]){fNContributorsHist[ba][rs]->Fill(avtx->GetNContributors());}

  // c) Centrality:
  if(fCentralityHist[ba]){fCentralityHist[ba]->Fill(fCentrality);}

  // d) Reference multiplicity (maintained centrally):
  AliAODHeader *aodheader = (AliAODHeader*)aAOD->GetHeader();
  if(!aodheader){cout<<__LINE__<<endl;exit(1);}
  if(fCentralMultiplicityHist[RefMultComb08][ba]){fCentralMultiplicityHist[RefMultComb08][ba]->Fill(aodheader->GetRefMultiplicityComb08());}

  // e) Remaining event distributions:
  if(fEventHistograms[ba][MagneticField]){fEventHistograms[ba][MagneticField]->Fill(aodheader->GetMagneticField());} // Solenoid Magnetic Field in kG
  if(fEventHistograms[ba][PrimaryVertex]){fEventHistograms[ba][PrimaryVertex]->Fill(0.44);} // here we only count # of events with valid pointer aAOD->GetPrimaryVertex()  
 } // if(aAOD)

 if(aMC)
 {
  // b) Vertex: TH1D *fVertexHist[2][2][3]; //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  //    TBI 20210802 this can be unified with AOD branch above
  AliAODVertex *avtx = (AliAODVertex*)aMC->GetPrimaryVertex();
  if(!avtx) return; 
  if(fVertexHist[ba][rs][X]){fVertexHist[ba][rs][X]->Fill(avtx->GetX());}
  if(fVertexHist[ba][rs][Y]){fVertexHist[ba][rs][Y]->Fill(avtx->GetY());}
  if(fVertexHist[ba][rs][Z]){fVertexHist[ba][rs][Z]->Fill(avtx->GetZ());}
  if(fNContributorsHist[ba][rs]){fNContributorsHist[ba][rs]->Fill(avtx->GetNContributors());}
 } // if(aMC)

} // AliAnalysisTaskMuPa::FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FillControlParticleHistograms(AliVParticle *vParticle, const Int_t ba, const Int_t rs)
{
 // Fill control histograms before particles cuts (ba = 0), or after (ba = 1). For reconstructed data (rs = 0), or simulated (rs = 1)
 // More detailed treatment of AOD track cuts can be found in Tasks/AliFlowTrackCuts.cxx
 
 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track;
 // b) Fill histograms for AOD track;
 // c) Fill histograms for AOD MC particle.

 //if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track:
 AliAODMCParticle *aodmcParticle = dynamic_cast<AliAODMCParticle*>(vParticle);
 AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(vParticle);
  
 if(rs == SIM && !aodmcParticle){cout<<__LINE__<<endl;exit(1);}
 if(rs == RECO && !aodTrack){cout<<__LINE__<<endl;exit(1);}

 // b) Fill histograms for AOD track:
 if(aodTrack)
 {
  // Kinematics:
  if(fKinematicsHist[ba][rs][PHI]){fKinematicsHist[ba][rs][PHI]->Fill(aodTrack->Phi());}
  if(fKinematicsHist[ba][rs][PT]){fKinematicsHist[ba][rs][PT]->Fill(aodTrack->Pt());}
  if(fKinematicsHist[ba][rs][ETA]){fKinematicsHist[ba][rs][ETA]->Fill(aodTrack->Eta());}
  if(fKinematicsHist[ba][rs][E]){fKinematicsHist[ba][rs][E]->Fill(aodTrack->E());}
  if(fKinematicsHist[ba][rs][CHARGE]){fKinematicsHist[ba][rs][CHARGE]->Fill(aodTrack->Charge());}

  // DCA:
  Float_t dcaXY = aodTrack->DCA(); 
  Float_t dcaZ = aodTrack->ZAtDCA(); 
  if(!aodTrack->TestBit(AliAODTrack::kIsDCA))
  {
   aodTrack->GetImpactParameters(dcaXY,dcaZ); 
  } 
  if(fDCAHist[ba][rs][0]){fDCAHist[ba][rs][0]->Fill(dcaXY);}
  if(fDCAHist[ba][rs][1]){fDCAHist[ba][rs][1]->Fill(dcaZ);}

  // Remaining:
  if(fParticleHist[ba][rs][TPCNcls]){fParticleHist[ba][rs][TPCNcls]->Fill(aodTrack->GetTPCNcls());}
  if(fParticleHist[ba][rs][TPCnclsS]){fParticleHist[ba][rs][TPCnclsS]->Fill(aodTrack->GetTPCnclsS());}
  if(fParticleHist[ba][rs][TPCnclsFractionShared])
  {
   if(TMath::Abs((Double_t)aodTrack->GetTPCncls())>0.)
   {
    fParticleHist[ba][rs][TPCnclsFractionShared]->Fill((Double_t)aodTrack->GetTPCnclsS()/(Double_t)aodTrack->GetTPCNcls()); // avoiding integer division truncation
   }
  }
  if(fParticleHist[ba][rs][TPCNCrossedRows]){fParticleHist[ba][rs][TPCNCrossedRows]->Fill(aodTrack->GetTPCNCrossedRows());}
  if(fParticleHist[ba][rs][TPCChi2perNDF])
  {
   if(aodTrack->GetTPCNcls()>0.)
   {
    // aodTrack->Chi2perNDF() is depreciated
    fParticleHist[ba][rs][TPCChi2perNDF]->Fill((Double_t)aodTrack->GetTPCchi2()/(Double_t)aodTrack->GetTPCNcls());
   }
  }
  if(fParticleHist[ba][rs][TPCFoundFraction])
  { 
   // "TPCFoundFraction" is defined as (Double_t)aodTrack->GetTPCNcls()/(Double_t)aodTrack->GetTPCNCrossedRows()
   fParticleHist[ba][rs][TPCFoundFraction]->Fill(aodTrack->GetTPCFoundFraction());
  }
  if(fParticleHist[ba][rs][Chi2TPCConstrainedVsGlobal]){fParticleHist[ba][rs][Chi2TPCConstrainedVsGlobal]->Fill(aodTrack->GetChi2TPCConstrainedVsGlobal());}
  if(fParticleHist[ba][rs][ITSNcls]){fParticleHist[ba][rs][ITSNcls]->Fill(aodTrack->GetITSNcls());}
  if(fParticleHist[ba][rs][ITSChi2perNDF]){if(TMath::Abs(aodTrack->GetITSNcls())>0.){fParticleHist[ba][rs][ITSChi2perNDF]->Fill(aodTrack->GetITSchi2()/aodTrack->GetITSNcls());}}
  if(fParticleHist[ba][rs][TPCNclsF]){fParticleHist[ba][rs][TPCNclsF]->Fill(aodTrack->GetTPCNclsF());}
  if(fParticleHist[ba][rs][HasPointOnITSLayer])
  {
   for(Int_t l=1;l<=6;l++) // loop over ITS layers
   {
    if(aodTrack->HasPointOnITSLayer(l)){fParticleHist[ba][rs][HasPointOnITSLayer]->Fill(l-0.5);} // fill in the center of the bin   
   }
  }
  if(fParticleHist[ba][rs][IsGlobalConstrained])
  {
   if(aodTrack->IsGlobalConstrained())
   {
    fParticleHist[ba][rs][IsGlobalConstrained]->Fill(0.44); // keep in sync. with the booking of this histogram
   } 
   else
   {
    fParticleHist[ba][rs][IsGlobalConstrained]->Fill(1.44); // keep in sync. with the booking of this histogram
   }  
  }
 } // if(aodTrack)

 // c) Fill histograms for AOD MC particle:
 if(aodmcParticle)
 {
  // Kinematics:
  if(fKinematicsHist[ba][rs][PHI]){fKinematicsHist[ba][rs][PHI]->Fill(aodmcParticle->Phi());}
  if(fKinematicsHist[ba][rs][PT]){fKinematicsHist[ba][rs][PT]->Fill(aodmcParticle->Pt());}
  if(fKinematicsHist[ba][rs][ETA]){fKinematicsHist[ba][rs][ETA]->Fill(aodmcParticle->Eta());}
  if(fKinematicsHist[ba][rs][E]){fKinematicsHist[ba][rs][E]->Fill(aodmcParticle->E());}
  if(fKinematicsHist[ba][rs][CHARGE]){fKinematicsHist[ba][rs][CHARGE]->Fill(aodmcParticle->Charge());} 
 } // if(aodmcParticle)

} // AliAnalysisTaskMuPa::FillControlParticleHistograms(AliVTrack *vTrack, const Int_t ba, const Int_t rs)

//=======================================================================================================================

Bool_t AliAnalysisTaskMuPa::SurvivesEventCuts(AliVEvent *ave)
{
 // Check if the current event survives event cuts. TBI 20210512 works only for 'reco' at the moment

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Only for QA count separately what each event cut is doing;
 // c) AOD.

 // a) Determine Ali{ESD,AOD}Event:
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) Only for QA count separately what each event cut is doing:
 if(fFillQAHistograms)
 {
  if(aAOD){this->EventCutCounter(aAOD); this->SequentialEventCutCounter(aAOD);}
 }

 // c) AOD:
 if(aAOD)
 {
  // Remark 0: Cut first on event-by-event quantities which we calculated in FilterEvent().
  //           At the moment, those are: fCentrality, fSelectedTracks

  // Remark 1: For each new cut you implement here, add corresponding support in EventCutCounter(...) and SequentialEventCutCounter(...)

  // Trigger:
  if(fUseTrigger)
  {
   UInt_t fSelectMask = fInputHandler->IsEventSelected();
   Bool_t bIsSelected = kFALSE;
   if(fTrigger.EqualTo("kAny"))
   {
    bIsSelected = fSelectMask & AliVEvent::kAny;
   }
   else if(fTrigger.EqualTo("kMB"))
   {
    bIsSelected = fSelectMask & AliVEvent::kMB;
   }
   else if(fTrigger.EqualTo("kINT7")) 
   {
    bIsSelected = fSelectMask & AliVEvent::kINT7;
   }
   else if(fTrigger.EqualTo("kCentral"))
   {
    bIsSelected = fSelectMask & AliVEvent::kCentral;
   }
   else if(fTrigger.EqualTo("kSemiCentral"))
   {
    bIsSelected = fSelectMask & AliVEvent::kSemiCentral;
   }
   else if(fTrigger.EqualTo("kCentral_kMB"))
   {
    bIsSelected = (fSelectMask & AliVEvent::kCentral) || (fSelectMask & AliVEvent::kMB);
   }
   else if(fTrigger.EqualTo("kSemiCentral_kMB"))
   {
    bIsSelected = (fSelectMask & AliVEvent::kSemiCentral) || (fSelectMask & AliVEvent::kMB);
   }
   else 
   {
    cout<<__LINE__<<endl; exit(1);
   }
   if(!bIsSelected){ return kFALSE; }
  } // if(fUseTrigger)

  // Centrality:
  if(fUseCentralityCuts)
  {
   if(fCentrality < fCentralityCuts[0]) return kFALSE;
   if(fCentrality > fCentralityCuts[1]) return kFALSE;
  } 

  // Selected tracks used to calculate Q-vectors, after all particle cuts have been applied:
  if(fUseSelectedTracksCuts)
  {
   if(fSelectedTracks < fSelectedTracksCuts[0]) return kFALSE;
   if(fSelectedTracks > fSelectedTracksCuts[1]) return kFALSE;
  }

  // Centrality correlations cuts:
  // TBI 20210727 this can be optimized further, but nevermind now, since by default it will be used in any case
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    if(fUseCentralityCorrelationsCuts[ce1][ce2])
    {
     Double_t centrality1 = ams->GetMultiplicityPercentile(sce[ce1].Data());
     Double_t centrality2 = ams->GetMultiplicityPercentile(sce[ce2].Data());
     if(centrality1 > 0. && centrality2 > 0.) 
     {
      switch(fCentralityCorrelationCutVersion)
      {
       case 0: // relative 
        if(TMath::Abs((centrality1-centrality2)/(centrality1+centrality2)) > fCentralityCorrelationsCuts[ce1][ce2]) return kFALSE;
       break;
       case 1: // absolute
        if(TMath::Abs((centrality1-centrality2)) > fCentralityCorrelationsCuts[ce1][ce2]) return kFALSE;
       break;
       default:
        cout<<__LINE__<<endl; exit(1);
       break; 
      } // switch(fCentralityCorrelationCutVersion)           
     } // if(centrality1 > 0. && centrality2 > 0.) 
    } 
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)

  // Vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) return kFALSE; 
  if(fUseNContributorsCuts)
  {
   if((Int_t)avtx->GetNContributors()<fNContributorsCuts[0]) return kFALSE;
   if((Int_t)avtx->GetNContributors()>fNContributorsCuts[1]) return kFALSE;
  }

  if(fUseMinVertexDistanceCut)
  {
   if(sqrt(pow(avtx->GetX(),2.) + pow(avtx->GetY(),2.) + pow(avtx->GetY(),2.)) < fMinVertexDistance)
   { 
    if(fQAAnomalousEvents){fQAAnomalousEvents->Fill(0.5);} // |vertex| = 0.
    return kFALSE;
   }
  }

  if(fUseVertexCuts[X])
  {
   if(avtx->GetX() < fVertexCuts[X][0]) return kFALSE;
   if(avtx->GetX() > fVertexCuts[X][1]) return kFALSE; 
  }
  if(fUseVertexCuts[Y])
  {
   if(avtx->GetY() < fVertexCuts[Y][0]) return kFALSE;
   if(avtx->GetY() > fVertexCuts[Y][1]) return kFALSE;
  }
  if(fUseVertexCuts[Z])
  {
   if(avtx->GetZ() < fVertexCuts[Z][0]) return kFALSE;
   if(avtx->GetZ() > fVertexCuts[Z][1]) return kFALSE;
  }

  // Centrality weights (flattening): 
  // Remark: since I am getting centrality weights from centrality distribution after the events cuts, flattening is applied here after all other event cuts:
  if(!ams){cout<<__LINE__<<endl;exit(1);} // pointer was obtained previously
  Double_t centralityWeight = -44.;
  for(Int_t ce=0;ce<gCentralityEstimators;ce++)
  {
   if(fUseCentralityWeights[ce])
   {
    centralityWeight = CentralityWeight(ams->GetMultiplicityPercentile(sce[ce].Data()),sce[ce].Data());
    if(gRandom->Uniform(0,1) > centralityWeight) return kFALSE; // yes, since centralityWeight is normalized probability (see CentralityWeight(...))
   } 
  }

  if(fUseGenericCorrelationsCuts[0])
  {
   if(avtx)
   {
    if(fSelectedTracks > fSelContrTreshold * (Int_t)avtx->GetNContributors()) return kFALSE;
   }
  }

 } // if(aAOD)

 return kTRUE;

} // AliAnalysisTaskMuPa::SurvivesEventCuts(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMuPa::EventCutCounter(AliVEvent *ave)
{
 // Simple QA function to summarize what each event cut is doing individually.

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) AOD.

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) AOD:
 if(aAOD)
 {
  // Remark: Keep in sync with TString secc[gQAEventCutCounter] = ...
 
  // Basic protection:
  if(!fQAEventCutCounter){return;}

  // Total number of events:
  fQAEventCutCounter->Fill(0.5);

  // Selected tracks for Q-vectors, after all event and particle cuts:
  if(fUseSelectedTracksCuts)
  { 
   if(fSelectedTracks < fSelectedTracksCuts[0]) fQAEventCutCounter->Fill(1.5);
   if(fSelectedTracks > fSelectedTracksCuts[1]) fQAEventCutCounter->Fill(2.5);
  }

  // Centrality:
  if(fUseCentralityCuts)
  { 
   if(fCentrality < fCentralityCuts[0]) fQAEventCutCounter->Fill(3.5);
   if(fCentrality > fCentralityCuts[1]) fQAEventCutCounter->Fill(4.5);
  }

  // Vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) fQAEventCutCounter->Fill(5.5);
  if(fUseNContributorsCuts)
  {
   if((Int_t)avtx->GetNContributors()<fNContributorsCuts[0]) fQAEventCutCounter->Fill(6.5);
   if((Int_t)avtx->GetNContributors()>fNContributorsCuts[1]) fQAEventCutCounter->Fill(7.5);
  }

  if(fUseMinVertexDistanceCut)
  {
   if(sqrt(pow(avtx->GetX(),2.) + pow(avtx->GetY(),2.) + pow(avtx->GetY(),2.)) < fMinVertexDistance)
   { 
    if(fQAAnomalousEvents){fQAAnomalousEvents->Fill(0.5);} // |vertex| = 0.
    fQAEventCutCounter->Fill(8.5);
   }
  }

  if(fUseVertexCuts[X])
  {
   if(avtx->GetX() < fVertexCuts[X][0]) fQAEventCutCounter->Fill(9.5);
   if(avtx->GetX() > fVertexCuts[X][1]) fQAEventCutCounter->Fill(10.5); 
  }
  if(fUseVertexCuts[Y])
  {
   if(avtx->GetY() < fVertexCuts[Y][0]) fQAEventCutCounter->Fill(11.5);
   if(avtx->GetY() > fVertexCuts[Y][1]) fQAEventCutCounter->Fill(12.5);
  }
  if(fUseVertexCuts[Z])
  {
   if(avtx->GetZ() < fVertexCuts[Z][0]) fQAEventCutCounter->Fill(13.5);
   if(avtx->GetZ() > fVertexCuts[Z][1]) fQAEventCutCounter->Fill(14.5);
  }

  // Centrality cuts:
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    if(fUseCentralityCorrelationsCuts[ce1][ce2])
    {
     Double_t centrality1 = ams->GetMultiplicityPercentile(sce[ce1].Data());
     Double_t centrality2 = ams->GetMultiplicityPercentile(sce[ce2].Data());
     if(centrality1+centrality2 > 0.) 
     {
      if(TMath::Abs((centrality1-centrality2)/(centrality1+centrality2)) > fCentralityCorrelationsCuts[ce1][ce2])
      {
       // Determine programatically which bin is that:
       for(Int_t b=1; b<=fQAEventCutCounter->GetNbinsX();b++)
       {
        if(TString(fQAEventCutCounter->GetXaxis()->GetBinLabel(b)).EqualTo(Form("centCorrCut[%d][%d]",ce1,ce2)))
        {
         //cout<<Form("%d %d",ce1,ce2)<<endl;
         //cout<<Form("%d",b)<<endl;
         //cout<<fQAEventCutCounter->GetXaxis()->GetBinLabel(b)<<endl;
         fQAEventCutCounter->Fill(b-0.5);
        } // if(TString(fQAEventCutCounter->GetBinLabel(b)).EqualTo(Form("centCorrCut[%d][%d]",ce1,ce2)))
       } // for(Int_t b=1; b<=fQAEventCutCounter->GetNBinsX();b++)
      } // if(TMath::Abs((centrality1-centrality2)/(centrality1+centrality2)) > fCentralityCorrelationsCuts[ce1][ce2]) 
     } // if(centrality1+centrality2 > 0.) 
    } // if(fUseCentralityCorrelationsCuts[ce1][ce2])
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)

  // Centrality weights (flattening): 
  // Remark: since I am getting centrality weights from centrality distribution after the events cuts, flattening is applied here after all other event cuts:
  if(!ams){cout<<__LINE__<<endl;exit(1);} // pointer was obtained previously
  Double_t centralityWeight = -44.;
  for(Int_t ce=0;ce<gCentralityEstimators;ce++)
  {
   if(fUseCentralityWeights[ce])
   {
    centralityWeight = CentralityWeight(ams->GetMultiplicityPercentile(sce[ce].Data()),sce[ce].Data());
    if(gRandom->Uniform(0,1) > centralityWeight) fQAEventCutCounter->Fill(21.5);
   } 
  }

  // Generic correlations:
  if(fUseGenericCorrelationsCuts[0])
  {
   if(avtx)
   {
    if(fSelectedTracks > fSelContrTreshold * (Int_t)avtx->GetNContributors()) fQAEventCutCounter->Fill(22.5);
   }
  }

 } // if(aAOD)

} // void AliAnalysisTaskMuPa::EventCutCounter(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMuPa::SequentialEventCutCounter(AliVEvent *ave)
{
 // Simple QA function to summarize what each event cut is doing individually.

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) AOD.

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) AOD:
 if(aAOD)
 {
  // Remark: Keep in sync with TString secc[gQASequentialEventCutCounter] = ...
 
  // Basic protection:
  if(!fQASequentialEventCutCounter){return;}

  // Total number of events:
  fQASequentialEventCutCounter->Fill(0.5);

  // Selected tracks for Q-vectors, after all event and particle cuts:
  if(fUseSelectedTracksCuts)
  { 
   if(fSelectedTracks < fSelectedTracksCuts[0]) { fQASequentialEventCutCounter->Fill(1.5); return; }
   if(fSelectedTracks > fSelectedTracksCuts[1]) { fQASequentialEventCutCounter->Fill(2.5); return; }
  }

  // Centrality:
  if(fUseCentralityCuts)
  { 
   if(fCentrality < fCentralityCuts[0]) { fQASequentialEventCutCounter->Fill(3.5); return; }
   if(fCentrality > fCentralityCuts[1]) { fQASequentialEventCutCounter->Fill(4.5); return; }
  }

  // Vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(!avtx) { fQASequentialEventCutCounter->Fill(5.5); return; }
  if(fUseNContributorsCuts)
  {
   if((Int_t)avtx->GetNContributors()<fNContributorsCuts[0]) { fQASequentialEventCutCounter->Fill(6.5); return; }
   if((Int_t)avtx->GetNContributors()>fNContributorsCuts[1]) { fQASequentialEventCutCounter->Fill(7.5); return; }
  }

  if(fUseMinVertexDistanceCut)
  {
   if(sqrt(pow(avtx->GetX(),2.) + pow(avtx->GetY(),2.) + pow(avtx->GetY(),2.)) < fMinVertexDistance)
   { 
    fQAAnomalousEvents->Fill(0.5); // |vertex| = 0.
    fQASequentialEventCutCounter->Fill(8.5); return; 
   }
  }

  if(fUseVertexCuts[X])
  {
   if(avtx->GetX() < fVertexCuts[X][0]) { fQASequentialEventCutCounter->Fill(9.5); return; }
   if(avtx->GetX() > fVertexCuts[X][1]) { fQASequentialEventCutCounter->Fill(10.5); return; }
  }
  if(fUseVertexCuts[Y])
  {
   if(avtx->GetY() < fVertexCuts[Y][0]) { fQASequentialEventCutCounter->Fill(11.5); return; }
   if(avtx->GetY() > fVertexCuts[Y][1]) { fQASequentialEventCutCounter->Fill(12.5); return; }
  }
  if(fUseVertexCuts[Z])
  {
   if(avtx->GetZ() < fVertexCuts[Z][0]) { fQASequentialEventCutCounter->Fill(13.5); return; }
   if(avtx->GetZ() > fVertexCuts[Z][1]) { fQASequentialEventCutCounter->Fill(14.5); return; }
  }

  // Centrality cuts:
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    if(fUseCentralityCorrelationsCuts[ce1][ce2])
    {
     Double_t centrality1 = ams->GetMultiplicityPercentile(sce[ce1].Data());
     Double_t centrality2 = ams->GetMultiplicityPercentile(sce[ce2].Data());
     if(centrality1+centrality2 > 0.) 
     {
      if(TMath::Abs((centrality1-centrality2)/(centrality1+centrality2)) > fCentralityCorrelationsCuts[ce1][ce2])
      {
       // Determine programatically which bin is that:
       for(Int_t b=1; b<=fQASequentialEventCutCounter->GetNbinsX();b++)
       {
        if(TString(fQASequentialEventCutCounter->GetXaxis()->GetBinLabel(b)).EqualTo(Form("centCorrCut[%d][%d]",ce1,ce2)))
        {
         //cout<<Form("%d %d",ce1,ce2)<<endl;
         //cout<<Form("%d",b)<<endl;
         //cout<<fQASequentialEventCutCounter->GetXaxis()->GetBinLabel(b)<<endl;
         fQASequentialEventCutCounter->Fill(b-0.5); return; 
        } // if(TString(fQASequentialEventCutCounter->GetBinLabel(b)).EqualTo(Form("centCorrCut[%d][%d]",ce1,ce2)))
       } // for(Int_t b=1; b<=fQASequentialEventCutCounter->GetNBinsX();b++)
      } // if(TMath::Abs((centrality1-centrality2)/(centrality1+centrality2)) > fCentralityCorrelationsCuts[ce1][ce2]) 
     } // if(centrality1+centrality2 > 0.) 
    } // if(fUseCentralityCorrelationsCuts[ce1][ce2])
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)

  // Centrality weights (flattening): 
  // Remark: since I am getting centrality weights from centrality distribution after the events cuts, flattening is applied here after all other event cuts:
  if(!ams){cout<<__LINE__<<endl;exit(1);} // pointer was obtained previously
  Double_t centralityWeight = -44.;
  for(Int_t ce=0;ce<gCentralityEstimators;ce++)
  {
   if(fUseCentralityWeights[ce])
   {
    centralityWeight = CentralityWeight(ams->GetMultiplicityPercentile(sce[ce].Data()),sce[ce].Data());
    if(gRandom->Uniform(0,1) > centralityWeight) { fQASequentialEventCutCounter->Fill(21.5); return; }
   } 
  }

  // Generic correlations:
  if(fUseGenericCorrelationsCuts[0])
  {
   if(avtx)
   {
    if(fSelectedTracks > fSelContrTreshold * (Int_t)avtx->GetNContributors()) { fQASequentialEventCutCounter->Fill(22.5); return; } 
   }
  }

  } // if(aAOD)

} // void AliAnalysisTaskMuPa::SequentialEventCutCounter(AliVEvent *ave)

//=======================================================================================================================

Bool_t AliAnalysisTaskMuPa::SurvivesParticleCuts(AliVParticle *vParticle)
{
 // Check if the current partice survives the specific track cuts (e.g. applied only on TPC-only tracks).

 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track;
 // b) Cut on AOD track;
 // c) Cut on AOD MC particle.

 //if(fVerbose){Green(__PRETTY_FUNCTION__);}

 if(fFillQAHistograms)
 {
  if(vParticle){ParticleCutCounter(vParticle);}
 }

 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track:
 AliAODMCParticle *aodmcParticle = dynamic_cast<AliAODMCParticle*>(vParticle);
 AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(vParticle);
 
 // b) Cut on AOD track:
 if(aodTrack)
 {
  // Filter bit:
  if(!aodTrack->TestFilterBit(fFilterBit)) return kFALSE;

  // Trivial cuts:
  if(fUseOnlyPrimaries)
  {
   if(aodTrack->GetType() != AliAODTrack::kPrimary) return kFALSE; // take into account only primaries
  }

  // Kinematics:
  if(fUseKinematicsCuts[PHI])
  {
   if(aodTrack->Phi() < fKinematicsCuts[PHI][0]) return kFALSE;
   if(aodTrack->Phi() >= fKinematicsCuts[PHI][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[PT])
  {
   if(aodTrack->Pt() < fKinematicsCuts[PT][0]) return kFALSE;
   if(aodTrack->Pt() >= fKinematicsCuts[PT][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[ETA])
  {
   if(aodTrack->Eta() < fKinematicsCuts[ETA][0]) return kFALSE;
   if(aodTrack->Eta() >= fKinematicsCuts[ETA][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[E])
  {
   if(aodTrack->E() < fKinematicsCuts[E][0]) return kFALSE;
   if(aodTrack->E() >= fKinematicsCuts[E][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[CHARGE])
  {
   if(0 == aodTrack->Charge()) return kFALSE; // take into account only charged particles 
   if(aodTrack->Charge() < fKinematicsCuts[CHARGE][0]) return kFALSE;
   if(aodTrack->Charge() >= fKinematicsCuts[CHARGE][1]) return kFALSE;
  }

  // DCA:
  if(fUseDCACuts[0]||fUseDCACuts[1])
  {
   Float_t dcaXY = aodTrack->DCA(); 
   Float_t dcaZ = aodTrack->ZAtDCA(); 
   if(!aodTrack->TestBit(AliAODTrack::kIsDCA))
   {
    aodTrack->GetImpactParameters(dcaXY,dcaZ); 
   } 
   // DCA xy:
   if(fUseDCACuts[0])
   {
    if(dcaXY < fDCACuts[0][0]) return kFALSE;
    if(dcaXY >= fDCACuts[0][1]) return kFALSE;
   } 
   // DCA z:
   if(fUseDCACuts[1])
   {
    if(dcaZ < fDCACuts[1][0]) return kFALSE;
    if(dcaZ >= fDCACuts[1][1]) return kFALSE;
   }
  } // if(fUseDCACuts[0]||fUseDCACuts[1])

  // Remaining particle cuts:
  if(fUseParticleCuts[TPCNcls])
  {
   if(aodTrack->GetTPCNcls() < fParticleCuts[TPCNcls][0]) return kFALSE;
   if(aodTrack->GetTPCNcls() >= fParticleCuts[TPCNcls][1]) return kFALSE;
  }  
  if(fUseParticleCuts[TPCnclsS])
  {
   if(aodTrack->GetTPCnclsS() < fParticleCuts[TPCnclsS][0]) return kFALSE;
   if(aodTrack->GetTPCnclsS() >= fParticleCuts[TPCnclsS][1]) return kFALSE;
  }
  if(fUseParticleCuts[TPCnclsFractionShared])
  {
   if(TMath::Abs((Double_t)aodTrack->GetTPCNcls())>0.)
   {
    if(!((Double_t)aodTrack->GetTPCnclsS()/(Double_t)aodTrack->GetTPCNcls() < fParticleCuts[TPCnclsFractionShared][0])) return kFALSE;
    // if((Double_t)aodTrack->GetTPCnclsS()/(Double_t)aodTrack->GetTPCNcls() >= fParticleCuts[TPCnclsFractionShared][1]) return kFALSE; // the upper limit is not needed
   }
  }
  if(fUseParticleCuts[TPCNCrossedRows])
  {  
   if(aodTrack->GetTPCNCrossedRows() < fParticleCuts[TPCNCrossedRows][0]) return kFALSE;
   if(aodTrack->GetTPCNCrossedRows() >= fParticleCuts[TPCNCrossedRows][1]) return kFALSE;
  }
  if(fUseParticleCuts[TPCChi2perNDF])
  {
   if(aodTrack->GetTPCNcls()>0.)
   {
    if(!((Double_t)aodTrack->GetTPCchi2()/(Double_t)aodTrack->GetTPCNcls() < fParticleCuts[TPCChi2perNDF][0])) return kFALSE;
    //if((Double_t)aodTrack->GetTPCchi2()/(Double_t)aodTrack->GetTPCNcls() >= fParticleCuts[TPCChi2perNDF][1]) return kFALSE; // open for the time being
   }
  } 
  if(fUseParticleCuts[TPCFoundFraction])
  {
   if(aodTrack->GetTPCFoundFraction() < fParticleCuts[TPCFoundFraction][0]) return kFALSE;
   if(aodTrack->GetTPCFoundFraction() >= fParticleCuts[TPCFoundFraction][1]) return kFALSE;
  }
  if(fUseParticleCuts[Chi2TPCConstrainedVsGlobal])
  {
   if(aodTrack->GetChi2TPCConstrainedVsGlobal() < fParticleCuts[Chi2TPCConstrainedVsGlobal][0]) return kFALSE;
   if(aodTrack->GetChi2TPCConstrainedVsGlobal() >= fParticleCuts[Chi2TPCConstrainedVsGlobal][1]) return kFALSE;
  }
  if(fUseParticleCuts[ITSNcls])
  {
   if(aodTrack->GetITSNcls() < fParticleCuts[ITSNcls][0]) return kFALSE;
   if(aodTrack->GetITSNcls() >= fParticleCuts[ITSNcls][1]) return kFALSE;
  }
  if(fUseParticleCuts[ITSChi2perNDF])
  {
   if(TMath::Abs(aodTrack->GetITSNcls())>0.)
   {
    if(aodTrack->GetITSchi2()/aodTrack->GetITSNcls() < fParticleCuts[ITSChi2perNDF][0]) return kFALSE;
    if(aodTrack->GetITSchi2()/aodTrack->GetITSNcls() >= fParticleCuts[ITSChi2perNDF][1]) return kFALSE;
   }
  }
  if(fUseParticleCuts[TPCNclsF])
  {
   if(aodTrack->GetTPCNclsF() < fParticleCuts[TPCNclsF][0]) return kFALSE;
   if(aodTrack->GetTPCNclsF() >= fParticleCuts[TPCNclsF][1]) return kFALSE;
  }
  if(fAtLeastOnePointInTheSPD)
  {
   if(!(aodTrack->HasPointOnITSLayer(1) || aodTrack->HasPointOnITSLayer(2))) return kFALSE;
  }
  if(fIgnoreGlobalConstrained) // kTRUE by default
  {
   if(aodTrack->IsGlobalConstrained()) return kFALSE;
  }

  // Toy NUA:
  if(fUseToyNUA[PHI] && aodTrack->Phi() >= fToyNUACuts[PHI][1] && aodTrack->Phi() < fToyNUACuts[PHI][2])
  {
   if(gRandom->Uniform(0,1) > fToyNUACuts[PHI][0]) return kFALSE; 
  }
  if(fUseToyNUA[PT] && aodTrack->Pt() >= fToyNUACuts[PT][1] && aodTrack->Pt() < fToyNUACuts[PT][2])
  {
   if(gRandom->Uniform(0,1) > fToyNUACuts[PT][0]) return kFALSE; 
  }
  if(fUseToyNUA[ETA] && aodTrack->Eta() >= fToyNUACuts[ETA][1] && aodTrack->Eta() < fToyNUACuts[ETA][2])
  {
   if(gRandom->Uniform(0,1) > fToyNUACuts[ETA][0]) return kFALSE; 
  }
  if(fUseToyNUA[E] && aodTrack->E() >= fToyNUACuts[E][1] && aodTrack->E() < fToyNUACuts[E][2])
  {
   if(gRandom->Uniform(0,1) > fToyNUACuts[E][0]) return kFALSE; 
  }
  if(fUseToyNUA[CHARGE] && aodTrack->Charge() >= fToyNUACuts[CHARGE][1] && aodTrack->Charge() < fToyNUACuts[CHARGE][2])
  {
   if(gRandom->Uniform(0,1) > fToyNUACuts[CHARGE][0]) return kFALSE; 
  }
 } // if(aodTrack)

 // c) Cut on AOD MC particle:
 if(aodmcParticle)
 {
  // Trivial cuts:
  if(fUseOnlyPrimaries)
  {
   if(!aodmcParticle->IsPrimary()) return kFALSE; // take into account only generated primaries 
   //if(!aodmcParticle->IsPhysicalPrimary()) return kFALSE; // take into account only what ALICE defines as primaries
  }

  // Kinematics:
  if(fUseKinematicsCuts[PHI])
  {
   if(aodmcParticle->Phi() < fKinematicsCuts[PHI][0]) return kFALSE;
   if(aodmcParticle->Phi() >= fKinematicsCuts[PHI][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[PT])
  {
   if(aodmcParticle->Pt() < fKinematicsCuts[PT][0]) return kFALSE;
   if(aodmcParticle->Pt() >= fKinematicsCuts[PT][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[ETA])
  {
   if(aodmcParticle->Eta() < fKinematicsCuts[ETA][0]) return kFALSE;
   if(aodmcParticle->Eta() >= fKinematicsCuts[ETA][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[E])
  {
   if(aodmcParticle->E() < fKinematicsCuts[E][0]) return kFALSE;
   if(aodmcParticle->E() >= fKinematicsCuts[E][1]) return kFALSE;
  }
  if(fUseKinematicsCuts[CHARGE])
  {
   if(0 == aodmcParticle->Charge()) return kFALSE;
   if(aodmcParticle->Charge()/3 < fKinematicsCuts[CHARGE][0]) return kFALSE; // yes, since aodmcParticle->Charge() = 3 if aodTrack->Charge() = 3
   if(aodmcParticle->Charge()/3 >= fKinematicsCuts[CHARGE][1]) return kFALSE;
  }

 } // if(aodmcParticle)
 
 return kTRUE;

} // Bool_t AliAnalysisTaskMuPa::SurvivesParticleCuts(AliVParticle *vParticle)

//=======================================================================================================================

void AliAnalysisTaskMuPa::ParticleCutCounter(AliVParticle *vParticle)
{
 // Used only in QA. Counter for each particle cut separately how many particles are rejected by the cut.

 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track;
 // b) AOD track;
 // c) AOD MC particle.

 // a) Determine Ali{AODMC,ESD,AOD}Particle/Track:
 AliAODMCParticle *aodmcParticle = dynamic_cast<AliAODMCParticle*>(vParticle);
 AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(vParticle);
 
 // b) Cut on AOD track:
 if(aodTrack)
 {
  // Basic protection:
  if(!fQAParticleCutCounter[RECO]){return;}

  // Filter bit:
  if(!aodTrack->TestFilterBit(fFilterBit)) fQAParticleCutCounter[RECO]->Fill(0.5);
 
  // Trivial cuts:
  if(fUseOnlyPrimaries)
  {
   if(aodTrack->GetType() != AliAODTrack::kPrimary) fQAParticleCutCounter[RECO]->Fill(1.5); // take into account only primaries
  }

  // Kinematics:
  if(fUseKinematicsCuts[PHI])
  {
   if(aodTrack->Phi() < fKinematicsCuts[PHI][0]) fQAParticleCutCounter[RECO]->Fill(2.5);
   if(aodTrack->Phi() >= fKinematicsCuts[PHI][1]) fQAParticleCutCounter[RECO]->Fill(3.5);
  }
  if(fUseKinematicsCuts[PT])
  {
   if(aodTrack->Pt() < fKinematicsCuts[PT][0]) fQAParticleCutCounter[RECO]->Fill(4.5);
   if(aodTrack->Pt() >= fKinematicsCuts[PT][1]) fQAParticleCutCounter[RECO]->Fill(5.5);
  }
  if(fUseKinematicsCuts[ETA])
  {
   if(aodTrack->Eta() < fKinematicsCuts[ETA][0]) fQAParticleCutCounter[RECO]->Fill(6.5);
   if(aodTrack->Eta() >= fKinematicsCuts[ETA][1]) fQAParticleCutCounter[RECO]->Fill(7.5);
  }
  if(fUseKinematicsCuts[E])
  {
   if(aodTrack->E() < fKinematicsCuts[E][0]) fQAParticleCutCounter[RECO]->Fill(8.5);
   if(aodTrack->E() >= fKinematicsCuts[E][1]) fQAParticleCutCounter[RECO]->Fill(9.5);
  }
  if(fUseKinematicsCuts[CHARGE])
  {
   if(0 == aodTrack->Charge()) fQAParticleCutCounter[RECO]->Fill(10.5); // take into account only charged particles 
   if(aodTrack->Charge() < fKinematicsCuts[CHARGE][0]) fQAParticleCutCounter[RECO]->Fill(11.5);
   if(aodTrack->Charge() >= fKinematicsCuts[CHARGE][1]) fQAParticleCutCounter[RECO]->Fill(12.5);
  }

  // DCA:
  if(fUseDCACuts[0]||fUseDCACuts[1])
  {
   Float_t dcaXY = aodTrack->DCA(); 
   Float_t dcaZ = aodTrack->ZAtDCA(); 
   if(!aodTrack->TestBit(AliAODTrack::kIsDCA))
   {
    aodTrack->GetImpactParameters(dcaXY,dcaZ); 
   } 
   // DCA xy:
   if(fUseDCACuts[0])
   {
    if(dcaXY < fDCACuts[0][0]) fQAParticleCutCounter[RECO]->Fill(13.5);
    if(dcaXY >= fDCACuts[0][1]) fQAParticleCutCounter[RECO]->Fill(14.5);
   } 
   // DCA z:
   if(fUseDCACuts[1])
   {
    if(dcaZ < fDCACuts[1][0]) fQAParticleCutCounter[RECO]->Fill(15.5);
    if(dcaZ >= fDCACuts[1][1]) fQAParticleCutCounter[RECO]->Fill(16.5);
   }
  } // if(fUseDCACuts[0]||fUseDCACuts[1])

  // Remaining particle cuts:
  if(fUseParticleCuts[TPCNcls])
  {
   if(aodTrack->GetTPCNcls() < fParticleCuts[TPCNcls][0]) fQAParticleCutCounter[RECO]->Fill(17.5);
   if(aodTrack->GetTPCNcls() >= fParticleCuts[TPCNcls][1]) fQAParticleCutCounter[RECO]->Fill(18.5);
  }
  if(fUseParticleCuts[TPCnclsS])
  {
   if(aodTrack->GetTPCnclsS() < fParticleCuts[TPCnclsS][0]) fQAParticleCutCounter[RECO]->Fill(19.5);
   if(aodTrack->GetTPCnclsS() >= fParticleCuts[TPCnclsS][1]) fQAParticleCutCounter[RECO]->Fill(20.5);
  }
  if(fUseParticleCuts[TPCnclsFractionShared])
  {
   if(TMath::Abs((Double_t)aodTrack->GetTPCNcls())>0.)
   {
    if(!((Double_t)aodTrack->GetTPCnclsS()/(Double_t)aodTrack->GetTPCNcls() < fParticleCuts[TPCnclsFractionShared][0])) fQAParticleCutCounter[RECO]->Fill(21.5);
    // if((Double_t)aodTrack->GetTPCnclsS()/(Double_t)aodTrack->GetTPCNcls() >= fParticleCuts[TPCnclsFractionShared][1]) fQAParticleCutCounter[RECO]->Fill(22.5); // the upper limit is not needed
   }
  }
  if(fUseParticleCuts[TPCNCrossedRows])
  {
   if(aodTrack->GetTPCNCrossedRows() < fParticleCuts[TPCNCrossedRows][0]) fQAParticleCutCounter[RECO]->Fill(23.5);
   if(aodTrack->GetTPCNCrossedRows() >= fParticleCuts[TPCNCrossedRows][1]) fQAParticleCutCounter[RECO]->Fill(24.5);
  }
  if(fUseParticleCuts[TPCChi2perNDF])
  {
   if(aodTrack->GetTPCNcls()>0.)
   {
    if(!((Double_t)aodTrack->GetTPCchi2()/(Double_t)aodTrack->GetTPCNcls() < fParticleCuts[TPCChi2perNDF][0])) fQAParticleCutCounter[RECO]->Fill(25.5);
    // if((Double_t)aodTrack->GetTPCchi2()/(Double_t)aodTrack->GetTPCNcls() >= fParticleCuts[TPCChi2perNDF][1]) fQAParticleCutCounter[RECO]->Fill(26.5); // open for the time being
   }
  } 
  if(fUseParticleCuts[TPCFoundFraction])
  {
   if(aodTrack->GetTPCFoundFraction() < fParticleCuts[TPCFoundFraction][0]) fQAParticleCutCounter[RECO]->Fill(27.5);
   if(aodTrack->GetTPCFoundFraction() >= fParticleCuts[TPCFoundFraction][1]) fQAParticleCutCounter[RECO]->Fill(28.5);
  }
  if(fUseParticleCuts[Chi2TPCConstrainedVsGlobal])
  {
   if(aodTrack->GetChi2TPCConstrainedVsGlobal() < fParticleCuts[Chi2TPCConstrainedVsGlobal][0]) fQAParticleCutCounter[RECO]->Fill(29.5);
   if(aodTrack->GetChi2TPCConstrainedVsGlobal() >= fParticleCuts[Chi2TPCConstrainedVsGlobal][1]) fQAParticleCutCounter[RECO]->Fill(30.5);
  }
  if(fUseParticleCuts[ITSNcls])
  {
   if(aodTrack->GetITSNcls() < fParticleCuts[ITSNcls][0]) fQAParticleCutCounter[RECO]->Fill(31.5);
   if(aodTrack->GetITSNcls() >= fParticleCuts[ITSNcls][1]) fQAParticleCutCounter[RECO]->Fill(32.5);
  }
  if(fUseParticleCuts[ITSChi2perNDF])
  {
   if(TMath::Abs(aodTrack->GetITSNcls())>0.)
   {
    if(aodTrack->GetITSchi2()/aodTrack->GetITSNcls() < fParticleCuts[ITSChi2perNDF][0]) fQAParticleCutCounter[RECO]->Fill(33.5);
    if(aodTrack->GetITSchi2()/aodTrack->GetITSNcls() >= fParticleCuts[ITSChi2perNDF][1]) fQAParticleCutCounter[RECO]->Fill(34.5);
   }
  }
  if(fUseParticleCuts[TPCNclsF])
  {
   if(aodTrack->GetTPCNclsF() < fParticleCuts[TPCNclsF][0]) fQAParticleCutCounter[RECO]->Fill(35.5);
   if(aodTrack->GetTPCNclsF() >= fParticleCuts[TPCNclsF][1]) fQAParticleCutCounter[RECO]->Fill(36.5);
  }
  if(fAtLeastOnePointInTheSPD)
  {
   if(!(aodTrack->HasPointOnITSLayer(1) || aodTrack->HasPointOnITSLayer(2))) fQAParticleCutCounter[RECO]->Fill(37.5);
  }
  if(fIgnoreGlobalConstrained) // kTRUE by default
  {
   if(aodTrack->IsGlobalConstrained()) fQAParticleCutCounter[RECO]->Fill(38.5);
  }
 } // if(aodTrack)

 // c) Cut on AOD MC particle:
 if(aodmcParticle)
 {
  // Basic protection:
  if(!fQAParticleCutCounter[SIM]){return;}

  // Trivial cuts:
  if(fUseOnlyPrimaries)
  {
   if(!aodmcParticle->IsPrimary())  fQAParticleCutCounter[SIM]->Fill(1.5); // take into account only generated primaries 
   //if(!aodmcParticle->IsPhysicalPrimary()) return kFALSE; // take into account only what ALICE defines as primaries
  }

  // Kinematics:
  if(fUseKinematicsCuts[PHI])
  {
   if(aodmcParticle->Phi() < fKinematicsCuts[PHI][0]) fQAParticleCutCounter[SIM]->Fill(2.5);
   if(aodmcParticle->Phi() >= fKinematicsCuts[PHI][1]) fQAParticleCutCounter[SIM]->Fill(3.5);
  }
  if(fUseKinematicsCuts[PT])
  {
   if(aodmcParticle->Pt() < fKinematicsCuts[PT][0]) fQAParticleCutCounter[SIM]->Fill(4.5);
   if(aodmcParticle->Pt() >= fKinematicsCuts[PT][1]) fQAParticleCutCounter[SIM]->Fill(5.5);
  }
  if(fUseKinematicsCuts[ETA])
  {
   if(aodmcParticle->Eta() < fKinematicsCuts[ETA][0]) fQAParticleCutCounter[SIM]->Fill(6.5);
   if(aodmcParticle->Eta() >= fKinematicsCuts[ETA][1]) fQAParticleCutCounter[SIM]->Fill(7.5);
  }
  if(fUseKinematicsCuts[E])
  {
   if(aodmcParticle->E() < fKinematicsCuts[E][0]) fQAParticleCutCounter[SIM]->Fill(8.5);
   if(aodmcParticle->E() >= fKinematicsCuts[E][1]) fQAParticleCutCounter[SIM]->Fill(9.5);
  }
  if(fUseKinematicsCuts[CHARGE])
  {
   if(0 == aodmcParticle->Charge()) fQAParticleCutCounter[SIM]->Fill(10.5);
   if(aodmcParticle->Charge()/3 < fKinematicsCuts[CHARGE][0]) fQAParticleCutCounter[SIM]->Fill(11.5); // yes, since aodmcParticle->Charge() = 3 if aodTrack->Charge() = 3
   if(aodmcParticle->Charge()/3 >= fKinematicsCuts[CHARGE][1]) fQAParticleCutCounter[SIM]->Fill(12.5);
  }

 } // if(aodmcParticle)

} // void AliAnalysisTaskMuPa::ParticleCutCounter(AliVParticle *vParticle)

//=======================================================================================================================

void AliAnalysisTaskMuPa::FillQAHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)
{
 // Fill all additional QA stuff here, that is already not filled in Control Event Histograms or Control Particle Histograms. 

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) AOD QA;
 // c) MC QA.

 //if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);


 // b) AOD QA:
 if(aAOD)
 {
  // Centrality: 
  AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
  if(!ams){cout<<__LINE__<<endl;exit(1);}
  TString sce[gCentralityEstimators] = {"V0M", "SPDTracklets", "CL0", "CL1"}; // keep this in sync with enum eCentralityEstimator in .h
  for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
  {
   if(fQACentralityHist[ce1][ba]){fQACentralityHist[ce1][ba]->Fill(ams->GetMultiplicityPercentile(sce[ce1].Data()));}
   for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++)
   {
    if(fQACentralityCorrHist[ce1][ce2][ba]){fQACentralityCorrHist[ce1][ce2][ba]->Fill(ams->GetMultiplicityPercentile(sce[ce1].Data()),ams->GetMultiplicityPercentile(sce[ce2].Data()));}
   } // for(Int_t ce2=ce1+1;ce2<gCentralityEstimators;ce2++) 
  } // for(Int_t ce1=0;ce1<gCentralityEstimators;ce1++)
   

  AliAODHeader *aodheader = (AliAODHeader*)aAOD->GetHeader();
  if(!aodheader){cout<<__LINE__<<endl;exit(1);}
  if(fQAMultiplicityCorrHist[ba])
  {
   fQAMultiplicityCorrHist[ba]->Fill(aodheader->GetRefMultiplicityComb08(),fSelectedTracks);
  }

  // *) Start loop over AOD tracks:
  Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  {
   AliAODTrack *aodTrack = NULL;
   if(!fUseFisherYates)
   {
    aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to "a track" (i.e. any track)
   }
   else
   {
    aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack((Int_t)fRandomIndices->GetAt(iTrack)));
   }
   if(!aodTrack){continue;}
   
   // Filter bit distribution:
   for(Int_t fb=0;fb<gFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
   {
    if(aodTrack->TestFilterBit(1<<fb))
    {
     if(fQAFilterBitScan){fQAFilterBitScan->Fill(fb);}
     if(fQAIDvsFilterBit){fQAIDvsFilterBit->Fill(fb,aodTrack->GetID());}
    }
   }

   // Filter bit kinematics:
   if(rs == RECO)
   {
    for(Int_t fb=0;fb<fQAFilterBits->GetSize();fb++)
    {
     if(aodTrack->TestFilterBit((Int_t)fQAFilterBits->GetAt(fb)))
     {
      if(fQAKinematicsFilterBits[fb][PHI]){fQAKinematicsFilterBits[fb][PHI]->Fill(aodTrack->Phi());}
      if(fQAKinematicsFilterBits[fb][PT]){fQAKinematicsFilterBits[fb][PT]->Fill(aodTrack->Pt());}
      if(fQAKinematicsFilterBits[fb][ETA]){fQAKinematicsFilterBits[fb][ETA]->Fill(aodTrack->Eta());}
      if(fQAKinematicsFilterBits[fb][E]){fQAKinematicsFilterBits[fb][E]->Fill(aodTrack->E());}
      if(fQAKinematicsFilterBits[fb][CHARGE]){fQAKinematicsFilterBits[fb][CHARGE]->Fill(aodTrack->Charge());}
     }
    }
   } // if(rs == RECO)
   
  } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks

  // Check for self-correlations in reconstructed data with two nested loops:
  if(fQACheckSelfCorrelations) 
  {
   for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++) // starting a loop over the first track
   {
    AliAODTrack *aodTrack1 = NULL;
    if(!fUseFisherYates)
    {
     aodTrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1)); // getting a pointer to "a track" (i.e. any track)
    }
    else
    {
     aodTrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack((Int_t)fRandomIndices->GetAt(iTrack1)));
    }
    if(!aodTrack1){continue;}
    if(!SurvivesParticleCuts(aodTrack1)){continue;} 
    for(Int_t iTrack2=iTrack1+1;iTrack2<nTracks;iTrack2++) // starting a loop over the second track
    {
     AliAODTrack *aodTrack2 = NULL;
     if(!fUseFisherYates)
     {
      aodTrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2)); // getting a pointer to "a track" (i.e. any track)
     }
     else
     {
      aodTrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack((Int_t)fRandomIndices->GetAt(iTrack2)));
     }
     if(!aodTrack2){continue;}
     if(!SurvivesParticleCuts(aodTrack2)){continue;} 
     if(fQASelfCorrelations[0]){fQASelfCorrelations[0]->Fill(aodTrack1->Phi()-aodTrack2->Phi());}
     if(fQASelfCorrelations[1]){fQASelfCorrelations[1]->Fill(aodTrack1->Pt()-aodTrack2->Pt());}  
     if(fQASelfCorrelations[2]){fQASelfCorrelations[2]->Fill(aodTrack1->Eta()-aodTrack2->Eta());}
    } // for(Int_t iTrack2=iTrack1+1;iTrack2<nTracks;iTrack2++) // starting a loop over the second track
   } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++) // starting a loop over the first track
  }

  // Trigger:
  if(rs == RECO && fQATrigger[ba])
  {
   UInt_t fSelectMask = fInputHandler->IsEventSelected();
   if(fSelectMask & AliVEvent::kAny)
   {
    fQATrigger[ba]->Fill(0.5);
   }
   if(fSelectMask & AliVEvent::kMB)
   {
    fQATrigger[ba]->Fill(1.5);
   }
   if(fSelectMask & AliVEvent::kINT7)
   {
    fQATrigger[ba]->Fill(2.5);
   }
   if(fSelectMask & AliVEvent::kCentral)
   {
    fQATrigger[ba]->Fill(3.5);
   }
   if(fSelectMask & AliVEvent::kSemiCentral)
   {
    fQATrigger[ba]->Fill(4.5);
   }
   if((fSelectMask & AliVEvent::kCentral) || (fSelectMask & AliVEvent::kMB))
   {
    fQATrigger[ba]->Fill(5.5);
   }
   if((fSelectMask & AliVEvent::kSemiCentral) || (fSelectMask & AliVEvent::kMB))
   {
    fQATrigger[ba]->Fill(6.5);
   }
  } // if(rs == RECO && fQATrigger[ba])

  // Generic correlations:
  if(rs == RECO)
  {
   // 0: fSelectedTracks vs. avtx->GetNContributors()
   AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
   if(avtx)
   {
    if(fQAGenericCorrHist[0][ba]){fQAGenericCorrHist[0][ba]->Fill(fSelectedTracks,(Int_t)avtx->GetNContributors());}
   }

   // 1: ...

  } // if(rs == RECO)

 } // if(aAOD)

 // c) MC QA:
 if(aMC)
 {
  // nothing is needed here at the moment
 } // if(aMC)

} // void AliAnalysisTaskMuPa::FillQAHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs)

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

 if(fUseFisherYates){cout<<__LINE__<<endl;exit(1);} // TBI 20210810 check and validate if also here Fisher-Yates needs to be applied

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
 // Generic expression <exp[i(n1*phi1)]>.

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

TComplex AliAnalysisTaskMuPa::Three(Int_t n1, Int_t n2, Int_t n3)
{
 // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.

 TComplex three = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
                - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3); 

 return three;

} // TComplex AliAnalysisTaskMuPa::Three(Int_t n1, Int_t n2, Int_t n3)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
{
 // Generic four-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.

 TComplex four = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
               - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
               + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
               + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);

 return four;

} // TComplex AliAnalysisTaskMuPa::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)
{
 // Generic five-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.

 TComplex five = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
               - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
               + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
               + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
               + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
               - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
               - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
               + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
               + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
               - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
               + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
               - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
               - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
               + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
               + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
               + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
               + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
               - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
               + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
               + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
               + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
               + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
               - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3) 
               - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
               - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
 
 return five;

} // TComplex AliAnalysisTaskMuPa::Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)
{
 // Generic six-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.

 TComplex six = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);

 return six;

} // TComplex AliAnalysisTaskMuPa::Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)
{
 // Generic seven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7)]>.

 Int_t harmonic[7] = {n1,n2,n3,n4,n5,n6,n7};

 TComplex seven = Recursion(7,harmonic); 

 return seven;

} // end of TComplex AliAnalysisTaskMuPa::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
 // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

 Int_t harmonic[8] = {n1,n2,n3,n4,n5,n6,n7,n8};

 TComplex eight = Recursion(8,harmonic); 

 return eight;

} // end of TComplex AliAnalysisTaskMuPa::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Nine(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9)
{
 // Generic nine-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9)]>.

 Int_t harmonic[9] = {n1,n2,n3,n4,n5,n6,n7,n8,n9};

 TComplex nine = Recursion(9,harmonic); 

 return nine;

} // end of TComplex AliAnalysisTaskMuPa::Nine(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Ten(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10)
{
 // Generic ten-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10)]>.

 Int_t harmonic[10] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10};

 TComplex ten = Recursion(10,harmonic); 

 return ten;

} // end of TComplex AliAnalysisTaskMuPa::Ten(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Eleven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11)
{
 // Generic eleven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11)]>.

 Int_t harmonic[11] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11};

 TComplex eleven = Recursion(11,harmonic); 

 return eleven;

} // end of TComplex AliAnalysisTaskMuPa::Eleven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Twelve(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11, Int_t n12)
{
 // Generic twelve-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11+n12*phi12)]>.

 Int_t harmonic[12] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12};

 TComplex twelve = Recursion(12,harmonic); 

 return twelve;

} // end of TComplex AliAnalysisTaskMuPa::Twelve(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11, Int_t n12)

//=======================================================================================================================

TComplex AliAnalysisTaskMuPa::Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliAnalysisTaskMuPa::Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip) 

//=======================================================================================================================

void AliAnalysisTaskMuPa::CalculateCorrelations() 
{
 // Calculate analytically multiparticle correlations from Q-vectors. 
 // By default, only correlations for which all harmonics are the same are evaluated. 

 // TBI 20210909 optimize that 2x I do not call e.g. Two(0,0).Re(), etc.

 for(Int_t h=1;h<=fMaxHarmonic;h++) // harmonic
 {
  // 2p:
  //cout<<"   => CalculateCorrelations(void), 2p .... "<<endl;
  if(fSelectedTracks<2){return;}
  TComplex twoC = Two(h,-h)/Two(0,0).Re(); // cos
  //TComplex twoS = Two(h,-h)/Two(0,0).Im(); // sin
  Double_t wTwo = Two(0,0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(2);
   harmonics->SetAt(h,0);
   harmonics->SetAt(-h,1);
   if(TMath::Abs(twoC.Re() - this->CalculateCustomNestedLoop(harmonics))>1.e-5)
   {          
    cout<<__LINE__<<endl; exit(1);
   } 
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 2-p, harmonic %d",h)<<endl;
    //cout<<Form("   value = %f",twoC.Re())<<endl;
   }
   delete harmonics;
  } // if(fCalculateCustomNestedLoop)
  if(fUseInternalValidation&&fRescaleWithTheoreticalInput&&TMath::Abs(fInternalValidationHarmonics->GetAt(h-1))>0.){twoC/=pow(fInternalValidationHarmonics->GetAt(h-1),2.);}
  // integrated:
  if(fCorrelationsPro[0][h-1][0]){fCorrelationsPro[0][h-1][0]->Fill(0.5,twoC,wTwo);}
  // vs. multiplicity:
  if(fCorrelationsPro[0][h-1][1]){fCorrelationsPro[0][h-1][1]->Fill(fSelectedTracks+0.5,twoC,wTwo);}
  // vs. centrality:
  if(fCorrelationsPro[0][h-1][2]){fCorrelationsPro[0][h-1][2]->Fill(fCentrality,twoC,wTwo);}

  // 4p:
  //cout<<"   => CalculateCorrelations(void), 4p .... "<<endl;
  if(fSelectedTracks<4){return;}
  TComplex fourC = Four(h,h,-h,-h)/Four(0,0,0,0).Re(); // cos
  //TComplex fourS = Four(h,h,-h,-h)/Four(0,0,0,0).Im(); // sin
  Double_t wFour = Four(0,0,0,0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(4);
   harmonics->SetAt(h,0);
   harmonics->SetAt(h,1);
   harmonics->SetAt(-h,2);
   harmonics->SetAt(-h,3);
   if(TMath::Abs(fourC.Re() - this->CalculateCustomNestedLoop(harmonics))>1.e-5)
   {          
    cout<<__LINE__<<endl; exit(1);
   }
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 4-p, harmonic %d",h)<<endl;
    //cout<<Form("   value = %f",fourC.Re())<<endl;
   }
   delete harmonics;
  } // if(fCalculateCustomNestedLoop)
  if(fUseInternalValidation&&fRescaleWithTheoreticalInput&&TMath::Abs(fInternalValidationHarmonics->GetAt(h-1))>0.){fourC/=pow(fInternalValidationHarmonics->GetAt(h-1),4.);}
  // integrated:
  if(fCorrelationsPro[1][h-1][0]){fCorrelationsPro[1][h-1][0]->Fill(0.5,fourC,wFour);}
  // vs. multiplicity:
  if(fCorrelationsPro[1][h-1][1]){fCorrelationsPro[1][h-1][1]->Fill(fSelectedTracks+0.5,fourC,wFour);}
  // vs. centrality:
  if(fCorrelationsPro[1][h-1][2]){fCorrelationsPro[1][h-1][2]->Fill(fCentrality,fourC,wFour);}

  // 6p:
  //cout<<"   => CalculateCorrelations(void), 6p .... "<<endl;
  if(fSelectedTracks<6){return;}
  TComplex sixC = Six(h,h,h,-h,-h,-h)/Six(0,0,0,0,0,0).Re(); // cos
  //TComplex sixS = Six(h,h,h,-h,-h,-h)/Six(0,0,0,0,0,0).Im(); // sin
  Double_t wSix = Six(0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(6);
   harmonics->SetAt(h,0);
   harmonics->SetAt(h,1);
   harmonics->SetAt(h,2);
   harmonics->SetAt(-h,3);
   harmonics->SetAt(-h,4);
   harmonics->SetAt(-h,5);
   if(TMath::Abs(sixC.Re() - this->CalculateCustomNestedLoop(harmonics))>1.e-5)
   {          
    cout<<__LINE__<<endl; exit(1);
   }
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 6-p, harmonic %d",h)<<endl;
   // cout<<Form("   value = %f",sixC.Re())<<endl;
   }
   delete harmonics;
  } // if(fCalculateCustomNestedLoop)
  if(fUseInternalValidation&&fRescaleWithTheoreticalInput&&TMath::Abs(fInternalValidationHarmonics->GetAt(h-1))>0.){sixC/=pow(fInternalValidationHarmonics->GetAt(h-1),6.);}
  // integrated:
  if(fCorrelationsPro[2][h-1][0]){fCorrelationsPro[2][h-1][0]->Fill(0.5,sixC,wSix);}
  // vs. multiplicity:
  if(fCorrelationsPro[2][h-1][1]){fCorrelationsPro[2][h-1][1]->Fill(fSelectedTracks+0.5,sixC,wSix);}
  // vs. centrality:
  if(fCorrelationsPro[2][h-1][2]){fCorrelationsPro[2][h-1][2]->Fill(fCentrality,sixC,wSix);}

  // 8p:
  //cout<<"   => CalculateCorrelations(void), 8p .... "<<endl;
  if(fSelectedTracks<8){return;}
  TComplex eightC = Eight(h,h,h,h,-h,-h,-h,-h)/Eight(0,0,0,0,0,0,0,0).Re(); // cos
  //TComplex eightS = Eight(h,h,h,h,-h,-h,-h,-h)/Eight(0,0,0,0,0,0).Im(); // sin
  Double_t wEight = Eight(0,0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(8);
   harmonics->SetAt(h,0);
   harmonics->SetAt(h,1);
   harmonics->SetAt(h,2);
   harmonics->SetAt(h,3);
   harmonics->SetAt(-h,4);
   harmonics->SetAt(-h,5);
   harmonics->SetAt(-h,6);
   harmonics->SetAt(-h,7);
   if(TMath::Abs(eightC.Re() - this->CalculateCustomNestedLoop(harmonics))>1.e-5)
   {          
    cout<<__LINE__<<endl; exit(1);
   }
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 8-p, harmonic %d",h)<<endl;
    //cout<<Form("   value = %f",eightC.Re())<<endl;
   }
   delete harmonics;
  } // if(fCalculateCustomNestedLoop)
  if(fUseInternalValidation&&fRescaleWithTheoreticalInput&&TMath::Abs(fInternalValidationHarmonics->GetAt(h-1))>0.){eightC/=pow(fInternalValidationHarmonics->GetAt(h-1),8.);}
  // integrated:
  if(fCorrelationsPro[3][h-1][0]){fCorrelationsPro[3][h-1][0]->Fill(0.5,eightC,wEight);}
  // vs. multiplicity:
  if(fCorrelationsPro[3][h-1][1]){fCorrelationsPro[3][h-1][1]->Fill(fSelectedTracks+0.5,eightC,wEight);}
  // vs. centrality:
  if(fCorrelationsPro[3][h-1][2]){fCorrelationsPro[3][h-1][2]->Fill(fCentrality,eightC,wEight);}

 } // for(Int_t h=1;h<=fMaxHarmonic;h++) // harmonic

} // void AliAnalysisTaskMuPa::CalculateCorrelations() 

//=======================================================================================================================

void AliAnalysisTaskMuPa::CalculateNestedLoops() 
{
 // Calculate correlations with nested loops. 

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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
    if(fNestedLoopsPro[0][h][0]){fNestedLoopsPro[0][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);}
    // fill cos, 2p, vs. M: 
    if(fNestedLoopsPro[0][h][1]){fNestedLoopsPro[0][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);}
    // fill cos, 2p, vs. centrality: 
    if(fNestedLoopsPro[0][h][2]){fNestedLoopsPro[0][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1-dPhi2)),dW1*dW2);}
   } // for(int h=1; h<=6; h++)
  } // for(int i2=0; i2<nTracks; ++i2)
 } // for(int i1=0; i1<nTracks; ++i1)

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
      if(fNestedLoopsPro[1][h][0]){fNestedLoopsPro[1][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);}
      // fill cos, 4p, all harmonics, vs. M: 
      if(fNestedLoopsPro[1][h][1]){fNestedLoopsPro[1][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);}
      // fill cos, 4p, all harmonics, vs. centrality: 
      if(fNestedLoopsPro[1][h][1]){fNestedLoopsPro[1][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2-dPhi3-dPhi4)),dW1*dW2*dW3*dW4);}
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
        if(fNestedLoopsPro[2][h][0]){fNestedLoopsPro[2][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);}
        // fill cos, 6p, all harmonics, vs. M: 
        if(fNestedLoopsPro[2][h][1]){fNestedLoopsPro[2][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);}
        // fill cos, 6p, all harmonics, vs. M: 
        if(fNestedLoopsPro[2][h][1]){fNestedLoopsPro[2][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3-dPhi4-dPhi5-dPhi6)),dW1*dW2*dW3*dW4*dW5*dW6);}
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
          if(fNestedLoopsPro[3][h][0]){fNestedLoopsPro[3][h][0]->Fill(0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);}
          // fill cos, 8p, all harmonics, vs. M: 
          if(fNestedLoopsPro[3][h][1]){fNestedLoopsPro[3][h][1]->Fill(fSelectedTracks+0.5,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);}
          // fill cos, 8p, all harmonics, vs. M: 
          if(fNestedLoopsPro[3][h][2]){fNestedLoopsPro[3][h][2]->Fill(fCentrality,TMath::Cos((h+1.)*(dPhi1+dPhi2+dPhi3+dPhi4-dPhi5-dPhi6-dPhi7-dPhi8)),dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8);}
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

Double_t AliAnalysisTaskMuPa::CalculateCustomNestedLoop(TArrayI *harmonics)
{
 // For the specified harmonics, get the correlation from nested loops.
 // Order of correlator is the number of harmonics, i.e. the number of elements in an array.

 // a) Determine the order of correlator;
 // b) Custom nested loop;
 // c) Return value. 

 if(!harmonics){cout<<__LINE__<<endl;exit(1);}

 // a) Determine the order of correlator;
 Int_t order = harmonics->GetSize();
 if(0==order||order>fMaxCorrelator){cout<<__LINE__<<endl;exit(1);}

 // b) Custom nested loop:
 TProfile *profile = new TProfile("profile","",1,0.,1.); // helper profile to get all averages automatically
 //profile->Sumw2();
 Double_t value = 0.; // cos of current multiplet
 Double_t weight = 1.; // weight of current multiplet
 for(int i1=0; i1<fSelectedTracks; i1++)
 {
  //if(i1==fSelectedTracks-1){cout<<"done with loop #1"<<endl;}
  Double_t dPhi1 = ftaNestedLoops[0]->GetAt(i1);
  Double_t dW1 = ftaNestedLoops[1]->GetAt(i1);
  if(1==order)
  {
   value = TMath::Cos(harmonics->GetAt(0)*dPhi1);   
   weight = dW1;
   profile->Fill(0.5,value,weight);
   continue;
  }
  for(int i2=0; i2<fSelectedTracks; i2++)
  {
   //if(i2==fSelectedTracks-1){cout<<"done with loop #2"<<endl;}
   if(i2==i1){continue;}
   Double_t dPhi2 = ftaNestedLoops[0]->GetAt(i2);
   Double_t dW2 = ftaNestedLoops[1]->GetAt(i2);
   if(2==order)
   {
    value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2);   
    weight = dW1*dW2;
    profile->Fill(0.5,value,weight);
    continue;
   }
   for(int i3=0; i3<fSelectedTracks; i3++)
   {
    //if(i3==fSelectedTracks-1){cout<<"done with loop #3"<<endl;}
    if(i3==i1||i3==i2){continue;}
    Double_t dPhi3 = ftaNestedLoops[0]->GetAt(i3);
    Double_t dW3 = ftaNestedLoops[1]->GetAt(i3);
    if(3==order)
    {
     value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3);   
     weight = dW1*dW2*dW3;
     profile->Fill(0.5,value,weight);
     continue;
    }
    for(int i4=0; i4<fSelectedTracks; i4++)
    {
     //if(i4==fSelectedTracks-1){cout<<"done with loop #4"<<endl;}
     if(i4==i1||i4==i2||i4==i3){continue;}
     Double_t dPhi4 = ftaNestedLoops[0]->GetAt(i4);
     Double_t dW4 = ftaNestedLoops[1]->GetAt(i4);
     if(4==order)
     {
      value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4);   
      weight = dW1*dW2*dW3*dW4;
      profile->Fill(0.5,value,weight);
      continue;
     }
     for(int i5=0; i5<fSelectedTracks; i5++)
     {
      //if(i5==fSelectedTracks-1){cout<<"done with loop #5"<<endl;}
      if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
      Double_t dPhi5 = ftaNestedLoops[0]->GetAt(i5);
      Double_t dW5 = ftaNestedLoops[1]->GetAt(i5);
      if(5==order)
      {
       value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5);   
       weight = dW1*dW2*dW3*dW4*dW5;
       profile->Fill(0.5,value,weight);
       continue;
      }
      for(int i6=0; i6<fSelectedTracks; i6++)
      {
       //if(i6==fSelectedTracks-1){cout<<"done with loop #6"<<endl;}
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
       Double_t dPhi6 = ftaNestedLoops[0]->GetAt(i6);
       Double_t dW6 = ftaNestedLoops[1]->GetAt(i6);
       if(6==order)
       {
        value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                         + harmonics->GetAt(5)*dPhi6);   
        weight = dW1*dW2*dW3*dW4*dW5*dW6;
        profile->Fill(0.5,value,weight);
        continue;
       }
       for(int i7=0; i7<fSelectedTracks; i7++)
       {
        //if(i7==fSelectedTracks-1){cout<<"done with loop #7"<<endl;}
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
        Double_t dPhi7 = ftaNestedLoops[0]->GetAt(i7);
        Double_t dW7 = ftaNestedLoops[1]->GetAt(i7);
        if(7==order)
        {
         value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                          + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7);   
         weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7;
         profile->Fill(0.5,value,weight);
         continue;
        }
        for(int i8=0; i8<fSelectedTracks; i8++)
        {
         //if(i8==fSelectedTracks-1){cout<<"done with loop #8"<<endl;}
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7){continue;}
         Double_t dPhi8 = ftaNestedLoops[0]->GetAt(i8);
         Double_t dW8 = ftaNestedLoops[1]->GetAt(i8);
         if(8==order)
         {
          value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                           + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7 + harmonics->GetAt(7)*dPhi8);   
          weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8;
          profile->Fill(0.5,value,weight);
          continue;
         }
         for(int i9=0; i9<fSelectedTracks; i9++)
         {
          //if(i9==fSelectedTracks-1){cout<<"done with loop #9"<<endl;}
          if(i9==i1||i9==i2||i9==i3||i9==i4||i9==i5||i9==i6||i9==i7||i9==i8){continue;}
          Double_t dPhi9 = ftaNestedLoops[0]->GetAt(i9);
          Double_t dW9 = ftaNestedLoops[1]->GetAt(i9);
          if(9==order)
          {
           value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                            + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7 + harmonics->GetAt(7)*dPhi8 + harmonics->GetAt(8)*dPhi9);   
           weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8*dW9;
           profile->Fill(0.5,value,weight);
           continue;
          }
          for(int i10=0; i10<fSelectedTracks; i10++)
          {
           //if(i10==fSelectedTracks-1){cout<<"done with loop #10"<<endl;}
           if(i10==i1||i10==i2||i10==i3||i10==i4||i10==i5||i10==i6||i10==i7||i10==i8||i10==i9){continue;}
           Double_t dPhi10 = ftaNestedLoops[0]->GetAt(i10);
           Double_t dW10 = ftaNestedLoops[1]->GetAt(i10);
           if(10==order)
           {
            value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                             + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7 + harmonics->GetAt(7)*dPhi8 + harmonics->GetAt(8)*dPhi9 + harmonics->GetAt(9)*dPhi10);   
            weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8*dW9*dW10;
            profile->Fill(0.5,value,weight);
            continue;
           }
           for(int i11=0; i11<fSelectedTracks; i11++)
           {
            //if(i11==fSelectedTracks-1){cout<<"done with loop #11"<<endl;}
            if(i11==i1||i11==i2||i11==i3||i11==i4||i11==i5||i11==i6||i11==i7||i11==i8||i11==i9||i11==i10){continue;}
            Double_t dPhi11 = ftaNestedLoops[0]->GetAt(i11);
            Double_t dW11 = ftaNestedLoops[1]->GetAt(i11);
            if(11==order)
            {
             value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                              + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7 + harmonics->GetAt(7)*dPhi8 + harmonics->GetAt(8)*dPhi9 + harmonics->GetAt(9)*dPhi10
                              + harmonics->GetAt(10)*dPhi11);   
             weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8*dW9*dW10*dW11;
             profile->Fill(0.5,value,weight);
             continue;
            }
            for(int i12=0; i12<fSelectedTracks; i12++)
            {
             //if(i12==fSelectedTracks-1){cout<<"done with loop #12"<<endl;}
             if(i12==i1||i12==i2||i12==i3||i12==i4||i12==i5||i12==i6||i12==i7||i12==i8||i12==i9||i12==i10||i12==i11){continue;}
             Double_t dPhi12 = ftaNestedLoops[0]->GetAt(i12);
             Double_t dW12 = ftaNestedLoops[1]->GetAt(i12);
             if(12==order)
             {
              value = TMath::Cos(harmonics->GetAt(0)*dPhi1 + harmonics->GetAt(1)*dPhi2 + harmonics->GetAt(2)*dPhi3 + harmonics->GetAt(3)*dPhi4 + harmonics->GetAt(4)*dPhi5
                               + harmonics->GetAt(5)*dPhi6 + harmonics->GetAt(6)*dPhi7 + harmonics->GetAt(7)*dPhi8 + harmonics->GetAt(8)*dPhi9 + harmonics->GetAt(9)*dPhi10
                               + harmonics->GetAt(10)*dPhi11 + harmonics->GetAt(11)*dPhi12);   
              weight = dW1*dW2*dW3*dW4*dW5*dW6*dW7*dW8*dW9*dW10*dW11*dW12;
              profile->Fill(0.5,value,weight);
              continue;
             }

             // ... it's easy to continue the above pattern here

            } // for(int i12=0; i12<fSelectedTracks; i12++)
           } // for(int i11=0; i11<fSelectedTracks; i11++)
          } // for(int i10=0; i10<fSelectedTracks; i10++)
         } // for(int i9=0; i9<fSelectedTracks; i9++)
        } // for(int i8=0; i8<fSelectedTracks; i8++)
       } // for(int i7=0; i7<fSelectedTracks; i7++)
      } // for(int i6=0; i6<fSelectedTracks; i6++)
     } // for(int i5=0; i5<fSelectedTracks; i5++)
    } // for(int i4=0; i4<fSelectedTracks; i4++)   
   } // for(int i3=0; i3<fSelectedTracks; i3++)
  } // for(int i2=0; i2<fSelectedTracks; i2++)
 } // for(int i1=0; i1<fSelectedTracks; i1++)

 // c) Return value:
 Double_t finalValue = profile->GetBinContent(1);
 delete profile;
 return finalValue;
 
} // Double_t AliAnalysisTaskMuPa::CalculateCustomNestedLoop(TArrayI *harmonics)

//=======================================================================================================================

void AliAnalysisTaskMuPa::ComparisonNestedLoopsVsCorrelations()
{
 // Make a ratio fNestedLoopsPro[....]/fCorrelationsPro[....]. If results are the same, these ratios must be 1.

 // a) Integrated comparison;
 // b) Comparison vs. multiplicity;
 // c) Comparison vs. centrality;

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

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
  cout<<Form("   ==== <<%d>>-particle correlations ====",2*(o+1))<<endl;
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
  cout<<Form("   ==== <<%d>>-particle correlations ====",2*(o+1))<<endl;
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
  cout<<Form("   ==== <<%d>>-particle correlations ====",2*(o+1))<<endl;
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

//=======================================================================================================================

Double_t AliAnalysisTaskMuPa::Weight(const Double_t &value, const char *variable) // value, [phi,pt,eta]
{
 // Determine particle weight. 

 // Basic protection:
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 if(!fWeightsHist[ppe]){cout<<__LINE__<<endl;exit(1);}

 Int_t bin = fWeightsHist[ppe]->FindBin(value);
 Double_t weight = 0.; 
 if(bin > fWeightsHist[ppe]->GetNbinsX())
 {
  weight = 0.; // we are in the overflow, ignore this particle TBI_20210524 is this really the correct procedure?
 } 
 else
 {
  weight = fWeightsHist[ppe]->GetBinContent(bin);
 }
 
 return weight;

} // AliAnalysisTaskMuPa::Weight(const Double_t &value, const char *variable) // value, [phi,pt,eta]

//=======================================================================================================================

Double_t AliAnalysisTaskMuPa::CentralityWeight(const Double_t &value, const char *estimator) // value, [V0M, SPDTracklets, CL0, CL1]
{
 // Determine centrality weight. 

 // Basic protection:
 if(!(TString(estimator).EqualTo("V0M") || TString(estimator).EqualTo("SPDTracklets") || TString(estimator).EqualTo("CL0") || TString(estimator).EqualTo("CL1"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ce = -1; // this part has to be in sync. with enum eCentralityEstimator
 if(TString(estimator).EqualTo("V0M")){ce=V0M;} 
 if(TString(estimator).EqualTo("SPDTracklets")){ce=SPDTracklets;} 
 if(TString(estimator).EqualTo("CL0")){ce=CL0;} 
 if(TString(estimator).EqualTo("CL1")){ce=CL1;} 

 if(!fCentralityWeightsHist[ce]){cout<<__LINE__<<endl;exit(1);}

 Int_t bin = fCentralityWeightsHist[ce]->FindBin(value);
 Double_t weight = 0.; 
 if(bin > fCentralityWeightsHist[ce]->GetNbinsX())
 {
  weight = 0.; // we are in the overflow, ignore this particle TBI_20210524 is this really the correct procedure?
 } 
 else
 {
  weight = fCentralityWeightsHist[ce]->GetBinContent(bin)*fCentralityWeightsHist[ce]->GetBinWidth(bin); // yes, since fCentralityWeightsHist is normalized p.d.f. (ensure that with the macro)
 }
  
 // In this context, it is assumed that centrality weight is a normalized probability (ensure that with the macro):
 if(weight < 0. || weight > 1.){Red(Form("\n weight = %f \n",weight));cout<<__LINE__<<endl;exit(1);}  

 return weight;

} // AliAnalysisTaskMuPa::Weight(const Double_t &value, const char *estimator) // value, [V0M, SPDTracklets, CL0, CL1]

//=======================================================================================================================

void AliAnalysisTaskMuPa::SetWeightsHist(TH1D* const hist, const char *variable)
{
 // Copy histogram holding weights from an external file to the corresponding data member. 
  
 // Basic protection:
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ppe=-1;
 if(TString(variable).EqualTo("phi")){ppe=0;} 
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 // Finally:
 hist->SetDirectory(0);
 fWeightsHist[ppe] = (TH1D*)hist->Clone();
 if(!fWeightsHist[ppe]){cout<<__LINE__<<endl; exit(1);}

 // Flag:
 fUseWeights[ppe] = kTRUE; 

} // void AliAnalysisTaskMuPa::SetWeightsHist(TH1D* const hwh, const char *type, const char *variable)

//=======================================================================================================================

void AliAnalysisTaskMuPa::SetCentralityWeightsHist(TH1D* const hist, const char *estimator)
{
 // Copy histogram holding weights from an external file to the corresponding data member. 
  
 // Basic protection:
 if(!(TString(estimator).EqualTo("V0M") || TString(estimator).EqualTo("SPDTracklets") || TString(estimator).EqualTo("CL0") || TString(estimator).EqualTo("CL1"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ce = -1; // this part has to be in sync. with enum eCentralityEstimator
 if(TString(estimator).EqualTo("V0M")){ce=V0M;} 
 if(TString(estimator).EqualTo("SPDTracklets")){ce=SPDTracklets;} 
 if(TString(estimator).EqualTo("CL0")){ce=CL0;} 
 if(TString(estimator).EqualTo("CL1")){ce=CL1;} 

 // Finally:
 hist->SetDirectory(0);
 fCentralityWeightsHist[ce] = (TH1D*)hist->Clone();
 if(!fCentralityWeightsHist[ce]){cout<<__LINE__<<endl; exit(1);}

 // Flag:
 fUseCentralityWeights[ce] = kTRUE; 

} // void AliAnalysisTaskMuPa::SetWeightsHist(TH1D* const hwh, const char *type, const char *variable)

//=======================================================================================================================

TH1D* AliAnalysisTaskMuPa::GetWeightsHist(const char *variable)
{
 // The standard getter. 
  
 // Basic protection:
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ppe=-1;
 if(TString(variable).EqualTo("phi")){ppe=0;} 
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 // Finally:
 return fWeightsHist[ppe];

} // TH1D* AliAnalysisTaskMuPa::GetWeightsHist(const char *variable)

//=======================================================================================================================

TH1D* AliAnalysisTaskMuPa::GetCentralityWeightsHist(const char *estimator)
{
 // The standard getter. 
  
 // Basic protection:
 if(!(TString(estimator).EqualTo("V0M") || TString(estimator).EqualTo("SPDTracklets") || TString(estimator).EqualTo("CL0") || TString(estimator).EqualTo("CL1"))){cout<<__LINE__<<endl;exit(1);}

 Int_t ce = -1; //  this part has to be in sync. with enum eCentralityEstimator
 if(TString(estimator).EqualTo("V0M")){ce=V0M;} 
 if(TString(estimator).EqualTo("SPDTracklets")){ce=SPDTracklets;} 
 if(TString(estimator).EqualTo("CL0")){ce=CL0;} 
 if(TString(estimator).EqualTo("CL1")){ce=CL1;} 

 // Finally:
 return fCentralityWeightsHist[ce];

} // TH1D* AliAnalysisTaskMuPa::GetWeightsHist(const char *estimator)

//=======================================================================================

TH1D *AliAnalysisTaskMuPa::GetHistogramWithWeights(const char *filePath, const char *variable)
{
 // Access from external ROOT file the desired histogram with particle weights. 
 // 'filePath' can be both abs and relative path (e.g. pwd)

 // a) Return value; 
 // b) Basic protection for arguments; 
 // c) Check if the external ROOT file exists at specified path; 
 // d) Access the external ROOT file and fetch the desired histogram with weights;
 // e) Close the external ROOT file. 

 // a) Return value:
 TH1D *hist = NULL; 

 // b) Basic protection for arguments:
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){cout<<__LINE__<<endl;exit(1);}

 // c) Check if the external ROOT file exists at specified path:
 if(gSystem->AccessPathName(filePath,kFileExists))
 {
  Red(Form("if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s",filePath)); 
  cout<<__LINE__<<endl;
  exit(1);
 }

 // d) Access the external ROOT file and fetch the desired histogram with weights:
 TFile *weightsFile = TFile::Open(filePath,"READ");
 if(!weightsFile){cout<<__LINE__<<endl;exit(1);}
 hist = (TH1D*)(weightsFile->Get(Form("%s_%s",variable,fTaskName.Data())));
 if(!hist){hist = (TH1D*)(weightsFile->Get(Form("%s",variable)));} // yes, for some simple tests I can have only histogram named e.g. 'phi'
 if(!hist){Red(Form("%s_%s",variable,fTaskName.Data())); cout<<__LINE__<<endl;exit(1);}
 hist->SetDirectory(0);
 hist->SetTitle(filePath);

 // e) Close the external ROOT file: 
 weightsFile->Close(); delete weightsFile;

 return hist;

} // TH1D *AliAnalysisTaskMuPa::GetHistogramWithWeights(const char *filePath, const char *variable)

//=======================================================================================

TH1D *AliAnalysisTaskMuPa::GetHistogramWithCentralityWeights(const char *filePath, const char *estimator)
{
 // Access from external ROOT file the desired histogram with particle weights. 
 // 'filePath' can be both abs and relative path (e.g. pwd)

 // a) Return value; 
 // b) Basic protection for arguments; 
 // c) Check if the external ROOT file exists at specified path; 
 // d) Access the external ROOT file and fetch the desired histogram with weights;
 // e) Close the external ROOT file. 

 // a) Return value:
 TH1D *hist = NULL; 

 // b) Basic protection for arguments:
 if(!(TString(estimator).EqualTo("V0M") || TString(estimator).EqualTo("SPDTracklets") || TString(estimator).EqualTo("CL0") || TString(estimator).EqualTo("CL1"))){cout<<__LINE__<<endl;exit(1);}

 // c) Check if the external ROOT file exists at specified path:
 if(gSystem->AccessPathName(filePath,kFileExists))
 {
  Red(Form("if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s",filePath));
  cout<<__LINE__<<endl;
  exit(1);
 }

 // d) Access the external ROOT file and fetch the desired histogram with weights:
 TFile *weightsFile = TFile::Open(filePath,"READ");
 if(!weightsFile){cout<<__LINE__<<endl;exit(1);}
 hist = (TH1D*)(weightsFile->Get(Form("%s_%s",estimator,fTaskName.Data())));
 if(!hist){hist = (TH1D*)(weightsFile->Get(Form("%s",estimator)));} // yes, for some simple tests I can have only histogram named e.g. 'phi'
 if(!hist){Red(Form("%s_%s",estimator,fTaskName.Data())); cout<<__LINE__<<endl;exit(1);}
 hist->SetDirectory(0);
 hist->SetTitle(filePath);

 // e) Close the external ROOT file: 
 weightsFile->Close(); delete weightsFile;

 return hist;

} // TH1D *AliAnalysisTaskMuPa::GetHistogramWithCentralityWeights(const char *filePath, const char *estimator)

//=======================================================================================

Bool_t AliAnalysisTaskMuPa::SpecifiedEvent(AliVEvent *ave)
{
 // Check if this is the event specified in a steering macro via the setter void SetProcessOnlySpecifiedEvent(Int_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period).

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Wait for specified event.

 if(fVerbose){Yellow(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 { 
  if(aAOD->GetRunNumber() != fRun) return kFALSE;
  else if(aAOD->GetBunchCrossNumber() != fBunchCross) return kFALSE;
  else if(aAOD->GetOrbitNumber() != fOrbit) return kFALSE;
  else if(aAOD->GetPeriodNumber() != fPeriod) return kFALSE;
 }

 return kTRUE;

} // void AliAnalysisTaskMuPa::SpecifiedEvent(Int_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)

//=======================================================================================

void AliAnalysisTaskMuPa::PrintEventInfo(AliVEvent *ave)
{
 // Print event metadata. Used for debugging. Enable via task->PrintEventInfo()
 
 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Wait for specified event.

 if(fVerbose){Yellow(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 //AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aAOD)
 {
  Yellow(Form("aAOD->GetRunNumber() = %d",aAOD->GetRunNumber()));
  Yellow(Form("aAOD->GetBunchCrossNumber() = %d",aAOD->GetBunchCrossNumber()));
  Yellow(Form("aAOD->GetOrbitNumber() = %d",aAOD->GetOrbitNumber()));
  Yellow(Form("aAOD->GetPeriodNumber() = %d",aAOD->GetPeriodNumber()));
 } // if(aAOD)

} // void AliAnalysisTaskMuPa::PrintEventInfo(AliVEvent *ave)

//=======================================================================================

void AliAnalysisTaskMuPa::Red(const char* text)
{ 
 cout<<"\n\033[1;31m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

void AliAnalysisTaskMuPa::Green(const char* text)
{ 
 cout<<"\n\033[1;32m"<<text<<"\033[0m\n"<<endl;
}
//=======================================================================================

void AliAnalysisTaskMuPa::Yellow(const char* text)
{ 
 cout<<"\n\033[1;33m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

void AliAnalysisTaskMuPa::Blue(const char* text)
{ 
 cout<<"\n\033[1;34m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

void AliAnalysisTaskMuPa::MakeLookUpTable(AliAODEvent *aAOD, AliMCEvent *aMC)
{
 // For Monte Carlo analysis, establish a look up table between reco and kine particles.
 // Algorithm: looping first over the reconstructed particles and building a look-up table which stores for each label index (kine) the corresponding track index (reco)
 
 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 if(!aAOD){cout<<__LINE__<<endl;exit(1);}
 if(!aMC){cout<<__LINE__<<endl;exit(1);}
 if(0 != fSimReco->GetSize()){fSimReco->Delete();} // yes, this method determines mapping from scratch each time

 AliAODTrack *aodTrack = NULL;
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all reconstructed tracks
 {
  aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(!aodTrack){continue;}

  if(!SurvivesParticleCuts(aodTrack)){continue;} // without this, I have problenm with TPC-only and Global tracks, since they have the same Monte Carlo label
  
  Int_t label = aodTrack->GetLabel(); // Monte Carlo label, i.e. index of the particle in AliMCEvent *aMC
  // The negative labels refer to fake tracks (see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TrackParametersMCTruth)
  // * The absolute value of this label is the index of the particle within the *
  // * MC stack. If the label is negative, this track was assigned a certain *
  // * number of clusters that did not in fact belong to this track. *
  // => In most analysis, fake tracks can be taken, use task->SetUseFakeTracks(kTRUE); to check the difference
  fSimReco->Add(label,iTrack); // "key" = label, "value" = iTrack

  /* 
  // Use this code snippet to cross-check that the mapping went allright:
  AliAODMCParticle *aodmcParticle = NULL;
  Int_t nTracksMC = aMC->GetNumberOfTracks(); // total number of Monte Carlo tracks
  for(Int_t iTrackMC=0;iTrackMC<nTracksMC;iTrackMC++) // starting a loop over all Monte Carlo tracks
  {
   if(TMath::Abs(label) == iTrackMC)
   {
    aodmcParticle = (AliAODMCParticle*)aMC->GetTrack(iTrackMC);    
    cout<<aodTrack->Phi()<<" "<<aodmcParticle->Phi()<<endl;
    cout<<aodTrack->Pt()<<" "<<aodmcParticle->Pt()<<endl;
    cout<<aodTrack->Charge()<<" "<<aodmcParticle->Charge()/3<<endl;
    cout<<endl;
    break;
   }
  }
  */

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all reconstructed tracks

} // void AliAnalysisTaskMuPa::MakeLookUpTable(AliAODEvent *aAOD, AliMCEvent *aMC)

//=======================================================================================

void AliAnalysisTaskMuPa::FillQAHistograms(AliAODEvent *aAOD, AliMCEvent *aMC)
{
 // Use this member function only when both AliAODEvent and AliMCEvent are needed. 
 // Alternatively, use FillQAHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs);

 // a) Check for self-correlations between simulated and reconstructed particles. I need only one loop here.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 if(!aAOD){cout<<__LINE__<<endl;exit(1);}
 if(!aMC){cout<<__LINE__<<endl;exit(1);}

 // a) Check for self-correlations between simulated and reconstructed particles. I need only one loop here:
 if(fQACheckSelfCorrelations)
 { 
  AliAODTrack *aodTrack = NULL;
  AliAODMCParticle *aodmcParticle = NULL;
  Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all reconstructed tracks
  {
   aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
   if(!aodTrack){continue;}
   if(!SurvivesParticleCuts(aodTrack)){continue;} // without this, I have problem with TPC-only and Global tracks, since they have the same Monte Carlo label

   // The corresponding Monte Carlo particle:
   aodmcParticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(aodTrack->GetLabel()));
   if(!aodmcParticle && fUseFakeTracks)
   {
    aodmcParticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(TMath::Abs(aodTrack->GetLabel())));
   }
   if(!aodmcParticle){cout<<__LINE__<<endl;exit(1);} // aod track doesn't have corresponding Monte Carlo particle
   if(fQASimRecoSelfCorrelations[0]){fQASimRecoSelfCorrelations[0]->Fill(aodTrack->Phi()-aodmcParticle->Phi());}
   if(fQASimRecoSelfCorrelations[1]){fQASimRecoSelfCorrelations[1]->Fill(aodTrack->Pt()-aodmcParticle->Pt());}  
   if(fQASimRecoSelfCorrelations[2]){fQASimRecoSelfCorrelations[2]->Fill(aodTrack->Eta()-aodmcParticle->Eta());}
  } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all reconstructed tracks
 } // if(fQACheckSelfCorrelations)

} // void AliAnalysisTaskMuPa::FillQAHistograms(AliAODEvent *aAOD, AliMCEvent *aMC)

//=======================================================================================

void AliAnalysisTaskMuPa::RandomIndices(AliVEvent *ave)
{
 // Randomize indices using Fisher-Yates algorithm. 

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Get total number of tracks;
 // c) Fisher-Yates algorithm.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) Get total number of tracks:
 Int_t nTracks = 0;
 if(aAOD)
 {
  nTracks = aAOD->GetNumberOfTracks();
 }
 else if(aMC)
 {
  nTracks = aMC->GetNumberOfTracks();
 }

 if(nTracks<1){return;}

 // c) Fisher-Yates algorithm:
 fRandomIndices = new TArrayI(nTracks);
 fRandomIndices->Reset(); 
 for(Int_t i=0;i<nTracks;i++)
 {
  fRandomIndices->AddAt(i,i);
 }
 for(Int_t i=nTracks-1;i>=1;i--)
 {
  Int_t j = gRandom->Integer(i+1);
  Int_t temp = fRandomIndices->GetAt(j);
  fRandomIndices->AddAt(fRandomIndices->GetAt(i),j);
  fRandomIndices->AddAt(temp,i);
 } // end of for(Int_t i=nTracks-1;i>=1;i--) 

} // void AliAnalysisTaskMuPa::RandomIndices(AliVEvent *ave)

//=======================================================================================

void AliAnalysisTaskMuPa::GenerateCorrelationsLabels()
{
 // Generate the labels of all correlations of interest.
 // Algorithm: a) If external file is specified, parse through it and fetch line by line all labels. FS is comma.
 //            b) If external file is NOT specified, determine the labels from the hardcoded loops below.

 Int_t counter[gMaxCorrelator] = {0}; // is this safe?
 for(Int_t o=0;o<gMaxCorrelator;o++){counter[o] = 0;} // now it's safe
  
 if(fFileWithLabels && !gSystem->AccessPathName(fFileWithLabels->Data(),kFileExists))
 {
  // external file exist locally, get all labels of interest from there:
  string line;
  ifstream myfile;
  myfile.open(fFileWithLabels->Data());
  Int_t order = -44;
  while (getline(myfile,line))
  { 
   order = TString(line).Tokenize(" ")->GetEntries();
   if(0 == order){continue;} // empty lines, or the label format which is not supported
   // 1-p => 0, 2-p => 1, etc.:
   fTest0Labels[order-1][counter[order-1]] = new TString(line); // okay...
   //cout<<fTest0Labels[order-1][counter[order-1]]->Data()<<endl;
   counter[order-1]++;
   //cout<<TString(line).Data()<<endl;
   //cout<<oa->GetEntries()<<endl;    
  }
  myfile.close();
 }
 else
 {
  // get all labels of interest systematically from loops below:
  // TBI 20210902 buggy, some entries are duplicated

  // 1p:
  // TBI 20210902 
  // fTest0Labels[0][index1p] = ... 

  // 2p:
  Int_t index2p = 0;
  for(Int_t h1=-fMaxHarmonic;h1<=fMaxHarmonic;h1++) 
  {  
   if(0==h1){continue;}
   for(Int_t h2=-fMaxHarmonic;h2<=fMaxHarmonic;h2++) 
   {
    if(0==h2){continue;}
    if(h2>h1){continue;} // eliminating trivial permutations
    if(abs(h2)>abs(h1)){continue;} // eliminating -1*(...) symmetry
    if(0 != h1+h2){continue;} // isotropy
    fTest0Labels[1][index2p] = new TString(Form("%d,%d",h1,h2));
    index2p++;
    //cout<<Form("%d,%d",h1,h2)<<endl;
   }
  }

  // 3p:
  Int_t index3p = 0;
  for(Int_t h1=-fMaxHarmonic;h1<=fMaxHarmonic;h1++) 
  { 
   if(0==h1){continue;}
   for(Int_t h2=-fMaxHarmonic;h2<=fMaxHarmonic;h2++) 
   {
    if(0==h2){continue;}
    if(h2>h1){continue;} // eliminating trivial permutations
    if(abs(h2)>abs(h1)){continue;} // eliminating -1*(...) symmetry
    for(Int_t h3=-fMaxHarmonic;h3<=fMaxHarmonic;h3++) 
    {
     if(0==h3){continue;}
     if(h3>h1||h3>h2){continue;} // eliminating trivial permutations
     if(abs(h3)>abs(h1)||abs(h2)>abs(h1)){continue;} // eliminating -1*(...) symmetry
     if(0 != h1+h2+h3){continue;} // isotropy
     fTest0Labels[2][index3p] = new TString(Form("%d,%d,%d",h1,h2,h3));
     index3p++;
     //cout<<Form("%d,%d,%d",h1,h2,h3)<<endl;
    }
   }
  }

  // 4p:
  Int_t index4p = 0;
  for(Int_t h1=-fMaxHarmonic;h1<=fMaxHarmonic;h1++) 
  { 
   if(0==h1){continue;}
   for(Int_t h2=-fMaxHarmonic;h2<=fMaxHarmonic;h2++) 
   {
    if(0==h2){continue;}
    if(h2>h1){continue;} // eliminating trivial permutations
    for(Int_t h3=-fMaxHarmonic;h3<=fMaxHarmonic;h3++) 
    {
     if(0==h3){continue;}
     if(h3>h1||h3>h2){continue;} // eliminating trivial permutations
     if(abs(h3)>abs(h1)){continue;} // eliminating -1*(...) symmetry
     for(Int_t h4=-fMaxHarmonic;h4<=fMaxHarmonic;h4++) 
     {
      if(0==h4){continue;}
      if(h4>h1||h4>h2||h4>h3){continue;} // eliminating trivial permutations
      if(abs(h4)>abs(h1)){continue;} // eliminating -1*(...) symmetry
      if(0 != h1+h2+h3+h4){continue;} // isotropy
      fTest0Labels[3][index4p] = new TString(Form("%d,%d,%d,%d",h1,h2,h3,h4));
      index4p++;
      //cout<<Form("%d,%d,%d,%d",h1,h2,h3,h4)<<endl;
     }
    }
   }
  }

  // 5p:
  Int_t index5p = 0;
  for(Int_t h1=-fMaxHarmonic;h1<=fMaxHarmonic;h1++) 
  { 
   if(0==h1){continue;}
   for(Int_t h2=-fMaxHarmonic;h2<=fMaxHarmonic;h2++) 
   { 
    if(0==h2){continue;}
    if(h2>h1){continue;} // eliminating trivial permutations
    for(Int_t h3=-fMaxHarmonic;h3<=fMaxHarmonic;h3++) 
    {
     if(0==h3){continue;}
     if(h3>h1||h3>h2){continue;} // eliminating trivial permutations
     if(abs(h3)>abs(h1)){continue;} // eliminating -1*(...) symmetry
     for(Int_t h4=-fMaxHarmonic;h4<=fMaxHarmonic;h4++) 
     {
      if(0==h4){continue;}
      if(h4>h1||h4>h2||h4>h3){continue;} // eliminating trivial permutations
      if(abs(h4)>abs(h1)){continue;} // eliminating -1*(...) symmetry
      for(Int_t h5=-fMaxHarmonic;h5<=fMaxHarmonic;h5++) 
      {
       if(0==h5){continue;}
       if(h5>h1||h5>h2||h5>h3||h5>h4){continue;} // eliminating trivial permutations
       if(abs(h5)>abs(h1)){continue;} // eliminating -1*(...) symmetry
       if(0 != h1+h2+h3+h4+h5){continue;} // isotropy
       fTest0Labels[4][index5p] = new TString(Form("%d,%d,%d,%d,%d",h1,h2,h3,h4,h5));
       index5p++;
       //cout<<Form("%d,%d,%d,%d,%d",h1,h2,h3,h4,h5)<<endl;
      }
     }
    }
   }
  }

  // 6p:
  Int_t index6p = 0;
  for(Int_t h1=-fMaxHarmonic;h1<=fMaxHarmonic;h1++) 
  { 
   if(0==h1){continue;}
   for(Int_t h2=-fMaxHarmonic;h2<=fMaxHarmonic;h2++) 
   { 
    if(0==h2){continue;}
    if(h2>h1){continue;} // eliminating trivial permutations
    for(Int_t h3=-fMaxHarmonic;h3<=fMaxHarmonic;h3++) 
    {
     if(0==h3){continue;}
     if(h3>h1||h3>h2){continue;} // eliminating trivial permutations
     if(abs(h3)>abs(h1)){continue;} // eliminating -1*(...) symmetry
     for(Int_t h4=-fMaxHarmonic;h4<=fMaxHarmonic;h4++) 
     {
      if(0==h4){continue;}
      if(h4>h1||h4>h2||h4>h3){continue;} // eliminating trivial permutations
      if(abs(h4)>abs(h1)){continue;} // eliminating -1*(...) symmetry
      for(Int_t h5=-fMaxHarmonic;h5<=fMaxHarmonic;h5++) 
      {
       if(0==h5){continue;}
       if(h5>h1||h5>h2||h5>h3||h5>h4){continue;} // eliminating trivial permutations
       if(abs(h5)>abs(h1)){continue;} // eliminating -1*(...) symmetry
       for(Int_t h6=-fMaxHarmonic;h6<=fMaxHarmonic;h6++) 
       {
        if(0==h6){continue;}
        if(h6>h1||h6>h2||h6>h3||h6>h4||h6>h5){continue;} // eliminating trivial permutations
        if(abs(h6)>abs(h1)){continue;} // eliminating -1*(...) symmetry
        if(0 != h1+h2+h3+h4+h5+h6){continue;} // isotropy
        fTest0Labels[5][index6p] = new TString(Form("%d,%d,%d,%d,%d,%d",h1,h2,h3,h4,h5,h6));
        index6p++;
        //cout<<Form("%d,%d,%d,%d,%d,%d",h1,h2,h3,h4,h5,h6)<<endl;
       }
      }
     }
    }
   }
  }

 } // else

} // void AliAnalysisTaskMuPa::GenerateCorrelationsLabels() 

//=======================================================================================

void AliAnalysisTaskMuPa::CalculateTest0()
{
 // Calculate Test0.

 Double_t correlation = 0.;
 Double_t weight = 0.;
 Int_t n[gMaxCorrelator] = {0}; // array holding harmonics

 for(Int_t mo=0;mo<gMaxCorrelator;mo++) 
 { 
  for(Int_t mi=0;mi<gMaxIndex;mi++) 
  { 
   if(fTest0Labels[mo][mi])
   {
    // Extract harmonics from TString, FS is " ": 
    for(Int_t h=0;h<=mo;h++)
    {
     n[h] = TString(fTest0Labels[mo][mi]->Tokenize(" ")->At(h)->GetName()).Atoi(); // okay...
    }

    switch(mo+1) // which order? yes, mo+1
    {
     case 1:
      correlation = One(n[0]).Re();
      weight = One(n[0]).Re();
     break;

     case 2: 
      correlation = Two(n[0],n[1]).Re();
      weight = Two(0,0).Re();
     break;

     case 3: 
      correlation = Three(n[0],n[1],n[2]).Re();
      weight = Three(0,0,0).Re();
     break;
 
     case 4: 
      correlation = Four(n[0],n[1],n[2],n[3]).Re();
      weight = Four(0,0,0,0).Re();
     break;

     case 5: 
      correlation = Five(n[0],n[1],n[2],n[3],n[4]).Re();
      weight = Five(0,0,0,0,0).Re();
     break;

     case 6: 
      correlation = Six(n[0],n[1],n[2],n[3],n[4],n[5]).Re();
      weight = Six(0,0,0,0,0,0).Re();
     break;

     case 7: 
      correlation = Seven(n[0],n[1],n[2],n[3],n[4],n[5],n[6]).Re();
      weight = Seven(0,0,0,0,0,0,0).Re();
     break;

     case 8: 
      correlation = Eight(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]).Re();
      weight = Eight(0,0,0,0,0,0,0,0).Re();
     break;

     case 9: 
      correlation = Nine(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8]).Re();
      weight = Nine(0,0,0,0,0,0,0,0,0).Re();
     break;

     case 10: 
      correlation = Ten(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9]).Re();
      weight = Ten(0,0,0,0,0,0,0,0,0,0).Re();
     break;

     case 11: 
      correlation = Eleven(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10]).Re();
      weight = Eleven(0,0,0,0,0,0,0,0,0,0,0).Re();
     break;

     case 12: 
      correlation = Twelve(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11]).Re();
      weight = Twelve(0,0,0,0,0,0,0,0,0,0,0,0).Re();
     break;

     default:
      cout<<fTest0Labels[mo][mi]->Data()<<endl;
      cout<<"not supported yet"<<endl;
      return; // TBI 20210907 or continue?
    } // switch(mo+1)

    // e-b-e sanity check:
    if(fCalculateCustomNestedLoop)
    {
     TArrayI *harmonics = new TArrayI(mo+1);
     for(Int_t i=0;i<mo+1;i++)
     {
      harmonics->SetAt(n[i],i);
     }
     if(!(weight>0.))
     {
      cout<<fTest0Labels[mo][mi]->Data()<<endl;   
      cout<<__LINE__<<endl; exit(1); 
     }
     if(TMath::Abs(correlation/weight - this->CalculateCustomNestedLoop(harmonics))>1.e-5)
     {       
      cout<<fTest0Labels[mo][mi]->Data()<<endl;   
      cout<<"correlation: "<<correlation/weight<<endl;   
      cout<<"custom loop: "<<this->CalculateCustomNestedLoop(harmonics)<<endl;   
      cout<<__LINE__<<endl; exit(1);
     }
     else
     {
      cout<<Form("=> e-b-e check with CustomNestedLoop is OK for %d-p Test0 corr. %s",mo+1,fTest0Labels[mo][mi]->Data())<<endl;
     }
     delete harmonics;
    } // if(fCalculateCustomNestedLoop)
 
    // Finally, fill:
    for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
    { 
     if(!(weight > 0.)){cout<<__LINE__<<endl;exit(1);}
     // integrated:
     if(fTest0Pro[mo][mi][0]){fTest0Pro[mo][mi][0]->Fill(0.5,correlation/weight,weight);}
     // vs. multiplicity:
     if(fTest0Pro[mo][mi][1]){fTest0Pro[mo][mi][1]->Fill(fSelectedTracks+0.5,correlation/weight,weight);}
     // vs. centrality:
     if(fTest0Pro[mo][mi][2]){fTest0Pro[mo][mi][2]->Fill(fCentrality,correlation/weight,weight);}
    } // for(Int_t v=0;v<3;v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
   } // if(fTest0Labels[mo][mi])
  } // for(Int_t mi=0;mi<gMaxIndex;mi++) 
 } // for(Int_t mo=0;mo<gMaxCorrelator;mo++) 

} // void AliAnalysisTaskMuPa::CalculateTest0()

//=======================================================================================

void AliAnalysisTaskMuPa::InternalValidation()
{
 // Internal validation against theoretical values in on-the-fly study for all implemented correlators. 

 // a) Configure Fourier like p.d.f. for azimuthal angles;
 // b) Loop over on-the-fly events.
 //    b0) Reset ebe quantities;
 //    b1) Determine multiplicity and reaction plane;
 //    b2) Loop over particles;
 //    b3) Calculate correlations;
 // c) Bail out directly from here when done;
 // d) Hasta la vista.

 Green(__PRETTY_FUNCTION__);

 // a) Configure Fourier like p.d.f. for azimuthal angles:
 TF1 *fPhiPDF = new TF1("fPhiPDF","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[3]*TMath::Cos(3.*(x-[0]))+2.*[4]*TMath::Cos(4.*(x-[0]))+2.*[5]*TMath::Cos(5.*(x-[0]))+2.*[6]*TMath::Cos(6.*(x-[0]))",0.,TMath::TwoPi());
 fPhiPDF->SetParName(0,"Reaction Plane");
 fPhiPDF->SetParameter(0,0.);
 for(Int_t h=1;h<=6;h++)
 {
  fPhiPDF->SetParName(h,Form("v_{%d}",h)); 
  if(h<=fInternalValidationHarmonics->GetSize())
  { 
   fPhiPDF->SetParameter(h,fInternalValidationHarmonics->GetAt(h-1));
  }
  else
  {
   fPhiPDF->SetParameter(h,0.);
  } 
 } // for(Int_t h=1;h<=6;h++)
 
 // b) Loop over on-the-fly events:
 Double_t step = 10.; // in percentage. Used only for the printout of progress
 TStopwatch watch;
 watch.Start();
 for(Int_t e=0;e<fnEventsInternalValidation;e++) 
 {
  if(1.*e/(fnEventsInternalValidation) > step/100.)
  {
   cout<<Form("Simulated %d%% events of requested %d",(Int_t)step,fnEventsInternalValidation)<<endl;
   watch.Print(); watch.Continue();
   step+=10.;
  } 

  // b0) Reset ebe quantities:
  this->ResetEventByEventQuantities();

  // b1) Determine multiplicity and reaction plane:
  Int_t nMult = gRandom->Uniform(fMultRangeInternalValidation[0],fMultRangeInternalValidation[1]);
  Double_t fReactionPlane = gRandom->Uniform(0.,TMath::TwoPi());
  fPhiPDF->SetParameter(0,fReactionPlane);

  // b2) Loop over particles:
  Double_t dPhi = 0.;
  for(Int_t p=0;p<nMult;p++) 
  {   
   // Particle angle:
   dPhi = fPhiPDF->GetRandom(); 
   // Fill Q-vector (simplified version, without weights):
   for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
   {
    for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
    {
     fQvector[h][wp] += TComplex(TMath::Cos(h*dPhi),TMath::Sin(h*dPhi));    
    } // for(Int_t wp=0;wp<fMaxCorrelator+1;wp++)
   } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)   

   // Nested loops containers: 
   if(fCalculateNestedLoops||fCalculateCustomNestedLoop)
   {
    if(ftaNestedLoops[0]){ftaNestedLoops[0]->AddAt(dPhi,p);} 
    if(ftaNestedLoops[1]){ftaNestedLoops[1]->AddAt(1.,p);} // yes, otherwise weights are automatically set to 0.
   }
  } // for(Int_t p=0;p<nMult;p++) 

  // b3) Calculate correlations:
  fSelectedTracks = nMult;
  fCentrality = gRandom->Uniform(0.,100.); // in any case it's meaningless in this exercise
  this->CalculateCorrelations();
  this->CalculateTest0();

  // b4) Optionally, cross-check with nested loops:
  if(fCalculateNestedLoops){this->CalculateNestedLoops();}

 } // for(Int_t e=0;e<fnEventsInternalValidation;e++) 

 // c) Bail out directly from here when done:
 if(fCalculateNestedLoops){this->ComparisonNestedLoopsVsCorrelations();}
 //    For the file name, I use again "AnalysisResults.root", not to bother with updating all scripts
 cout<<Form("\nInternal validation is over after %d events on-the-fly.\nDumping results in the file %s ....",fnEventsInternalValidation,"AnalysisResults.root")<<endl;
 sleep(2);
 TFile *f = new TFile("AnalysisResults.root","recreate");
 fBaseList->Write(fBaseList->GetName(),TObject::kSingleKey);
 f->Close();
 cout<<"Dumped!\n"<<endl;

 // d) Hasta la vista:
 exit(1);

} // void AliAnalysisTaskMuPa::InternalValidation()

















