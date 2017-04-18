/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// --- ROOT system ---
#include <TFile.h>
#include <TGeoManager.h>
#include <TStreamerInfo.h>

// ---- ANALYSIS system ----
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
//#include "AliTriggerAnalysis.h"
#include "AliESDVZERO.h"
#include "AliVCaloCells.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliMultSelection.h"

// ---- Detectors ----
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

// ---- CaloTrackCorr ---
#include "AliCalorimeterUtils.h"
#include "AliCaloTrackReader.h"

// ---- Jets ----
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"

/// \cond CLASSIMP
ClassImp(AliCaloTrackReader) ;
/// \endcond

//________________________________________
/// Constructor. Initialize parameters.
//________________________________________
AliCaloTrackReader::AliCaloTrackReader() :
TObject(),                   fEventNumber(-1), //fCurrentFileName(""),
fDataType(0),                fDebug(0),
fFiducialCut(0x0),           fCheckFidCut(kFALSE),
fComparePtHardAndJetPt(0),   fPtHardAndJetPtFactor(0),
fComparePtHardAndClusterPt(0),fPtHardAndClusterPtFactor(0),
fCTSPtMin(0),                fEMCALPtMin(0),                  fPHOSPtMin(0),
fCTSPtMax(0),                fEMCALPtMax(0),                  fPHOSPtMax(0),
fEMCALBadChMinDist(0),       fPHOSBadChMinDist (0),           
fEMCALNCellsCut(0),          fPHOSNCellsCut(0),
fUseEMCALTimeCut(1),         fUseParamTimeCut(0),
fUseTrackTimeCut(0),         fAccessTrackTOF(0),
fEMCALTimeCutMin(-10000),    fEMCALTimeCutMax(10000),
fEMCALParamTimeCutMin(),     fEMCALParamTimeCutMax(),
fTrackTimeCutMin(-10000),    fTrackTimeCutMax(10000),
fUseTrackDCACut(0),
fAODBranchList(0x0),
fCTSTracks(0x0),             fEMCALClusters(0x0),
fDCALClusters(0x0),          fPHOSClusters(0x0),
fEMCALCells(0x0),            fPHOSCells(0x0),
fInputEvent(0x0),            fOutputEvent(0x0),fMC(0x0),
fFillCTS(0),                 fFillEMCAL(0),
fFillDCAL(0),                fFillPHOS(0),
fFillEMCALCells(0),          fFillPHOSCells(0),
fRecalculateClusters(kFALSE),fCorrectELinearity(kTRUE),
fSelectEmbeddedClusters(kFALSE),
fSmearShowerShape(0),        fSmearShowerShapeWidth(0),       fRandom(),
fSmearingFunction(0),        fSmearNLMMin(0),                 fSmearNLMMax(0),
fTrackStatus(0),             fSelectSPDHitTracks(0),
fTrackMultNPtCut(0),         fTrackMultEtaCut(0.9),
fReadStack(kFALSE),          fReadAODMCParticles(kFALSE),
fDeltaAODFileName(""),       fFiredTriggerClassName(""),

fEventTriggerMask(0),        fMixEventTriggerMask(0),         fEventTriggerAtSE(0),
fEventTrigMinBias(0),        fEventTrigCentral(0),
fEventTrigSemiCentral(0),    fEventTrigEMCALL0(0),
fEventTrigEMCALL1Gamma1(0),  fEventTrigEMCALL1Gamma2(0),
fEventTrigEMCALL1Jet1(0),    fEventTrigEMCALL1Jet2(0),
fBitEGA(0),                  fBitEJE(0),

fEventType(-1),
fTaskName(""),               fCaloUtils(0x0),
fWeightUtils(0x0),           fEventWeight(1),
fMixedEvent(NULL),           fNMixedEvent(0),                 fVertex(NULL),
fListMixedTracksEvents(),    fListMixedCaloEvents(),
fLastMixedTracksEvent(-1),   fLastMixedCaloEvent(-1),
fWriteOutputDeltaAOD(kFALSE),
fEMCALClustersListName(""),  fZvtxCut(0.),
fAcceptFastCluster(kFALSE),  fRemoveLEDEvents(kFALSE),
//Trigger rejection
fRemoveBadTriggerEvents(0),  fTriggerPatchClusterMatch(1),
fTriggerPatchTimeWindow(),   fTriggerL0EventThreshold(0),
fTriggerL1EventThreshold(0), fTriggerL1EventThresholdFix(0),
fTriggerClusterBC(0),        fTriggerClusterIndex(0),         fTriggerClusterId(0),
fIsExoticEvent(0),           fIsBadCellEvent(0),              fIsBadMaxCellEvent(0),
fIsTriggerMatch(0),          fIsTriggerMatchOpenCut(),
fTriggerClusterTimeRecal(kTRUE), fRemoveUnMatchedTriggers(kTRUE),
fDoPileUpEventRejection(kFALSE), fDoV0ANDEventSelection(kFALSE),
fDoVertexBCEventSelection(kFALSE),
fDoRejectNoTrackEvents(kFALSE),
fUseEventsWithPrimaryVertex(kFALSE),
//fTriggerAnalysis (0x0),
fTimeStampEventSelect(0),
fTimeStampEventFracMin(0),   fTimeStampEventFracMax(0),
fTimeStampRunMin(0),         fTimeStampRunMax(0),
fNPileUpClusters(-1),        fNNonPileUpClusters(-1),         fNPileUpClustersCut(3),
fVertexBC(-200),             fRecalculateVertexBC(0),
fUseAliCentrality(0),        fCentralityClass(""),            fCentralityOpt(0),
fEventPlaneMethod(""),
fFillInputNonStandardJetBranch(kFALSE),
fNonStandardJets(new TClonesArray("AliAODJet",100)),          fInputNonStandardJetBranchName("jets"),
fFillInputBackgroundJetBranch(kFALSE), 
fBackgroundJets(0x0),fInputBackgroundJetBranchName("jets"),
fAcceptEventsWithBit(0),     fRejectEventsWithBit(0),         fRejectEMCalTriggerEventsWith2Tresholds(0),
fMomentum(),                 fOutputContainer(0x0),           fEnergyHistogramNbins(0),
fhNEventsAfterCut(0),        fNMCGenerToAccept(0),            fMCGenerEventHeaderToAccept("")
{
  for(Int_t i = 0; i < 8; i++) fhEMCALClusterCutsE [i]= 0x0 ;    
  for(Int_t i = 0; i < 7; i++) fhPHOSClusterCutsE  [i]= 0x0 ;  
  for(Int_t i = 0; i < 6; i++) fhCTSTrackCutsPt    [i]= 0x0 ;    
  for(Int_t j = 0; j < 5; j++) { fMCGenerToAccept  [j] =  ""; fMCGenerIndexToAccept[j] = -1; }
  
  InitParameters();
}

//_______________________________________
/// Destructor.
//_______________________________________
AliCaloTrackReader::~AliCaloTrackReader()
{
  DeletePointers();
}

//_______________________________________
/// Destructor. Called by the destructors  
/// of this class and derived classes.
//_______________________________________
void AliCaloTrackReader::DeletePointers()
{  
  delete fFiducialCut ;
	
  if(fAODBranchList)
  {
    fAODBranchList->Delete();
    delete fAODBranchList ;
  }
  
  if(fCTSTracks)
  {
    if(fDataType!=kMC)fCTSTracks->Clear() ;
    else              fCTSTracks->Delete() ;
    delete fCTSTracks ;
  }
  
  if(fEMCALClusters)
  {
    if(fDataType!=kMC)fEMCALClusters->Clear("C") ;
    else              fEMCALClusters->Delete() ;
    delete fEMCALClusters ;
  }
  
  if(fDCALClusters)
  {
    if(fDataType!=kMC)fDCALClusters->Clear("C") ;
    else              fDCALClusters->Delete() ;
    delete fDCALClusters ;
  }
  
  if(fPHOSClusters)
  {
    if(fDataType!=kMC)fPHOSClusters->Clear("C") ;
    else              fPHOSClusters->Delete() ;
    delete fPHOSClusters ;
  }
  
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++)
    {
      delete [] fVertex[i] ;
      
    }
    delete [] fVertex ;
  }
  
  //delete fTriggerAnalysis;
  
  if(fNonStandardJets)
  {
    if(fDataType!=kMC) fNonStandardJets->Clear("C") ;
    else               fNonStandardJets->Delete() ;
    delete fNonStandardJets ;
  }
  delete fBackgroundJets ;

  fRejectEventsWithBit.Reset();
  fAcceptEventsWithBit.Reset();
  
  if ( fWeightUtils ) delete fWeightUtils ;
    
  //  Pointers not owned, done by the analysis frame
  //  if(fInputEvent)  delete fInputEvent ;
  //  if(fOutputEvent) delete fOutputEvent ;
  //  if(fMC)          delete fMC ;
  //  Pointer not owned, deleted by maker
  //  if (fCaloUtils) delete fCaloUtils ;
}

//____________________________________________________________
/// Accept track if DCA is smaller than function.
/// \param pt of track.
/// \param dca of track.
//____________________________________________________________
Bool_t  AliCaloTrackReader::AcceptDCA(Float_t pt, Float_t dca)
{  
  Float_t cut = fTrackDCACut[0]+fTrackDCACut[1]/TMath::Power(pt,fTrackDCACut[2]);
  
  if(TMath::Abs(dca) < cut)
    return kTRUE;
  else
    return kFALSE;
}

//_____________________________________________________
/// Accept events that pass the physics selection
/// depending on an array of trigger bits set during the configuration.
//_____________________________________________________
Bool_t  AliCaloTrackReader::AcceptEventWithTriggerBit()
{  
  Int_t nAccept = fAcceptEventsWithBit.GetSize();
  
  //printf("N accept %d\n", nAccept);
  
  if( nAccept <= 0 )
    return kTRUE ; // accept the event
  
  UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  for(Int_t ibit = 0; ibit < nAccept; ibit++)
  {
    Bool_t accept = (trigFired & fAcceptEventsWithBit.At(ibit));
    
    //printf("accept %d, ibit %d, bit %d \n",accept, ibit,fAcceptEventsWithBit.At(ibit));
    if(accept) return kTRUE ; // accept the event
  }
  
  return kFALSE ; // reject the event
}

//_____________________________________________________
/// Reject particles/clusters depending on the generator 
/// of origin of the MC label.
///
/// \param mcLabel label index of particle originating the cluster or track or mc stack
//_____________________________________________________
Bool_t  AliCaloTrackReader::AcceptParticleMCLabel(Int_t mcLabel) const
{
  if( !fMC || fNMCGenerToAccept <= 0 ) return kTRUE;
  
  TString genName;
  Int_t genIndex;
  genIndex = GetCocktailGeneratorAndIndex(mcLabel, genName);
  //fMC->GetCocktailGenerator(mcLabel,genName);
  
  Bool_t generOK = kFALSE;
  for(Int_t ig = 0; ig < fNMCGenerToAccept; ig++) 
  {
    if ( fMCGenerToAccept[ig].Contains(genName) ) generOK = kTRUE;
    
    if ( generOK && fMCGenerIndexToAccept[ig] >= 0 && fMCGenerToAccept[ig] != genIndex) generOK = kFALSE;
  }
  
  if ( !generOK ) AliDebug(1, Form("skip label %d, gen %s",mcLabel,genName.Data()) );

  return generOK;
}

//_____________________________________________________________________
/// Get the name of the generator that generated a given primary particle 
/// Copy of AliMCEvent::GetCocktailGeneratorAndIndex(), modified to get the 
/// the generator index in the cocktail
///
/// \param index: mc label index
/// \param nameGen: cocktail generator name for this index
/// \return cocktail generator index
//_____________________________________________________________________
Int_t AliCaloTrackReader::GetCocktailGeneratorAndIndex(Int_t index, TString & nameGen) const
{
  //method that gives the generator for a given particle with label index (or that of the corresponding primary)
  AliVParticle* mcpart0 = (AliVParticle*) GetMC()->GetTrack(index);
  Int_t genIndex = -1;
  
  if(!mcpart0)
  {
    printf("AliMCEvent-BREAK: No valid AliMCParticle at label %i\n",index);
    return -1;
  }
  
  nameGen = GetGeneratorNameAndIndex(index,genIndex);
  
  if(nameGen.Contains("nococktailheader") ) return -1;
  
  Int_t lab=index;
  
  while(nameGen.IsWhitespace())
  {
    AliVParticle* mcpart = (AliVParticle*) GetMC()->GetTrack(lab);
    
    if(!mcpart)
    {
      printf("AliMCEvent-BREAK: No valid AliMCParticle at label %i\n",lab);
      break;
    }
    
    Int_t mother=0;
    mother = mcpart->GetMother();
    
    if(mother<0)
    {
      printf("AliMCEvent - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    
    AliVParticle* mcmom = (AliVParticle*) GetMC()->GetTrack(mother);
    if(!mcmom)
    {
      printf("AliMCEvent-BREAK: No valid AliMCParticle mother at label %i\n",mother);
      break;
    }
    
    lab=mother;
    
    nameGen = GetGeneratorNameAndIndex(mother,genIndex);
  }
  
  return genIndex;
}
//_____________________________________________________________________
/// Get the name of the generator that generated a given primary particle 
/// Copy of AliMCEvent::GetGenerator(), modified to get the 
/// the generator index in the cocktail
///
/// \param index: mc label index
/// \param genIndex: cocktail generator name for this index
/// \return cocktail generator name string
//_____________________________________________________________________
TString AliCaloTrackReader::GetGeneratorNameAndIndex(Int_t index, Int_t & genIndex) const
{
  Int_t nsumpart = GetMC()->GetNumberOfPrimaries();
  
  genIndex = -1;
  
  TList* lh = GetMC()->GetCocktailList();
  if(!lh)
  { 
    TString noheader="nococktailheader";
    return noheader;
  }
  
  Int_t nh = lh->GetEntries();
  
  for (Int_t i = nh-1; i >= 0; i--)
  {
    AliGenEventHeader* gh = (AliGenEventHeader*)lh->At(i);
    
    TString genname = gh->GetName();
    
    Int_t npart=gh->NProduced();
    
    if (i == 0) npart = nsumpart;
    
    if(index < nsumpart && index >= (nsumpart-npart)) 
    { 
      genIndex = i ;
      return genname;
    }
    
    nsumpart-=npart;
  }
  
  TString empty="";
  return empty;
}


//_____________________________________________________
/// Reject events that pass the physics selection
/// depending on an array of trigger bits set during the configuration.
//_____________________________________________________
Bool_t  AliCaloTrackReader::RejectEventWithTriggerBit()
{
  Int_t nReject = fRejectEventsWithBit.GetSize();
  
  //printf("N reject %d\n", nReject);
  
  if( nReject <= 0 )
    return kTRUE ; // accept the event
  
  UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  for(Int_t ibit = 0; ibit < nReject; ibit++)
  {
    Bool_t reject = (trigFired & fRejectEventsWithBit.At(ibit));
    
    //printf("reject %d, ibit %d, bit %d \n",reject, ibit,fRejectEventsWithBit.At(ibit));
    if(reject) return kFALSE ; // reject the event
  }
  
  return kTRUE ; // accept the event
}

//_____________________________________________
/// Do different selection of the event
/// depending on trigger name, event type, 
/// goodness of the EMCal trigger ...
//_____________________________________________
Bool_t AliCaloTrackReader::CheckEventTriggers()
{  
  //-----------------------------------------------------------
  // Reject events depending on the trigger name 
  //-----------------------------------------------------------
  
  AliDebug(1,Form("FiredTriggerClass <%s>, selected class <%s>, compare name %d",
                  GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data(),
                  GetFiredTriggerClasses().Contains(fFiredTriggerClassName)));
  
  if ( fFiredTriggerClassName != "" )
  {
    if ( !GetFiredTriggerClasses().Contains(fFiredTriggerClassName) ) 
      return kFALSE;
    else 
      AliDebug(1,"Accepted triggered event");
  }
  
  fhNEventsAfterCut->Fill(1.5);
  
  //-----------------------------------------------------------
  // Reject events depending on the event species type
  //-----------------------------------------------------------
  
  // Event types:
  //	  kStartOfRun =       1,    // START_OF_RUN
  //	  kEndOfRun =         2,    // END_OF_RUN
  //	  kStartOfRunFiles =  3,    // START_OF_RUN_FILES
  //	  kEndOfRunFiles =    4,    // END_OF_RUN_FILES
  //	  kStartOfBurst =     5,    // START_OF_BURST
  //	  kEndOfBurst =       6,    // END_OF_BURST
  //	  kPhysicsEvent =     7,    // PHYSICS_EVENT
  //	  kCalibrationEvent = 8,    // CALIBRATION_EVENT
  //	  kFormatError =      9,    // EVENT_FORMAT_ERROR
  //	  kStartOfData =      10,   // START_OF_DATA
  //	  kEndOfData =        11,   // END_OF_DATA
  //	  kSystemSoftwareTriggerEvent   = 12, // SYSTEM_SOFTWARE_TRIGGER_EVENT
  //	  kDetectorSoftwareTriggerEvent = 13  // DETECTOR_SOFTWARE_TRIGGER_EVENT

  Int_t eventType = 0;
  if(fInputEvent->GetHeader())
    eventType = ((AliVHeader*)fInputEvent->GetHeader())->GetEventType();
  
  AliDebug(3,Form("Event type %d",eventType));
  
  // Select only Physics events in data, eventType = 7, 
  // usually done by physics selection
  // MC not set, eventType =0, I think, so do not apply a value to fEventType by default
  // LED events have eventType = 8, implemented a selection cut in the past
  // but not really useful for EMCal data since LED are not reconstructed (unless wrongly assumed as physics)
  
  if ( fEventType >= 0 && eventType != fEventType ) return kFALSE ;
  
  AliDebug(1,"Pass event species selection");

  fhNEventsAfterCut->Fill(2.5);
  
  //-----------------------------------------------------------------
  // In case of mixing analysis, select here the trigger of the event
  //-----------------------------------------------------------------
  
  UInt_t isTrigger = kFALSE;
  UInt_t isMB      = kFALSE;
  
  if(!fEventTriggerAtSE)
  {
    // In case of mixing analysis, accept MB events, not only Trigger
    // Track and cluster arrays filled for MB in order to create the pool in the corresponding analysis
    // via de method in the base class FillMixedEventPool()
    
    AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
    
    if(!inputHandler) return kFALSE ;  // to content coverity
    
    isTrigger = inputHandler->IsEventSelected() & fEventTriggerMask;
    isMB      = inputHandler->IsEventSelected() & fMixEventTriggerMask;
    
    if(!isTrigger && !isMB) return kFALSE;
    
    //printf("Selected triggered event : %s\n",GetFiredTriggerClasses().Data());
    AliDebug(0,"Pass uninteresting triggered events rejection in case of mixing analysis");  
    
    fhNEventsAfterCut->Fill(3.5);
  }

  
  //-------------------------------------------------------------------------------------
  // Reject or accept events depending on the trigger bit
  //-------------------------------------------------------------------------------------
  
  Bool_t okA = AcceptEventWithTriggerBit();
  Bool_t okR = RejectEventWithTriggerBit();
  
  //printf("AliCaloTrackReader::FillInputEvent() - Accept event? %d, Reject event %d? \n",okA,okR);
  
  if(!okA || !okR) return kFALSE;
  
  AliDebug(1,"Pass event bit rejection");
  
  fhNEventsAfterCut->Fill(4.5);
  
  //----------------------------------------------------------------------
  // Do not count events that were likely triggered by an exotic cluster
  // or out BC cluster
  //----------------------------------------------------------------------
  
  // Set a bit with the event kind, MB, L0, L1 ...
  SetEventTriggerBit();
  
  // In case of Mixing, avoid checking the triggers in the min bias events
  if(!fEventTriggerAtSE && (isMB && !isTrigger)) return kTRUE;
  
  if( (IsEventEMCALL1() || IsEventEMCALL0())  &&  fTriggerPatchClusterMatch)
  {
    if(fRejectEMCalTriggerEventsWith2Tresholds)
    {
      // Reject triggered events when there is coincidence on both EMCal trigger thresholds,
      // but the requested trigger is the low trigger threshold
      if(IsEventEMCALL1Jet1  () && IsEventEMCALL1Jet2  () && fFiredTriggerClassName.Contains("EJ2")) return kFALSE;
      if(IsEventEMCALL1Gamma1() && IsEventEMCALL1Gamma2() && fFiredTriggerClassName.Contains("EG2")) return kFALSE;
    }
    
    //Get Patches that triggered
    TArrayI patches = GetTriggerPatches(fTriggerPatchTimeWindow[0],fTriggerPatchTimeWindow[1]);
    
    MatchTriggerCluster(patches);
    
    patches.Reset();
    
    // If requested, remove badly triggeed events, but only when the EMCal trigger bit is set
    if(fRemoveBadTriggerEvents)
    {
     AliDebug(1,Form("ACCEPT triggered event? \n exotic? %d - bad cell %d - bad Max cell %d - BC %d  - Matched %d\n",
                     fIsExoticEvent,fIsBadCellEvent, fIsBadMaxCellEvent, fTriggerClusterBC,fIsTriggerMatch));
      if     (fIsExoticEvent)         return kFALSE;
      else if(fIsBadCellEvent)        return kFALSE;
      else if(fRemoveUnMatchedTriggers && !fIsTriggerMatch) return kFALSE ;
      else if(fTriggerClusterBC != 0) return kFALSE;
      AliDebug(1,Form("\t *** YES for %s",GetFiredTriggerClasses().Data()));
    }
    
    AliDebug(1,"Pass EMCal triggered event rejection \n"); 
    
    fhNEventsAfterCut->Fill(5.5);
  }
  
  //-------------------------------------------------------------------------------------
  // Select events only fired by a certain trigger configuration if it is provided
  //-------------------------------------------------------------------------------------
  
  if (GetFiredTriggerClasses().Contains("FAST")  && !GetFiredTriggerClasses().Contains("ALL") && !fAcceptFastCluster)
  {
    AliDebug(1,Form("Do not count events from fast cluster, trigger name %s\n",fFiredTriggerClassName.Data()));
    return kFALSE;
    
    fhNEventsAfterCut->Fill(6.5);
  }
  
  //-------------------------------------------------------------------------------------
  // Reject event if large clusters with large energy
  // Use only for LHC11a data for the moment, and if input is clusterizer V1 or V1+unfolding
  // If clusterzer NxN or V2 it does not help
  //-------------------------------------------------------------------------------------
  
  //Int_t run = fInputEvent->GetRunNumber();
  
  if ( fRemoveLEDEvents )
  {
    Bool_t reject = RejectLEDEvents();

    if(reject) return kFALSE;
    
    AliDebug(1,"Pass LED event rejection");
    
    fhNEventsAfterCut->Fill(7.5);
  } // Remove LED events

  // All selection criteria passed, accept the event
  return kTRUE ;
}

//________________________________________________
/// Check the MC PYTHIA event, if the requested 
/// pT-hard is much smaller than the jet pT, then,
/// there can be a problem in the tails of the 
/// distributions and the event should be rejected.
//________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndJetPt()
{  
  //printf("AliCaloTrackReader::ComparePtHardAndJetPt() - GenHeaderName : %s\n",GetGenEventHeader()->ClassName());
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    TParticle * jet =  0;
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Int_t nTriggerJets =  pygeh->NTriggerJets();
    Float_t ptHard = pygeh->GetPtHard();
    
    AliDebug(1,Form("Njets: %d, pT Hard %f",nTriggerJets, ptHard));
    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      pygeh->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      
      AliDebug(1,Form("jet %d; pycell jet pT %f",ijet, jet->Pt()));
      
      //Compare jet pT and pt Hard
      if(jet->Pt() > fPtHardAndJetPtFactor * ptHard)
      {
        AliInfo(Form("Reject jet event with : pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n",
                     ptHard, jet->Pt(), fPtHardAndJetPtFactor));
        return kFALSE;
      }
    }
    
    if(jet) delete jet;
  }
  
  return kTRUE ;
}

//____________________________________________________
/// Check the MC PYTHIA event, if the requested 
/// pT-hard is smaller than the calorimeter cluster E, 
/// there can be a problem in the tails of the 
/// distributions and the event should be rejected.
//____________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndClusterPt()
{  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Float_t ptHard = pygeh->GetPtHard();
    
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = fInputEvent->GetCaloCluster(iclus) ;
      Float_t ecluster = clus->E();
      
      if(ecluster > fPtHardAndClusterPtFactor*ptHard)
      {
        AliInfo(Form("Reject : ecluster %2.2f, calo %d, factor %2.2f, ptHard %f",ecluster,clus->GetType(),fPtHardAndClusterPtFactor,ptHard));
        
        return kFALSE;
      }
    }
    
  }
  
  return kTRUE ;
}

//___________________________________________________
/// Fill the output list of initialized control histograms.
/// Cluster or track spectra histograms, depending on different selection cuts.
//___________________________________________________
TList * AliCaloTrackReader::GetCreateControlHistograms()
{  
  
  fhNEventsAfterCut = new TH1I("hNEventsAfterCut", "Number of analyzed events", 19, 0, 19) ;
  //fhNEventsAfterCut->SetXTitle("Selection");
  fhNEventsAfterCut->SetYTitle("# events");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(1 ,"1=Input");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(2 ,"2=Trigger string");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(3 ,"3=Event Type");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(4 ,"4=Mixing Event");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(5 ,"5=Trigger Bit");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(6 ,"6=Good EMC Trigger");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(7 ,"7=!Fast Cluster");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(8 ,"8=!LED");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(9 ,"9=Time stamp"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(10,"10=Z vertex"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(11,"11=Primary vertex"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(12,"12=Pile-up"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(13,"13=V0AND"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(14,"14=Centrality"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(15,"15=GenHeader"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(16,"16=PtHard-Jet");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(17,"17=PtHard-Cluster"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(18,"18=Track multi."); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(19,"19=TOF BC"); 
  fOutputContainer->Add(fhNEventsAfterCut);

  
  if(fFillEMCAL)
  {
    for(Int_t i = 0; i < 8; i++)
    {
      TString names[] = {"NoCut", "Corrected", "GoodCluster", "NonLinearity", 
        "EnergyAndFidutial", "Time", "NCells", "BadDist"};
      
      fhEMCALClusterCutsE[i] = new TH1F(Form("hEMCALReaderClusterCuts_%d_%s",i,names[i].Data()),
                                        Form("EMCal %d, %s",i,names[i].Data()),   
                                        fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]);
      fhEMCALClusterCutsE[i]->SetYTitle("# clusters");
      fhEMCALClusterCutsE[i]->SetXTitle("#it{E} (GeV)");
      fOutputContainer->Add(fhEMCALClusterCutsE[i]);
    }
  }
  
  if(fFillPHOS)
  {
    for(Int_t i = 0; i < 7; i++)
    {
      TString names[] = {"NoCut", "ExcludeCPV", "BorderCut", "FiducialCut", "EnergyCut", "NCells", "BadDist"};
      
      fhPHOSClusterCutsE[i] = new TH1F(Form("hPHOSReaderClusterCuts_%d_%s",i,names[i].Data()),
                                       Form("PHOS Cut %d, %s",i,names[i].Data()), 
                                       fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
      fhPHOSClusterCutsE[i]->SetYTitle("# clusters");
      fhPHOSClusterCutsE[i]->SetXTitle("#it{E} (GeV)");
      fOutputContainer->Add(fhPHOSClusterCutsE[i]);
    }
  }
  
  if(fFillCTS)
  {
    for(Int_t i = 0; i < 6; i++)
    {
      TString names[] = {"NoCut", "Status", "ESD_AOD", "TOF", "DCA","PtAcceptanceMult"};

      fhCTSTrackCutsPt[i] = new TH1F(Form("hCTSReaderClusterCuts_%d_%s",i,names[i].Data()),
                                     Form("CTS Cut %d, %s",i,names[i].Data()), 
                                     fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
      fhCTSTrackCutsPt[i]->SetYTitle("# tracks");
      fhCTSTrackCutsPt[i]->SetXTitle("#it{p}_{T} (GeV)");
      fOutputContainer->Add(fhCTSTrackCutsPt[i]);
    }
  }
  
  return fOutputContainer ;
}

//_____________________________________________________
/// Save parameters used for analysis in a string.
//_____________________________________________________
TObjString *  AliCaloTrackReader::GetListOfParameters()
{
  TString parList ; //this will be list of parameters used for this analysis.

  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- Reader ---:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Data type: %d; zvertex cut %2.2f; EMC cluster name: <%s> ",fDataType, fZvtxCut, fEMCALClustersListName.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Use detector: EMC %d, DCA %d, PHOS %d, CTS %d, EMCcells %d, PHOScells %d ; ",
           fFillEMCAL,fFillDCAL,fFillPHOS,fFillCTS,fFillEMCALCells,fFillPHOSCells) ;
  snprintf(onePar,buffersize,"E-pT window: EMC (%2.1f,%2.1f), PHOS (%2.1f,%2.1f), CTS (%2.1f,%2.1f); ",
           fEMCALPtMin,fEMCALPtMax,fPHOSPtMin,fPHOSPtMax,fCTSPtMin,fCTSPtMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Dist to bad channel: EMC > %2.1f, PHOS > %2.1f; ",fEMCALBadChMinDist,fPHOSBadChMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"N cells: EMC > %d, PHOS > %d; ",fEMCALNCellsCut,fPHOSNCellsCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"EMC time cut single window (%2.2f,%2.2f); ",fEMCALTimeCutMin,fEMCALTimeCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Check: calo fid cut %d; ",fCheckFidCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Track: status %d, SPD hit %d; ",(Int_t) fTrackStatus, fSelectSPDHitTracks) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"multip. eta cut %1.1f; npt cuts %d;",fTrackMultEtaCut, fTrackMultNPtCut) ;
  parList+=onePar ;
  
  if(fUseTrackDCACut)
  {
    snprintf(onePar,buffersize,"DCA cut ON, param (%2.4f,%2.4f,%2.4f); ",fTrackDCACut[0],fTrackDCACut[1],fTrackDCACut[2]) ;
    parList+=onePar ;
  }
  
  snprintf(onePar,buffersize,"Recalculate Clusters = %d, E linearity = %d; ",fRecalculateClusters, fCorrectELinearity) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"SE trigger sel. %d, not? trigger Mask? %d, MB Trigger Mask for mixed %d; ",
           fEventTriggerAtSE, fEventTriggerMask,fMixEventTriggerMask);
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"Select fired trigger %s; Remove Bad trigger event %d, unmatched %d; Accept fastcluster %d, Reject LED %d ",
          fFiredTriggerClassName.Data(), fRemoveBadTriggerEvents, fRemoveUnMatchedTriggers, fAcceptFastCluster, fRemoveLEDEvents);
  parList+=onePar ;
  
  if(fNMCGenerToAccept)
  {
    snprintf(onePar,buffersize,"Accept only labels from: ");
    parList+=onePar ;
    for(Int_t i = 0; i< fNMCGenerToAccept; i++)
      parList+=(fMCGenerToAccept[i]+" ") ;
         
    snprintf(onePar,buffersize,"; ");
    parList+=onePar ;
  }
  
  if(fSmearShowerShape)
  {
    snprintf(onePar,buffersize,"EMC M02 smear ON, function %d, param %2.4f, %d<=NLM<=%d; ",
             fSmearingFunction,fSmearShowerShapeWidth,fSmearNLMMin,fSmearNLMMax) ;
    parList+=onePar ;
  }
  
  if(fComparePtHardAndJetPt)
  {
    snprintf(onePar,buffersize,"jet pt / pt hard < %2.1f; ",fPtHardAndJetPtFactor);
    parList+=onePar ;
  }
  
  if(fComparePtHardAndClusterPt)
  {
    snprintf(onePar,buffersize,"cluster pt / pt hard < %2.2f",fPtHardAndClusterPtFactor);
    parList+=onePar ;
  }
  
  snprintf(onePar,buffersize,"Centrality: Class %s, Option %d, Bin [%d,%d]; New centrality %d; Event plane method %s; ", 
           fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1],fUseAliCentrality,fEventPlaneMethod.Data()) ;
  parList+=onePar ;
  
  return new TObjString(parList) ;
}


//____________________________________________
/// \return pointer to stack (AliStack)
//____________________________________________
AliStack* AliCaloTrackReader::GetStack() const
{
  if(fMC)
    return fMC->Stack();
  else
  {
    AliDebug(1,"Stack is not available");
    return 0x0 ;
  }
}

//______________________________________________
/// \return pointer to header (AliHeader)
//______________________________________________
AliHeader* AliCaloTrackReader::GetHeader() const
{
  if(fMC)
  {
    return fMC->Header();
  }
  else
  {
    AliInfo("Header is not available");
    return 0x0 ;
  }
}

//______________________________________________________________
/// \return pointer to Generated event header (AliGenEventHeader)
//______________________________________________________________
AliGenEventHeader* AliCaloTrackReader::GetGenEventHeader() const
{
  if     (ReadStack() && fMC)
  {
    AliGenEventHeader * eventHeader = fMC->GenEventHeader();
    
    if(fMCGenerEventHeaderToAccept=="") return eventHeader ;
        
    AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    
    if(!cocktail) return 0x0 ;
    
    TList *genHeaders = cocktail->GetHeaders();
    
    Int_t nGenerators = genHeaders->GetEntries();
    
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      //printf("ESD Event header %d %s\n",igen,name.Data());

      if(name.Contains(fMCGenerEventHeaderToAccept,TString::kIgnoreCase)) 
        return eventHeader2 ;
    }

    return 0x0;
    
  }
  else if(ReadAODMCParticles() && GetAODMCHeader())
  {
    Int_t nGenerators = GetAODMCHeader()->GetNCocktailHeaders();

    if( nGenerators <= 0)        return 0x0;
    
    if(fMCGenerEventHeaderToAccept=="") return GetAODMCHeader()->GetCocktailHeader(0);
        
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader = GetAODMCHeader()->GetCocktailHeader(igen) ;
      TString name = eventHeader->GetName();
      //printf("AOD Event header %d %s\n",igen,name.Data());
      
      if(name.Contains(fMCGenerEventHeaderToAccept,TString::kIgnoreCase))
        return eventHeader ;
    }
    
    return 0x0;
        
  }
  else
  {
    return 0x0;
  }
}

//_________________________________________________________
/// \return list of particles in AOD, 
/// Implemented in AliCaloTrackAODReader.
//_________________________________________________________
TClonesArray* AliCaloTrackReader::GetAODMCParticles() const
{  
  AliInfo("Input are not AODs");
  
  return NULL ;
}

//________________________________________________________
/// \return MC header in AOD. 
/// Implemented in AliCaloTrackAODReader.
//________________________________________________________
AliAODMCHeader* AliCaloTrackReader::GetAODMCHeader() const
{  
  AliInfo("Input are not AODs");
  
  return NULL ;
}

//___________________________________________________________
/// \return vertex Bunch Cross Number.
/// In old AODs BC not stored, recalculate it here,
/// loop over the global track and select those which have 
/// small DCA to primary vertex (e.g. primary).
/// If at least one of these primaries has valid BC != 0, then 
/// this vertex is a pile-up candidate.
/// Execute after CTS filtering.
//___________________________________________________________
Int_t AliCaloTrackReader::GetVertexBC(const AliVVertex * vtx)
{  
  Int_t vertexBC=vtx->GetBC();
  
  if(!fRecalculateVertexBC) return vertexBC;
  
  // Value not available, recalculate it.
  
  Double_t bz  = fInputEvent->GetMagneticField();
  Bool_t   bc0 = kFALSE;
  Int_t    ntr = GetCTSTracks()->GetEntriesFast();
  //printf("N Tracks %d\n",ntr);
  
  for(Int_t i = 0 ; i < ntr ; i++)
  {
    AliVTrack * track =  (AliVTrack*) (GetCTSTracks()->At(i));
    
    //Check if has TOF info, if not skip
    ULong_t status  = track->GetStatus();
    Bool_t  okTOF   = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
    vertexBC        = track->GetTOFBunchCrossing(bz);
    Float_t pt      = track->Pt();
    
    if(!okTOF) continue;
    
    // Get DCA x, y
    Double_t dca[2]   = {1e6,1e6};
    Double_t covar[3] = {1e6,1e6,1e6};
    track->PropagateToDCA(vtx,bz,100.,dca,covar);
    
    if(AcceptDCA(pt,dca[0]))
    {
      if     (vertexBC !=0 && fVertexBC != AliVTrack::kTOFBCNA) return vertexBC;
      else if(vertexBC == 0)                                    bc0 = kTRUE;
    }
  }
  
  if( bc0 ) vertexBC = 0 ;
  else      vertexBC = AliVTrack::kTOFBCNA ;
  
  return vertexBC;
}

//_____________________________
/// 
/// Return track ID, different for ESD and AODs.
/// See AliCaloTrackAODReader for AOD correspondent
///
/// \return track ID
/// \param track: pointer to track
//_____________________________
Int_t AliCaloTrackReader::GetTrackID(AliVTrack* track) 
{
  return track->GetID();
}

//_____________________________
/// Init the reader. 
/// Method to be called in AliAnaCaloTrackCorrMaker.
//_____________________________
void AliCaloTrackReader::Init()
{  
  if(fReadStack && fReadAODMCParticles)
  {
    AliInfo("Cannot access stack and mcparticles at the same time, change them");
    fReadStack          = kFALSE;
    fReadAODMCParticles = kFALSE;
  }

  // Activate debug level in AliAnaWeights
  if( fWeightUtils->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(fWeightUtils->ClassName(), fWeightUtils->GetDebug());
}

//_______________________________________
/// Initialize the parameters with default.
//_______________________________________
void AliCaloTrackReader::InitParameters()
{
  fDataType   = kESD ;
  fCTSPtMin   = 0.1 ;
  fEMCALPtMin = 0.1 ;
  fPHOSPtMin  = 0.1 ;
  fCTSPtMax   = 1000. ;
  fEMCALPtMax = 1000. ;
  fPHOSPtMax  = 1000. ;
  
  fEMCALBadChMinDist = 0; // open, 2; // standard       
  fPHOSBadChMinDist  = 0; // open, 2; // standard   
  
  fEMCALNCellsCut    = 0; // open, 1; // standard          
  fPHOSNCellsCut     = 0; // open, 2; // standard
  
  //Track DCA cuts
  // dca_xy cut = 0.0105+0.0350/TMath::Power(pt,1.1);
  fTrackDCACut[0] = 0.0105;
  fTrackDCACut[1] = 0.0350;
  fTrackDCACut[2] = 1.1;
  
  //Do not filter the detectors input by default.
  fFillEMCAL      = kFALSE;
  fFillDCAL       = kFALSE;
  fFillPHOS       = kFALSE;
  fFillCTS        = kFALSE;
  fFillEMCALCells = kFALSE;
  fFillPHOSCells  = kFALSE;
  
  fReadStack             = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fReadAODMCParticles    = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fDeltaAODFileName      = "deltaAODPartCorr.root";
  fFiredTriggerClassName = "";
  fEventTriggerMask      = AliVEvent::kAny;
  fMixEventTriggerMask   = AliVEvent::kAnyINT;
  fEventTriggerAtSE      = kTRUE; // Use only events that pass event selection at SE base class
  
  fAcceptFastCluster = kTRUE;
  fEventType         = -1;
  
  //We want tracks fitted in the detectors:
  //fTrackStatus=AliESDtrack::kTPCrefit;
  //fTrackStatus|=AliESDtrack::kITSrefit;
  fTrackStatus     = 0;
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0;
  fV0Mul[0] = 0;   fV0Mul[1] = 0;
  
  fZvtxCut   = 10.;
  
  fNMixedEvent = 1;
  
  fPtHardAndJetPtFactor     = 7.;
  fPtHardAndClusterPtFactor = 1.;
  
  //Centrality
  fUseAliCentrality = kFALSE;
  fCentralityClass  = "V0M";
  fCentralityOpt    = 100;
  fCentralityBin[0] = fCentralityBin[1]=-1;
  
  fEventPlaneMethod = "V0";
  
  // Allocate memory (not sure this is the right place)
  fCTSTracks       = new TObjArray();
  fEMCALClusters   = new TObjArray();
  fDCALClusters    = new TObjArray();
  fPHOSClusters    = new TObjArray();
  //fTriggerAnalysis = new AliTriggerAnalysis;
  fAODBranchList   = new TList ;
  fOutputContainer = new TList ;
  
  fEnergyHistogramNbins    = 200;
  fEnergyHistogramLimit[0] = 0  ;
  fEnergyHistogramLimit[1] = 100;
  
  fPileUpParamSPD[0] = 3   ; fPileUpParamSPD[1] = 0.8 ;
  fPileUpParamSPD[2] = 3.0 ; fPileUpParamSPD[3] = 2.0 ; fPileUpParamSPD[4] = 5.0;
  
  // Parametrized time cut (LHC11d)
  fEMCALParamTimeCutMin[0] =-5; fEMCALParamTimeCutMin[1] =-1 ; fEMCALParamTimeCutMin[2] = 3.5 ; fEMCALParamTimeCutMin[3] = 1.  ;
  fEMCALParamTimeCutMax[0] = 5; fEMCALParamTimeCutMax[1] = 50; fEMCALParamTimeCutMax[2] = 0.45; fEMCALParamTimeCutMax[3] = 1.25;
  
  // Parametrized time cut (LHC11c)
  //fEMCALParamTimeCutMin[0] =-5;   fEMCALParamTimeCutMin[1] =-1 ; fEMCALParamTimeCutMin[2] = 1.87; fEMCALParamTimeCutMin[3] = 0.4;
  //fEMCALParamTimeCutMax[0] = 3.5; fEMCALParamTimeCutMax[1] = 50; fEMCALParamTimeCutMax[2] = 0.15; fEMCALParamTimeCutMax[3] = 1.6;
  
  fTimeStampRunMin = -1;
  fTimeStampRunMax = 1e12;
  fTimeStampEventFracMin = -1;
  fTimeStampEventFracMax = 2;
  
  for(Int_t i = 0; i < 19; i++)
  {
    fEMCalBCEvent   [i] = 0;
    fEMCalBCEventCut[i] = 0;
    fTrackBCEvent   [i] = 0;
    fTrackBCEventCut[i] = 0;
  }
  
  // Trigger match-rejection
  fTriggerPatchTimeWindow[0] = 8;
  fTriggerPatchTimeWindow[1] = 9;
  
  fTriggerClusterBC        = -10000 ;
  fTriggerL0EventThreshold = -1;
  fTriggerL1EventThreshold = -1;
  fTriggerClusterIndex     = -1;
  fTriggerClusterId        = -1;
  
  //Jets
  fInputNonStandardJetBranchName = "jets";
  fFillInputNonStandardJetBranch = kFALSE;
  if(!fNonStandardJets) fNonStandardJets = new TClonesArray("AliAODJet",100);
  fInputBackgroundJetBranchName = "jets";
  fFillInputBackgroundJetBranch = kFALSE; 
  if(!fBackgroundJets) fBackgroundJets = new AliAODJetEventBackground();

  fSmearShowerShapeWidth = 0.005;
  fSmearNLMMin = 1;
  fSmearNLMMax = 1;
  
  fWeightUtils = new AliAnaWeights() ;
  fEventWeight = 1 ;
    
  fTrackMultNPtCut = 8;
  fTrackMultPtCut[0] = 0.15; fTrackMultPtCut[1] = 0.5;  fTrackMultPtCut[2] = 1.0; 
  fTrackMultPtCut[3] = 2.0 ; fTrackMultPtCut[4] = 4.0;  fTrackMultPtCut[5] = 6.0;  
  fTrackMultPtCut[6] = 8.0 ; fTrackMultPtCut[7] = 10.;  
  fTrackMultPtCut[8] = 15.0; fTrackMultPtCut[9] = 20.;  
}

//__________________________________________________________________________
/// Select the cluster depending on a time window, either a simple
/// range or a parametrization depending on the energy.
//__________________________________________________________________________
Bool_t AliCaloTrackReader::IsInTimeWindow(Double_t tof, Float_t energy) const
{  
  // Parametrized cut depending on E
  if(fUseParamTimeCut)
  {
    Float_t minCut= fEMCALParamTimeCutMin[0]+fEMCALParamTimeCutMin[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMin[2])/fEMCALParamTimeCutMin[3]);
    Float_t maxCut= fEMCALParamTimeCutMax[0]+fEMCALParamTimeCutMax[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMax[2])/fEMCALParamTimeCutMax[3]);
    //printf("tof %f, minCut %f, maxCut %f\n",tof,minCut,maxCut);
    if( tof < minCut || tof > maxCut )  return kFALSE ;
  }
  
  //In any case, the time should to be larger than the fixed window ...
  if( tof < fEMCALTimeCutMin  || tof > fEMCALTimeCutMax )  return kFALSE ;
  
  return kTRUE ;
}

//________________________________________________
/// Check if event is from pile-up determined by SPD.
/// Default values: (3, 0.8, 3., 2., 5.).
//________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPD() const
{
  return fInputEvent->IsPileupFromSPD((Int_t) fPileUpParamSPD[0] , fPileUpParamSPD[1] ,
                                      fPileUpParamSPD[2] , fPileUpParamSPD[3] , fPileUpParamSPD[4] );
  //printf("Param : %d, %2.2f, %2.2f, %2.2f, %2.2f\n",(Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
}

//__________________________________________________
/// Check if event is from pile-up determined by EMCal
//__________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCal() const
{
  if(fNPileUpClusters > fNPileUpClustersCut) return kTRUE ;
  else                                       return kFALSE;
}

//________________________________________________________
/// Check if event is from pile-up determined by SPD and EMCal.
//________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndEMCal() const
{
  if( IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//_______________________________________________________
/// Check if event is from pile-up determined by SPD or EMCal.
//_______________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDOrEMCal() const
{
  if( IsPileUpFromSPD() || IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
/// Check if event is from pile-up determined by SPD and not by EMCal.
//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndNotEMCal() const
{
  if( IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
/// Check if event is from pile-up determined by EMCal, not by SPD.
//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCalAndNotSPD() const
{
  if( !IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                           return kFALSE;
}

//______________________________________________________________
/// Check if event not from pile-up determined neither by SPD nor by EMCal.
//______________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromNotSPDAndNotEMCal() const
{
  if( !IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                            return kFALSE;
}

//___________________________________________________________________________________
/// Event and tracks/clusters filtering method. Main steps:
/// * Accept/reject the event looking to the triggers, vertex, pile-up, time stamps, 
/// PYTHIA pT hard, centrality etc.
/// * Filter the tracks and calorimeter clusters, even correct the clusters if requested
/// and put them in lists.
///
/// Called by the analysis maker.
//___________________________________________________________________________________
Bool_t AliCaloTrackReader::FillInputEvent(Int_t iEntry, const char * /*curFileName*/)
{  
  fEventNumber         = iEntry;
  fTriggerClusterIndex = -1;
  fTriggerClusterId    = -1;
  fIsTriggerMatch      = kFALSE;
  fTriggerClusterBC    = -10000;
  fIsExoticEvent       = kFALSE;
  fIsBadCellEvent      = kFALSE;
  fIsBadMaxCellEvent   = kFALSE;
  
  fIsTriggerMatchOpenCut[0] = kFALSE ;
  fIsTriggerMatchOpenCut[1] = kFALSE ;
  fIsTriggerMatchOpenCut[2] = kFALSE ;
  
  //fCurrentFileName = TString(currentFileName);
  if(!fInputEvent)
  {
    AliInfo("Input event not available, skip event analysis");
    return kFALSE;
  }
  
  fhNEventsAfterCut->Fill(0.5);
  
  //-----------------------------------------------
  // Select the event depending on the trigger type
  // and other event characteristics
  // like the goodness of the EMCal trigger
  //-----------------------------------------------
  
  Bool_t accept = CheckEventTriggers();
  if(!accept) return kFALSE;
  
  AliDebug(1,"Pass Event trigger selection");
  

  //------------------------------------------------------
  // Event rejection depending on time stamp
  //------------------------------------------------------
  
  if(fDataType==kESD && fTimeStampEventSelect)
  {
    AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
    if(esd)
    {
      Int_t timeStamp = esd->GetTimeStamp();
      Float_t timeStampFrac = 1.*(timeStamp-fTimeStampRunMin) / (fTimeStampRunMax-fTimeStampRunMin);
      
      //printf("stamp0 %d, max0 %d, frac %f\n", timeStamp-fTimeStampRunMin,fTimeStampRunMax-fTimeStampRunMin, timeStampFrac);
      
      if(timeStampFrac < fTimeStampEventFracMin || timeStampFrac > fTimeStampEventFracMax) return kFALSE;
    }
    
    AliDebug(1,"Pass Time Stamp rejection");
    fhNEventsAfterCut->Fill(8.5);
  }

  //------------------------------------------------------
  // Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------
  
  FillVertexArray();
  
  //Reject events with Z vertex too large, only for SE analysis, if not, cut on the analysis code
  if(!GetMixedEvent() && TMath::Abs(fVertex[0][2]) > fZvtxCut) return kFALSE;
  
  fhNEventsAfterCut->Fill(9.5);

  if(fUseEventsWithPrimaryVertex)
  {
    if( !CheckForPrimaryVertex() )              return kFALSE; // algorithm in ESD/AOD Readers
    if( TMath::Abs(fVertex[0][0] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][1] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][2] ) < 1.e-6    ) return kFALSE;
  }
  
  AliDebug(1,"Pass primary vertex rejection");
  
  fhNEventsAfterCut->Fill(10.5);

  //printf("Reader : IsPileUp %d, Multi %d\n",IsPileUpFromSPD(),fInputEvent->IsPileupFromSPDInMultBins());
  
  if(fDoPileUpEventRejection)
  {
    // Do not analyze events with pileup
    Bool_t bPileup = IsPileUpFromSPD();
    //IsPileupFromSPDInMultBins() // method to try
    //printf("pile-up %d, %d, %2.2f, %2.2f, %2.2f, %2.2f\n",bPileup, (Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
    if(bPileup) return kFALSE;
    
    AliDebug(1,"Pass Pile-Up event rejection");
    
    fhNEventsAfterCut->Fill(11.5);
  }
  
  if(fDoV0ANDEventSelection)
  {
    AliVVZERO* v0 = fInputEvent->GetVZEROData();

    Bool_t bV0AND = ((v0->GetV0ADecision()==1) && (v0->GetV0CDecision()==1));
    //bV0AND = fTriggerAnalysis->IsOfflineTriggerFired((AliESDEvent*)fInputEvent, AliTriggerAnalysis::kV0AND);
    //printf("V0AND event? %d\n",bV0AND);

    if(!bV0AND)
    {
      AliDebug(1,"Reject event by V0AND");
      return kFALSE;
    }
    
    AliDebug(1,"Pass V0AND event rejection");
    
    fhNEventsAfterCut->Fill(12.5);
  }

  //------------------------------------------------------
  // Check if there is a centrality value, PbPb analysis,
  // and if a centrality bin selection is requested
  //------------------------------------------------------
  
  //If we need a centrality bin, we select only those events in the corresponding bin.
  Int_t cen = -1;
  if ( fCentralityBin[0] >= 0 && fCentralityBin[1] >= 0 )
  {
    cen = GetEventCentrality();
      
    if(cen > fCentralityBin[1] || cen <= fCentralityBin[0]) return kFALSE; //reject events out of bin.
    
    AliDebug(1,"Pass centrality rejection");
    
    fhNEventsAfterCut->Fill(13.5);
  }
  
  //----------------------------------------------------------------
  // Reject the event if the event header name is not
  // the one requested among the possible generators.
  // Needed in case of cocktail MC generation with multiple options.
  //----------------------------------------------------------------
  if(fMCGenerEventHeaderToAccept!="" && !GetGenEventHeader()) 
    return kFALSE;
  
  fhNEventsAfterCut->Fill(14.5);

  //---------------------------------------------------------------------------
  // In case of analysis of events with jets, skip those with jet pt > 5 pt hard
  // To be used on for MC data in pT hard bins
  //---------------------------------------------------------------------------
  
  if(fComparePtHardAndJetPt)
  {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
    AliDebug(1,"Pass Pt Hard - Jet rejection");
    fhNEventsAfterCut->Fill(15.5);
  }
  
  if(fComparePtHardAndClusterPt)
  {
    if(!ComparePtHardAndClusterPt()) return kFALSE ;
    AliDebug(1,"Pass Pt Hard - Cluster rejection");
    fhNEventsAfterCut->Fill(16.5);
  }
  
  //------------------------------------------------------------------
  // Recover the weight assigned to the event, if provided
  // right now only for pT-hard bins and centrality depedent weights
  //------------------------------------------------------------------
  if ( fWeightUtils->IsWeightSettingOn() )
  {
    fWeightUtils->SetCentrality(cen);
    
    fWeightUtils->SetPythiaEventHeader(((AliGenPythiaEventHeader*)GetGenEventHeader()));
      
    fEventWeight = fWeightUtils->GetWeight();
  }
  
  //-------------------------------------------------------
  // Get the main vertex BC, in case not available
  // it is calculated in FillCTS checking the BC of tracks
  //------------------------------------------------------
  fVertexBC = fInputEvent->GetPrimaryVertex()->GetBC();
  
  //-----------------------------------------------
  // Fill the arrays with cluster/tracks/cells data
  //-----------------------------------------------
  
  if(fFillCTS)
  {
    FillInputCTS();
    //Accept events with at least one track
    if(fTrackMult[0] == 0 && fDoRejectNoTrackEvents) return kFALSE ;
    
    fhNEventsAfterCut->Fill(17.5);
    
    AliDebug(1,"Pass rejection of null track events");
  }
  
  if(fDoVertexBCEventSelection)
  {
    if(fVertexBC != 0 && fVertexBC != AliVTrack::kTOFBCNA) return kFALSE ;
    
    AliDebug(1,"Pass rejection of events with vertex at BC!=0");
    
    fhNEventsAfterCut->Fill(18.5);
  }
  
  if(fFillEMCALCells)
    FillInputEMCALCells();
  
  if(fFillPHOSCells)
    FillInputPHOSCells();
  
  if(fFillEMCAL || fFillDCAL)
    FillInputEMCAL();
  
  if(fFillPHOS)
    FillInputPHOS();
  
  FillInputVZERO();
  
  //one specified jet branch
  if(fFillInputNonStandardJetBranch)
    FillInputNonStandardJets();
  if(fFillInputBackgroundJetBranch)
    FillInputBackgroundJets();

  AliDebug(1,"Event accepted for analysis");

  return kTRUE ;
}

//__________________________________________________
/// \return Current event centrality bin. 
/// Different percentile options and centrality class can be requested.
//__________________________________________________
Int_t AliCaloTrackReader::GetEventCentrality() const
{  
  if(fUseAliCentrality)
  {
    if ( !GetCentrality() ) return -1;
    
    AliDebug(1,Form("Cent. Percentile: V0M %2.2f, CL0 %2.2f, CL1 %2.2f; selected class %s", 
                    GetCentrality()->GetCentralityPercentile("V0M"), 
                    GetCentrality()->GetCentralityPercentile("CL0"), 
                    GetCentrality()->GetCentralityPercentile("CL1"), 
                    fCentralityClass.Data()));
    
    if     (fCentralityOpt == 100) return (Int_t) GetCentrality()->GetCentralityPercentile(fCentralityClass); // 100 bins max
    else if(fCentralityOpt ==  10) return GetCentrality()->GetCentralityClass10(fCentralityClass);// 10 bins max
    else if(fCentralityOpt ==  20) return GetCentrality()->GetCentralityClass5(fCentralityClass); // 20 bins max
    else
    {
      AliInfo(Form("Unknown centrality option %d, use 10, 20 or 100\n",fCentralityOpt));
      return -1;
    }
  }
  else
  {
    if ( !GetMultSelCen() ) return -1;
    
    AliDebug(1,Form("Mult. Percentile: V0M %2.2f, CL0 %2.2f, CL1 %2.2f; selected class %s", 
                    GetMultSelCen()->GetMultiplicityPercentile("V0M",1), 
                    GetMultSelCen()->GetMultiplicityPercentile("CL0",1), 
                    GetMultSelCen()->GetMultiplicityPercentile("CL1",1), 
                    fCentralityClass.Data()));
    
    return GetMultSelCen()->GetMultiplicityPercentile(fCentralityClass, kTRUE); // returns centrality only for events used in calibration
    
    // equivalent to
    //GetMultSelCen()->GetMultiplicityPercentile("V0M", kFALSE); // returns centrality for any event
    //Int_t    qual = GetMultSelCen()->GetEvSelCode(); if (qual ! = 0) cent = qual;
  }
}

//_____________________________________________________
/// \return Current event plane angle.
/// Different methods options can be requested.
//_____________________________________________________
Double_t AliCaloTrackReader::GetEventPlaneAngle() const
{  
  if( !GetEventPlane() ) return -1000;
  
  Float_t ep =  GetEventPlane()->GetEventplane(GetEventPlaneMethod(), GetInputEvent());
  
  if(GetEventPlaneMethod()=="Q" && (ep < 0 || ep > TMath::Pi()))
  {
    AliDebug(1,Form("Bad EP for <Q> method : %f\n",ep));
    return -1000;
  }
  else if(GetEventPlaneMethod().Contains("V0")  )
  {
    if((ep > TMath::Pi()/2 || ep < -TMath::Pi()/2))
    {
      AliDebug(1,Form("Bad EP for <%s> method : %f\n",GetEventPlaneMethod().Data(), ep));
      return -1000;
    }
    
    ep+=TMath::Pi()/2; // put same range as for <Q> method
  }
  
  AliDebug(3,Form("Event plane angle %f",ep));
    
//    if(fDebug > 0 )
//    {
//      if     (ep > TMath::Pi()) printf("AliCaloTrackReader::GetEventPlaneAngle() - Too large angle = %f\n",ep);
//      else if(ep < 0          ) printf("AliCaloTrackReader::GetEventPlaneAngle() - Negative angle = %f\n" ,ep);
//    }
    
  return ep;
}

//__________________________________________________________
/// \return Vertex position to be used for single event analysis.
//__________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3]) const
{
  vertex[0] = fVertex[0][0];
  vertex[1] = fVertex[0][1];
  vertex[2] = fVertex[0][2];
}

//__________________________________________________________________________
/// \return Vertex position for mixed event, recover the vertex in a particular event.
//__________________________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3], Int_t evtIndex) const
{  
  vertex[0] = fVertex[evtIndex][0];  
  vertex[1] = fVertex[evtIndex][1];  
  vertex[2] = fVertex[evtIndex][2];
}

//________________________________________
/// Fill data member fVertex. 
/// In case of Mixed event, multiple vertices.
//________________________________________
void AliCaloTrackReader::FillVertexArray()
{  
  // Delete previous vertex
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++)
    {
      delete [] fVertex[i] ;
    }
    delete [] fVertex ;
  }
  
  fVertex = new Double_t*[fNMixedEvent] ;
  for (Int_t i = 0; i < fNMixedEvent; i++)
  {
    fVertex[i] = new Double_t[3] ;
    fVertex[i][0] = 0.0 ;
    fVertex[i][1] = 0.0 ;
    fVertex[i][2] = 0.0 ;
  }
  
  if (!fMixedEvent)
  { // Single event analysis
    if(fDataType!=kMC)
    {
      
      if(fInputEvent->GetPrimaryVertex())
      {
        fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex[0]);
      }
      else
      {
        AliWarning("NULL primary vertex");
        fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
      }//Primary vertex pointer do not exist
      
    } else
    {// MC read event
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
    
    AliDebug(1,Form("Single Event Vertex : %f,%f,%f\n",fVertex[0][0],fVertex[0][1],fVertex[0][2]));
    
  } else
  { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++)
    {
      if (fMixedEvent->GetVertexOfEvent(iev))
        fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
      else
         AliWarning("No vertex found");
      
      AliDebug(1,Form("Multi Event %d Vertex : %f,%f,%f",iev,fVertex[iev][0],fVertex[iev][1],fVertex[iev][2]));
    }
  }
}

//_____________________________________
/// Fill the array with Central Tracking System (CTS) 
/// filtered tracks. To select the tracks, kinematic cuts, DCA, 
/// re-fit status and timing cuts are applied. 
/// Other more ESD/AOD dependent cuts are applied in *SelectTrack()* method,
/// see AliCaloTrackAODReader and AliCaloTrackESDReader.
//_____________________________________
void AliCaloTrackReader::FillInputCTS()
{  
  AliDebug(1,"Begin");
  
  Double_t pTrack[3] = {0,0,0};
  
  Int_t nTracks = fInputEvent->GetNumberOfTracks() ;
  Int_t nstatus = 0;
  Double_t bz   = GetInputEvent()->GetMagneticField();
  
  for(Int_t i = 0; i < 19; i++)
  {
    fTrackBCEvent   [i] = 0;
    fTrackBCEventCut[i] = 0;
  }
  
  for(Int_t iptCut = 0; iptCut < fTrackMultNPtCut; iptCut++ )
  {
    fTrackMult [iptCut] = 0;
    fTrackSumPt[iptCut] = 0;
  }
  
  Bool_t   bc0  = kFALSE;
  if(fRecalculateVertexBC) fVertexBC = AliVTrack::kTOFBCNA;
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++)
  {////////////// track loop
    AliVTrack * track = (AliVTrack*)fInputEvent->GetTrack(itrack) ; // retrieve track from esd
    
    if ( !AcceptParticleMCLabel( TMath::Abs(track->GetLabel()) ) ) continue ;
    
    fhCTSTrackCutsPt[0]->Fill(track->Pt());
    
    //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
    ULong_t status = track->GetStatus();
    
    if (fTrackStatus && !((status & fTrackStatus) == fTrackStatus))
      continue ;

    fhCTSTrackCutsPt[1]->Fill(track->Pt());
    
    nstatus++;

    //-------------------------
    // Select the tracks depending on cuts of AOD or ESD
    if(!SelectTrack(track, pTrack)) continue ;
    
    fhCTSTrackCutsPt[2]->Fill(track->Pt());
    
    //-------------------------
    // TOF cuts
    Bool_t okTOF  = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
    Double_t tof  = -1000;
    Int_t trackBC = -1000 ;
    
    if(fAccessTrackTOF)
    {
      if(okTOF)
      {
        trackBC = track->GetTOFBunchCrossing(bz);
        SetTrackEventBC(trackBC+9);
        
        tof = track->GetTOFsignal()*1e-3;
        
        //After selecting tracks with small DCA, pointing to vertex, set vertex BC depeding on tracks BC
        if(fRecalculateVertexBC)
        {
          if     (trackBC != 0 && trackBC != AliVTrack::kTOFBCNA) fVertexBC = trackBC;
          else if(trackBC == 0)                                   bc0       = kTRUE;
        }
        
        //In any case, the time should to be larger than the fixed window ...
        if( fUseTrackTimeCut && (trackBC !=0 || tof < fTrackTimeCutMin  || tof > fTrackTimeCutMax) )
        {
          //printf("Remove track time %f and bc = %d\n",tof,trackBC);
          continue ;
        }
        //else printf("Accept track time %f and bc = %d\n",tof,trackBC);
      }
    }
    
    fhCTSTrackCutsPt[3]->Fill(track->Pt());

    //---------------------
    // DCA cuts
    //
    fMomentum.SetPxPyPzE(pTrack[0],pTrack[1],pTrack[2],0);
        
    if(fUseTrackDCACut)
    {      
      Float_t dcaTPC =-999;
      //In case of AODs, TPC tracks cannot be propagated back to primary vertex,
      if( fDataType == kAOD ) dcaTPC = ((AliAODTrack*) track)->DCA();

      //normal way to get the dca, cut on dca_xy
      if(dcaTPC==-999)
      {
        Double_t dca[2]   = {1e6,1e6};
        Double_t covar[3] = {1e6,1e6,1e6};
        Bool_t okDCA = track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),bz,100.,dca,covar);
        if( okDCA) okDCA = AcceptDCA(fMomentum.Pt(),dca[0]);
        if(!okDCA)
        {
          //printf("AliCaloTrackReader::FillInputCTS() - Reject track pt %2.2f, dca_xy %2.4f\n",fMomentum.Pt(),dca[0]);
          continue ;
        }
      }
    }// DCA cuts
    
    fhCTSTrackCutsPt[4]->Fill(track->Pt());

    //-------------------------
    // Kinematic/acceptance cuts
    //
    // Count the tracks in eta < 0.9 and different pT cuts
    Float_t ptTrack = fMomentum.Pt();
    if(TMath::Abs(track->Eta())< fTrackMultEtaCut) 
    {
      for(Int_t iptCut = 0; iptCut < fTrackMultNPtCut; iptCut++ )
      {
        if(ptTrack > fTrackMultPtCut[iptCut]) 
        {
          fTrackMult [iptCut]++;
          fTrackSumPt[iptCut]+=ptTrack;
        }
      }
    }
    
    if(fCTSPtMin > ptTrack || fCTSPtMax < ptTrack) continue ;
    
    // Check effect of cuts on track BC
    if(fAccessTrackTOF && okTOF) SetTrackEventBCcut(trackBC+9);
    
    if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kCTS)) continue;
    
    fhCTSTrackCutsPt[5]->Fill(track->Pt());
        
    // ------------------------------
    // Add selected tracks to array
    AliDebug(2,Form("Selected tracks pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));
        
    fCTSTracks->Add(track);
    
    // TODO, check if remove
    if (fMixedEvent)  track->SetID(itrack);
    
  }// track loop
	
  if( fRecalculateVertexBC && (fVertexBC == 0 || fVertexBC == AliVTrack::kTOFBCNA))
  {
    if( bc0 ) fVertexBC = 0 ;
    else      fVertexBC = AliVTrack::kTOFBCNA ;
  }
  
  AliDebug(1,Form("AOD entries %d, input tracks %d, pass status %d, multipliticy %d", fCTSTracks->GetEntriesFast(), nTracks, nstatus, fTrackMult[0]));//fCTSTracksNormalInputEntries);
}

//_______________________________________________________________________________
/// Correct, if requested, and select here the EMCal cluster.
/// If selected add it to the EMCal clusters array.
/// The actions taken are:
///
///   * If requested, recalibrate and recalculate most of the cluster parameters (careful not to be done if tender applied or other tasks executed after)
///   * Select clusters without bad channels, exotic channels or close to borders
///   * If requested, correct cluster non linearity (careful not to be done if tender applied or other tasks executed after)
///   * Select clusters within an energy window and passing fiducial cuts
///   * Select clusters within a time window
///   * Select clusters with a minimum number of cells and not too close to a bad channel
///   * Smear the shower shape, to be done only for MC
///   * Besides, some counters on the number of clusters with time in different BC are stored
///
/// \param clus: AliVCaloCluster pointer
/// \param iclus: cluster index, only needed in case of mixing frame (not used recently)
///
/// Method called by *FillInputEMCAL()*
//_______________________________________________________________________________
void AliCaloTrackReader::FillInputEMCALAlgorithm(AliVCluster * clus, Int_t iclus)
{
  // Accept clusters with the proper label
  if ( clus->GetLabel() >= 0 )  // -1 corresponds to noisy MC
  { 
    if ( !AcceptParticleMCLabel(clus->GetLabel()) ) return ;
  }
  
  // TODO, not sure if needed anymore
  Int_t vindex = 0 ;
  if (fMixedEvent)
    vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
    
  clus->GetMomentum(fMomentum, fVertex[vindex]);

  // No correction/cut applied yet
  fhEMCALClusterCutsE[0]->Fill(clus->E());

  //if( (fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10 )
  AliDebug(2,Form("Input cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                  fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));

  //---------------------------
  // Embedding case
  // TO BE REVISED
  if(fSelectEmbeddedClusters)
  {
    if(clus->GetNLabels()==0 || clus->GetLabel() < 0) return;
    //else printf("Embedded cluster,  %d, n label %d label %d  \n",iclus,clus->GetNLabels(),clus->GetLabel());
  }

  //--------------------------------------
  // Apply some corrections in the cluster
  //
  if(fRecalculateClusters)
  {
    //Recalibrate the cluster energy
    if(GetCaloUtils()->IsRecalibrationOn())
    {
      Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, GetEMCALCells());
      
      clus->SetE(energy);
      //printf("Recalibrated Energy %f\n",clus->E());
      
      GetCaloUtils()->RecalculateClusterShowerShapeParameters(GetEMCALCells(),clus);
      GetCaloUtils()->RecalculateClusterPID(clus);
      
    } // recalculate E
    
    //Recalculate distance to bad channels, if new list of bad channels provided
    GetCaloUtils()->RecalculateClusterDistanceToBadChannel(GetEMCALCells(),clus);
    
    //Recalculate cluster position
    if(GetCaloUtils()->IsRecalculationOfClusterPositionOn())
    {
      GetCaloUtils()->RecalculateClusterPosition(GetEMCALCells(),clus);
      //clus->GetPosition(pos);
      //printf("After  Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
    }
    
    // Recalculate TOF
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn())
    {
      Double_t tof      = clus->GetTOF();
      Float_t  frac     =-1;
      Int_t    absIdMax = GetCaloUtils()->GetMaxEnergyCell(fEMCALCells, clus,frac);
      
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
      
      //additional L1 phase shift
      if(GetCaloUtils()->GetEMCALRecoUtils()->IsL1PhaseInTimeRecalibrationOn())
      {
        GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absIdMax), 
                                                                        fInputEvent->GetBunchCrossNumber(), tof);
      }

      clus->SetTOF(tof);
      
    }// Time recalibration
  }
  
  // Check effect of corrections
  fhEMCALClusterCutsE[1]->Fill(clus->E());

  //-----------------------------------------------------------------
  // Reject clusters with bad channels, close to borders and exotic
  //
  Bool_t goodCluster = GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,
                                                                          GetCaloUtils()->GetEMCALGeometry(),
                                                                          GetEMCALCells(),fInputEvent->GetBunchCrossNumber());
  
  if(!goodCluster)
  {
    //if( (fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10 )
    AliDebug(2,Form("Bad cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));

    return;
  }
    
  // Check effect of bad cluster removal 
  fhEMCALClusterCutsE[2]->Fill(clus->E());

  //Float_t pos[3];
  //clus->GetPosition(pos);
  //printf("Before Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
  
  //--------------------------------------
  // Correct non linearity or smear energy
  //
  if(fCorrectELinearity && GetCaloUtils()->IsCorrectionOfClusterEnergyOn())
  {
    GetCaloUtils()->CorrectClusterEnergy(clus) ;
    
    //if( (fDebug > 5 && fMomentum.E() > 0.1) || fDebug > 10 )
    AliDebug(5,Form("Correct Non Lin: Old E %3.2f, New E %3.2f",
                    fMomentum.E(),clus->E()));

    // In case of MC analysis, to match resolution/calibration in real data
    // Not needed anymore, just leave for MC studies on systematics
    if( GetCaloUtils()->GetEMCALRecoUtils()->IsClusterEnergySmeared() )
    {
      Float_t rdmEnergy = GetCaloUtils()->GetEMCALRecoUtils()->SmearClusterEnergy(clus);
      
      //if( (fDebug > 5 && fMomentum.E() > 0.1) || fDebug > 10 )
      AliDebug(5,Form("Smear energy: Old E %3.2f, New E %3.2f",clus->E(),rdmEnergy));
    
      clus->SetE(rdmEnergy);
    }
  }
  
  clus->GetMomentum(fMomentum, fVertex[vindex]);

  // Check effect linearity correction, energy smearing
  fhEMCALClusterCutsE[3]->Fill(clus->E());

  // Check the event BC depending on EMCal clustr before final cuts
  Double_t tof = clus->GetTOF()*1e9;
  
  Int_t bc = TMath::Nint(tof/50) + 9;
  //printf("tof %2.2f, bc+5=%d\n",tof,bc);
  
  SetEMCalEventBC(bc);
  
  //--------------------------------------
  // Apply some kinematical/acceptance cuts
  //
  if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) return ;

  // Select cluster fiducial region
  //
  Bool_t bEMCAL = kFALSE;
  Bool_t bDCAL  = kFALSE;
  if(fCheckFidCut)
  {
    if(fFillEMCAL && fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kEMCAL)) bEMCAL = kTRUE ;
    if(fFillDCAL  && fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kDCAL )) bDCAL  = kTRUE ;
  }
  else
  {
    bEMCAL = kTRUE;
  }

  //---------------------------------------------------------------------
  // Mask all cells in collumns facing ALICE thick material if requested
  //
  if(GetCaloUtils()->GetNMaskCellColumns())
  {
    Int_t absId   = -1;
    Int_t iSupMod = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Bool_t shared = kFALSE;
    GetCaloUtils()->GetEMCALRecoUtils()->GetMaxEnergyCell(GetCaloUtils()->GetEMCALGeometry(), GetEMCALCells(),clus,absId,iSupMod,ieta,iphi,shared);
    
    AliDebug(2,Form("Masked collumn: cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));
    
    
    if(GetCaloUtils()->MaskFrameCluster(iSupMod, ieta)) return;
  }
  
  // Check effect of energy and fiducial cuts
  fhEMCALClusterCutsE[4]->Fill(clus->E());
  
  
  //------------------------------------------
  // Apply time cut, count EMCal BC before cut
  //
  SetEMCalEventBCcut(bc);
  
  if(!IsInTimeWindow(tof,clus->E()))
  {
    fNPileUpClusters++ ;
    if(fUseEMCALTimeCut) 
    {
      AliDebug(2,Form("Out of time window E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f, time %e",
                      fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta(),tof));

      return ;
    }
  }
  else
    fNNonPileUpClusters++;
    
  // Check effect of time cut
  fhEMCALClusterCutsE[5]->Fill(clus->E());
  
  
  //----------------------------------------------------
  // Apply N cells cut
  //
  if(clus->GetNCells() <= fEMCALNCellsCut && fDataType != AliCaloTrackReader::kMC) return ;

  // Check effect of n cells cut
  fhEMCALClusterCutsE[6]->Fill(clus->E());

  //----------------------------------------------------
  // Apply distance to bad channel cut
  //
  Double_t distBad = clus->GetDistanceToBadChannel() ; //Distance to bad channel
  
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  
  if(distBad < fEMCALBadChMinDist) return  ;
  
  // Check effect distance to bad channel cut
  fhEMCALClusterCutsE[7]->Fill(clus->E());

  //----------------------------------------------------
  // Smear the SS to try to match data and simulations,
  // do it only for simulations.
  //
  Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(clus, GetEMCALCells()); 
  // Int_t nMaxima = clus->GetNExMax(); // For Run2
  if( fSmearShowerShape  && clus->GetNCells() > 2 && 
      nMaxima >= fSmearNLMMin && nMaxima <= fSmearNLMMax )
  {
    AliDebug(2,Form("Smear shower shape - Original: %2.4f\n", clus->GetM02()));
    if(fSmearingFunction == kSmearingLandau)
    {
      clus->SetM02( clus->GetM02() + fRandom.Landau(0, fSmearShowerShapeWidth) );
    }
    else if(fSmearingFunction == kSmearingLandauShift)
    {
      if(iclus%3 == 0 && clus->GetM02() > 0.1) clus->SetM02( clus->GetM02() + fRandom.Landau(0.05, fSmearShowerShapeWidth) );     //fSmearShowerShapeWidth = 0.035
    }
    else if (fSmearingFunction == kNoSmearing)
    {
      clus->SetM02( clus->GetM02() );
    }
    //clus->SetM02( fRandom.Landau(clus->GetM02(), fSmearShowerShapeWidth) );
    AliDebug(2,Form("Width %2.4f         Smeared : %2.4f\n", fSmearShowerShapeWidth,clus->GetM02()));
  }

  //--------------------------------------------------------
  // Fill the corresponding array with the selected clusters
  // Usually just filling EMCal array with upper or lower clusters is enough, 
  // but maybe we want to do EMCal-DCal correlations.
  
  //if((fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10)
  AliDebug(2,Form("Selected clusters (EMCAL%d, DCAL%d), E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                  bEMCAL,bDCAL,fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));

  
  if     (bEMCAL) fEMCALClusters->Add(clus);
  else if(bDCAL ) fDCALClusters ->Add(clus);
  
  // TODO, not sure if needed anymore
  if (fMixedEvent)
    clus->SetID(iclus) ;
}

//_______________________________________
/// Fill the array with EMCAL clusters. 
/// Source of clusters can be different, 
/// external branch, output of some other 
/// analysis task or the standard.
/// External branch is requested when providing 
/// its name in *fEMCALClustersListName*.
//_______________________________________
void AliCaloTrackReader::FillInputEMCAL()
{  
  AliDebug(1,"Begin");
  
  // First recalibrate cells, time or energy
  //  if(GetCaloUtils()->IsRecalibrationOn())
  //    GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCells(GetCaloUtils()->GetEMCALGeometry(),
  //                                                          GetEMCALCells(),
  //                                                          fInputEvent->GetBunchCrossNumber());
  
  fNPileUpClusters    = 0; // Init counter
  fNNonPileUpClusters = 0; // Init counter
  for(Int_t i = 0; i < 19; i++)
  {
    fEMCalBCEvent   [i] = 0;
    fEMCalBCEventCut[i] = 0;
  }
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  if(fEMCALClustersListName=="")
  {
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = 0;
      if ( (clus = fInputEvent->GetCaloCluster(iclus)) )
      {
        if (clus->IsEMCAL())
        {
          FillInputEMCALAlgorithm(clus, iclus);
        }//EMCAL cluster
      }// cluster exists
    }// cluster loop
    
    //Recalculate track matching
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent,0x0,fMC);
    
  }//Get the clusters from the input event
  else
  {
    TClonesArray * clusterList = 0x0;
    
    if      (fInputEvent->FindListObject(fEMCALClustersListName))
    {
      clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(fEMCALClustersListName));
    }
    else if(fOutputEvent)
    {
      clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
    }
    
    if(!clusterList)
    {
      AliWarning(Form("Wrong name of list with clusters?  <%s>",fEMCALClustersListName.Data()));
      return;
    }
    
    Int_t nclusters = clusterList->GetEntriesFast();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
      //printf("E %f\n",clus->E());
      if (clus) FillInputEMCALAlgorithm(clus, iclus);
      else      AliWarning("Null cluster in list!");
    }// cluster loop
    
    // Recalculate the pile-up time, in case long time clusters removed during clusterization
    //printf("Input event INIT : Pile-up clusters %d, NO pile-up %d\n",fNPileUpClusters,fNNonPileUpClusters);
    
    fNPileUpClusters    = 0; // Init counter
    fNNonPileUpClusters = 0; // Init counter
    for(Int_t i = 0; i < 19; i++)
    {
      fEMCalBCEvent   [i] = 0;
      fEMCalBCEventCut[i] = 0;
    }
    
    for (Int_t iclus =  0; iclus < fInputEvent->GetNumberOfCaloClusters(); iclus++)
    {
      AliVCluster * clus = 0;
      
      if ( (clus = fInputEvent->GetCaloCluster(iclus)) )
      {
        if (clus->IsEMCAL())
        {
          
          Float_t  frac     =-1;
          Int_t    absIdMax = GetCaloUtils()->GetMaxEnergyCell(fEMCALCells, clus,frac);
          Double_t tof = clus->GetTOF();
          GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
          //additional L1 phase shift
          if(GetCaloUtils()->GetEMCALRecoUtils()->IsL1PhaseInTimeRecalibrationOn())
          {
            GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absIdMax), fInputEvent->GetBunchCrossNumber(), tof);
          }

          tof*=1e9;
          
          //printf("Input event cluster : AbsIdMax %d, E %2.2f, time %2.2f \n", absIdMax,clus->E(),tof);
          
          //Reject clusters with bad channels, close to borders and exotic;
          if(!GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,GetCaloUtils()->GetEMCALGeometry(),GetEMCALCells(),fInputEvent->GetBunchCrossNumber()))  continue;
          
          Int_t bc = TMath::Nint(tof/50) + 9;
          SetEMCalEventBC(bc);
          
          if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) continue ;
          
          clus->GetMomentum(fMomentum, fVertex[0]);
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kEMCAL)) return ;
          
          SetEMCalEventBCcut(bc);
          
          if(!IsInTimeWindow(tof,clus->E()))
            fNPileUpClusters++ ;
          else
            fNNonPileUpClusters++;
          
        }
      }
    }
    
    //printf("Input event : Pile-up clusters %d, NO pile-up %d\n",fNPileUpClusters,fNNonPileUpClusters);
    
    // Recalculate track matching, not necessary if already done in the reclusterization task.
    // in case it was not done ...
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent,clusterList,fMC);
    
  }
  
  AliDebug(1,Form("AOD entries %d, n pile-up clusters %d, n non pile-up %d", fEMCALClusters->GetEntriesFast(),fNPileUpClusters,fNNonPileUpClusters));
}

//_______________________________________
/// Fill the array with PHOS filtered clusters. 
//_______________________________________
void AliCaloTrackReader::FillInputPHOS()
{  
  AliDebug(1,"Begin");
  
  // Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
  TString genName;
  for (Int_t iclus = 0; iclus < nclusters; iclus++)
  {
    AliVCluster * clus = fInputEvent->GetCaloCluster(iclus) ;
    if ( !clus ) continue ;
        
    if ( !clus->IsPHOS() ) continue ;
        
    if(clus->GetLabel() >=0 ) // -1 corresponds to noisy MC
    {
      if ( !AcceptParticleMCLabel(clus->GetLabel()) ) continue ;
    }
    
    fhPHOSClusterCutsE[0]->Fill(clus->E());
    
    // Skip CPV input
    if( clus->GetType() == AliVCluster::kPHOSCharged ) continue ;
    
    fhPHOSClusterCutsE[1]->Fill(clus->E());
    
    //---------------------------------------------
    // Old feature, try to rely on PHOS tender
    //
    if(fRecalculateClusters)
    {
      // Recalibrate the cluster energy
      if(GetCaloUtils()->IsRecalibrationOn())
      {
        Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliAODCaloCells*)GetPHOSCells());
        clus->SetE(energy);
      }
    }
    
    //----------------------------------------------------------------------------------
    // Check if the cluster contains any bad channel and if close to calorimeter borders    
    //
    // Old feature, try to rely on PHOS tender
    if( GetCaloUtils()->ClusterContainsBadChannel(kPHOS,clus->GetCellsAbsId(), clus->GetNCells()))
      continue;
    
    if(!GetCaloUtils()->CheckCellFiducialRegion(clus, fInputEvent->GetPHOSCells()))
      continue;
 
    // TODO, add exotic cut???
    
    fhPHOSClusterCutsE[2]->Fill(clus->E());

    // TODO Dead code? remove?
    Int_t vindex = 0 ;
    if (fMixedEvent)
      vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
    
    clus->GetMomentum(fMomentum, fVertex[vindex]);
    
    //----------------------------------------------------------------------------------
    // Remove clusters close to borders
    //
    if (fCheckFidCut && !fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kPHOS) ) 
      continue ;
    
    fhPHOSClusterCutsE[3]->Fill(clus->E());
    
    //----------------------------------------------------------------------------------
    // Remove clusters with too low energy
    //
    if (fPHOSPtMin > fMomentum.E() || fPHOSPtMax < fMomentum.E() ) 
      continue ;
    
    fhPHOSClusterCutsE[4]->Fill(clus->E());

    //----------------------------------------------------
    // Apply N cells cut
    //
    if(clus->GetNCells() <= fPHOSNCellsCut && fDataType != AliCaloTrackReader::kMC) return ;
    
    // Check effect of n cells cut
    fhPHOSClusterCutsE[5]->Fill(clus->E());
    
    //----------------------------------------------------
    // Apply distance to bad channel cut
    //
    Double_t distBad = clus->GetDistanceToBadChannel() ; //Distance to bad channel
    
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    
    if(distBad < fPHOSBadChMinDist) return  ;
    
    // Check effect distance to bad channel cut
    fhPHOSClusterCutsE[6]->Fill(clus->E());

    // TODO, add time cut

    //----------------------------------------------------------------------------------
    // Add selected clusters to array
    //
    //if(fDebug > 2 && fMomentum.E() > 0.1)
    AliDebug(2,Form("Selected clusters E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));
    
    fPHOSClusters->Add(clus);

    // TODO Dead code? remove?
    if (fMixedEvent)
      clus->SetID(iclus) ;    
    
  } // esd/aod cluster loop
  
  AliDebug(1,Form("AOD entries %d",fPHOSClusters->GetEntriesFast())) ;  
}

//____________________________________________
/// Connects the array with EMCAL cells and the pointer.
//____________________________________________
void AliCaloTrackReader::FillInputEMCALCells()
{  
  fEMCALCells = fInputEvent->GetEMCALCells();
}

//___________________________________________
/// Connects the array with PHOS cells and the pointer.
//___________________________________________
void AliCaloTrackReader::FillInputPHOSCells()
{  
  fPHOSCells = fInputEvent->GetPHOSCells();
}

//_______________________________________
/// Fill VZERO information in data member, 
/// add all the channels information.
//_______________________________________
void AliCaloTrackReader::FillInputVZERO()
{
  AliVVZERO* v0 = fInputEvent->GetVZEROData();
  //printf("Init V0: ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  
  if (v0)
  {
    AliESDVZERO* esdV0 = dynamic_cast<AliESDVZERO*> (v0);
    for (Int_t i = 0; i < 32; i++)
    {
      if(esdV0)
      {//Only available in ESDs
        fV0ADC[0] += (Int_t)esdV0->GetAdcV0C(i);
        fV0ADC[1] += (Int_t)esdV0->GetAdcV0A(i);
      }
      
      fV0Mul[0] += (Int_t)v0->GetMultiplicityV0C(i);
      fV0Mul[1] += (Int_t)v0->GetMultiplicityV0A(i);
    }
    
    AliDebug(1,Form("ADC (%d,%d), Multiplicity (%d,%d)",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]));
  }
  else
  {
    AliDebug(1,"Cannot retrieve V0 ESD! Run w/ null V0 charges");
  }
}

//_________________________________________________
/// Fill array with non standard jets
///
/// Author: Adam T. Matyja
//_________________________________________________
void AliCaloTrackReader::FillInputNonStandardJets()
{  
  AliDebug(2,"Begin");
  
  //
  //check if branch name is given
  if(!fInputNonStandardJetBranchName.Length())
  {
    fInputEvent->Print();
    AliFatal("No non-standard jet branch name specified. Specify among existing ones.");
  }
  
  fNonStandardJets = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fInputNonStandardJetBranchName.Data()));
  
  if(!fNonStandardJets)
  {
    //check if jet branch exist; exit if not
    fInputEvent->Print();

    AliFatal(Form("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fInputNonStandardJetBranchName.Data()));
  }
  else
  {
    AliDebug(1,Form("AOD input jets %d", fNonStandardJets->GetEntriesFast()));
  }
}

//_________________________________________________
/// Fill array with Background jets
///
/// Author: Adam T. Matyja
//_________________________________________________
void AliCaloTrackReader::FillInputBackgroundJets()
{
  AliDebug(1,"Begin");
  //
  //check if branch name is given
  if(!fInputBackgroundJetBranchName.Length())
  {
    fInputEvent->Print();
    
    AliFatal("No background jet branch name specified. Specify among existing ones.");
  }
  
  fBackgroundJets = (AliAODJetEventBackground*)(fInputEvent->FindListObject(fInputBackgroundJetBranchName.Data()));
  
  if(!fBackgroundJets)
  {
    //check if jet branch exist; exit if not
    fInputEvent->Print();

    AliFatal(Form("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fInputBackgroundJetBranchName.Data()));
  }
  else
  {
    AliDebug(1,"FillInputBackgroundJets");
    fBackgroundJets->Print("");
  }
}

//____________________________________________________________________
/// Recover the patches that triggered, either L0 or L1.
///
/// \param tmin: minimum L0 time bin cut
/// \param tmax: maximum L0 time bin cut
/// \return TArrayI, array with patches index
//____________________________________________________________________
TArrayI AliCaloTrackReader::GetTriggerPatches(Int_t tmin, Int_t tmax )
{
  // init some variables
  Int_t  trigtimes[30], globCol, globRow,ntimes, i;
  Int_t  absId  = -1; //[100];
  Int_t  nPatch = 0;
  
  TArrayI patches(0);
  
  // get object pointer
  AliVCaloTrigger *caloTrigger = GetInputEvent()->GetCaloTrigger( "EMCAL" );
  
  if(!caloTrigger) 
  {
    AliError("Trigger patches input (AliVCaloTrigger) not available in data!");
    return patches;
  }
  
  //printf("CaloTrigger Entries %d\n",caloTrigger->GetEntries() );
  
  // class is not empty
  if( caloTrigger->GetEntries() > 0 )
  {
    // must reset before usage, or the class will fail
    caloTrigger->Reset();

    // go throuth the trigger channels
    while( caloTrigger->Next() )
    {
      // get position in global 2x2 tower coordinates
      caloTrigger->GetPosition( globCol, globRow );

      //L0
      if(IsEventEMCALL0())
      {
        // get dimension of time arrays
        caloTrigger->GetNL0Times( ntimes );
        
        // no L0s in this channel
        // presence of the channel in the iterator still does not guarantee that L0 was produced!!
        if( ntimes < 1 )
          continue;
        
        // get timing array
        caloTrigger->GetL0Times( trigtimes );
        //printf("Get L0 patch : n times %d - trigger time window %d - %d\n",ntimes, tmin,tmax);
        
        // go through the array
        for( i = 0; i < ntimes; i++ )
        {
          // check if in cut - 8,9 shall be accepted in 2011
          if( trigtimes[i] >= tmin && trigtimes[i] <= tmax )
          {
            //printf("Accepted trigger time %d \n",trigtimes[i]);
            //if(nTrig > 99) continue;
            GetCaloUtils()->GetEMCALGeometry()->GetAbsFastORIndexFromPositionInEMCAL(globCol,globRow, absId);
            //printf("pass the time cut globCol %d, globRow %d absId %d\n",globCol,globRow, absId);
            patches.Set(nPatch+1);
            patches.AddAt(absId,nPatch++);
          }
        } // trigger time array
      }//L0
      else if(IsEventEMCALL1()) // L1
      {
        Int_t bit = 0;
        caloTrigger->GetTriggerBits(bit);
        
        Int_t sum = 0;
        caloTrigger->GetL1TimeSum(sum);
        //fBitEGA-=2;
        Bool_t isEGA1 = ((bit >>  fBitEGA   ) & 0x1) && IsEventEMCALL1Gamma1() ;
        Bool_t isEGA2 = ((bit >> (fBitEGA+1)) & 0x1) && IsEventEMCALL1Gamma2() ;
        Bool_t isEJE1 = ((bit >>  fBitEJE   ) & 0x1) && IsEventEMCALL1Jet1  () ;
        Bool_t isEJE2 = ((bit >> (fBitEJE+1)) & 0x1) && IsEventEMCALL1Jet2  () ;
        
        //if((bit>> fBitEGA   )&0x1) printf("Trig Bit %d - bit %d - EG1 %d - EG2 %d\n",fBitEGA  ,bit,IsEventEMCALL1Gamma1(),IsEventEMCALL1Gamma2());
        //if((bit>>(fBitEGA+1))&0x1) printf("Trig Bit %d - bit %d - EG1 %d - EG2 %d\n",fBitEGA+1,bit,IsEventEMCALL1Gamma1(),IsEventEMCALL1Gamma2());
        
        if(!isEGA1 && !isEJE1 && !isEGA2 && !isEJE2) continue;
        
        Int_t patchsize = -1;
        if      (isEGA1 || isEGA2) patchsize =  2;
        else if (isEJE1 || isEJE2) patchsize = 16;
        
        //printf("**** Get L1 Patch: Bit %x, sum %d, patchsize %d, EGA1 %d, EGA2 %d, EJE1 %d, EJE2 %d, EGA bit %d, EJE bit %d, Trigger Gamma %d, Trigger Jet %d\n",
        //       bit,sum,patchsize,isEGA1,isEGA2,isEJE1,isEJE2,fBitEGA,fBitEJE,IsEventEMCALL1Gamma(),IsEventEMCALL1Jet());

        
        // add 2x2 (EGA) or 16x16 (EJE) patches
        for(Int_t irow=0; irow < patchsize; irow++)
        {
          for(Int_t icol=0; icol < patchsize; icol++)
          {
            GetCaloUtils()->GetEMCALGeometry()->GetAbsFastORIndexFromPositionInEMCAL(globCol+icol,globRow+irow, absId);
            //printf("pass the time cut globCol %d, globRow %d absId %d\n",globCol,globRow, absId);
            patches.Set(nPatch+1);
            patches.AddAt(absId,nPatch++);
          }
        }
        
      } // L1
      
    } // trigger iterator
  } // go through triggers
  
  if(patches.GetSize()<=0) AliInfo(Form("No patch found! for triggers: %s and selected <%s>",
                                        GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data()));
  //else                     printf(">>>>> N patches %d, test %d,first %d, last %d\n",patches.GetSize(), nPatch, patches.At(0), patches.At(patches.GetSize()-1));
                 
  return patches;
}

//____________________________________________________________
/// Finds the cluster that triggered.
/// It compares the cells of the trigger patches and 
/// high energy clusters.
//____________________________________________________________
void  AliCaloTrackReader::MatchTriggerCluster(TArrayI patches)
{
  // Init info from previous event
  fTriggerClusterIndex = -1;
  fTriggerClusterId    = -1;
  fTriggerClusterBC    = -10000;
  fIsExoticEvent       = kFALSE;
  fIsBadCellEvent      = kFALSE;
  fIsBadMaxCellEvent   = kFALSE;
  
  fIsTriggerMatch           = kFALSE;
  fIsTriggerMatchOpenCut[0] = kFALSE;
  fIsTriggerMatchOpenCut[1] = kFALSE;
  fIsTriggerMatchOpenCut[2] = kFALSE;
  
  // Do only analysis for triggered events
  if(!IsEventEMCALL1() && !IsEventEMCALL0())
  {
    fTriggerClusterBC = 0;
    return;
  }
  
  //printf("***** Try to match trigger to cluster %d **** L0 %d, L1 %d\n",fTriggerPatchClusterMatch,IsEventEMCALL0(),IsEventEMCALL1());
  
  //Recover the list of clusters
  TClonesArray * clusterList = 0;
  if      (fInputEvent->FindListObject(fEMCALClustersListName))
  {
    clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(fEMCALClustersListName));
  }
  else if(fOutputEvent)
  {
    clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
  }
  
  // Get number of clusters and of trigger patches
  Int_t   nclusters   = fInputEvent->GetNumberOfCaloClusters();
  if(clusterList)
    nclusters    = clusterList->GetEntriesFast();
  
  Int_t   nPatch      = patches.GetSize();
  Float_t exoDiffTime = GetCaloUtils()->GetEMCALRecoUtils()->GetExoticCellDiffTimeCut();
  
  //Init some variables used in the cluster loop
  Float_t tofPatchMax = 100000;
  Float_t ePatchMax   =-1;
  
  Float_t tofMax      = 100000;
  Float_t eMax        =-1;
  
  Int_t   clusMax     =-1;
  Int_t   idclusMax   =-1;
  Bool_t  badClMax    = kFALSE;
  Bool_t  badCeMax    = kFALSE;
  Bool_t  exoMax      = kFALSE;
//  Int_t   absIdMaxTrig= -1;
  Int_t   absIdMaxMax = -1;
  
  Int_t   nOfHighECl  = 0 ;
  
  //
  // Check what is the trigger threshold
  // set minimu energym of candidate for trigger cluster
  //
  SetEMCALTriggerThresholds();
  
  Float_t triggerThreshold = fTriggerL1EventThreshold;
  if(IsEventEMCALL0()) triggerThreshold = fTriggerL0EventThreshold;
  Float_t minE = triggerThreshold / 2.;

  // This method is not really suitable for JET trigger
  // but in case, reduce the energy cut since we do not trigger on high energy particle
  if(IsEventEMCALL1Jet() || minE < 1) minE = 1;

  AliDebug(1,Form("IsL1Trigger %d, IsL1JetTrigger? %d, IsL0Trigger %d, L1 threshold %2.1f, L0 threshold %2.1f, Min cluster E %2.2f",IsEventEMCALL1Jet(), IsEventEMCALL1(), IsEventEMCALL0(), fTriggerL1EventThreshold,fTriggerL0EventThreshold,minE));  
  
  //
  // Loop on the clusters, check if there is any that falls into one of the patches
  //
  for (Int_t iclus =  0; iclus <  nclusters; iclus++)
  {
    AliVCluster * clus = 0;
    if(clusterList) clus = (AliVCluster*) clusterList->At(iclus);
    else            clus = fInputEvent->GetCaloCluster(iclus);
    
    if ( !clus )            continue ;
    
    if ( !clus->IsEMCAL() ) continue ;
    
    //Skip clusters with too low energy to be triggering
    if ( clus->E() < minE ) continue ;
    
    Float_t  frac       = -1;
    Int_t    absIdMax   = GetCaloUtils()->GetMaxEnergyCell(fInputEvent->GetEMCALCells(), clus,frac);
    
    Bool_t   badCluster = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),
                                                                                         clus->GetCellsAbsId(),clus->GetNCells());
    UShort_t cellMax[]  = {(UShort_t) absIdMax};
    Bool_t   badCell    = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),cellMax,1);
    
    // if cell is bad, it can happen that time calibration is not available,
    // when calculating if it is exotic, this can make it to be exotic by default
    // open it temporarily for this cluster
    if(badCell)
      GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(10000000);
    
    Bool_t   exotic     = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCluster(clus, fInputEvent->GetEMCALCells());
    
    // Set back the time cut on exotics
    if(badCell)
      GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(exoDiffTime);
    
    // Energy threshold for exotic Ecross typically at 4 GeV,
    // for lower energy, check that there are more than 1 cell in the cluster
    if(!exotic && clus->GetNCells() < 2) exotic = kTRUE;
    
    Float_t  energy     = clus->E();
    Int_t    idclus     = clus->GetID();
    
    Double_t tof        = clus->GetTOF();
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn() && fTriggerClusterTimeRecal){
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
      //additional L1 phase shift
      if(GetCaloUtils()->GetEMCALRecoUtils()->IsL1PhaseInTimeRecalibrationOn()){
	GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absIdMax), fInputEvent->GetBunchCrossNumber(), tof);
      }
    }
    tof *=1.e9;
    
    //printf("cluster %d, ID %d, E %2.2f, tof %2.2f, AbsId max %d, exotic %d, bad Cluster %d, bad Cell %d\n",
    //       iclus,idclus, energy,tof,absIdMax, exotic, badCluster,badCell);
    
    // Find the highest energy cluster, avobe trigger threshold
    // in the event in case no match to trigger is found
    if( energy > eMax )
    {
      tofMax       = tof;
      eMax         = energy;
      badClMax     = badCluster;
      badCeMax     = badCell;
      exoMax       = exotic;
      clusMax      = iclus;
      idclusMax    = idclus;
      absIdMaxMax  = absIdMax;
    }
    
    // count the good clusters in the event avobe the trigger threshold
    // to check the exotic events
    if(!badCluster && !exotic)
      nOfHighECl++;
    
    // Find match to trigger
    if(fTriggerPatchClusterMatch && nPatch>0)
    {
      for(Int_t iabsId =0; iabsId < nPatch; iabsId++)
      {
        Int_t absIDCell[4];
        GetCaloUtils()->GetEMCALGeometry()->GetCellIndexFromFastORIndex(patches.At(iabsId), absIDCell);
        //if(tof > 75 ) printf("E %2.2f TOF %2.2f Trigger patch %d, cells : %d, %d, %d, %d\n",
        //                     clus->E(),tof,patches.At(iabsId), absIDCell[0],absIDCell[1],absIDCell[2],absIDCell[3]);
        
        for(Int_t ipatch = 0; ipatch < 4; ipatch++)
        {
          if(absIdMax == absIDCell[ipatch])
          {
            //printf("*** Patches : absId %d, E %2.1f, tof %f \n",absIdMax,clus->E(), tof);
            if(energy > ePatchMax)
            {
              tofPatchMax          = tof;
              ePatchMax            = energy;
              fIsBadCellEvent      = badCluster;
              fIsBadMaxCellEvent   = badCell;
              fIsExoticEvent       = exotic;
              fTriggerClusterIndex = iclus;
              fTriggerClusterId    = idclus;
              fIsTriggerMatch      = kTRUE;
//              absIdMaxTrig         = absIdMax;
            }
          }
        }// cell patch loop
      }// trigger patch loop
    } // Do trigger patch matching
    
  }// Cluster loop
  
  // If there was no match, assign as trigger
  // the highest energy cluster in the event
  if(!fIsTriggerMatch)
  {
    tofPatchMax          = tofMax;
    ePatchMax            = eMax;
    fIsBadCellEvent      = badClMax;
    fIsBadMaxCellEvent   = badCeMax;
    fIsExoticEvent       = exoMax;
    fTriggerClusterIndex = clusMax;
    fTriggerClusterId    = idclusMax;
  }
  
  Double_t tofPatchMaxUS = TMath::Abs(tofPatchMax);
  
  if     (tofPatchMaxUS < 28 ) fTriggerClusterBC = 0 ;
  else if(tofPatchMaxUS < 75 ) fTriggerClusterBC = 1 ;
  else if(tofPatchMaxUS < 125) fTriggerClusterBC = 2 ;
  else if(tofPatchMaxUS < 175) fTriggerClusterBC = 3 ;
  else if(tofPatchMaxUS < 225) fTriggerClusterBC = 4 ;
  else if(tofPatchMaxUS < 275) fTriggerClusterBC = 5 ;
  else
  {
    //printf("AliCaloTrackReader::MatchTriggerCluster() - Large BC - tof %2.3f - Index %d\n",tofPatchMaxUS,fTriggerClusterIndex);
    if(fTriggerClusterIndex >= 0) fTriggerClusterBC = 6 ;
    else
    {
      fTriggerClusterIndex = -2;
      fTriggerClusterId    = -2;
    }
  }
  
  if(tofPatchMax < 0) fTriggerClusterBC*=-1;
  
  
  //  printf("AliCaloTrackReader::MatchTriggerCluster(TArrayI patches) - Trigger cluster: index %d, ID %d, E = %2.2f, tof = %2.2f (BC = %d), bad cluster? %d, bad cell? %d, exotic? %d, patch match? %d, n High E cluster %d, absId Max %d\n",
  //         fTriggerClusterIndex, fTriggerClusterId,ePatchMax, tofPatchMax,
  //         fTriggerClusterBC, fIsBadCellEvent,fIsBadMaxCellEvent,fIsExoticEvent, fIsTriggerMatch, nOfHighECl,absIdMaxMax);
  //
  //  if(!fIsTriggerMatch) printf("\t highest energy cluster:  index %d, ID %d, E = %2.2f, tof = %2.2f, bad cluster? %d, bad cell? %d, exotic? %d\n",
  //                              clusMax, idclusMax, eMax,tofMax, badClMax, badCeMax,exoMax);
  
  //Redo matching but open cuts
  if(!fIsTriggerMatch && fTriggerClusterId >= 0)
  {
    // Open time patch time
    TArrayI patchOpen = GetTriggerPatches(7,10);
    
    Int_t patchAbsIdOpenTime = -1;
    for(Int_t iabsId =0; iabsId < patchOpen.GetSize(); iabsId++)
    {
      Int_t absIDCell[4];
      patchAbsIdOpenTime = patchOpen.At(iabsId);
      GetCaloUtils()->GetEMCALGeometry()->GetCellIndexFromFastORIndex(patchAbsIdOpenTime, absIDCell);
      //if(tof > 75 ) printf("E %2.2f TOF %2.2f Trigger patch %d, cells : %d, %d, %d, %d\n",
      //                     clus->E(),tof,patches.At(iabsId), absIDCell[0],absIDCell[1],absIDCell[2],absIDCell[3]);
      
      for(Int_t ipatch = 0; ipatch < 4; ipatch++)
      {
        if(absIdMaxMax == absIDCell[ipatch])
        {
          fIsTriggerMatchOpenCut[0] = kTRUE;
          break;
        }
      }// cell patch loop
    }// trigger patch loop
    
    // Check neighbour patches
    Int_t patchAbsId = -1;
    Int_t globalCol  = -1;
    Int_t globalRow  = -1;
    GetCaloUtils()->GetEMCALGeometry()->GetFastORIndexFromCellIndex(absIdMaxMax, patchAbsId);
    GetCaloUtils()->GetEMCALGeometry()->GetPositionInEMCALFromAbsFastORIndex(patchAbsId,globalCol,globalRow);
    
    // Check patches with strict time cut
    Int_t patchAbsIdNeigh = -1;
    for(Int_t icol = globalCol-1; icol <= globalCol+1; icol++)
    {
      if(icol < 0 || icol > 47) continue;
      
      for(Int_t irow = globalRow; irow <= globalRow+1; irow++)
      {
        if(irow < 0 || irow > 63) continue;
        
        GetCaloUtils()->GetEMCALGeometry()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, patchAbsIdNeigh);
        
        if ( patchAbsIdNeigh < 0 ) continue;
        
        for(Int_t iabsId =0; iabsId < patches.GetSize(); iabsId++)
        {
          if(patchAbsIdNeigh == patches.At(iabsId))
          {
            fIsTriggerMatchOpenCut[1] = kTRUE;
            break;
          }
        }// trigger patch loop
        
      }// row
    }// col
    
    // Check patches with open time cut
    Int_t patchAbsIdNeighOpenTime = -1;
    for(Int_t icol = globalCol-1; icol <= globalCol+1; icol++)
    {
      if(icol < 0 || icol > 47) continue;
      
      for(Int_t irow = globalRow; irow <= globalRow+1; irow++)
      {
        if(irow < 0 || irow > 63) continue;
        
        GetCaloUtils()->GetEMCALGeometry()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, patchAbsIdNeighOpenTime);
        
        if ( patchAbsIdNeighOpenTime < 0 ) continue;
        
        for(Int_t iabsId =0; iabsId < patchOpen.GetSize(); iabsId++)
        {
          if(patchAbsIdNeighOpenTime == patchOpen.At(iabsId))
          {
            fIsTriggerMatchOpenCut[2] = kTRUE;
            break;
          }
        }// trigger patch loop
        
      }// row
    }// col
    
    //    printf("No match, new match: Open time %d-%d, open Neigh %d-%d, both open %d-%d\n",fIsTriggerMatchOpenCut[0],patchAbsIdOpenTime,
    //           fIsTriggerMatchOpenCut[1],patchAbsIdNeigh,
    //           fIsTriggerMatchOpenCut[2],patchAbsIdNeighOpenTime);
    
    patchOpen.Reset();
    
  }// No trigger match found
  //printf("Trigger BC %d, Id %d, Index %d\n",fTriggerClusterBC,fTriggerClusterId,fTriggerClusterIndex);
}

//_________________________________________________________
/// Recover the EMCal L1 trigger threshold from data.
/// Set the EMCal L0 threshold depending on the run number.
/// Set the threshold only if requested. 
//_________________________________________________________
void AliCaloTrackReader::SetEMCALTriggerThresholds()
{
  if(!fTriggerL1EventThresholdFix)
  {
    // get object pointer
    AliVCaloTrigger *caloTrigger = GetInputEvent()->GetCaloTrigger( "EMCAL" );

    if ( fBitEGA == 6 )
    {
      if     (IsEventEMCALL1Gamma1()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(1);
      else if(IsEventEMCALL1Gamma2()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(3);
      else if(IsEventEMCALL1Jet1  ()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(0);
      else if(IsEventEMCALL1Jet2  ()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(2);
      
      //      printf("L1 trigger Threshold Jet1 %f, Gamma1 %f, Jet2 %f, Gamma2 %f\n",
      //             0.07874*caloTrigger->GetL1Threshold(0),
      //             0.07874*caloTrigger->GetL1Threshold(1),
      //             0.07874*caloTrigger->GetL1Threshold(2),
      //             0.07874*caloTrigger->GetL1Threshold(3));
    }
    else
    {
      // Old AOD data format, in such case, order of thresholds different!!!
      if     (IsEventEMCALL1Gamma1()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(0);
      else if(IsEventEMCALL1Gamma2()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(2);
      else if(IsEventEMCALL1Jet1  ()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(1);
      else if(IsEventEMCALL1Jet2  ()) fTriggerL1EventThreshold =  0.07874*caloTrigger->GetL1Threshold(3);
      
      //      printf("L1 trigger Threshold Jet1 %f, Gamma1 %f, Jet2 %f, Gamma2 %f\n",
      //             0.07874*caloTrigger->GetL1Threshold(1),
      //             0.07874*caloTrigger->GetL1Threshold(0),
      //             0.07874*caloTrigger->GetL1Threshold(3),
      //             0.07874*caloTrigger->GetL1Threshold(2));
    }
  }
  
  // Set L0 threshold, if not set by user
  if( IsEventEMCALL0() && fTriggerL0EventThreshold < 0)
  { 
    // Revise for periods > LHC11d 
    Int_t runNumber = fInputEvent->GetRunNumber();
    if     (runNumber < 146861) fTriggerL0EventThreshold = 3. ;
    else if(runNumber < 154000) fTriggerL0EventThreshold = 4. ;
    else if(runNumber < 165000) fTriggerL0EventThreshold = 5.5;
    else                        fTriggerL0EventThreshold = 2  ;
  }  
}

//________________________________________________________
/// Print some relevant parameters set for the analysis
//________________________________________________________
void AliCaloTrackReader::Print(const Option_t * opt) const
{  
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n",    GetName(), GetTitle() ) ;
  printf("Task name      : %s\n",          fTaskName.Data()) ;
  printf("Data type      : %d\n",          fDataType) ;
  printf("CTS Min pT     : %2.1f GeV/c\n", fCTSPtMin) ;
  printf("EMCAL Min pT   : %2.1f GeV/c\n", fEMCALPtMin) ;
  printf("PHOS Min pT    : %2.1f GeV/c\n", fPHOSPtMin) ;
  printf("CTS Max pT     : %2.1f GeV/c\n", fCTSPtMax) ;
  printf("EMCAL Max pT   : %2.1f GeV/c\n", fEMCALPtMax) ;
  printf("PHOS Max pT    : %2.1f GeV/c\n", fPHOSPtMax) ;
  printf("EMCAL Bad Dist > %2.1f \n"     , fEMCALBadChMinDist) ;
  printf("PHOS  Bad Dist > %2.1f \n"     , fPHOSBadChMinDist) ;
  printf("EMCAL N cells  > %d \n"        , fEMCALNCellsCut) ;
  printf("PHOS  N cells  > %d \n"        , fPHOSNCellsCut) ;
  printf("EMCAL Time Cut: %3.1f < TOF  < %3.1f\n", fEMCALTimeCutMin, fEMCALTimeCutMax);
  printf("Use CTS         =     %d\n",     fFillCTS) ;
  printf("Use EMCAL       =     %d\n",     fFillEMCAL) ;
  printf("Use DCAL        =     %d\n",     fFillDCAL)  ;
  printf("Use PHOS        =     %d\n",     fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n",     fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n",     fFillPHOSCells) ;
  printf("Track status    =     %d\n", (Int_t) fTrackStatus) ;

  printf("Track Mult Eta Cut =  %2.2f\n",  fTrackMultEtaCut) ;

  printf("Track Mult Pt Cuts:") ;
  for(Int_t icut = 0; icut < fTrackMultNPtCut; icut++) printf(" %2.2f GeV;",fTrackMultPtCut[icut]);
  printf("    \n") ;
 
  printf("Write delta AOD =     %d\n",     fWriteOutputDeltaAOD) ;
  printf("Recalculate Clusters = %d, E linearity = %d\n",    fRecalculateClusters, fCorrectELinearity) ;
  
  printf("Use Triggers selected in SE base class %d; If not what Trigger Mask? %d; MB Trigger Mask for mixed %d \n",
         fEventTriggerAtSE, fEventTriggerMask,fMixEventTriggerMask);
  
  if(fComparePtHardAndJetPt)
    printf("Compare jet pt and pt hard to accept event, factor = %2.2f",fPtHardAndJetPtFactor);
  
  if(fComparePtHardAndClusterPt)
    printf("Compare cluster pt and pt hard to accept event, factor = %2.2f",fPtHardAndClusterPtFactor);
  
  printf("Read Kine from, stack? %d, AOD ? %d \n", fReadStack, fReadAODMCParticles) ;
  printf("Delta AOD File Name =     %s\n", fDeltaAODFileName.Data()) ;
  printf("Centrality: Class %s, Option %d, Bin [%d,%d] \n", fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1]) ;
  
  printf("    \n") ;
}

//__________________________________________
/// LED Events in period LHC11a contaminated 
/// EMCAL clusters sample, simple method
/// to reject such events.
//__________________________________________
Bool_t  AliCaloTrackReader::RejectLEDEvents()
{
  // Count number of cells with energy larger than 0.1 in SM3, cut on this number
  Int_t ncellsSM3 = 0;
  for(Int_t icell = 0; icell < fInputEvent->GetEMCALCells()->GetNumberOfCells(); icell++)
  {
    Int_t absID = fInputEvent->GetEMCALCells()->GetCellNumber(icell);
    Int_t sm    = GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absID);
    if(fInputEvent->GetEMCALCells()->GetAmplitude(icell) > 0.1 && sm==3) ncellsSM3++;
  }
  
  Int_t ncellcut = 21;
  if(GetFiredTriggerClasses().Contains("EMC")) ncellcut = 35;
  
  if(ncellsSM3 >= ncellcut)
  {
    AliDebug(1,Form("Reject event with ncells in SM3 %d, cut %d, trig %s",
                    ncellsSM3,ncellcut,GetFiredTriggerClasses().Data()));
    return kTRUE;
  }
  
  return kFALSE;
}

//_________________________________________________________
/// MC label for Cells not remapped after ESD filtering, do it here.
/// Needed for old MC/data productions done with AliRoot older 
/// than v5-02-Rev09 (more or less, not sure)
//_________________________________________________________
void AliCaloTrackReader::RemapMCLabelForAODs(Int_t & label)
{  
  if(label < 0) return ;
  
  AliAODEvent  * evt = dynamic_cast<AliAODEvent*> (fInputEvent) ;
  if(!evt) return ;
  
  TClonesArray * arr = dynamic_cast<TClonesArray*>(evt->FindListObject("mcparticles")) ;
  if(!arr) return ;
  
  if(label < arr->GetEntriesFast())
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(label));
    if(!particle) return ;
    
    if(label == particle->Label()) return ; // label already OK
    //else printf("AliCaloTrackReader::RemapMCLabelForAODs() - Label  %d - AOD stack %d \n",label, particle->Label());
  }
  //else printf("AliCaloTrackReader::RemapMCLabelForAODs() - Label  %d > AOD labels %d \n",label, arr->GetEntriesFast());
  
  // loop on the particles list and check if there is one with the same label
  for(Int_t ind = 0; ind < arr->GetEntriesFast(); ind++ )
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(ind));
    if(!particle) continue ;
    
    if(label == particle->Label())
    {
      label = ind;
      //printf("AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs() - New Label Index  %d \n",label);
      return;
    }
  }
  
  label = -1;
  
  //printf("AliCaloTrackReader::RemapMCLabelForAODs() - Label not found set to -1 \n");
}

//___________________________________
/// Reset lists, called in AliAnaCaloTrackCorrMaker.
//___________________________________
void AliCaloTrackReader::ResetLists()
{  
  if(fCTSTracks)       fCTSTracks     -> Clear();
  if(fEMCALClusters)   fEMCALClusters -> Clear("C");
  if(fPHOSClusters)    fPHOSClusters  -> Clear("C");
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0;
  fV0Mul[0] = 0;   fV0Mul[1] = 0;
  
  if(fNonStandardJets) fNonStandardJets -> Clear("C");
  fBackgroundJets->Reset();
}

//___________________________________________
/// Tag event depending on trigger name.
/// Set also the L1 bit defining the EGA or EJE triggers.
/// depending on the trigger class version, if not set by user.
//___________________________________________
void AliCaloTrackReader::SetEventTriggerBit()
{	
  fEventTrigMinBias       = kFALSE;
  fEventTrigCentral       = kFALSE;
  fEventTrigSemiCentral   = kFALSE;
  fEventTrigEMCALL0       = kFALSE;
  fEventTrigEMCALL1Gamma1 = kFALSE;
  fEventTrigEMCALL1Gamma2 = kFALSE;
  fEventTrigEMCALL1Jet1   = kFALSE;
  fEventTrigEMCALL1Jet2   = kFALSE;
  
  AliDebug(1,Form("Select trigger mask bit %d - Trigger Event %s - Select <%s>",
                  fEventTriggerMask,GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data()));
  
  if(fEventTriggerMask <=0 )// in case no mask set
  {
    // EMC triggered event? Which type?
    if( GetFiredTriggerClasses().Contains("-B-") || GetFiredTriggerClasses().Contains("-S-") || GetFiredTriggerClasses().Contains("-I-") )
    {
      if     ( GetFiredTriggerClasses().Contains("EGA" ) ||
              GetFiredTriggerClasses().Contains("EG1" )   )
      {
        fEventTrigEMCALL1Gamma1 = kTRUE;
        if( GetFiredTriggerClasses().Contains("EG1" ) && !fFiredTriggerClassName.Contains("EG1") ) fEventTrigEMCALL1Gamma1 = kFALSE;
      }
      else if( GetFiredTriggerClasses().Contains("EG2" )   )
      {
        fEventTrigEMCALL1Gamma2   = kTRUE;
        if( !fFiredTriggerClassName.Contains("EG2") ) fEventTrigEMCALL1Gamma2 = kFALSE;
      }
      else if( GetFiredTriggerClasses().Contains("EJE" ) ||
              GetFiredTriggerClasses().Contains("EJ1" )   )
      {
        fEventTrigEMCALL1Jet1   = kTRUE;
        if( GetFiredTriggerClasses().Contains("EJ1" ) && !fFiredTriggerClassName.Contains("EJ1") )
          fEventTrigEMCALL1Jet1 = kFALSE;
      }
      else if( GetFiredTriggerClasses().Contains("EJ2" )   )
      {
        fEventTrigEMCALL1Jet2   = kTRUE;
        if( !fFiredTriggerClassName.Contains("EJ2") ) fEventTrigEMCALL1Jet2 = kFALSE;
      }
      else if( GetFiredTriggerClasses().Contains("CEMC") &&
              !GetFiredTriggerClasses().Contains("EGA" ) &&
              !GetFiredTriggerClasses().Contains("EJE" ) &&
              !GetFiredTriggerClasses().Contains("EG1" ) &&
              !GetFiredTriggerClasses().Contains("EJ1" ) &&
              !GetFiredTriggerClasses().Contains("EG2" ) &&
              !GetFiredTriggerClasses().Contains("EJ2" )    )
      {
        fEventTrigEMCALL0 = kTRUE;
      }
      
      //Min bias event trigger?
      if     (GetFiredTriggerClasses().Contains("CCENT_R2-B-NOPF-ALLNOTRD"))
      {
        fEventTrigCentral     = kTRUE;
      }
      else if(GetFiredTriggerClasses().Contains("CSEMI_R1-B-NOPF-ALLNOTRD"))
      {
        fEventTrigSemiCentral = kTRUE;
      }
			else if((GetFiredTriggerClasses().Contains("CINT") || GetFiredTriggerClasses().Contains("CPBI2_B1") ) &&
              GetFiredTriggerClasses().Contains("-NOPF-ALLNOTRD") )
      {
			  fEventTrigMinBias = kTRUE;
      }
    }
  }
  else
	{
	  // EMC L1 Gamma
	  if     ( fEventTriggerMask & AliVEvent::kEMCEGA      )
    {
      //printf("EGA trigger bit\n");
      if     (GetFiredTriggerClasses().Contains("EG"))
      {
        if     (GetFiredTriggerClasses().Contains("EGA")) fEventTrigEMCALL1Gamma1 = kTRUE;
        else
        {
          if(GetFiredTriggerClasses().Contains("EG1")) fEventTrigEMCALL1Gamma1 = kTRUE;
          if(GetFiredTriggerClasses().Contains("EG2")) fEventTrigEMCALL1Gamma2 = kTRUE;
        }
      }
    }
	  // EMC L1 Jet
	  else if( fEventTriggerMask & AliVEvent::kEMCEJE      )
    {
      //printf("EGA trigger bit\n");
      if     (GetFiredTriggerClasses().Contains("EJ"))
      {
        if     (GetFiredTriggerClasses().Contains("EJE")) fEventTrigEMCALL1Jet1 = kTRUE;
        else
        {
          if(GetFiredTriggerClasses().Contains("EJ1")) fEventTrigEMCALL1Jet1 = kTRUE;
          if(GetFiredTriggerClasses().Contains("EJ2")) fEventTrigEMCALL1Jet2 = kTRUE;
        }
      }
    }
		// EMC L0
	  else if((fEventTriggerMask & AliVEvent::kEMC7) ||
            (fEventTriggerMask & AliVEvent::kEMC1)       )
    {
      //printf("L0 trigger bit\n");
	    fEventTrigEMCALL0 = kTRUE;
    }
	  // Min Bias Pb-Pb
	  else if( fEventTriggerMask & AliVEvent::kCentral     )
    {
      //printf("MB semi central trigger bit\n");
	    fEventTrigSemiCentral = kTRUE;
    }
	  // Min Bias Pb-Pb
	  else if( fEventTriggerMask & AliVEvent::kSemiCentral )
    {
      //printf("MB central trigger bit\n");
	    fEventTrigCentral = kTRUE;
    }
	  // Min Bias pp, PbPb, pPb
	  else if((fEventTriggerMask & AliVEvent::kMB  ) ||
            (fEventTriggerMask & AliVEvent::kINT7) ||
            (fEventTriggerMask & AliVEvent::kINT8) ||
            (fEventTriggerMask & AliVEvent::kAnyINT) )
    {
      //printf("MB trigger bit\n");
	    fEventTrigMinBias = kTRUE;
    }
	}
  
  AliDebug(1,Form("Event bits: \n \t MB   %d, Cen  %d, Sem  %d, L0   %d, L1G1 %d, L1G2 %d, L1J1 %d, L1J2 %d",
                  fEventTrigMinBias,      fEventTrigCentral,       fEventTrigSemiCentral,
                  fEventTrigEMCALL0 ,     fEventTrigEMCALL1Gamma1, fEventTrigEMCALL1Gamma2,
                  fEventTrigEMCALL1Jet1 , fEventTrigEMCALL1Jet2));
  
  // L1 trigger bit
  if( fBitEGA == 0 && fBitEJE == 0 )
  {
    // Init the trigger bit once, correct depending on AliESD(AOD)CaloTrigger header version
    
    // Simpler way to do it ...
//    if(  fInputEvent->GetRunNumber() < 172000 )
//      reader->SetEventTriggerL1Bit(4,5); // current LHC11 data
//    else
//      reader->SetEventTriggerL1Bit(6,8); // LHC12-13 data
    
    // Old values
    fBitEGA = 4;
    fBitEJE = 5;
        
    TFile* file = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
    
    const TList *clist = file->GetStreamerInfoCache();
    
    if(clist)
    {
      TStreamerInfo *cinfo = (TStreamerInfo*)clist->FindObject("AliESDCaloTrigger");
      Int_t verid = 5; // newer ESD header version
      if(!cinfo)
      {
        cinfo = (TStreamerInfo*)clist->FindObject("AliAODCaloTrigger");
        verid = 3; // newer AOD header version
      }
      
      if(cinfo)
	    {
	      Int_t classversionid = cinfo->GetClassVersion();
	      //printf("********* Header class version %d *********** \n",classversionid);
	      
        if (classversionid >= verid)
        {
          fBitEGA = 6;
          fBitEJE = 8;
        }
	    }  else AliInfo("AliCaloTrackReader()::SetEventTriggerBit() - Streamer info for trigger class not available, bit not changed");
    } else AliInfo("AliCaloTrackReader::SetEventTriggerBit() -  Streamer list not available!, bit not changed");
    
  } // set once the EJE, EGA trigger bit
}

//____________________________________________________________
/// Define here the input event and mixed event.
/// Called in AliAnaCaloTrackCorrMaker.
//____________________________________________________________
void AliCaloTrackReader::SetInputEvent(AliVEvent* const input)
{
  fInputEvent  = input;
  fMixedEvent = dynamic_cast<AliMixedEvent*>(GetInputEvent()) ;
  if (fMixedEvent)
    fNMixedEvent = fMixedEvent->GetNumberOfEvents() ;
  
  //Delete previous vertex
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++)
    {
      delete [] fVertex[i] ;
    }
    delete [] fVertex ;
  }
  
  fVertex = new Double_t*[fNMixedEvent] ;
  for (Int_t i = 0; i < fNMixedEvent; i++)
  {
    fVertex[i] = new Double_t[3] ;
    fVertex[i][0] = 0.0 ;
    fVertex[i][1] = 0.0 ;
    fVertex[i][2] = 0.0 ;
  }
}


