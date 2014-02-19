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

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and
// Central Barrel Tracking detectors (CTS).
// Not all MC particles/tracks/clusters are kept, some kinematical/fiducial restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader : Fills Kinematics data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS)
//-- Author: Gustavo Conesa (LNF-INFN)
//////////////////////////////////////////////////////////////////////////////


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
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliTriggerAnalysis.h"
#include "AliESDVZERO.h"
#include "AliVCaloCells.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"

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

ClassImp(AliCaloTrackReader)


//________________________________________
AliCaloTrackReader::AliCaloTrackReader() :
TObject(),                   fEventNumber(-1), //fCurrentFileName(""),
fDataType(0),                fDebug(0),
fFiducialCut(0x0),           fCheckFidCut(kFALSE),
fComparePtHardAndJetPt(0),   fPtHardAndJetPtFactor(0),
fComparePtHardAndClusterPt(0),fPtHardAndClusterPtFactor(0),
fCTSPtMin(0),                fEMCALPtMin(0),                  fPHOSPtMin(0),
fCTSPtMax(0),                fEMCALPtMax(0),                  fPHOSPtMax(0),
fUseEMCALTimeCut(1),         fUseParamTimeCut(0),             fUseTrackTimeCut(0),
fEMCALTimeCutMin(-10000),    fEMCALTimeCutMax(10000),
fEMCALParamTimeCutMin(),     fEMCALParamTimeCutMax(),
fTrackTimeCutMin(-10000),    fTrackTimeCutMax(10000),
fUseTrackDCACut(0),
fAODBranchList(0x0),
fCTSTracks(0x0),             fEMCALClusters(0x0),             fPHOSClusters(0x0),
fEMCALCells(0x0),            fPHOSCells(0x0),
fInputEvent(0x0),            fOutputEvent(0x0),fMC(0x0),
fFillCTS(0),                 fFillEMCAL(0),                   fFillPHOS(0),
fFillEMCALCells(0),          fFillPHOSCells(0),
fRecalculateClusters(kFALSE),fCorrectELinearity(kTRUE),
fSelectEmbeddedClusters(kFALSE),
fTrackStatus(0),             fTrackFilterMask(0),             fTrackFilterMaskComplementary(0),
fESDtrackCuts(0),            fESDtrackComplementaryCuts(0),   fConstrainTrack(kFALSE),
fSelectHybridTracks(0),      fSelectPrimaryTracks(0),
fSelectSPDHitTracks(0),      fSelectFractionTPCSharedClusters(0), fCutTPCSharedClustersFraction(0),
fTrackMult(0),               fTrackMultEtaCut(0.9),
fReadStack(kFALSE),          fReadAODMCParticles(kFALSE),
fDeltaAODFileName(""),       fFiredTriggerClassName(""),

fEventTriggerMask(0),        fMixEventTriggerMask(0),         fEventTriggerAtSE(0),
fEventTrigMinBias(0),        fEventTrigCentral(0),
fEventTrigSemiCentral(0),    fEventTrigEMCALL0(0),
fEventTrigEMCALL1Gamma1(0),  fEventTrigEMCALL1Gamma2(0),
fEventTrigEMCALL1Jet1(0),    fEventTrigEMCALL1Jet2(0),
fBitEGA(0),                  fBitEJE(0),

fAnaLED(kFALSE),
fTaskName(""),               fCaloUtils(0x0),
fMixedEvent(NULL),           fNMixedEvent(0),                 fVertex(NULL),
fListMixedTracksEvents(),    fListMixedCaloEvents(),
fLastMixedTracksEvent(-1),   fLastMixedCaloEvent(-1),
fWriteOutputDeltaAOD(kFALSE),fOldAOD(kFALSE),
fEMCALClustersListName(""),  fZvtxCut(0.),
fAcceptFastCluster(kFALSE),  fRemoveLEDEvents(kTRUE),
//Trigger rejection
fRemoveBadTriggerEvents(0),  fTriggerPatchClusterMatch(0),
fTriggerPatchTimeWindow(),   fTriggerL0EventThreshold(0),
fTriggerL1EventThreshold(0), fTriggerL1EventThresholdFix(0),
fTriggerClusterBC(0),        fTriggerClusterIndex(0),         fTriggerClusterId(0),
fIsExoticEvent(0),           fIsBadCellEvent(0),              fIsBadMaxCellEvent(0),
fIsTriggerMatch(0),          fIsTriggerMatchOpenCut(),
fTriggerClusterTimeRecal(kTRUE), fRemoveUnMatchedTriggers(kTRUE),
fDoEventSelection(kFALSE),   fDoV0ANDEventSelection(kFALSE),
fDoVertexBCEventSelection(kFALSE),
fDoRejectNoTrackEvents(kFALSE),
fUseEventsWithPrimaryVertex(kFALSE),
fTriggerAnalysis (0x0),      fTimeStampEventSelect(0),
fTimeStampEventFracMin(0),   fTimeStampEventFracMax(0),
fTimeStampRunMin(0),         fTimeStampRunMax(0),
fNPileUpClusters(-1),        fNNonPileUpClusters(-1),         fNPileUpClustersCut(3),
fVertexBC(-200),             fRecalculateVertexBC(0),
fCentralityClass(""),        fCentralityOpt(0),
fEventPlaneMethod(""),
fAcceptOnlyHIJINGLabels(0),  fNMCProducedMin(0), fNMCProducedMax(0),
fFillInputNonStandardJetBranch(kFALSE),
fNonStandardJets(new TClonesArray("AliAODJet",100)),fInputNonStandardJetBranchName("jets"),
fFillInputBackgroundJetBranch(kFALSE), 
fBackgroundJets(0x0),fInputBackgroundJetBranchName("jets"),
fAcceptEventsWithBit(0),     fRejectEventsWithBit(0), fRejectEMCalTriggerEventsWith2Tresholds(0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//_______________________________________
AliCaloTrackReader::~AliCaloTrackReader()
{
  //Dtor
  
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
  
  delete fESDtrackCuts;
  delete fESDtrackComplementaryCuts;
  delete fTriggerAnalysis;
  
  if(fNonStandardJets)
  {
    if(fDataType!=kMC) fNonStandardJets->Clear("C") ;
    else               fNonStandardJets->Delete() ;
    delete fNonStandardJets ;
  }
  delete fBackgroundJets ;

  fRejectEventsWithBit.Reset();
  fAcceptEventsWithBit.Reset();
  
  //  Pointers not owned, done by the analysis frame
  //  if(fInputEvent)  delete fInputEvent ;
  //  if(fOutputEvent) delete fOutputEvent ;
  //  if(fMC)          delete fMC ;
  //  Pointer not owned, deleted by maker
  //  if (fCaloUtils) delete fCaloUtils ;
  
}

//________________________________________________________________________
Bool_t  AliCaloTrackReader::AcceptDCA(Float_t pt, Float_t dca)
{
  // Accept track if DCA is smaller than function
  
  Float_t cut = fTrackDCACut[0]+fTrackDCACut[1]/TMath::Power(pt,fTrackDCACut[2]);
  
  if(TMath::Abs(dca) < cut)
    return kTRUE;
  else
    return kFALSE;
  
}

//_____________________________________________________
Bool_t  AliCaloTrackReader::AcceptEventWithTriggerBit()
{
  // Accept events that pass the physics selection
  // depending on an array of trigger bits set during the configuration
  
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
Bool_t  AliCaloTrackReader::RejectEventWithTriggerBit()
{
  // Reject events that pass the physics selection
  // depending on an array of trigger bits set during the configuration

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
Bool_t AliCaloTrackReader::CheckEventTriggers()
{
  // Do different selection of the event
  // depending on trigger name, event type, goodness of the EMCal trigger ...
  
  //-----------------------------------------------------------
  // Reject events depending on the trigger name and event type
  //-----------------------------------------------------------
  
  Int_t eventType = 0;
  if(fInputEvent->GetHeader())
  eventType = ((AliVHeader*)fInputEvent->GetHeader())->GetEventType();
  
  if( fFiredTriggerClassName  !="" && !fAnaLED)
  {
    //printf("Event type %d\n",eventType);
    if(eventType!=7)
    return kFALSE; //Only physics event, do not use for simulated events!!!
    
    if(fDebug > 0)
      printf("AliCaloTrackReader::CheckEventTriggers() - FiredTriggerClass <%s>, selected class <%s>, compare name %d\n",
             GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data(), GetFiredTriggerClasses().Contains(fFiredTriggerClassName));
    
    if( !GetFiredTriggerClasses().Contains(fFiredTriggerClassName) ) return kFALSE;
    else if(fDebug > 0) printf("AliCaloTrackReader::CheckEventTriggers() - Accepted triggered event\n");
  }
  else if(fAnaLED)
  {
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
    
	  if(eventType!=7 && fDebug > 1 )printf("AliCaloTrackReader::CheckEventTriggers() - DO LED, Event Type <%d>, 8 Calibration \n",  eventType);
	  if(eventType!=8)return kFALSE;
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::CheckEventTriggers() - Pass Trigger name rejection \n");

  //-------------------------------------------------------------------------------------
  // Reject or accept events depending on the trigger bit
  //-------------------------------------------------------------------------------------
  
  Bool_t okA = AcceptEventWithTriggerBit();
  Bool_t okR = RejectEventWithTriggerBit();
  
  //printf("AliCaloTrackReader::FillInputEvent() - Accept event? %d, Reject event %d? \n",okA,okR);
  
  if(!okA || !okR) return kFALSE;
  
  if(fDebug > 0) printf("AliCaloTrackReader::CheckEventTriggers() - Pass event bit rejection \n");
  
  //----------------------------------------------------------------------
  // Do not count events that were likely triggered by an exotic cluster
  // or out BC cluster
  //----------------------------------------------------------------------
  
  // Set a bit with the event kind, MB, L0, L1 ...
  SetEventTriggerBit();
  
  if( IsEventEMCALL1() || IsEventEMCALL0()  )
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
     if(fDebug > 0)
      printf("AliCaloTrackReader::CheckEventTriggers() - ACCEPT triggered event? \n exotic? %d - bad cell %d - bad Max cell %d - BC %d  - Matched %d\n",
             fIsExoticEvent,fIsBadCellEvent, fIsBadMaxCellEvent, fTriggerClusterBC,fIsTriggerMatch);
      if     (fIsExoticEvent)         return kFALSE;
      else if(fIsBadCellEvent)        return kFALSE;
      else if(fRemoveUnMatchedTriggers && !fIsTriggerMatch) return kFALSE ;
      else if(fTriggerClusterBC != 0) return kFALSE;
      if(fDebug > 0) printf("\t *** YES for %s\n",GetFiredTriggerClasses().Data());
    }
    
    if(fDebug > 0) printf("AliCaloTrackReader::CheckEventTriggers() - Pass EMCal triggered event rejection \n");
  }
  
  //-------------------------------------------------------------------------------------
  //Select events only fired by a certain trigger configuration if it is provided
  //-------------------------------------------------------------------------------------
  if (GetFiredTriggerClasses().Contains("FAST")  && !GetFiredTriggerClasses().Contains("ALL") && !fAcceptFastCluster)
  {
    if(fDebug > 0)  printf("AliCaloTrackReader::CheckEventTriggers() - Do not count events from fast cluster, trigger name %s\n",fFiredTriggerClassName.Data());
    return kFALSE;
  }
  
  //-------------------------------------------------------------------------------------
  // Reject event if large clusters with large energy
  // Use only for LHC11a data for the moment, and if input is clusterizer V1 or V1+unfolding
  // If clusterzer NxN or V2 it does not help
  //-------------------------------------------------------------------------------------
  Int_t run = fInputEvent->GetRunNumber();
  if( fRemoveLEDEvents && run > 146857  && run < 146861 )
  {
    Bool_t reject = RejectLEDEvents();

    if(reject) return kFALSE;
    
    if(fDebug > 0) printf("AliCaloTrackReader::CheckEventTriggers() - Pass LED event rejection \n");

  }// Remove LED events
  
  return kTRUE;
}

//________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndJetPt()
{
  // Check the event, if the requested ptHard is much smaller than the jet pT, then there is a problem.
  // Only for PYTHIA.
  
  //printf("AliCaloTrackReader::ComparePtHardAndJetPt() - GenHeaderName : %s\n",GetGenEventHeader()->ClassName());
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    TParticle * jet =  0;
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Int_t nTriggerJets =  pygeh->NTriggerJets();
    Float_t ptHard = pygeh->GetPtHard();
    
    if(fDebug > 1)
      printf("AliCaloTrackReader::ComparePtHardAndJetPt() - Njets: %d, pT Hard %f\n",nTriggerJets, ptHard);
    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      pygeh->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      
      if(fDebug > 1)
        printf("AliCaloTrackReader::ComparePtHardAndJetPt() - jet %d; pycell jet pT %f\n",ijet, jet->Pt());
      
      //Compare jet pT and pt Hard
      if(jet->Pt() > fPtHardAndJetPtFactor * ptHard)
      {
        printf("AliCaloTrackReader::ComparePtHardAndJetPt() - Reject jet event with : pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n",
               ptHard, jet->Pt(), fPtHardAndJetPtFactor);
        return kFALSE;
      }
    }
    
    if(jet) delete jet;
  }
  
  return kTRUE ;
  
}

//____________________________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndClusterPt()
{
  // Check the event, if the requested ptHard is smaller than the calorimeter cluster E, then there is a problem.
  // Only for PYTHIA.
  
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
        printf("AliCaloTrackReader::ComparePtHardAndClusterPt() - Reject : ecluster %2.2f, calo %d, factor %2.2f, ptHard %f\n",ecluster,clus->GetType(),fPtHardAndClusterPtFactor,ptHard);
        
        return kFALSE;
      }
    }
    
  }
  
  return kTRUE ;
  
}

//____________________________________________
AliStack* AliCaloTrackReader::GetStack() const
{
  //Return pointer to stack
  if(fMC)
    return fMC->Stack();
  else
  {
    if(fDebug > 1) printf("AliCaloTrackReader::GetStack() - Stack is not available\n");
    return 0x0 ;
  }
}

//______________________________________________
AliHeader* AliCaloTrackReader::GetHeader() const
{
  //Return pointer to header
  if(fMC)
  {
    return fMC->Header();
  }
  else
  {
    printf("AliCaloTrackReader::Header is not available\n");
    return 0x0 ;
  }
}

//____________________________________________________
void AliCaloTrackReader::SetGeneratorMinMaxParticles()
{
  // In case of access only to hijing particles in cocktail
  // get the min and max labels
  // TODO: Check when generator is not the first one ...
  
  fNMCProducedMin = 0;
  fNMCProducedMax = 0;
  
  if     (ReadStack() && fMC)
  {
    AliGenEventHeader * eventHeader = fMC->GenEventHeader();
    
    if(!fAcceptOnlyHIJINGLabels) return ;
    
    // TODO Check if it works from here ...
    
    AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    
    if(!cocktail) return ;
    
    TList *genHeaders = cocktail->GetHeaders();
    
    Int_t nGenerators = genHeaders->GetEntries();
    //printf("N generators %d \n", nGenerators);
    
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle());
      
      fNMCProducedMin = fNMCProducedMax;
      fNMCProducedMax+= eventHeader2->NProduced();
      
			if(name.Contains("Hijing",TString::kIgnoreCase)) return ;
    }
        
  }
  else if(ReadAODMCParticles() && GetAODMCHeader())
  {
    Int_t nGenerators = GetAODMCHeader()->GetNCocktailHeaders();
    //printf("AliCaloTrackReader::GetGenEventHeader() - N headers %d\n",nGenerators);
    
    if( nGenerators <= 0)        return ;
    
    if(!fAcceptOnlyHIJINGLabels) return ;
    
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader = GetAODMCHeader()->GetCocktailHeader(igen) ;
      TString name = eventHeader->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader->ClassName(), name.Data(), eventHeader->GetTitle());
      
      fNMCProducedMin = fNMCProducedMax;
      fNMCProducedMax+= eventHeader->NProduced();
      
			if(name.Contains("Hijing",TString::kIgnoreCase)) return ;
    }
        
  }
}


//______________________________________________________________
AliGenEventHeader* AliCaloTrackReader::GetGenEventHeader() const
{
  // Return pointer to Generated event header
  // If requested and cocktail, search for the hijing generator
  
  if     (ReadStack() && fMC)
  {
    AliGenEventHeader * eventHeader = fMC->GenEventHeader();
    
    if(!fAcceptOnlyHIJINGLabels) return eventHeader ;
    
    // TODO Check if it works from here ...
    
    AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    
    if(!cocktail) return 0x0 ;
    
    TList *genHeaders = cocktail->GetHeaders();
    
    Int_t nGenerators = genHeaders->GetEntries();
    //printf("N generators %d \n", nGenerators);
    
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle());
      
			if(name.Contains("Hijing",TString::kIgnoreCase)) return eventHeader2 ;
    }

    return 0x0;
    
  }
  else if(ReadAODMCParticles() && GetAODMCHeader())
  {
    Int_t nGenerators = GetAODMCHeader()->GetNCocktailHeaders();
    //printf("AliCaloTrackReader::GetGenEventHeader() - N headers %d\n",nGenerators);

    if( nGenerators <= 0)        return 0x0;
    
    if(!fAcceptOnlyHIJINGLabels) return GetAODMCHeader()->GetCocktailHeader(0);
        
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader = GetAODMCHeader()->GetCocktailHeader(igen) ;
      TString name = eventHeader->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader->ClassName(), name.Data(), eventHeader->GetTitle());
      
			if(name.Contains("Hijing",TString::kIgnoreCase)) return eventHeader ;
    }
    
    return 0x0;
        
  }
  else
  {
    //printf("AliCaloTrackReader::GetGenEventHeader() - MC header not available! \n");
    return 0x0;
  }
}

//____________________________________________________________________
TClonesArray* AliCaloTrackReader::GetAODMCParticles() const
{
  //Return list of particles in AOD. Do it for the corresponding input event.
  
  TClonesArray * rv = NULL ;
  if(fDataType == kAOD)
  {
    //Normal input AOD
    AliAODEvent * evt = dynamic_cast<AliAODEvent*> (fInputEvent) ;
    if(evt)
      rv = (TClonesArray*)evt->FindListObject("mcparticles");
    else
      printf("AliCaloTrackReader::GetAODMCParticles() - Null AOD event \n");
  }
  else
  {
    printf("AliCaloTrackReader::GetAODMCParticles() - Input are not AODs\n");
  }
  
  return rv ;
}

//________________________________________________________
AliAODMCHeader* AliCaloTrackReader::GetAODMCHeader() const
{
  //Return MC header in AOD. Do it for the corresponding input event.
  
  AliAODMCHeader *mch = NULL;
  
  if(fDataType == kAOD)
  {
    AliAODEvent * aod = dynamic_cast<AliAODEvent*> (fInputEvent);
    if(aod) mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
  }
  else
  {
    printf("AliCaloTrackReader::GetAODMCHeader() - Input are not AODs\n");
  }
  
  return mch;
}

//___________________________________________________________
Int_t AliCaloTrackReader::GetVertexBC(const AliVVertex * vtx)
{
  // Get the vertex BC
  
  Int_t vertexBC=vtx->GetBC();
  if(!fRecalculateVertexBC) return vertexBC;
  
  // In old AODs BC not stored, recalculate it
  // loop over the global track and select those which have small DCA to primary vertex (e.g. primary).
  // If at least one of these primaries has valid BC != 0, then this vertex is a pile-up candidate.
  // Execute after CTS
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
void AliCaloTrackReader::Init()
{
  //Init reader. Method to be called in AliAnaPartCorrMaker
  
  //printf(" AliCaloTrackReader::Init() %p \n",gGeoManager);
  
  if(fReadStack && fReadAODMCParticles)
  {
    printf("AliCaloTrackReader::Init() - Cannot access stack and mcparticles at the same time, change them \n");
    fReadStack          = kFALSE;
    fReadAODMCParticles = kFALSE;
  }

  if(!fESDtrackCuts)
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //initialize with TPC only tracks
	
}

//_______________________________________
void AliCaloTrackReader::InitParameters()
{
  //Initialize the parameters of the analysis.
  fDataType   = kESD ;
  fCTSPtMin   = 0.1 ;
  fEMCALPtMin = 0.1 ;
  fPHOSPtMin  = 0.1 ;
  fCTSPtMax   = 1000. ;
  fEMCALPtMax = 1000. ;
  fPHOSPtMax  = 1000. ;
  
  //Track DCA cuts
  // dca_xy cut = 0.0105+0.0350/TMath::Power(pt,1.1);
  fTrackDCACut[0] = 0.0105;
  fTrackDCACut[1] = 0.0350;
  fTrackDCACut[2] = 1.1;
  
  //Do not filter the detectors input by default.
  fFillEMCAL      = kFALSE;
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
  fAnaLED            = kFALSE;
  
  //We want tracks fitted in the detectors:
  //fTrackStatus=AliESDtrack::kTPCrefit;
  //fTrackStatus|=AliESDtrack::kITSrefit;
  fTrackStatus     = 0;
  fTrackFilterMask = 128; //For AODs, but what is the difference between fTrackStatus and fTrackFilterMask?
  fTrackFilterMaskComplementary = 0; // in case of hybrid tracks, without using the standard method
  
  fSelectFractionTPCSharedClusters = kTRUE;
  fCutTPCSharedClustersFraction = 0.4,
  
  fESDtrackCuts = 0;
  fESDtrackComplementaryCuts = 0;
  
  fConstrainTrack = kFALSE ; // constrain tracks to vertex
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0;
  fV0Mul[0] = 0;   fV0Mul[1] = 0;
  
  fZvtxCut   = 10.;
  
  fNMixedEvent = 1;
  
  fPtHardAndJetPtFactor     = 7.;
  fPtHardAndClusterPtFactor = 1.;
  
  //Centrality
  fCentralityClass  = "V0M";
  fCentralityOpt    = 10;
  fCentralityBin[0] = fCentralityBin[1]=-1;
  
  fEventPlaneMethod = "V0";
  
  // Allocate memory (not sure this is the right place)
  fCTSTracks       = new TObjArray();
  fEMCALClusters   = new TObjArray();
  fPHOSClusters    = new TObjArray();
  fTriggerAnalysis = new AliTriggerAnalysis;
  fAODBranchList   = new TList ;
  
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
  fTriggerL0EventThreshold =  2.;
  fTriggerL1EventThreshold =  4.;
  fTriggerClusterIndex     = -1;
  fTriggerClusterId        = -1;
  
  //Jets
  fInputNonStandardJetBranchName = "jets";
  fFillInputNonStandardJetBranch = kFALSE;
  if(!fNonStandardJets) fNonStandardJets = new TClonesArray("AliAODJet",100);
  fInputBackgroundJetBranchName = "jets";
  fFillInputBackgroundJetBranch = kFALSE; 
  if(!fBackgroundJets) fBackgroundJets = new AliAODJetEventBackground();

}

//___________________________________________________________________
Bool_t AliCaloTrackReader::IsEMCALCluster(AliVCluster* cluster) const
{
  // Check if it is a cluster from EMCAL. For old AODs cluster type has
  // different number and need to patch here
  
  if(fDataType==kAOD && fOldAOD)
  {
    if (cluster->GetType() == 2) return kTRUE;
    else                         return kFALSE;
  }
  else
  {
    return cluster->IsEMCAL();
  }
  
}

//___________________________________________________________________
Bool_t AliCaloTrackReader::IsPHOSCluster(AliVCluster * cluster) const
{
  //Check if it is a cluster from PHOS.For old AODs cluster type has
  // different number and need to patch here
  
  if(fDataType==kAOD && fOldAOD)
  {
    Int_t type = cluster->GetType();
    if (type == 0 || type == 1) return kTRUE;
    else                        return kFALSE;
  }
  else
  {
    return cluster->IsPHOS();
  }
  
}

//________________________________________________________________________
Bool_t AliCaloTrackReader::IsHIJINGLabel(Int_t label)
{
 
  // Find if cluster/track was generated by HIJING
  
  AliGenHijingEventHeader*  hijingHeader =  dynamic_cast<AliGenHijingEventHeader *> (GetGenEventHeader());
  
  //printf("header %p, label %d\n",hijingHeader,label);
  
  if(!hijingHeader || label < 0 ) return kFALSE;
  
  
  //printf("pass a), N produced %d\n",nproduced);
  
  if(label >= fNMCProducedMin && label < fNMCProducedMax)
  {
    //printf(" accept!, label is smaller than produced, N %d\n",nproduced);

    return kTRUE;
  }
  
  if(ReadStack())
  {
    if(!GetStack()) return kFALSE;
    
    Int_t nprimaries = GetStack()->GetNtrack();
    
    if(label > nprimaries) return kFALSE;
    
    TParticle * mom = GetStack()->Particle(label);
    
    Int_t iMom = label;
    Int_t iParent = mom->GetFirstMother();
    while(iParent!=-1)
    {
      if(iParent >= fNMCProducedMin && iParent < fNMCProducedMax)
      {
        //printf("\t accept, mother is %d \n",iParent)
        return kTRUE;
      }
      
      iMom = iParent;
      mom = GetStack()->Particle(iMom);
      iParent = mom->GetFirstMother();
    }
    
    return kFALSE ;
    
  } // ESD
  else
  {
    TClonesArray* mcparticles = GetAODMCParticles();
    
    if(!mcparticles) return kFALSE;
    
    Int_t nprimaries = mcparticles->GetEntriesFast();
    
    if(label > nprimaries) return kFALSE;
    
    //printf("pass b) N primaries %d \n",nprimaries);
    
    if(label >= fNMCProducedMin && label < fNMCProducedMax)
    {
      return kTRUE;
    }
    
    // Find grand parent, check if produced in the good range
    AliAODMCParticle * mom = (AliAODMCParticle *) mcparticles->At(label);
    
    Int_t iMom = label;
    Int_t iParent = mom->GetMother();
    while(iParent!=-1)
    {
      if(iParent >= fNMCProducedMin && iParent < fNMCProducedMax)
      {
        //printf("\t accept, mother is %d, with nProduced %d \n",iParent, nproduced);
        return kTRUE;
      }
      
      iMom = iParent;
      mom = (AliAODMCParticle *) mcparticles->At(iMom);
      iParent = mom->GetMother();

    }
    
    //printf("pass c), no match found \n");
    
    return kFALSE ;
    
  }//AOD
}

//__________________________________________________________________________
Bool_t AliCaloTrackReader::IsInTimeWindow(Double_t tof, Float_t energy) const
{
  // Cluster time selection window
  
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
Bool_t AliCaloTrackReader::IsPileUpFromSPD() const
{
  // Check if event is from pile-up determined by SPD
  // Default values: (3, 0.8, 3., 2., 5.)
  return fInputEvent->IsPileupFromSPD((Int_t) fPileUpParamSPD[0] , fPileUpParamSPD[1] ,
                                      fPileUpParamSPD[2] , fPileUpParamSPD[3] , fPileUpParamSPD[4] );
  //printf("Param : %d, %2.2f, %2.2f, %2.2f, %2.2f\n",(Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
  
}

//__________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCal() const
{
  // Check if event is from pile-up determined by EMCal
  if(fNPileUpClusters > fNPileUpClustersCut) return kTRUE ;
  else                                       return kFALSE;
}

//________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndEMCal() const
{
  // Check if event is from pile-up determined by SPD and EMCal
  if( IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//_______________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDOrEMCal() const
{
  // Check if event is from pile-up determined by SPD or EMCal
  if( IsPileUpFromSPD() || IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndNotEMCal() const
{
  // Check if event is from pile-up determined by SPD and not by EMCal
  if( IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCalAndNotSPD() const
{
  // Check if event is from pile-up determined by EMCal, not by SPD
  if( !IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                           return kFALSE;
}

//______________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromNotSPDAndNotEMCal() const
{
  // Check if event not from pile-up determined neither by SPD nor by EMCal
  if( !IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                            return kFALSE;
}

//___________________________________________________________________________________
Bool_t AliCaloTrackReader::FillInputEvent(Int_t iEntry, const char * /*curFileName*/)
{
  //Fill the event counter and input lists that are needed, called by the analysis maker.
  
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
    if(fDebug >= 0) printf("AliCaloTrackReader::FillInputEvent() - Input event not available, skip event analysis\n");
    return kFALSE;
  }
  
  Bool_t accept = CheckEventTriggers();
  if(!accept) return kFALSE;
  
  //---------------------------------------------------------------------------
  // In case of analysis of events with jets, skip those with jet pt > 5 pt hard
  // To be used on for MC data in pT hard bins
  //---------------------------------------------------------------------------
  if(fComparePtHardAndJetPt)
  {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
  }
  
  if(fComparePtHardAndClusterPt)
  {
    if(!ComparePtHardAndClusterPt()) return kFALSE ;
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass Pt Hard rejection \n");
  
  //Fill Vertex array
  FillVertexArray();
  //Reject events with Z vertex too large, only for SE analysis, if not, cut on the analysis code
  if(!GetMixedEvent() && TMath::Abs(fVertex[0][2]) > fZvtxCut) return kFALSE;
  
  //------------------------------------------------------
  //Event rejection depending on vertex, pileup, v0and
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
    //printf("\t accept time stamp\n");
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass Time Stamp rejection \n");
  
  //------------------------------------------------------
  //Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------
  
  if(fUseEventsWithPrimaryVertex)
  {
    if( !CheckForPrimaryVertex() )              return kFALSE;
    if( TMath::Abs(fVertex[0][0] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][1] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][2] ) < 1.e-6    ) return kFALSE;
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass primary vertex rejection \n");
  
  //printf("Reader : IsPileUp %d, Multi %d\n",IsPileUpFromSPD(),fInputEvent->IsPileupFromSPDInMultBins());
  
  if(fDoEventSelection)
  {
    // Do not analyze events with pileup
    Bool_t bPileup = IsPileUpFromSPD();
    //IsPileupFromSPDInMultBins() // method to try
    //printf("pile-up %d, %d, %2.2f, %2.2f, %2.2f, %2.2f\n",bPileup, (Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
    if(bPileup) return kFALSE;
    
    if(fDoV0ANDEventSelection)
    {
      Bool_t bV0AND = kTRUE;
      AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
      if(esd)
        bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0AND);
      //else bV0AND = //FIXME FOR AODs
      if(!bV0AND) return kFALSE;
    }
  }// Event selection/AliceSoft/AliRoot/trunk/PWG/CaloTrackCorrBase/AliCaloTrackReader.h
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass Pile-Up, V0AND event rejection \n");

  //------------------------------------------------------
  
  //Check if there is a centrality value, PbPb analysis, and if a centrality bin selection is requested
  //If we need a centrality bin, we select only those events in the corresponding bin.
  if(GetCentrality() && fCentralityBin[0]>=0 && fCentralityBin[1]>=0 && fCentralityOpt==100)
  {
    Int_t cen = GetEventCentrality();
    if(cen > fCentralityBin[1] || cen < fCentralityBin[0]) return kFALSE; //reject events out of bin.
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass centrality rejection \n");

  
  //Fill the arrays with cluster/tracks/cells data
  
  if(!fEventTriggerAtSE)
  {
    // In case of mixing analysis, accept MB events, not only Trigger
    // Track and cluster arrays filled for MB in order to create the pool in the corresponding analysis
    // via de method in the base class FillMixedEventPool()
    
    AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
    
    if(!inputHandler) return kFALSE ;  // to content coverity
    
    UInt_t isTrigger = inputHandler->IsEventSelected() & fEventTriggerMask;
    UInt_t isMB      = inputHandler->IsEventSelected() & fMixEventTriggerMask;
    
    if(!isTrigger && !isMB) return kFALSE;
    
    //printf("Selected triggered event : %s\n",GetFiredTriggerClasses().Data());
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass uninteresting triggered events rejection in case of mixing analysis \n");

  
  // Get the main vertex BC, in case not available
  // it is calculated in FillCTS checking the BC of tracks
  // with DCA small (if cut applied, if open)
  fVertexBC=fInputEvent->GetPrimaryVertex()->GetBC();
  
  if(fAcceptOnlyHIJINGLabels) SetGeneratorMinMaxParticles();
  
  //printf("N min %d, N max %d\n",fNMCProducedMin,fNMCProducedMax);
  
  if(fFillCTS)
  {
    FillInputCTS();
    //Accept events with at least one track
    if(fTrackMult == 0 && fDoRejectNoTrackEvents) return kFALSE ;
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass rejection of null track events \n");

  
  if(fDoVertexBCEventSelection)
  {
    if(fVertexBC!=0 && fVertexBC!=AliVTrack::kTOFBCNA) return kFALSE ;
  }
  
  if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent()-Pass rejection of events with vertex at BC!=0 \n");

  
  if(fFillEMCALCells)
    FillInputEMCALCells();
  
  if(fFillPHOSCells)
    FillInputPHOSCells();
  
  if(fFillEMCAL)
    FillInputEMCAL();
  
  if(fFillPHOS)
    FillInputPHOS();
  
  //FillInputVZERO();
  
  //one specified jet branch
  if(fFillInputNonStandardJetBranch)
    FillInputNonStandardJets();
  if(fFillInputBackgroundJetBranch)
    FillInputBackgroundJets();

  
  return kTRUE ;
}

//__________________________________________________
Int_t AliCaloTrackReader::GetEventCentrality() const
{
  //Return current event centrality
  
  if(GetCentrality())
  {
    if     (fCentralityOpt==100) return (Int_t) GetCentrality()->GetCentralityPercentile(fCentralityClass); // 100 bins max
    else if(fCentralityOpt==10)  return GetCentrality()->GetCentralityClass10(fCentralityClass);// 10 bins max
    else if(fCentralityOpt==20)  return GetCentrality()->GetCentralityClass5(fCentralityClass); // 20 bins max
    else
    {
      printf("AliCaloTrackReader::GetEventCentrality() - Unknown centrality option %d, use 10, 20 or 100\n",fCentralityOpt);
      return -1;
    }
  }
  else return -1;
  
}

//_____________________________________________________
Double_t AliCaloTrackReader::GetEventPlaneAngle() const
{
  //Return current event centrality
  
  if(GetEventPlane())
  {
    Float_t ep =  GetEventPlane()->GetEventplane(GetEventPlaneMethod(), GetInputEvent());
    
    if(GetEventPlaneMethod()=="Q" && (ep < 0 || ep > TMath::Pi()))
    {
      if(fDebug > 0 ) printf("AliCaloTrackReader::GetEventPlaneAngle() -  Bad EP for <Q> method : %f\n",ep);
      return -1000;
    }
    else if(GetEventPlaneMethod().Contains("V0")  )
    {
      if((ep > TMath::Pi()/2 || ep < -TMath::Pi()/2))
      {
        if(fDebug > 0 ) printf("AliCaloTrackReader::GetEventPlaneAngle() -  Bad EP for <%s> method : %f\n",GetEventPlaneMethod().Data(), ep);
        return -1000;
      }
      
      ep+=TMath::Pi()/2; // put same range as for <Q> method
      
    }
    
    //printf("AliCaloTrackReader::GetEventPlaneAngle() = %f\n",ep);
    if(fDebug > 0 )
    {
      if     (ep > TMath::Pi()) printf("AliCaloTrackReader::GetEventPlaneAngle() - Too large angle = %f\n",ep);
      else if(ep < 0          ) printf("AliCaloTrackReader::GetEventPlaneAngle() - Negative angle = %f\n" ,ep);
    }
    
    return ep;
  }
  else
  {
    if(fDataType!=kMC && fDebug > 0) printf("AliCaloTrackReader::GetEventPlaneAngle() -  No EP pointer\n");
    return -1000;
  }
  
}

//__________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3]) const
{
  //Return vertex position to be used for single event analysis
  vertex[0]=fVertex[0][0];
  vertex[1]=fVertex[0][1];
  vertex[2]=fVertex[0][2];
}

//__________________________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3], Int_t evtIndex) const
{
  //Return vertex position for mixed event, recover the vertex in a particular event.
  
  vertex[0]=fVertex[evtIndex][0];  vertex[1]=fVertex[evtIndex][1];  vertex[2]=fVertex[evtIndex][2];
  
}

//________________________________________
void AliCaloTrackReader::FillVertexArray()
{
  
  //Fill data member with vertex
  //In case of Mixed event, multiple vertices
  
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
  
  if (!fMixedEvent)
  { //Single event analysis
    if(fDataType!=kMC)
    {
      
      if(fInputEvent->GetPrimaryVertex())
      {
        fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex[0]);
      }
      else
      {
        printf("AliCaloTrackReader::FillVertexArray() - NULL primary vertex\n");
        fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
      }//Primary vertex pointer do not exist
      
    } else
    {//MC read event
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
    
    if(fDebug > 1)
      printf("AliCaloTrackReader::FillVertexArray() - Single Event Vertex : %f,%f,%f\n",fVertex[0][0],fVertex[0][1],fVertex[0][2]);
    
  } else
  { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++)
    {
      if (fMixedEvent->GetVertexOfEvent(iev))
        fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
      else
      { // no vertex found !!!!
        AliWarning("No vertex found");
      }
      
      if(fDebug > 1)
        printf("AliCaloTrackReader::FillVertexArray() - Multi Event %d Vertex : %f,%f,%f\n",iev,fVertex[iev][0],fVertex[iev][1],fVertex[iev][2]);
      
    }
  }
  
}

//_____________________________________
void AliCaloTrackReader::FillInputCTS()
{
  //Return array with Central Tracking System (CTS) tracks
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS()\n");
  
  Double_t pTrack[3] = {0,0,0};
  
  Int_t nTracks = fInputEvent->GetNumberOfTracks() ;
  fTrackMult    = 0;
  Int_t nstatus = 0;
  Double_t bz   = GetInputEvent()->GetMagneticField();
  
  for(Int_t i = 0; i < 19; i++)
  {
    fTrackBCEvent   [i] = 0;
    fTrackBCEventCut[i] = 0;
  }
  
  Bool_t   bc0  = kFALSE;
  if(fRecalculateVertexBC) fVertexBC=AliVTrack::kTOFBCNA;
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++)
  {////////////// track loop
    AliVTrack * track = (AliVTrack*)fInputEvent->GetTrack(itrack) ; // retrieve track from esd
    
    //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
    ULong_t status = track->GetStatus();
    
    if (fTrackStatus && !((status & fTrackStatus) == fTrackStatus))
      continue ;
    
    nstatus++;
    
    Float_t dcaTPC =-999;
    
    if     (fDataType==kESD)
    {
      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
      
      if(esdTrack)
      {
        if(fESDtrackCuts->AcceptTrack(esdTrack))
        {
          track->GetPxPyPz(pTrack) ;
          
          if(fConstrainTrack)
          {
            if(esdTrack->GetConstrainedParam())
            {
              const AliExternalTrackParam* constrainParam = esdTrack->GetConstrainedParam();
              esdTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
              esdTrack->GetConstrainedPxPyPz(pTrack);
            }
            else continue;
            
          } // use constrained tracks
          
          if(fSelectSPDHitTracks)
          {//Not much sense to use with TPC only or Hybrid tracks
            if(!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1)) continue ;
          }
        }
        // Complementary track to global : Hybrids (make sure that the previous selection is for Global)
        else  if(fESDtrackComplementaryCuts && fESDtrackComplementaryCuts->AcceptTrack(esdTrack))
        {
          // constrain the track
          if(esdTrack->GetConstrainedParam())
          {
            esdTrack->Set(esdTrack->GetConstrainedParam()->GetX(),esdTrack->GetConstrainedParam()->GetAlpha(),esdTrack->GetConstrainedParam()->GetParameter(),esdTrack->GetConstrainedParam()->GetCovariance());
            
            track->GetPxPyPz(pTrack) ;
            
          }
          else continue;
        }
        else continue;
      }
    } // ESD
    else if(fDataType==kAOD)
    {
      AliAODTrack *aodtrack = dynamic_cast <AliAODTrack*>(track);
      
      if(aodtrack)
      {
        if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS():AOD track type: %d (primary %d), hybrid? %d \n",
                               aodtrack->GetType(),AliAODTrack::kPrimary,
                               aodtrack->IsHybridGlobalConstrainedGlobal());
        
        if (fSelectHybridTracks && fTrackFilterMaskComplementary == 0)
        {
          if (!aodtrack->IsHybridGlobalConstrainedGlobal())       continue ;
        }
        else
        {
          Bool_t accept = aodtrack->TestFilterBit(fTrackFilterMask);
          
          if(!fSelectHybridTracks && !accept) continue ;
          
          if(fSelectHybridTracks)
          {
            Bool_t acceptcomplement = aodtrack->TestFilterBit(fTrackFilterMaskComplementary);
            if (!accept && !acceptcomplement) continue ;
          }
        }
        
        if(fSelectSPDHitTracks)
        {//Not much sense to use with TPC only or Hybrid tracks
          if(!aodtrack->HasPointOnITSLayer(0) && !aodtrack->HasPointOnITSLayer(1)) continue ;
        }
        
        if ( fSelectFractionTPCSharedClusters )
        {
          Double_t frac = Double_t(aodtrack->GetTPCnclsS()) / Double_t(aodtrack->GetTPCncls());
          if (frac > fCutTPCSharedClustersFraction)
          {
             if (fDebug > 2 )printf("\t Reject track, shared cluster fraction %f > %f\n",frac, fCutTPCSharedClustersFraction);
            continue ;
          }
        }
        
        if ( fSelectPrimaryTracks )
        {
          if ( aodtrack->GetType()!= AliAODTrack::kPrimary )
          {
             if (fDebug > 2 ) printf("\t Remove not primary track\n");
            continue ;
          }
        }

        if (fDebug > 2 ) printf("\t accepted track! \n");
        
        //In case of AODs, TPC tracks cannot be propagated back to primary vertex,
        // info stored here
        dcaTPC = aodtrack->DCA();
        
        track->GetPxPyPz(pTrack) ;
        
      } // aod track exists
      else continue ;
      
    } // AOD
    
    TLorentzVector momentum(pTrack[0],pTrack[1],pTrack[2],0);
    
    Bool_t okTOF  = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
    Double_t tof  = -1000;
    Int_t trackBC = -1000 ;
    
    if(okTOF)
    {
      trackBC = track->GetTOFBunchCrossing(bz);
      SetTrackEventBC(trackBC+9);
      
      tof = track->GetTOFsignal()*1e-3;
    }
    
    if(fUseTrackDCACut)
    {
      //normal way to get the dca, cut on dca_xy
      if(dcaTPC==-999)
      {
        Double_t dca[2]   = {1e6,1e6};
        Double_t covar[3] = {1e6,1e6,1e6};
        Bool_t okDCA = track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),bz,100.,dca,covar);
        if( okDCA) okDCA = AcceptDCA(momentum.Pt(),dca[0]);
        if(!okDCA)
        {
          //printf("AliCaloTrackReader::FillInputCTS() - Reject track pt %2.2f, dca_xy %2.4f, BC %d\n",momentum.Pt(),dca[0],trackBC);
          continue ;
        }
      }
    }// DCA cuts
    
    if(okTOF)
    {
      //SetTrackEventBCcut(bc);
      SetTrackEventBCcut(trackBC+9);
      
      //After selecting tracks with small DCA, pointing to vertex, set vertex BC depeding on tracks BC
      if(fRecalculateVertexBC)
      {
        if     (trackBC !=0 && trackBC != AliVTrack::kTOFBCNA) fVertexBC = trackBC;
        else if(trackBC == 0)                                  bc0 = kTRUE;
      }
      
      //In any case, the time should to be larger than the fixed window ...
      if( fUseTrackTimeCut && (trackBC!=0 || tof < fTrackTimeCutMin  || tof > fTrackTimeCutMax) )
      {
        //printf("Remove track time %f and bc = %d\n",tof,trackBC);
        continue ;
      }
      //else printf("Accept track time %f and bc = %d\n",tof,trackBC);
      
    }

    //Count the tracks in eta < 0.9
    //printf("Eta %f cut  %f\n",TMath::Abs(track->Eta()),fTrackMultEtaCut);
    if(TMath::Abs(track->Eta())< fTrackMultEtaCut) fTrackMult++;
    
    if(fCTSPtMin > momentum.Pt() || fCTSPtMax < momentum.Pt()) continue ;
    
    if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"CTS")) continue;
    
    if(fDebug > 2 && momentum.Pt() > 0.1)
      printf("AliCaloTrackReader::FillInputCTS() - Selected tracks pt %3.2f, phi %3.2f, eta %3.2f\n",
             momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    
    if (fMixedEvent)  track->SetID(itrack);
    
    if(fAcceptOnlyHIJINGLabels && !IsHIJINGLabel(track->GetLabel())) continue ;
    
    fCTSTracks->Add(track);
    
  }// track loop
	
  if(fVertexBC ==0 || fVertexBC == AliVTrack::kTOFBCNA)
  {
    if( bc0 ) fVertexBC = 0 ;
    else      fVertexBC = AliVTrack::kTOFBCNA ;
  }
  
  
  if(fDebug > 1)
    printf("AliCaloTrackReader::FillInputCTS()   - aod entries %d, input tracks %d, pass status %d, multipliticy %d\n", fCTSTracks->GetEntriesFast(), nTracks, nstatus, fTrackMult);//fCTSTracksNormalInputEntries);
  
}

//_______________________________________________________________________________
void AliCaloTrackReader::FillInputEMCALAlgorithm(AliVCluster * clus, Int_t iclus)
{
  //Fill the EMCAL data in the array, do it
  
  Int_t vindex = 0 ;
  if (fMixedEvent)
    vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
  
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
      
      if(fDataType==AliCaloTrackReader::kESD)
      {
        tof = fEMCALCells->GetCellTime(absIdMax);
      }
      
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
      
      clus->SetTOF(tof);
      
    }// Time recalibration
  }
  
  //Reject clusters with bad channels, close to borders and exotic;
  if(!GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,GetCaloUtils()->GetEMCALGeometry(),GetEMCALCells(),fInputEvent->GetBunchCrossNumber())) return;
  
  //Mask all cells in collumns facing ALICE thick material if requested
  if(GetCaloUtils()->GetNMaskCellColumns())
  {
    Int_t absId   = -1;
    Int_t iSupMod = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Bool_t shared = kFALSE;
    GetCaloUtils()->GetEMCALRecoUtils()->GetMaxEnergyCell(GetCaloUtils()->GetEMCALGeometry(), GetEMCALCells(),clus,absId,iSupMod,ieta,iphi,shared);
    if(GetCaloUtils()->MaskFrameCluster(iSupMod, ieta)) return;
  }
  
  if(fSelectEmbeddedClusters)
  {
    if(clus->GetNLabels()==0 || clus->GetLabel() < 0) return;
    //else printf("Embedded cluster,  %d, n label %d label %d  \n",iclus,clus->GetNLabels(),clus->GetLabel());
  }
  
  //Float_t pos[3];
  //clus->GetPosition(pos);
  //printf("Before Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
  
  //Correct non linearity
  if(fCorrectELinearity && GetCaloUtils()->IsCorrectionOfClusterEnergyOn())
  {
    GetCaloUtils()->CorrectClusterEnergy(clus) ;
    
    //In case of MC analysis, to match resolution/calibration in real data
    Float_t rdmEnergy = GetCaloUtils()->GetEMCALRecoUtils()->SmearClusterEnergy(clus);
    // printf("\t Energy %f, smeared %f\n", clus->E(),rdmEnergy);
    clus->SetE(rdmEnergy);
  }
    
  Double_t tof = clus->GetTOF()*1e9;
  
  Int_t bc = TMath::Nint(tof/50) + 9;
  //printf("tof %2.2f, bc+5=%d\n",tof,bc);
  
  SetEMCalEventBC(bc);
  
  if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) return ;
  
  TLorentzVector momentum ;
  
  clus->GetMomentum(momentum, fVertex[vindex]);
  
  if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) return ;
  
  SetEMCalEventBCcut(bc);
  
  if(!IsInTimeWindow(tof,clus->E()))
  {
    fNPileUpClusters++ ;
    if(fUseEMCALTimeCut) return ;
  }
  else
    fNNonPileUpClusters++;
  
  if(fDebug > 2 && momentum.E() > 0.1)
    printf("AliCaloTrackReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
           momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
  
  if (fMixedEvent)
    clus->SetID(iclus) ;
  
  if(fAcceptOnlyHIJINGLabels && !IsHIJINGLabel( clus->GetLabel() )) return ;
  
//  if(fAcceptOnlyHIJINGLabels)
//  {
//    printf("Accept label %d?\n",clus->GetLabel());
//
//    if( !IsHIJINGLabel( clus->GetLabel() ) ) { printf("\t Reject label\n") ; return ; }
//    else                                       printf("\t Accept label\n") ;
//  }
  
  fEMCALClusters->Add(clus);
  
}

//_______________________________________
void AliCaloTrackReader::FillInputEMCAL()
{
  //Return array with EMCAL clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputEMCAL()\n");
  
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
        if (IsEMCALCluster(clus))
        {
          FillInputEMCALAlgorithm(clus, iclus);
        }//EMCAL cluster
      }// cluster exists
    }// cluster loop
    
    //Recalculate track matching
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent);
    
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
      printf("AliCaloTrackReader::FillInputEMCAL() - Wrong name of list with clusters?  <%s>\n",fEMCALClustersListName.Data());
      return;
    }
    
    Int_t nclusters = clusterList->GetEntriesFast();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
      //printf("E %f\n",clus->E());
      if (clus) FillInputEMCALAlgorithm(clus, iclus);
      else printf("AliCaloTrackReader::FillInputEMCAL() - Null cluster in list!\n");
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
        if (IsEMCALCluster(clus))
        {
          
          Float_t  frac     =-1;
          Int_t    absIdMax = GetCaloUtils()->GetMaxEnergyCell(fEMCALCells, clus,frac);
          Double_t tof = clus->GetTOF();
          GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
          tof*=1e9;
          
          //printf("Input event cluster : AbsIdMax %d, E %2.2f, time %2.2f \n", absIdMax,clus->E(),tof);
          
          //Reject clusters with bad channels, close to borders and exotic;
          if(!GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,GetCaloUtils()->GetEMCALGeometry(),GetEMCALCells(),fInputEvent->GetBunchCrossNumber()))  continue;
          
          Int_t bc = TMath::Nint(tof/50) + 9;
          SetEMCalEventBC(bc);
          
          if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) continue ;
          
          TLorentzVector momentum ;
          
          clus->GetMomentum(momentum, fVertex[0]);
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) return ;
          
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
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent,clusterList);
    
  }
  
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod entries %d, n pile-up clusters %d, n non pile-up %d \n",  fEMCALClusters->GetEntriesFast(),fNPileUpClusters,fNNonPileUpClusters);
  
}

//______________________________________
void AliCaloTrackReader::FillInputPHOS()
{
  //Return array with PHOS clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputPHOS()\n");
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
  for (Int_t iclus = 0; iclus < nclusters; iclus++)
  {
    AliVCluster * clus = 0;
    if ( (clus = fInputEvent->GetCaloCluster(iclus)) )
    {
      if (IsPHOSCluster(clus))
      {
        //Check if the cluster contains any bad channel and if close to calorimeter borders
        Int_t vindex = 0 ;
        if (fMixedEvent)
          vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
        if( GetCaloUtils()->ClusterContainsBadChannel("PHOS",clus->GetCellsAbsId(), clus->GetNCells()))
          continue;
        if(!GetCaloUtils()->CheckCellFiducialRegion(clus, fInputEvent->GetPHOSCells(), fInputEvent, vindex))
          continue;
        
        if(fRecalculateClusters)
        {
          //Recalibrate the cluster energy
          if(GetCaloUtils()->IsRecalibrationOn())
          {
            Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliAODCaloCells*)GetPHOSCells());
            clus->SetE(energy);
          }
        }
        
        TLorentzVector momentum ;
        
        clus->GetMomentum(momentum, fVertex[vindex]);
        
        if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"PHOS")) continue;
        
        if(fPHOSPtMin > momentum.E() || fPHOSPtMax < momentum.E())          continue;
        
        if(fDebug > 2 && momentum.E() > 0.1)
          printf("AliCaloTrackReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
                 momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
        
        
        if (fMixedEvent)
        {
          clus->SetID(iclus) ;
        }
        
        if(fAcceptOnlyHIJINGLabels && !IsHIJINGLabel(clus->GetLabel())) continue ;
        
        fPHOSClusters->Add(clus);
        
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod entries %d\n",  fPHOSClusters->GetEntriesFast());
  
}

//____________________________________________
void AliCaloTrackReader::FillInputEMCALCells()
{
  //Return array with EMCAL cells in aod format
  
  fEMCALCells = fInputEvent->GetEMCALCells();
  
}

//___________________________________________
void AliCaloTrackReader::FillInputPHOSCells()
{
  //Return array with PHOS cells in aod format
  
  fPHOSCells = fInputEvent->GetPHOSCells();
  
}

//_______________________________________
void AliCaloTrackReader::FillInputVZERO()
{
  //Fill VZERO information in data member, add all the channels information.
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
    if(fDebug > 0)
      printf("AliCaloTrackReader::FillInputVZERO() - ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  }
  else
  {
    if(fDebug > 0)
      printf("AliCaloTrackReader::FillInputVZERO() - Cannot retrieve V0 ESD! Run w/ null V0 charges\n ");
  }
}

//_________________________________________________
void AliCaloTrackReader::FillInputNonStandardJets()
{
  //
  //fill array with non standard jets
  //
  // Adam T. Matyja
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputNonStandardJets()\n");
  //
  //check if branch name is given
  if(!fInputNonStandardJetBranchName.Length())
  {
    Printf("No non-standard jet branch name specified. Specify among existing ones.");
    fInputEvent->Print();
    abort();
  }
  
  fNonStandardJets = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fInputNonStandardJetBranchName.Data()));
  
  if(!fNonStandardJets)
  {
    //check if jet branch exist; exit if not
    Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fInputNonStandardJetBranchName.Data());
    fInputEvent->Print();
    abort();
  }
  else
  {
    if(fDebug > 1)
      printf("AliCaloTrackReader::FillInputNonStandardJets() - aod input jets %d\n", fNonStandardJets->GetEntriesFast() );
  }
  
}

//_________________________________________________
void AliCaloTrackReader::FillInputBackgroundJets()
{
  //
  //fill array with Background jets
  //
  // Adam T. Matyja
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputBackgroundJets()\n");
  //
  //check if branch name is given
  if(!fInputBackgroundJetBranchName.Length())
  {
    Printf("No background jet branch name specified. Specify among existing ones.");
    fInputEvent->Print();
    abort();
  }
  
  fBackgroundJets = (AliAODJetEventBackground*)(fInputEvent->FindListObject(fInputBackgroundJetBranchName.Data()));
  
  if(!fBackgroundJets)
  {
    //check if jet branch exist; exit if not
    Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fInputBackgroundJetBranchName.Data());
    fInputEvent->Print();
    abort();
  }
  else
  {
    if(fDebug > 1){
      printf("AliCaloTrackReader::FillInputBackgroundJets()\n");
      fBackgroundJets->Print("");
    }
  }
  
}


//________________________________________________
Bool_t AliCaloTrackReader::CheckForPrimaryVertex()
{
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //Only for ESDs ...
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (fInputEvent);
  if(!event) return kTRUE;
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0)
  {
    return kTRUE;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1)
  {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0)
    {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1)
    {
      //      cout<<"bad vertex type::"<< event->GetPrimaryVertex()->GetName() << endl;
      return kFALSE;
    }
  }
  
  return kFALSE;
  
}

//________________________________________________________________________________
TArrayI AliCaloTrackReader::GetTriggerPatches(Int_t tmin, Int_t tmax )
{
  // Select the patches that triggered
  // Depend on L0 or L1
	
  // some variables
  Int_t  trigtimes[30], globCol, globRow,ntimes, i;
  Int_t  absId  = -1; //[100];
  Int_t  nPatch = 0;
	
  TArrayI patches(0);
  
  // get object pointer
  AliVCaloTrigger *caloTrigger = GetInputEvent()->GetCaloTrigger( "EMCAL" );

  // Recover the threshold of the event that triggered, only possible for L1
  if(!fTriggerL1EventThresholdFix)
  {
    if(fBitEGA==6)
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
  
  if(patches.GetSize()<=0) printf("AliCaloTrackReader::GetTriggerPatches() - No patch found! for triggers: %s and selected <%s>\n",
                                  GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data());
  //else                     printf(">>>>> N patches %d, test %d,first %d, last %d\n",patches.GetSize(), nPatch, patches.At(0), patches.At(patches.GetSize()-1));
                 
  return patches;
}

//______________________________________________________________________
void  AliCaloTrackReader::MatchTriggerCluster(TArrayI patches)
{
  // Finds the cluster that triggered
  
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
  
  Float_t triggerThreshold = fTriggerL1EventThreshold;
  if(IsEventEMCALL0()) triggerThreshold = fTriggerL0EventThreshold;
  //printf("Threshold %f\n",triggerThreshold);
  Float_t minE = triggerThreshold / 2.;

  // This method is not really suitable for JET trigger
  // but in case, reduce the energy cut since we do not trigger on high energy particle
  if(IsEventEMCALL1Jet() || minE < 1) minE = 1;

  //printf("Min trigger Energy threshold %f\n",minE);
  
  // Loop on the clusters, check if there is any that falls into one of the patches
  for (Int_t iclus =  0; iclus <  nclusters; iclus++)
  {
    AliVCluster * clus = 0;
    if(clusterList) clus = (AliVCluster*) clusterList->At(iclus);
    else            clus = fInputEvent->GetCaloCluster(iclus);
    
    if ( !clus )                continue ;
    
    if ( !IsEMCALCluster(clus)) continue ;
    
    //Skip clusters with too low energy to be triggering
    if ( clus->E() < minE )    continue ;
    
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
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn() && fTriggerClusterTimeRecal)
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
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

//________________________________________________________
void AliCaloTrackReader::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
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
  printf("EMCAL Time Cut: %3.1f < TOF  < %3.1f\n", fEMCALTimeCutMin, fEMCALTimeCutMax);
  printf("Use CTS         =     %d\n",     fFillCTS) ;
  printf("Use EMCAL       =     %d\n",     fFillEMCAL) ;
  printf("Use PHOS        =     %d\n",     fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n",     fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n",     fFillPHOSCells) ;
  printf("Track status    =     %d\n", (Int_t) fTrackStatus) ;
  printf("AODs Track filter mask  =  %d or hybrid %d (if filter bit comp %d), select : SPD hit %d, primary %d\n",
         (Int_t) fTrackFilterMask, fSelectHybridTracks, (Int_t) fTrackFilterMaskComplementary, fSelectSPDHitTracks,fSelectPrimaryTracks) ;
  printf("Track Mult Eta Cut =  %d\n", (Int_t) fTrackMultEtaCut) ;
  printf("Write delta AOD =     %d\n",     fWriteOutputDeltaAOD) ;
  printf("Recalculate Clusters = %d, E linearity = %d\n",    fRecalculateClusters, fCorrectELinearity) ;
  
  printf("Use Triggers selected in SE base class %d; If not what Trigger Mask? %d; MB Trigger Mask for mixed %d \n",
         fEventTriggerAtSE, fEventTriggerMask,fMixEventTriggerMask);
  
  if(fComparePtHardAndClusterPt)
    printf("Compare jet pt and pt hard to accept event, factor = %2.2f",fPtHardAndJetPtFactor);
  
  if(fComparePtHardAndClusterPt)
    printf("Compare cluster pt and pt hard to accept event, factor = %2.2f",fPtHardAndClusterPtFactor);
  
  printf("Read Kine from, stack? %d, AOD ? %d \n", fReadStack, fReadAODMCParticles) ;
  printf("Delta AOD File Name =     %s\n", fDeltaAODFileName.Data()) ;
  printf("Centrality: Class %s, Option %d, Bin [%d,%d] \n", fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1]) ;
  
  printf("    \n") ;
  
}

//__________________________________________
Bool_t  AliCaloTrackReader::RejectLEDEvents()
{
  // LED Events in period LHC11a contaminated sample, simple method
  // to reject such events
  
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
    if(fDebug > 0)
      printf(" AliCaloTrackReader::FillInputEvent() - reject event with ncells in SM3 %d, cut %d, trig %s\n",
             ncellsSM3,ncellcut,GetFiredTriggerClasses().Data());
    return kTRUE;
  }
  
  return kFALSE;
  
}

//_________________________________________________________
void AliCaloTrackReader::RemapMCLabelForAODs(Int_t & label)
{
  // MC label for Cells not remapped after ESD filtering, do it here.
  
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
void AliCaloTrackReader::ResetLists()
{
  //  Reset lists, called by the analysis maker
  
  if(fCTSTracks)       fCTSTracks     -> Clear();
  if(fEMCALClusters)   fEMCALClusters -> Clear("C");
  if(fPHOSClusters)    fPHOSClusters  -> Clear("C");
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0;
  fV0Mul[0] = 0;   fV0Mul[1] = 0;
  
  if(fNonStandardJets) fNonStandardJets -> Clear("C");
  fBackgroundJets->Reset();

}

//___________________________________________
void AliCaloTrackReader::SetEventTriggerBit()
{
  // Tag event depeding on trigger name
	
  fEventTrigMinBias       = kFALSE;
  fEventTrigCentral       = kFALSE;
  fEventTrigSemiCentral   = kFALSE;
  fEventTrigEMCALL0       = kFALSE;
  fEventTrigEMCALL1Gamma1 = kFALSE;
  fEventTrigEMCALL1Gamma2 = kFALSE;
  fEventTrigEMCALL1Jet1   = kFALSE;
  fEventTrigEMCALL1Jet2   = kFALSE;
  
  if(fDebug > 0)
    printf("AliCaloTrackReader::SetEventTriggerBit() - Select trigger mask bit %d - Trigger Event %s - Select <%s>\n",
           fEventTriggerMask,GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data());
  
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
  
  if(fDebug > 0 )
    printf("AliCaloTrackReader::SetEventTriggerBit() - Event bits: \n \t MB   %d, Cen  %d, Sem  %d, L0   %d, L1G1 %d, L1G2 %d, L1J1 %d, L1J2 %d \n",
           fEventTrigMinBias,      fEventTrigCentral,       fEventTrigSemiCentral,
           fEventTrigEMCALL0 ,     fEventTrigEMCALL1Gamma1, fEventTrigEMCALL1Gamma2,
           fEventTrigEMCALL1Jet1 , fEventTrigEMCALL1Jet2);
  
  if(fBitEGA == 0 && fBitEJE ==0)
  {
    // Init the trigger bit once, correct depending on AliESDAODCaloTrigger header version
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
	    }  else printf("AliCaloTrackReader()::SetEventTriggerBit() - Streamer info for trigger class not available, bit not changed\n");
    } else printf("AliCaloTrackReader::SetEventTriggerBit() -  Streamer list not available!, bit not changed\n");
    
  } // set once the EJE, EGA trigger bit
  
}

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

//____________________________________________________________
void  AliCaloTrackReader::SetTrackCuts(AliESDtrackCuts * cuts)
{
  // Set Track cuts
  
  if(fESDtrackCuts) delete fESDtrackCuts ;
  
  fESDtrackCuts = cuts ;
  
}

//_________________________________________________________________________
void  AliCaloTrackReader::SetTrackComplementaryCuts(AliESDtrackCuts * cuts)
{
  // Set Track cuts for complementary tracks (hybrids)
  
  if(fESDtrackComplementaryCuts) delete fESDtrackComplementaryCuts ;
  
  fESDtrackComplementaryCuts = cuts ;
  
}


