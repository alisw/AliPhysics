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
#include <TObjString.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TStreamerInfo.h>

// ---- ANALYSIS system ----
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenEventHeader.h"
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
#include "AliMultSelection.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

// ---- Detectors ----
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEmcalTriggerDecisionContainer.h"

// ---- CaloTrackCorr ---
#include "AliCalorimeterUtils.h"
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"

// ---- Jets ----
#include "AliAODJet.h"
//#include "AliAODJetEventBackground.h"

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
fComparePtHardAndPromptPhotonPt(0),fPtHardAndPromptPhotonPtFactor(0),
fCTSPtMin(0),                fEMCALPtMin(0),                  fPHOSPtMin(0),
fCTSPtMax(0),                fEMCALPtMax(0),                  fPHOSPtMax(0),
fEMCALBadChMinDist(0),       fPHOSBadChMinDist (0),           
fEMCALNCellsCut(0),          fEMCALNCellsCutEnDepEnMin(0),
fEMCALNCellsCutEnDepConstant(0), fEMCALNCellsCutEnDepSlope(0),      
fPHOSNCellsCut(0),
fEMCALHighEnergyNdiffCut(0), fEMCALMinCellEnNdiffCut(0),
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
fInputEvent(0x0),            fOutputEvent(0x0),               fMC(0x0),
fSelectEmbeddedClusters(kFALSE),
fFillCTS(0),                 fFillEMCAL(0),
fFillDCAL(0),                fFillPHOS(0),
fFillEMCALCells(0),          fFillPHOSCells(0),
fRecalculateClusters(kFALSE),fCorrectELinearity(kTRUE),
fScaleEPerSM(kFALSE),       
fSmearShowerShape(0),        fSmearShowerShapeWidth(0),       fRandom(),
fSmearingFunction(0),        fSmearNLMMin(0),                 fSmearNLMMax(0),
fTrackStatus(0),             fSelectSPDHitTracks(0),          
fSelectMinITSclusters(0),    fSelectMaxChi2PerITScluster(10000),
fSelectMinTPCclusters(0),    fSelectMaxChi2PerTPCcluster(10000),
fTrackMultNPtCut(0),         fTrackMultEtaCut(0.9),
fDeltaAODFileName(""),       fFiredTriggerClassName(""),

fEventTriggerMaskInput(0),   fEventTriggerMask(0),        fMixEventTriggerMask(0),         
fEventTriggerAtSE(0),
fEventTrigMinBias(0),        fEventTrigCentral(0),
fEventTrigSemiCentral(0),    fEventTrigEMCALL0(0),
fEventTrigEMCALL1Gamma1(0),  fEventTrigEMCALL1Gamma2(0),
fEventTrigEMCALL1Jet1(0),    fEventTrigEMCALL1Jet2(0),
fEventTrigDCALL0(0),
fEventTrigDCALL1Gamma1(0),   fEventTrigDCALL1Gamma2(0),
fEventTrigDCALL1Jet1(0),     fEventTrigDCALL1Jet2(0),

fEventTrigMinBiasCaloOnly(0),        fEventTrigEMCALL0CaloOnly(0),
fEventTrigEMCALL1Gamma1CaloOnly(0),  fEventTrigEMCALL1Gamma2CaloOnly(0),
fEventTrigEMCALL1Jet1CaloOnly(0),    fEventTrigEMCALL1Jet2CaloOnly(0),
fEventTrigDCALL0CaloOnly(0),
fEventTrigDCALL1Gamma1CaloOnly(0),   fEventTrigDCALL1Gamma2CaloOnly(0),
fEventTrigDCALL1Jet1CaloOnly(0),     fEventTrigDCALL1Jet2CaloOnly(0),

fBitEGA(0),                  fBitEJE(0),

fEventType(-1),
fTaskName(""),               fCaloUtils(0x0),                 fMCUtils(0x0), 
fWeightUtils(0x0),           fEventWeight(1),
fMixedEvent(NULL),           fNMixedEvent(0),                 fVertex(NULL),
fEventCuts(1),               fUseEventCutsClass(kFALSE),
fListMixedTracksEvents(),    fListMixedCaloEvents(),
fLastMixedTracksEvent(-1),   fLastMixedCaloEvent(-1),
fWriteOutputDeltaAOD(kFALSE),
fEMCALClustersListName(""),  fEMCALCellsListName(""),  
fZvtxCut(0.),
fAcceptFastCluster(kFALSE),  fRemoveLEDEvents(0),
fLEDHighEnergyCutSM(0),      fLEDHighNCellsCutSM(0), 
fLEDLowEnergyCutSM3(0),      fLEDLowNCellsCutSM3(0),
fLEDMinCellEnergy(0),        fLEDMaxCellEnergy(0),
fRemoveLEDStripEvents(0),    fLEDEventMaxNumberOfStrips(0),
fLEDLowEnergyCutSM3Strip(0), fLEDLowNCellsCutSM3Strip(0),

//Trigger rejection
fRemoveBadTriggerEventsFromEMCalTriggerMaker(0),
fEMCalTriggerMakerDecissionContainerName(0),
fRemoveBadTriggerEvents(0),  fTriggerPatchClusterMatch(0),
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
fTimeStampEventCTPBCCorrExclude(0),
fTimeStampEventCTPBCCorrMin(0), fTimeStampEventCTPBCCorrMax(0),
fNPileUpClusters(-1),        fNNonPileUpClusters(-1),         fNPileUpClustersCut(3),
fVertexBC(-200),             fRecalculateVertexBC(0),
fUseAliCentrality(0),        fMultWithEventSel(0),
fCentralityClass(""),        fCentralityOpt(0),
fEventPlaneMethod(""),
fFillInputNonStandardJetBranch(kFALSE),
fNonStandardJets(new TClonesArray("AliAODJet",100)),          fInputNonStandardJetBranchName("jets"),
fFillInputBackgroundJetBranch(kFALSE), 
//fBackgroundJets(0x0),
fBackgroundJets(new TClonesArray("AliAODJet",100)),
fInputBackgroundJetBranchName("jets"),
fAcceptEventsWithBit(0),     fRejectEventsWithBit(0),         
fRejectEMCalTriggerEventsL1HighWithL1Low(0),
fRemoveCentralityTriggerOutliers(0),
fMomentum(),                 fParRun(kFALSE),                 fCurrentParIndex(0),
fOutputContainer(0x0),       fhEMCALClusterEtaPhi(0),         fhEMCALClusterEtaPhiFidCut(0),     
fhEMCALClusterDisToBadE(0),  fhEMCALClusterTimeE(0),      
fhEMCALClusterBadTrigger(0), fhCentralityBadTrigger(0),       fhEMCALClusterCentralityBadTrigger(0),
fhEMCALNSumEnCellsPerSM(0),    fhEMCALNSumEnCellsPerSMAfter(0), fhEMCALNSumEnCellsPerSMAfterStripCut(0),
fhEMCALNSumEnCellsPerStrip(0), fhEMCALNSumEnCellsPerStripAfter(0),
fhPtHardPtJetPtRatio(0),     fhPtHardPromptPhotonPtRatio(0),
fhPtHardEnClusterRatio(0),   fhPtHardEnClusterCenRatio(0),
fEnergyHistogramNbins(0),    fHistoCentDependent(0),          fHistoPtDependent(0),
fhNEventsAfterCut(0),        fNMCGenerToAccept(0),            fMCGenerEventHeaderToAccept(""),
fGenEventHeader(0),          fGenPythiaEventHeader(0),        fCheckPythiaEventHeader(1),
fAcceptMCPromptPhotonOnly(0),fRejectMCFragmentationPhoton(0)
{
  for(Int_t i = 0; i < 9; i++) fhEMCALClusterCutsE   [i]= 0x0 ;
  for(Int_t i = 0; i < 9; i++) fhEMCALClusterCutsECen[i]= 0x0 ;
  for(Int_t i = 0; i < 9; i++) fhEMCALClusterCutsESignal   [i]= 0x0 ;
  for(Int_t i = 0; i < 9; i++) fhEMCALClusterCutsECenSignal[i]= 0x0 ;
  for(Int_t i = 0; i < 7; i++) fhPHOSClusterCutsE    [i]= 0x0 ;
  for(Int_t i = 0; i < 6; i++) fhCTSTrackCutsPt      [i]= 0x0 ;
  for(Int_t i = 0; i < 6; i++) fhCTSTrackCutsPtSignal[i]= 0x0 ;
  for(Int_t i = 0; i < 6; i++) fhCTSTrackCutsPtCen   [i]= 0x0 ;
  for(Int_t i = 0; i < 6; i++) fhCTSTrackCutsPtCenSignal[i]= 0x0 ;
  for(Int_t j = 0; j < 5; j++) { fMCGenerToAccept    [j] =  ""; fMCGenerIndexToAccept[j] = -1; }
  
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

  if(fBackgroundJets)
  {
    if(fDataType!=kMC) fBackgroundJets->Clear("C") ;
    else               fBackgroundJets->Delete() ;
    delete fBackgroundJets;
  }
  //  delete fBackgroundJets ;

  fRejectEventsWithBit.Reset();
  fAcceptEventsWithBit.Reset();
  
  if ( fWeightUtils ) delete fWeightUtils ;
    
  if ( fMCUtils     ) delete fMCUtils ; 

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
Bool_t  AliCaloTrackReader::AcceptEventWithTriggerBit(UInt_t trigFired)
{  
  Int_t nAccept = fAcceptEventsWithBit.GetSize();
  
  //printf("N accept %d\n", nAccept);
  
  if( nAccept <= 0 )
    return kTRUE ; // accept the event
  
  //UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
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
  
  if ( !mcpart0 )
  {
    AliWarning(Form("AliMCEvent-BREAK: No valid AliMCParticle at label %i",index));
    return -1;
  }
  
  nameGen = GetGeneratorNameAndIndex(index,genIndex);
  
  if(nameGen.Contains("nococktailheader") ) return -1;
  
  Int_t lab=index;
  
  while(nameGen.IsWhitespace())
  {
    AliVParticle* mcpart = (AliVParticle*) GetMC()->GetTrack(lab);
    
    if ( !mcpart )
    {
      AliWarning(Form("AliMCEvent-BREAK: No valid AliMCParticle at label %i",lab));
      break;
    }
    
    Int_t mother=0;
    mother = mcpart->GetMother();
    
    if ( mother<0 )
    {
      AliWarning("AliMCEvent - BREAK: Reached primary particle without valid mother");
      break;
    }
    
    AliVParticle* mcmom = (AliVParticle*) GetMC()->GetTrack(mother);
    if ( !mcmom )
    {
      AliWarning(Form("AliMCEvent-BREAK: No valid AliMCParticle mother at label %i",mother));
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
Bool_t  AliCaloTrackReader::RejectEventWithTriggerBit(UInt_t trigFired)
{
  Int_t nReject = fRejectEventsWithBit.GetSize();
  
  //printf("N reject %d\n", nReject);
  
  if( nReject <= 0 )
    return kTRUE ; // accept the event
  
  //UInt_t trigFired = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
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

  fhNEventsAfterCut->Fill(1.5);
  
  //-----------------------------------------------------------------
  // In case of mixing analysis, select here the trigger of the event
  //-----------------------------------------------------------------
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if ( !inputHandler ) return kFALSE ;  // to content coverity

  Bool_t isTrigger = kFALSE;
  Bool_t isMB      = kFALSE;
  fEventTriggerMaskInput  = inputHandler->IsEventSelected();
//fEventTriggerMaskInput = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  if ( !fEventTriggerAtSE )
  {
    // In case of mixing analysis, accept MB events, not only Trigger
    // Track and cluster arrays filled for MB in order to create the pool in the corresponding analysis
    // via de method in the base class FillMixedEventPool()
    
    isTrigger = fEventTriggerMaskInput & fEventTriggerMask;
    isMB      = fEventTriggerMaskInput & fMixEventTriggerMask;
    
    if ( !isTrigger && !isMB ) return kFALSE;
    
    //printf("Selected triggered event : %s\n",GetFiredTriggerClasses().Data());
    AliDebug(1,"Pass uninteresting triggered events rejection in case of mixing analysis");  
    
    fhNEventsAfterCut->Fill(2.5);
  }

  //-----------------------------------------------------------
  // Reject events depending on the trigger name 
  // Careful!, if a special MB event string is selected but the option
  // to select events via the mask in the reader is done, it will not 
  // be taken into account.
  //-----------------------------------------------------------
  
  AliDebug(1,Form("FiredTriggerClass <%s>, selected class <%s>, compare name %d",
                  GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data(),
                  GetFiredTriggerClasses().Contains(fFiredTriggerClassName)));  
  
  if ( fFiredTriggerClassName != "" && !isMB )
  {
    if ( !GetFiredTriggerClasses().Contains(fFiredTriggerClassName) ) 
      return kFALSE;
    
    AliDebug(1,"Accepted triggered event");
    
    fhNEventsAfterCut->Fill(3.5);
  }
  
  //-------------------------------------------------------------------------------------
  // Reject or accept events depending on the trigger bit
  //-------------------------------------------------------------------------------------
  
  Bool_t okA = AcceptEventWithTriggerBit(fEventTriggerMaskInput);
  Bool_t okR = RejectEventWithTriggerBit(fEventTriggerMaskInput);
  
  //printf("AliCaloTrackReader::FillInputEvent() - Accept event? %d, Reject event %d? \n",okA,okR);
  
  if ( !okA || !okR ) return kFALSE;
  
  AliDebug(1,"Pass event bit rejection");
  
  fhNEventsAfterCut->Fill(4.5);
  
  //----------------------------------------------------------------------
  // Do not count events that were likely triggered by an exotic cluster
  // or out BC cluster
  //----------------------------------------------------------------------
  
  // Set a bit with the event kind, MB, L0, L1 ...
  SetEventTriggerBit(fEventTriggerMaskInput);
  
  // In case of Mixing, avoid checking the triggers in the min bias events
  if ( !fEventTriggerAtSE && (isMB && !isTrigger) ) return kTRUE;
  
  // Reject triggered events when there is coincidence on both EMCal/DCal L1 high and low trigger thresholds,
  // but the requested trigger is the high trigger threshold
  // Check trigger string selection set in ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
  //
  if ( fRejectEMCalTriggerEventsL1HighWithL1Low )
  {    
    if ( IsEventEMCALL1Jet1  () && IsEventEMCALL1Jet2  () && fFiredTriggerClassName.Contains("J1") ) return kFALSE;
    if ( IsEventEMCALL1Gamma1() && IsEventEMCALL1Gamma2() && fFiredTriggerClassName.Contains("G1") ) return kFALSE;
    if ( IsEventDCALL1Jet1   () && IsEventDCALL1Jet2   () && fFiredTriggerClassName.Contains("J1") ) return kFALSE;
    if ( IsEventDCALL1Gamma1 () && IsEventDCALL1Gamma2 () && fFiredTriggerClassName.Contains("G1") ) return kFALSE;
    
    // Not sure if coincidences with kCaloOnly are possible but just in case
    if ( fEventTriggerMaskInput & AliVEvent::kCaloOnly )
    {
      if ( IsEventEMCALL1Jet1CaloOnly  () && IsEventEMCALL1Jet2CaloOnly  () && fFiredTriggerClassName.Contains("J1") ) return kFALSE;
      if ( IsEventEMCALL1Gamma1CaloOnly() && IsEventEMCALL1Gamma2CaloOnly() && fFiredTriggerClassName.Contains("G1") ) return kFALSE;
      if ( IsEventDCALL1Jet1CaloOnly   () && IsEventDCALL1Jet2CaloOnly   () && fFiredTriggerClassName.Contains("J1") ) return kFALSE;
      if ( IsEventDCALL1Gamma1CaloOnly () && IsEventDCALL1Gamma2CaloOnly () && fFiredTriggerClassName.Contains("G1") ) return kFALSE;
      
      // Coincidence L0-L2
      if ( IsEventDCALL0CaloOnly() && IsEventDCALL1Gamma2CaloOnly() && fFiredTriggerClassName.Contains("G2") ) return kFALSE;
      if ( IsEventDCALL0CaloOnly() && IsEventDCALL1Gamma1CaloOnly() && fFiredTriggerClassName.Contains("G1") ) return kFALSE;
    }
    
     fhNEventsAfterCut->Fill(5.5);
  }
  
  // Reject events from centrality triggers with centrality out of expected range
  //
  if ( fRemoveCentralityTriggerOutliers  )
  {
    Float_t centrality = GetEventCentralityF();
//    printf("Check outliers for cent %2.1f, central? %d, semicentral? %d; mb %d; run %d\n",
//           centrality, fEventTrigCentral, fEventTrigSemiCentral,fEventTrigMinBias, fInputEvent->GetRunNumber());

    // In case of OR of all MB triggers, do not discard events considered as pure MB
    Bool_t checkMBcent = kTRUE;
    if ( ( (fEventTriggerMask & AliVEvent::kCentral) || (fEventTriggerMask & AliVEvent::kSemiCentral) ) && 
         ( (fEventTriggerMask & AliVEvent::kMB)      || (fEventTriggerMask & AliVEvent::kINT7)        )    )
    {
      if ( fEventTrigMinBias ) checkMBcent = kFALSE;
    }
    
    if ( checkMBcent )
    {
      if ( fEventTrigSemiCentral && (fEventTriggerMask & AliVEvent::kSemiCentral) ) 
      {
        Int_t centMin = 0; // LHC11h
        Int_t centMax = 50;
        if ( fInputEvent->GetRunNumber() > 295274 ) 
        {
          centMin = 30; // LHC18qr
        }
        
        if ( centrality < centMin ) 
        {
          // Do not skip good central events when central mask
          if (  ( (fEventTriggerMask & AliVEvent::kCentral) && fEventTrigCentral && centrality >= 10) || !fEventTrigCentral )
          {
            //printf("%s\n",GetFiredTriggerClasses().Data());
            AliInfo(Form("Skip semi-central event with centrality %2.1f, out of [%d,%d]",
                         centrality, centMin, centMax));
            return kFALSE;
          }
        }
        else if  ( centrality >= centMax  ) 
        {
          AliInfo(Form("Skip semi-central event with centrality %2.1f, out of [%d,%d]",
                       centrality, centMin, centMax));
          return kFALSE;
        }
        
      }
      
      if ( fEventTrigCentral && centrality >= 10  && (fEventTriggerMask & AliVEvent::kCentral) ) 
      {
        printf("%s\n",GetFiredTriggerClasses().Data());
        AliInfo(Form("Skip central event with centrality %2.1f",centrality));
        return kFALSE;
      }
      
    }
    
    if ( (fEventTriggerMask & AliVEvent::kEMCEGA) || (fEventTriggerMask & AliVEvent::kCaloOnly) )
    {
      if  ( fEventTrigEMCALL1Gamma2 || fEventTrigEMCALL1Gamma2CaloOnly || 
            fEventTrigDCALL1Gamma2  || fEventTrigDCALL1Gamma2CaloOnly    )
      {
        if ( centrality < 50 && fInputEvent->GetRunNumber() > 295274 && fFiredTriggerClassName.Contains("G2"))
        {
          //printf("%s\n",GetFiredTriggerClasses().Data());
          AliInfo(Form("Skip L1-G2 event with centrality %2.1f",centrality));
          return kFALSE;
        }
      } // L1-Low threshold
      
      if ( (fEventTrigEMCALL1Gamma1 || fEventTrigEMCALL1Gamma1CaloOnly || 
            fEventTrigDCALL1Gamma1  || fEventTrigDCALL1Gamma1CaloOnly)   )
      {
        if ( centrality > 50 && fInputEvent->GetRunNumber() > 295274 && fFiredTriggerClassName.Contains("G1") )
        {
          //printf("%s\n",GetFiredTriggerClasses().Data());
          AliInfo(Form("Skip L1-G1 event with centrality %2.1f",centrality));
          return kFALSE;
        }
      } // L1-High threshold
    } // EMCal triggers
    
  } //  fRemoveCentralityTriggerOutliers
  
  // Match triggers
  //
  if ( fTriggerPatchClusterMatch &&
      ( IsEventEMCALL1() || IsEventEMCALL0() || IsEventDCALL1() || IsEventDCALL0() ) )
  {
    //Get Patches that triggered
    TArrayI patches = GetTriggerPatches(fTriggerPatchTimeWindow[0],fTriggerPatchTimeWindow[1]);
    
    MatchTriggerCluster(patches);
    
    patches.Reset();
    
    // If requested, remove badly triggered events, but only when the EMCal trigger bit is set
    if ( fRemoveBadTriggerEvents )
    {
      AliDebug(1,Form("ACCEPT triggered event? \n exotic? %d - bad cell %d - bad Max cell %d - BC %d  - Matched %d\n",
                      fIsExoticEvent,fIsBadCellEvent, fIsBadMaxCellEvent, fTriggerClusterBC,fIsTriggerMatch));
      
      if     (fIsExoticEvent)         return kFALSE;
      else if(fIsBadCellEvent)        return kFALSE;
      else if(fRemoveUnMatchedTriggers && !fIsTriggerMatch) return kFALSE ;
      else if(fTriggerClusterBC != 0) return kFALSE;
      
      AliDebug(1,Form("\t *** YES for %s",GetFiredTriggerClasses().Data()));
    }
    
    AliDebug(1,"Pass EMCal triggered event rejection"); 
    
    fhNEventsAfterCut->Fill(6.5);
  }
  else if ( fRemoveBadTriggerEventsFromEMCalTriggerMaker )
  {
    auto trgsel = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(GetInputEvent()->FindListObject(fEMCalTriggerMakerDecissionContainerName));
    if ( trgsel )
    {
      AliDebug(1,Form("Trigger Maker, check decision: EG1 %d, EG2 %d, DG1 %d, DG2 %d, EGA %d; EMCL0 %d, DMCL0 %d request %s",
                      trgsel->IsEventSelected("EG1"),trgsel->IsEventSelected("EG2"),
                      trgsel->IsEventSelected("DG1"),trgsel->IsEventSelected("DG2"),
                      trgsel->IsEventSelected("EGA"),
                      trgsel->IsEventSelected("EMCL0"),trgsel->IsEventSelected("DMCL0"),
                      fFiredTriggerClassName.Data()));
      
      Bool_t reject = kFALSE;
      // Check trigger string selection set in ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
      if      ( fFiredTriggerClassName.Contains("EG1") && !trgsel->IsEventSelected("EG1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("EGA") && !trgsel->IsEventSelected("EGA") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("DG1") && !trgsel->IsEventSelected("DG1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("EG2") && !trgsel->IsEventSelected("EG2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("DG2") && !trgsel->IsEventSelected("DG2") ) reject = kTRUE;
      
      else if ( fFiredTriggerClassName.Contains("EMC") && !trgsel->IsEventSelected("EMCL0") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("DMC") && !trgsel->IsEventSelected("DMCL0") ) reject = kTRUE;

      else if ( fFiredTriggerClassName.Contains("EJ1") && !trgsel->IsEventSelected("EJ1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("DJ1") && !trgsel->IsEventSelected("DJ1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("EJ2") && !trgsel->IsEventSelected("EJ2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName.Contains("DJ2") && !trgsel->IsEventSelected("DJ2") ) reject = kTRUE;
      
      else if ( fFiredTriggerClassName == "EG" && !trgsel->IsEventSelected("EG1") && !trgsel->IsEventSelected("EG2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "DG" && !trgsel->IsEventSelected("DG1") && !trgsel->IsEventSelected("DG2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "EJ" && !trgsel->IsEventSelected("EJ1") && !trgsel->IsEventSelected("EJ2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "DJ" && !trgsel->IsEventSelected("DJ1") && !trgsel->IsEventSelected("DJ2") ) reject = kTRUE;

      else if ( fFiredTriggerClassName == "G1" && !trgsel->IsEventSelected("EG1") && !trgsel->IsEventSelected("DG1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "G2" && !trgsel->IsEventSelected("EG2") && !trgsel->IsEventSelected("DG2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "J1" && !trgsel->IsEventSelected("EJ1") && !trgsel->IsEventSelected("DJ1") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "J2" && !trgsel->IsEventSelected("EJ2") && !trgsel->IsEventSelected("DJ2") ) reject = kTRUE;
      
      else if ( fFiredTriggerClassName == "MC" && !trgsel->IsEventSelected("EMCL0") && !trgsel->IsEventSelected("DMCL0") ) reject = kTRUE;
      
      else if ( fFiredTriggerClassName == "G"  && !trgsel->IsEventSelected("EG1") && !trgsel->IsEventSelected("DG1") && 
                                                  !trgsel->IsEventSelected("EG2") && !trgsel->IsEventSelected("DG2") ) reject = kTRUE;
      else if ( fFiredTriggerClassName == "J"  && !trgsel->IsEventSelected("EJ1") && !trgsel->IsEventSelected("DJ1") && 
                                                  !trgsel->IsEventSelected("EJ2") && !trgsel->IsEventSelected("DJ2") ) reject = kTRUE;
      
      if ( reject ) 
      {
        if ( fFillEMCAL )
        {
          fhCentralityBadTrigger->Fill(GetEventCentrality());       
          
          TClonesArray * clusterList = 0x0;
          if ( fEMCALClustersListName == "" )
            clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject("caloClusters"));
          if      (fInputEvent->FindListObject(fEMCALClustersListName))
            clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(fEMCALClustersListName));
          else if ( fOutputEvent )
            clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
          
          if ( clusterList )
          {
            Int_t nclusters = clusterList->GetEntriesFast();
            for (Int_t iclus =  0; iclus <  nclusters; iclus++)
            {
              AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
              
              if ( !clus )            continue;
              if ( !clus->IsEMCAL() ) continue;
              
              //printf("E %f\n",clus->E());
              
              fhEMCALClusterBadTrigger->Fill(clus->E()); 
              fhEMCALClusterCentralityBadTrigger->Fill(clus->E(), GetEventCentrality());
            } // cluster loop
          } // clusterList
          else AliError("No cluster list");
        }
      
        AliInfo(Form("Trigger Maker, event rejected! EG1 %d, EG2 %d, DG1 %d, DG2 %d, EGA %d, EMCL0 %d, DMCL0 %d; request %s",
                     trgsel->IsEventSelected("EG1")  , trgsel->IsEventSelected("EG2"),
                     trgsel->IsEventSelected("DG1")  , trgsel->IsEventSelected("DG2"),
                     trgsel->IsEventSelected("EGA")  ,
                     trgsel->IsEventSelected("EMCL0"), trgsel->IsEventSelected("DMCL0"),
                     fFiredTriggerClassName.Data()));
        
        return kFALSE;
      }
      
      AliDebug(1,"Pass EMCal triggered event rejection"); 
      
      fhNEventsAfterCut->Fill(6.5);
    }
    //else AliError("Trigger decision container not found, select event");
    
  } // fRemoveBadTriggerEventsFromEMCalTriggerMaker 
  
  //-------------------------------------------------------------------------------------
  // Select events only fired by a certain trigger configuration if it is provided
  //-------------------------------------------------------------------------------------
  
  if (GetFiredTriggerClasses().Contains("FAST")  && !GetFiredTriggerClasses().Contains("ALL") && !fAcceptFastCluster)
  {
    AliDebug(1,Form("Do not count events from fast cluster, trigger name %s\n",fFiredTriggerClassName.Data()));
    return kFALSE;
    
    fhNEventsAfterCut->Fill(7.5);
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
    
    fhNEventsAfterCut->Fill(8.5);
  } // Remove LED events

  // All selection criteria passed, accept the event
  return kTRUE ;
}

//________________________________________________
/// Check the MC PYTHIA event, if the requested 
/// pT-hard is much smaller than the jet pT, then,
/// there can be a problem in the tails of the 
/// distributions and the event should be rejected.
/// Do this only for pythia gamma-jet or jet-jet events
///
/// \param process pythia process from AliMCAnalysisUtils::GetPythiaEventHeader()
/// \param processName Jet-Jet or Gamma-Jet processes from AliMCAnalysisUtils::GetPythiaEventHeader()
//________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndJetPt(Int_t process, TString processName)
{  
  if ( !fGenEventHeader ) 
  {
    AliError("Skip event, event header is not available!");
    return kFALSE;
  }
  
  if ( fGenPythiaEventHeader )
  {  
    // Do this check only for jet-jet and gamma-jet productions
    
    if ( !processName.Contains("Jet") ) return kTRUE;
    
    Int_t nTriggerJets =  fGenPythiaEventHeader->NTriggerJets();
    Float_t ptHard = fGenPythiaEventHeader->GetPtHard();
    TParticle * jet =  0;

    AliDebug(1,Form("Njets: %d, pT Hard %f",nTriggerJets, ptHard));
    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      fGenPythiaEventHeader->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      
      AliDebug(1,Form("jet %d; pycell jet pT %f",ijet, jet->Pt()));

      if ( ptHard > 0 )
        fhPtHardPtJetPtRatio->Fill(jet->Pt()/ptHard);

      // Compare jet pT and pt Hard
      if ( jet->Pt() > fPtHardAndJetPtFactor * ptHard )
      {
        AliInfo(Form("Reject jet event with : process %d <%s>, pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f",
                     process, processName.Data(), ptHard, jet->Pt(), fPtHardAndJetPtFactor));

        return kFALSE;
      }
    } // jet loop
    
    if(jet) delete jet;
  } // pythia header
  
  return kTRUE ;
}

//____________________________________________________
/// Check the MC PYTHIA event, if the requested 
/// pT-hard is smaller than the calorimeter cluster E, 
/// there can be a problem in the tails of the 
/// distributions and the event should be rejected.
/// Do this only for pythia gamma-jet events
///
/// \param process pythia process from AliMCAnalysisUtils::GetPythiaEventHeader()
/// \param processName Jet-Jet or Gamma-Jet processes from AliMCAnalysisUtils::GetPythiaEventHeader()
//____________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndClusterPt(Int_t process, TString processName)
{ 
  if ( !fGenEventHeader ) 
  {
    AliError("Skip event, event header is not available!");
    return kFALSE;
  }

  if ( fGenPythiaEventHeader )
  {  
    // Do this check only for gamma-jet productions
    if ( processName != "Gamma-Jet" ) return kTRUE;
    
    Float_t ptHard = fGenPythiaEventHeader->GetPtHard();
    
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    if ( fEMCALClustersListName == "" )
    {
      for (Int_t iclus =  0; iclus <  nclusters; iclus++)
      {
        AliVCluster * clus = fInputEvent->GetCaloCluster(iclus) ;
        Float_t ecluster = clus->E();

        if ( ptHard > 0 )
        {
          fhPtHardEnClusterRatio->Fill(ecluster/ptHard);
          if ( fHistoCentDependent )
            fhPtHardEnClusterCenRatio->Fill(ecluster/ptHard, GetEventCentrality());
        }

        if ( ecluster > fPtHardAndClusterPtFactor * ptHard )
        {
          AliInfo(Form("Reject : process %d <%s>, ecluster %2.2f, calo %d, factor %2.2f, ptHard %f",
                       process, processName.Data(), ecluster ,clus->GetType(), fPtHardAndClusterPtFactor,ptHard));
          
          return kFALSE;
        }
      } // cluster loop
    }
    else
    {
      TClonesArray * clusterList = 0x0;
      
      if      ( fInputEvent->FindListObject(fEMCALClustersListName) )
      {
        clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(fEMCALClustersListName));
      }
      else if ( fOutputEvent )
      {
        clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
      }
      
      if ( !clusterList )
      {
        AliWarning(Form("Wrong name of list with clusters?  <%s>",fEMCALClustersListName.Data()));
        return kTRUE;
      }
      
      Int_t nclusters = clusterList->GetEntriesFast();
      for (Int_t iclus =  0; iclus <  nclusters; iclus++)
      {
        AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
        
        Float_t ecluster = clus->E();

        if ( ptHard > 0 )
        {
          fhPtHardEnClusterRatio->Fill(ecluster/ptHard);
          if ( fHistoCentDependent )
            fhPtHardEnClusterCenRatio->Fill(ecluster/ptHard, GetEventCentrality());
        }

        if ( ecluster > fPtHardAndClusterPtFactor * ptHard )
        {
          AliInfo(Form("Reject : process %d <%s>, ecluster %2.2f, calo %d, factor %2.2f, ptHard %f",
                       process, processName.Data(), ecluster ,clus->GetType(), fPtHardAndClusterPtFactor,ptHard));

          return kFALSE;
        }
      } // cluster loop
    }
  } // pythia header
  
  return kTRUE ;
}

//____________________________________________________
/// Check the MC PYTHIA event, if the requested
/// pT-hard is smaller than the generated prompt photon
/// there can be a problem in the tails of the
/// distributions and the event should be rejected.
/// Do this only for pythia gamma-jet events
///
/// \param process pythia process from AliMCAnalysisUtils::GetPythiaEventHeader()
/// \param processName Jet-Jet or Gamma-Jet processes from AliMCAnalysisUtils::GetPythiaEventHeader()
//____________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndPromptPhotonPt(Int_t process, TString processName)
{
  if ( !fGenEventHeader )
  {
    AliError("Skip event, event header is not available!");
    return kFALSE;
  }

  if ( !fGenPythiaEventHeader ) return kTRUE;

  // Do this check only for gamma-jet productions
  if ( processName != "Gamma-Jet" ) return kTRUE;

  Float_t ptHard = fGenPythiaEventHeader->GetPtHard();

  // Loop on pythia generated particles
  Int_t   firstParticle = 0 ;

  // Loop only over likely final particles not partons
  if ( GetGenPythiaEventHeader() )
    firstParticle = GetMCAnalysisUtils()->GetPythiaMaxPartParent();

  AliVParticle * primary = 0;

  Int_t    nprim     = GetMC()->GetNumberOfTracks();
  for(Int_t i = firstParticle ; i < nprim; i++)
  {
    if ( !AcceptParticleMCLabel( i ) ) continue ;

    primary = GetMC()->GetTrack(i) ;
    if ( !primary )
    {
      AliWarning("primaries pointer not available!!");
      continue;
    }

    // Select prompt photon
    if ( primary->PdgCode()      != 22 ) continue;
    if ( primary->MCStatusCode() != 1  ) continue;

    // Get tag of this particle photon from fragmentation, decay, prompt ...
    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(i, GetMC(),
                                                  GetNameOfMCEventHederGeneratorToAccept(),
                                                  primary->E()); // Not used, should be cluster
    if ( !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt) )
    {
      continue;
    }

    if ( ptHard > 0 )
      fhPtHardPromptPhotonPtRatio->Fill(primary->Pt()/ptHard);

    if ( primary->Pt() > fPtHardAndPromptPhotonPtFactor * ptHard )
    {
      AliInfo(Form("Reject : process %d <%s>, prompt photon %2.2f, factor %2.2f, ptHard %f",
                   process, processName.Data(), primary->Pt(), fPtHardAndPromptPhotonPtFactor,ptHard));

      return kFALSE;
    }
  } // cluster loop

  return kTRUE ;
}


//___________________________________________________
/// Fill the output list of initialized control histograms.
/// Cluster or track spectra histograms, depending on different selection cuts.
//___________________________________________________
TList * AliCaloTrackReader::GetCreateControlHistograms()
{  
  fhNEventsAfterCut = new TH1I("hNEventsAfterCut", "Number of analyzed events", 22, 0, 22) ;
  //fhNEventsAfterCut->SetXTitle("Selection");
  fhNEventsAfterCut->SetYTitle("# events");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(1 ,"1=Input");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(2 ,"2=Event Type");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(3 ,"3=Mixing Event");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(4 ,"4=Trigger string");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(5 ,"5=Trigger Bit");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(6 ,"6=L1 no L2"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(7 ,"7=Good EMC Trigger");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(8 ,"8=!Fast Cluster");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(9 ,"9=!LED");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(10,"10=Time stamp"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(11,"11=Primary vertex"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(12,"12=Null 3 vertex"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(13,"13=Z vertex window"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(14,"14=Pile-up"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(15,"15=V0AND"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(16,"16=Centrality"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(17,"17=GenHeader"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(18,"18=PtHard-Jet");
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(19,"19=PtHard-Cluster"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(20,"20=N Track>0"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(21,"21=TOF BC"); 
  fhNEventsAfterCut->GetXaxis()->SetBinLabel(22,"22=AliEventCuts"); 
  fOutputContainer->Add(fhNEventsAfterCut);

  if ( fFillEMCAL )
  {
    TString names[] =
    { "NoCut", "Corrected", "GoodCluster", "NonLinearity",
      "EnergyAndFidutial", "NCells", "BadDist", "Time","NcellsDiff" } ;
    
    for(Int_t i = 0; i < 9; i++)
    {  
      if ( names[i] == "Corrected"    && !fSelectEmbeddedClusters  && !fRecalculateClusters ) continue;
      if ( names[i] == "NonLinearity" && !fScaleEPerSM && !fCorrectELinearity ) continue;
      if ( names[i] == "BadDist"      && fEMCALBadChMinDist <= 0  )             continue;
      if ( names[i] == "NCells"       && fEMCALNCellsCut    <= 0  )             continue;
      if ( names[i] == "Time"         && !fUseEMCALTimeCut        )             continue;
      if ( names[i] == "NcellsDiff"   && (fEMCALHighEnergyNdiffCut > 200 || fEMCALHighEnergyNdiffCut < 40) ) continue;
      
      if ( !fHistoCentDependent )
      {
        fhEMCALClusterCutsE[i] = new TH1F
        (Form("hEMCALReaderClusterCuts_%d_%s",i,names[i].Data()),
         Form("EMCal %d, %s",i,names[i].Data()),
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]);
        fhEMCALClusterCutsE[i]->SetYTitle("# clusters");
        fhEMCALClusterCutsE[i]->SetXTitle("#it{E} (GeV)");
        if ( fHistoPtDependent )
          fhEMCALClusterCutsE[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fOutputContainer->Add(fhEMCALClusterCutsE[i]);
        
        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )//&&
           // !fSelectEmbeddedClusters && !fAcceptMCPromptPhotonOnly )
        {
          fhEMCALClusterCutsESignal[i] = new TH1F
          (Form("hEMCALReaderClusterCutsSignal_%d_%s",i,names[i].Data()),
           Form("EMCal %d, %s, embedded signal",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]);
          fhEMCALClusterCutsESignal[i]->SetYTitle("# clusters");
          fhEMCALClusterCutsESignal[i]->SetXTitle("#it{E} (GeV)");
          if ( fHistoPtDependent )
            fhEMCALClusterCutsESignal[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fOutputContainer->Add(fhEMCALClusterCutsESignal[i]);
        }
      }
      else
      {
        fhEMCALClusterCutsECen[i] = new TH2F
        (Form("hEMCALReaderClusterCutsCen_%d_%s",i,names[i].Data()),
         Form("EMCal %d, %s",i,names[i].Data()),
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
         100,0,100);
        fhEMCALClusterCutsECen[i]->SetZTitle("# clusters");
        fhEMCALClusterCutsECen[i]->SetXTitle("#it{E} (GeV)");
        if ( fHistoPtDependent )
          fhEMCALClusterCutsECen[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhEMCALClusterCutsECen[i]->SetYTitle("Centrality (%)");
        fOutputContainer->Add(fhEMCALClusterCutsECen[i]);
        
        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] ) //&&
            //!fSelectEmbeddedClusters && !fAcceptMCPromptPhotonOnly )
        {
          fhEMCALClusterCutsECenSignal[i] = new TH2F
          (Form("hEMCALReaderClusterCutsCenSignal_%d_%s",i,names[i].Data()),
           Form("EMCal %d, %s, embedded signal",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
           100,0,100);
          fhEMCALClusterCutsECenSignal[i]->SetZTitle("# clusters");
          fhEMCALClusterCutsECenSignal[i]->SetXTitle("#it{E} (GeV)");
          if ( fHistoPtDependent )
            fhEMCALClusterCutsECenSignal[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhEMCALClusterCutsECenSignal[i]->SetYTitle("Centrality (%)");
          fOutputContainer->Add(fhEMCALClusterCutsECenSignal[i]);
        }
      }
    }
    
    fhEMCALClusterTimeE  = new TH2F 
    ("hEMCALReaderTimeE","#it{time}_{cluster} vs #it{E}_{cluster} after cuts", 125,0,250,2000,-1000,1000);
    fhEMCALClusterTimeE->SetXTitle("#it{E}_{cluster} (GeV)");
    if ( fHistoPtDependent )
       fhEMCALClusterTimeE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhEMCALClusterTimeE->SetYTitle("#it{time}_{cluster} (ns)");
    fOutputContainer->Add(fhEMCALClusterTimeE);

    fhEMCALClusterDisToBadE  = new TH2F 
    ("hEMCALReaderDistToBadE","Distance to bad cell vs #it{E}_{cluster}", 50,0,50,100,0,20);
    fhEMCALClusterDisToBadE->SetXTitle("#it{E}_{cluster} (GeV)");
    if ( fHistoPtDependent ) 
       fhEMCALClusterDisToBadE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhEMCALClusterDisToBadE->SetYTitle("Distance to bad cell");
    fOutputContainer->Add(fhEMCALClusterDisToBadE);
    
    fhEMCALClusterEtaPhi  = new TH2F 
    ("hEMCALReaderEtaPhi","#eta vs #varphi",80,-2, 2,100, 0,10);
    // Very open limits to check problems
    fhEMCALClusterEtaPhi->SetXTitle("#eta");
    fhEMCALClusterEtaPhi->SetYTitle("#varphi (rad)");
    fOutputContainer->Add(fhEMCALClusterEtaPhi);    
    
    fhEMCALClusterEtaPhiFidCut  = new TH2F 
    ("hEMCALReaderEtaPhiFidCut","#eta vs #varphi after fidutial cut",80,-2, 2,100, 0,10);
    fhEMCALClusterEtaPhiFidCut->SetXTitle("#eta");
    fhEMCALClusterEtaPhiFidCut->SetYTitle("#varphi (rad)");
    fOutputContainer->Add(fhEMCALClusterEtaPhiFidCut);
    
    if ( fRemoveLEDEvents > 1 )
    {
      fhEMCALNSumEnCellsPerSM = new TH2F 
      ("hEMCALNSumEnCellsPerSM","Total number of cells and energy in any SM",144,0,1152,250,0,5000);
      fhEMCALNSumEnCellsPerSM->SetXTitle("#it{n}_{cells}^{SM}");
      fhEMCALNSumEnCellsPerSM->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
      fOutputContainer->Add(fhEMCALNSumEnCellsPerSM);
      
      fhEMCALNSumEnCellsPerSMAfter = new TH2F 
      ("hEMCALNSumEnCellsPerSMAfter","Total number of cells and energy in any SM, after LED SM event rejection",144,0,1152,250,0,5000);
      fhEMCALNSumEnCellsPerSMAfter->SetXTitle("#it{n}_{cells}^{SM}");
      fhEMCALNSumEnCellsPerSMAfter->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
      fOutputContainer->Add(fhEMCALNSumEnCellsPerSMAfter);
      
      if ( fRemoveLEDStripEvents )
      {
        fhEMCALNSumEnCellsPerSMAfterStripCut = new TH2F 
        ("hEMCALNSumEnCellsPerSMAfterStripCut","Total number of cells and energy in any SM, after LED SM and strip event rejection ",144,0,1152,250,0,5000);
        fhEMCALNSumEnCellsPerSMAfterStripCut->SetXTitle("#it{n}_{cells}^{SM}");
        fhEMCALNSumEnCellsPerSMAfterStripCut->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
        fOutputContainer->Add(fhEMCALNSumEnCellsPerSMAfterStripCut);
      }
    }
    
    if ( fRemoveLEDStripEvents && fRemoveLEDEvents > 0 )
    {
      fhEMCALNSumEnCellsPerStrip = new TH2F 
      ("hEMCALNSumEnCellsPerStrip","Total number of cells and energy in any strip, after LED SM event rejection",48,0,48,100,0,500);
      fhEMCALNSumEnCellsPerStrip->SetXTitle("#it{n}_{cells}^{strip}");
      fhEMCALNSumEnCellsPerStrip->SetYTitle("#Sigma #it{E}_{cells}^{strip} (GeV)");
      fOutputContainer->Add(fhEMCALNSumEnCellsPerStrip);
      
      fhEMCALNSumEnCellsPerStripAfter = new TH2F 
      ("hEMCALNSumEnCellsPerStripAfter","Total number of cells and energy in any strip, after LED SM event rejection",48,0,48,100,0,500);
      fhEMCALNSumEnCellsPerStripAfter->SetXTitle("#it{n}_{cells}^{strip}");
      fhEMCALNSumEnCellsPerStripAfter->SetYTitle("#Sigma #it{E}_{cells}^{strip} (GeV)");
      fOutputContainer->Add(fhEMCALNSumEnCellsPerStripAfter);
    }
    
    if ( fRemoveBadTriggerEventsFromEMCalTriggerMaker )
    {
      fhEMCALClusterBadTrigger = new TH1F
      ("hEMCALReaderClusterBadTrigger","Clusters in rejected triggered events",   
      fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]);
      fhEMCALClusterBadTrigger->SetYTitle("# clusters");
      fhEMCALClusterBadTrigger->SetXTitle("#it{E} (GeV)");
      fOutputContainer->Add(fhEMCALClusterBadTrigger);
      
      fhCentralityBadTrigger = new TH1F
      ("hCentralityBadTrigger","Rejected triggered events",   
      100, 0, 100);
      fhCentralityBadTrigger->SetYTitle("# eveents");
      fhCentralityBadTrigger->SetXTitle("centrality");
      fOutputContainer->Add(fhCentralityBadTrigger);
      
      fhEMCALClusterCentralityBadTrigger = new TH2F
      ("hEMCALReaderClusterCentralityBadTrigger","Clusters vs centrality in rejected triggered events",   
      fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],20,0,100);
      fhEMCALClusterCentralityBadTrigger->SetYTitle("centrality");
      fhEMCALClusterCentralityBadTrigger->SetXTitle("#it{E} (GeV)");
      fOutputContainer->Add(fhEMCALClusterCentralityBadTrigger);
    }
  }
  
  if ( fFillPHOS )
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
  
  if ( fFillCTS )
  {
    TString names[] = {"NoCut", "Status", "ESD_AOD", "TOF", "DCA","PtAcceptanceMult"};

    for(Int_t i = 0; i < 6; i++)
    {
      if ( names[i].Contains("Acceptance")  && !fFiducialCut    ) continue;
      if ( names[i] == "Status"             && !fTrackStatus    ) continue;
      if ( names[i] == "DCA"                && !fUseTrackDCACut ) continue;
      if ( names[i] == "TOF"                && !fAccessTrackTOF ) continue;
     
      if ( !fHistoCentDependent )
      {
        fhCTSTrackCutsPt[i] = new TH1F
        (Form("hCTSReaderTrackCuts_%d_%s",i,names[i].Data()),
         Form("CTS Cut %d, %s",i,names[i].Data()), 
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
        fhCTSTrackCutsPt[i]->SetYTitle("# tracks");
        fhCTSTrackCutsPt[i]->SetXTitle("#it{p}_{T} (GeV)");
        fOutputContainer->Add(fhCTSTrackCutsPt[i]);

        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
        {
          fhCTSTrackCutsPtSignal[i] = new TH1F
          (Form("hCTSReaderTrackCutsSignal_%d_%s",i,names[i].Data()),
           Form("CTS Cut %d, %s",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
          fhCTSTrackCutsPtSignal[i]->SetYTitle("# tracks");
          fhCTSTrackCutsPtSignal[i]->SetXTitle("#it{p}_{T} (GeV)");
          fOutputContainer->Add(fhCTSTrackCutsPtSignal[i]);
        }
      }
      else
      {
        fhCTSTrackCutsPtCen[i] = new TH2F
        (Form("hCTSReaderTrackCutsCen_%d_%s",i,names[i].Data()),
         Form("CTS Cut %d, %s",i,names[i].Data()), 
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
         100, 0, 100) ;
        fhCTSTrackCutsPtCen[i]->SetZTitle("# tracks");
        fhCTSTrackCutsPtCen[i]->SetXTitle("#it{p}_{T} (GeV)");
        fhCTSTrackCutsPtCen[i]->SetYTitle("Centrality (%)");
        fOutputContainer->Add(fhCTSTrackCutsPtCen[i]);

        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
        {
          fhCTSTrackCutsPtCenSignal[i] = new TH2F
          (Form("hCTSReaderTrackCutsCenSignal_%d_%s",i,names[i].Data()),
           Form("CTS Cut %d, %s",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
           100, 0, 100) ;
          fhCTSTrackCutsPtCenSignal[i]->SetZTitle("# tracks");
          fhCTSTrackCutsPtCenSignal[i]->SetXTitle("#it{p}_{T} (GeV)");
          fhCTSTrackCutsPtCenSignal[i]->SetYTitle("Centrality (%)");
          fOutputContainer->Add(fhCTSTrackCutsPtCenSignal[i]);
        }
      }
    }
  }
  
  if ( fComparePtHardAndJetPt )
  {
    fhPtHardPtJetPtRatio = new TH1F
    ("hPtHardPtJetPtRatio","Generated jet #it{p}_{T} / #it{p}_{T}^{hard}",100,0,10);
    fhPtHardPtJetPtRatio->SetYTitle("# events");
    fhPtHardPtJetPtRatio->SetXTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{hard}");
    fOutputContainer->Add(fhPtHardPtJetPtRatio);
  }

  if ( fComparePtHardAndPromptPhotonPt )
   {
     fhPtHardPromptPhotonPtRatio = new TH1F
     ("hPtHardPtPromptPhotonPtRatio","Generated prompt #gamma #it{p}_{T} / #it{p}_{T}^{hard}",100,0,10);
     fhPtHardPromptPhotonPtRatio->SetYTitle("# events");
     fhPtHardPromptPhotonPtRatio->SetXTitle("#it{p}_{T}^{prompt #gamma} / #it{p}_{T}^{hard}");
     fOutputContainer->Add(fhPtHardPromptPhotonPtRatio);
   }

  if ( fComparePtHardAndClusterPt )
  {
    fhPtHardEnClusterRatio = new TH1F
    ("hPtHardEnClusterRatio","Cluster energy / #it{p}_{T}^{hard}",100,0,10);
    fhPtHardEnClusterRatio->SetYTitle("# events");
    fhPtHardEnClusterRatio->SetXTitle("#it{E}_{cluster} / #it{p}_{T}^{hard}");
    fOutputContainer->Add(fhPtHardEnClusterRatio);

    if ( fHistoCentDependent )
    {
      fhPtHardEnClusterCenRatio = new TH2F
      ("hPtHardEnClusterCenRatio","Cluster energy / #it{p}_{T}^{hard}",100,0,10,20,0,100);
      fhPtHardEnClusterCenRatio->SetYTitle("Centrality (%)");
      fhPtHardEnClusterCenRatio->SetXTitle("#it{E}_{cluster} / #it{p}_{T}^{hard}");
      fOutputContainer->Add(fhPtHardEnClusterCenRatio);
    }
  }

  if ( fUseEventCutsClass )
    fEventCuts.AddQAplotsToList(fOutputContainer); 
  
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
  if ( fEMCALNCellsCutEnDepEnMin < 200 )
  {
    snprintf(onePar,buffersize,"For Emin>%2.1f, ncell>%2.2f+ %2.2fE; ",
             fEMCALNCellsCutEnDepEnMin,fEMCALNCellsCutEnDepConstant,fEMCALNCellsCutEnDepSlope) ;
    parList+=onePar ;
  }
  snprintf(onePar,buffersize,"N cells diff TCard: E > %2.2f, Ecell > %1.2f; ",fEMCALHighEnergyNdiffCut,fEMCALMinCellEnNdiffCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"EMC time cut single window (%2.2f,%2.2f); ",fEMCALTimeCutMin,fEMCALTimeCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Check: calo fid cut %d; ",fCheckFidCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Track: status %d, SPD hit %d; ITS cluster >= %d; ITS chi2 > %2.1f; TPC cluster >= %d; TPC chi2 > %2.1f ",
           (Int_t) fTrackStatus,  fSelectSPDHitTracks,
           fSelectMinITSclusters, fSelectMaxChi2PerITScluster,
           fSelectMinTPCclusters, fSelectMaxChi2PerTPCcluster) ;
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
  
  snprintf(onePar,buffersize,"Select fired trigger %s; Remove Bad trigger event %d, unmatched %d; Accept fastcluster %d; Trigger maker: bad %d, name %s",
          fFiredTriggerClassName.Data(), fRemoveBadTriggerEvents, fRemoveUnMatchedTriggers, fAcceptFastCluster,
           fRemoveBadTriggerEventsFromEMCalTriggerMaker, fEMCalTriggerMakerDecissionContainerName.Data());
  parList+=onePar ;
  
  if ( fRemoveLEDEvents > 0 )
  {
    snprintf(onePar,buffersize,"Remove LED %d, %2.1f < Ecell < %1.2f : SM - nCell > %d, Sum E > %2.0f; SM3 - nCell < %d, Sum E < %2.0f;",
             fRemoveLEDEvents, fLEDMinCellEnergy, fLEDMaxCellEnergy, 
             fLEDHighNCellsCutSM, fLEDHighEnergyCutSM, fLEDLowNCellsCutSM3, fLEDLowEnergyCutSM3);
    parList+=onePar ;
    
    if ( fRemoveLEDStripEvents > 0 )
    {
      snprintf(onePar,buffersize,"Remove LED strip? %d, with n strip > %d: "
               "Full SM, nCell > %d, Sum E > %2.0f; "
               "1/3 SM, nCell > %d, Sum E > %2.0f; "
               "SM3, nCell < %d, Sum E < %2.0f;",
               fRemoveLEDStripEvents    , fLEDEventMaxNumberOfStrips, 
               fLEDHighNCellsCutStrip[0], fLEDHighEnergyCutStrip[0], 
               fLEDHighNCellsCutStrip[1], fLEDHighEnergyCutStrip[1], 
               fLEDLowNCellsCutSM3Strip , fLEDLowEnergyCutSM3Strip);
      parList+=onePar ;
    }
  }
  
  
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
  
  if ( fComparePtHardAndJetPt )
  {
    snprintf(onePar,buffersize,"jet pt / pt hard < %2.1f; ",fPtHardAndJetPtFactor);
    parList+=onePar ;
  }
  
  if ( fComparePtHardAndPromptPhotonPt )
  {
    snprintf(onePar,buffersize,"prompt photon pt / pt hard < %2.1f; ",fPtHardAndPromptPhotonPtFactor);
    parList+=onePar ;
  }

  if ( fComparePtHardAndClusterPt )
  {
    snprintf(onePar,buffersize,"cluster pt / pt hard < %2.2f",fPtHardAndClusterPtFactor);
    parList+=onePar ;
  }
  
  snprintf(onePar,buffersize,"Centrality: Class %s, Option %d, Bin [%d,%d]; New centrality %d; Mult. PS %d; Event plane method %s; ", 
           fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1],
           fUseAliCentrality,fMultWithEventSel,fEventPlaneMethod.Data()) ;
  parList+=onePar ;
  
  return new TObjString(parList) ;
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
  // Activate debug level in reader
  if( fDebug >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(this->ClassName(),fDebug);
  
  // Activate debug level in AliAnaWeights
  if( fWeightUtils->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(fWeightUtils->ClassName(), fWeightUtils->GetDebug());
  
  // Activate debug level in AliMCAnalysisUtils
  if( GetMCAnalysisUtils()->GetDebug() >= 0 )
  (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetMCAnalysisUtils()->ClassName(),GetMCAnalysisUtils()->GetDebug());
  
  //printf("Debug levels: Reader %d, Neutral Sel %d, Iso %d\n",fDebug,GetMCAnalysisUtils()->GetDebug(),fWeightUtils->GetDebug());

  if ( fAcceptMCPromptPhotonOnly && fRejectMCFragmentationPhoton )
    AliFatal("Prompt and Frag photon filtering cannot be activated at the same time!");
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
  
  fEMCALNCellsCutEnDepEnMin    = 10000 ;  
  fEMCALNCellsCutEnDepConstant = 4.9;    // for emin=40 GeV
  fEMCALNCellsCutEnDepSlope    = 0.04;   // for emin=40 GeV
  
  // For clusters with energy above fEMCALHighEnergyNdiffCut, count cells
  // with E > fEMCALMinCellEnNdiffCut in different T-Card than highest energy cell.
  // Reject if 0.
  fEMCALHighEnergyNdiffCut = 50.;
  fEMCALMinCellEnNdiffCut  = 0.;
  
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
  
  fDeltaAODFileName      = "deltaAODPartCorr.root";
  fFiredTriggerClassName = "";
  fEventTriggerMask      = AliVEvent::kAny;
  fMixEventTriggerMask   = AliVEvent::kAnyINT;
  fEventTriggerAtSE      = kTRUE; // Use only events that pass event selection at SE base class
 
  fEmbeddedEvent[0] = kFALSE;
  fEmbeddedEvent[1] = kFALSE;
  
  fAcceptFastCluster = kTRUE;
  fEventType         = -1;
  
  fRemoveLEDEvents = 0;
  fLEDHighEnergyCutSM = 500.; fLEDHighNCellsCutSM = 100;
  fLEDLowEnergyCutSM3 = 20. ; fLEDLowNCellsCutSM3 = 20;
  fLEDMinCellEnergy = 0.5;
  fLEDMaxCellEnergy = 15.;
  
  fRemoveLEDStripEvents     = 0  ;
  fLEDEventMaxNumberOfStrips= 0  ; 
  fLEDHighEnergyCutStrip[0] = 80 ; fLEDHighEnergyCutStrip[1] = 55 ; 
  fLEDHighNCellsCutStrip[0] = 24 ; fLEDHighNCellsCutStrip[1] = 15 ;
  fLEDLowEnergyCutSM3Strip  = 100; // open
  fLEDLowNCellsCutSM3Strip  = 100; // open
  
  //We want tracks fitted in the detectors:
  //fTrackStatus=AliESDtrack::kTPCrefit;
  //fTrackStatus|=AliESDtrack::kITSrefit;
  fTrackStatus     = 0;
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0;
  fV0Mul[0] = 0;   fV0Mul[1] = 0;
  
  fZvtxCut   = 10.;
  
  fNMixedEvent = 1;
  
  fPtHardAndJetPtFactor     = 4. ;
  fPtHardAndClusterPtFactor = 1.5;
  fPtHardAndPromptPhotonPtFactor = 2.;
  
  //Centrality
  fUseAliCentrality = kFALSE;
  fMultWithEventSel = kTRUE;
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
  
  fEnergyHistogramNbins    = 500 ;
  fEnergyHistogramLimit[0] = 0.  ;
  fEnergyHistogramLimit[1] = 250.;
  
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
  
  fTimeStampEventCTPBCCorrMin = -1;
  fTimeStampEventCTPBCCorrMax = 1e12;
  
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
  
  fEMCalTriggerMakerDecissionContainerName = "EmcalTriggerDecision";
  
  //Jets
  fInputNonStandardJetBranchName = "jets";
  fFillInputNonStandardJetBranch = kFALSE;
  //if(!fNonStandardJets) fNonStandardJets = new TClonesArray("AliAODJet",100);
  if(!fNonStandardJets) fNonStandardJets = new TClonesArray("AliEmcalJet",100);
  fInputBackgroundJetBranchName = "jets";
  fFillInputBackgroundJetBranch = kFALSE; 
  //if(!fBackgroundJets) fBackgroundJets = new AliAODJetEventBackground();
  if(!fBackgroundJets) fBackgroundJets = new TClonesArray("AliEmcalJet",100);

  fSmearShowerShapeWidth = 0.005;
  fSmearNLMMin = 1;
  fSmearNLMMax = 1;
  
  fMCUtils     = new AliMCAnalysisUtils() ;

  fWeightUtils = new AliAnaWeights() ;
  fEventWeight = 1 ;
      
  fTrackMultNPtCut = 8;
  fTrackMultPtCut[0] = 0.15; fTrackMultPtCut[1] = 0.5;  fTrackMultPtCut[2] = 1.0; 
  fTrackMultPtCut[3] = 2.0 ; fTrackMultPtCut[4] = 4.0;  fTrackMultPtCut[5] = 6.0;  
  fTrackMultPtCut[6] = 8.0 ; fTrackMultPtCut[7] = 10.;  
  fTrackMultPtCut[8] = 15.0; fTrackMultPtCut[9] = 20.;  
  
  for(Int_t ism = 0; ism < 22; ism++) fScaleFactorPerSM[ism] = 1. ;    
}

//__________________________________________________________________________
/// Select the cluster depending on a time window, either a simple
/// range or a parametrization depending on the energy.
//__________________________________________________________________________
Bool_t AliCaloTrackReader::IsInTimeWindow(Double_t tof, Float_t energy) const
{  
  // Parametrized cut depending on E
  if ( fUseParamTimeCut )
  {
    Float_t minCut= fEMCALParamTimeCutMin[0]+fEMCALParamTimeCutMin[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMin[2])/fEMCALParamTimeCutMin[3]);
    Float_t maxCut= fEMCALParamTimeCutMax[0]+fEMCALParamTimeCutMax[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMax[2])/fEMCALParamTimeCutMax[3]);
    //printf("tof %f, minCut %f, maxCut %f\n",tof,minCut,maxCut);
    if ( tof < minCut || tof > maxCut )  return kFALSE ;
  }
  
  //In any case, the time should to be larger than the fixed window ...
  if ( tof < fEMCALTimeCutMin  || tof > fEMCALTimeCutMax )  return kFALSE ;
  
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

  fCurrentParIndex = 0;
  if ( IsParRun() )
  {
    ULong64_t globalEventID = 
    (ULong64_t)fInputEvent->GetBunchCrossNumber() + 
    (ULong64_t)fInputEvent->GetOrbitNumber () * (ULong64_t)3564 + 
    (ULong64_t)fInputEvent->GetPeriodNumber() * (ULong64_t)59793994260;
    
    for(Short_t ipar=0;ipar<GetCaloUtils()->GetEMCALRecoUtils()->GetNPars();ipar++)
    {
      if(globalEventID >= GetCaloUtils()->GetEMCALRecoUtils()->GetGlobalIDPar(ipar)) 
      {
        fCurrentParIndex++;
      }
    }
  }
  
  GetCaloUtils()->GetEMCALRecoUtils()->SetCurrentParNumber(fCurrentParIndex);
  
  //fCurrentFileName = TString(currentFileName);
  if ( !fInputEvent )
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
  
  if ( fTimeStampEventSelect )
  {
    Int_t timeStamp = fInputEvent->GetTimeStamp();
    Float_t timeStampFrac = 1.*(timeStamp-fTimeStampRunMin) / (fTimeStampRunMax-fTimeStampRunMin);
    
    //printf("stamp0 %d, max0 %d, frac %f\n", timeStamp-fTimeStampRunMin,fTimeStampRunMax-fTimeStampRunMin, timeStampFrac);
    
    if(timeStampFrac < fTimeStampEventFracMin || timeStampFrac > fTimeStampEventFracMax) return kFALSE;
    
    AliDebug(1,"Pass Time Stamp rejection");
    
    fhNEventsAfterCut->Fill(9.5);
  }

  if ( fDataType==kESD && fTimeStampEventCTPBCCorrExclude )
  {
    AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
    if(esd)
    {
      Int_t timeStamp = esd->GetTimeStampCTPBCCorr();
      
      if(timeStamp > fTimeStampEventCTPBCCorrMin && timeStamp <= fTimeStampEventCTPBCCorrMax) return kFALSE;
    }
    
    AliDebug(1,"Pass Time Stamp CTPBCCorr rejection");
    
    fhNEventsAfterCut->Fill(9.5);
  }

  
  //------------------------------------------------------
  // Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------
  
  FillVertexArray();
  
  if ( fUseEventsWithPrimaryVertex )
  {
    if ( !CheckForPrimaryVertex() )             return kFALSE; // algorithm in ESD/AOD Readers

    fhNEventsAfterCut->Fill(10.5);

    if( TMath::Abs(fVertex[0][0] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][1] ) < 1.e-6 &&
        TMath::Abs(fVertex[0][2] ) < 1.e-6    ) return kFALSE;

    AliDebug(1,"Pass primary vertex/null rejection");
    
    fhNEventsAfterCut->Fill(11.5);
  }

  //Reject events with Z vertex too large, only for SE analysis, if not, cut on the analysis code
  if(!GetMixedEvent() && TMath::Abs(fVertex[0][2]) > fZvtxCut) return kFALSE;
  
  fhNEventsAfterCut->Fill(12.5);

  AliDebug(1,"Pass z vertex rejection");

  
  //printf("Reader : IsPileUp %d, Multi %d\n",IsPileUpFromSPD(),fInputEvent->IsPileupFromSPDInMultBins());
  
  if ( fDoPileUpEventRejection )
  {
    // Do not analyze events with pileup
    Bool_t bPileup = IsPileUpFromSPD();
    //IsPileupFromSPDInMultBins() // method to try
    //printf("pile-up %d, %d, %2.2f, %2.2f, %2.2f, %2.2f\n",bPileup, (Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
    if(bPileup) return kFALSE;
    
    AliDebug(1,"Pass Pile-Up event rejection");
    
    fhNEventsAfterCut->Fill(13.5);
  }
  
  if ( fDoV0ANDEventSelection )
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
    
    fhNEventsAfterCut->Fill(14.5);
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
    
    AliDebug(1,Form("Centrality %d in [%d,%d]?", cen, fCentralityBin[0], fCentralityBin[1]));

    if ( cen >= fCentralityBin[1] || cen <  fCentralityBin[0]   )  return kFALSE; //reject events out of bin.
    
    AliDebug(1,"Pass centrality rejection");
    
    fhNEventsAfterCut->Fill(15.5);
  }

  //----------------------------------------------------------------
  // MC events selections
  //----------------------------------------------------------------
  if ( GetMC() )
  {
    //----------------------------------------------------------------
    // Get the event headers
    //----------------------------------------------------------------
    
    // Main header
    // Init it first to 0 to tell the method to recover it.
    fGenEventHeader = 0;
    fGenEventHeader = GetGenEventHeader();

    if ( fGenEventHeader )
    {
      AliDebug(1,Form("Selected event header class <%s>, name <%s>; cocktail %p",
                      fGenEventHeader->ClassName(),
                      fGenEventHeader->GetName(),
                      GetMC()->GetCocktailList()));
    }
   
    //----------------------------------------------------------------
    // Reject the event if the event header name is not
    // the one requested among the possible generators.
    // Needed in case of cocktail MC generation with multiple options.
    //----------------------------------------------------------------
    if ( fMCGenerEventHeaderToAccept!="" ) 
    {
      if ( !fGenEventHeader ) return kFALSE;
      
      AliDebug(1,"Pass Event header selection");
      
      fhNEventsAfterCut->Fill(16.5);
    }
    
    // Pythia header
    TString pyGenName       = ""; 
    TString pyProcessName   = "";  
    Int_t   pyProcess       = 0;
    Int_t   pyFirstGenPart  = 0; 
    Int_t   pythiaVersion   = 0;
    
    // Init it first to 0 to tell the method to recover it.
    if ( fCheckPythiaEventHeader )
      fGenPythiaEventHeader = 
      GetMCAnalysisUtils()->GetPythiaEventHeader(GetMC(),fMCGenerEventHeaderToAccept,
                                                 pyGenName,pyProcessName,pyProcess,pyFirstGenPart,pythiaVersion);

    if ( pyProcessName != "Gamma-Jet" && fAcceptMCPromptPhotonOnly    ) 
      AliFatal("Not a pythia gamma-jet process, set reader->SwitchOffMCPromptPhotonsSelection()");
    
    if ( pyProcessName != "Jet-Jet"   && fRejectMCFragmentationPhoton ) 
      AliFatal("Not a pythia jet-jet process, set reader->SwitchOffMCFragmentationPhotonsRejection()");
    
    if ( fGenPythiaEventHeader )
    {
      AliDebug(2,Form("Pythia v%d name <%s>, process %d <%s>, first generated particle %d",
                   pythiaVersion, pyGenName.Data(), pyProcess, pyProcessName.Data(), pyFirstGenPart));
      
      //---------------------------------------------------------------------------
      // In case of analysis of events with jets, skip those with jet pt > 5 pt hard
      // To be used on for MC data in pT hard bins
      //---------------------------------------------------------------------------
      
      if ( fComparePtHardAndJetPt )
      {
        if ( !ComparePtHardAndJetPt(pyProcess, pyProcessName) ) return kFALSE ;
        
        AliDebug(1,"Pass Pt Hard - Jet rejection");
        
        fhNEventsAfterCut->Fill(17.5);
      }
      
      if ( fComparePtHardAndClusterPt )
      {
        if ( !ComparePtHardAndClusterPt(pyProcess, pyProcessName) ) return kFALSE ;
        
        AliDebug(1,"Pass Pt Hard - Cluster rejection");
        
        if ( !fComparePtHardAndPromptPhotonPt ) // avoid double counting in next filling
          fhNEventsAfterCut->Fill(18.5);
      }

      if ( fComparePtHardAndPromptPhotonPt )
      {
        if ( !ComparePtHardAndPromptPhotonPt(pyProcess, pyProcessName) ) return kFALSE ;

        AliDebug(1,"Pass Pt Hard - Prompt photon rejection");

        fhNEventsAfterCut->Fill(18.5);
      }

    } // pythia header
  } // MC
  
  //------------------------------------------------------------------
  // Recover the weight assigned to the event, if provided
  // right now only for pT-hard bins and centrality dependent weights
  //------------------------------------------------------------------
  if ( fWeightUtils->IsWeightSettingOn() )
  {
    fWeightUtils->SetCentrality(cen);
    
    fWeightUtils->SetPythiaEventHeader(fGenPythiaEventHeader);
      
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
  
  if ( fFillCTS )
  {
    FillInputCTS();
    
    // Accept events with at least one track
    if ( fTrackMult[0] == 0 && fDoRejectNoTrackEvents ) 
      return kFALSE ;
    
    AliDebug(1,"Pass rejection of null track events");

    fhNEventsAfterCut->Fill(19.5);    
  }
  
  if ( fDoVertexBCEventSelection )
  {
    if ( fVertexBC != 0 && fVertexBC != AliVTrack::kTOFBCNA ) 
      return kFALSE ;
    
    AliDebug(1,"Pass rejection of events with vertex at BC!=0");
    
    fhNEventsAfterCut->Fill(20.5);
  }
  
  //-----------------------------------
  // Check the event cuts implemented in AliEventCuts class
  // Right after other event cuts applied here 
  // Some might be also applied in AliEventCuts
  //-----------------------------------
  if ( fUseEventCutsClass )
  {
    Bool_t accept = fEventCuts.AcceptEvent(GetInputEvent());
    
    if ( !accept ) return kFALSE;
    
    AliDebug(1,"Pass AliEventCuts!");
    
    fhNEventsAfterCut->Fill(21.5);
  }
  
  //-----------------------------------
  // Get and filter calorimeter data
  //-----------------------------------
  
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
/// \return  AliCentrality pointer object
//__________________________________________________
AliCentrality*    AliCaloTrackReader::GetCentrality() const 
{ 
  if ( fDataType == kMC ) return 0x0; 
  
  AliVEvent * event = NULL; 

  // In case of analysis of pure MC event used in embedding 
  // get bkg PbPb event since cuts for embedded event are based on 
  // data centrality and not on pp simu with no centrality
  if ( fEmbeddedEvent[1] ) 
    event = ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetEvent();
  else 
    event = fInputEvent;
  
  return event->GetCentrality() ;
} 
//__________________________________________________
/// \return  AliMultiplicity pointer object
//__________________________________________________
AliMultSelection* AliCaloTrackReader::GetMultSelCen() const 
{ 
  if ( fDataType == kMC ) return 0x0; 
  
  AliVEvent * event = NULL; 

  // In case of analysis of pure MC event used in embedding 
  // get bkg PbPb event since cuts for embedded event are based on 
  // data centrality and not on pp simu with no centrality
  if ( fEmbeddedEvent[1] ) 
    event = ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetEvent();
  else 
    event = fInputEvent;
  
  return (AliMultSelection * ) event->FindListObject("MultSelection") ; 
} 


//__________________________________________________
/// \return Current event centrality bin. 
/// Different percentile options and centrality class can be requested.
//__________________________________________________
Float_t AliCaloTrackReader::GetEventCentralityF() const
{  
  if(fUseAliCentrality)
  {
    if ( !GetCentrality() ) return -1;
    
    AliDebug(1,Form("Cent. Percentile: V0M %2.2f, CL0 %2.2f, CL1 %2.2f; selected class %s", 
                    GetCentrality()->GetCentralityPercentile("V0M"), 
                    GetCentrality()->GetCentralityPercentile("CL0"), 
                    GetCentrality()->GetCentralityPercentile("CL1"), 
                    fCentralityClass.Data()));
    
    if     (fCentralityOpt == 100) return GetCentrality()->GetCentralityPercentile(fCentralityClass); // 100 bins max
    else if(fCentralityOpt ==  10) return GetCentrality()->GetCentralityClass10(fCentralityClass);// 10 bins max
    else if(fCentralityOpt ==  20) return GetCentrality()->GetCentralityClass5(fCentralityClass); // 20 bins max
    else
    {
      AliInfo(Form("Unknown centrality option %d, use 10, 20 or 100",fCentralityOpt));
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
    
    return GetMultSelCen()->GetMultiplicityPercentile(fCentralityClass, fMultWithEventSel); // returns centrality only for events used in calibration
    
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
    AliDebug(1,Form("Bad EP for <Q> method : %f",ep));
    return -1000;
  }
  else if(GetEventPlaneMethod().Contains("V0")  )
  {
    if((ep > TMath::Pi()/2 || ep < -TMath::Pi()/2))
    {
      AliDebug(1,Form("Bad EP for <%s> method : %f",GetEventPlaneMethod().Data(), ep));
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
  
  if ( !fMixedEvent )
  { // Single event analysis
    if ( fDataType != kMC )
    {
      AliVEvent * event = NULL;
      // In case of analysis of pure MC event used in embedding 
      // get bkg PbPb event since cuts for embedded event are based on 
      // data vertex and not on pp simu vertex
      if ( fEmbeddedEvent[1] ) // Input event is MC
        event = ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetEvent();
      else 
        event = fInputEvent;
      
      if ( event->GetPrimaryVertex() )
      {
        event->GetPrimaryVertex()->GetXYZ(fVertex[0]);
      }
      else
      {
        AliWarning("NULL primary vertex");
        fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
      }//Primary vertex pointer do not exist
      
    } 
    else
    {// MC read event
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
    
    AliDebug(1,Form("Single Event Vertex : %f,%f,%f",
                    fVertex[0][0],fVertex[0][1],fVertex[0][2]));
  } 
  else
  { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++)
    {
      if ( fMixedEvent->GetVertexOfEvent(iev) )
        fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
      else
        AliWarning("No vertex found");
      
      AliDebug(1,Form("Multi Event %d Vertex : %f,%f,%f",
                      iev,fVertex[iev][0],fVertex[iev][1],fVertex[iev][2]));
    }
  }
}

//_____________________________________
/// Fill the array with Central Tracking System (CTS) 
/// filtered tracks. Filtering done in FillInputCTSSelectTracks()
//_____________________________________
void AliCaloTrackReader::FillInputCTS()
{  
  AliDebug(1,"Begin");
    
  Int_t nTracks = fInputEvent->GetNumberOfTracks() ;
  
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
  if ( fRecalculateVertexBC ) fVertexBC = AliVTrack::kTOFBCNA;
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++)
  {
    AliVTrack * track = (AliVTrack*)fInputEvent->GetTrack(itrack) ; 
    
    FillInputCTSSelectTrack(track, itrack, bc0);
  }
  
  // Add embedded tracks from external event
  // Only if input event is just data and not already external
  if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
  {
    AliVEvent * externalEvent = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent();
    if ( externalEvent )
    {
      for (Int_t jtrack =  0; jtrack < externalEvent->GetNumberOfTracks() ; jtrack++)
      {
        AliVTrack * extTrack = (AliVTrack*) externalEvent->GetTrack(jtrack) ; 
        
        FillInputCTSSelectTrack(extTrack, jtrack, bc0);
      } // track loop
    } 
    else 
      AliInfo(Form("No external event for embed mc %d embed data %d",fEmbeddedEvent[0], fEmbeddedEvent[1]));
  }
  
  if( fRecalculateVertexBC && (fVertexBC == 0 || fVertexBC == AliVTrack::kTOFBCNA))
  {
    if ( bc0 ) fVertexBC = 0 ;
    else       fVertexBC = AliVTrack::kTOFBCNA ;
  }
  
  AliDebug(1,Form("CTS entries %d, input tracks %d, multipliticy %d", 
                  fCTSTracks->GetEntriesFast(), nTracks, fTrackMult[0]));
}

//_______________________________________________________________________________
/// Select the tracks based on kinematic cuts, DCA, 
/// re-fit status and timing cuts are applied. 
/// 
/// Other more ESD/AOD dependent cuts are applied in *SelectTrack()* method,
/// see AliCaloTrackAODReader and AliCaloTrackESDReader.
///
/// \param track: AliVTrack pointer
/// \param iclus: track index, only needed in case of mixing frame (not used recently)
/// \param bc0: bunch crossing bool, at least one vertex track at bc=0
///
/// Method called by *FillInputCTS()*
//_______________________________________________________________________________
void AliCaloTrackReader::FillInputCTSSelectTrack(AliVTrack * track, Int_t itrack,
                                                 Bool_t & bc0)
{
  if ( !AcceptParticleMCLabel( TMath::Abs(track->GetLabel()) ) ) return ;

  Int_t cen = GetEventCentrality();
  Bool_t fillEmbedSignalTrack = kFALSE;
  if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1]  && track->GetLabel() >=0 )
    fillEmbedSignalTrack = kTRUE;

  if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [0]->Fill(track->Pt());
  else                        fhCTSTrackCutsPtCen[0]->Fill(track->Pt(),cen);
    
  if ( fillEmbedSignalTrack )
  {
    if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [0]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCenSignal[0]->Fill(track->Pt(),cen);
  }

  //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
  ULong_t status = track->GetStatus();
  
  if ( fTrackStatus && !((status & fTrackStatus) == fTrackStatus) )
    return ;
  
  if ( fTrackStatus )
  { 
    if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [1]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCen[1]->Fill(track->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [1]->Fill(track->Pt());
      else                        fhCTSTrackCutsPtCenSignal[1]->Fill(track->Pt(),cen);
    }
  }
  
  //-------------------------
  // Select the tracks depending on cuts of AOD or ESD
  Double_t pTrack[3] = {0,0,0};
  if ( !SelectTrack(track, pTrack) ) return ;
  
  if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [2]->Fill(track->Pt());
  else                        fhCTSTrackCutsPtCen[2]->Fill(track->Pt(),cen);
  
  if ( fillEmbedSignalTrack )
  {
    if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [2]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCenSignal[2]->Fill(track->Pt(),cen);
  }

  //-------------------------
  // TOF cuts
  Bool_t okTOF  = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
  Double_t tof  = -1000;
  Int_t trackBC = -1000 ;

  if ( fAccessTrackTOF )
  {
    if ( okTOF )
    {
      trackBC = track->GetTOFBunchCrossing(GetInputEvent()->GetMagneticField());
      SetTrackEventBC(trackBC+9);
      
      tof = track->GetTOFsignal()*1e-3;
      
      // After selecting tracks with small DCA, pointing to vertex, set vertex BC depeding on tracks BC
      if ( fRecalculateVertexBC )
      {
        if     (trackBC != 0 && trackBC != AliVTrack::kTOFBCNA) fVertexBC = trackBC;
        else if(trackBC == 0)                                   bc0       = kTRUE;
      }
      
      //In any case, the time should to be larger than the fixed window ...
      if( fUseTrackTimeCut && (trackBC !=0 || tof < fTrackTimeCutMin  || tof > fTrackTimeCutMax) )
      {
        //printf("Remove track time %f and bc = %d\n",tof,trackBC);
        return ;
      }
      //else printf("Accept track time %f and bc = %d\n",tof,trackBC);
    }
    
    if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [3]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCen[3]->Fill(track->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [3]->Fill(track->Pt());
      else                        fhCTSTrackCutsPtCenSignal[3]->Fill(track->Pt(),cen);
    }
  }
  
  //---------------------
  // DCA cuts
  //
  fMomentum.SetPxPyPzE(pTrack[0],pTrack[1],pTrack[2],0);
  
  if ( fUseTrackDCACut )
  {      
    Float_t dcaTPC =-999;
    //In case of AODs, TPC tracks cannot be propagated back to primary vertex,
    if( fDataType == kAOD ) dcaTPC = ((AliAODTrack*) track)->DCA();
    
    //normal way to get the dca, cut on dca_xy
    if ( dcaTPC==-999 )
    {
      Double_t dca[2]   = {1e6,1e6};
      Double_t covar[3] = {1e6,1e6,1e6};
      Bool_t okDCA = track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),GetInputEvent()->GetMagneticField(),100.,dca,covar);
      if (  okDCA ) okDCA = AcceptDCA(fMomentum.Pt(),dca[0]);
      if ( !okDCA )
      {
        //printf("AliCaloTrackReader::FillInputCTS() - Reject track pt %2.2f, dca_xy %2.4f\n",fMomentum.Pt(),dca[0]);
        return ;
      }
    }
    
    if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [4]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCen[4]->Fill(track->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [4]->Fill(track->Pt());
      else                        fhCTSTrackCutsPtCenSignal[4]->Fill(track->Pt(),cen);
    }
  }// DCA cuts
  
  //-------------------------
  // Kinematic/acceptance cuts
  //
  // Count the tracks in eta < 0.9 and different pT cuts
  Float_t ptTrack = fMomentum.Pt();
  if ( TMath::Abs(track->Eta()) <  fTrackMultEtaCut ) 
  {
    for(Int_t iptCut = 0; iptCut < fTrackMultNPtCut; iptCut++ )
    {
      if ( ptTrack > fTrackMultPtCut[iptCut] ) 
      {
        fTrackMult [iptCut]++;
        fTrackSumPt[iptCut]+=ptTrack;
      }
    }
  }
  
  if ( fCTSPtMin > ptTrack || fCTSPtMax < ptTrack ) return ;
  
  // Check effect of cuts on track BC
  if ( fAccessTrackTOF && okTOF ) SetTrackEventBCcut(trackBC+9);
  
  if ( fCheckFidCut ) 
  {
    if ( !fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kCTS) ) return;
    
    if ( !fHistoCentDependent ) fhCTSTrackCutsPt   [5]->Fill(track->Pt());
    else                        fhCTSTrackCutsPtCen[5]->Fill(track->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSTrackCutsPtSignal   [5]->Fill(track->Pt());
      else                        fhCTSTrackCutsPtCenSignal[5]->Fill(track->Pt(),cen);
    }
  }
  
  // ------------------------------
  // Add selected tracks to array
  AliDebug(2,Form("Selected tracks pt %3.2f, phi %3.2f deg, eta %3.2f",
                  fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));
  
  fCTSTracks->Add(track);
  
  // TODO, check if remove
  if ( fMixedEvent )  track->SetID(itrack);
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
void AliCaloTrackReader::FillInputEMCALSelectCluster(AliVCluster * clus, Int_t iclus)
{
  // Accept clusters with the proper label, only applicable for MC
  //
  Int_t mclabel = clus->GetLabel();
  if ( mclabel >= 0 )  // -1 corresponds to noisy MC
  { 
    if ( !AcceptParticleMCLabel(clus->GetLabel()) ) return ;
  }
  
//  // If requested, accept only prompt photon clusters or
//  // reject fragmentation photon clusters
//  //
//  if ( fMC && (fAcceptMCPromptPhotonOnly || fRejectMCFragmentationPhoton) )
//  {
//    if ( mclabel < 0 && fAcceptMCPromptPhotonOnly ) return ;
//
//    Int_t tag = 0;
//    if ( mclabel >= 0 )
//      tag = GetMCAnalysisUtils()->CheckOrigin(mclabel, GetMC(),
//                                              GetNameOfMCEventHederGeneratorToAccept(),
//                                              clus->E());
//    if ( fAcceptMCPromptPhotonOnly &&
//        !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt       ) ) return ;
//
//    if ( fRejectMCFragmentationPhoton &&
//         (GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCFragmentation) ||
//          GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCISR)) ) return ;
//  }

  // TODO, not sure if needed anymore
  Int_t vindex = 0 ;
  if (fMixedEvent)
    vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
    
  clus->GetMomentum(fMomentum, fVertex[vindex]);
  Float_t energyOrMom = clus->E();
  if ( fHistoPtDependent ) energyOrMom = fMomentum.Pt();
  
  Int_t cen = GetEventCentrality();
  Bool_t fillEmbedSignalCluster = kFALSE;
  if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] && !fSelectEmbeddedClusters && // && !fAcceptMCPromptPhotonOnly
       clus->GetNLabels() > 0 && clus->GetLabel() >=0 )
    fillEmbedSignalCluster = kTRUE;
  
  // No correction/cut applied yet
  if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [0]->Fill(energyOrMom);
  else                        fhEMCALClusterCutsECen[0]->Fill(energyOrMom,cen);
  
  if ( fillEmbedSignalCluster )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [0]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECenSignal[0]->Fill(energyOrMom,cen);
  }
  
  // Get the maximum cell energy, its SM number and its col, row location, needed in 
  // different places of this method, although not active by default, one can consider
  // deactivate this and only activate it when requiered.
  Int_t absIdMax= -1;
  Int_t iSupMod = -1;
  Int_t iphiMax = -1;
  Int_t ietaMax = -1;
  Bool_t shared = kFALSE;
  GetCaloUtils()->GetEMCALRecoUtils()->GetMaxEnergyCell(GetCaloUtils()->GetEMCALGeometry(), 
                                                        GetEMCALCells(),clus,absIdMax,iSupMod,
                                                        ietaMax,iphiMax,shared);
  
  //if( (fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10 )
  AliDebug(2,Form("Input cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f, nCells %d, SM %d",
                  fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta(), clus->GetNCells(), iSupMod));

  //---------------------------
  // Embedding case
  if ( fSelectEmbeddedClusters )
  {
    if ( clus->GetNLabels()==0 || clus->GetLabel() < 0 ) return;
    //else printf("Embedded cluster,  %d, n label %d label %d  \n",iclus,clus->GetNLabels(),clus->GetLabel());
  }

  //--------------------------------------
  // Apply some corrections in the cluster
  //
  if ( fRecalculateClusters )
  {
    // Recalibrate the cluster energy
    if ( GetCaloUtils()->IsRecalibrationOn() )
    {
      Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, GetEMCALCells());
      
      clus->SetE(energy);
      //printf("Recalibrated Energy %f\n",clus->E());
      
      GetCaloUtils()->RecalculateClusterShowerShapeParameters(GetEMCALCells(),clus);
      GetCaloUtils()->RecalculateClusterPID(clus);
      
      clus->GetMomentum(fMomentum, fVertex[vindex]);
      energyOrMom = clus->E();
      if ( fHistoPtDependent ) energyOrMom = fMomentum.Pt();
    } // recalculate E
    
    //Recalculate distance to bad channels, if new list of bad channels provided
    GetCaloUtils()->RecalculateClusterDistanceToBadChannel(GetEMCALCells(),clus);
    
    //Recalculate cluster position
    if ( GetCaloUtils()->IsRecalculationOfClusterPositionOn() )
    {
      GetCaloUtils()->RecalculateClusterPosition(GetEMCALCells(),clus);
      //clus->GetPosition(pos);
      //printf("After  Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
    }
    
    // Recalculate TOF
    if ( GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn() )
    {
      Double_t tof      = clus->GetTOF();
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
      
      //additional L1 phase shift
      if(GetCaloUtils()->GetEMCALRecoUtils()->IsL1PhaseInTimeRecalibrationOn())
      {
        GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(iSupMod,fInputEvent->GetBunchCrossNumber(), tof, fCurrentParIndex);
      }

      clus->SetTOF(tof);
      
    }// Time recalibration
  }
  
  // Check effect of corrections
  if ( fSelectEmbeddedClusters  || fRecalculateClusters )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [1]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[1]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
     {
       if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [1]->Fill(energyOrMom);
       else                        fhEMCALClusterCutsECenSignal[1]->Fill(energyOrMom,cen);
     }
  }
  
  //-----------------------------------------------------------------
  // Reject clusters with bad channels, close to borders and exotic
  //
  Bool_t goodCluster = GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,
                                                                          GetCaloUtils()->GetEMCALGeometry(),
                                                                          GetEMCALCells(),fInputEvent->GetBunchCrossNumber());
  
  if ( !goodCluster )
  {
    //if( (fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10 )
    AliDebug(1,Form("Bad cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));

    return;
  }
    
  // Check effect of bad cluster removal 
  if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [2]->Fill(energyOrMom);
  else                        fhEMCALClusterCutsECen[2]->Fill(energyOrMom,cen);
  
  if ( fillEmbedSignalCluster )
   {
     if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [2]->Fill(energyOrMom);
     else                        fhEMCALClusterCutsECenSignal[2]->Fill(energyOrMom,cen);
   }
  
  //Float_t pos[3];
  //clus->GetPosition(pos);
  //printf("Before Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
  
  //-----------------------------------------------------
  // Correct non linearity, smear energy or scale per SM
  //
  if ( fCorrectELinearity && GetCaloUtils()->IsCorrectionOfClusterEnergyOn() )
  {
    GetCaloUtils()->CorrectClusterEnergy(clus) ;
    
    AliDebug(5,Form("Correct Non Lin: Old E %3.2f, New E %3.2f",
                    fMomentum.E(),clus->E()));
  }
  
  // In case of MC analysis, to match resolution/calibration in real data
  // Not needed anymore, just leave for MC studies on systematics
  if( GetCaloUtils()->GetEMCALRecoUtils()->IsClusterEnergySmeared() )
  {
    Float_t rdmEnergy = GetCaloUtils()->GetEMCALRecoUtils()->SmearClusterEnergy(clus);
    
    AliDebug(5,Form("Smear energy: Old E %3.2f, New E %3.2f",clus->E(),rdmEnergy));
    
    clus->SetE(rdmEnergy);
  }

  // In case of uncalibrated data, or non final non linearity for MC
  // or calibration of MC for SMs behind TRD, apply a global scale factor per SM  
  if ( fScaleEPerSM && iSupMod < 22 && iSupMod >=0)
  {
    Float_t scale = fScaleFactorPerSM[iSupMod];
    
    AliDebug(5,Form("Scale energy for SM %d: Old E %3.2f, scale factor %1.5f",iSupMod,clus->E(),scale));
    
    clus->SetE(clus->E()*scale);
  }  
  
  clus->GetMomentum(fMomentum, fVertex[vindex]);
  
  energyOrMom = clus->E();
  if ( fHistoPtDependent ) energyOrMom = fMomentum.Pt();
  
  fhEMCALClusterEtaPhi->Fill(fMomentum.Eta(),GetPhi(fMomentum.Phi()));
  
  // Check effect linearity correction, energy smearing
  if ( fScaleEPerSM ||  fCorrectELinearity )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [3]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[3]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [3]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[3]->Fill(energyOrMom,cen);
    }
  }
  
  // Check the event BC depending on EMCal clustr before final cuts
  Double_t tof = clus->GetTOF()*1e9;
  
  Int_t bc = TMath::Nint(tof/50) + 9;
  //printf("tof %2.2f, bc+5=%d\n",tof,bc);
  
  SetEMCalEventBC(bc);
  
  //--------------------------------------
  // Apply some kinematical/acceptance cuts
  //
  if ( fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E() ) 
  {
    AliDebug(2,Form("Cluster E out of range, %2.2f < %2.2f < %2.2f",fEMCALPtMin,clus->E(),fEMCALPtMax));
    return ;
  }
  
  // Select cluster fiducial region
  //
  Bool_t bEMCAL = kFALSE;
  Bool_t bDCAL  = kFALSE;
  if ( fCheckFidCut )
  {
    if ( fFillEMCAL && fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kEMCAL) ) bEMCAL = kTRUE ;
    if ( fFillDCAL  && fFiducialCut->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kDCAL ) ) bDCAL  = kTRUE ;
  }
  else
  {
    bEMCAL = kTRUE;
  }
  
  //---------------------------------------------------------------------
  // Mask all cells in collumns facing ALICE thick material if requested
  //
  if ( GetCaloUtils()->GetNMaskCellColumns() )
  {
    AliDebug(2,Form("Masked collumn: cluster E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                    fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));
    
    if ( GetCaloUtils()->MaskFrameCluster(iSupMod, ietaMax) )
    {
      AliDebug(2,"Mask cluster");
      return;
    }
  }
  
  // Check effect of energy and fiducial cuts  
  if ( bEMCAL || bDCAL ) 
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [4]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[4]->Fill(energyOrMom,cen);

    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [4]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[4]->Fill(energyOrMom,cen);
    }
    
    fhEMCALClusterEtaPhiFidCut->Fill(fMomentum.Eta(),GetPhi(fMomentum.Phi()));
  }
  else 
  {
    AliDebug(2,"Cluster not on EMCal or DCal selected region");
    return ;
  }
  
  //----------------------------------------------------
  // Apply N cells cut
  //
  if ( fDataType != AliCaloTrackReader::kMC )
  {
    Int_t nCells =  clus->GetNCells();
    if (  nCells <= fEMCALNCellsCut ) 
    {
      AliDebug(2,Form("Reject cluster with n cells %d < %d",nCells, fEMCALNCellsCut));
      return ;
    }
    
    // Energy dependent n cell cut
    if ( clus->E() > fEMCALNCellsCutEnDepEnMin )
    {
      Float_t nCellsCut = fEMCALNCellsCutEnDepConstant + fEMCALNCellsCutEnDepSlope*clus->E();
      if ( nCells <= nCellsCut )
      {
        AliDebug(2,Form("Reject cluster with n cells %d < %2.1f for E %2.1f > %2.1f",
                        nCells, nCellsCut, clus->E(), fEMCALNCellsCutEnDepEnMin));
        return ;
      }
    }
  }
  
  // Check effect of n cells cut
  if ( fEMCALNCellsCut > 0 )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [5]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[5]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [5]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[5]->Fill(energyOrMom,cen);
    }
  }
  
  //----------------------------------------------------
  // Apply distance to bad channel cut
  //
  Double_t distBad = clus->GetDistanceToBadChannel() ; //Distance to bad channel
  fhEMCALClusterDisToBadE->Fill(energyOrMom,distBad);
  
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  
  if ( distBad < fEMCALBadChMinDist ) 
  {
    AliDebug(2, Form("Cluster close to bad, dist %2.2f < %2.2f",distBad,fEMCALBadChMinDist));
    return  ;
  }
  
  // Check effect distance to bad channel cut
  if ( fEMCALBadChMinDist > 0 )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [6]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[6]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [6]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[6]->Fill(energyOrMom,cen);
    }
  }
  //------------------------------------------
  // Apply time cut, count EMCal BC before cut
  //
  SetEMCalEventBCcut(bc);

  // Shift time in case of no calibration with rough factor
  Double_t tofShift = tof;
  //if(tof > 400) tofShift-=615;
  fhEMCALClusterTimeE->Fill(energyOrMom,tofShift);
  
  if ( !IsInTimeWindow(tof,energyOrMom) )
  {
    fNPileUpClusters++ ;
    if ( fUseEMCALTimeCut ) 
    {
      AliDebug(2,Form("Out of time window E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f, time %e",
                      fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta(),tof));
      
      return ;
    }
  }
  else
    fNNonPileUpClusters++;
  
  // Check effect of time cut
  if ( fUseEMCALTimeCut )
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [7]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[7]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [7]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[7]->Fill(energyOrMom,cen);
    }
  }
  
  //----------------------------------------
  // Apply cut on number of cells in different T-Card
  // than highest energy cell. At high E it should be more than 0
  // if not, sign of exotic cluster 
  
  Int_t   nDiff = 0, nSame = 0;
  Float_t eDiff = 0, eSame = 0;
  GetCaloUtils()->GetEnergyAndNumberOfCellsInTCard(clus, absIdMax, GetEMCALCells(), 
                                                   nDiff, nSame, eDiff, eSame,
                                                   fEMCALMinCellEnNdiffCut); 
  
  if ( nDiff == 0 && clus->E() > fEMCALHighEnergyNdiffCut )
  {
    AliInfo(Form("** Reader: Reject cluster with E = %2.1f (min %2.1f) and n cells in diff TCard = %d, for Ecell min = %1.2f; m02 %2.2f, ncells %d",
           clus->E(),fEMCALHighEnergyNdiffCut,nDiff,fEMCALMinCellEnNdiffCut,clus->GetM02(),clus->GetNCells()));
    return;
  }
  
  if ( fEMCALHighEnergyNdiffCut >= 40 &&  fEMCALHighEnergyNdiffCut <= 200)
  {
    if ( !fHistoCentDependent ) fhEMCALClusterCutsE   [8]->Fill(energyOrMom);
    else                        fhEMCALClusterCutsECen[8]->Fill(energyOrMom,cen);
    
    if ( fillEmbedSignalCluster )
    {
      if ( !fHistoCentDependent ) fhEMCALClusterCutsESignal   [8]->Fill(energyOrMom);
      else                        fhEMCALClusterCutsECenSignal[8]->Fill(energyOrMom,cen);
    }
  }
  
  //----------------------------------------------------
  // Smear the SS to try to match data and simulations,
  // do it only for simulations.
  //
  if ( fSmearShowerShape  && clus->GetNCells() > 2 )
  {
    Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(clus, GetEMCALCells()); 
//  Int_t nMaxima = clus->GetNExMax(); // For Run2

    if (  nMaxima >= fSmearNLMMin && nMaxima <= fSmearNLMMax )
    {
      AliDebug(2,Form("Smear shower shape - Original: %2.4f", clus->GetM02()));
      if(fSmearingFunction == kSmearingLandau)
      {
        clus->SetM02( clus->GetM02() + fRandom.Landau(0, fSmearShowerShapeWidth) );
      }
      else if ( fSmearingFunction == kSmearingLandauShift )
      {
        if(iclus%3 == 0 && clus->GetM02() > 0.1) clus->SetM02( clus->GetM02() + fRandom.Landau(0.05, fSmearShowerShapeWidth) );     //fSmearShowerShapeWidth = 0.035
      }
      else if (fSmearingFunction == kNoSmearing)
      {
        clus->SetM02( clus->GetM02() );
      }
      //clus->SetM02( fRandom.Landau(clus->GetM02(), fSmearShowerShapeWidth) );
      AliDebug(2,Form("Width %2.4f         Smeared : %2.4f", fSmearShowerShapeWidth,clus->GetM02()));
    }
  }
  
  //--------------------------------------------------------
  // Fill the corresponding array with the selected clusters
  // Usually just filling EMCal array with upper or lower clusters is enough, 
  // but maybe we want to do EMCal-DCal correlations.
  
  //if((fDebug > 2 && fMomentum.E() > 0.1) || fDebug > 10)
  AliDebug(2,Form("Selected clusters (EMCAL%d, DCAL%d), E %3.2f, pt %3.2f, phi %3.2f deg, eta %3.2f",
                  bEMCAL,bDCAL,fMomentum.E(),fMomentum.Pt(),RadToDeg(GetPhi(fMomentum.Phi())),fMomentum.Eta()));

  
  if      ( bEMCAL ) fEMCALClusters->Add(clus);
  else if ( bDCAL  ) fDCALClusters ->Add(clus);
  
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
  if ( fEMCALClustersListName == "" )
  {
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = fInputEvent->GetCaloCluster(iclus);
      
      if ( clus && clus->IsEMCAL() )
      {
        FillInputEMCALSelectCluster(clus, iclus);
      }//EMCAL cluster exists
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
    
    if ( !clusterList) 
    {
      AliWarning(Form("Wrong name of list with clusters?  <%s>",fEMCALClustersListName.Data()));
      
      if ( fInputEvent )
      {
        Int_t nInput = fInputEvent->GetList()->GetEntries();
        printf("\t N branches input event %d\n",nInput);
        for(Int_t i = 0; i < nInput; i++) 
          printf("\t %d %s\n",i,fInputEvent->GetList()->At(i)->GetName());
      }
      
      if ( fOutputEvent )
      {
        Int_t nOutput = fOutputEvent->GetList()->GetEntries();
        printf("\t N branches output event %d\n",nOutput);
        for(Int_t j = 0; j < nOutput; j++) 
          printf("\t %d %s\n",j,fOutputEvent->GetList()->At(j)->GetName());
      }
      return;
    }
    
    Int_t nclusters = clusterList->GetEntriesFast();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
    {
      AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
      //printf("E %f\n",clus->E());
      if (clus) FillInputEMCALSelectCluster(clus, iclus);
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
            GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absIdMax), fInputEvent->GetBunchCrossNumber(), tof, fCurrentParIndex);
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
  
  AliDebug(1,Form("EMCal selected clusters %d", 
                  fEMCALClusters->GetEntriesFast()));
  AliDebug(2,Form("\t n pile-up clusters %d, n non pile-up %d", 
                  fNPileUpClusters,fNNonPileUpClusters));
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
  
  AliDebug(1,Form("PHOS selected clusters %d",fPHOSClusters->GetEntriesFast())) ;  
}

//____________________________________________
/// Connects the array with EMCAL cells and the pointer.
//____________________________________________
void AliCaloTrackReader::FillInputEMCALCells()
{  
  if(fEMCALCellsListName.Length() == 0)
    fEMCALCells = fInputEvent->GetEMCALCells();
  else
    fEMCALCells = (AliVCaloCells*) fInputEvent->FindListObject(fEMCALCellsListName);
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
    if(GetDebug()>3) fNonStandardJets->Print("");
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
  
  fBackgroundJets = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fInputBackgroundJetBranchName.Data()));
    
  if(!fBackgroundJets)
  {
    //check if jet branch exist; exit if not
    fInputEvent->Print();

    AliFatal(Form("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fInputBackgroundJetBranchName.Data()));
  }
  else
  {
    AliDebug(1,"FillInputBackgroundJets");
    if(GetDebug()>3) fBackgroundJets->Print("");
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
	GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absIdMax), fInputEvent->GetBunchCrossNumber(), tof, fCurrentParIndex);
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
    if     (runNumber < 146861) fTriggerL0EventThreshold = 3. ;  // LHC11a
    else if(runNumber < 154000) fTriggerL0EventThreshold = 4. ;  // LHC11b,c
    else if(runNumber < 165000) fTriggerL0EventThreshold = 5.5;  // LHC11c,d,e
    else if(runNumber < 194000) fTriggerL0EventThreshold = 2  ;  // LHC12
    else if(runNumber < 197400) fTriggerL0EventThreshold = 3  ;  // LHC13def 
    else if(runNumber < 197400) fTriggerL0EventThreshold = 2  ;  // LHC13g 
    else if(runNumber < 244300) fTriggerL0EventThreshold = 5  ;  // LHC15 in, phys 1, 5 in phys2 
    else if(runNumber < 266400) fTriggerL0EventThreshold = 2.5;  // LHC16ir 
    else                        fTriggerL0EventThreshold = 3.5;  // LHC16s 
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
  if ( fEMCALNCellsCutEnDepEnMin < 200 )
   {
     printf("\t For Emin>%2.1f, ncell>%2.2f+ %2.2fE; \n",
              fEMCALNCellsCutEnDepEnMin,fEMCALNCellsCutEnDepConstant,fEMCALNCellsCutEnDepSlope) ;
   }
  printf("PHOS  N cells  > %d \n"        , fPHOSNCellsCut) ;
  printf("EMCAL Reject cluster N cells in diff = 0, for E>%2.2f and E cell > %1.2f \n", fEMCALHighEnergyNdiffCut,fEMCALMinCellEnNdiffCut) ;
  printf("EMCAL Time Cut: %3.1f < TOF  < %3.1f\n", fEMCALTimeCutMin, fEMCALTimeCutMax);
  printf("Use CTS         =     %d\n",     fFillCTS) ;
  printf("Use EMCAL       =     %d\n",     fFillEMCAL) ;
  printf("Use DCAL        =     %d\n",     fFillDCAL)  ;
  printf("Use PHOS        =     %d\n",     fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n",     fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n",     fFillPHOSCells) ;
  printf("Track status    =     %d\n", (Int_t) fTrackStatus) ;
  printf("Track SPD hit %d; ITS cluster >= %d; ITS chi2 < %2.1f; TPC cluster >= %d; TPC chi2 < %2.1f\n",
         fSelectSPDHitTracks,
         fSelectMinITSclusters, fSelectMaxChi2PerITScluster,
         fSelectMinTPCclusters, fSelectMaxChi2PerTPCcluster) ;
  printf("Track Mult Eta Cut =  %2.2f\n",  fTrackMultEtaCut) ;

  printf("Track Mult Pt Cuts:") ;
  for(Int_t icut = 0; icut < fTrackMultNPtCut; icut++) printf(" %2.2f GeV;",fTrackMultPtCut[icut]);
  printf("    \n") ;
 
  printf("Write delta AOD =     %d\n",     fWriteOutputDeltaAOD) ;
  printf("Recalculate Clusters = %d, E linearity = %d\n",    fRecalculateClusters, fCorrectELinearity) ;
  
  printf("Use Triggers selected in SE base class %d; If not what Trigger Mask? %d; MB Trigger Mask for mixed %d; \n",
         fEventTriggerAtSE, fEventTriggerMask,fMixEventTriggerMask);

  printf("Reject L1-G1 with L1-G2 %d; n bits accepted %d, n bits rejected %d; reject centrality trigger outliers %d \n",
         fRejectEMCalTriggerEventsL1HighWithL1Low,fAcceptEventsWithBit.GetSize(),
         fRejectEventsWithBit.GetSize(),fRemoveCentralityTriggerOutliers);
  
  if ( fComparePtHardAndJetPt )
    printf("Compare jet pt and pt hard to accept event, factor = %2.2f\n",fPtHardAndJetPtFactor);
  
  if ( fComparePtHardAndClusterPt )
    printf("Compare cluster pt and pt hard to accept event, factor = %2.2f\n",fPtHardAndClusterPtFactor);
  
  if ( fComparePtHardAndPromptPhotonPt )
    printf("Compare prompt photon pt and pt hard to accept event, factor = %2.2f\n",fPtHardAndPromptPhotonPtFactor);
  
  if ( fRemoveLEDEvents > 0 )
  {
    printf("Remove LED events %d, %2.1f < Ecell < %1.2f:\n", fRemoveLEDEvents, fLEDMinCellEnergy, fLEDMaxCellEnergy  );
    printf("\t SM - nCell >= %d - Sum E >= %2.0f; \n", fLEDHighNCellsCutSM, fLEDHighEnergyCutSM);
    printf("\t SM3: nCell <= %d - Sum E <= %2.0f  \n", fLEDLowNCellsCutSM3, fLEDLowEnergyCutSM3);
    
    if ( fRemoveLEDStripEvents > 0 )
    {
      printf("Remove LED strip? %d, with n strip > %d\n: "
               "\t Full SM, nCell > %d, Sum E > %2.0f;\n "
               "\t  1/3 SM, nCell > %d, Sum E > %2.0f;\n "
               "\t     SM3, nCell < %d, Sum E < %2.0f\n",
               fRemoveLEDStripEvents    , fLEDEventMaxNumberOfStrips, 
               fLEDHighNCellsCutStrip[0], fLEDHighEnergyCutStrip[0], 
               fLEDHighNCellsCutStrip[1], fLEDHighEnergyCutStrip[1], 
               fLEDLowNCellsCutSM3Strip , fLEDLowEnergyCutSM3Strip);
    }
  }
  
  printf("Delta AOD File Name =     %s\n", fDeltaAODFileName.Data()) ;
  printf("Centrality: Class %s, Option %d, Bin [%d,%d] \n",
         fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1]) ;

  printf("Accept only prompt photon clusters %d; Reject fragmentation photon clusters %d\n",
         fAcceptMCPromptPhotonOnly,fRejectMCFragmentationPhoton);

  printf("Accept clusters from N=%d generators\n",fNMCGenerToAccept);
  for(Int_t igen = 0; igen <= fNMCGenerToAccept; igen++ )
  printf("\t igen %d %s, index %d\n",
         igen, fMCGenerToAccept[igen].Data(),fMCGenerIndexToAccept[igen]);
  printf("Accept event header %s, Check Pythia event header %d\n",
         fMCGenerEventHeaderToAccept.Data(),fCheckPythiaEventHeader);

  //printf("    \n") ;
}

//__________________________________________
/// LED Events in period LHC11a contaminated 
/// EMCAL clusters sample, simple method
/// to reject such events. 
/// For period LHC11a only SM3 and sometimes SM4 gave problems, fRemoveLEDEvents=1 handles it
/// For testing, a generalization for all SMs is introduced for fRemoveLEDEvents>1
//__________________________________________
Bool_t  AliCaloTrackReader::RejectLEDEvents()
{
  // For LHC11a
  // Count number of cells with energy larger than 0.1 in SM3, cut on this number
  Int_t   ncellsSM[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Float_t ecellsSM[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  // LHC11a case
  if ( fRemoveLEDEvents == 1 )
  {
    for(Int_t icell = 0; icell < fInputEvent->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      Int_t absID = fInputEvent->GetEMCALCells()->GetCellNumber(icell);
      Int_t sm    = GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absID);
      
      if ( fInputEvent->GetEMCALCells()->GetAmplitude(icell) > 0.1 && sm == 3) ncellsSM[3]++;
    }
    
    Int_t ncellcut = 21;
    if ( GetFiredTriggerClasses().Contains("EMC") ) ncellcut = 35;
    
    if ( ncellsSM[3] >= ncellcut )
    {
      AliDebug(1,Form("Reject event with ncells in SM3 %d, cut %d, trig %s",
                      ncellsSM[3],ncellcut,GetFiredTriggerClasses().Data()));
      return kTRUE;
    }
  }
  // Run2 and general case
  else if ( fRemoveLEDEvents > 1 )
  {
    for(Int_t icell = 0; icell < fInputEvent->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      Int_t absID = fInputEvent->GetEMCALCells()->GetCellNumber(icell);
      Int_t sm    = GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absID);
      Float_t amp = fInputEvent->GetEMCALCells()->GetAmplitude(icell);
      
      if ( amp >= fLEDMinCellEnergy && amp <= fLEDMaxCellEnergy ) 
      {
        ncellsSM[sm]++;
        ecellsSM[sm]+=amp;
      }
    }
    
    for(Int_t ism = 0; ism < 20; ism++)
      fhEMCALNSumEnCellsPerSM->Fill(ncellsSM[ism],ecellsSM[ism]);
    
    // Run2
    if ( fRemoveLEDEvents == 2 ) 
    {
      // if there is some activity in SM3, accept the event
      if ( ncellsSM[3] <= fLEDLowNCellsCutSM3 || ecellsSM[3] <= fLEDLowEnergyCutSM3 ) 
      {      
        for(Int_t ism = 0; ism < 20; ism++)
        {
          if ( ism == 3 ) continue;
          
          if ( ncellsSM[ism] >=  fLEDHighNCellsCutSM )
          {
            printf("Reject event because of SM%d: ",ism);
            for(Int_t jsm = 0; jsm < 20; jsm++){
              if ( ncellsSM[jsm] > 0 ) printf("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]);}
            return kTRUE;
          }
          
          if ( ecellsSM[ism] >= fLEDHighEnergyCutSM )
          {
            printf("Reject event because of SM%d: ",ism);
            for(Int_t jsm = 0; jsm < 20; jsm++) {
              if ( ncellsSM[jsm] > 0 ) printf("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]);}
            return kTRUE;
          }
        }
      } // SM3 activity
      
    } // fRemoveLEDEvents == 2
    
    // General case without condition on SM3 low activity
    else //fRemoveLEDEvents > 2
    {      
      for(Int_t ism = 0; ism <  GetCaloUtils()->GetEMCALGeometry()->GetNumberOfSuperModules(); ism++)
      {
        if ( ncellsSM[ism] >=  fLEDHighNCellsCutSM )
        {
          printf("Reject event because of SM%d: ",ism);
          for(Int_t jsm = 0; jsm < 20; jsm++){
            if ( ncellsSM[jsm] > 0 ) printf("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]);}
          return kTRUE;
        }
        
        if ( ecellsSM[ism] >= fLEDHighEnergyCutSM )
        {
          printf("Reject event because of SM%d: ",ism);
          for(Int_t jsm = 0; jsm < 20; jsm++) {
            if ( ncellsSM[jsm] > 0 ) printf("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]);}
          return kTRUE;
        }
      } // SM loop
    }  //fRemoveLEDEvents > 2
    
    for(Int_t ism = 0; ism < 20; ism++)
      fhEMCALNSumEnCellsPerSMAfter->Fill(ncellsSM[ism],ecellsSM[ism]);
      
  } // fRemoveLEDEvents > 1
  
  // Check activity inside all the strips 
  // (24 strips, 2x24 cells in full SM, 2x8 cellsin 1/3 SM), 
  // n cells (max 48) and sum of cells energy
  if ( fRemoveLEDStripEvents )
  {    
    Float_t amp1   = 0., amp2   = 0. ;
    Int_t   absId1 = -1, absId2 = -1 ;
    Float_t enCellsStrip[20][24];
    Int_t    nCellsStrip[20][24];
    
    for (Int_t ism = 0; ism < 20; ism++)
    {      
      for (Int_t ieta = 0; ieta < 48; ieta=ieta+2)
      {
        enCellsStrip[ism][ieta/2] = 0.; 
        nCellsStrip [ism][ieta/2] = 0 ; 
        
        for (Int_t iphi = 0; iphi < 24; iphi++)
        {
          absId1 = GetCaloUtils()->GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(ism, iphi, ieta);
          if ( absId1 < 0 || absId1 > 17664 ) continue;
          
          absId2 = GetCaloUtils()->GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(ism, iphi, ieta+1);
          if ( absId2 < 0 || absId2 > 17664 ) continue;   
          
          amp1 = fInputEvent->GetEMCALCells()->GetCellAmplitude(absId1);
          amp2 = fInputEvent->GetEMCALCells()->GetCellAmplitude(absId2);
          
          if ( amp1 > fLEDMinCellEnergy && amp1 < fLEDMaxCellEnergy )
          {            
             nCellsStrip[ism][ieta/2]++;
            enCellsStrip[ism][ieta/2]+=amp1;
//            printf("Reader: cell1 %d amp %f; SM %d, strip %d, n cells %d, sum E %f \n",
//                   absId1, amp1,
//                   ism,ieta,
//                    nCellsStrip[ism][ieta/2],
//                   enCellsStrip[ism][ieta/2] 
//                   );
          }
          
          if ( amp2 > fLEDMinCellEnergy && amp2 < fLEDMaxCellEnergy )
          {            
             nCellsStrip[ism][ieta/2]++;
            enCellsStrip[ism][ieta/2]+=amp2;
//            printf("Reader: cell2 %d amp %f; SM %d, strip %d, n cells %d, sum E %f \n",
//                   absId2, amp2,
//                   ism,ieta,
//                    nCellsStrip[ism][ieta/2],
//                   enCellsStrip[ism][ieta/2] 
//                   );
          }
        }// iphi
        
        fhEMCALNSumEnCellsPerStrip->Fill(nCellsStrip[ism][ieta/2],enCellsStrip[ism][ieta/2]);
        
      } // ieta
    } // ism 
    
    // Count per event over event cut
    // Low activity on SM3 for emin = 0.5
    Bool_t bSM3StripsLowActivity = kTRUE;
    for (Int_t ieta = 0; ieta < 24; ieta++)
    {
      if ( enCellsStrip[3][ieta] > fLEDLowEnergyCutSM3Strip || 
            nCellsStrip[3][ieta] > fLEDLowNCellsCutSM3Strip   ) 
        bSM3StripsLowActivity = kFALSE;
    }
    
    // Count number or active strips, depending on cuts
    //
    Int_t   nStrips = 0;
    
    if ( bSM3StripsLowActivity )
    {
      Int_t   maxNCells = 24;
      Float_t maxECells = 80;
      for (Int_t ism = 0; ism < 20; ism++)
      {
        if ( ism == 3 ) continue ;

        maxNCells = fLEDHighNCellsCutStrip[0];
        maxECells = fLEDHighEnergyCutStrip[0];
        if ( ism == 10 || ism == 11 || 
             ism == 18 || ism == 19   ) 
        {
          maxNCells = fLEDHighNCellsCutStrip[1];
          maxECells = fLEDHighEnergyCutStrip[1];
        }
  
        for (Int_t ieta = 0; ieta < 24; ieta++)
        {
//          if ( nCellsStrip[ism][ieta] > 0 || enCellsStrip[ism][ieta] > 0 )
//            printf("Reader: SM %d, strip %d, n cells %d < %d , sum E %f < %f \n",
//                   ism,ieta,
//                   nCellsStrip[ism][ieta],maxNCells,
//                   enCellsStrip[ism][ieta],maxECells 
//                   );

          if( enCellsStrip[ism][ieta] >= maxECells || 
               nCellsStrip[ism][ieta] >= maxNCells   )
          {
            nStrips++;
          }
        } // ieta
      } // ism
    } // bSM03
    
    //printf("Reader: Number of strips %d\n",nStrips);

    if ( nStrips > fLEDEventMaxNumberOfStrips ) return kTRUE;
    
    for (Int_t ism = 0; ism < 20; ism++)
    {
      if ( fRemoveLEDEvents > 1 )
        fhEMCALNSumEnCellsPerSMAfterStripCut->Fill(ncellsSM[ism],ecellsSM[ism]);
      
      for (Int_t ieta = 0; ieta < 24; ieta++)
      {
        fhEMCALNSumEnCellsPerStripAfter->Fill(nCellsStrip[ism][ieta],enCellsStrip[ism][ieta]);
      }
    }
    
  } // remove strip LED events
  
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
  if(fBackgroundJets)  fBackgroundJets  -> Clear("C");
  //fBackgroundJets->Reset();
}

//___________________________________________
/// Tag event depending on trigger name.
/// Set also the L1 bit defining the EGA or EJE triggers.
/// depending on the trigger class version, if not set by user.
//___________________________________________
void AliCaloTrackReader::SetEventTriggerBit(UInt_t mask)
{	
  fEventTrigMinBias       = kFALSE;
  fEventTrigCentral       = kFALSE;
  fEventTrigSemiCentral   = kFALSE;
  
  fEventTrigEMCALL0       = kFALSE;
  fEventTrigEMCALL1Gamma1 = kFALSE;
  fEventTrigEMCALL1Gamma2 = kFALSE;
  fEventTrigEMCALL1Jet1   = kFALSE;
  fEventTrigEMCALL1Jet2   = kFALSE;

  fEventTrigDCALL0        = kFALSE;
  fEventTrigDCALL1Gamma1  = kFALSE;
  fEventTrigDCALL1Gamma2  = kFALSE;
  fEventTrigDCALL1Jet1    = kFALSE;
  fEventTrigDCALL1Jet2    = kFALSE;
  
  fEventTrigMinBiasCaloOnly       = kFALSE;
  fEventTrigEMCALL0CaloOnly       = kFALSE;
  fEventTrigEMCALL1Gamma1CaloOnly = kFALSE;
  fEventTrigEMCALL1Gamma2CaloOnly = kFALSE;
  fEventTrigEMCALL1Jet1CaloOnly   = kFALSE;
  fEventTrigEMCALL1Jet2CaloOnly   = kFALSE;
  
  fEventTrigDCALL0CaloOnly        = kFALSE;
  fEventTrigDCALL1Gamma1CaloOnly  = kFALSE;
  fEventTrigDCALL1Gamma2CaloOnly  = kFALSE;
  fEventTrigDCALL1Jet1CaloOnly    = kFALSE;
  fEventTrigDCALL1Jet2CaloOnly    = kFALSE;
  
  AliDebug(1,Form("Select trigger mask bit %d - Trigger Event %s - Select <%s>",
                  fEventTriggerMask,GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data()));
  
  if ( fEventTriggerMask <=0 )// in case no mask set
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
	  // EMC/DMC L1 Gamma
	  if     ( mask & AliVEvent::kEMCEGA )
    {
      //printf("EGA trigger bit\n");
      // EMCal
      if     (GetFiredTriggerClasses().Contains("EG"))
      {
        if     (GetFiredTriggerClasses().Contains("EGA")) fEventTrigEMCALL1Gamma1 = kTRUE;
        else
        {
          if(GetFiredTriggerClasses().Contains("EG1")) fEventTrigEMCALL1Gamma1 = kTRUE;
          if(GetFiredTriggerClasses().Contains("EG2")) fEventTrigEMCALL1Gamma2 = kTRUE;
        }
      }
      
      // DCal
      if     (GetFiredTriggerClasses().Contains("DG"))
      {
        if(GetFiredTriggerClasses().Contains("DG1")) fEventTrigDCALL1Gamma1 = kTRUE;
        if(GetFiredTriggerClasses().Contains("DG2")) fEventTrigDCALL1Gamma2 = kTRUE;
      }
    }
    
	  // EMC L1 Jet
    if ( mask & AliVEvent::kEMCEJE )
    {
      //printf("EGA trigger bit\n");
      // EMCal
      if     (GetFiredTriggerClasses().Contains("EJ"))
      {
        if     (GetFiredTriggerClasses().Contains("EJE")) fEventTrigEMCALL1Jet1 = kTRUE;
        else
        {
          if(GetFiredTriggerClasses().Contains("EJ1")) fEventTrigEMCALL1Jet1 = kTRUE;
          if(GetFiredTriggerClasses().Contains("EJ2")) fEventTrigEMCALL1Jet2 = kTRUE;
        }
      }
      
      // DCal
      if     (GetFiredTriggerClasses().Contains("DJ"))
      {
        if(GetFiredTriggerClasses().Contains("DJ1")) fEventTrigDCALL1Jet1 = kTRUE;
        if(GetFiredTriggerClasses().Contains("DJ2")) fEventTrigDCALL1Jet2 = kTRUE;
      }
    }
		
    // EMC L0
    if( ( mask & AliVEvent::kEMC7 ) ||
        ( mask & AliVEvent::kEMC1 )       )
    {
      //printf("L0 trigger bit\n");
      if      ( GetFiredTriggerClasses().Contains("EMC") ) fEventTrigEMCALL0 = kTRUE;
      else if ( GetFiredTriggerClasses().Contains("DMC") ) fEventTrigDCALL0  = kTRUE;
    }
    
    //------------
    // kCaloOnly
    if ( mask & AliVEvent::kCaloOnly )
    {
      // EMC/DMC L1 Gamma
      if ( GetFiredTriggerClasses().Contains("EG") )
      {
        if ( GetFiredTriggerClasses().Contains("EG1") ) fEventTrigEMCALL1Gamma1CaloOnly = kTRUE;
        if ( GetFiredTriggerClasses().Contains("EG2") ) fEventTrigEMCALL1Gamma2CaloOnly = kTRUE;
      }
      
      // DCal L1 Gamma
      if ( GetFiredTriggerClasses().Contains("DG") )
      {
        if ( GetFiredTriggerClasses().Contains("DG1") ) fEventTrigDCALL1Gamma1CaloOnly = kTRUE;
        if ( GetFiredTriggerClasses().Contains("DG2") ) fEventTrigDCALL1Gamma2CaloOnly = kTRUE;
      }
      
      // EMC L1 Jet
      if ( GetFiredTriggerClasses().Contains("EJ") )
      {
        if ( GetFiredTriggerClasses().Contains("EJ1") ) fEventTrigEMCALL1Jet1CaloOnly = kTRUE;
        if ( GetFiredTriggerClasses().Contains("EJ2") ) fEventTrigEMCALL1Jet2CaloOnly = kTRUE;
      }
      
      // DCal L1 Jet
      if ( GetFiredTriggerClasses().Contains("DJ") )
      {
        if ( GetFiredTriggerClasses().Contains("DJ1") ) fEventTrigDCALL1Jet1CaloOnly = kTRUE;
        if ( GetFiredTriggerClasses().Contains("DJ2") ) fEventTrigDCALL1Jet2CaloOnly = kTRUE;
      }
      
      if ( GetFiredTriggerClasses().Contains("CDMC7PER") )
      {
        fEventTrigDCALL0CaloOnly = kTRUE;
      }
      
      if ( GetFiredTriggerClasses().Contains("CINT7-B-NOPF-CALOPLUS") )
      {
        fEventTrigMinBiasCaloOnly = kTRUE;
      }
    }
    //------------
	  
    // Min Bias Pb-Pb
    if ( mask & AliVEvent::kCentral )
    {
      //printf("MB central trigger bit\n");
	    fEventTrigCentral = kTRUE;
    }
	  
    // Min Bias Pb-Pb
    if ( mask & AliVEvent::kSemiCentral )
    {
      //printf("MB semi central trigger bit\n");
	    fEventTrigSemiCentral = kTRUE;
    }
	  
    // Min Bias pp, PbPb, pPb
    if ( (mask & AliVEvent::kMB  ) ||
         (mask & AliVEvent::kINT7)    )
    {
      //printf("MB trigger bit\n");
	    fEventTrigMinBias = kTRUE;
    }
	}
  
  AliDebug(1,Form("Event bits: MB       %d, Cen    %d, Sem    %d, CaloMB %d\n"
                  "            L0 EMC   %d, L1-EG1 %d, L1-EG2 %d, L1-EJ1 %d, L1-EJ2 %d,\n"
                  "            L0 DMC   %d, L1-DG1 %d, L1-DG2 %d, L1-DJ1 %d, L1-DJ2 %d,\n"
                  "kCaloOnly:  L0 EMC   %d, L1-EG1 %d, L1-EG2 %d, L1-EJ1 %d, L1-EJ2 %d,\n"
                  "            L0 DMC   %d, L1-DG1 %d, L1-DG2 %d, L1-DJ1 %d, L1-DJ2 %d;\n",
                  fEventTrigMinBias, fEventTrigCentral      , fEventTrigSemiCentral  , fEventTrigMinBiasCaloOnly,
                  fEventTrigEMCALL0, fEventTrigEMCALL1Gamma1, fEventTrigEMCALL1Gamma2, fEventTrigEMCALL1Jet1    , fEventTrigEMCALL1Jet2,
                  fEventTrigDCALL0 , fEventTrigDCALL1Gamma1 , fEventTrigDCALL1Gamma2 , fEventTrigDCALL1Jet1     , fEventTrigDCALL1Jet2 ,
                  fEventTrigEMCALL0CaloOnly, fEventTrigEMCALL1Gamma1CaloOnly, fEventTrigEMCALL1Gamma2CaloOnly, fEventTrigEMCALL1Jet1CaloOnly, fEventTrigEMCALL1Jet2CaloOnly,
                  fEventTrigDCALL0CaloOnly , fEventTrigDCALL1Gamma1CaloOnly , fEventTrigDCALL1Gamma2CaloOnly , fEventTrigDCALL1Jet1CaloOnly , fEventTrigDCALL1Jet2CaloOnly  )  );
  
  
  // L1 trigger bit
  if ( fBitEGA == 0 && fBitEJE == 0 )
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
/// Called in ESD/AOD readers
//____________________________________________________________
void AliCaloTrackReader::SetInputEvent(AliVEvent * input)
{
  fInputEvent  = input;
  
  fMixedEvent = dynamic_cast<AliMixedEvent*>(GetInputEvent()) ;
  if ( fMixedEvent )
    fNMixedEvent = fMixedEvent->GetNumberOfEvents() ;
  
  // Delete previous vertex
  if ( fVertex )
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
  
  // Recover embedded event instead
  if ( fEmbeddedEvent[1] ) 
    fInputEvent = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent();
}

//____________________________________________________________
/// Set the MC event used in analysis.
/// In case of embedding, recover it from external event
//____________________________________________________________
void AliCaloTrackReader::SetMC(AliMCEvent * mc)              
{ 
  fMC = mc ; // MCEvent(); Set in the main steering task.
  
  // In case of embedding get it from external event
  //
  if ( fEmbeddedEvent[0] )
  {
    fMC = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalMCEvent();
    if ( !fMC ) 
    {
      AliWarning("Embedded MC event not found\n");
    }
  } // embedded
}
