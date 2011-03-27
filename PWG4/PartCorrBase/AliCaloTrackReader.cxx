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
/* $Id:  $ */

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors (CTS).
// Not all MC particles/tracks/clusters are kept, some kinematical/fiducial restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader: Fills Kinematics data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS) 
//              
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()

//-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TFile.h"
#include "TRandom3.h"

//---- ANALYSIS system ----
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliTriggerAnalysis.h"
#include "AliESDVZERO.h"

//---- PartCorr/EMCAL ---
#include "AliEMCALRecoUtils.h"
#include "AliCaloTrackReader.h"

ClassImp(AliCaloTrackReader)
  
  
//____________________________________________________________________________
  AliCaloTrackReader::AliCaloTrackReader() : 
    TObject(), fEventNumber(-1), //fCurrentFileName(""),
    fDataType(0), fDebug(0), 
    fFiducialCut(0x0), fCheckFidCut(kFALSE), fComparePtHardAndJetPt(kFALSE), fPtHardAndJetPtFactor(7),
    fCTSPtMin(0), fEMCALPtMin(0),fPHOSPtMin(0), fAODBranchList(new TList ),
    fCTSTracks(new TObjArray()), fEMCALClusters(new TObjArray()), fPHOSClusters(new TObjArray()),
    fEMCALCells(0x0), fPHOSCells(0x0),
    fInputEvent(0x0), fOutputEvent(0x0),fMC(0x0),
    fFillCTS(0),fFillEMCAL(0),fFillPHOS(0),
    fFillEMCALCells(0),fFillPHOSCells(0),  fSelectEmbeddedClusters(kFALSE),
    fRemoveSuspiciousClusters(kFALSE), fSmearClusterEnergy(kFALSE), fRandom(),
//    fSecondInputAODTree(0x0), fSecondInputAODEvent(0x0),
//    fSecondInputFileName(""),fSecondInputFirstEvent(0), 
//    fCTSTracksNormalInputEntries(0), fEMCALClustersNormalInputEntries(0), 
//    fPHOSClustersNormalInputEntries(0), 
    fTrackStatus(0),   fESDtrackCuts(0), fTrackMult(0), fTrackMultEtaCut(0.8),
    fReadStack(kFALSE), fReadAODMCParticles(kFALSE), 
    fDeltaAODFileName("deltaAODPartCorr.root"),fFiredTriggerClassName(""),
    fAnaLED(kFALSE),fTaskName(""),fCaloUtils(0x0), 
    fMixedEvent(NULL), fNMixedEvent(1), fVertex(NULL), 
    fWriteOutputDeltaAOD(kFALSE),fOldAOD(kFALSE),fCaloFilterPatch(kFALSE),
    fEMCALClustersListName(""),fZvtxCut(0.), 
    fDoEventSelection(kFALSE),   fDoV0ANDEventSelection(kFALSE),
    fTriggerAnalysis (new AliTriggerAnalysis), fCentralityClass("V0M"),fCentralityOpt(10)
   
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//_________________________________
AliCaloTrackReader::~AliCaloTrackReader() {
  //Dtor
  
  if(fFiducialCut) delete fFiducialCut ;
	
  if(fAODBranchList){
    fAODBranchList->Delete();
    delete fAODBranchList ;
  }  
  
  if(fCTSTracks){
    if(fDataType!=kMC)fCTSTracks->Clear() ; 
    else              fCTSTracks->Delete() ; 
    delete fCTSTracks ;
  }
  
  if(fEMCALClusters){
    if(fDataType!=kMC)fEMCALClusters->Clear("C") ; 
    else              fEMCALClusters->Delete() ; 
    delete fEMCALClusters ;
  }
  
  if(fPHOSClusters){
    if(fDataType!=kMC)fPHOSClusters->Clear("C") ; 
    else              fPHOSClusters->Delete() ; 
    delete fPHOSClusters ;
  }
  
//  if(fEMCALCells){
//    delete fEMCALCells ;
//  }
//  
//  if(fPHOSCells){
//    delete fPHOSCells ;
//  }

  if(fVertex){
    for (Int_t i = 0; i < fNMixedEvent; i++) {
      delete [] fVertex[i] ;

    }
    delete [] fVertex ;
	}

  if(fESDtrackCuts)    delete fESDtrackCuts;
  if(fTriggerAnalysis) delete fTriggerAnalysis;
  
//  Pointers not owned, done by the analysis frame
//  if(fInputEvent)  delete fInputEvent ;
//  if(fOutputEvent) delete fOutputEvent ;
//  if(fMC)          delete fMC ;  
	
//  if(fSecondInputAODTree){
//    fSecondInputAODTree->Clear();
//    delete fSecondInputAODTree;
//  }
//	
//  if(fSecondInputAODEvent) delete fSecondInputAODEvent ;
	
  //  Pointer not owned, deleted by maker
  //if (fCaloUtils) delete fCaloUtils ;

}

//_________________________________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndJetPt(){
  // Check the event, if the requested ptHard is much larger than the jet pT, then there is a problem.
  // Only for PYTHIA.
  if(!fReadStack) return kTRUE; //Information not filtered to AOD
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader")){
    TParticle * jet =  0;
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Int_t nTriggerJets =  pygeh->NTriggerJets();
    Float_t ptHard = pygeh->GetPtHard();
    
    //if(fDebug > 1) printf("AliMCAnalysisUtils::PythiaEventHeader: Njets: %d, pT Hard %f\n",nTriggerJets, ptHard);
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
      pygeh->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      //Compare jet pT and pt Hard
      //if(fDebug > 1) printf("AliMCAnalysisUtils:: %d pycell jet pT %f\n",ijet, jet->Pt());
      if(jet->Pt() > fPtHardAndJetPtFactor * ptHard) {
        printf("AliMCAnalysisUtils::PythiaEventHeader: Njets: %d, pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n",
               nTriggerJets, ptHard, jet->Pt(), fPtHardAndJetPtFactor);
        return kFALSE;
      }
    }
    if(jet) delete jet; 
  }
  
  return kTRUE ;
  
}

//____________________________________________________________________________
AliStack* AliCaloTrackReader::GetStack() const {
  //Return pointer to stack
  if(fMC)
    return fMC->Stack();
  else{
    if(fDebug > 1) printf("AliCaloTrackReader::GetStack() - Stack is not available\n"); 
    return 0x0 ;
  }
}

//____________________________________________________________________________
AliHeader* AliCaloTrackReader::GetHeader() const {
  //Return pointer to header
  if(fMC)
    return fMC->Header();
  else{
    printf("AliCaloTrackReader::Header is not available\n"); 
    return 0x0 ;
  }
}
//____________________________________________________________________________
AliGenEventHeader* AliCaloTrackReader::GetGenEventHeader() const {
  //Return pointer to Generated event header
  if(fMC)
    return fMC->GenEventHeader();
  else{
    printf("AliCaloTrackReader::GenEventHeader is not available\n"); 
    return 0x0 ;
  }
}

//____________________________________________________________________________
TClonesArray* AliCaloTrackReader::GetAODMCParticles(Int_t input) const {
  //Return list of particles in AOD. Do it for the corresponding input event.
  
  TClonesArray * rv = NULL ; 
  if(fDataType == kAOD){
 
    if(input == 0){
      //Normal input AOD
      AliAODEvent * evt = dynamic_cast<AliAODEvent*> (fInputEvent) ;
      if(evt)
        rv = (TClonesArray*)evt->FindListObject("mcparticles");
      else  
        printf("AliCaloTrackReader::GetAODMCParticles() - wrong AOD input index? %d, or non existing tree? \n",input); 
      
    } //else if(input == 1 && fSecondInputAODEvent){ //Second input AOD   
//      rv = (TClonesArray*) fSecondInputAODEvent->FindListObject("mcparticles");	
//    } 
    
  } else {
      printf("AliCaloTrackReader::GetAODMCParticles() - Input are not AODs\n"); 
  }
  
  return rv ; 
}

//____________________________________________________________________________
AliAODMCHeader* AliCaloTrackReader::GetAODMCHeader(Int_t input) const {
  //Return MC header in AOD. Do it for the corresponding input event.
  AliAODMCHeader *mch = NULL;
  if(fDataType == kAOD){
    //Normal input AOD
    if(input == 0) {
      mch = (AliAODMCHeader*)((AliAODEvent*)fInputEvent)->FindListObject("mcheader");
    }
    //		//Second input AOD
    //		else if(input == 1){ 
    //       mch = (AliAODMCHeader*) fSecondInputAODEvent->FindListObject("mcheader");
    //  }
    else {
      printf("AliCaloTrackReader::GetAODMCHeader() - wrong AOD input index, %d\n",input);
    }
  }
  else {
    printf("AliCaloTrackReader::GetAODMCHeader() - Input are not AODs\n");
  }
  
  return mch;
}

//_______________________________________________________________
void AliCaloTrackReader::Init()
{
  //Init reader. Method to be called in AliAnaPartCorrMaker
  
  //Get the file with second input events if the filename is given
  //Get the tree and connect the AODEvent. Only with AODs
  
  if(fReadStack && fReadAODMCParticles){
    printf("AliCaloTrackReader::Init() - Cannot access stack and mcparticles at the same time, change them \n");
    fReadStack = kFALSE;
    fReadAODMCParticles = kFALSE;
  }
  
//	if(fSecondInputFileName!=""){
//		if(fDataType == kAOD){
//			TFile * input2   = new TFile(fSecondInputFileName,"read");
//			printf("AliCaloTrackReader::Init() - Second input file opened: %s, size %d \n", input2->GetName(), (Int_t) input2->GetSize());
//			fSecondInputAODTree = (TTree*) input2->Get("aodTree");
//			if(fSecondInputAODTree) printf("AliCaloTrackReader::Init() - Second input tree opened: %s, entries %d \n", 
//										   fSecondInputAODTree->GetName(), (Int_t) fSecondInputAODTree->GetEntries());
//			else{
//			 printf("AliCaloTrackReader::Init() - Second input tree not available, STOP \n");
//			 abort();
//			}
//			fSecondInputAODEvent = new AliAODEvent;
//			fSecondInputAODEvent->ReadFromTree(fSecondInputAODTree);
//			if(fSecondInputFirstEvent >= fSecondInputAODTree->GetEntriesFast()){
//				printf("AliCaloTrackReader::Init() - Requested first event of second input %d, is larger than number of events %d, STOP\n", 
//					   fSecondInputFirstEvent, (Int_t) fSecondInputAODTree->GetEntriesFast());
//				abort();
//			}
//		}
//		else printf("AliCaloTrackReader::Init() - Second input not added, reader is not AOD\n");
//	}
}

//_______________________________________________________________
void AliCaloTrackReader::InitParameters()
{
  //Initialize the parameters of the analysis.
  fDataType   = kESD ;
  fCTSPtMin   = 0.2 ;
  fEMCALPtMin = 0.2 ;
  fPHOSPtMin  = 0.2 ;

  //Do not filter the detectors input by default.
  fFillEMCAL      = kFALSE;
  fFillPHOS       = kFALSE;
  fFillCTS        = kFALSE;
  fFillEMCALCells = kFALSE;
  fFillPHOSCells  = kFALSE;

  //fSecondInputFileName   = "" ;
  //fSecondInputFirstEvent = 0 ;
  fReadStack             = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fReadAODMCParticles    = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fDeltaAODFileName      = "deltaAODPartCorr.root";
  fFiredTriggerClassName      = "";
	 	
  fAnaLED = kFALSE;
  
  //We want tracks fitted in the detectors:
  //fTrackStatus=AliESDtrack::kTPCrefit;
  //fTrackStatus|=AliESDtrack::kITSrefit;  
  
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

  fV0ADC[0] = 0;   fV0ADC[1] = 0; 
  fV0Mul[0] = 0;   fV0Mul[1] = 0; 

  fZvtxCut   = 10.;
  
  //Centrality
  fCentralityBin[0]=fCentralityBin[1]=-1;
  
  //Cluster smearing
  fSmearClusterEnergy   = kFALSE;
  fSmearClusterParam[0] = 0.07; // * sqrt E term
  fSmearClusterParam[1] = 0.02; // * E term
  fSmearClusterParam[2] = 0.00; // constant
  
}

//________________________________________________________________
void AliCaloTrackReader::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  printf("Task name      : %s\n", fTaskName.Data()) ;
  printf("Data type      : %d\n", fDataType) ;
  printf("CTS Min pT     : %2.1f GeV/c\n", fCTSPtMin) ;
  printf("EMCAL Min pT   : %2.1f GeV/c\n", fEMCALPtMin) ;
  printf("PHOS Min pT    : %2.1f GeV/c\n", fPHOSPtMin) ;
  printf("Use CTS         =     %d\n", fFillCTS) ;
  printf("Use EMCAL       =     %d\n", fFillEMCAL) ;
  printf("Use PHOS        =     %d\n", fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n", fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n", fFillPHOSCells) ;
  printf("Track status    =     %d\n", (Int_t) fTrackStatus) ;
  printf("Track Mult Eta Cut =  %d\n", (Int_t) fTrackMultEtaCut) ;
  printf("Write delta AOD =     %d\n", fWriteOutputDeltaAOD) ;

  if(fComparePtHardAndJetPt)
	  printf("Compare jet pt and pt hard to accept event, factor = %2.2f",fPtHardAndJetPtFactor);
	
//  if(fSecondInputFileName!="") {
//	  printf("Second Input File Name     =     %s\n", fSecondInputFileName.Data()) ;
//	  printf("Second Input First Event   =     %d\n", fSecondInputFirstEvent) ;
//  }
	
  printf("Read Kine from, stack? %d, AOD ? %d \n", fReadStack, fReadAODMCParticles) ;
  printf("Delta AOD File Name =     %s\n", fDeltaAODFileName.Data()) ;
  printf("Centrality: Class %s, Option %d, Bin [%d,%d] \n", fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1]) ;

  printf("    \n") ;
  
} 

//___________________________________________________
Bool_t AliCaloTrackReader::FillInputEvent(const Int_t iEntry, const char * /*currentFileName*/) {
  //Fill the event counter and input lists that are needed, called by the analysis maker.

  fEventNumber = iEntry;
  //fCurrentFileName = TString(currentFileName);
  if(!fInputEvent) {
	  if(fDebug >= 0) printf("AliCaloTrackReader::FillInputEvent() - Input event not available, skip event analysis\n");
	  return kFALSE;
  }
  //Select events only fired by a certain trigger configuration if it is provided
  Int_t eventType = 0;
  if(fInputEvent->GetHeader())
	  eventType = ((AliVHeader*)fInputEvent->GetHeader())->GetEventType();
  if( fFiredTriggerClassName  !="" && !fAnaLED){
    if(eventType!=7)
      return kFALSE; //Only physics event, do not use for simulated events!!!
    if(fDebug > 0) 
      printf("AliCaloTrackReader::FillInputEvent() - FiredTriggerClass <%s>, selected class <%s>, compare name %d\n",
	     GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data(), GetFiredTriggerClasses().Contains(fFiredTriggerClassName));
    if( !GetFiredTriggerClasses().Contains(fFiredTriggerClassName) ) return kFALSE;
  }
  else if(fAnaLED){
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
	 
	  if(eventType!=7 && fDebug > 1 )printf("AliCaloTrackReader::FillInputEvent() - DO LED, Event Type <%d>, 8 Calibration \n",  eventType);
	  if(eventType!=8)return kFALSE;
  }
		
  //In case of analysis of events with jets, skip those with jet pt > 5 pt hard	
  if(fComparePtHardAndJetPt && GetStack()) {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
  }

  //In case of mixing events with other AOD file	
 // if(fDataType == kAOD && fSecondInputAODTree){
//	 
//    if(fDebug > 1) 
//      printf("AliCaloTrackReader::FillInputEvent() - Get event %d from second input AOD file \n", iEntry+fSecondInputFirstEvent);
//    if(fSecondInputAODTree->GetEntriesFast() <= iEntry+fSecondInputFirstEvent) {
//      if(fSecondInputAODTree->GetEntriesFast() == iEntry+fSecondInputFirstEvent) 
//			 printf("AliCaloTrackReader::FillInputEvent() - Skip events from event %d, no more events in second AOD file \n", iEntry);
//      return kFALSE;
//    }
//    
//    //Get the Event
//    Int_t nbytes = fSecondInputAODTree->GetEvent(iEntry+fSecondInputFirstEvent);
//    if ( nbytes == 0 ) {//If nothing in AOD
//      printf("AliCaloTrackReader::FillInputEvent() - Nothing in Second AOD input, STOP\n");
//      abort() ; 
//    }
//    
//  }
	
  //Fill Vertex array
  FillVertexArray();
  //Reject events with Z vertex too large, only for SE analysis, if not, cut on the analysis code
  if(!GetMixedEvent() && TMath::Abs(fVertex[0][2]) > fZvtxCut) return kFALSE;  
  
  //------------------------------------------------------
  //Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------
  if(fDoEventSelection){
    if(!fCaloFilterPatch){
      //Do not analyze events with pileup
      Bool_t bPileup = fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.); //Default values, if not it does not compile
      //Bool_t bPileup = event->IsPileupFromSPD(); 
      if(bPileup) return kFALSE;
      
      if(fDoV0ANDEventSelection){
        Bool_t bV0AND = kTRUE; 
        AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
        if(esd) 
          bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0AND);
        //else bV0AND = //FIXME FOR AODs
        if(!bV0AND) return kFALSE;
      }
      
      if(!CheckForPrimaryVertex()) return kFALSE;
    }//CaloFilter patch
    else{ 
      if(fInputEvent->GetNumberOfCaloClusters() > 0) {
        AliVCluster * calo = fInputEvent->GetCaloCluster(0);
        if(calo->GetNLabels() == 4){
          Int_t * selection = calo->GetLabels();
          Bool_t bPileup = selection[0];
          if(bPileup) return kFALSE;
          
          Bool_t bGoodV = selection[1]; 
          if(!bGoodV) return kFALSE;
          
          if(fDoV0ANDEventSelection){
            Bool_t bV0AND = selection[2]; 
            if(!bV0AND) return kFALSE;
          }
          
          fTrackMult = selection[3];
          if(fTrackMult == 0) return kFALSE;
        } else {
          //First filtered AODs, track multiplicity stored there.  
          fTrackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
          if(fTrackMult == 0) return kFALSE;          
        }
      }//at least one cluster
      else {
        //printf("AliCaloTrackReader::FillInputEvent() - No clusters in event\n");
        //Remove events with  vertex (0,0,0), bad vertex reconstruction
        if(TMath::Abs(fVertex[0][0]) < 1.e-6 && TMath::Abs(fVertex[0][1]) < 1.e-6 && TMath::Abs(fVertex[0][2]) < 1.e-6) return kFALSE;
        
        //First filtered AODs, track multiplicity stored there.  
        fTrackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
        if(fTrackMult == 0) return kFALSE;
      }// no cluster
    }// CaloFileter patch
  }// Event selection
  //------------------------------------------------------

  //Check if there is a centrality value, PbPb analysis, and if a centrality bin selection is requested
  //If we need a centrality bin, we select only those events in the corresponding bin.
  if(GetCentrality() && fCentralityBin[0]>=0 && fCentralityBin[1]>=0 && fCentralityOpt==100){
    Int_t cen = GetEventCentrality();
    if(cen > fCentralityBin[1] || cen < fCentralityBin[0]) return kFALSE; //reject events out of bin.
  }
  
  //Fill the arrays with cluster/tracks/cells data
   if(fFillEMCALCells) 
    FillInputEMCALCells();
  if(fFillPHOSCells)  
    FillInputPHOSCells();
	
  if(fFillCTS){   
    FillInputCTS();
    //Accept events with at least one track
    if(fTrackMult == 0 && fDoEventSelection) return kFALSE;
  }
  
  if(fFillEMCAL) 
    FillInputEMCAL();
  if(fFillPHOS)  
    FillInputPHOS();

  FillInputVZERO();
	
  return kTRUE ;
}

//__________________________________________________
void AliCaloTrackReader::ResetLists() {
  //  Reset lists, called by the analysis maker 

  if(fCTSTracks)       fCTSTracks     -> Clear();
  if(fEMCALClusters)   fEMCALClusters   -> Clear("C");
  if(fPHOSClusters)    fPHOSClusters    -> Clear("C");
//  if(fEMCALCells) fEMCALCells -> Clear("");
//  if(fPHOSCells)  fPHOSCells  -> Clear("");

  fV0ADC[0] = 0;   fV0ADC[1] = 0; 
  fV0Mul[0] = 0;   fV0Mul[1] = 0; 

}

//____________________________________________________________________________
void AliCaloTrackReader::SetInputEvent(AliVEvent* const input)  
{
  fInputEvent  = input;
  fMixedEvent = dynamic_cast<AliMixedEvent*>(GetInputEvent()) ; 
  if (fMixedEvent) {
    fNMixedEvent = fMixedEvent->GetNumberOfEvents() ; 
  }

  //Delete previous vertex
  if(fVertex){
    for (Int_t i = 0; i < fNMixedEvent; i++) {
      delete [] fVertex[i] ; 
    }
    delete [] fVertex ;
  }
  
  fVertex = new Double_t*[fNMixedEvent] ; 
  for (Int_t i = 0; i < fNMixedEvent; i++) {
    fVertex[i] = new Double_t[3] ; 
    fVertex[i][0] = 0.0 ; 
    fVertex[i][1] = 0.0 ; 
    fVertex[i][2] = 0.0 ; 
  }
}

//__________________________________________________
Int_t AliCaloTrackReader::GetEventCentrality() const {
  //Return current event centrality
  
  if(GetCentrality()){
    if(fCentralityOpt==100)     return (Int_t) GetCentrality()->GetCentralityPercentile(fCentralityClass);
    else if(fCentralityOpt==10) return GetCentrality()->GetCentralityClass10(fCentralityClass); 
    else if(fCentralityOpt==5)  return GetCentrality()->GetCentralityClass5(fCentralityClass);
    else {
      printf("AliAnaPartCorrBaseClass::Unknown centrality option %d, use 5, 10 or 100\n",fCentralityOpt);
      return 0;
    } 
  }
  else return 0;
  
}

//____________________________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3]) const {
  //Return vertex position to be used for single event analysis
  vertex[0]=fVertex[0][0];  
  vertex[1]=fVertex[0][1];  
  vertex[2]=fVertex[0][2];
}

//____________________________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3], const Int_t evtIndex) const {
  //Return vertex position for mixed event, recover the vertex in a particular event.
  
  //Int_t evtIndex = 0; // for single events only one vertex stored in position 0, default value
  //if (fMixedEvent && clusterID >=0) {
  //  evtIndex=GetMixedEvent()->EventIndexForCaloCluster(clusterID) ; 
  //}
  
  vertex[0]=fVertex[evtIndex][0];  vertex[1]=fVertex[evtIndex][1];  vertex[2]=fVertex[evtIndex][2];
  
}
//


//____________________________________________________________________________
void AliCaloTrackReader::FillVertexArray() {
  
  //Fill data member with vertex
  //In case of Mixed event, multiple vertices
  
  //Delete previous vertex
  if(fVertex){
    for (Int_t i = 0; i < fNMixedEvent; i++) {
      delete [] fVertex[i] ; 
    }
    delete [] fVertex ;  
  }
  
  fVertex = new Double_t*[fNMixedEvent] ; 
  for (Int_t i = 0; i < fNMixedEvent; i++) {
    fVertex[i] = new Double_t[3] ; 
    fVertex[i][0] = 0.0 ; 
    fVertex[i][1] = 0.0 ; 
    fVertex[i][2] = 0.0 ; 
  }          
  
  if (!fMixedEvent) { //Single event analysis
    if(fDataType!=kMC){

      if(fInputEvent->GetPrimaryVertex()){
        fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex[0]); 
      }
      else {
        printf("AliCaloTrackReader::FillVertexArray() - NULL primary vertex\n");
        fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
      }//Primary vertex pointer do not exist
      
    } else {//MC read event 
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
      
    if(fDebug > 1)
      printf("AliCaloTrackReader::FillVertexArray() - Single Event Vertex : %f,%f,%f\n",fVertex[0][0],fVertex[0][1],fVertex[0][2]);

  } else { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++) {
      if (fMixedEvent->GetVertexOfEvent(iev))
        fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
      else { // no vertex found !!!!
        AliWarning("No vertex found");
      }

      if(fDebug > 1)
        printf("AliCaloTrackReader::FillVertexArray() - Multi Event %d Vertex : %f,%f,%f\n",iev,fVertex[iev][0],fVertex[iev][1],fVertex[iev][2]);

    }
  }
  
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputCTS() {
  //Return array with Central Tracking System (CTS) tracks
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS()\n");
  
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  Double_t p[3];
  fTrackMult = 0;
  Int_t nstatus = 0;
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliVTrack * track = (AliVTrack*)fInputEvent->GetTrack(itrack) ; // retrieve track from esd

    //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
    if (fTrackStatus && !((track->GetStatus() & fTrackStatus) == fTrackStatus)) 
      continue ;
    
    nstatus++;
    
    if(fDataType==kESD && !fESDtrackCuts->AcceptTrack((AliESDtrack*)track)) continue;
    
    // Track filter selection
    //if (fTrackFilter) {
	  //  selectInfo = fTrackFilter->IsSelected(esdTrack);
	  //  if (!selectInfo && !(esd->GetPrimaryVertex())->UsesTrack(esdTrack->GetID())) continue;
   // }
    
    //Count the tracks in eta < 0.9
    //printf("Eta %f cut  %f\n",TMath::Abs(track->Eta()),fTrackMultEtaCut);
    if(TMath::Abs(track->Eta())< fTrackMultEtaCut) fTrackMult++;
    
    track->GetPxPyPz(p) ;
    TLorentzVector momentum(p[0],p[1],p[2],0);
    
    if(fCTSPtMin < momentum.Pt()){
      
      if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"CTS")) 
        continue;
      
      if(fDebug > 2 && momentum.Pt() > 0.1) 
        printf("AliCaloTrackReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
               momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
      
      if (fMixedEvent) {
        track->SetID(itrack);
      }
      
      fCTSTracks->Add(track);        
       
    }//Pt and Fiducial cut passed. 
  }// track loop
	
  //fCTSTracksNormalInputEntries = fCTSTracks->GetEntriesFast();
  if(fDebug > 1) 
    printf("AliCaloTrackReader::FillInputCTS()   - aod entries %d, input tracks %d, pass status %d, multipliticy %d\n", fCTSTracks->GetEntriesFast(), nTracks, nstatus, fTrackMult);//fCTSTracksNormalInputEntries);
  
    //  //If second input event available, add the clusters.
    //  if(fSecondInputAODTree && fSecondInputAODEvent){
    //	  nTracks   = fSecondInputAODEvent->GetNumberOfTracks() ;
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputCTS()   - Add second input tracks, entries %d\n", nTracks) ;
    //	  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    //		  AliAODTrack * track = ((AliAODEvent*)fSecondInputAODEvent)->GetTrack(itrack) ; // retrieve track from esd
    //		  
    //		  //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
    //		  if (fTrackStatus && !((track->GetStatus() & fTrackStatus) == fTrackStatus)) continue ;
    //		  
    //		  track->GetPxPyPz(p) ;
    //		  TLorentzVector momentum(p[0],p[1],p[2],0);
    //		  
    //		  if(fCTSPtMin < momentum.Pt()){
    //
    //			  if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"CTS")) continue;
    //
    //			  if(fDebug > 2 && momentum.Pt() > 0.1) printf("AliCaloTrackReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
    //								       momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    //			  
    //			  fCTSTracks->Add(track);
    //			  
    //		  }//Pt and Fiducial cut passed. 
    //	  }// track loop
    //	  
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputCTS()   - aod normal entries %d, after second input %d\n", fCTSTracksNormalInputEntries, fCTSTracks->GetEntriesFast());
    //  }	//second input loop
    //	
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputEMCALAlgorithm(AliVCluster * clus, const Int_t iclus) {
  //Fill the EMCAL data in the array, do it
  
  Int_t vindex = 0 ;  
  if (fMixedEvent) 
    vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
  
  //Check if the cluster contains any bad channel and if close to calorimeter borders
  if(GetCaloUtils()->ClusterContainsBadChannel("EMCAL",clus->GetCellsAbsId(), clus->GetNCells())) 
    return;
  if(!GetCaloUtils()->CheckCellFiducialRegion(clus, (AliVCaloCells*)fInputEvent->GetEMCALCells(), fInputEvent, vindex)) 
    return;
  
  //Remove suspicious clusters
  if(fRemoveSuspiciousClusters){
    Int_t ncells      = clus->GetNCells();
    Float_t energy    = clus->E();
    Float_t minNCells = 1+energy/3;//-x*x*0.0033
    if(ncells < minNCells) {
      //if(energy > 2)printf("AliCaloTrackReader::FillInputEMCALAlgorithm() - Remove cluster: e %2.2f, Ncells %d, min Ncells %2.1f\n",energy,ncells,minNCells);
      return;
    }
//    else {
//      if(energy > 2)printf("AliCaloTrackReader::FillInputEMCALAlgorithm() - Keep cluster: e %2.2f, Ncells %d, min Ncells %2.1f\n",energy,ncells,minNCells);
//    }
  }//Suspicious
  
  if(fSelectEmbeddedClusters){
    if(clus->GetNLabels()==0) return;
    //else printf("Embedded cluster,  %d, n label %d label %d  \n",iclus,clus->GetNLabels(),clus->GetLabel());
  }
  
  TLorentzVector momentum ;
  
  clus->GetMomentum(momentum, fVertex[vindex]);      
  
  if(fEMCALPtMin < momentum.Pt()){
    
    if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) 
      return;
    
    if(fDebug > 2 && momentum.E() > 0.1) 
      printf("AliCaloTrackReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
             momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    
    //Float_t pos[3];
    //clus->GetPosition(pos);
    //printf("Before Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
    
    //Recalibrate the cluster energy 
    if(GetCaloUtils()->IsRecalibrationOn()) {
      Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, GetEMCALCells());
      clus->SetE(energy);
      //printf("Recalibrated Energy %f\n",clus->E());  
      GetCaloUtils()->RecalculateClusterShowerShapeParameters(GetEMCALCells(),clus);
      GetCaloUtils()->RecalculateClusterPID(clus);
      
    }
    
    //Recalculate distance to bad channels, if new list of bad channels provided
    GetCaloUtils()->RecalculateClusterDistanceToBadChannel(GetEMCALCells(),clus);
    
    //Recalculate cluster position
    if(GetCaloUtils()->IsRecalculationOfClusterPositionOn()){
      GetCaloUtils()->RecalculateClusterPosition(GetEMCALCells(),clus); 
      //clus->GetPosition(pos);
      //printf("After  Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
    }
    
    //Correct non linearity
    if(GetCaloUtils()->IsCorrectionOfClusterEnergyOn()){
      GetCaloUtils()->CorrectClusterEnergy(clus) ;
      //printf("Linearity Corrected Energy %f\n",clus->E());  
    }          
    
    //In case of MC analysis, to match resolution/calibration in real data
    if(fSmearClusterEnergy){
      Float_t energy    = clus->E();
      Float_t rdmEnergy = fRandom.Gaus(energy,fSmearClusterParam[0]*TMath::Sqrt(energy)+
                                       fSmearClusterParam[1]*energy+fSmearClusterParam[2]);
      clus->SetE(rdmEnergy);
      if(fDebug > 2) printf("\t Energy %f, smeared %f\n", energy, clus->E());
    }
    
    if (fMixedEvent) 
      clus->SetID(iclus) ; 
    
    fEMCALClusters->Add(clus);	
  }
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputEMCAL()\n");
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  if(fEMCALClustersListName==""){
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++) {
      AliVCluster * clus = 0;
      if ( (clus = fInputEvent->GetCaloCluster(iclus)) ) {
        if (IsEMCALCluster(clus)){          
          FillInputEMCALAlgorithm(clus, iclus);
        }//EMCAL cluster
      }// cluster exists
    }// cluster loop
    
    //Recalculate track matching
    if(fDataType==kESD)GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent);
    
  }//Get the clusters from the input event
  else {
    TClonesArray * clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
    if(!clusterList){
      printf("AliCaloTrackReader::FillInputEMCAL() - Wrong name of list with clusters? <%s>\n",fEMCALClustersListName.Data());
      return;
    }
    Int_t nclusters = clusterList->GetEntriesFast();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++) {
      AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
      //printf("E %f\n",clus->E());
      if (clus) FillInputEMCALAlgorithm(clus, iclus);
      else printf("AliCaloTrackReader::FillInputEMCAL() - Null cluster in list!\n");
      
    }// cluster loop
    
    //Recalculate track matching, not necessary, already done in the reclusterization task
    //GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent,clusterList);
    
  }
    
  //fEMCALClustersNormalInputEntries = fEMCALClusters->GetEntriesFast();
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod entries %d\n",  fEMCALClusters->GetEntriesFast());//fEMCALClustersNormalInputEntries);
  
    //If second input event available, add the clusters.
    //  if(fSecondInputAODTree && fSecondInputAODEvent){
    //	  GetSecondInputAODVertex(v);
    //	  nclusters = ((AliAODEvent*)fSecondInputAODEvent)->GetNumberOfCaloClusters();
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - Add second input clusters, entries %d\n", nclusters) ;
    //		for (Int_t iclus =  0; iclus < nclusters; iclus++) {
    //			AliAODCaloCluster * clus = 0;
    //			if ( (clus = ((AliAODEvent*)fSecondInputAODEvent)->GetCaloCluster(iclus)) ) {
    //				if (clus->IsEMCAL()){
    //					TLorentzVector momentum ;
    //					clus->GetMomentum(momentum, v);      
    //					
    //					if(fEMCALPtMin < momentum.Pt()){
    //
    //						if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) continue;
    //
    //						if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
    //																	momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    //					  fEMCALClusters->Add(clus);	
    //					}//Pt and Fiducial cut passed.
    //				}//EMCAL cluster
    //			}// cluster exists
    //		}// cluster loop
    //		
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod normal entries %d, after second input %d\n", fEMCALClustersNormalInputEntries, fEMCALClusters->GetEntriesFast());
    //
    //	} //second input loop
}

  //____________________________________________________________________________
void AliCaloTrackReader::FillInputPHOS() {
  //Return array with PHOS clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputPHOS()\n");
	  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
  for (Int_t iclus = 0; iclus < nclusters; iclus++) {
    AliVCluster * clus = 0;
    if ( (clus = fInputEvent->GetCaloCluster(iclus)) ) {
      if (IsPHOSCluster(clus)){
        //Check if the cluster contains any bad channel and if close to calorimeter borders
        Int_t vindex = 0 ;  
        if (fMixedEvent) 
          vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
        if( GetCaloUtils()->ClusterContainsBadChannel("PHOS",clus->GetCellsAbsId(), clus->GetNCells())) 
          continue;
        if(!GetCaloUtils()->CheckCellFiducialRegion(clus, fInputEvent->GetPHOSCells(), fInputEvent, vindex)) 
          continue;
        
        TLorentzVector momentum ;
        
        clus->GetMomentum(momentum, fVertex[vindex]);      
        
        if(fPHOSPtMin < momentum.Pt()){
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"PHOS")) 
            continue;
          
          if(fDebug > 2 && momentum.E() > 0.1) 
            printf("AliCaloTrackReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
                   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
          
            //Recalibrate the cluster energy 
          if(GetCaloUtils()->IsRecalibrationOn()) {
            Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliAODCaloCells*)GetPHOSCells());
            clus->SetE(energy);
          }
          
          if (fMixedEvent) {
            clus->SetID(iclus) ; 
          }              
          
          fPHOSClusters->Add(clus);	
          
        }//Pt and Fiducial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  //fPHOSClustersNormalInputEntries = fPHOSClusters->GetEntriesFast() ;
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod entries %d\n",  fPHOSClusters->GetEntriesFast());//fPHOSClustersNormalInputEntries);
  
    //If second input event available, add the clusters.
    //  if(fSecondInputAODTree && fSecondInputAODEvent){  
    //	  GetSecondInputAODVertex(v);
    //	  nclusters = ((AliAODEvent*)fSecondInputAODEvent)->GetNumberOfCaloClusters();
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - Add second input clusters, entries %d\n", nclusters);
    //		for (Int_t iclus =  0; iclus < nclusters; iclus++) {
    //			AliAODCaloCluster * clus = 0;
    //			if ( (clus = ((AliAODEvent*)fSecondInputAODEvent)->GetCaloCluster(iclus)) ) {
    //				if (clus->IsPHOS()){
    //					TLorentzVector momentum ;
    //					clus->GetMomentum(momentum, v);      
    //					
    //					if(fPHOSPtMin < momentum.Pt()){
    //
    //						if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"PHOS")) continue;
    //
    //						if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
    //																	momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    //						fPHOSClusters->Add(clus);	
    //					}//Pt and Fiducial cut passed.
    //				}//PHOS cluster
    //			}// cluster exists
    //		}// cluster loop
    //		if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod normal entries %d, after second input %d\n", fPHOSClustersNormalInputEntries, fPHOSClusters->GetEntriesFast());
    //  }	//second input loop
  
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputEMCALCells() {
    //Return array with EMCAL cells in aod format
  
  fEMCALCells = fInputEvent->GetEMCALCells(); 
  
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputPHOSCells() {
    //Return array with PHOS cells in aod format
  
  fPHOSCells = fInputEvent->GetPHOSCells(); 
  
}

//____________________________________________________________________________
void AliCaloTrackReader::FillInputVZERO(){
  //Fill VZERO information in data member, add all the channels information.
  AliVVZERO* v0 = fInputEvent->GetVZEROData();
  //printf("Init V0: ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  
  if (v0) 
  {
    AliESDVZERO* esdV0 = dynamic_cast<AliESDVZERO*> (v0);
    for (Int_t i = 0; i < 32; i++)
    {
      if(esdV0){//Only available in ESDs
        fV0ADC[0] += (Int_t)esdV0->GetAdcV0C(i);
        fV0ADC[1] += (Int_t)esdV0->GetAdcV0A(i);
      }
      fV0Mul[0] += (Int_t)v0->GetMultiplicityV0C(i);
      fV0Mul[1] += (Int_t)v0->GetMultiplicityV0A(i);
    }
    if(fDebug > 0)
      printf("V0: ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  }
  else
  {
    if(fDebug > 0)
      printf("Cannot retrieve V0 ESD! Run w/ null V0 charges\n ");
  }
}


//____________________________________________________________________________
Bool_t AliCaloTrackReader::IsEMCALCluster(AliVCluster* cluster) const {
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

//____________________________________________________________________________
Bool_t AliCaloTrackReader::IsPHOSCluster(AliVCluster * cluster) const {
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

//____________________________________________________________________________
Bool_t AliCaloTrackReader::CheckForPrimaryVertex(){
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //Only for ESDs ...
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (fInputEvent);
  if(!event) return kFALSE;
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0) {
    return kTRUE;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1) {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0) {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1) {
      //      cout<<"bad vertex type::"<< event->GetPrimaryVertex()->GetName() << endl;
      return kFALSE;
    }
  }
  
  return kFALSE;  
  
}



