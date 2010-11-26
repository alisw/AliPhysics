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

//---- ANALYSIS system ----
#include "AliCaloTrackReader.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
#include "AliESDtrack.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDtrackCuts.h"

ClassImp(AliCaloTrackReader)
  
  
//____________________________________________________________________________
  AliCaloTrackReader::AliCaloTrackReader() : 
    TObject(), fEventNumber(-1), fCurrentFileName(""),fDataType(0), fDebug(0), 
    fFiducialCut(0x0), fCheckFidCut(kFALSE), fComparePtHardAndJetPt(kFALSE), fPtHardAndJetPtFactor(7),
    fCTSPtMin(0), fEMCALPtMin(0),fPHOSPtMin(0), fAODBranchList(new TList ),
    fAODCTS(new TObjArray()), fAODEMCAL(new TObjArray()), fAODPHOS(new TObjArray()),
    fEMCALCells(0x0), fPHOSCells(0x0),
    fInputEvent(0x0), fOutputEvent(0x0),fMC(0x0),
    fFillCTS(0),fFillEMCAL(0),fFillPHOS(0),
    fFillEMCALCells(0),fFillPHOSCells(0), 
//    fSecondInputAODTree(0x0), fSecondInputAODEvent(0x0),
//    fSecondInputFileName(""),fSecondInputFirstEvent(0), 
//    fAODCTSNormalInputEntries(0), fAODEMCALNormalInputEntries(0), 
//    fAODPHOSNormalInputEntries(0), 
    fTrackStatus(0),   fESDtrackCuts(0), fTrackMult(0), fTrackMultEtaCut(0.9),
    fReadStack(kFALSE), fReadAODMCParticles(kFALSE), 
    fDeltaAODFileName("deltaAODPartCorr.root"),fFiredTriggerClassName(""),
    fAnaLED(kFALSE),fTaskName(""),fCaloUtils(0x0), 
    fMixedEvent(NULL), fNMixedEvent(1), fVertex(NULL), 
    fWriteOutputDeltaAOD(kFALSE),fOldAOD(kFALSE){
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
  
  if(fAODCTS){
    fAODCTS->Clear() ; 
    delete fAODCTS ;
  }
  
  if(fAODEMCAL){
    fAODEMCAL->Clear() ; 
    delete fAODEMCAL ;
  }
  
  if(fAODPHOS){
    fAODPHOS->Clear() ; 
    delete fAODPHOS ;
  }
  
  if(fEMCALCells){
    fEMCALCells->Clear() ; 
    delete fEMCALCells ;
  }
  
  if(fPHOSCells){
    fPHOSCells->Clear() ; 
    delete fPHOSCells ;
  }

  if(fVertex){
    for (Int_t i = 0; i < fNMixedEvent; i++) {
      delete [] fVertex[i] ;

    }
    delete [] fVertex ;
	}

  if(fESDtrackCuts)   delete fESDtrackCuts;
  
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
  printf("    \n") ;
} 

//___________________________________________________
Bool_t AliCaloTrackReader::FillInputEvent(const Int_t iEntry, const char * currentFileName) {
  //Fill the event counter and input lists that are needed, called by the analysis maker.

  fEventNumber = iEntry;
  fCurrentFileName = TString(currentFileName);
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
  
  //Fill the arrays with cluster/tracks/cells data
   if(fFillEMCALCells) 
    FillInputEMCALCells();
  if(fFillPHOSCells)  
    FillInputPHOSCells();
	
  if(fFillCTS)   
    FillInputCTS();
  if(fFillEMCAL) 
    FillInputEMCAL();
  if(fFillPHOS)  
    FillInputPHOS();

	
  return kTRUE ;
}

//__________________________________________________
void AliCaloTrackReader::ResetLists() {
  //  Reset lists, called by the analysis maker 

  if(fAODCTS)     fAODCTS     -> Clear();
  if(fAODEMCAL)   fAODEMCAL   -> Clear();
  if(fAODPHOS)    fAODPHOS    -> Clear();
  if(fEMCALCells) fEMCALCells -> Clear();
  if(fPHOSCells)  fPHOSCells  -> Clear();
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
    
    if(fDataType!=kMC)fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex[0]); 
    else { 
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
      
    if(fDebug > 1)
      printf("AliCaloTrackReader::FillVertexArray() - Single Event Vertex : %f,%f,%f\n",fVertex[0][0],fVertex[0][1],fVertex[0][2]);

  } else { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++) {
      fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
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
      
      fAODCTS->Add(track);        
       
    }//Pt and Fiducial cut passed. 
  }// track loop
	
  //fAODCTSNormalInputEntries = fAODCTS->GetEntriesFast();
  if(fDebug > 1) 
    printf("AliCaloTrackReader::FillInputCTS()   - aod entries %d, input tracks %d, pass status %d, multipliticy %d\n", fAODCTS->GetEntriesFast(), nTracks, nstatus, fTrackMult);//fAODCTSNormalInputEntries);
  
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
    //			  fAODCTS->Add(track);
    //			  
    //		  }//Pt and Fiducial cut passed. 
    //	  }// track loop
    //	  
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputCTS()   - aod normal entries %d, after second input %d\n", fAODCTSNormalInputEntries, fAODCTS->GetEntriesFast());
    //  }	//second input loop
    //	
}

  //____________________________________________________________________________
void AliCaloTrackReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputEMCAL()\n");
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
  for (Int_t iclus =  0; iclus <  nclusters; iclus++) {
    AliVCluster * clus = 0;
    if ( (clus = fInputEvent->GetCaloCluster(iclus)) ) {
      if (IsEMCALCluster(clus)){
        
        //Check if the cluster contains any bad channel and if close to calorimeter borders
        Int_t vindex = 0 ;  
        if (fMixedEvent) 
          vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
        
        if(GetCaloUtils()->ClusterContainsBadChannel("EMCAL",clus->GetCellsAbsId(), clus->GetNCells())) 
          continue;
        if(!GetCaloUtils()->CheckCellFiducialRegion(clus, (AliVCaloCells*)fInputEvent->GetEMCALCells(), fInputEvent, vindex)) 
          continue;
        
        TLorentzVector momentum ;
        
        clus->GetMomentum(momentum, fVertex[vindex]);      
        
        if(fEMCALPtMin < momentum.Pt()){
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) 
            continue;
          
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
          
          //Correct non linearity
          if(GetCaloUtils()->IsCorrectionOfClusterEnergyOn()){
            GetCaloUtils()->CorrectClusterEnergy(clus) ;
            //printf("Linearity Corrected Energy %f\n",clus->E());  
          }
          
          //Recalculate cluster position
          if(GetCaloUtils()->IsRecalculationOfClusterPositionOn()){
            GetCaloUtils()->RecalculateClusterPosition(GetEMCALCells(),clus); 
            //clus->GetPosition(pos);
            //printf("After  Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
          }
          
          if (fMixedEvent) {
            clus->SetID(iclus) ; 
          }
            
          fAODEMCAL->Add(clus);	
          
        }//Pt and Fiducial cut passed.
      }//EMCAL cluster
    }// cluster exists
  }// cluster loop
  
  //Recalculate track matching
  GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent);
  
  //fAODEMCALNormalInputEntries = fAODEMCAL->GetEntriesFast();
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod entries %d\n",  fAODEMCAL->GetEntriesFast());//fAODEMCALNormalInputEntries);
  
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
    //					  fAODEMCAL->Add(clus);	
    //					}//Pt and Fiducial cut passed.
    //				}//EMCAL cluster
    //			}// cluster exists
    //		}// cluster loop
    //		
    //	  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod normal entries %d, after second input %d\n", fAODEMCALNormalInputEntries, fAODEMCAL->GetEntriesFast());
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
          
          fAODPHOS->Add(clus);	
          
        }//Pt and Fiducial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  //fAODPHOSNormalInputEntries = fAODPHOS->GetEntriesFast() ;
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod entries %d\n",  fAODPHOS->GetEntriesFast());//fAODPHOSNormalInputEntries);
  
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
    //						fAODPHOS->Add(clus);	
    //					}//Pt and Fiducial cut passed.
    //				}//PHOS cluster
    //			}// cluster exists
    //		}// cluster loop
    //		if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod normal entries %d, after second input %d\n", fAODPHOSNormalInputEntries, fAODPHOS->GetEntriesFast());
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

