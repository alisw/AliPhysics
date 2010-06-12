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
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS) 
//                
//-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TFile.h"

//---- ANALYSIS system ----
#include "AliCaloTrackReader.h"
#include "AliFiducialCut.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODEvent.h"

ClassImp(AliCaloTrackReader)
  
  
//____________________________________________________________________________
  AliCaloTrackReader::AliCaloTrackReader() : 
    TObject(), fEventNumber(-1), fCurrentFileName(""),fDataType(0), fDebug(0), 
    fFiducialCut(0x0), fComparePtHardAndJetPt(kFALSE), fPtHardAndJetPtFactor(7),
    fCTSPtMin(0), fEMCALPtMin(0),fPHOSPtMin(0),
    fAODCTS(new TObjArray()), fAODEMCAL(new TObjArray()), fAODPHOS(new TObjArray()),
    fEMCALCells(0x0), fPHOSCells(0x0),
    fInputEvent(0x0), fOutputEvent(0x0),fMC(0x0),
    fFillCTS(0),fFillEMCAL(0),fFillPHOS(0),
    fFillEMCALCells(0),fFillPHOSCells(0), 
    fSecondInputAODTree(0x0), fSecondInputAODEvent(0x0),
    fSecondInputFileName(""),fSecondInputFirstEvent(0), 
    fAODCTSNormalInputEntries(0), fAODEMCALNormalInputEntries(0), 
    fAODPHOSNormalInputEntries(0), fTrackStatus(0), 
    fReadStack(kFALSE), fReadAODMCParticles(kFALSE), 
    fCleanOutputStdAOD(kFALSE), fDeltaAODFileName("deltaAODPartCorr.root"),fFiredTriggerClassName(""),
    fAnaLED(kFALSE),fTaskName(""),fCaloUtils(0x0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliCaloTrackReader::AliCaloTrackReader(const AliCaloTrackReader & reader) :   
  TObject(reader), fEventNumber(reader.fEventNumber), fCurrentFileName(reader.fCurrentFileName), 
  fDataType(reader.fDataType), fDebug(reader.fDebug),
  fFiducialCut(reader.fFiducialCut),
  fComparePtHardAndJetPt(reader.fComparePtHardAndJetPt),
  fPtHardAndJetPtFactor(reader.fPtHardAndJetPtFactor),
  fCTSPtMin(reader.fCTSPtMin), fEMCALPtMin(reader.fEMCALPtMin),fPHOSPtMin(reader.fPHOSPtMin), 
  fAODCTS(new TObjArray(*reader.fAODCTS)),  
  fAODEMCAL(new TObjArray(*reader.fAODEMCAL)),
  fAODPHOS(new TObjArray(*reader.fAODPHOS)),
  fEMCALCells(new TNamed(*reader.fEMCALCells)),
  fPHOSCells(new TNamed(*reader.fPHOSCells)),
  fInputEvent(reader.fInputEvent), fOutputEvent(reader.fOutputEvent), fMC(reader.fMC),
  fFillCTS(reader.fFillCTS),fFillEMCAL(reader.fFillEMCAL),fFillPHOS(reader.fFillPHOS),
  fFillEMCALCells(reader.fFillEMCALCells),fFillPHOSCells(reader.fFillPHOSCells),
  fSecondInputAODTree(reader.fSecondInputAODTree), 
  fSecondInputAODEvent(reader.fSecondInputAODEvent),
  fSecondInputFileName(reader.fSecondInputFileName), 
  fSecondInputFirstEvent(reader.fSecondInputFirstEvent),
  fAODCTSNormalInputEntries(reader.fAODCTSNormalInputEntries), 
  fAODEMCALNormalInputEntries(reader.fAODEMCALNormalInputEntries), 
  fAODPHOSNormalInputEntries(reader.fAODPHOSNormalInputEntries),
  fTrackStatus(reader.fTrackStatus),
  fReadStack(reader.fReadStack), fReadAODMCParticles(reader.fReadAODMCParticles),
  fCleanOutputStdAOD(reader.fCleanOutputStdAOD), fDeltaAODFileName(reader.fDeltaAODFileName),
  fFiredTriggerClassName(reader.fFiredTriggerClassName),
  fAnaLED(reader.fAnaLED), 
  fTaskName(reader.fTaskName),
  fCaloUtils(new AliCalorimeterUtils(*reader.fCaloUtils))
{
  // cpy ctor  
}

//_________________________________________________________________________
//AliCaloTrackReader & AliCaloTrackReader::operator = (const AliCaloTrackReader & source)
//{
//  // assignment operator
//  
//  if(&source == this) return *this;
//  
//  fDataType    = source.fDataType ;
//  fDebug       = source.fDebug ;
//  fEventNumber = source.fEventNumber ;
//  fCurrentFileName = source.fCurrentFileName ;
//  fFiducialCut = source.fFiducialCut;
//	
//  fComparePtHardAndJetPt = source.fComparePtHardAndJetPt;
//  fPtHardAndJetPtFactor  = source.fPtHardAndJetPtFactor;
//	
//  fCTSPtMin    = source.fCTSPtMin ;
//  fEMCALPtMin  = source.fEMCALPtMin ;
//  fPHOSPtMin   = source.fPHOSPtMin ; 
//  
//  fAODCTS     = new TObjArray(*source.fAODCTS) ;
//  fAODEMCAL   = new TObjArray(*source.fAODEMCAL) ;
//  fAODPHOS    = new TObjArray(*source.fAODPHOS) ;
//  fEMCALCells = new TNamed(*source.fEMCALCells) ;
//  fPHOSCells  = new TNamed(*source.fPHOSCells) ;
//
//  fInputEvent  = source.fInputEvent;
//  fOutputEvent = source.fOutputEvent;
//  fMC          = source.fMC;
//  
//  fFillCTS        = source.fFillCTS;
//  fFillEMCAL      = source.fFillEMCAL;
//  fFillPHOS       = source.fFillPHOS;
//  fFillEMCALCells = source.fFillEMCALCells;
//  fFillPHOSCells  = source.fFillPHOSCells;
//
//  fSecondInputAODTree    = source.fSecondInputAODTree;
//  fSecondInputAODEvent   = source.fSecondInputAODEvent;
//  fSecondInputFileName   = source.fSecondInputFileName;
//  fSecondInputFirstEvent = source.fSecondInputFirstEvent;
//
//  fAODCTSNormalInputEntries   = source.fAODCTSNormalInputEntries; 
//  fAODEMCALNormalInputEntries = source.fAODEMCALNormalInputEntries; 
//  fAODPHOSNormalInputEntries  = source.fAODPHOSNormalInputEntries;
//	
//  fTrackStatus        = source.fTrackStatus;
//  fReadStack          = source.fReadStack;
//  fReadAODMCParticles = source.fReadAODMCParticles;	
//	
//  fCleanOutputStdAOD  = source.fCleanOutputStdAOD;
//  fDeltaAODFileName   = source.fDeltaAODFileName;
//	
//  fFiredTriggerClassName = source.fFiredTriggerClassName  ;
//	
//  return *this;
//  
//}

//_________________________________
AliCaloTrackReader::~AliCaloTrackReader() {
  //Dtor
  
  if(fFiducialCut) delete fFiducialCut ;
	
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
		TParticle * jet =  new TParticle;
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
	if(fDataType == kAOD){
	 //Normal input AOD
	 if(input == 0) return (TClonesArray*)((AliAODEvent*)fInputEvent)->FindListObject("mcparticles");
	  //Second input AOD
	 else if(input == 1 && fSecondInputAODEvent) return (TClonesArray*) fSecondInputAODEvent->FindListObject("mcparticles");	
	 else {
	     printf("AliCaloTrackReader::GetAODMCParticles() - wrong AOD input index? %d, or non existing tree? \n",input); 
		 return 0x0;
	 }
	}
	else {
		printf("AliCaloTrackReader::GetAODMCParticles() - Input are not AODs\n"); 
		return 0x0;
	}
}

//____________________________________________________________________________
AliAODMCHeader* AliCaloTrackReader::GetAODMCHeader(Int_t input) const {
	//Return MC header in AOD. Do it for the corresponding input event.
	if(fDataType == kAOD){
		//Normal input AOD
		if(input == 0) return (AliAODMCHeader*)((AliAODEvent*)fInputEvent)->FindListObject("mcheader");
		//Second input AOD
		else if(input == 1) return  (AliAODMCHeader*) fSecondInputAODEvent->FindListObject("mcheader");	
		else {
			printf("AliCaloTrackReader::GetAODMCHeader() - wrong AOD input index, %d\n",input);
			return 0x0;
		}
	}
	else {
		printf("AliCaloTrackReader::GetAODMCHeader() - Input are not AODs\n");
		return 0x0;
	}
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
	
	if(fSecondInputFileName!=""){
		if(fDataType == kAOD){
			TFile * input2   = new TFile(fSecondInputFileName,"read");
			printf("AliCaloTrackReader::Init() - Second input file opened: %s, size %d \n", input2->GetName(), (Int_t) input2->GetSize());
			fSecondInputAODTree = (TTree*) input2->Get("aodTree");
			if(fSecondInputAODTree) printf("AliCaloTrackReader::Init() - Second input tree opened: %s, entries %d \n", 
										   fSecondInputAODTree->GetName(), (Int_t) fSecondInputAODTree->GetEntries());
			else{
			 printf("AliCaloTrackReader::Init() - Second input tree not available, STOP \n");
			 abort();
			}
			fSecondInputAODEvent = new AliAODEvent;
			fSecondInputAODEvent->ReadFromTree(fSecondInputAODTree);
			if(fSecondInputFirstEvent >= fSecondInputAODTree->GetEntriesFast()){
				printf("AliCaloTrackReader::Init() - Requested first event of second input %d, is larger than number of events %d, STOP\n", 
					   fSecondInputFirstEvent, (Int_t) fSecondInputAODTree->GetEntriesFast());
				abort();
			}
		}
		else printf("AliCaloTrackReader::Init() - Second input not added, reader is not AOD\n");
	}
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

  fFiducialCut           = new AliFiducialCut();
  fSecondInputFileName   = "" ;
  fSecondInputFirstEvent = 0 ;
  fReadStack             = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fReadAODMCParticles    = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fCleanOutputStdAOD     = kFALSE; // Clean the standard clusters/tracks?
  fDeltaAODFileName      = "deltaAODPartCorr.root";
  fFiredTriggerClassName      = "";
	 	
  fAnaLED = kFALSE;
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
  if(fComparePtHardAndJetPt)
	  printf("Compare jet pt and pt hard to accept event, factor = %2.2f",fPtHardAndJetPtFactor);
	
  if(fSecondInputFileName!="") {
	  printf("Second Input File Name     =     %s\n", fSecondInputFileName.Data()) ;
	  printf("Second Input First Event   =     %d\n", fSecondInputFirstEvent) ;
  }
	
  printf("Read Kine from, stack? %d, AOD ? %d \n", fReadStack, fReadAODMCParticles) ;
  printf("Clean std AOD       =     %d\n", fCleanOutputStdAOD) ;
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
	if(eventType!=7)return kFALSE; //Only physics event, do not use for simulated events!!!
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
		
  if(fOutputEvent && (fDataType != kAOD) && ((fOutputEvent->GetCaloClusters())->GetEntriesFast()!=0 ||(fOutputEvent->GetTracks())->GetEntriesFast()!=0)){
    if (fFillCTS || fFillEMCAL || fFillPHOS) {
      printf("AliCaloTrackReader::AODCaloClusters or AODTracks already filled by the filter, do not use the ESD reader, use the AOD reader, STOP\n");
		  abort();
    }
  }

  //In case of analysis of events with jets, skip those with jet pt > 5 pt hard	
  if(fComparePtHardAndJetPt && GetStack()) {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
  }

  //In case of mixing events with other AOD file	
  if(fDataType == kAOD && fSecondInputAODTree){
	 
    if(fDebug > 1) 
      printf("AliCaloTrackReader::FillInputEvent() - Get event %d from second input AOD file \n", iEntry+fSecondInputFirstEvent);
    if(fSecondInputAODTree->GetEntriesFast() <= iEntry+fSecondInputFirstEvent) {
      if(fSecondInputAODTree->GetEntriesFast() == iEntry+fSecondInputFirstEvent) 
			 printf("AliCaloTrackReader::FillInputEvent() - Skip events from event %d, no more events in second AOD file \n", iEntry);
      return kFALSE;
    }
    
    //Get the Event
    Int_t nbytes = fSecondInputAODTree->GetEvent(iEntry+fSecondInputFirstEvent);
    if ( nbytes == 0 ) {//If nothing in AOD
      printf("AliCaloTrackReader::FillInputEvent() - Nothing in Second AOD input, STOP\n");
      abort() ; 
    }
    
  }
	
  if(fFillEMCALCells) FillInputEMCALCells();
  if(fFillPHOSCells)  FillInputPHOSCells();
	
  if(fFillCTS)   FillInputCTS();
  if(fFillEMCAL) FillInputEMCAL();
  if(fFillPHOS)  FillInputPHOS();

	
  return kTRUE ;
}

//__________________________________________________
void AliCaloTrackReader::ResetLists() {
  //  Reset lists, called by the analysis maker 

  if(fAODCTS)   fAODCTS -> Clear();
  if(fAODEMCAL) fAODEMCAL -> Clear();
  if(fAODPHOS)  fAODPHOS -> Clear();
  if(fEMCALCells) fEMCALCells -> Clear();
  if(fPHOSCells)  fPHOSCells -> Clear();
  if(fCleanOutputStdAOD && fOutputEvent ){
    //Only keep copied tracks and clusters if requested
    fOutputEvent->GetTracks()      ->Clear();
    fOutputEvent->GetCaloClusters()->Clear();
  }

}
