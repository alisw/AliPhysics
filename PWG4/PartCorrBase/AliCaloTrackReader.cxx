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
#include "TGeoManager.h"

//---- ANALYSIS system ----
#include "AliCaloTrackReader.h"
#include "AliFiducialCut.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"

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
    fEMCALGeoName("EMCAL_COMPLETE"),fPHOSGeoName("PHOSgeo"), 
    fEMCALGeo(0x0), fPHOSGeo(0x0), fEMCALGeoMatrixSet(kFALSE), fPHOSGeoMatrixSet(kFALSE), fAnaLED(kFALSE),
    fRemoveBadChannels(kFALSE),fEMCALBadChannelMap(new TObjArray()),
    fPHOSBadChannelMap(new TObjArray()), fTaskName("")
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliCaloTrackReader::AliCaloTrackReader(const AliCaloTrackReader & g) :   
  TObject(g), fEventNumber(g.fEventNumber), fCurrentFileName(g.fCurrentFileName), 
  fDataType(g.fDataType), fDebug(g.fDebug),
  fFiducialCut(g.fFiducialCut),
  fComparePtHardAndJetPt(g.fComparePtHardAndJetPt),
  fPtHardAndJetPtFactor(g.fPtHardAndJetPtFactor),
  fCTSPtMin(g.fCTSPtMin), fEMCALPtMin(g.fEMCALPtMin),fPHOSPtMin(g.fPHOSPtMin), 
  fAODCTS(new TObjArray(*g.fAODCTS)),  
  fAODEMCAL(new TObjArray(*g.fAODEMCAL)),
  fAODPHOS(new TObjArray(*g.fAODPHOS)),
  fEMCALCells(new TNamed(*g.fEMCALCells)),
  fPHOSCells(new TNamed(*g.fPHOSCells)),
  fInputEvent(g.fInputEvent), fOutputEvent(g.fOutputEvent), fMC(g.fMC),
  fFillCTS(g.fFillCTS),fFillEMCAL(g.fFillEMCAL),fFillPHOS(g.fFillPHOS),
  fFillEMCALCells(g.fFillEMCALCells),fFillPHOSCells(g.fFillPHOSCells),
  fSecondInputAODTree(g.fSecondInputAODTree), 
  fSecondInputAODEvent(g.fSecondInputAODEvent),
  fSecondInputFileName(g.fSecondInputFileName), 
  fSecondInputFirstEvent(g.fSecondInputFirstEvent),
  fAODCTSNormalInputEntries(g.fAODCTSNormalInputEntries), 
  fAODEMCALNormalInputEntries(g.fAODEMCALNormalInputEntries), 
  fAODPHOSNormalInputEntries(g.fAODPHOSNormalInputEntries),
  fTrackStatus(g.fTrackStatus),
  fReadStack(g.fReadStack), fReadAODMCParticles(g.fReadAODMCParticles),
  fCleanOutputStdAOD(g.fCleanOutputStdAOD), fDeltaAODFileName(g.fDeltaAODFileName),
  fFiredTriggerClassName(g.fFiredTriggerClassName),
  fEMCALGeoName(g.fEMCALGeoName),           fPHOSGeoName(g.fPHOSGeoName),
  fEMCALGeo(new AliEMCALGeoUtils(*g.fEMCALGeo)), fPHOSGeo(new AliPHOSGeoUtils(*g.fPHOSGeo)),
  fEMCALGeoMatrixSet(g.fEMCALGeoMatrixSet), fPHOSGeoMatrixSet(g.fPHOSGeoMatrixSet),
  fAnaLED(g.fAnaLED),  fRemoveBadChannels(g.fRemoveBadChannels),
  fEMCALBadChannelMap(new TObjArray(*g.fEMCALBadChannelMap)),
  fPHOSBadChannelMap(new TObjArray(*g.fPHOSBadChannelMap)),
  fTaskName(g.fTaskName)
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
//  fEMCALGeoName      = source.fEMCALGeoName ; 
//  fPHOSGeoName       = source.fPHOSGeoName ; 
//  fEMCALGeo          = new AliEMCALGeoUtils(*source.fEMCALGeo);  
//  fPHOSGeo           = new AliPHOSGeoUtils(*source.fPHOSGeo);
//  fEMCALGeoMatrixSet = source.fEMCALGeoMatrixSet; 
//  fPHOSGeoMatrixSet  = source.fPHOSGeoMatrixSet;
//  fAnaLED            = source.fAnaLED;
//  fRemoveBadChannels = source.fRemoveBadChannels;
//  fEMCALBadChannelMap= source.fEMCALBadChannelMap;
//  fPHOSBadChannelMap = source.fPHOSBadChannelMap;
//
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

  if(fInputEvent)  delete fInputEvent ;
  if(fOutputEvent) delete fOutputEvent ;
  if(fMC)          delete fMC ;  
	
  if(fSecondInputAODTree){
    fSecondInputAODTree->Clear();
    delete fSecondInputAODTree;
  }
	
  if(fSecondInputAODEvent) delete fSecondInputAODEvent ;

  if(fPHOSGeo)  delete fPHOSGeo  ;
  if(fEMCALGeo) delete fEMCALGeo ;
	
  if(fEMCALBadChannelMap) { 
    fEMCALBadChannelMap->Clear();
    delete  fEMCALBadChannelMap;
  }
  if(fPHOSBadChannelMap) { 
    fPHOSBadChannelMap->Clear();
    delete  fPHOSBadChannelMap;
  }

  //fEMCALBadChannelMap. Delete();
  //fPHOSBadChannelMap. Delete();
	
}

//_________________________________________________________________________________________________________
Bool_t AliCaloTrackReader::ClusterContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells){
	// Check that in the cluster cells, there is no bad channel of those stored 
	// in fEMCALBadChannelMap or fPHOSBadChannelMap
	
	if (!fRemoveBadChannels) return kFALSE;
	
	if(calorimeter == "EMCAL" && !fEMCALBadChannelMap->GetEntries()) return kFALSE;
	if(calorimeter == "PHOS"  && !fPHOSBadChannelMap ->GetEntries()) return kFALSE;
	printf("hello\n");
	Int_t icol = -1;
	Int_t irow = -1;
	Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
	
		//Get the column and row
		if(calorimeter == "EMCAL"){
			Int_t iTower = -1, iIphi = -1, iIeta = -1; 
			fEMCALGeo->GetCellIndex(cellList[iCell],imod,iTower,iIphi,iIeta); 
			if(fEMCALBadChannelMap->GetEntries() <= imod) continue;
			fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);			
			if(GetEMCALChannelStatus(imod, icol, irow))return kTRUE;
		}
		else if(calorimeter=="PHOS"){
			Int_t    relId[4];
			fPHOSGeo->AbsToRelNumbering(cellList[iCell],relId);
			irow = relId[2];
			icol = relId[3];
			imod = relId[0]-1;
			if(fPHOSBadChannelMap->GetEntries() <= imod)continue;
			if(GetPHOSChannelStatus(imod, icol, irow)) return kTRUE;
		}
		else return kFALSE;
		
	}// cell cluster loop
	printf("hello 2\n");
	return kFALSE;

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

	//printf("TaskName in reader: %s\n",fTaskName.Data());
	fEMCALBadChannelMap->SetName(Form("EMCALBadMap_%s",fTaskName.Data()));
	fPHOSBadChannelMap->SetName(Form("PHOSBadMap_%s",fTaskName.Data()));
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
  fEMCALGeoName = "EMCAL_COMPLETE";
  fPHOSGeoName  = "PHOSgeo";
	
  if(gGeoManager) {// geoManager was set
	if(fDebug > 2)printf("AliCaloTrackReader::InitParameters() - Geometry manager available\n");
	fEMCALGeoMatrixSet = kTRUE;	 
	fPHOSGeoMatrixSet  = kTRUE;	 
  }
  else{
	fEMCALGeoMatrixSet = kFALSE;
	fPHOSGeoMatrixSet  = kFALSE;
  }
	
  fAnaLED = kFALSE;
	
  fRemoveBadChannels = kFALSE;
}

//________________________________________________________________
void AliCaloTrackReader::InitEMCALBadChannelStatusMap(){
  //Init EMCAL bad channels map
   if(fDebug > 0 )printf("AliCaloTrackReader::InitEMCALBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  for (int i = 0; i < 12; i++) fEMCALBadChannelMap->Add(new TH2I(Form("EMCALBadChannelMap_SM%d_%s",i,fTaskName.Data()),Form("EMCALBadChannelMap_SM%d",i),  48, 0, 48, 24, 0, 24));
  
  fEMCALBadChannelMap->SetOwner(kTRUE);
  fEMCALBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
void AliCaloTrackReader::InitPHOSBadChannelStatusMap(){
  //Init PHOS bad channels map
  if(fDebug > 0 )printf("AliCaloTrackReader::InitPHOSBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOSBadChannelMap_Mod%d_%s",i,fTaskName.Data()),Form("PHOSBadChannelMap_Mod%d",i), 56, 0, 56, 64, 0, 64));
    
  fPHOSBadChannelMap->SetOwner(kTRUE);
  fPHOSBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
void AliCaloTrackReader::InitEMCALGeometry()
{
	//Initialize EMCAL geometry if it did not exist previously
	if (!fEMCALGeo){
		fEMCALGeo = new AliEMCALGeoUtils(fEMCALGeoName); 
		if (!gGeoManager && fDebug > 0) printf("AliCaloTrackReader::InitEMCALGeometry() - Careful!, gGeoManager not loaded, load misalign matrices\n");
	}
}

//________________________________________________________________
void AliCaloTrackReader::InitPHOSGeometry()
{
	//Initialize PHOS geometry if it did not exist previously
	if (!fPHOSGeo){
		fPHOSGeo = new AliPHOSGeoUtils(fPHOSGeoName); 
		if (!gGeoManager && fDebug > 0) printf("AliCaloTrackReader::InitPHOSGeometry() - Careful!, gGeoManager not loaded, load misalign matrices\n");
	}	
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
  printf("Remove Clusters with bad channels? %d\n",fRemoveBadChannels);
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
	
  //Get the EMCAL transformation geometry matrices from ESD 
  if (!gGeoManager && fEMCALGeo) {//&& !fEMCALGeoMatrixSet) {
    if(fDebug > 1) 
      printf(" AliCaloTrackReader::FillInputEvent() - Load EMCAL misalignment matrices. \n");
    if(!strcmp(fInputEvent->GetName(),"AliESDEvent"))  {
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ 
	if(((AliESDEvent*)fInputEvent)->GetEMCALMatrix(mod)) {
	  //printf("EMCAL: mod %d, matrix %p\n",mod, ((AliESDEvent*)fInputEvent)->GetEMCALMatrix(mod));
	  fEMCALGeo->SetMisalMatrix(((AliESDEvent*)fInputEvent)->GetEMCALMatrix(mod),mod) ;
	  fEMCALGeoMatrixSet = kTRUE;//At least one, so good
	}
      }// loop over super modules	
    }//ESD as input
    else {
      if(fDebug > 1)
	printf("AliCaloTrackReader::FillInputEvent() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
	  }//AOD as input
  }//EMCAL geo && no geoManager
  
  //Get the PHOS transformation geometry matrices from ESD 
  if (!gGeoManager && fPHOSGeo && !fPHOSGeoMatrixSet) {
    if(fDebug > 1) 
      printf(" AliCaloTrackReader::FillInputEvent() - Load PHOS misalignment matrices. \n");
    if(!strcmp(fInputEvent->GetName(),"AliESDEvent"))  {
      for(Int_t mod=0; mod < 5; mod++){ 
	if(((AliESDEvent*)fInputEvent)->GetPHOSMatrix(mod)) {
	  //printf("PHOS: mod %d, matrix %p\n",mod, ((AliESDEvent*)fInputEvent)->GetPHOSMatrix(mod));
	  fPHOSGeo->SetMisalMatrix(((AliESDEvent*)fInputEvent)->GetPHOSMatrix(mod),mod) ;
	  fPHOSGeoMatrixSet  = kTRUE; //At least one so good
	}
      }// loop over modules	
    }//ESD as input
    else {
      if(fDebug > 1) 
	printf("AliCaloTrackReader::FillInputEvent() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
    }//AOD as input
  }//PHOS geo	and  geoManager was not set
	
  if(fFillCTS)   FillInputCTS();
  if(fFillEMCAL) FillInputEMCAL();
  if(fFillPHOS)  FillInputPHOS();
  if(fFillEMCALCells) FillInputEMCALCells();
  if(fFillPHOSCells)  FillInputPHOSCells();
	
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
