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

/* $Id:$ */

////////////////////////////////////////////////////////////////////////
//
// Analysis class to Correct Underlying Event studies
// (requires input AOD)
//
// Different analysis are performed according to input.
//
// Run on MC:
// - fraction of diffractive events after different cuts
// - tracking efficiency and contamination
// 
// Run on DATA:
// - fraction of events after different cuts
// - vertex reconstruction efficiency
//
// vallero@physi.uni-heidelberg.de
// 
////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAnalyseUE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCorrectionsUE.h"
#include "AliAODHandler.h"
#include "AliESDHandler.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHistogramsUE.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisHelperJetTasks.h"

ClassImp( AliAnalysisTaskCorrectionsUE)

// Define global pointer
AliAnalysisTaskCorrectionsUE* AliAnalysisTaskCorrectionsUE::fgTaskCorrectionsUE=NULL;

//____________________________________________________________________
AliAnalysisTaskCorrectionsUE:: AliAnalysisTaskCorrectionsUE(const char* name):
AliAnalysisTask(name,""),
fAnaUE(0x0),
fAOD(0x0), 
fESDHandler(0x0),
fAODBranch("jets"),
fCFManager(0x0),
fDebug(0),
fHistosUE(0x0),
fInputHandler(0x0),
fListOfHistos(0x0),  
fMcHandler(0x0),
fMcEvent(0x0),
fBinsPtInHist(30),     
fIsNorm2Area(kTRUE),
fMaxJetPtInHist(300.), 
fMinJetPtInHist(0.),
fConstrainDistance(kTRUE),
fMinDistance(0.2),
fSimulateChJetPt(kFALSE),
fUseAliStack(kTRUE),
fUseMCParticleBranch(kFALSE),
fnTracksVertex(3),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(5.),
fAnaType(1),         
fConePosition(1),
fConeRadius(0.7),
fFilterBit(0xFF),
fJetsOnFly(kFALSE),
fRegionType(1),
fUseChargeHadrons(kFALSE),
fUseChPartJet(kFALSE),
fUsePositiveCharge(kTRUE),
fUseSingleCharge(kFALSE),
fOrdering(1),
fChJetPtMin(5.0),
fJet1EtaCut(0.2),
fJet2DeltaPhiCut(2.616),    // 150 degrees
fJet2RatioPtCut(0.8),
fJet3PtCut(15.),
fTrackEtaCut(0.9),
fTrackPtCut(0.),
fAvgTrials(1)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
  DefineOutput(1, AliCFContainer::Class());

}

//____________________________________________________________________
AliAnalysisTaskCorrectionsUE:: AliAnalysisTaskCorrectionsUE(const AliAnalysisTaskCorrectionsUE & original):
AliAnalysisTask(),
fAnaUE(original.fAnaUE),
fAOD(original.fAOD),            
fESDHandler(original.fESDHandler),            
fAODBranch(original.fAODBranch),
fCFManager(original.fCFManager),
fDebug(original.fDebug),
fHistosUE(original.fHistosUE),
fInputHandler(original.fInputHandler),
fListOfHistos(original.fListOfHistos),  
fMcHandler(original.fMcHandler),
fMcEvent(original.fMcEvent),
fBinsPtInHist(original.fBinsPtInHist),     
fIsNorm2Area(original.fIsNorm2Area),
fMaxJetPtInHist(original.fMaxJetPtInHist), 
fMinJetPtInHist(original.fMinJetPtInHist),
fConstrainDistance(original.fConstrainDistance),
fMinDistance(original.fMinDistance),
fSimulateChJetPt(original.fSimulateChJetPt),
fUseAliStack(original.fUseAliStack),
fUseMCParticleBranch(original.fUseMCParticleBranch),
fnTracksVertex(original.fnTracksVertex),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(original.fZVertex),
fAnaType(original.fAnaType),         
fConePosition(original.fConePosition),
fConeRadius(original.fConeRadius),
fFilterBit(original.fFilterBit),
fJetsOnFly(original.fJetsOnFly),
fRegionType(original.fRegionType),
fUseChargeHadrons(original.fUseChargeHadrons),
fUseChPartJet(original.fUseChPartJet),
fUsePositiveCharge(original.fUsePositiveCharge),
fUseSingleCharge(original.fUseSingleCharge),
fOrdering(original.fOrdering),
fChJetPtMin(original.fChJetPtMin),
fJet1EtaCut(original.fJet1EtaCut),
fJet2DeltaPhiCut(original.fJet2DeltaPhiCut),    // 150 degrees
fJet2RatioPtCut(original.fJet2RatioPtCut),
fJet3PtCut(original.fJet3PtCut),
fTrackEtaCut(original.fTrackEtaCut),
fTrackPtCut(original.fTrackPtCut),
fAvgTrials(original.fAvgTrials)
{
  // Copy constructor
}


//______________________________________________________________
AliAnalysisTaskCorrectionsUE & AliAnalysisTaskCorrectionsUE::operator = (const AliAnalysisTaskCorrectionsUE & /*source*/)
{
  // assignment operator
  return *this;
}

//______________________________________________________________
Bool_t AliAnalysisTaskCorrectionsUE::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // Copy from AliAnalysisTaskJFSystematics
  fAvgTrials = 1;
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t trials  = 1;
  if(tree){
  	TFile *curfile = tree->GetCurrentFile();
  	if (!curfile) {
  		Error("Notify","No current file");
  		return kFALSE;
  		}
  	
	AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,trials); 
        fHistosUE->GetXsec()->Fill("<#sigma>",xsection);

  	// construct average trials
  	Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  	if(trials>=nEntries && nEntries>0.)fAvgTrials = trials/nEntries;
  	}
  
  return kTRUE;

}

//____________________________________________________________________
void AliAnalysisTaskCorrectionsUE::ConnectInputData(Option_t* /*option*/)
{
  // Connect the input data  
  
  // We need AODs with tracks and jets.
  // Since AODs can either be connected to the InputEventHandler 
  // or to the OutputEventHandler ( the AOD is created by a previous task in the train )
  // we need to check where it is and get the pointer to AODEvent in the right way
  
  // Delta AODs are also accepted
  
 
  if (fDebug > 1) AliInfo("ConnectInputData() ");
  
  //Get the input handler
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  if( handler && handler->InheritsFrom("AliAODInputHandler") ) { // input AOD
	fAOD = ((AliAODInputHandler*)handler)->GetEvent();
  	if(!fJetsOnFly){
		if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODInputHandler");
		}else{
		if (fDebug > 1) AliInfo(" ==== Tracks from AliAODInputHandler / Jets on-the-fly");
		}
  	} else {  //output AOD
  	handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
  	if( handler && handler->InheritsFrom("AliAODHandler") ) {
  		fAOD = ((AliAODHandler*)handler)->GetAOD();
  		if (!fJetsOnFly){
  			if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODHandler");
  			} else {
  			if (fDebug > 1) AliInfo(" ==== Tracks from AliAODHandler / Jets on-the-fly");
  			}
  		}else {
  		AliFatal("I can't get any AOD Event Handler");
  		return;
  		}
  	}

  //Get ESD input handler
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   AliESDEvent *esdEvent = (AliESDEvent*)inputHandler->GetEvent();
   if (!esdEvent && fDebug > 1) {
  	AliInfo("********************** No ESD event: cannot retrive DCA values from AOD !!! ");
  }else fAnaUE->SetESDEvent(esdEvent);

  //Get MC handler 
  fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
 
  //Initialize AliAnalysisUE class
  fAnaUE->Initialize(fAnaType, fAOD, fConeRadius, fDebug, fFilterBit, fJet1EtaCut, fJet2DeltaPhiCut, fJet2RatioPtCut, fJet3PtCut, fOrdering, fRegionType, fSimulateChJetPt, fTrackEtaCut, fTrackPtCut, fUseChargeHadrons, fUseChPartJet, fUsePositiveCharge, fUseSingleCharge, fHistosUE);
 
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::CreateOutputObjects()
{
  // Create the output container
  
  if (fDebug > 1) AliInfo("CreateOutPutData()");
   
  // Create pointer to AliAnalysisUE, a class implementing the main analysis algorithms
  fAnaUE = new AliAnalyseUE();
  if (!fAnaUE){
     AliError("UE analysis class not initialized!!!");
     return;
  }
  // Create pointer to AliHistogramsUE, a class handling histogram creation/filling
  fHistosUE = new AliHistogramsUE();
  if (!fHistosUE){
     AliError("UE histograms class not initialized!!!");
     return;
  }

  // Create list of output histograms
  if (fListOfHistos != NULL){
	delete fListOfHistos;
        fListOfHistos = NULL;
  	}
  if (!fListOfHistos){
  	fListOfHistos = new TList();
  	fListOfHistos->SetOwner(kTRUE); 
  	}
  
  
  //Initialize output histograms
  fHistosUE->CreateHistogramsCorrections(fListOfHistos,fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, fTrackEtaCut);
  AddSettingsTree();

  //Configure the CF manager
  if (fCFManager != NULL){
	delete fCFManager;
        fCFManager = NULL;
  	}
  if (!fCFManager){
        fCFManager = new AliCFManager();
  	}
  fHistosUE->CreateCorrectionsContainer(fCFManager,fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, fTrackEtaCut,fJet1EtaCut);

  //Post outputs
  PostData(0,fListOfHistos);
  PostData(1,fCFManager->GetEventContainer());

}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::Exec(Option_t */*option*/)
{

  Bool_t flag=kTRUE;

  //Determine the corrections
  if (fMcHandler){
  	fMcEvent = fMcHandler->MCEvent();
       	if ( fDebug > 3 ) AliInfo( " Processing MC event..." );
  }else{
       	if ( fDebug > 3 ) AliInfo( " Processing DATA event..." );
	}
  

  flag = EvaluateCorrections(); 

  
  if (flag){ //executed if event passes trigger+vertex selection
  	//Fetch the pythia header info and get the trials
  	Float_t nTrials = 1;
  	if (fMcHandler && fMcEvent) {
  		AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(fMcEvent);
  		if(pythiaGenHeader) nTrials = pythiaGenHeader->Trials();
  		}
  	fHistosUE->GetTrials()->Fill("#sum{ntrials}",fAvgTrials);
 	PostData(0,fListOfHistos);
 	PostData(1,fCFManager->GetEventContainer());
 	}
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UECorrectionsSettings","Analysis Settings in UE corrections");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fConeRadius", &fConeRadius,"Rad/D");
  settingsTree->Branch("fJet1EtaCut", &fJet1EtaCut, "LeadJetEtaCut/D");
  settingsTree->Branch("fJet2DeltaPhiCut", &fJet2DeltaPhiCut, "DeltaPhi/D");
  settingsTree->Branch("fJet2RatioPtCut", &fJet2RatioPtCut, "Jet2Ratio/D");
  settingsTree->Branch("fJet3PtCut", &fJet3PtCut, "Jet3PtCut/D");
  settingsTree->Branch("fTrackPtCut", &fTrackPtCut, "TrackPtCut/D");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fAnaType", &fAnaType, "Ana/I");        
  settingsTree->Branch("fRegionType", &fRegionType,"Reg/I");
  settingsTree->Branch("fOrdering", &fOrdering,"OrderMeth/I");
  settingsTree->Branch("fUseChPartJet", &fUseChPartJet,"UseChPart/O");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fUseSingleCharge", &fUseSingleCharge,"UseSingleCh/O");
  settingsTree->Branch("fUsePositiveCharge", &fUsePositiveCharge,"UsePositiveCh/O");
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
} 


//____________________________________________________________________
Bool_t  AliAnalysisTaskCorrectionsUE::EvaluateCorrections()
{
  
  /////////////////////////////////////////////////////////////////////////
  // 
  // EVENT SELECTION:
  // CF containers are filled to get the number of entries after every cut.
  // *** Cuts: ***
  // 0 - triggered
  // 1 - physics selection
  // 2 - vertex selection
  // 3 - event topology
  // 4 - leading track pT cut
  // 5 - leading track correctly identified
  // *** Variables: ***
  // 0 - leading track pT (reco)
  // 1 - leading track eta (reco)
  // 2 - process type:
  //	1: non diffractive
  //	2: double diffractive
  //	4: single diffractive
  // 3 - leading track pT (true)
  // 4 - leading track eta (true)
  // 5 - delta eta (reco-true)
  // 6 - delta phi (reco-true)
  // 7 - radius (reco-true)
  //
  // TRACK-LEVEL CORRECTIONS:
  // Fill histograms similar to AliAnalysisTaskUE.
  //
  /////////////////////////////////////////////////////////////////////////

  Double_t containerInput[8];// relevant variables (see above)  
  
  //PROCESS TYPE (ND,SD,DD)
  AliAnalysisHelperJetTasks::MCProcessType eventId = AliAnalysisHelperJetTasks::kInvalidProcess;
  if (fMcHandler && fMcEvent) {
  		AliGenEventHeader* genHeader = fMcEvent->GenEventHeader();
		eventId = AliAnalysisHelperJetTasks::GetPythiaEventProcessType(genHeader,kFALSE);
		if (eventId<0){
		eventId = AliAnalysisHelperJetTasks::GetDPMjetEventProcessType(genHeader,kFALSE);
		}
		if (eventId<0 && fDebug>1)AliInfo("No Pythia or Phojet Header retrived!"); 
  }else if (fDebug>1) AliInfo("No MC handler or Event!"); 

  //Initialize container inputs 
  for (Int_t i =0; i<8; i++){  
  	containerInput[i]=-999.;
  	}
  
  //Assign process type
  if (eventId == 1 ) containerInput[2]=1.; //Non diffractive
  if (eventId == 2 ) containerInput[2]=2.; //Double diffractive
  if (eventId == 4 ) containerInput[2]=4.; //Single diffractive
  
  // Execute analysis for current event ******************************
  
  // Get jets and order by pT
  TVector3 jetVect[3];
  *jetVect = fAnaUE->GetOrderedClusters(fAODBranch, fUseChPartJet, fChJetPtMin );
  
  //now define leading track pT and eta
  containerInput[0]=jetVect[0].Pt();
  containerInput[1]=jetVect[0].Eta();

  fCFManager->GetEventContainer()->Fill(containerInput,kCFStepTriggered);  //fill CF container 
  
  // Physics selection ************************************************
  fInputHandler = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!fAnaUE->TriggerSelection(fInputHandler)) return kFALSE;
  
  fCFManager->GetEventContainer()->Fill(containerInput, kCFStepPhysSelect); //fill CF container
  
  // Event selection (vertex) *****************************************
  
  if(!fAnaUE->VertexSelection(fAOD,fnTracksVertex,fZVertex)) return kFALSE;
  //if(!fAnaUE->VertexSelectionOld(fAOD)) return kFALSE; // temporary to compare with old task and to have same cuts for MC !!!

  fCFManager->GetEventContainer()->Fill(containerInput, kCFStepVertexSelect); //fill CF container
   
   // Select events according to analysis type *************************
   // (in the leading track analysis it should not happen that there are no "jets")
   if( jetVect[0].Pt() < 0. ) {
   	if( fDebug > 1 ) AliInfo("\n   Skipping Event, not jet found...");
      	return kTRUE;
    } else {
	if (fDebug >1 ) AliInfo(Form("\n   Pt Leading Jet = %6.1f eta=%5.3f ",  jetVect[0].Pt(), jetVect[0].Eta() ));
 	}
  
  
  if ( ! (fAnaUE->AnaTypeSelection(jetVect))) return kTRUE;
   
   
   // For the events selected check the real MC leading particle
   // (only when running on MC)
   TVector3 jetVectMC[3];
   if (fMcEvent){
   	*jetVectMC = fAnaUE->GetLeadingTracksMC(fMcEvent); 
   	if (fAnaUE->AnaTypeSelection(jetVectMC) ){
   		//now define leading track pT and eta
   		containerInput[3]=jetVectMC[0].Pt();
   		containerInput[4]=jetVectMC[0].Eta();
   		//Check distance between real and reco leading-particle
        	containerInput[5] = jetVect[0].Eta()-jetVectMC[0].Eta();
		containerInput[6] = TMath::Abs(jetVect[0].Phi()-jetVectMC[0].Phi());
		if (containerInput[6] >= 2.*TMath::Pi())containerInput[6] -= 2.*TMath::Pi(); 
   		containerInput[7] = sqrt(TMath::Power(containerInput[5],2.) + TMath::Power(containerInput[6],2.));
		} 
   	}  
   
   fCFManager->GetEventContainer()->Fill(containerInput, kCFStepAnaTopology); //fill CF container
   
   //Cut on the leading track pT 
   if (!(jetVect[0].Pt() >= 1.)) return kTRUE;
   // (only when running on MC)   
   if (fMcEvent){ 
   	if (!(jetVectMC[0].Pt()>= 1.)){
		containerInput[3]=containerInput[4]=-999.;
		containerInput[5]=containerInput[6]=containerInput[7]=-999.;
		}
	}	
   
   fCFManager->GetEventContainer()->Fill(containerInput, kCFStepLtPtCut1); //fill CF container

   // Check if leading track is correctly identified
   // (only when running on MC)
   if (fMcEvent){
   	Int_t labelLt = fAnaUE->GetLtLabel();
   	Int_t labelLtMC = fAnaUE->GetLtMCLabel();
   	if (labelLt == labelLtMC){
        	fCFManager->GetEventContainer()->Fill(containerInput, kCFStepLtCorrect); //fill CF container
		} 
   	}

   // For track efficiency and contamination
   // (only when running on MC)
   if (fMcEvent){
   	//Run once on reco...
   	fAnaUE->FindMaxMinRegions( jetVect,  fConePosition, 0, 0 );//normal track cut 
   	fAnaUE->FindMaxMinRegions( jetVect,  fConePosition, 0, 1 );//for efficiency: cut on pT and eta are set on MC true  
   	//and once on MC true
        fAnaUE->FindMaxMinRegions( jetVect,  fConePosition, 1, 0 );
  }else{
  	// For d0 distribution, runs only on real data
   	fAnaUE->FindMaxMinRegions( jetVect,  fConePosition, 0, 0 );//run on real data 
  	}


   return kTRUE;
}

//____________________________________________________________________
AliAnalysisTaskCorrectionsUE* AliAnalysisTaskCorrectionsUE::Instance()
{ 
  
  //Create instance
  if (fgTaskCorrectionsUE) {
	return fgTaskCorrectionsUE;
  } else {
	fgTaskCorrectionsUE = new AliAnalysisTaskCorrectionsUE();
	return fgTaskCorrectionsUE;
  	}
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  
  if (fDebug >1) AliAnalysisHelperJetTasks::PrintDirectorySize("PWG4_JetTasksOutput.root");
  
  if (!gROOT->IsBatch()){
  	fListOfHistos = dynamic_cast<TList*> (GetOutputData(0));
  	if (!fListOfHistos){
		AliError("Histogram List is not available");
		return;
  		}

  } else {
        AliInfo(" Batch mode, not histograms will be shown...");
  }

  
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::WriteSettings()
{ 

}
