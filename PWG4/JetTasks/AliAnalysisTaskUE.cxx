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

//class AliAODInputHandler;
//class AliAODHandler;
//class AliLog;
//class TROOT;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TVector3.h>
#include <TH1F.h>

#include "AliAnalyseUE.h"
#include "AliAnalysisTaskUE.h"
#include "AliHistogramsUE.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"

////////////////////////////////////////////////////////////////////////
//
// Analysis class for Underlying Event studies
//
// Look for correlations on the tranverse regions to 
// the leading charged jet
//
// This class needs as input AOD with track and Jets.
// The output is a list of histograms
//
// AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files, so this case class should be 
// in the train after the Jet finder
//
//    Arian.Abrahantes.Quintana@cern.ch 
//    Ernesto.Lopez.Torres@cern.ch
// 
////////////////////////////////////////////////////////////////////////


ClassImp( AliAnalysisTaskUE)

// Define global pointer
AliAnalysisTaskUE* AliAnalysisTaskUE::fgTaskUE=NULL;

//____________________________________________________________________
AliAnalysisTaskUE:: AliAnalysisTaskUE(const char* name):
AliAnalysisTask(name,""),
fAnaUE(0x0),
fAOD(0x0),            
fAODBranch("jets"),
fDebug(0),
fListOfHistos(0x0),  
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

}

//____________________________________________________________________
AliAnalysisTaskUE:: AliAnalysisTaskUE(const AliAnalysisTaskUE & original):
AliAnalysisTask(),
fAnaUE(original.fAnaUE),
fAOD(original.fAOD),            
fAODBranch(original.fAODBranch),
fDebug(original.fDebug),
fListOfHistos(original.fListOfHistos),  
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
AliAnalysisTaskUE & AliAnalysisTaskUE::operator = (const AliAnalysisTaskUE & /*source*/)
{
  // assignment operator
  return *this;
}

//______________________________________________________________
Bool_t AliAnalysisTaskUE::Notify()
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
        
	fAnaUE->FillXsec("<#sigma>",xsection);

  	// construct average trials
  	Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  	if(trials>=nEntries && nEntries>0.)fAvgTrials = trials/nEntries;
  	}
  
  return kTRUE;
}

//____________________________________________________________________
void AliAnalysisTaskUE::ConnectInputData(Option_t* /*option*/)
{
  // Connect the input data  
  
  // We need AOD with tracks and jets.
  // Since AOD can be either connected to the InputEventHandler (input chain fron AOD files)
  // or to the OutputEventHandler ( AOD is create by a previus task in the train )
  // we need to check where it is and get the pointer to AODEvent in the right way
  
  // Delta AODs are also accepted
  
 
  if (fDebug > 1) AliInfo("ConnectInputData() ");
  
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) { // input AOD
  	fAOD = ((AliAODInputHandler*)handler)->GetEvent();
  	if(!fJetsOnFly){
		if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODInputHandler");
		}else{
		AliInfo(" ==== Tracks from AliAODInputHandler / Jets on-the-fly");
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

   fAnaUE->Initialize( *this );
}

//____________________________________________________________________
void  AliAnalysisTaskUE::CreateOutputObjects()
{
  // Create the output container
  
  if (fDebug > 1) AliInfo("CreateOutPutData()");
   
  //Initialize AliAnalysisUE, a class implementing the main analysis algorithms
  fAnaUE = new AliAnalyseUE();
  
  //Initialize output histograms
  fAnaUE->CreateHistograms(fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist);
  fListOfHistos = (TList*)fAnaUE->GetHistograms()->Clone();
  fListOfHistos->SetOwner(kTRUE); 
  PostData(0,fListOfHistos);

}

//____________________________________________________________________
void  AliAnalysisTaskUE::Exec(Option_t */*option*/)
{
  // Trigger selection ************************************************
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
         ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if( !inputHandler->InheritsFrom("AliAODInputHandler") ) { // input AOD
  	if (inputHandler->IsEventSelected()) {
  		if (fDebug > 1) AliInfo(" Trigger Selection: event ACCEPTED ... ");
  	} else {
  		if (fDebug > 1) AliInfo(" Trigger Selection: event REJECTED ... ");
  		return;
  		}
        }                                
  // Event selection (vertex) *****************************************
   
  if(!fAnaUE->VertexSelection(fAOD,fnTracksVertex,fZVertex)) return;
  
  // Execute analysis for current event ******************************
  
  if ( fDebug > 3 ) AliInfo( " Processing event..." );

  // fetch the pythia header info and get the trials
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  Float_t nTrials = 1;
  if (mcHandler) {
  	AliMCEvent* mcEvent = mcHandler->MCEvent();
  	if (mcEvent) {
  		AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
  		if(pythiaGenHeader){
  			nTrials = pythiaGenHeader->Trials();
  			}
  		}
  	}
  fAnaUE->FillTrials("#sum{ntrials}",fAvgTrials);
  
  //analyse the event
  AnalyseUE();
  fListOfHistos = (TList*)fAnaUE->GetHistograms()->Clone();
  PostData(0,fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskUE::AnalyseUE()
{

  
   // Get jets and order by pT
   TVector3 jetVect[3];
   *jetVect = fAnaUE->GetOrderedClusters(fAODBranch, fUseChPartJet, fChJetPtMin);
  
   if( jetVect[0].Pt() < 0. ) {
   	if( fDebug > 1 ) AliInfo("\n   Skipping Event, not jet found...");
      	return;
    	} else {
		if (fDebug >1 ) AliInfo(Form("\n   Pt Leading Jet = %6.1f eta=%5.3f ",  jetVect[0].Pt(), jetVect[0].Eta() ));
 		}

   // Select events according to analysis type
   if ( ! (fAnaUE->AnaTypeSelection(jetVect))) return;

  // Find max and min regions with real tracks
  if (!fUseMCParticleBranch){
  	// Primary vertex distribution
    	AliAODVertex* vertex = (AliAODVertex*)fAOD->GetPrimaryVertex();
	fAnaUE->FillVertex(vertex->GetNContributors());
	
	// Fill leading "jet" histogram
  	fAnaUE->FillLeadingJet(jetVect[0].Pt());

  	fAnaUE->FindMaxMinRegions( jetVect,  fConePosition );

  }else { 
    
    // this is the part we only use when we have MC information
    // More than a test for values of it also resumes the reconstruction efficiency of jets
    // As commented bellow if available for the data, we try to pair reconstructed jets with simulated ones
    // afterwards we kept angular variables of MC jet to perform UE analysis over MC particles
    // TODO: Handle Multiple jet environment. 06/2009 just suited for inclusive jet condition ( fAnaType = 1 ) 
    
      AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!mcHandler) {
        Printf("ERROR: Could not retrieve MC event handler");
        return;
      }
      
      AliMCEvent* mcEvent = mcHandler->MCEvent();
      if (!mcEvent) {
        Printf("ERROR: Could not retrieve MC event");
        return;
      }
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    if(!pythiaGenHeader){
      return;
    }
    
    fAnaUE->AnalyseMC(jetVect,mcEvent,pythiaGenHeader, fConePosition, fUseAliStack, fConstrainDistance, fMinDistance);

  }

  fAnaUE->FillRegions(fIsNorm2Area, jetVect);
  
}

//____________________________________________________________________
AliAnalysisTaskUE* AliAnalysisTaskUE::Instance()
{ 
  
  //Create instance
  if (fgTaskUE) {
	return fgTaskUE;
  } else {
	fgTaskUE = new AliAnalysisTaskUE();
	return fgTaskUE;
  	}
}

//____________________________________________________________________
void  AliAnalysisTaskUE::Terminate(Option_t */*option*/)
{
  
  // Terminate analysis

  fListOfHistos = dynamic_cast<TList*> (GetOutputData(0));
  if (!fListOfHistos){
	AliError("Histogram List is not available");
	return;
  }
  if (!gROOT->IsBatch()){
  //call class AliHistogramsUE
  AliHistogramsUE *histos=new AliHistogramsUE(fListOfHistos);
  histos->DrawUE(0);
  } else {
        AliInfo(" Batch mode, not histograms will be shown...");
  }

  if( fDebug > 1 ) 
    AliInfo("End analysis");
 
}

void  AliAnalysisTaskUE::WriteSettings()
{ 

  //Print analysis settings on screen
  if (fDebug > 5){
    AliInfo(" All Analysis Settings in Saved Tree");
    fAnaUE->WriteSettings();
  }
}
