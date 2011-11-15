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
/* $Id$ */

//_________________________________________________________________________
// Analysis task that executes the analysis classes
// that depend on the PartCorr frame, frame for Particle identification and correlations.
// Specially designed for calorimeters but also can be used for charged tracks
// Input of this task is a configuration file that contains all the settings of the analyis
//
// -- Author: Gustavo Conesa (INFN-LNF)

#include <cstdlib>

// --- Root ---
#include <TROOT.h>
#include <TInterpreter.h>
#include <TClonesArray.h>
//#include <Riostream.h>
//#include <TObjectTable.h>

// --- Analysis ---
#include "AliAnalysisTaskParticleCorrelationM.h"
#include "AliAnaPartCorrMaker.h"
#include "AliCaloTrackReader.h"
#include "AliPDG.h"
#include "AliAnalysisManager.h"
#include "AliMultiEventInputHandler.h"
#include "AliMixedEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisDataSlot.h"

ClassImp(AliAnalysisTaskParticleCorrelationM)

////////////////////////////////////////////////////////////////////////
AliAnalysisTaskParticleCorrelationM::AliAnalysisTaskParticleCorrelationM():
  AliAnalysisTaskME(),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), fCuts(0x0), fInputEvent(NULL)
{
  // Default constructor
}

//_____________________________________________________
AliAnalysisTaskParticleCorrelationM::AliAnalysisTaskParticleCorrelationM(const char* name):
  AliAnalysisTaskME(name),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), fCuts(0x0), fInputEvent(NULL)
{
  // Default constructor
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // will contain cuts or local params
}

//_____________________________________________________
AliAnalysisTaskParticleCorrelationM::~AliAnalysisTaskParticleCorrelationM() 
{
  // Remove all pointers
	
//	printf("********** Delete Task\n");
//  // Do not delete it here, already done somewhere else, need to understand where.
//  if(fOutputContainer && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
//    fOutputContainer->Clear() ; 
//    delete fOutputContainer ;
//  }

  delete fInputEvent ; 
  if(fAna) delete fAna;

//  printf("********** Task deleted \n");
}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelationM::UserCreateOutputObjects()
{
  // Create the output container
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::UserCreateOutputObjects() - Begin\n");
  
  //Get list of aod arrays, add each aod array to analysis frame 
  TClonesArray *array = 0;
  TList * list = fAna->FillAndGetAODBranchList();
  if (DebugLevel() >= 1) printf("AliAnalysisTaskParticleCorrelationM::UserCreateOutputObjects() - n AOD branches %d\n",list->GetEntries());

  //Put the delta AODs in output file, std or delta
  if((fAna->GetReader())->WriteDeltaAODToFile()){
    TString deltaAODName = (fAna->GetReader())->GetDeltaAODFileName();
    for(Int_t iaod = 0; iaod < list->GetEntries(); iaod++){
      array = (TClonesArray*) list->At(iaod);
      if(deltaAODName!="") AddAODBranch("TClonesArray", &array, deltaAODName);//Put it in DeltaAOD file
      else AddAODBranch("TClonesArray", &array);//Put it in standard AOD file
    } 
	}
	
  //Histograms container
  OpenFile(1);
  fOutputContainer = fAna->GetOutputContainer();
  
  if (DebugLevel() >= 1) printf("AliAnalysisTaskParticleCorrelationM::UserCreateOutputObjects() - n histograms %d\n",fOutputContainer->GetEntries());
  
  fOutputContainer->SetOwner(kTRUE);
  
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::UserCreateOutputObjects() - End\n");
 
  PostData(1,fOutputContainer);
	
}
//_____________________________________________________
void AliAnalysisTaskParticleCorrelationM::LocalInit()
{
	// Local Initialization
	
	//Call the Init to initialize the configuration of the analysis
	Init();
	
	// Create cuts/param objects and publish to slot
	fCuts = fAna->GetListOfAnalysisCuts();
  fCuts ->SetOwner(kTRUE);

	// Post Data
	PostData(2, fCuts);
	
}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelationM::Init()
{
  // Initialization
 
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::Init() - Begin\n");
  
  fInputEvent = new AliMixedEvent() ; 
  
  // Call configuration file if specified
  
  if (fConfigName.Length()) {
    printf("AliAnalysisTaskParticleCorrelationM::Init() - ### Configuration file is %s.C ###\n", fConfigName.Data());
	gROOT->LoadMacro(fConfigName+".C");
	fAna = (AliAnaPartCorrMaker*) gInterpreter->ProcessLine("ConfigAnalysis()");
  }
  
  if(!fAna) {
	printf("AliAnalysisTaskParticleCorrelationM::Init() - Analysis maker pointer not initialized, no analysis specified, STOP !\n");
	abort();
  }
  
  // Add different generator particles to PDG Data Base 
  // to avoid problems when reading MC generator particles
  AliPDG::AddParticlesToPdgDataBase();

  //Set in the reader the name of the task in case is needed
  (fAna->GetReader())->SetTaskName(GetName());
	
  // Initialise analysis
  fAna->Init();

  //Delta AOD
  if((fAna->GetReader())->GetDeltaAODFileName()!="")
	  AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile((fAna->GetReader())->GetDeltaAODFileName());

  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::Init() - End\n");
  
}


//_____________________________________________________
void AliAnalysisTaskParticleCorrelationM::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::UserExec() - Begin\n");

   //Get the type of data, check if type is correct
  Int_t  datatype = fAna->GetReader()->GetDataType();
  if(datatype != AliCaloTrackReader::kESD && datatype != AliCaloTrackReader::kAOD &&
     datatype != AliCaloTrackReader::kMC){
    printf("AliAnalysisTaskParticleCorrelationM::UserExec() - Wrong type of data\n");
    return ;
  }

  Int_t nev = fInputHandler->GetBufferSize();
  fInputEvent->Reset();
  
  for (Int_t iev = 0; iev < nev; iev++) {
    if (datatype == AliCaloTrackReader::kESD) {
      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(GetEvent(iev));
      fInputEvent->AddEvent(esd);
    } else if (datatype == AliCaloTrackReader::kAOD) {
      AliAODEvent* aod = dynamic_cast<AliAODEvent*>(GetEvent(iev));
      fInputEvent->AddEvent(aod);
    } else {
      AliFatal("need to implement mixed event for MC") ; 
    }
   }
  fInputEvent->Init();

  fAna->GetReader()->SetInputOutputMCEvent(fInputEvent, AODEvent(), MCEvent());

  //Process event
  fAna->ProcessEvent((Int_t) Entry(), CurrentFileName());
  //printf("AliAnalysisTaskParticleCorrelationM::Current Event %d; Current File Name : %s\n",(Int_t) Entry(), CurrentFileName());
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelationM::UserExec() - End\n");
    
  PostData(1, fOutputContainer);
	
  AliAnalysisDataSlot *out0 = GetOutputSlot(0);
  if (out0 && out0->IsConnected()) PostData(0, fTreeA);  
  
  //gObjectTable->Print();

}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelationM::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  // Do some plots

  // Get merged histograms from the output container
  // Propagate histagrams to maker
  fAna->Terminate((TList*)GetOutputData(1));
	
}

