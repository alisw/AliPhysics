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
/* $Id: $ */

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
#include "AliAnalysisTaskParticleCorrelation.h"
#include "AliAnaPartCorrMaker.h"
#include "AliCaloTrackReader.h"
#include "AliPDG.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisTaskParticleCorrelation)

////////////////////////////////////////////////////////////////////////
AliAnalysisTaskParticleCorrelation::AliAnalysisTaskParticleCorrelation():
  AliAnalysisTaskSE(),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), fCuts(0x0)
{
  // Default constructor
}

//_____________________________________________________
AliAnalysisTaskParticleCorrelation::AliAnalysisTaskParticleCorrelation(const char* name):
  AliAnalysisTaskSE(name),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), fCuts(0x0)
{
  // Default constructor
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // will contain cuts or local params
}

//_____________________________________________________
AliAnalysisTaskParticleCorrelation::~AliAnalysisTaskParticleCorrelation() 
{
  // Remove all pointers
	
//	printf("********** Delete Task\n");
//  // Do not delete it here, already done somewhere else, need to understand where.
//  if(fOutputContainer){
//    fOutputContainer->Clear() ; 
//    delete fOutputContainer ;
//  }

  if(fAna) delete fAna;

//  printf("********** Task deleted \n");
}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelation::UserCreateOutputObjects()
{
  // Create the output container
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::UserCreateOutputObjects() - Begin\n");
  
  //Get list of aod arrays, add each aod array to analysis frame 
  TClonesArray *array = 0;
  TList * list = fAna->GetAODBranchList();
  TString deltaAODName = (fAna->GetReader())->GetDeltaAODFileName();
  for(Int_t iaod = 0; iaod < list->GetEntries(); iaod++){
    array = (TClonesArray*) list->At(iaod);
	if(deltaAODName!="") AddAODBranch("TClonesArray", &array, deltaAODName);//Put it in DeltaAOD file
	else AddAODBranch("TClonesArray", &array);//Put it in standard AOD file
  } 
	
  //Histograms container
  OpenFile(1);
  fOutputContainer = fAna->GetOutputContainer();
  fOutputContainer->SetOwner(kTRUE);
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::UserCreateOutputObjects() - End\n");
 
  PostData(1,fOutputContainer);
	
}
//_____________________________________________________
void AliAnalysisTaskParticleCorrelation::LocalInit()
{
	// Local Initialization
	
	//Call the Init to initialize the configuration of the analysis
	Init();
	
	// Create cuts/param objects and publish to slot
	fCuts = fAna->GetListOfAnalysisCuts();
	
	// Post Data
	PostData(2, fCuts);
	
}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelation::Init()
{
  // Initialization
 
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::Init() - Begin\n");
  
  // Call configuration file if specified
  
  if (fConfigName.Length()) {
    printf("AliAnalysisTaskParticleCorrelation::Init() - ### Configuration file is %s.C ###\n", fConfigName.Data());
	gROOT->LoadMacro(fConfigName+".C");
	fAna = (AliAnaPartCorrMaker*) gInterpreter->ProcessLine("ConfigAnalysis()");
  }
  
  if(!fAna) {
	printf("AliAnalysisTaskParticleCorrelation::Init() - Analysis maker pointer not initialized, no analysis specified, STOP !\n");
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

  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::Init() - End\n");
  
}


//_____________________________________________________
void AliAnalysisTaskParticleCorrelation::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::UserExec() - Begin\n");

   //Get the type of data, check if type is correct
  Int_t  datatype = fAna->GetReader()->GetDataType();
  if(datatype != AliCaloTrackReader::kESD && datatype != AliCaloTrackReader::kAOD &&
     datatype != AliCaloTrackReader::kMC){
    printf("AliAnalysisTaskParticleCorrelation::UserExec() - Wrong type of data\n");
    return ;
  }

  fAna->GetReader()->SetInputOutputMCEvent(InputEvent(), AODEvent(), MCEvent());

  //Process event
  fAna->ProcessEvent((Int_t) Entry(), CurrentFileName());
  //printf("AliAnalysisTaskParticleCorrelation::Current Event %d; Current File Name : %s\n",(Int_t) Entry(), CurrentFileName());
  if (DebugLevel() > 1) printf("AliAnalysisTaskParticleCorrelation::UserExec() - End\n");
	
  PostData(1, fOutputContainer);
	
  //gObjectTable->Print();

  
}

//_____________________________________________________
void AliAnalysisTaskParticleCorrelation::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  // Do some plots

  // Get merged histograms from the output container
  // Propagate histagrams to maker
  fAna->Terminate((TList*)GetOutputData(1));
	
}

