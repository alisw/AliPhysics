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

// root
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <Riostream.h>

// analysis
#include "AliAnalysisTaskGamma.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAnaGamma.h"
#include "AliGammaReader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliStack.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskGamma)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskGamma::AliAnalysisTaskGamma():
    fAna(0x0),
    fChain(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeG(0x0),
    fOutputContainer(0x0),
    fConfigName(0)
{
  // Default constructor
}

//_____________________________________________________
AliAnalysisTaskGamma::AliAnalysisTaskGamma(const char* name):
    AliAnalysisTask(name, "AnalysisTaskGamma"),
    fAna(0x0),
    fChain(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeG(0x0),
    fOutputContainer(0x0),
    fConfigName("ConfigGammaAnalysis")
{
  // Default constructor
 
  DefineInput (0, TChain::Class());
  DefineOutput(0, TTree::Class());
  DefineOutput(1, TList::Class());

}

//_____________________________________________________
AliAnalysisTaskGamma::~AliAnalysisTaskGamma() 
{
  // Remove all pointers
 
  if(fOutputContainer){
    fOutputContainer->Clear() ; 
    delete fOutputContainer ;
  }
  
  if(fTreeG) delete fTreeG ; 

}

//_____________________________________________________
void AliAnalysisTaskGamma::CreateOutputObjects()
{
  // Create the output container
  
  //AODs
//  OpenFile(0);
  AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  fAOD   = handler->GetAOD();
  fTreeG = handler->GetTree();
  fAna->ConnectAOD(fAOD);

  //Histograms container
  OpenFile(1);
  fOutputContainer = fAna->GetOutputContainer();
  
}

//_____________________________________________________
void AliAnalysisTaskGamma::Init()
{
  // Initialization
  AliDebug(1,"Begin");
  
  // Call configuration file

  if(fConfigName == ""){
    fConfigName="ConfigGammaAnalysis";
  }
 
  AliInfo(Form("### Configuration file is %s.C ###", fConfigName.Data()));
  gROOT->LoadMacro(fConfigName+".C");
  fAna = (AliAnaGamma*) gInterpreter->ProcessLine("ConfigGammaAnalysis()");
  
  if(!fAna)
    AliFatal("Analysis pointer not initialized, abort analysis!");
  
  // Initialise Gamma Analysis
  fAna->Init();
  
  //In case of MC analysis
/*
  // NOT ALLOWED TO SET MC HANDLER FROM WITHIN AN ANALYSIS TASK !!! (MG)
  Int_t  datatype = fAna->GetReader()->GetDataType();
  if(datatype == AliGammaReader::kMC || datatype == AliGammaReader::kMCData ){
    AliMCEventHandler * mc = new AliMCEventHandler();
    (AliAnalysisManager::GetAnalysisManager())->SetMCtruthEventHandler(mc);
  }
*/  
  AliDebug(1,"End");
  
}

//_____________________________________________________
void AliAnalysisTaskGamma::ConnectInputData(Option_t */*option*/)
{
  // Connect the input data
  //
  AliDebug(1,"ConnectInputData() ");
  AliESDInputHandler* esdH = (AliESDInputHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  fESD = esdH->GetEvent();
  fChain = (TChain*)GetInputData(0);
}

//_____________________________________________________
void AliAnalysisTaskGamma::Exec(Option_t */*option*/)
{
  // Execute analysis for current event
  //

  //Get the type of data, check if type is correct
  Int_t  datatype = fAna->GetReader()->GetDataType();
  if(datatype != AliGammaReader::kData && 
     datatype != AliGammaReader::kMC && 
     datatype != AliGammaReader::kMCData){
    AliFatal("Wrong type of data");
    return ;
  }

  //Get MC data
  AliStack* stack = 0x0; 
  if(datatype == AliGammaReader::kMC || datatype == AliGammaReader::kMCData ){
    AliMCEventHandler*    mctruth = (AliMCEventHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    
    if(mctruth)
      stack = mctruth->MCEvent()->Stack();
    
  }
  
  //Get Event
  Long64_t ientry = fChain->GetReadEntry();
  if ( !((ientry)%100) ) 
    AliInfo(Form("Analysing event # %5d\n", (Int_t) ientry));
  
  //Pass ESD pointer to analysis      
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  } 
  
  fAna->SetData(fESD);
  
  //In case of montecarlo analysis, pass the stack also.
  if((datatype == AliGammaReader::kMC || datatype == AliGammaReader::kMCData ) && stack)
    fAna -> SetKine(stack);
  
  //Process event
  fAna->ProcessEvent(ientry);
  
  PostData(0, fTreeG); 
  PostData(1, fOutputContainer);
  
}

//_____________________________________________________
void AliAnalysisTaskGamma::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(1,"Do nothing in Terminate");
 
}

