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
    AliAnalysisTaskSE(),
    fAna(0x0),
    fOutputContainer(0x0),
    fConfigName(0)
{
  // Default constructor
}

//_____________________________________________________
AliAnalysisTaskGamma::AliAnalysisTaskGamma(const char* name):
    AliAnalysisTaskSE(name),
    fAna(0x0),
    fOutputContainer(0x0),
    fConfigName("ConfigGammaAnalysis")
{
  // Default constructor

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

}

//_____________________________________________________
void AliAnalysisTaskGamma::UserCreateOutputObjects()
{
  // Create the output container
  
  //AODs

  fAna->ConnectAOD(AODEvent());

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
  
  AliDebug(1,"End");
  
}


//_____________________________________________________
void AliAnalysisTaskGamma::UserExec(Option_t */*option*/)
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
  
  fAna->SetData(InputEvent());
  
  //In case of montecarlo analysis, pass the stack also.
  if((datatype == AliGammaReader::kMC || datatype == AliGammaReader::kMCData ) && MCEvent())
    fAna -> SetKine(MCEvent()->Stack());
  
  //Process event
  fAna->ProcessEvent();
  
  PostData(1, fOutputContainer);
  
}

//_____________________________________________________
void AliAnalysisTaskGamma::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(1,"Do nothing in Terminate");
  //fAna->Terminate();
}

