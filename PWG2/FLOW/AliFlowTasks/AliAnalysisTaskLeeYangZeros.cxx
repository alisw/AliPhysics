/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

#include "Riostream.h" //needed as include
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TProfile.h"

class AliAnalysisTaskSE;
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliFlowLYZConstants.h"   
#include "AliAnalysisTaskLeeYangZeros.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZHist1.h"
#include "AliFlowLYZHist2.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"

// AliAnalysisTaskLeeYangZeros:
// analysis task for Lee Yang Zeros method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLeeYangZeros)

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name, Bool_t firstrun) : 
  AliAnalysisTaskSE(name), 
  fEvent(0),
  fLyz(0),
  fFirstRunFile(0),
  fListHistos(NULL),
  fFirstRunLYZ(firstrun), //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE)       //set boolean for use sum to initial value
{
  // Constructor
  cout<<"AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, AliFlowEventSimple::Class());
  if (!firstrun) DefineInput(1, TList::Class()); //for second loop 
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());  
   
} 

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros() :
  AliAnalysisTaskSE(),
  fEvent(0),
  fLyz(0),
  fFirstRunFile(0),
  fListHistos(NULL),
  fFirstRunLYZ(kTRUE), //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE)    //set boolean for use sum to initial value
{
  // Constructor
  cout<<"AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros()"<<endl;

}

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::~AliAnalysisTaskLeeYangZeros()
{

  //destructor

}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::UserCreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskLeeYangZeros::CreateOutputObjects()"<<endl;

  
  //Analyser
  fLyz = new AliFlowAnalysisWithLeeYangZeros() ;
  fLyz -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fLyz -> SetUseSum(GetUseSumLYZ());       //set use sum true or false

  // Get data from input slot 1
  if (GetNinputs() == 2) {                   //if there are two input slots
    TList* pFirstRunList = (TList*)GetInputData(1);
    if (pFirstRunList) {
      fLyz -> SetFirstRunList(pFirstRunList);
    } else { cout<<"No first run List!"<<endl; exit(0); }
  }
  
  fLyz -> Init();

  if (fLyz->GetHistList()) {
    fListHistos = fLyz->GetHistList();
    //    fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }

 PostData(1,fListHistos);
  
}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent) {
    fLyz->Make(fEvent);
  }
  else {
    cout << "Warning no input data!!!" << endl; }
  
  PostData(1,fListHistos); 
  
}      

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Terminate(Option_t *) 
{
  // Called once at the end of the query
  
  AliFlowAnalysisWithLeeYangZeros* fLyzTerm = new AliFlowAnalysisWithLeeYangZeros() ;
  fLyzTerm -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fLyzTerm -> SetUseSum(GetUseSumLYZ());       //set use sum true or false
   
  fListHistos = (TList*)GetOutputData(1);
  
  if(fListHistos) 
  {
   fLyzTerm -> GetOutputHistograms(fListHistos);
   fLyzTerm -> Finish();
	PostData(1,fListHistos);
  } else 
    {
     cout << "histogram list pointer in Lee-Yang Zeros is empty in AliAnalysisTaskLYZ::Terminate ()" << endl;
    } 

  
}
