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
#include "TFile.h" //needed as include
#include "TList.h"


class AliAnalysisTask;
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"

#include "AliAnalysisTaskScalarProduct.h"
#include "AliFlowAnalysisWithScalarProduct.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

// AliAnalysisTaskScalarProduct:
//
// analysis task for Scalar Product Method
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskScalarProduct)

//________________________________________________________________________
AliAnalysisTaskScalarProduct::AliAnalysisTaskScalarProduct(const char *name) : 
  AliAnalysisTask(name, ""), 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL)
  {
  // Constructor
  cout<<"AliAnalysisTaskScalarProduct::AliAnalysisTaskScalarProduct(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  DefineInput(0, AliFlowEventSimple::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
   
}

//________________________________________________________________________
AliAnalysisTaskScalarProduct::AliAnalysisTaskScalarProduct() : 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL)
  {
  // Constructor
  cout<<"AliAnalysisTaskScalarProduct::AliAnalysisTaskScalarProduct()"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskScalarProduct::~AliAnalysisTaskScalarProduct()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  //  if (ListHistos) {
  //    delete fListHistos;
  //    fListHistos = NULL;
  //  }
}

//________________________________________________________________________
void AliAnalysisTaskScalarProduct::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskScalarProduct::ConnectInputData(Option_t *)"<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskScalarProduct::CreateOutputObjects() 
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskScalarProduct::CreateOutputObjects()"<<endl;
  
  //Analyser
  fSP  = new AliFlowAnalysisWithScalarProduct() ;
  fSP-> Init();
  

  if (fSP->GetHistList()) {
    //fSP->GetHistList()->Print();
    fListHistos = fSP->GetHistList();
    //fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }
}

//________________________________________________________________________
void AliAnalysisTaskScalarProduct::Exec(Option_t *) 
{
  // Main loop
  // Called for each event


  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent){
    fSP->Make(fEvent);
  }
  else {
    cout << "Warning no input data!!!" << endl;
  }
    
  //fListHistos->Print();	
  PostData(0,fListHistos);
  
} 

//________________________________________________________________________
void AliAnalysisTaskScalarProduct::Terminate(Option_t *) 
{
  // Called once at the end of the query -- do not call in case of CAF
  //  fSP->Finish();
  //  PostData(0,fListHistos);

  fListHistos = (TList*)GetOutputData(0);
  // cout << "histgram list in Terminate" << endl;
  if (fListHistos)  {
    //    fListHistos->Print();
  }	
  else { cout << "histgram list pointer is empty" << endl; }

}
