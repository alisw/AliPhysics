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

class TFile;
class TList;
class AliAnalysisTask;

#include "TProfile.h"  //needed as include
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
  // Called once at the end of the query
  AliFlowAnalysisWithScalarProduct* fSPTerm = new AliFlowAnalysisWithScalarProduct() ;
  fListHistos = (TList*)GetOutputData(0);
  if (fListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHist = dynamic_cast<AliFlowCommonHist*> 
      (fListHistos->FindObject("AliFlowCommonHistSP"));
    AliFlowCommonHistResults *pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
      (fListHistos->FindObject("AliFlowCommonHistResultsSP"));
    TProfile* pHistProQaQb     = dynamic_cast<TProfile*>(fListHistos->FindObject("Flow_QaQb_SP"));
    TProfile* pHistProUQetaRP  = dynamic_cast<TProfile*>(fListHistos->FindObject("Flow_UQetaRP_SP"));
    TProfile* pHistProUQetaPOI = dynamic_cast<TProfile*>(fListHistos->FindObject("Flow_UQetaPOI_SP"));
    TProfile* pHistProUQPtRP   = dynamic_cast<TProfile*>(fListHistos->FindObject("Flow_UQPtRP_SP"));
    TProfile* pHistProUQPtPOI  = dynamic_cast<TProfile*>(fListHistos->FindObject("Flow_UQPtPOI_SP"));
    if (pCommonHist && pCommonHistResults && pHistProQaQb && pHistProUQetaRP
	&& pHistProUQetaPOI && pHistProUQPtRP && pHistProUQPtPOI) {
      fSPTerm -> SetCommonHists(pCommonHist);
      fSPTerm -> SetCommonHistsRes(pCommonHistResults);
      fSPTerm -> SetHistProQaQb(pHistProQaQb);
      fSPTerm -> SetHistProUQetaRP(pHistProUQetaRP);
      fSPTerm -> SetHistProUQetaPOI(pHistProUQetaPOI);
      fSPTerm -> SetHistProUQPtRP(pHistProUQPtRP);
      fSPTerm -> SetHistProUQPtPOI(pHistProUQPtPOI);
      fSPTerm -> Finish();
      PostData(0,fListHistos);
    }
  }
  else { cout << "histgram list pointer is empty in Scalar Product" << endl; }

}
