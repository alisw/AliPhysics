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

///////////////////////////////////////////////
// AliAnalysisTaskMSP:
//
// analysis task for Scalar Product Method
//
// Author: Paul Kuijer (Paul.Kuijer@nikhef.nl)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"

#include "AliAnalysisTaskMSP.h"
#include "AliFlowAnalysisWithMSP.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "AliLog.h"

using std::endl;
using std::cout;
ClassImp(AliAnalysisTaskMSP)

//________________________________________________________________________
AliAnalysisTaskMSP::AliAnalysisTaskMSP(const char *name, Bool_t usePhiWeights) : 
  AliAnalysisTaskSE(name), 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL),
  fMinimalBook(kFALSE),
  fUsePhiWeights(usePhiWeights),
  fListWeights(NULL),
  fApplyCorrectionForNUA(kFALSE),
  fHarmonic(2)
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskMSP::AliAnalysisTaskMSP(const char *name)");

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  DefineInput(0, AliFlowEventSimple::Class());
  // Input slot #1 is needed for the weights input file
  if(usePhiWeights) {
    DefineInput(1, TList::Class()); }
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMSP::AliAnalysisTaskMSP() : 
  AliAnalysisTaskSE(), 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL),
  fMinimalBook(kFALSE),
  fUsePhiWeights(kFALSE),
  fListWeights(NULL),
  fApplyCorrectionForNUA(kFALSE),
  fHarmonic(0)
  {
  // Constructor
    AliDebug(2,"AliAnalysisTaskMSP::AliAnalysisTaskMSP()");
}

//________________________________________________________________________
AliAnalysisTaskMSP::~AliAnalysisTaskMSP()
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
void AliAnalysisTaskMSP::UserCreateOutputObjects() 
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskMSP::CreateOutputObjects()");
  
  //Analyser
  fSP = new AliFlowAnalysisWithMSP();
  fSP->UseCommonConstants(kTRUE);
  fSP->EnableCommonHistograms(!fMinimalBook);

  //set the allowed relative difference in the subevent multiplicities
  //fSP->SetRelDiffMsub(fRelDiffMsub); 
    
  //apply automatic correction for non-uniform acceptance:
  if (fApplyCorrectionForNUA) {
    AliDebug(2,"Corrections for non-uniform acceptance applied in the MSP method");
  }
  fSP->EnableNUA(fApplyCorrectionForNUA);
  // harmonic: 
  fSP->SetHarmonic(fHarmonic);

  /* TODO: Phi weigths not implemented yet:
  //for using phi weights:
  if(fUsePhiWeights) {
    //pass the flag to the analysis class:
    fSP->SetUsePhiWeights(fUsePhiWeights);
    //get data from input slot #1 which is used for weights:
    if(GetNinputs()==2) {                   
      fListWeights = (TList*)GetInputData(1); 
    }
    //pass the list with weights to the analysis class:
    if(fListWeights) fSP->SetWeightsList(fListWeights);
  }
  */
  
  fSP-> Init();

  if( !(fListHistos = fSP->ListHistograms()) ) {
     Printf("ERROR: Could not retrieve histogram list"); 
  }else{
     fListHistos->SetOwner(kTRUE); // TODO: required by analysis manager but will clash with destructor of AliAnalysisWithMSP
  }

  PostData(1,fListHistos);
}

//________________________________________________________________________
void AliAnalysisTaskMSP::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent){
    fSP->Make(fEvent);
  }
  else {
    cout << "Warning no input data for Scalar Product task!!!" << endl;
  }
    
  //fListHistos->Print();	
  PostData(1,fListHistos);
  
} 

//________________________________________________________________________
void AliAnalysisTaskMSP::Terminate(Option_t *) 
{
  // Called once at the end of the query
  // TODO: check if it is really required to do Finish here. The results should be calculated after merging, not here. Finish does not do much else.
  
   fListHistos = (TList*)GetOutputData(1);
   if (fListHistos) {
      // Create a temporary analysis object to calculate the results because the original is not available on CAF.
      AliFlowAnalysisWithMSP* fMSPTerm = new AliFlowAnalysisWithMSP(fListHistos) ;   // Load histograms
      fMSPTerm->EnableNUA(fApplyCorrectionForNUA);
      fMSPTerm -> Finish();                       // Create result histograms
      fListHistos=fMSPTerm->ListHistograms();     // Get new list including results
      fListHistos->SetOwner(kTRUE); // TODO: Required by analysis manager, but will clash with destructor of this object
      PostData(1,fListHistos);
   } else { 
      cout << "histgram list pointer is empty in Scalar Product" << endl; 
   }
}
