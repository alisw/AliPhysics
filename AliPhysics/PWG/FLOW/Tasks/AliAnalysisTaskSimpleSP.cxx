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
// AliAnalysisTaskSimpleSP:
//
// analysis task for Scalar Product Method
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"

#include "AliAnalysisTaskSimpleSP.h"
#include "AliFlowAnalysisWithSimpleSP.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "AliLog.h"

using std::endl;
using std::cout;
ClassImp(AliAnalysisTaskSimpleSP)

//________________________________________________________________________
AliAnalysisTaskSimpleSP::AliAnalysisTaskSimpleSP(const char *name, Bool_t usePhiWeights) : 
  AliAnalysisTaskSE(name), 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL),
  fMinimalBook(kFALSE),
  fUsePhiWeights(usePhiWeights),
  fListWeights(NULL),
  fRelDiffMsub(1.0),
  fApplyCorrectionForNUA(kFALSE),
  fHarmonic(2),
  fWeights(kFALSE),
  fScaling(kTRUE),
  fNormalizationType(1),
  fV0SanityCheck(0),
  fTotalQvector(NULL) 
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskSimpleSP::AliAnalysisTaskSimpleSP(const char *name)");

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  DefineInput(0, AliFlowEventSimple::Class());
  // Input slot #1 is needed for the weights input file
  if(usePhiWeights) {
    DefineInput(1, TList::Class()); }
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());

  // Total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
  fTotalQvector = new TString("QaQb");
}

//________________________________________________________________________
AliAnalysisTaskSimpleSP::AliAnalysisTaskSimpleSP() : 
  AliAnalysisTaskSE(), 
  fEvent(NULL),
  fSP(NULL),
  fListHistos(NULL),
  fMinimalBook(kFALSE),
  fUsePhiWeights(kFALSE),
  fListWeights(NULL),
  fRelDiffMsub(1.0),
  fApplyCorrectionForNUA(kFALSE),
  fHarmonic(2),
  fWeights(kFALSE),
  fScaling(kTRUE),
  fNormalizationType(1),
  fV0SanityCheck(0),
  fTotalQvector(NULL) 
  {
  // Constructor
    AliDebug(2,"AliAnalysisTaskSimpleSP::AliAnalysisTaskSimpleSP()");
}

//________________________________________________________________________
AliAnalysisTaskSimpleSP::~AliAnalysisTaskSimpleSP()
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
void AliAnalysisTaskSimpleSP::UserCreateOutputObjects() 
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskSimpleSP::CreateOutputObjects()");
  
  //Analyser
  fSP = new AliFlowAnalysisWithSimpleSP();
  fSP->SetBookOnlyBasicCCH(fMinimalBook);

  //set the allowed relative difference in the subevent multiplicities
  //fSP->SetRelDiffMsub(fRelDiffMsub); 
    
  //apply automatic correction for non-uniform acceptance:
  if (fApplyCorrectionForNUA) {
    AliDebug(2,"Corrections for non-uniform acceptance applied in the Scalar Product method");
  }
  fSP->SetApplyCorrectionForNUA(fApplyCorrectionForNUA);
  // harmonic: 
  fSP->SetHarmonic(fHarmonic);
  fSP->SetWeights(fWeights);
  fSP->SetScaling(fScaling);
  fSP->SetNormalizationType( fNormalizationType );
  fSP->SetV0SanityCheck(fV0SanityCheck);
  // total Q-vector:
  Int_t totalQ = 0;
  if( fTotalQvector->Contains("Qa") ) totalQ += 1;
  if( fTotalQvector->Contains("Qb") ) totalQ += 2;
  fSP->SetTotalQvector( totalQ );
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
  
  fSP-> Init();

  if (fSP->GetHistList()) {
    fListHistos = fSP->GetHistList();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }

  PostData(1,fListHistos);

}

//________________________________________________________________________
void AliAnalysisTaskSimpleSP::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent){
    fSP->ProcessEvent(fEvent);
  }
  else {
    cout << "Warning no input data for Scalar Product task!!!" << endl;
  }
    
  //fListHistos->Print();	
  PostData(1,fListHistos);
  
} 

//________________________________________________________________________
void AliAnalysisTaskSimpleSP::Terminate(Option_t *) 
{
  // Called once at the end of the query
  AliFlowAnalysisWithSimpleSP* fSPTerm = new AliFlowAnalysisWithSimpleSP() ;
  fListHistos = (TList*)GetOutputData(1);
  if (fListHistos) {
      fSPTerm -> GetOutputHistograms(fListHistos);
      fSPTerm -> Finish();
      PostData(1,fListHistos);
    }
    
  else { cout << "histgram list pointer is empty in Scalar Product" << endl; }

}
