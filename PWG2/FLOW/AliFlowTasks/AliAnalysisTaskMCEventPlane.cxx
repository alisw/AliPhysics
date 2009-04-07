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
#include "TProfile.h"
#include "TProfile2D.h"
#include "TList.h"


class AliAnalysisTask;
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"


#include "AliAnalysisTaskMCEventPlane.h"
#include "AliFlowAnalysisWithMCEventPlane.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

// AliAnalysisTaskMCEventPlane:
//
// analysis task for Monte Carlo Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskMCEventPlane)

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name) : 
  AliAnalysisTask(name, ""), 
  fEvent(NULL),
  fMc(NULL),
  fListHistos(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, AliFlowEventSimple::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane() : 
  fEvent(NULL),
  fMc(NULL),
  fListHistos(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane()"<<endl;

}

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::~AliAnalysisTaskMCEventPlane()
{

  //destructor

}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskMCEventPlane::ConnectInputData(Option_t *)"<<endl;

}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskMCEventPlane::CreateOutputObjects()"<<endl;

  //Analyser
  fMc  = new AliFlowAnalysisWithMCEventPlane() ;
      
  fMc-> Init();

  if (fMc->GetHistList()) {
    //fMc->GetHistList()->Print();
    fListHistos = fMc->GetHistList();
    //fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }

}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent){
    fMc->Make(fEvent);
  }
  else {
    cout << "Warning no input data!!!" << endl;
  }

  PostData(0,fListHistos); 
}      


//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  AliFlowAnalysisWithMCEventPlane* fMcTerm = new AliFlowAnalysisWithMCEventPlane() ;

  //Get output data
  fListHistos = (TList*)GetOutputData(0);
  // cout << "histogram list in Terminate" << endl;
  if (fListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHists = dynamic_cast<AliFlowCommonHist*> 
      (fListHistos->FindObject("AliFlowCommonHistMCEP"));
    AliFlowCommonHistResults *pCommonHistResults = 
      dynamic_cast<AliFlowCommonHistResults*> 
      (fListHistos->FindObject("AliFlowCommonHistResultsMCEP"));

    TProfile *pHistProIntFlow = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_V_MCEP")); 
      
    TProfile2D *pHistProDiffFlowPtEtaRP = dynamic_cast<TProfile2D*> 
      (fListHistos->FindObject("FlowPro_VPtEtaRP_MCEP")); 
                               
    TProfile *pHistProDiffFlowPtRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPtRP_MCEP")); 
     
    TProfile *pHistProDiffFlowEtaRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VetaRP_MCEP"));
 
    TProfile2D *pHistProDiffFlowPtEtaPOI = dynamic_cast<TProfile2D*> 
      (fListHistos->FindObject("FlowPro_VPtEtaPOI_MCEP")); 
          
    TProfile *pHistProDiffFlowPtPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPtPOI_MCEP")); 
     
    TProfile *pHistProDiffFlowEtaPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VetaPOI_MCEP"));                             

    if (pCommonHists && pCommonHistResults && pHistProIntFlow && 
	pHistProDiffFlowPtRP && pHistProDiffFlowEtaRP && 
	pHistProDiffFlowPtPOI && pHistProDiffFlowEtaPOI) {
      fMcTerm->SetCommonHists(pCommonHists);
      fMcTerm->SetCommonHistsRes(pCommonHistResults);
      fMcTerm->SetHistProIntFlow(pHistProIntFlow);
      fMcTerm->SetHistProDiffFlowPtEtaRP(pHistProDiffFlowPtEtaRP);
      fMcTerm->SetHistProDiffFlowPtRP(pHistProDiffFlowPtRP);      
      fMcTerm->SetHistProDiffFlowEtaRP(pHistProDiffFlowEtaRP);  
      fMcTerm->SetHistProDiffFlowPtEtaPOI(pHistProDiffFlowPtEtaPOI);
      fMcTerm->SetHistProDiffFlowPtPOI(pHistProDiffFlowPtPOI);      
      fMcTerm->SetHistProDiffFlowEtaPOI(pHistProDiffFlowEtaPOI);          
      fMcTerm->Finish();
      PostData(0,fListHistos);
    } else {
      cout<<"WARNING: Histograms needed to run Finish() are not accessable!"<<endl;  }
    
    //fListHistos->Print();
  } else { cout << "histogram list pointer is empty" << endl;}
}

