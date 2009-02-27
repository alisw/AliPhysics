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
#include "TProfile.h"
#include "TFile.h"
#include "TList.h"

class AliAnalysisTask;
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskLYZEventPlane.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZEventPlane.h"
#include "AliFlowAnalysisWithLYZEventPlane.h"

// AliAnalysisTaskLYZEventPlane:
//
// analysis task for Lee Yang Zeros Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLYZEventPlane)

//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name) : 
  AliAnalysisTask(name, ""), 
  fEvent(NULL), 
  fLyzEp(NULL),
  fLyz(NULL),
  fListHistos(NULL),
  fSecondRunFile(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  DefineInput(0, AliFlowEventSimple::Class());
  DefineInput(1, TList::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane() : 
  fEvent(NULL), 
  fLyzEp(NULL),
  fLyz(NULL),
  fListHistos(NULL),
  fSecondRunFile(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane()"<<endl;
}


//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::~AliAnalysisTaskLYZEventPlane() 
{
  //destructor

}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskLYZEventPlane::ConnectInputData(Option_t *)"<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskLYZEventPlane::CreateOutputObjects()"<<endl;
  
  //lee yang zeros event plane
  fLyzEp = new AliFlowLYZEventPlane() ;
  //Analyser
  fLyz = new AliFlowAnalysisWithLYZEventPlane() ;
     
  // Get data from input slot
  TList* pSecondRunList = (TList*)GetInputData(1);
  if (pSecondRunList) {
    fLyzEp -> SetSecondRunList(pSecondRunList);
    fLyz -> SetSecondRunList(pSecondRunList);
  } else { cout<<"No Second run List!"<<endl; exit(0); }

  fLyzEp-> Init();
  fLyz-> Init();

  if (fLyz->GetHistList()) {
    fListHistos = fLyz->GetHistList();
    //fListHistos->Print();
  }
  else { cout<<"ERROR: Could not retrieve histogram list"<<endl;}


}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent) {
    fLyz->Make(fEvent,fLyzEp);
  }
  else {
    cout << "Warning no input data!!!" << endl;}
    
  PostData(0,fListHistos);
  
}      

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  AliFlowAnalysisWithLYZEventPlane* fLyzTerm = new AliFlowAnalysisWithLYZEventPlane() ;
  fListHistos = (TList*)GetOutputData(0);
  //cout << "histogram list in Terminate" << endl;
  if (fListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHist = dynamic_cast<AliFlowCommonHist*> 
      (fListHistos->FindObject("AliFlowCommonHistLYZEP"));
    AliFlowCommonHistResults *pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
      (fListHistos->FindObject("AliFlowCommonHistResultsLYZEP"));

    TProfile* pHistProR0theta = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("First_FlowPro_r0theta_LYZ"));

    TProfile* pHistProVetaRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VetaRP_LYZEP"));
    TProfile* pHistProVetaPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VetaPOI_LYZEP"));
    TProfile* pHistProVPtRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPtRP_LYZEP"));
    TProfile* pHistProVPtPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPtPOI_LYZEP"));

    TH1F* pHistQsumforChi = dynamic_cast<TH1F*> 
      (fListHistos->FindObject("Flow_QsumforChi_LYZEP"));

    if (pCommonHist && pCommonHistResults && pHistProR0theta &&
	pHistProVetaRP && pHistProVetaPOI && pHistProVPtRP && 
	pHistProVPtPOI && pHistQsumforChi ) {
      fLyzTerm -> SetCommonHists(pCommonHist);
      fLyzTerm -> SetCommonHistsRes(pCommonHistResults);
      fLyzTerm -> SetFirstr0theta(pHistProR0theta);
      fLyzTerm -> SetHistProVetaRP(pHistProVetaRP);
      fLyzTerm -> SetHistProVetaPOI(pHistProVetaPOI);
      fLyzTerm -> SetHistProVPtRP(pHistProVPtRP);
      fLyzTerm -> SetHistProVPtPOI(pHistProVPtPOI);
      fLyzTerm -> SetHistQsumforChi(pHistQsumforChi);
      fLyzTerm -> Finish();
      PostData(0,fListHistos);
    } else { 
      cout<<"WARNING: Histograms needed to run Finish() are not accessable!"<<endl; 
    }

    //fListHistos->Print(); 
  } else { cout << "histogram list pointer is empty" << endl;}

  //cout<<".....finished LYZ EventPlane"<<endl;  
}


