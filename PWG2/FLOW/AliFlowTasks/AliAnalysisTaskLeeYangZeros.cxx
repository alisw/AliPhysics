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

class AliAnalysisTask;
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
  AliAnalysisTask(name, ""), 
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
  DefineOutput(0, TList::Class());  
   
} 

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros() :  
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
void AliAnalysisTaskLeeYangZeros::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskLeeYangZeros::ConnectInputData(Option_t *)"<<endl;
 
}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::CreateOutputObjects() 
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
  
}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  if (fEvent) {
    fLyz->Make(fEvent);
  }
  else {
    cout << "Warning no input data!!!" << endl; }
  
  PostData(0,fListHistos); //here for CAF
  
}      

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Terminate(Option_t *) 
{
  // Called once at the end of the query

  const Int_t iNtheta = AliFlowLYZConstants::kTheta;

  AliFlowAnalysisWithLeeYangZeros* fLyzTerm = new AliFlowAnalysisWithLeeYangZeros() ;
  fLyzTerm -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fLyzTerm -> SetUseSum(GetUseSumLYZ());       //set use sum true or false
   
  fListHistos = (TList*)GetOutputData(0);
  //cout << "histogram list in Terminate" << endl;

  if (fListHistos) {

    //define histograms for first and second run
    AliFlowCommonHist *pCommonHist = NULL;
    AliFlowCommonHistResults *pCommonHistResults = NULL;
    TProfile* pHistProVtheta = NULL;
    TProfile* pHistProReDenom = NULL;
    TProfile* pHistProImDenom = NULL;
    TProfile* pHistProReDtheta = NULL;
    TProfile* pHistProImDtheta = NULL;
    TProfile* pHistProVetaRP = NULL;
    TProfile* pHistProVetaPOI = NULL;
    TProfile* pHistProVPtRP  = NULL;
    TProfile* pHistProVPtPOI  = NULL;
    AliFlowLYZHist1 *pLYZHist1[iNtheta] = {NULL};      //array of pointers to AliFlowLYZHist1
    AliFlowLYZHist2 *pLYZHist2RP[iNtheta] = {NULL};    //array of pointers to AliFlowLYZHist2
    AliFlowLYZHist2 *pLYZHist2POI[iNtheta] = {NULL};   //array of pointers to AliFlowLYZHist2

    if (GetFirstRunLYZ()) { //first run
      //Get the common histograms from the output list
      pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	(fListHistos->FindObject("AliFlowCommonHistLYZ1"));
      pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	(fListHistos->FindObject("AliFlowCommonHistResultsLYZ1"));
    }
    else { //second run
      //Get the common histograms from the output list
      pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	(fListHistos->FindObject("AliFlowCommonHistLYZ2"));
      pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	(fListHistos->FindObject("AliFlowCommonHistResultsLYZ2"));
    }

    TProfile* pHistProR0theta = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("First_FlowPro_r0theta_LYZ"));

    TH1F* pHistQsumforChi = dynamic_cast<TH1F*> 
      (fListHistos->FindObject("Flow_QsumforChi_LYZ"));

    
    if (GetFirstRunLYZ()) { //for firstrun
      //Get the histograms from the output list
      for(Int_t theta = 0;theta<iNtheta;theta++){
	TString name = "AliFlowLYZHist1_"; 
	name += theta;
	pLYZHist1[theta] = dynamic_cast<AliFlowLYZHist1*> 
	  (fListHistos->FindObject(name));
      }
      pHistProVtheta = dynamic_cast<TProfile*> 
	  (fListHistos->FindObject("First_FlowPro_Vtheta_LYZ"));

      //Set the histogram pointers and call Finish()
      if (pCommonHist && pCommonHistResults && pLYZHist1[0] && 
	  pHistProVtheta && pHistProR0theta && pHistQsumforChi ) {
	fLyzTerm->SetCommonHists(pCommonHist);
	fLyzTerm->SetCommonHistsRes(pCommonHistResults);
	fLyzTerm->SetHist1(pLYZHist1);
	fLyzTerm->SetHistProVtheta(pHistProVtheta);
	fLyzTerm->SetHistProR0theta(pHistProR0theta);
	fLyzTerm->SetHistQsumforChi(pHistQsumforChi);
	fLyzTerm->Finish();
	PostData(0,fListHistos);
      } else { 
	cout<<"WARNING: Histograms needed to run Finish() firstrun are not accessable!"<<endl; 
      }
    } else { //for second run
      //Get the histograms from the output list
      for(Int_t theta = 0;theta<iNtheta;theta++){
	TString nameRP = "AliFlowLYZHist2RP_"; 
	nameRP += theta;
	pLYZHist2RP[theta] = dynamic_cast<AliFlowLYZHist2*> 
	  (fListHistos->FindObject(nameRP));
	TString namePOI = "AliFlowLYZHist2POI_"; 
	namePOI += theta;
	pLYZHist2POI[theta] = dynamic_cast<AliFlowLYZHist2*> 
	  (fListHistos->FindObject(namePOI));

      }
      
      pHistProReDenom = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ReDenom_LYZ"));
      pHistProImDenom = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ImDenom_LYZ"));

      pHistProReDtheta = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ReDtheta_LYZ"));
      pHistProImDtheta = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ImDtheta_LYZ"));

      pHistProVetaRP = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_VetaRP_LYZ"));
      pHistProVetaPOI = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_VetaPOI_LYZ"));
      pHistProVPtRP = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_VPtRP_LYZ"));
      pHistProVPtPOI = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_VPtPOI_LYZ"));


      //Set the histogram pointers and call Finish()
      if (pCommonHist && pCommonHistResults && pLYZHist2RP[0] && pLYZHist2POI[0] && 
	  pHistProR0theta && pHistProReDenom && pHistProImDenom && pHistProVetaRP && 
	  pHistProVetaPOI && pHistProVPtRP && pHistProVPtPOI) {
	fLyzTerm->SetCommonHists(pCommonHist);
	fLyzTerm->SetCommonHistsRes(pCommonHistResults);
	fLyzTerm->SetHist2RP(pLYZHist2RP);
	fLyzTerm->SetHist2POI(pLYZHist2POI);
	fLyzTerm->SetHistProR0theta(pHistProR0theta);
	fLyzTerm->SetHistProReDenom(pHistProReDenom);
	fLyzTerm->SetHistProImDenom(pHistProImDenom);
	fLyzTerm->SetHistProReDtheta(pHistProReDtheta);
	fLyzTerm->SetHistProImDtheta(pHistProImDtheta);
	fLyzTerm->SetHistProVetaRP(pHistProVetaRP);
	fLyzTerm->SetHistProVetaPOI(pHistProVetaPOI);
	fLyzTerm->SetHistProVPtRP(pHistProVPtRP);
	fLyzTerm->SetHistProVPtPOI(pHistProVPtPOI);
	fLyzTerm->SetHistQsumforChi(pHistQsumforChi);
	fLyzTerm->Finish();
	PostData(0,fListHistos);
      } else { 
	cout<<"WARNING: Histograms needed to run Finish() secondrun are not accessable!"<<endl; 
      }
    }
          
    //    fListHistos->Print(); 
  }	
  else { cout << "histogram list pointer in Lee-Yang Zeros is empty" << endl;}

  //cout<<".....finished"<<endl;
}
