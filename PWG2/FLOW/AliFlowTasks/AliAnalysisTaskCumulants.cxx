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

/**************************************
 * analysis task for cumulant method  * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/
 
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCumulants.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliCumulantsFunctions.h"

ClassImp(AliAnalysisTaskCumulants)

//================================================================================================================

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fEvent(NULL),
 fGFCA(NULL), // Generating Function Cumulant (GFCA) analysis object
 fListHistos(NULL),
 fUseWeights(useWeights),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with a TChain
 DefineInput(0, AliFlowEventSimple::Class());
 
 // Input slot #1 is needed for the weights 
 if(useWeights)
 {
  DefineInput(1, TList::Class());   
 }
  
 // Output slot #0 writes into a TList container
 DefineOutput(0, TList::Class());   
}

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants():
 fEvent(NULL),
 fGFCA(NULL), // Generating Function Cumulant (GFCA) analysis object
 fListHistos(NULL),
 fUseWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // dummy constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskCumulants::ConnectInputData(Option_t *) 
{
 // connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskCumulants::ConnectInputData(Option_t *)"<<endl;
}

//================================================================================================================

void AliAnalysisTaskCumulants::CreateOutputObjects() 
{
 // called at every worker node to initialize
 cout<<"AliAnalysisTaskCumulants::CreateOutputObjects()"<<endl;

 // analyser
 fGFCA = new AliFlowAnalysisWithCumulants();
 fGFCA->Init();
 
 //weights:
 if(fUseWeights)
 {
  //pass the flags to class:
  if(fUsePhiWeights) fGFCA->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fGFCA->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fGFCA->SetUseEtaWeights(fUseEtaWeights);
  //get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fListWeights = (TList*)GetInputData(1); 
  }
  //pass the list with weights to class:
  if(fListWeights) fGFCA->SetWeightsList(fListWeights);
 }

 if(fGFCA->GetHistList()) 
 {
  fListHistos = fGFCA->GetHistList();
  //fListHistos->Print();
 }
 else
 {
  Printf(" ERROR: Could not retrieve histogram list (GFCA, Task::COO)"); 
 }
}

//================================================================================================================

void AliAnalysisTaskCumulants::Exec(Option_t *) 
{
 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // generating function cumulants
 if(fEvent) 
 {
  fGFCA->Make(fEvent);
 }else 
  {
   cout<<" WARNING: No input data (GFCA, Task::E) !!!"<<endl;
   cout<<endl;
  }
  
 PostData(0,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskCumulants::Terminate(Option_t *) 
{  
 //accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();
 
 fGFCA = new AliFlowAnalysisWithCumulants();  
 
 if(fListHistos)
 {
  fGFCA->GetOutputHistograms(fListHistos);
  fGFCA->Finish();
 } else
   {
    cout<<" WARNING: histogram list pointer is empty (GFC, Task::T)"<<endl;
    cout<<endl;
   }
}





















