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
 * analysis task for Q-cumulants      * 
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
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TBits.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskQCumulants.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

ClassImp(AliAnalysisTaskQCumulants)

//================================================================================================================

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fEvent(NULL),
 fQCA(NULL), // Q-cumulant Analysis (QCA) object
 fListHistos(NULL),
 fUseWeights(useWeights),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name)"<<endl;
 
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

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(): 
 fEvent(NULL),
 fQCA(NULL),//Q-cumulant Analysis (QCA) object
 fListHistos(NULL),
 fUseWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fListWeights(NULL)
{
 // dummy constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::ConnectInputData(Option_t *) 
{
 // connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskQCumulants::ConnectInputData(Option_t *)"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::CreateOutputObjects() 
{
 // called at every worker node to initialize
 cout<<"AliAnalysisTaskQCumulants::CreateOutputObjects()"<<endl;

 // analyser
 fQCA = new AliFlowAnalysisWithQCumulants();
 fQCA->Init();
 
 //weights:
 if(fUseWeights)
 {
  //pass the flags to class:
  if(fUsePhiWeights) fQCA->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fQCA->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fQCA->SetUseEtaWeights(fUseEtaWeights);
  //get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fListWeights = (TList*)GetInputData(1); 
  }
  //pass the list with weights to class:
  if(fListWeights) fQCA->SetWeightsList(fListWeights);
 }
 
 if(fQCA->GetHistList()) 
 {
  fListHistos = fQCA->GetHistList();
  //fListHistos->Print();
 }
 else 
 {
  Printf(" ERROR: Could not retrieve histogram list (QC, Task::COO)"); 
 }

 //PostData(0,fListHistos);
 
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Exec(Option_t *) 
{
  // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Q-cumulants
 if(fEvent) 
 {
  fQCA->Make(fEvent);
 }else 
  {
   cout<<" WARNING: No input data (QC, Task::E) !!!"<<endl;
   cout<<endl;
  }
  
 PostData(0,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Terminate(Option_t *) 
{
 //accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();
 
 fQCA = new AliFlowAnalysisWithQCumulants(); 
 
 if (fListHistos) 
 {
  fQCA->GetOutputHistograms(fListHistos);
  fQCA->Finish();
 }
 else
 {
  cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate())"<<endl;
  cout<<endl;
 }
}





















