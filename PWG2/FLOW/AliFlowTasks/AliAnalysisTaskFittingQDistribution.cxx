/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *f
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

/**************************************
 *    analysis task for fitting       * 
 *         q-distribution             *
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

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskFittingQDistribution.h"
#include "AliFittingQDistribution.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHistResults.h"
#include "AliFittingFunctionsForQDistribution.h"

ClassImp(AliAnalysisTaskFittingQDistribution)

//================================================================================================================

AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(const char *name, Bool_t useWeights): 
  AliAnalysisTask(name,""), 
  fEvent(NULL),
  fFQDA(NULL),//Fitting Q_Distribution Analysis (FQDA) object
  fListHistos(NULL),
  fUseWeights(useWeights),
  fUsePhiWeights(kFALSE),
  fListWeights(NULL)
{
  //constructor
  cout<<"AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(const char *name)"<<endl;
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, AliFlowEventSimple::Class());
  
  // Input slot #1 is needed for the weights 
  if(useWeights) {
    DefineInput(1, TList::Class());   
  }
  
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(): 
  fEvent(NULL),
  fFQDA(NULL),//Fitting q-distribution Analysis (FQDA) object
  fListHistos(NULL),  
  fUseWeights(kFALSE),
  fUsePhiWeights(kFALSE),
  fListWeights(NULL)
{
 //dummy constructor
 cout<<"AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::ConnectInputData(Option_t *) 
{
 //connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskFittingQDistribution::ConnectInputData(Option_t *)"<<endl;

}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::CreateOutputObjects() 
{
  //called at every worker node to initialize
  cout<<"AliAnalysisTaskFittingQDistribution::CreateOutputObjects()"<<endl;
  
  
  //analyser
  fFQDA = new AliFittingQDistribution();
  fFQDA->Init();
  
  //weights:
  if(fUseWeights) {
    //pass the flags to class:
    if(fUsePhiWeights) fFQDA->SetUsePhiWeights(fUsePhiWeights);
    //get data from input slot #1 which is used for weights:
    if(GetNinputs()==2) {                   
      fListWeights = (TList*)GetInputData(1); 
    }
    //pass the list with weights to class:
    if(fListWeights) fFQDA->SetWeightsList(fListWeights);
  }
  
  if(fFQDA->GetHistList()) {
    fListHistos = fFQDA->GetHistList();
    //fListHistos->Print();
  }
  else {
    Printf("ERROR: Could not retrieve histogram list"); 
  }
  
}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::Exec(Option_t *) 
{
  //main loop (called for each event)
  fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

  //fitting q-distribution 
  if (fEvent) {
    fFQDA->Make(fEvent);
  }
  else {
    cout << "Warning no input data!!!" << endl;
  }
  PostData(0,fListHistos); 
}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::Terminate(Option_t *) 
{  
  //accessing the output list
  fListHistos = (TList*)GetOutputData(0);
  //fListHistos->Print();
  
  fFQDA = new AliFittingQDistribution();
  
  if(fListHistos) 
  {	     
   fFQDA->GetOutputHistograms(fListHistos);
   fFQDA->Finish();  
  } else 
    {
     cout<<"histogram list pointer is empty"<<endl;
    }
}

//================================================================================================================



















