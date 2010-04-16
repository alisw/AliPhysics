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

/**********************************
 * analysis task for nested loops * 
 *                                * 
 * authors: Naomi van der Kolk    *
 *           (kolk@nikhef.nl)     *  
 *          Raimond Snellings     *
 *           (snelling@nikhef.nl) * 
 *          Ante Bilandzic        *
 *           (anteb@nikhef.nl)    * 
 * *******************************/
 
class TFile;
class TList;
class AliAnalysisTaskSE; 
 
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskNestedLoops.h"
#include "AliFlowAnalysisWithNestedLoops.h"

ClassImp(AliAnalysisTaskNestedLoops)

//================================================================================================================

AliAnalysisTaskNestedLoops::AliAnalysisTaskNestedLoops(const char *name, Bool_t useParticleWeights): 
AliAnalysisTaskSE(name), 
fEvent(NULL),
fNL(NULL), 
fListHistos(NULL),
fUseParticleWeights(useParticleWeights),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskNestedLoops::AliAnalysisTaskNestedLoops(const char *name, Bool_t useParticleWeights)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with an AliFlowEventSimple
 DefineInput(0, AliFlowEventSimple::Class());  
 // Input slot #1 is needed for the weights input file:
 if(useParticleWeights)
 {
  DefineInput(1, TList::Class());   
 }  
 // Output slot #0 is reserved              
 // Output slot #1 writes into a TList container
 DefineOutput(1, TList::Class());  
}

AliAnalysisTaskNestedLoops::AliAnalysisTaskNestedLoops(): 
AliAnalysisTaskSE(),
fEvent(NULL),
fNL(NULL),
fListHistos(NULL),
fUseParticleWeights(kFALSE),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL)
{
 // Dummy constructor
 cout<<"AliAnalysisTaskNestedLoops::AliAnalysisTaskNestedLoops()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskNestedLoops::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize
 cout<<"AliAnalysisTaskNestedLoops::UserCreateOutputObjects()"<<endl;

 // Analyser:
 fNL = new AliFlowAnalysisWithNestedLoops();

 if(fUseParticleWeights)
 {
  // Pass the flags to class:
  if(fUsePhiWeights) fNL->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fNL->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fNL->SetUseEtaWeights(fUseEtaWeights);
  // Get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fWeightsList = (TList*)GetInputData(1); 
  }
  // Pass the list with weights to class:
  if(fWeightsList) fNL->SetWeightsList(fWeightsList);
 }
 
 fNL->Init();
 
 if(fNL->GetHistList()) 
 {
  fListHistos = fNL->GetHistList();
  // fListHistos->Print();
 } else 
   {
    Printf("ERROR: Could not retrieve histogram list (NL, Task::UserCreateOutputObjects()) !!!!"); 
   }
  
 PostData(1,fListHistos);
 
} // end of void AliAnalysisTaskNestedLoops::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskNestedLoops::UserExec(Option_t *) 
{
 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Nested Loops:
 if(fEvent) 
 {
  fNL->Make(fEvent);
 } else 
   {
    cout<<"WARNING: No input data (NL, Task::UserExec()) !!!!"<<endl;
    cout<<endl;
   }
  
 PostData(1,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskNestedLoops::Terminate(Option_t *) 
{
 //accessing the merged output list: 
 fListHistos = (TList*)GetOutputData(1);
 
 fNL = new AliFlowAnalysisWithNestedLoops(); 
 
 if(fListHistos) 
 {
  fNL->GetOutputHistograms(fListHistos);
  fNL->Finish();
  PostData(1,fListHistos);
 } else
   {
    cout<<" WARNING: histogram list pointer is empty (NL, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
    
} // end of void AliAnalysisTaskNestedLoops::Terminate(Option_t *)





















