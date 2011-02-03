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
 
class TFile;
class TList;
class AliAnalysisTaskSE; 
 
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCumulants.h"
#include "AliFlowAnalysisWithCumulants.h"

ClassImp(AliAnalysisTaskCumulants)

//================================================================================================================

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name, Bool_t useWeights): 
AliAnalysisTaskSE(name), 
fEvent(NULL),
fGFC(NULL),
fListHistos(NULL),
fHarmonic(2),  
fMultiple(1),
fCalculateVsMultiplicity(kFALSE),
fnBinsMult(10000),  
fMinMult(0.),   
fMaxMult(10000.),
fUseWeights(useWeights),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL),
fTuneParameters(kFALSE)
{
 // Constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name, Bool_t useWeights)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with an AliFlowEventSimple
 DefineInput(0, AliFlowEventSimple::Class());
 
 // Input slot #1 is needed for the weights input files 
 if(useWeights) 
 {
  DefineInput(1, TList::Class());   
 }
 // Output slot #0 is reserved
 // Output slot #1 writes into a TList container
 DefineOutput(1, TList::Class());   
 
 // Initilize arrays:
 for(Int_t r=0;r<10;r++)
 {
  fTuningR0[r] = 0.;
 }
} // end of constructor

AliAnalysisTaskCumulants::AliAnalysisTaskCumulants():
AliAnalysisTaskSE(),
fEvent(NULL),
fGFC(NULL),
fListHistos(NULL),
fHarmonic(0),  
fMultiple(0),
fCalculateVsMultiplicity(kFALSE),
fnBinsMult(0),  
fMinMult(0.),   
fMaxMult(0.),
fUseWeights(kFALSE),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL),
fTuneParameters(kFALSE)
{
 // Dummy constructor
 cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants()"<<endl;

 // Initilize arrays:
 for(Int_t r=0;r<10;r++)
 {
  fTuningR0[r] = 0.;
 }
}

//================================================================================================================

void AliAnalysisTaskCumulants::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize
 cout<<"AliAnalysisTaskCumulants::UserCreateOutputObjects()"<<endl;

 // Analyser:
 fGFC = new AliFlowAnalysisWithCumulants(); 
 fGFC->SetHarmonic(fHarmonic);
 
 // Calculation vs multiplicity:
 if(fCalculateVsMultiplicity)
 {
  fGFC->SetCalculateVsMultiplicity(fCalculateVsMultiplicity);
  fGFC->SetnBinsMult(fnBinsMult);
  fGFC->SetMinMult(fMinMult);
  fGFC->SetMaxMult(fMaxMult);
 }
 
 // Weights:
 if(fUseWeights)
 {
  // Pass the flags to class:
  if(fUsePhiWeights) fGFC->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fGFC->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fGFC->SetUseEtaWeights(fUseEtaWeights);
  // Get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fWeightsList = (TList*)GetInputData(1); 
  }
  // Pass the list with weights to class:
  if(fWeightsList) fGFC->SetWeightsList(fWeightsList);
 }

 // Tuning:
 if(fTuneParameters)
 {
  fGFC->SetTuneParameters(fTuneParameters); 
  for(Int_t r=0;r<10;r++) {fGFC->SetTuningR0(fTuningR0[r],r);}
 }
 
 fGFC->Init();

 if(fGFC->GetHistList()) 
 {
  fListHistos = fGFC->GetHistList();
  //fListHistos->Print();
 } else
   {
    Printf("ERROR: Could not retrieve histogram list (GFC, Task::UserCreateOutputObjects()) !!!!"); 
   }

 PostData(1,fListHistos);
   
} // end of void AliAnalysisTaskCumulants::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskCumulants::UserExec(Option_t *) 
{
 // Main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Generating function cumulants (GFC):
 if(fEvent) 
 {
  fGFC->Make(fEvent);
 } else 
   {
    cout<<"WARNING: No input data (GFC, Task::UserExec()) !!!!"<<endl;
   }
  
 PostData(1,fListHistos);
 
} // end of void AliAnalysisTaskCumulants::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskCumulants::Terminate(Option_t *) 
{  
 // Accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(1);
 //fListHistos->Print();
 
 fGFC = new AliFlowAnalysisWithCumulants();  
 
 if(fListHistos)
 {
  fGFC->GetOutputHistograms(fListHistos);
  fGFC->Finish();
  PostData(1,fListHistos);
 } else
   {
    cout<<"WARNING: histogram list pointer is empty (GFC, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
   
} // end of void AliAnalysisTaskCumulants::Terminate(Option_t *) 





















