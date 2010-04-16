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
 
class TFile;
class TList;
class AliAnalysisTaskSE; 
 
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskQCumulants.h"
#include "AliFlowAnalysisWithQCumulants.h"

ClassImp(AliAnalysisTaskQCumulants)

//================================================================================================================

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fEvent(NULL),
 fQC(NULL), 
 fListHistos(NULL),
 fHarmonic(2),  
 fApplyCorrectionForNUA(kTRUE), 
 fCalculate2DFlow(kFALSE),
 fStoreDistributions(kFALSE),
 fUseParticleWeights(useParticleWeights),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fWeightsList(NULL),
 fMultiplicityWeight(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t useParticleWeights)"<<endl;
 
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
 
 // Event weights:
 fMultiplicityWeight = new TString("combinations");
}

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(): 
 AliAnalysisTaskSE(),
 fEvent(NULL),
 fQC(NULL),
 fListHistos(NULL),
 fHarmonic(0),  
 fApplyCorrectionForNUA(kFALSE), 
 fCalculate2DFlow(kFALSE),
 fStoreDistributions(kFALSE),
 fUseParticleWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fWeightsList(NULL),
 fMultiplicityWeight(NULL)
{
 // Dummy constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize
 cout<<"AliAnalysisTaskQCumulants::UserCreateOutputObjects()"<<endl;

 // Analyser:
 fQC = new AliFlowAnalysisWithQCumulants();
 
 // Common:
 fQC->SetHarmonic(fHarmonic);
 fQC->SetApplyCorrectionForNUA(fApplyCorrectionForNUA);
 fQC->SetCalculate2DFlow(fCalculate2DFlow);
 fQC->SetStoreDistributions(fStoreDistributions);
 // Particle weights:
 if(fUseParticleWeights)
 {
  // Pass the flags to class:
  if(fUsePhiWeights) fQC->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fQC->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fQC->SetUseEtaWeights(fUseEtaWeights);
  // Get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fWeightsList = (TList*)GetInputData(1); 
  }
  // Pass the list with weights to class:
  if(fWeightsList) fQC->SetWeightsList(fWeightsList);
 }
 // Event weights:
 if(!(strcmp(fMultiplicityWeight->Data(),"combinations")==0)) // default is "combinations"
 {
  fQC->SetMultiplicityWeight(fMultiplicityWeight->Data());
 }

 fQC->Init();
 
 if(fQC->GetHistList()) 
 {
  fListHistos = fQC->GetHistList();
  // fListHistos->Print();
 } else 
   {
    Printf("ERROR: Could not retrieve histogram list (QC, Task::UserCreateOutputObjects()) !!!!"); 
   }
 
 PostData(1,fListHistos);
  
} // end of void AliAnalysisTaskQCumulants::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskQCumulants::UserExec(Option_t *) 
{
 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Q-cumulants
 if(fEvent) 
 {
  fQC->Make(fEvent);
 } else 
   {
    cout<<"WARNING: No input data (QC, Task::UserExec()) !!!!"<<endl;
    cout<<endl;
   }
  
 PostData(1,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Terminate(Option_t *) 
{
 //accessing the merged output list: 
 fListHistos = (TList*)GetOutputData(1);
 
 fQC = new AliFlowAnalysisWithQCumulants(); 
 
 if(fListHistos) 
 {
  fQC->GetOutputHistograms(fListHistos);
  fQC->Finish();
  PostData(1,fListHistos);
 } else
   {
    cout<<" WARNING: histogram list pointer is empty (QC, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
    
} // end of void AliAnalysisTaskQCumulants::Terminate(Option_t *)





















