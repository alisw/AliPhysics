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
 * analysis task for mixed harmomics  * 
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
#include "AliAnalysisTaskMixedHarmonics.h"
#include "AliFlowAnalysisWithMixedHarmonics.h"

ClassImp(AliAnalysisTaskMixedHarmonics)

//================================================================================================================

AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(const char *name, Bool_t useParticleWeights): 
AliAnalysisTaskSE(name), 
fEvent(NULL),
fMH(NULL), 
fListHistos(NULL),
fCorrelatorInteger(1),
fNoOfMultipicityBins(10),
fMultipicityBinWidth(2.),
fMinMultiplicity(3.),
fCorrectForDetectorEffects(kTRUE),
fUseParticleWeights(useParticleWeights),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(const char *name, Bool_t useParticleWeights)"<<endl;
 
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

AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(): 
AliAnalysisTaskSE(),
fEvent(NULL),
fMH(NULL),
fListHistos(NULL),
fCorrelatorInteger(0),
fNoOfMultipicityBins(0),
fMultipicityBinWidth(0),
fMinMultiplicity(0),
fCorrectForDetectorEffects(kFALSE),
fUseParticleWeights(kFALSE),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fWeightsList(NULL)
{
 // Dummy constructor
 cout<<"AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize
 cout<<"AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects()"<<endl;

 // Analyser:
 fMH = new AliFlowAnalysisWithMixedHarmonics();
  
 // Common:
 fMH->SetCorrelatorInteger(fCorrelatorInteger);
 fMH->SetNoOfMultipicityBins(fNoOfMultipicityBins);
 fMH->SetMultipicityBinWidth(fMultipicityBinWidth);
 fMH->SetMinMultiplicity(fMinMultiplicity);
 fMH->SetCorrectForDetectorEffects(fCorrectForDetectorEffects);
 if(fUseParticleWeights)
 {
  // Pass the flags to class:
  if(fUsePhiWeights) fMH->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fMH->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fMH->SetUseEtaWeights(fUseEtaWeights);
  // Get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fWeightsList = (TList*)GetInputData(1); 
  }
  // Pass the list with weights to class:
  if(fWeightsList) fMH->SetWeightsList(fWeightsList);
 }
 
 fMH->Init();
 
 if(fMH->GetHistList()) 
 {
  fListHistos = fMH->GetHistList();
  // fListHistos->Print();
 } else 
   {
    Printf("ERROR: Could not retrieve histogram list (MH, Task::UserCreateOutputObjects()) !!!!"); 
   }
 
 PostData(1,fListHistos);
  
} // end of void AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::UserExec(Option_t *) 
{
 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // Mixed Harmonics:
 if(fEvent) 
 {
  fMH->Make(fEvent);
 } else 
   {
    cout<<"WARNING: No input data (MH, Task::UserExec()) !!!!"<<endl;
    cout<<endl;
   }
  
 PostData(1,fListHistos);
}

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::Terminate(Option_t *) 
{
 //accessing the merged output list: 
 fListHistos = (TList*)GetOutputData(1);
 
 fMH = new AliFlowAnalysisWithMixedHarmonics(); 
 
 if(fListHistos) 
 {
  fMH->GetOutputHistograms(fListHistos);
  fMH->Finish();
  PostData(1,fListHistos);
 } else
   {
    cout<<" WARNING: histogram list pointer is empty (MH, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
    
} // end of void AliAnalysisTaskMixedHarmonics::Terminate(Option_t *)





















