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

/****************************************
 * analysis task for flow analysis with *
 *     multi-particle correlations      * 
 *                                      * 
 * author: Ante Bilandzic               *
 *         (abilandzic@gmail.com)       * 
 ***************************************/
  
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskMultiparticleCorrelations.h"
#include "AliFlowAnalysisWithMultiparticleCorrelations.h"
#include "AliLog.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMultiparticleCorrelations)

//================================================================================================================

AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fEvent(NULL),
 fMPC(NULL), 
 fHistList(NULL),
 fUseInternalFlags(kFALSE),
 fMinNoRPs(-44),
 fMaxNoRPs(-44),
 fExactNoRPs(-44),
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),
 fFillMultCorrelationsHist(kFALSE),
 fCalculateQvector(kFALSE),
 fCalculateDiffQvectors(kFALSE),
 fCalculateCorrelations(kFALSE),
 fCalculateIsotropic(kFALSE),
 fCalculateSame(kFALSE),
 fSkipZeroHarmonics(kFALSE),
 fCalculateSameIsotropic(kFALSE),
 fCalculateAll(kFALSE),
 fDontGoBeyond(0),
 fCalculateOnlyForHarmonicQC(kFALSE),
 fCalculateOnlyForSC(kFALSE),
 fCalculateOnlyCos(kFALSE),
 fCalculateOnlySin(kFALSE),
 fCalculateEbECumulants(kFALSE),
 fCrossCheckWithNestedLoops(kFALSE),
 fCrossCheckDiffWithNestedLoops(kFALSE),
 fCalculateStandardCandles(kFALSE),
 fPropagateErrorSC(kTRUE),
 fCalculateQcumulants(kFALSE),
 fHarmonicQC(2),
 fPropagateErrorQC(kTRUE),
 fCalculateDiffCorrelations(kFALSE)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights)");

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

  // Initialize all arrays:
  for(Int_t rp=0;rp<2;rp++) // [RP,POI]
  {
   for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
   {
    fUseWeights[rp][ppe] = kFALSE;
    fWeightsHist[rp][ppe] = NULL; 
   }
  }
 
} // AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations(): 
 AliAnalysisTaskSE(),
 fEvent(NULL),
 fMPC(NULL),
 fHistList(NULL),
 fUseInternalFlags(kFALSE),
 fMinNoRPs(-44),
 fMaxNoRPs(-44),
 fExactNoRPs(-44),
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),
 fFillMultCorrelationsHist(kFALSE),
 fCalculateQvector(kFALSE),
 fCalculateDiffQvectors(kFALSE),
 fCalculateCorrelations(kFALSE),
 fCalculateIsotropic(kFALSE),
 fCalculateSame(kFALSE),
 fSkipZeroHarmonics(kFALSE),
 fCalculateSameIsotropic(kFALSE),
 fCalculateAll(kFALSE),
 fDontGoBeyond(0),
 fCalculateOnlyForHarmonicQC(kFALSE),
 fCalculateOnlyForSC(kFALSE),
 fCalculateOnlyCos(kFALSE),
 fCalculateOnlySin(kFALSE),
 fCalculateEbECumulants(kFALSE),
 fCrossCheckWithNestedLoops(kFALSE),
 fCrossCheckDiffWithNestedLoops(kFALSE),
 fCalculateStandardCandles(kFALSE),
 fPropagateErrorSC(kFALSE),
 fCalculateQcumulants(kFALSE),
 fHarmonicQC(0),
 fPropagateErrorQC(kFALSE),
 fCalculateDiffCorrelations(kFALSE)
 {
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations()");

  // Initialize all arrays:
  for(Int_t rp=0;rp<2;rp++) // [RP,POI]
  {
   for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
   {
    fUseWeights[rp][ppe] = kFALSE;
    fWeightsHist[rp][ppe] = NULL; 
   }
  }

} // AliAnalysisTaskMultiparticleCorrelations::AliAnalysisTaskMultiparticleCorrelations():

//================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.
  
 AliDebug(2,"AliAnalysisTaskMultiparticleCorrelations::UserCreateOutputObjects()");
 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::UserCreateOutputObjects()";

 // Analyser:
 fMPC = new AliFlowAnalysisWithMultiparticleCorrelations();
 
 // Setters:
 if(fUseInternalFlags){fMPC->SetMinNoRPs(fMinNoRPs);}
 if(fUseInternalFlags){fMPC->SetMaxNoRPs(fMaxNoRPs);}
 if(fUseInternalFlags){fMPC->SetExactNoRPs(fExactNoRPs);}
 fMPC->SetFillControlHistograms(fFillControlHistograms);
 fMPC->SetFillKinematicsHist(fFillKinematicsHist);
 fMPC->SetFillMultDistributionsHist(fFillMultDistributionsHist);
 fMPC->SetFillMultCorrelationsHist(fFillMultCorrelationsHist);
 fMPC->SetCalculateQvector(fCalculateQvector);
 fMPC->SetCalculateDiffQvectors(fCalculateDiffQvectors);
 fMPC->SetCalculateCorrelations(fCalculateCorrelations);
 fMPC->SetCalculateIsotropic(fCalculateIsotropic);
 fMPC->SetCalculateSame(fCalculateSame);
 fMPC->SetSkipZeroHarmonics(fSkipZeroHarmonics);
 fMPC->SetCalculateSameIsotropic(fCalculateSameIsotropic);
 fMPC->SetCalculateAll(fCalculateAll);
 fMPC->SetDontGoBeyond(fDontGoBeyond);
 fMPC->SetCalculateOnlyForHarmonicQC(fCalculateOnlyForHarmonicQC);
 fMPC->SetCalculateOnlyForSC(fCalculateOnlyForSC);
 fMPC->SetCalculateOnlyCos(fCalculateOnlyCos);
 fMPC->SetCalculateOnlySin(fCalculateOnlySin);
 fMPC->SetCalculateEbECumulants(fCalculateEbECumulants);
 fMPC->SetCrossCheckWithNestedLoops(fCrossCheckWithNestedLoops);
 fMPC->SetCrossCheckDiffWithNestedLoops(fCrossCheckDiffWithNestedLoops);
 fMPC->SetCalculateStandardCandles(fCalculateStandardCandles);
 fMPC->SetPropagateErrorSC(fPropagateErrorSC);
 fMPC->SetCalculateQcumulants(fCalculateQcumulants);
 fMPC->SetHarmonicQC(fHarmonicQC);
 fMPC->SetPropagateErrorQC(fPropagateErrorQC);
 fMPC->SetCalculateDiffCorrelations(fCalculateDiffCorrelations);
 // Weights:
 TString type[2] = {"RP","POI"};
 TString variable[3] = {"phi","pt","eta"};
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   if(fUseWeights[rp][ppe])
   {
    if(!fWeightsHist[rp][ppe]){Fatal(sMethodName.Data(),"fWeightsHist[%d][%d]",rp,ppe);}
    else{fMPC->SetWeightsHist(fWeightsHist[rp][ppe],type[rp].Data(),variable[ppe].Data());}
   }
  } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
 } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]

 // Initialize:
 fMPC->Init();
 if(fMPC->GetHistList()) 
 {
  fHistList = fMPC->GetHistList();
  // fHistList->Print();
 } else 
   {
    Printf("ERROR: Could not retrieve histogram list (MPC, Task::UserCreateOutputObjects()) !!!!"); 
   }
 
 PostData(1,fHistList);
  
} // void AliAnalysisTaskMultiparticleCorrelations::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 // It's time for multi-particle correlations:
 if(fEvent) 
 {
  fMPC->Make(fEvent);
 } else 
   {
    cout<<" WARNING: No input data (MPC, Task::UserExec()) !!!!"<<endl;
    cout<<endl;
   }
  
 PostData(1,fHistList);

} // void AliAnalysisTaskMultiparticleCorrelations::UserExec(Option_t *) 

//================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::Terminate(Option_t *) 
{
 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);
 
 fMPC = new AliFlowAnalysisWithMultiparticleCorrelations(); 
 
 if(fHistList) 
 {
  fMPC->GetOutputHistograms(fHistList);
  fMPC->Finish();
  PostData(1,fHistList);
 } else
   {
    cout<<" WARNING: fHistList is NULL (MPC, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
    
} // end of void AliAnalysisTaskMultiparticleCorrelations::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetWeightsHist(TH1D* const hist, const char *type, const char *variable)
{
 // Pass histogram holding weights from an external file to the corresponding data member. 
 
 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetWeightsHist(TH1D* const hist, const char *type, const char *variable)";
 
 // Basic protection:
 if(!hist){Fatal(sMethodName.Data(),"hist");}
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI"))){Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);}
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);}

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 // Finally:
 hist->SetDirectory(0);
 fWeightsHist[rp][ppe] = (TH1D*)hist->Clone();
 if(!fWeightsHist[rp][ppe]){Fatal(sMethodName.Data(),"fWeightsHist[%d][%d]",rp,ppe);}

 fUseWeights[rp][ppe] = kTRUE; 

} // void AliAnalysisTaskMultiparticleCorrelations::SetWeightsHist(TH1D* const hwh, const char *type, const char *variable)

//================================================================================================================



