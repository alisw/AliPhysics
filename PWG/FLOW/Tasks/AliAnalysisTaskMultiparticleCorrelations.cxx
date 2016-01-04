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
 fAnalysisTag(""),
 fDumpThePoints(kFALSE),
 fMaxNoEventsPerFile(100),
 fSelectRandomlyRPs(kFALSE),
 fnSelectedRandomlyRPs(-44),
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),
 fFillMultCorrelationsHist(kFALSE),
 fSkipSomeIntervals(kFALSE),
 fCalculateQvector(kFALSE),
 fCalculateDiffQvectors(kFALSE),
 fProduction(""),
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
 fCalculateDiffCorrelations(kFALSE),
 fCalculateDiffCos(kTRUE),
 fCalculateDiffSin(kFALSE),
 fCalculateDiffCorrelationsVsPt(kTRUE),
 fUseDefaultBinning(kTRUE),
 fnDiffBins(-44),
 fRangesDiffBins(NULL),
 fCalculateSymmetryPlanes(kFALSE),
 fCalculateEtaGaps(kFALSE)
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
  // For nested loops arrays:
  fCrossCheckDiffCSCOBN[0] = 0; // cos/sin
  fCrossCheckDiffCSCOBN[1] = 2; // correlator order
  fCrossCheckDiffCSCOBN[2] = 4; // bin number 
  // Initialize default vaulues for fDontFill[3]:
  fDontFill[0] = kFALSE;
  fDontFill[1] = kFALSE;
  fDontFill[2] = kFALSE;
  // Initialize default binning values for fKinematicsHist[2][3]:
  // nBins:
  fnBins[0][0] = 360;  // [RP][phi]
  fnBins[0][1] = 1000; // [RP][pt]
  fnBins[0][2] = 1000; // [RP][eta]
  fnBins[1][0] = 360;  // [POI][phi]
  fnBins[1][1] = 1000; // [POI][pt]
  fnBins[1][2] = 1000; // [POI][eta]
  // Min:
  fMin[0][0] = 0.;  // [RP][phi]
  fMin[0][1] = 0.;  // [RP][pt]
  fMin[0][2] = -1.; // [RP][eta]
  fMin[1][0] = 0.;  // [POI][phi]
  fMin[1][1] = 0.;  // [POI][pt]
  fMin[1][2] = -1.; // [POI][eta]
  // Max:
  fMax[0][0] = TMath::TwoPi(); // [RP][phi]
  fMax[0][1] = 10.;            // [RP][pt]
  fMax[0][2] = 1.;             // [RP][eta]
  fMax[1][0] = TMath::TwoPi(); // [POI][phi]
  fMax[1][1] = 10.;            // [POI][pt]
  fMax[1][2] = 1.;             // [POI][eta]
  // Initialize default binning values for fMultCorrelationsHist[3]:
  // nBins:
  fnBinsMult[0] = 3000; // [RP]
  fnBinsMult[1] = 3000; // [POI]
  fnBinsMult[2] = 3000; // [REF]
  // Min:
  fMinMult[0] = 0.; // [RP]
  fMinMult[1] = 0.; // [POI]
  fMinMult[2] = 0.; // [REF]
  // Max:
  fMaxMult[0] = 3000.; // [RP]
  fMaxMult[1] = 3000.; // [POI]
  fMaxMult[2] = 3000.; // [REF]
  // Initialize default rp, phi and eta intervals to be skipped:
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   for(Int_t i=0;i<10;i++) // interval boundaries, 5 intervals, 10 boundaries TBI
   { 
    fSkip[ppe][i] = -44.;
   } 
  } // for(Int_t ppe=0;ppe<3;ppe++)

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
 fAnalysisTag(""),
 fDumpThePoints(kFALSE),
 fMaxNoEventsPerFile(0),
 fSelectRandomlyRPs(kFALSE),
 fnSelectedRandomlyRPs(-44),
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),
 fFillMultCorrelationsHist(kFALSE),
 fSkipSomeIntervals(kFALSE),
 fCalculateQvector(kFALSE),
 fCalculateDiffQvectors(kFALSE),
 fProduction(""),
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
 fCalculateDiffCorrelations(kFALSE),
 fCalculateDiffCos(kTRUE),
 fCalculateDiffSin(kFALSE),
 fCalculateDiffCorrelationsVsPt(kTRUE),
 fUseDefaultBinning(kTRUE),
 fnDiffBins(-44),
 fRangesDiffBins(NULL),
 fCalculateSymmetryPlanes(kFALSE),
 fCalculateEtaGaps(kFALSE)
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
  // For nested loops arrays:
  fCrossCheckDiffCSCOBN[0] = 0; // cos/sin
  fCrossCheckDiffCSCOBN[1] = 2; // correlator order
  fCrossCheckDiffCSCOBN[2] = 4; // bin number 
  // Initialize default vaulues for fDontFill[3]:
  fDontFill[0] = kFALSE;
  fDontFill[1] = kFALSE;
  fDontFill[2] = kFALSE;
  // Initialize default binning values for fKinematicsHist[2][3]:
  // nBins:
  fnBins[0][0] = 360;  // [RP][phi]
  fnBins[0][1] = 1000; // [RP][pt]
  fnBins[0][2] = 1000; // [RP][eta]
  fnBins[1][0] = 360;  // [POI][phi]
  fnBins[1][1] = 1000; // [POI][pt]
  fnBins[1][2] = 1000; // [POI][eta]
  // Min:
  fMin[0][0] = 0.;  // [RP][phi]
  fMin[0][1] = 0.;  // [RP][pt]
  fMin[0][2] = -1.; // [RP][eta]
  fMin[1][0] = 0.;  // [POI][phi]
  fMin[1][1] = 0.;  // [POI][pt]
  fMin[1][2] = -1.; // [POI][eta]
  // Max:
  fMax[0][0] = TMath::TwoPi(); // [RP][phi]
  fMax[0][1] = 10.;            // [RP][pt]
  fMax[0][2] = 1.;             // [RP][eta]
  fMax[1][0] = TMath::TwoPi(); // [POI][phi]
  fMax[1][1] = 10.;            // [POI][pt]
  fMax[1][2] = 1.;             // [POI][eta]
  // Initialize default binning values for fMultCorrelationsHist[3]:
  // nBins:
  fnBinsMult[0] = 3000; // [RP]
  fnBinsMult[1] = 3000; // [POI]
  fnBinsMult[2] = 3000; // [REF]
  // Min:
  fMinMult[0] = 0.; // [RP]
  fMinMult[1] = 0.; // [POI]
  fMinMult[2] = 0.; // [REF]
  // Max:
  fMaxMult[0] = 3000.; // [RP]
  fMaxMult[1] = 3000.; // [POI]
  fMaxMult[2] = 3000.; // [REF]
  // Initialize default rp, phi and eta intervals to be skipped:
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   for(Int_t i=0;i<10;i++) // interval boundaries, 10 boundaries at max
   { 
    fSkip[ppe][i] = -44.;
   } 
  } // for(Int_t ppe=0;ppe<3;ppe++)

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
 fMPC->SetAnalysisTag(fAnalysisTag.Data());
 fMPC->SetDumpThePoints(fDumpThePoints,fMaxNoEventsPerFile);
 if(fSelectRandomlyRPs){fMPC->SetSelectRandomlyRPs(fnSelectedRandomlyRPs);}
 fMPC->SetFillControlHistograms(fFillControlHistograms);
 if(fDontFill[0]){fMPC->SetDontFill("RP");}
 if(fDontFill[1]){fMPC->SetDontFill("POI");}
 if(fDontFill[2]){fMPC->SetDontFill("REF");}
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
 fMPC->SetCrossCheckDiffCSCOBN(fCrossCheckDiffCSCOBN[0],fCrossCheckDiffCSCOBN[1],fCrossCheckDiffCSCOBN[2]);  
 fMPC->SetCalculateStandardCandles(fCalculateStandardCandles);
 fMPC->SetPropagateErrorSC(fPropagateErrorSC);
 fMPC->SetCalculateQcumulants(fCalculateQcumulants);
 fMPC->SetHarmonicQC(fHarmonicQC);
 fMPC->SetPropagateErrorQC(fPropagateErrorQC);
 fMPC->SetCalculateDiffCorrelations(fCalculateDiffCorrelations);
 fMPC->SetCalculateDiffCos(fCalculateDiffCos);
 fMPC->SetCalculateDiffSin(fCalculateDiffSin);
 fMPC->SetCalculateDiffCorrelationsVsPt(fCalculateDiffCorrelationsVsPt);
 fMPC->SetUseDefaultBinning(fUseDefaultBinning);
 fMPC->SetnDiffBins(fnDiffBins);
 fMPC->SetRangesDiffBins(fRangesDiffBins);
 fMPC->SetCalculateSymmetryPlanes(fCalculateSymmetryPlanes);
 fMPC->SetCalculateEtaGaps(fCalculateEtaGaps);

 // Weights:
 TString type[2] = {"RP","POI"};
 TString variable[3] = {"phi","pt","eta"}; 
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   if(fUseWeights[rp][ppe])
   {
    if(!fWeightsHist[rp][ppe])
    {
     if(fProduction.EqualTo("")){Fatal(sMethodName.Data(),"fProduction is empty, for one reason or another...");}
     fWeightsHist[rp][ppe] = GetHistogramWithWeights(TString(Form("%s/%s",gSystem->pwd(),"weights.root")).Data(),TString(this->fName).Data(),type[rp].Data(),variable[ppe].Data(),fProduction.Data());
    }
    if(!fWeightsHist[rp][ppe])
    {
     Fatal(sMethodName.Data(),"fWeightsHist[%d][%d]",rp,ppe);
    } else{fMPC->SetWeightsHist(fWeightsHist[rp][ppe],type[rp].Data(),variable[ppe].Data());}
   } // if(fUseWeights[rp][ppe])
  } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
 } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]

 // Control histos:
 // Kinematics:
 //TString typeKine[2] = {"RP","POI"};
 //TString variable[3] = {"phi","pt","eta"};
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fMPC->SetnBins(type[rp].Data(),variable[ppe].Data(),fnBins[rp][ppe]);
   fMPC->SetMin(type[rp].Data(),variable[ppe].Data(),fMin[rp][ppe]);
   fMPC->SetMax(type[rp].Data(),variable[ppe].Data(),fMax[rp][ppe]);
  } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
 } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 // Multiplicites:
 TString typeMult[3] = {"RP","POI","REF"};
 for(Int_t rpr=0;rpr<3;rpr++) // [RP,POI,REF]
 {  
  fMPC->SetnBinsMult(typeMult[rpr].Data(),fnBinsMult[rpr]);
  fMPC->SetMinMult(typeMult[rpr].Data(),fMinMult[rpr]);
  fMPC->SetMaxMult(typeMult[rpr].Data(),fMaxMult[rpr]);
 } // for(Int_t rpr=0;rpr<3;rpr++) // [RP,POI,REF]

 // Intervals to skip:
 if(fSkipSomeIntervals)
 {
  // TBI some things are clearly hardwired here:
  fMPC->SetIntervalsToSkip("Phi",10,fSkip[0]);
  fMPC->SetIntervalsToSkip("Pt",10,fSkip[1]);
  fMPC->SetIntervalsToSkip("Eta",10,fSkip[2]);
 }

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

void AliAnalysisTaskMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)
{
 // Set number of bins for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fnBins[rp][ppe] = nBins;

} // void AliAnalysisTaskMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)
{
 // Set min bin range for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fMin[rp][ppe] = min;

} // void AliAnalysisTaskMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t max)
{
 // Set max bin range for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t max)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fMax[rp][ppe] = max;

} // void AliAnalysisTaskMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t min)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)
{
 // Set number of bins for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fnBinsMult[rpr] = nBinsMult;

} // void AliAnalysisTaskMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)
{
 // Set min bin range for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fMinMult[rpr] = minMult;

} // void AliAnalysisTaskMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetMaxMult(const char *type, Double_t maxMult)
{
 // Set max bin range for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetMaxMult(const char *type, Double_t maxMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fMaxMult[rpr] = maxMult;

} // void AliAnalysisTaskMultiparticleCorrelations::SetMaxMult(const char *type, Double_t minMult)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t nBoundaries, Double_t *boundaries)
{
 // Set all pt, phi and eta intervals to be skipped. 

 // Example usage in the steering macro (before Init()):
 //  Double_t skip[4] = {-0.1,0.2,0.8,0.9};
 //  taskMPC->SetIntervalsToSkip("Eta",4,skip);
 
 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t n, Double_t *boundaries)";
 
 // Basic protection:
 if(!(TString(ppe).EqualTo("Phi") || TString(ppe).EqualTo("Pt") || TString(ppe).EqualTo("Eta")))
 {
  cout<<"Well, could you perhaps try to use only Phi, Pt or Eta here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(ppe).EqualTo... type = %s ",ppe);
 }
  
 if(nBoundaries>10)
 {
  cout<<"Maximum number of boundaries is hardwired to be 10 at the moment, sorry..."<<endl;
  Fatal(sMethodName.Data(),"nBoundaries = %d ",nBoundaries);
 }

 fSkipSomeIntervals = kTRUE;

 Int_t index = -44;
 if(TString(ppe).EqualTo("Phi"))
 {
  index = 0;
 } 
 else if(TString(ppe).EqualTo("Pt"))
 {
  index = 1;
 }
 else
 {
  index = 2;
 }

 for(Int_t b=0;b<nBoundaries;b++) // boundaries
 {
  fSkip[index][b] = boundaries[b];
 }

} // void AliAnalysisTaskMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t n, Double_t *boundaries)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)
{
 // Set harmonics for all differential correlators.

 TString sMethodName = "void AliAnalysisTaskMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)";

 // TBI to be finalized, along the same lines as it was done in void AliFlowAnalysisWithMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)
 order = harmonics[0]; // TBI eliminating temporarily the warnings
 if(order>44) harmonics[0] = -44; // TBI eliminating temporarily the warnings, the implementation of this method has to be finalized one day...

} // void AliAnalysisTaskMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)

//=======================================================================================================================




