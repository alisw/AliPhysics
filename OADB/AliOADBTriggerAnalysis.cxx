//-------------------------------------------------------------------------
//     OADB container for trigger analysis configuration (cut ranges.. ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include "AliOADBTriggerAnalysis.h"
#include "AliLog.h"
#include "TBrowser.h"
#include "TFolder.h"
#include <iostream>
#include "TObjString.h"
#include "TObjArray.h"

using namespace std;

ClassImp(AliOADBTriggerAnalysis)

AliOADBTriggerAnalysis::AliOADBTriggerAnalysis(TString name) : TNamed(name.Data(), "OADB object storing trigger analysis settings"),
fFMDLowCut(0.2),
fFMDHitCut(0.5),
fZDCCutRefSum(-568.5),
fZDCCutRefDelta(-2.1),
fZDCCutSigmaSum(3.25),
fZDCCutSigmaDelta(2.25),
fZDCCutRefSumCorr(-65.5),
fZDCCutRefDeltaCorr(-2.1),
fZDCCutSigmaSumCorr(6.0),
fZDCCutSigmaDeltaCorr(1.2),
fZDCCutZNATimeCorrMin(0.0),
fZDCCutZNATimeCorrMax(2.0),
fZDCCutZNCTimeCorrMin(0.0),
fZDCCutZNCTimeCorrMax(5.0),
fSPDClsVsTklA(65),
fSPDClsVsTklB(4),
fV0MOnVsOfA(-145.),
fV0MOnVsOfB(7.2),
fSPDOnVsOfA(-4.16),
fSPDOnVsOfB(0.84),
fVtxMinContributors(3),
fVtxMinZdist(0.8),
fVtxNSigmaZdist(3.),
fVtxNSigmaDiamXY(2.),
fVtxNSigmaDiamZ(5.),
fV0CasymA(-25),
fV0CasymB(0.15),
fNBCsPast(4),
fNBCsFuture(7),
fVHMBBAflags(-1),
fVHMBBCflags(-1),
fVHMBGAflags(-1),
fVHMBGCflags(-1),
fVIRBBAflags(10),
fVIRBBCflags(10),
fVIRBGAflags(33),
fVIRBGCflags(33),
fV0MOnThreshold(-1),
fV0MOfThreshold(400),
fSPDGFOThreshold(2),
fSH1OuterThreshold(-1),
fSH2OuterThreshold(-1),
fTRDptHSE(3.),
fTRDpidHSE(144),
fTRDptHQU(2.),
fTRDpidHQU(164.),
fTRDptHEE(3.),
fTRDpidHEE(144),
fTRDminSectorHEE(6),
fTRDmaxSectorHEE(8),
fTRDptHJT(3.),
fTRDnHJT(3)
{
}

AliOADBTriggerAnalysis::~AliOADBTriggerAnalysis(){
}
  
void AliOADBTriggerAnalysis::Browse(TBrowser *b){
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly
  TObjArray* parArray = new TObjArray(50);
  parArray->Add(new TObjString(Form("FMDLowCut            %f", fFMDLowCut           )));
  parArray->Add(new TObjString(Form("FMDHitCut            %f", fFMDHitCut           )));
  parArray->Add(new TObjString(Form("ZDCCutRefSum         %f", fZDCCutRefSum        )));
  parArray->Add(new TObjString(Form("ZDCCutRefDelta       %f", fZDCCutRefDelta      )));
  parArray->Add(new TObjString(Form("ZDCCutSigmaSum       %f", fZDCCutSigmaSum      )));
  parArray->Add(new TObjString(Form("ZDCCutSigmaDelta     %f", fZDCCutSigmaDelta    )));
  parArray->Add(new TObjString(Form("ZDCCutRefSumCorr     %f", fZDCCutRefSumCorr    )));
  parArray->Add(new TObjString(Form("ZDCCutRefDeltaCorr   %f", fZDCCutRefDeltaCorr  )));
  parArray->Add(new TObjString(Form("ZDCCutSigmaSumCorr   %f", fZDCCutSigmaSumCorr  )));
  parArray->Add(new TObjString(Form("ZDCCutSigmaDeltaCorr %f", fZDCCutSigmaDeltaCorr)));
  parArray->Add(new TObjString(Form("ZDCCutZNATimeCorrMin %f", fZDCCutZNATimeCorrMin)));
  parArray->Add(new TObjString(Form("ZDCCutZNATimeCorrMax %f", fZDCCutZNATimeCorrMax)));
  parArray->Add(new TObjString(Form("ZDCCutZNCTimeCorrMin %f", fZDCCutZNCTimeCorrMin)));
  parArray->Add(new TObjString(Form("ZDCCutZNCTimeCorrMax %f", fZDCCutZNCTimeCorrMax)));
  parArray->Add(new TObjString(Form("SPDClsVsTklA         %f", fSPDClsVsTklA        )));
  parArray->Add(new TObjString(Form("SPDClsVsTklB         %f", fSPDClsVsTklB        )));
  parArray->Add(new TObjString(Form("V0MOnVsOfA           %f", fV0MOnVsOfA          )));
  parArray->Add(new TObjString(Form("V0MOnVsOfB           %f", fV0MOnVsOfB          )));
  parArray->Add(new TObjString(Form("SPDOnVsOfA           %f", fSPDOnVsOfA          )));
  parArray->Add(new TObjString(Form("SPDOnVsOfB           %f", fSPDOnVsOfB          )));
  parArray->Add(new TObjString(Form("VtxMinContributors   %i", fVtxMinContributors  )));
  parArray->Add(new TObjString(Form("VtxMinZdist          %f", fVtxMinZdist         )));
  parArray->Add(new TObjString(Form("VtxNSigmaZdist       %f", fVtxNSigmaZdist      )));
  parArray->Add(new TObjString(Form("VtxNSigmaDiamXY      %f", fVtxNSigmaDiamXY     )));
  parArray->Add(new TObjString(Form("VtxNSigmaDiamZ       %f", fVtxNSigmaDiamZ      )));
  parArray->Add(new TObjString(Form("V0CasymA             %f", fV0CasymA            )));
  parArray->Add(new TObjString(Form("V0CasymB             %f", fV0CasymB            )));
  parArray->Add(new TObjString(Form("NBCsPast             %i", fNBCsPast            )));
  parArray->Add(new TObjString(Form("NBCsFuture           %i", fNBCsFuture          )));
  parArray->Add(new TObjString(Form("VHMBBAflags          %i", fVHMBBAflags         )));
  parArray->Add(new TObjString(Form("VHMBBCflags          %i", fVHMBBCflags         )));
  parArray->Add(new TObjString(Form("VHMBGAflags          %i", fVHMBGAflags         )));
  parArray->Add(new TObjString(Form("VHMBGCflags          %i", fVHMBGCflags         )));
  parArray->Add(new TObjString(Form("VIRBBAflags          %i", fVIRBBAflags         )));
  parArray->Add(new TObjString(Form("VIRBBCflags          %i", fVIRBBCflags         )));
  parArray->Add(new TObjString(Form("VIRBGAflags          %i", fVIRBGAflags         )));
  parArray->Add(new TObjString(Form("VIRBGCflags          %i", fVIRBGCflags         )));
  parArray->Add(new TObjString(Form("V0MOnThreshold       %i", fV0MOnThreshold      )));
  parArray->Add(new TObjString(Form("V0MOfThreshold       %f", fV0MOfThreshold      )));
  parArray->Add(new TObjString(Form("SPDGFOThreshold      %i", fSPDGFOThreshold     )));
  parArray->Add(new TObjString(Form("SH1OuterThreshold    %i", fSH1OuterThreshold   )));
  parArray->Add(new TObjString(Form("SH2OuterThreshold    %i", fSH2OuterThreshold   )));
  parArray->Add(new TObjString(Form("TRDptHSE             %f", fTRDptHSE            )));
  parArray->Add(new TObjString(Form("TRDpidHSE            %i", fTRDpidHSE           )));
  parArray->Add(new TObjString(Form("TRDptHQU             %f", fTRDptHQU            )));
  parArray->Add(new TObjString(Form("TRDpidHQU            %i", fTRDpidHQU           )));
  parArray->Add(new TObjString(Form("TRDptHEE             %f", fTRDptHEE            )));
  parArray->Add(new TObjString(Form("TRDpidHEE            %i", fTRDpidHEE           )));
  parArray->Add(new TObjString(Form("TRDminSectorHEE      %f", fTRDminSectorHEE     )));
  parArray->Add(new TObjString(Form("TRDmaxSectorHEE      %f", fTRDmaxSectorHEE     )));
  parArray->Add(new TObjString(Form("TRDptHJT             %f", fTRDptHJT            )));
  parArray->Add(new TObjString(Form("TRDnHJT              %i", fTRDnHJT             )));
  
  if (b) b->Add(parArray);
  else TObject::Browse(b);
}

void AliOADBTriggerAnalysis::Print(Option_t* option) const {
  // Print Class contents
  printf("Physics selection cut settings:\n");
  printf(Form("  FMDLowCut            %f\n", fFMDLowCut           ));
  printf(Form("  FMDHitCut            %f\n", fFMDHitCut           ));
  printf(Form("  ZDCCutRefSum         %f\n", fZDCCutRefSum        ));
  printf(Form("  ZDCCutRefDelta       %f\n", fZDCCutRefDelta      ));
  printf(Form("  ZDCCutSigmaSum       %f\n", fZDCCutSigmaSum      ));
  printf(Form("  ZDCCutSigmaDelta     %f\n", fZDCCutSigmaDelta    ));
  printf(Form("  ZDCCutRefSumCorr     %f\n", fZDCCutRefSumCorr    ));
  printf(Form("  ZDCCutRefDeltaCorr   %f\n", fZDCCutRefDeltaCorr  ));
  printf(Form("  ZDCCutSigmaSumCorr   %f\n", fZDCCutSigmaSumCorr  ));
  printf(Form("  ZDCCutSigmaDeltaCorr %f\n", fZDCCutSigmaDeltaCorr));
  printf(Form("  ZDCCutZNATimeCorrMin %f\n", fZDCCutZNATimeCorrMin));
  printf(Form("  ZDCCutZNATimeCorrMax %f\n", fZDCCutZNATimeCorrMax));
  printf(Form("  ZDCCutZNCTimeCorrMin %f\n", fZDCCutZNCTimeCorrMin));
  printf(Form("  ZDCCutZNCTimeCorrMax %f\n", fZDCCutZNCTimeCorrMax));
  printf(Form("  SPDClsVsTklA         %f\n", fSPDClsVsTklA        ));
  printf(Form("  SPDClsVsTklB         %f\n", fSPDClsVsTklB        ));
  printf(Form("  V0MOnVsOfA           %f\n", fV0MOnVsOfA          ));
  printf(Form("  V0MOnVsOfB           %f\n", fV0MOnVsOfB          ));
  printf(Form("  SPDOnVsOfA           %f\n", fSPDOnVsOfA          ));
  printf(Form("  SPDOnVsOfB           %f\n", fSPDOnVsOfB          ));
  printf(Form("  VtxMinContributors   %i\n", fVtxMinContributors  ));
  printf(Form("  VtxMinZdist          %f\n", fVtxMinZdist         ));
  printf(Form("  VtxNSigmaZdist       %f\n", fVtxNSigmaZdist      ));
  printf(Form("  VtxNSigmaDiamXY      %f\n", fVtxNSigmaDiamXY     ));
  printf(Form("  VtxNSigmaDiamZ       %f\n", fVtxNSigmaDiamZ      ));
  printf(Form("  V0CasymA             %f\n", fV0CasymA            ));
  printf(Form("  V0CasymB             %f\n", fV0CasymB            ));
  printf(Form("  NBCsPast             %i\n", fNBCsPast            ));
  printf(Form("  NBCsFuture           %i\n", fNBCsFuture          ));
  printf(Form("  VHMBBAflags          %i\n", fVHMBBAflags         ));
  printf(Form("  VHMBBCflags          %i\n", fVHMBBCflags         ));
  printf(Form("  VHMBGAflags          %i\n", fVHMBGAflags         ));
  printf(Form("  VHMBGCflags          %i\n", fVHMBGCflags         ));
  printf(Form("  VIRBBAflags          %i\n", fVIRBBAflags         ));
  printf(Form("  VIRBBCflags          %i\n", fVIRBBCflags         ));
  printf(Form("  VIRBGAflags          %i\n", fVIRBGAflags         ));
  printf(Form("  VIRBGCflags          %i\n", fVIRBGCflags         ));
  printf(Form("  V0MOnThreshold       %i\n", fV0MOnThreshold      ));
  printf(Form("  V0MOfThreshold       %f\n", fV0MOfThreshold      ));
  printf(Form("  SPDGFOThreshold      %i\n", fSPDGFOThreshold     ));
  printf(Form("  SH1OuterThreshold    %i\n", fSH1OuterThreshold   ));
  printf(Form("  SH2OuterThreshold    %i\n", fSH2OuterThreshold   ));
  printf(Form("  TRDptHSE             %f\n", fTRDptHSE            ));
  printf(Form("  TRDpidHSE            %i\n", fTRDpidHSE           ));
  printf(Form("  TRDptHQU             %f\n", fTRDptHQU            ));
  printf(Form("  TRDpidHQU            %i\n", fTRDpidHQU           ));
  printf(Form("  TRDptHEE             %f\n", fTRDptHEE            ));
  printf(Form("  TRDpidHEE            %i\n", fTRDpidHEE           ));
  printf(Form("  TRDminSectorHEE      %f\n", fTRDminSectorHEE     ));
  printf(Form("  TRDmaxSectorHEE      %f\n", fTRDmaxSectorHEE     ));
  printf(Form("  TRDptHJT             %f\n", fTRDptHJT            ));
  printf(Form("  TRDnHJT              %i\n", fTRDnHJT             ));
}

