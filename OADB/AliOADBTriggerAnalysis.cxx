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
fZDCCutRefSumCorr(0.),
fZDCCutRefDeltaCorr(0.),
fZDCCutSigmaSumCorr(2.),
fZDCCutSigmaDeltaCorr(2.),
fZDCCutZNATimeCorrMin(5.0),
fZDCCutZNATimeCorrMax(100.0),
fZDCCutZNCTimeCorrMin(5.0),
fZDCCutZNCTimeCorrMax(100.0),
fSPDClsVsTklA(65),
fSPDClsVsTklB(4),
fV0C012vsTklA(150),
fV0C012vsTklB(20),
fV0MOnVsOfA(-59.56),
fV0MOnVsOfB(5.22),
fSPDOnVsOfA(-5.62),
fSPDOnVsOfB(0.85),
fVtxMinContributors(5),
fVtxMinZdist(0.8),
fVtxNSigmaZdist(3.),
fVtxNSigmaDiamXY(2.),
fVtxNSigmaDiamZ(5.),
fV0CasymA(-25),
fV0CasymB(0.15),
fNBCsPast(4),
fNBCsFuture(7),
fVIRBBAflags(10),
fVIRBBCflags(10),
fVIRBGAflags(33),
fVIRBGCflags(33),
fVHMBBAflags(-1),
fVHMBBCflags(-1),
fVHMBGAflags(33),
fVHMBGCflags(33),
fV0MOnThreshold(-1),
fV0MOfThreshold(-1),
fSPDGFOThreshold(2),
fSH1OuterThreshold(-1),
fSH2OuterThreshold(-1),
fTklThreshold(-1),
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
  parArray->SetName("parameters");
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
  parArray->Add(new TObjString(Form("V0C012vsTklA         %f", fV0C012vsTklA        )));
  parArray->Add(new TObjString(Form("V0C012vsTklB         %f", fV0C012vsTklB        )));
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
  parArray->Add(new TObjString(Form("TklThreshold         %i", fTklThreshold        )));
  parArray->Add(new TObjString(Form("TRDptHSE             %f", fTRDptHSE            )));
  parArray->Add(new TObjString(Form("TRDpidHSE            %i", fTRDpidHSE           )));
  parArray->Add(new TObjString(Form("TRDptHQU             %f", fTRDptHQU            )));
  parArray->Add(new TObjString(Form("TRDpidHQU            %i", fTRDpidHQU           )));
  parArray->Add(new TObjString(Form("TRDptHEE             %f", fTRDptHEE            )));
  parArray->Add(new TObjString(Form("TRDpidHEE            %i", fTRDpidHEE           )));
  parArray->Add(new TObjString(Form("TRDminSectorHEE      %i", fTRDminSectorHEE     )));
  parArray->Add(new TObjString(Form("TRDmaxSectorHEE      %i", fTRDmaxSectorHEE     )));
  parArray->Add(new TObjString(Form("TRDptHJT             %f", fTRDptHJT            )));
  parArray->Add(new TObjString(Form("TRDnHJT              %i", fTRDnHJT             )));
  
  if (b) b->Add(parArray);
  else TObject::Browse(b);
}

void AliOADBTriggerAnalysis::Print(Option_t* option) const {
  // Print Class contents
  printf("Physics selection cut settings:\n");
  printf("  FMDLowCut            %f\n", fFMDLowCut           );
  printf("  FMDHitCut            %f\n", fFMDHitCut           );
  printf("  ZDCCutRefSum         %f\n", fZDCCutRefSum        );
  printf("  ZDCCutRefDelta       %f\n", fZDCCutRefDelta      );
  printf("  ZDCCutSigmaSum       %f\n", fZDCCutSigmaSum      );
  printf("  ZDCCutSigmaDelta     %f\n", fZDCCutSigmaDelta    );
  printf("  ZDCCutRefSumCorr     %f\n", fZDCCutRefSumCorr    );
  printf("  ZDCCutRefDeltaCorr   %f\n", fZDCCutRefDeltaCorr  );
  printf("  ZDCCutSigmaSumCorr   %f\n", fZDCCutSigmaSumCorr  );
  printf("  ZDCCutSigmaDeltaCorr %f\n", fZDCCutSigmaDeltaCorr);
  printf("  ZDCCutZNATimeCorrMin %f\n", fZDCCutZNATimeCorrMin);
  printf("  ZDCCutZNATimeCorrMax %f\n", fZDCCutZNATimeCorrMax);
  printf("  ZDCCutZNCTimeCorrMin %f\n", fZDCCutZNCTimeCorrMin);
  printf("  ZDCCutZNCTimeCorrMax %f\n", fZDCCutZNCTimeCorrMax);
  printf("  SPDClsVsTklA         %f\n", fSPDClsVsTklA        );
  printf("  SPDClsVsTklB         %f\n", fSPDClsVsTklB        );
  printf("  V0C012vsTklA         %f\n", fV0C012vsTklA        );
  printf("  V0C012vsTklB         %f\n", fV0C012vsTklB        );
  printf("  V0MOnVsOfA           %f\n", fV0MOnVsOfA          );
  printf("  V0MOnVsOfB           %f\n", fV0MOnVsOfB          );
  printf("  SPDOnVsOfA           %f\n", fSPDOnVsOfA          );
  printf("  SPDOnVsOfB           %f\n", fSPDOnVsOfB          );
  printf("  VtxMinContributors   %i\n", fVtxMinContributors  );
  printf("  VtxMinZdist          %f\n", fVtxMinZdist         );
  printf("  VtxNSigmaZdist       %f\n", fVtxNSigmaZdist      );
  printf("  VtxNSigmaDiamXY      %f\n", fVtxNSigmaDiamXY     );
  printf("  VtxNSigmaDiamZ       %f\n", fVtxNSigmaDiamZ      );
  printf("  V0CasymA             %f\n", fV0CasymA            );
  printf("  V0CasymB             %f\n", fV0CasymB            );
  printf("  NBCsPast             %i\n", fNBCsPast            );
  printf("  NBCsFuture           %i\n", fNBCsFuture          );
  printf("  VHMBBAflags          %i\n", fVHMBBAflags         );
  printf("  VHMBBCflags          %i\n", fVHMBBCflags         );
  printf("  VHMBGAflags          %i\n", fVHMBGAflags         );
  printf("  VHMBGCflags          %i\n", fVHMBGCflags         );
  printf("  VIRBBAflags          %i\n", fVIRBBAflags         );
  printf("  VIRBBCflags          %i\n", fVIRBBCflags         );
  printf("  VIRBGAflags          %i\n", fVIRBGAflags         );
  printf("  VIRBGCflags          %i\n", fVIRBGCflags         );
  printf("  V0MOnThreshold       %i\n", fV0MOnThreshold      );
  printf("  V0MOfThreshold       %f\n", fV0MOfThreshold      );
  printf("  SPDGFOThreshold      %i\n", fSPDGFOThreshold     );
  printf("  SH1OuterThreshold    %i\n", fSH1OuterThreshold   );
  printf("  SH2OuterThreshold    %i\n", fSH2OuterThreshold   );
  printf("  TklThreshold         %i\n", fTklThreshold        );
  printf("  TRDptHSE             %f\n", fTRDptHSE            );
  printf("  TRDpidHSE            %i\n", fTRDpidHSE           );
  printf("  TRDptHQU             %f\n", fTRDptHQU            );
  printf("  TRDpidHQU            %i\n", fTRDpidHQU           );
  printf("  TRDptHEE             %f\n", fTRDptHEE            );
  printf("  TRDpidHEE            %i\n", fTRDpidHEE           );
  printf("  TRDminSectorHEE      %i\n", fTRDminSectorHEE     );
  printf("  TRDmaxSectorHEE      %i\n", fTRDmaxSectorHEE     );
  printf("  TRDptHJT             %f\n", fTRDptHJT            );
  printf("  TRDnHJT              %i\n", fTRDnHJT             );
}

