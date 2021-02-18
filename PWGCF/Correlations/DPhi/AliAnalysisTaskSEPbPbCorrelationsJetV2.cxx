#include "AliLog.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TAxis.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"

#include "AliAnalysisTaskSEPbPbCorrelationsJetV2.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"

#include "AliCodeTimer.h"
#include "AliOADBContainer.h"

#include "AliEventPoolManager.h"
#include "AliTHn.h"

ClassImp(AliAnalysisTaskSEPbPbCorrelationsJetV2)
ClassImp(AliBasicParticleST)

 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::cont = NULL;
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2trs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2trs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx3trs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy3trs[14] = { NULL };

TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ2V0A[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ3V0A[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ2V0C[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ3V0C[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ2trk[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ3trk[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ2V0AC[90] = { NULL };
TSpline3 *AliAnalysisTaskSEPbPbCorrelationsJetV2::fSplQ3V0AC[90] = { NULL };

//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::AliAnalysisTaskSEPbPbCorrelationsJetV2() : 
  AliAnalysisTaskSE(), 
  fAOD(0x0),
  fPoolMgr(0x0),
  ftrk(0x0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fNbinsAssocPt(1),
  fNbinsZvtx(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fPtTrigAxis(0x0),
  fPtAssocAxis(0x0),
  fEtaAxis(0x0),
  fZvtxAxis(0x0),
  fRemovePileup(kFALSE),
  fRemovePileup2(kFALSE),
  fRemovePileup3(kFALSE),
  fUseRes(kTRUE),
  fAverageRes(kTRUE),
  fN1(0),
  fN2(-1),
  fLowESEcut(0.),
  fHighESEcut(100.),
  fESEdet(0),
  fEP(kFALSE),
  fMinHardPt(0.),

  fHistACv2(0x0),
  fHistATv2(0x0),
  fHistCTv2(0x0),
  fHistACv3(0x0),
  fHistATv3(0x0),
  fHistCTv3(0x0),
  fHistACv2LowESE(0x0),
  fHistATv2LowESE(0x0),
  fHistCTv2LowESE(0x0),
  fHistACv2HighESE(0x0),
  fHistATv2HighESE(0x0),
  fHistCTv2HighESE(0x0),
  fHistACv3LowESE(0x0),
  fHistATv3LowESE(0x0),
  fHistCTv3LowESE(0x0),
  fHistACv3HighESE(0x0),
  fHistATv3HighESE(0x0),
  fHistCTv3HighESE(0x0),
  fHistACv2vsESE(0x0),
  fHistATv2vsESE(0x0),
  fHistCTv2vsESE(0x0),
  fHistACv3vsESE(0x0),
  fHistATv3vsESE(0x0),
  fHistCTv3vsESE(0x0),
  fHistQxT2(0x0),
  fHistQyT2(0x0),
  fHistQxT3(0x0),
  fHistQyT3(0x0),

  fHistdPhidEtaPt(0x0), 
  fHistdPhidEtaPt_Mixed(0x0), 
  fHistdPhidEtaPt_SS(0x0), 
  fHistdPhidEtaPt_Mixed_SS(0x0), 

  fHistv2ESEvsCent(0x0),
  fHistv3ESEvsCent(0x0), 
  
  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fHistNtrVsCent(0x0),
  fHistV0CVsCent(0x0),
  fHistESDvsTPC(0x0),

  fHistCentVsZ(0x0),
  fHistCentVsZMixed(0x0),
  
  fCentMethod("V0M"),
  fOutputList(0x0),
  fOutputList1(0x0),
  fRunN(-1),
  fv0mult(-1.),
  fv0multonline(-1.),
  fitsmult(-1.),
  nITSTrkls(-1),
  fv0mpercentile(-1.),
  fv0meqpercentile(-1.),
  fv0apercentile(-1.),
  fv0cpercentile(-1.),
  fcl0percentile(-1.),
  fcl1percentile(-1.),
  fzvtx(0.),
  multtrkcut(-1),
  Qx2trkcut(0.),
  Qy2trkcut(0.),
  Qx3trkcut(0.),
  Qy3trkcut(0.),

  sumMa(-1.),
  sumMc(-1.),
  Qxa2(0.),
  Qya2(0.),
  Qxc2(0.),
  Qyc2(0.),
  Qxa3(0.),
  Qya3(0.),
  Qxc3(0.),
  Qyc3(0.),
   
  Qya2Cor(0.),
  Qxa2Cor(0.),
  Qya3Cor(0.),
  Qxa3Cor(0.),
  Qyc2Cor(0.),
  Qxc2Cor(0.),
  Qyc3Cor(0.),
  Qxc3Cor(0.),
  
  Qytr2Cor(0.),
  Qxtr2Cor(0.),
  Qytr3Cor(0.),
  Qxtr3Cor(0.),

  fResACv2(0x0),
  fResATv2(0x0),
  fResCTv2(0x0),
  fResACv2LowESE(0x0),
  fResATv2LowESE(0x0),
  fResCTv2LowESE(0x0),
  fResACv2HighESE(0x0),
  fResATv2HighESE(0x0),
  fResCTv2HighESE(0x0),
  fResACv3(0x0),
  fResATv3(0x0),
  fResCTv3(0x0),
  fResACv3LowESE(0x0),
  fResATv3LowESE(0x0),
  fResCTv3LowESE(0x0),
  fResACv3HighESE(0x0),
  fResATv3HighESE(0x0),
  fResCTv3HighESE(0x0),

  fLowCenCut(0x0),
  fHighCenCut(0x0),
  fV0MultOfOnCut(0x0),
  
  fMultV0(0x0)
  
{

  // Default constructor
  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for(Int_t iZvtx = 0; iZvtx<fNMaxBinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
      fHistSP2A[iCent][iZvtx][index] = NULL;
      fHistSP2C[iCent][iZvtx][index] = NULL;
      fHistSP2T[iCent][iZvtx][index] = NULL;
      fHistSP2ALowESE[iCent][iZvtx][index] = NULL;
      fHistSP2CLowESE[iCent][iZvtx][index] = NULL;
      fHistSP2TLowESE[iCent][iZvtx][index] = NULL;
      fHistSP2AHighESE[iCent][iZvtx][index] = NULL;
      fHistSP2CHighESE[iCent][iZvtx][index] = NULL;
      fHistSP2THighESE[iCent][iZvtx][index] = NULL;
      fHistSP3A[iCent][iZvtx][index] = NULL;
      fHistSP3C[iCent][iZvtx][index] = NULL;
      fHistSP3T[iCent][iZvtx][index] = NULL;
      fHistSP3ALowESE[iCent][iZvtx][index] = NULL;
      fHistSP3CLowESE[iCent][iZvtx][index] = NULL;
      fHistSP3TLowESE[iCent][iZvtx][index] = NULL;
      fHistSP3AHighESE[iCent][iZvtx][index] = NULL;
      fHistSP3CHighESE[iCent][iZvtx][index] = NULL;
      fHistSP3THighESE[iCent][iZvtx][index] = NULL;

      fHistSP2AvsESE[iCent][iZvtx][index] = NULL;
      fHistSP2CvsESE[iCent][iZvtx][index] = NULL;
      fHistSP2TvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3AvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3CvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3TvsESE[iCent][iZvtx][index] = NULL;
      }
      for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt) {
	for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
//          fHistdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP3AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaMixed[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaSSMixed[iCent][iZvtx][ipt][jpt] = NULL;
	}
      }
    }
  }

/*
  for(Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++)
  {
   for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt) 
   {
    for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt) 
    {
     fHistdPhidEta[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaSS[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaMixed[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaSSMixed[iCent][ipt][jpt] = NULL;
    }
   }
  }
*/

  for(Int_t i = 0; i < 14; ++i) {
    fQx2mV0A[i] = fQy2mV0A[i] = fQx2sV0A[i] = fQy2sV0A[i] = NULL;
    fQx2mV0C[i] = fQy2mV0C[i] = fQx2sV0C[i] = fQy2sV0C[i] = NULL;
    fQx2mTrk[i] = fQy2mTrk[i] = fQx2sTrk[i] = fQy2sTrk[i] = NULL;
    fQx3mV0A[i] = fQy3mV0A[i] = fQx3sV0A[i] = fQy3sV0A[i] = NULL;
    fQx3mV0C[i] = fQy3mV0C[i] = fQx3sV0C[i] = fQy3sV0C[i] = NULL;
    fQx3mTrk[i] = fQy3mTrk[i] = fQx3sTrk[i] = fQy3sTrk[i] = NULL;
  }

  // for(Int_t i = 0; i < 90; ++i) {
  //   fSplQ2V0A[i] = fSplQ3V0A[i] = fSplQ2V0C[i] = fSplQ3V0C[i] = fSplQ2trk[i] = fSplQ3trk[i] = fSplQ2V0AC[i] = fSplQ3V0AC[i] = NULL;
  // }

  for(Int_t i = 0; i < 5; ++i) {
    fResACv2vsESE[i] = fResATv2vsESE[i] = fResCTv2vsESE[i] = NULL;
    fResACv3vsESE[i] = fResATv3vsESE[i] = fResCTv3vsESE[i] = NULL;
  }

}


//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::AliAnalysisTaskSEPbPbCorrelationsJetV2(const char *name) : 
  AliAnalysisTaskSE(name), 
  fAOD(0x0),
  fPoolMgr(0x0),
  ftrk(0x0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fNbinsAssocPt(1),
  fNbinsZvtx(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fPtTrigAxis(0x0),
  fPtAssocAxis(0x0),
  fEtaAxis(0x0),
  fZvtxAxis(0x0),
  fRemovePileup(kFALSE),
  fRemovePileup2(kFALSE),
  fRemovePileup3(kFALSE),
  fUseRes(kTRUE),
  fAverageRes(kTRUE),
  fN1(0),
  fN2(-1),
  fLowESEcut(0.),
  fHighESEcut(100.),
  fESEdet(0),
  fEP(kFALSE),
  fMinHardPt(0.),

  fHistACv2(0x0),
  fHistATv2(0x0),
  fHistCTv2(0x0),
  fHistACv3(0x0),
  fHistATv3(0x0),
  fHistCTv3(0x0),
  fHistACv2LowESE(0x0),
  fHistATv2LowESE(0x0),
  fHistCTv2LowESE(0x0),
  fHistACv2HighESE(0x0),
  fHistATv2HighESE(0x0),
  fHistCTv2HighESE(0x0),
  fHistACv3LowESE(0x0),
  fHistATv3LowESE(0x0),
  fHistCTv3LowESE(0x0),
  fHistACv3HighESE(0x0),
  fHistATv3HighESE(0x0),
  fHistCTv3HighESE(0x0),
  fHistACv2vsESE(0x0),
  fHistATv2vsESE(0x0),
  fHistCTv2vsESE(0x0),
  fHistACv3vsESE(0x0),
  fHistATv3vsESE(0x0),
  fHistCTv3vsESE(0x0),
  fHistQxT2(0x0),
  fHistQyT2(0x0),
  fHistQxT3(0x0),
  fHistQyT3(0x0),
 
  fHistdPhidEtaPt(0x0),
  fHistdPhidEtaPt_Mixed(0x0),
  fHistdPhidEtaPt_SS(0x0),
  fHistdPhidEtaPt_Mixed_SS(0x0),
 
  fHistv2ESEvsCent(0x0),
  fHistv3ESEvsCent(0x0), 
  
  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fHistNtrVsCent(0x0),
  fHistV0CVsCent(0x0),
  fHistESDvsTPC(0x0),

  fHistCentVsZ(0x0),
  fHistCentVsZMixed(0x0),
  
  fCentMethod("V0M"),
  fOutputList(0x0),
  fOutputList1(0x0),
  fRunN(-1),
  fv0mult(-1.),
  fv0multonline(-1.),
  fitsmult(-1.),
  nITSTrkls(-1),
  fv0mpercentile(-1.),
  fv0meqpercentile(-1.),
  fv0apercentile(-1.),
  fv0cpercentile(-1.),
  fcl0percentile(-1.),
  fcl1percentile(-1.),
  fzvtx(0.),
  multtrkcut(-1),
  Qx2trkcut(0.),
  Qy2trkcut(0.),
  Qx3trkcut(0.),
  Qy3trkcut(0.),

  sumMa(-1.),
  sumMc(-1.),
  Qxa2(0.),
  Qya2(0.),
  Qxc2(0.),
  Qyc2(0.),
  Qxa3(0.),
  Qya3(0.),
  Qxc3(0.),
  Qyc3(0.),
   
  Qya2Cor(0.),
  Qxa2Cor(0.),
  Qya3Cor(0.),
  Qxa3Cor(0.),
  Qyc2Cor(0.),
  Qxc2Cor(0.),
  Qyc3Cor(0.),
  Qxc3Cor(0.),
  
  Qytr2Cor(0.),
  Qxtr2Cor(0.),
  Qytr3Cor(0.),
  Qxtr3Cor(0.),

  fResACv2(0x0),
  fResATv2(0x0),
  fResCTv2(0x0),
  fResACv2LowESE(0x0),
  fResATv2LowESE(0x0),
  fResCTv2LowESE(0x0),
  fResACv2HighESE(0x0),
  fResATv2HighESE(0x0),
  fResCTv2HighESE(0x0),
  fResACv3(0x0),
  fResATv3(0x0),
  fResCTv3(0x0),
  fResACv3LowESE(0x0),
  fResATv3LowESE(0x0),
  fResCTv3LowESE(0x0),
  fResACv3HighESE(0x0),
  fResATv3HighESE(0x0),
  fResCTv3HighESE(0x0),

  fLowCenCut(0x0),
  fHighCenCut(0x0),
  fV0MultOfOnCut(0x0),
  
  fMultV0(0x0)

{

  // Constructor
  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for(Int_t iZvtx = 0; iZvtx<fNMaxBinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
      fHistSP2A[iCent][iZvtx][index] = NULL;
      fHistSP2C[iCent][iZvtx][index] = NULL;
      fHistSP2T[iCent][iZvtx][index] = NULL;
      fHistSP2ALowESE[iCent][iZvtx][index] = NULL;
      fHistSP2CLowESE[iCent][iZvtx][index] = NULL;
      fHistSP2TLowESE[iCent][iZvtx][index] = NULL;
      fHistSP2AHighESE[iCent][iZvtx][index] = NULL;
      fHistSP2CHighESE[iCent][iZvtx][index] = NULL;
      fHistSP2THighESE[iCent][iZvtx][index] = NULL;
      fHistSP3A[iCent][iZvtx][index] = NULL;
      fHistSP3C[iCent][iZvtx][index] = NULL;
      fHistSP3T[iCent][iZvtx][index] = NULL;
      fHistSP3ALowESE[iCent][iZvtx][index] = NULL;
      fHistSP3CLowESE[iCent][iZvtx][index] = NULL;
      fHistSP3TLowESE[iCent][iZvtx][index] = NULL;
      fHistSP3AHighESE[iCent][iZvtx][index] = NULL;
      fHistSP3CHighESE[iCent][iZvtx][index] = NULL;
      fHistSP3THighESE[iCent][iZvtx][index] = NULL;

      fHistSP2AvsESE[iCent][iZvtx][index] = NULL;
      fHistSP2CvsESE[iCent][iZvtx][index] = NULL;
      fHistSP2TvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3AvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3CvsESE[iCent][iZvtx][index] = NULL;
      fHistSP3TvsESE[iCent][iZvtx][index] = NULL;
      }
      for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt) {
	for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP3AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaMixed[iCent][iZvtx][ipt][jpt] = NULL;
//	  fHistdPhidEtaSSMixed[iCent][iZvtx][ipt][jpt] = NULL;
	}
      }
    }
  }
/*
  for(Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++)
  {
   for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt)
   {
    for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt)
    {
     fHistdPhidEta[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaSS[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaMixed[iCent][ipt][jpt] = NULL;
     fHistdPhidEtaSSMixed[iCent][ipt][jpt] = NULL;
    }
   }
  }
*/


  for(Int_t i = 0; i < 14; ++i) {
    fQx2mV0A[i] = fQy2mV0A[i] = fQx2sV0A[i] = fQy2sV0A[i] = NULL;
    fQx2mV0C[i] = fQy2mV0C[i] = fQx2sV0C[i] = fQy2sV0C[i] = NULL;
    fQx2mTrk[i] = fQy2mTrk[i] = fQx2sTrk[i] = fQy2sTrk[i] = NULL;
    fQx3mV0A[i] = fQy3mV0A[i] = fQx3sV0A[i] = fQy3sV0A[i] = NULL;
    fQx3mV0C[i] = fQy3mV0C[i] = fQx3sV0C[i] = fQy3sV0C[i] = NULL;
    fQx3mTrk[i] = fQy3mTrk[i] = fQx3sTrk[i] = fQy3sTrk[i] = NULL;
  }

  // for(Int_t i = 0; i < 90; ++i) {
  //   fSplQ2V0A[i] = fSplQ3V0A[i] = fSplQ2V0C[i] = fSplQ3V0C[i] = fSplQ2trk[i] = fSplQ3trk[i] = fSplQ2V0AC[i] = fSplQ3V0AC[i] = NULL;
  // }    

  for(Int_t i = 0; i < 5; ++i) {
    fResACv2vsESE[i] = fResATv2vsESE[i] = fResCTv2vsESE[i] = NULL;
    fResACv3vsESE[i] = fResATv3vsESE[i] = fResCTv3vsESE[i] = NULL;
  }

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::~AliAnalysisTaskSEPbPbCorrelationsJetV2() {
  
  if (fOutputList  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
    delete fOutputList;

  if (fOutputList1  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    delete fOutputList1;
}

//___________________________________________________________________________
void AliAnalysisTaskSEPbPbCorrelationsJetV2::NotifyRun()
{
  /// Notify run
  printf("Loading calibration for run %d\n",fInputHandler->GetEvent()->GetRunNumber());
  OpenInfoCalbration(fInputHandler->GetEvent()->GetRunNumber());
}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::UserCreateOutputObjects() {
  
  //=====================================================================
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fOutputList1 = new TList();
  fOutputList1->SetOwner(kTRUE);

    for (Int_t iCent=0; iCent<fNbinsCent; iCent++) {
      for(Int_t iZvtx = 0; iZvtx < fNbinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
	fHistSP2A[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2A_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
	fOutputList -> Add(fHistSP2A[iCent][iZvtx][index]);

	fHistSP2C[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2C_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2C[iCent][iZvtx][index]);
	
	fHistSP2T[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2T_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2T[iCent][iZvtx][index]);

	fHistSP2ALowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2ALowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2ALowESE[iCent][iZvtx][index]);
	
	fHistSP2CLowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2CLowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2CLowESE[iCent][iZvtx][index]);

	fHistSP2TLowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2TLowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2TLowESE[iCent][iZvtx][index]);

	fHistSP2AHighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2AHighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2AHighESE[iCent][iZvtx][index]);
	
	fHistSP2CHighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2CHighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2CHighESE[iCent][iZvtx][index]);

	fHistSP2THighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2THighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP2THighESE[iCent][iZvtx][index]);

	fHistSP3A[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3A_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3A[iCent][iZvtx][index]);
	
	fHistSP3C[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3C_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3C[iCent][iZvtx][index]);

	fHistSP3T[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3T_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3T[iCent][iZvtx][index]);

	fHistSP3ALowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3ALowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3ALowESE[iCent][iZvtx][index]);
	
	fHistSP3CLowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3CLowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3CLowESE[iCent][iZvtx][index]);

	fHistSP3TLowESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3TLowESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							    Form("%d-%d %%, %.1f<z<%.1f",
								 Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								 Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								 fZvtxAxis->GetBinLowEdge(iZvtx+1),
								 fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							    fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3TLowESE[iCent][iZvtx][index]);

	fHistSP3AHighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3AHighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3AHighESE[iCent][iZvtx][index]);
	
	fHistSP3CHighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3CHighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3CHighESE[iCent][iZvtx][index]);

	fHistSP3THighESE[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP3THighESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray());
//	fOutputList -> Add(fHistSP3THighESE[iCent][iZvtx][index]);

	// Histos versus ESE
	fHistSP2AvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP2AvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP2AvsESE[iCent][iZvtx][index]);
	
	fHistSP2CvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP2CvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP2CvsESE[iCent][iZvtx][index]);
	
	fHistSP2TvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP2TvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP2TvsESE[iCent][iZvtx][index]);
	////////////////////

	fHistSP3AvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP3AvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP3AvsESE[iCent][iZvtx][index]);
	
	fHistSP3CvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP3CvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP3CvsESE[iCent][iZvtx][index]);
	
	fHistSP3TvsESE[iCent][iZvtx][index] = new TProfile2D(Form("fHist%sSP3TvsESE_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
							     Form("%d-%d %%, %.1f<z<%.1f",
								  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								  fZvtxAxis->GetBinLowEdge(iZvtx+1),
								  fZvtxAxis->GetBinUpEdge(iZvtx+1)),
							     fNbinsPt,fPtAxis->GetXbins()->GetArray(),
							     5,0,100);
//	fOutputList -> Add(fHistSP3TvsESE[iCent][iZvtx][index]);
      }
      for(Int_t iPt = 0; iPt < fNbinsPtTrig; ++iPt) {
	for(Int_t jPt = 0; jPt < fNbinsAssocPt; ++jPt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2AdPhidEta_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								Form("%d-%d %%, %.1f<z<%.1f",
								     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								     fZvtxAxis->GetBinLowEdge(iZvtx+1),
								     fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
	  fOutputList1 -> Add(fHistSP2AdPhidEta[iCent][iZvtx][iPt][jPt]);
	  
	  if (0) { // for output file size reasons
	  fHistSP2CdPhidEta[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2CdPhidEta_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								Form("%d-%d %%, %.1f<z<%.1f",
								     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								     fZvtxAxis->GetBinLowEdge(iZvtx+1),
								     fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
//	  fOutputList -> Add(fHistSP2CdPhidEta[iCent][iZvtx][iPt][jPt]);
	  
	  fHistSP2TdPhidEta[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2TdPhidEta_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								Form("%d-%d %%, %.1f<z<%.1f",
								     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								     fZvtxAxis->GetBinLowEdge(iZvtx+1),
								     fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
//	  fOutputList -> Add(fHistSP2TdPhidEta[iCent][iZvtx][iPt][jPt]);
	  }
	  	  
	  // same-sign track pairs
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2AdPhidEtaSS_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								  Form("%d-%d %%, %.1f<z<%.1f",
								       Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								       fZvtxAxis->GetBinLowEdge(iZvtx+1),
								       fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								  24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
	  fOutputList1 -> Add(fHistSP2AdPhidEtaSS[iCent][iZvtx][iPt][jPt]);
	  
	  if (0) { // for output file size reasons
	  fHistSP2CdPhidEtaSS[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2CdPhidEtaSS_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								  Form("%d-%d %%, %.1f<z<%.1f",
								       Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								       fZvtxAxis->GetBinLowEdge(iZvtx+1),
								       fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								  24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
//	  fOutputList -> Add(fHistSP2CdPhidEtaSS[iCent][iZvtx][iPt][jPt]);

	  fHistSP2TdPhidEtaSS[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2TdPhidEtaSS_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								  Form("%d-%d %%, %.1f<z<%.1f",
								       Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								       fZvtxAxis->GetBinLowEdge(iZvtx+1),
								       fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								  24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
//	  fOutputList -> Add(fHistSP2TdPhidEtaSS[iCent][iZvtx][iPt][jPt]);
	  }	
	  

	  fHistSP3AdPhidEtaSS[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP3AdPhidEtaSS_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								  Form("%d-%d %%, %.1f<z<%.1f",
								       Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								       fZvtxAxis->GetBinLowEdge(iZvtx+1),
								       fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								  24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
//	  fOutputList -> Add(fHistSP3AdPhidEtaSS[iCent][iZvtx][iPt][jPt]);
	}
      }
     }
    }

 const Int_t nbins_dPhidEtaPt[] = {24, 24, 20, fNbinsPtTrig, fNbinsAssocPt, fNbinsCent};
 const Int_t nVar = sizeof(nbins_dPhidEtaPt) / sizeof(Int_t); 

 const TArrayD *aBin_Trig(fPtTrigAxis->GetXbins());
 const Int_t nBin_Trig(aBin_Trig->GetSize() - 1);
 Double_t dBin_Trig[nBin_Trig+1];

 const TArrayD *aBin_Asso(fPtAssocAxis->GetXbins());
 const Int_t nBin_Asso(aBin_Asso->GetSize() - 1);
 Double_t dBin_Asso[nBin_Asso+1];

 const TArrayD *aBin_Cen(fCentAxis->GetXbins());
 const Int_t nBin_Cen(aBin_Cen->GetSize() - 1);
 Double_t dBin_Cen[nBin_Cen+1];

 for (Int_t i=0; i<=nBin_Trig; ++i) dBin_Trig[i] = (*aBin_Trig)[i];
 for (Int_t i=0; i<=nBin_Asso; ++i) dBin_Asso[i] = (*aBin_Asso)[i];
 for (Int_t i=0; i<=nBin_Cen;  ++i) dBin_Cen[i]  = (*aBin_Cen)[i];


 fHistdPhidEtaPt = new AliTHn("fHistdPhidEtaPt", "fHistdPhidEtaPt", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt->SetBinLimits(5, dBin_Cen);
 
 fHistdPhidEtaPt->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt->SetVarTitle(4, "asso p_{T} GeV/c");
 fHistdPhidEtaPt->SetVarTitle(5, "centrality");

 fHistdPhidEtaPt_SS = new AliTHn("fHistdPhidEtaPt_SS", "fHistdPhidEtaPt_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_SS->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_SS->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_SS->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_SS->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_SS->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_SS->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_SS->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_SS->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_SS->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_SS->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_SS->SetVarTitle(4, "asso p_{T} GeV/c");
 fHistdPhidEtaPt_SS->SetVarTitle(5, "centrality");
 
 fHistdPhidEtaPt_Mixed = new AliTHn("fHistdPhidEtaPt_Mixed", "fHistdPhidEtaPt_Mixed", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_Mixed->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_Mixed->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_Mixed->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_Mixed->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_Mixed->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_Mixed->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_Mixed->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_Mixed->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_Mixed->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_Mixed->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_Mixed->SetVarTitle(4, "asso p_{T} GeV/c"); 
 fHistdPhidEtaPt_Mixed->SetVarTitle(5, "centrality"); 

 fHistdPhidEtaPt_Mixed_SS = new AliTHn("fHistdPhidEtaPt_Mixed_SS", "fHistdPhidEtaPt_Mixed_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(4, "asso p_{T} GeV/c"); 
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(5, "centrality"); 

 
 const Int_t nbins_dTrig[] = {fNbinsPtTrig, 20, fNbinsCent};
 const Int_t nVar_Trig = sizeof(nbins_dTrig) / sizeof(Int_t);

 fHistTrig = new AliTHn("fHistTrig", "fHistTrig", 1, nVar_Trig, nbins_dTrig);
 fHistTrig->SetBinLimits(0, dBin_Trig);
 fHistTrig->SetBinLimits(1, -10, 10);
 fHistTrig->SetBinLimits(2, dBin_Cen);
 fHistTrig->SetVarTitle(0, "leading p_{T} GeV/c");
 fHistTrig->SetVarTitle(1, "Vz");
 fHistTrig->SetVarTitle(2, "centrality");


 fOutputList1->Add(fHistdPhidEtaPt); 
 fOutputList1->Add(fHistdPhidEtaPt_SS); 
 fOutputList1->Add(fHistdPhidEtaPt_Mixed); 
 fOutputList1->Add(fHistdPhidEtaPt_Mixed_SS); 
 fOutputList1->Add(fHistTrig);


  fHistACv2 = new TProfile("fHistACv2", "; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{V0C}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistACv2);
  fHistATv2 = new TProfile("fHistATv2", "; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistATv2);
  fHistCTv2 = new TProfile("fHistCTv2", "; centrality percentile; <Q^{2}_{V0C}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistCTv2);
  fHistACv3 = new TProfile("fHistACv3", "; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{V0C}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistACv3);
  fHistATv3 = new TProfile("fHistATv3", "; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistATv3);
  fHistCTv3 = new TProfile("fHistCTv3", "; centrality percentile; <Q^{3}_{V0C}.Q*^{3}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistCTv3);

  fHistACv2LowESE = new TProfile("fHistACv2LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{V0C}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistACv2LowESE);
  fHistATv2LowESE = new TProfile("fHistATv2LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{Trk}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistATv2LowESE);
  fHistCTv2LowESE = new TProfile("fHistCTv2LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{2}_{V0C}.Q*^{2}_{Trk}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistCTv2LowESE);

  fHistACv2HighESE = new TProfile("fHistACv2HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{V0C}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistACv2HighESE);
  fHistATv2HighESE = new TProfile("fHistATv2HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{Trk}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistATv2HighESE);
  fHistCTv2HighESE = new TProfile("fHistCTv2HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{2}_{V0C}.Q*^{2}_{Trk}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistCTv2HighESE);

  fHistACv3LowESE = new TProfile("fHistACv3LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{V0C}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistACv3LowESE);
  fHistATv3LowESE = new TProfile("fHistATv3LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{Trk}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistATv3LowESE);
  fHistCTv3LowESE = new TProfile("fHistCTv3LowESE", Form("%.1f%% Low ESE; centrality percentile; <Q^{3}_{V0C}.Q*^{3}_{Trk}>",fLowESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistCTv3LowESE);

  fHistACv3HighESE = new TProfile("fHistACv3HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{V0C}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistACv3HighESE);
  fHistATv3HighESE = new TProfile("fHistATv3HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{3}_{V0A}.Q*^{3}_{Trk}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistATv3HighESE);
  fHistCTv3HighESE = new TProfile("fHistCTv3HighESE", Form("%.1f%% High ESE; centrality percentile; <Q^{3}_{V0C}.Q*^{3}_{Trk}>",fHighESEcut), fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
//  fOutputList->Add(fHistCTv3HighESE);

  fHistACv2vsESE = new TProfile2D("fHistACv2vsESE", "<Q^{2}_{V0A}.Q*^{2}_{V0C}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistACv2vsESE);
  fHistATv2vsESE = new TProfile2D("fHistATv2vsESE", "<Q^{2}_{V0A}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistATv2vsESE);
  fHistCTv2vsESE = new TProfile2D("fHistCTv2vsESE", "<Q^{2}_{V0C}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistCTv2vsESE);
  fHistACv3vsESE = new TProfile2D("fHistACv3vsESE", "<Q^{3}_{V0A}.Q*^{3}_{V0C}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistACv3vsESE);
  fHistATv3vsESE = new TProfile2D("fHistATv3vsESE", "<Q^{3}_{V0A}.Q*^{3}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistATv3vsESE);
  fHistCTv3vsESE = new TProfile2D("fHistCTv3vsESE", "<Q^{3}_{V0C}.Q*^{3}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), 5, 0, 100);
//  fOutputList->Add(fHistCTv3vsESE);

  fHistQxT2 = new TProfile("fHistQxT2","Qx2 SPD",fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistQxT2);
  fHistQyT2 = new TProfile("fHistQyT2","Qy2 SPD",fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistQyT2);
  fHistQxT3 = new TProfile("fHistQxT3","Qx3 SPD",fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistQxT3);
  fHistQyT3 = new TProfile("fHistQyT3","Qy3 SPD",fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistQyT3);

  fHistv2ESEvsCent = new TH2D("fHistv2ESEvsCent","",
			      fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			      10,0,100);
//  fOutputList->Add(fHistv2ESEvsCent);
  fHistv3ESEvsCent = new TH2D("fHistv3ESEvsCent","",
			      fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			      10,0,100);
//  fOutputList->Add(fHistv3ESEvsCent);

  TString fFuncv2FN;
  if (fESEdet == 0)
    fFuncv2FN = "V0A";
  else if (fESEdet == 1)
    fFuncv2FN = "V0C";
  else if (fESEdet == 2)
    fFuncv2FN = "SPD";
  else if (fESEdet == 3)
    fFuncv2FN = "V0AC";
  
  fHistV0Multiplicity = new TH2D("fHistV0Multiplicity", "V0 Multiplicity (Online vs Offline)",
				 500, 0, 50000,
				 500, 0, 50000);
  fHistV0Multiplicity -> SetXTitle("Multiplicity Offline");
  fHistV0Multiplicity -> SetYTitle("Multiplicity Online");
  fHistV0Multiplicity -> Sumw2();

  fHistITSMultiplicity = new TH1D("fHistITSMultiplicity", "ITS Multiplicity", 5000, 0, 25000);
  fHistITSMultiplicity -> SetXTitle("N_{Clusters}");
  fHistITSMultiplicity -> Sumw2();

  fHistCentrality = new TH1D("fHistCentrality", Form("%s_Centrality",fCentMethod.Data()), 400, -100, 300);
  fHistCentrality -> SetXTitle("Centrality  [%]");
  fHistCentrality -> Sumw2();


  fPileup1_Before = new TH2D("fPileup1_Before", "Pile_up_before (Online vs cut)",
                                 500, 0, 50000,
                                 500, 0, 50000);
  fPileup1_Before -> SetXTitle("Multiplicity Online Cut");
  fPileup1_Before -> SetYTitle("Multiplicity Online");
  fPileup1_Before -> Sumw2();

  fPileup1_After = new TH2D("fPileup1_After", "Pile_up_After (Online vs cut)",
                                 500, 0, 50000,
                                 500, 0, 50000);
  fPileup1_After -> SetXTitle("Multiplicity Online Cut");
  fPileup1_After -> SetYTitle("Multiplicity Online");
  fPileup1_After -> Sumw2();

  fPileup2_Before = new TH2D("fPileup2_Before", "ITS Multiplicity vs Num Tracklets", 5000, 0, 25000, 2000, 0 , 10000);
  fPileup2_Before -> SetXTitle("N_{Clusters}");
  fPileup2_Before -> SetYTitle("N_{Tracklets}");
  fPileup2_Before -> Sumw2();
  
  fPileup2_After = new TH2D("fPileup2_After", "ITS Multiplicity vs Num Tracklets", 5000, 0, 25000, 2000, 0 , 10000);
  fPileup2_After -> SetXTitle("N_{Clusters}");
  fPileup2_After -> SetYTitle("N_{Tracklets}");
  fPileup2_After -> Sumw2();

  fPileup3_Before_Low = new TH2D("fPileup3_Before_Low", "CL0 vs Low Event Cut", 100, 0, 100, 100, 0 , 100);
  fPileup3_Before_Low -> SetXTitle("Centrality  [%]");
  fPileup3_Before_Low -> SetYTitle("Centrality  [%]");
  fPileup3_Before_Low -> Sumw2();

  fPileup3_After_Low = new TH2D("fPileup3_After_Low", "CL0 vs Low Event Cut", 100, 0, 100, 100, 0 , 100);
  fPileup3_After_Low -> SetXTitle("Centrality  [%]");
  fPileup3_After_Low -> SetYTitle("Centrality  [%]");
  fPileup3_After_Low -> Sumw2();

  fPileup3_Before_High = new TH2D("fPileup3_Before_High", "CL0 vs High Event Cut", 100, 0, 100, 100, 0 , 100);
  fPileup3_Before_High -> SetXTitle("Centrality  [%]");
  fPileup3_Before_High -> SetYTitle("Centrality  [%]");
  fPileup3_Before_High -> Sumw2();

  fPileup3_After_High = new TH2D("fPileup3_After_High", "CL0 vs High Event Cut", 100, 0, 100, 100, 0 , 100);
  fPileup3_After_High -> SetXTitle("Centrality  [%]");
  fPileup3_After_High -> SetYTitle("Centrality  [%]");
  fPileup3_After_High -> Sumw2();


  fOutputList -> Add(fHistV0Multiplicity);
  fOutputList -> Add(fHistITSMultiplicity);
  fOutputList -> Add(fHistCentrality);

  fOutputList -> Add(fPileup1_Before);
  fOutputList -> Add(fPileup1_After);

  fOutputList -> Add(fPileup2_Before);
  fOutputList -> Add(fPileup2_After);

  fOutputList -> Add(fPileup3_Before_Low);
  fOutputList -> Add(fPileup3_After_Low);

  fOutputList -> Add(fPileup3_Before_High);
  fOutputList -> Add(fPileup3_After_High);



  fHistEvStat = new TH1D("fHistEvStat","Event cuts statistics",25,-0.5,24.5);
  fHistEvStat->SetXTitle("Cut index");
  fOutputList->Add(fHistEvStat);

  fHistNtrVsCent = new TH2D("fHistNtrVsCent","Number of tracks vs centrality",
			    fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			    2000,0.,20000.);
  fOutputList->Add(fHistNtrVsCent);

  fHistV0CVsCent = new TH2D("fHistV0CVsCent","V0C vs centrality",
			    fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			    500,0.,50000.);
  fOutputList->Add(fHistV0CVsCent);

  fHistESDvsTPC = new TH2D("fHistESDvsTPC","ESD vs TPC multiplicity",
			   500,0.,20000.,
			   500,0.,20000.);
  fOutputList->Add(fHistESDvsTPC);

  fHistCentVsZ = new TH2D("fHistCentVsZ","",20,-10,10,20,0,100);
  fOutputList->Add(fHistCentVsZ);
  fHistCentVsZMixed = new TH2D("fHistCentVsZMixed","Mixed events",20,-10,10,20,0,100);
  fOutputList->Add(fHistCentVsZMixed);
 
  PostData(1, fOutputList); 
  PostData(2, fOutputList1); 

  //================================================================
  if (fESEdet == 2) {
    TFile *fSpl2 = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2TrklNoEtaCutRun2.root");
    if (!fSpl2)
      fSpl2 = TFile::Open("./calibSpq2TrklNoEtaCutRun2.root");
  
    if (!fSpl2) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl2 = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2TrklNoEtaCutRun2.root");
    }

    if (!fSpl2){
      printf("Spline file q2 cannot be opened \n");
      return;
    }

    TFile *fSpl3 = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3TrklNoEtaCutRun2.root");
    if (!fSpl3)
      fSpl3 = TFile::Open("./calibSpq3TrklNoEtaCutRun2.root");
  
    if (!fSpl3) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
    
      fSpl3 = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3TrklNoEtaCutRun2.root");
    }

    if (!fSpl3){
      printf("Spline file q3 cannot be opened \n");
      return;
    }
    for (Int_t isp = 0; isp < 90; isp++) {
      if (!fSplQ2trk[isp]) fSplQ2trk[isp] = (TSpline3*)fSpl2->Get(Form("sp_q2Trkl_%d", isp));
      if (!fSplQ3trk[isp]) fSplQ3trk[isp] = (TSpline3*)fSpl3->Get(Form("sp_q3Trkl_%d", isp));
    }
  }

  if (fESEdet == 0) {
    TFile *fSpl2V0A = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0ARun2.root");
    if (!fSpl2V0A)
      fSpl2V0A = TFile::Open("./calibSpq2V0ARun2.root");
    
    if (!fSpl2V0A) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl2V0A = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0ARun2.root");
    }
    
    if (!fSpl2V0A){
      printf("Spline file V0A q2 cannot be opened \n");
      return;
    }
    
    TFile *fSpl3V0A = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0ARun2.root");
    if (!fSpl3V0A)
      fSpl3V0A = TFile::Open("./calibSpq3V0ARun2.root");
  
    if (!fSpl3V0A) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl3V0A = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0ARun2.root");
    }

    if (!fSpl3V0A){
      printf("Spline file V0A q3 cannot be opened \n");
      return;
    }
  
    for (Int_t isp = 0; isp < 90; isp++) {
      if (!fSplQ2V0A[isp]) fSplQ2V0A[isp] = (TSpline3*)fSpl2V0A->Get(Form("sp_q2V0A_%d", isp));
      if (!fSplQ3V0A[isp]) fSplQ3V0A[isp] = (TSpline3*)fSpl3V0A->Get(Form("sp_q3V0A_%d", isp));
    }
  }
  if (fESEdet == 1) {
    TFile *fSpl2V0C = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0CRun2.root");
    if (!fSpl2V0C)
      fSpl2V0C = TFile::Open("./calibSpq2V0CRun2.root");
    
    if (!fSpl2V0C) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl2V0C = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0CRun2.root");
    }
    
    if (!fSpl2V0C){
      printf("Spline file V0C q2 cannot be opened \n");
      return;
    }
    
    TFile *fSpl3V0C = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0CRun2.root");
    if (!fSpl3V0C)
      fSpl3V0C = TFile::Open("./calibSpq3V0CRun2.root");
  
    if (!fSpl3V0C) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl3V0C = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0CRun2.root");
    }

    if (!fSpl3V0C){
      printf("Spline file V0C q3 cannot be opened \n");
      return;
    }
  
    for (Int_t isp = 0; isp < 90; isp++) {
      if (!fSplQ2V0C[isp]) fSplQ2V0C[isp] = (TSpline3*)fSpl2V0C->Get(Form("sp_q2V0C_%d", isp));
      if (!fSplQ3V0C[isp]) fSplQ3V0C[isp] = (TSpline3*)fSpl3V0C->Get(Form("sp_q3V0C_%d", isp));
    }
  }
  if (fESEdet == 3) {
    TFile *fSpl2V0AC = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0CombACRun2New.root");
    if (!fSpl2V0AC)
      fSpl2V0AC = TFile::Open("./calibSpq2V0CombACRun2New.root");
    
    if (!fSpl2V0AC) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl2V0AC = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq2V0CombACRun2New.root");
    }
    
    if (!fSpl2V0AC){
      printf("Spline file V0AC q2 cannot be opened \n");
      return;
    }
    
    TFile *fSpl3V0AC = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0CombACRun2New.root");
    if (!fSpl3V0AC)
      fSpl3V0AC = TFile::Open("./calibSpq3V0CombACRun2New.root");
  
    if (!fSpl3V0AC) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }
      
      fSpl3V0AC = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibSpq3V0CombACRun2New.root");
    }

    if (!fSpl3V0AC){
      printf("Spline file V0AC q3 cannot be opened \n");
      return;
    }
  
    for (Int_t isp = 0; isp < 90; isp++) {
      if (!fSplQ2V0AC[isp]) fSplQ2V0AC[isp] = (TSpline3*)fSpl2V0AC->Get(Form("sp_q2V0AC_%d", isp));
      if (!fSplQ3V0AC[isp]) fSplQ3V0AC[isp] = (TSpline3*)fSpl3V0AC->Get(Form("sp_q3V0AC_%d", isp));
    }
  }
  
  if (fUseRes) { // read resolution curves
    TString qVectResFN;
    if (fESEdet == 0)
      qVectResFN = "V0A";
    else if (fESEdet == 1)
      qVectResFN = "V0C";
    else if (fESEdet == 2)
      qVectResFN = "SPD";
    else if (fESEdet == 3)
      qVectResFN = "V0AC";

    if (fEP) qVectResFN += "_EP";
    
    TFile *fFileRes = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/qVectResMuons_%s.root",qVectResFN.Data()));
    if (!fFileRes)
      fFileRes = TFile::Open(Form("./qVectResMuons_%s.root",qVectResFN.Data()));

    if (!fFileRes) {
      if (!gGrid) {
	TGrid::Connect("alien://");
      }

      fFileRes = TFile::Open(Form("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/qVectResMuons_%s.root",qVectResFN.Data()));
    }

    if (fAverageRes) {
      fResACv2 = (TF1*)fFileRes->Get("fResACv2");
      fResATv2 = (TF1*)fFileRes->Get("fResATv2");
      fResCTv2 = (TF1*)fFileRes->Get("fResCTv2");
      fResACv2LowESE = (TF1*)fFileRes->Get("fResACv2LowESE");
      fResATv2LowESE = (TF1*)fFileRes->Get("fResATv2LowESE");
      fResCTv2LowESE = (TF1*)fFileRes->Get("fResCTv2LowESE");
      fResACv2HighESE = (TF1*)fFileRes->Get("fResACv2HighESE");
      fResATv2HighESE = (TF1*)fFileRes->Get("fResATv2HighESE");
      fResCTv2HighESE = (TF1*)fFileRes->Get("fResCTv2HighESE");
      fResACv3 = (TF1*)fFileRes->Get("fResACv3");
      fResATv3 = (TF1*)fFileRes->Get("fResATv3");
      fResCTv3 = (TF1*)fFileRes->Get("fResCTv3");
      fResACv3LowESE = (TF1*)fFileRes->Get("fResACv3LowESE");
      fResATv3LowESE = (TF1*)fFileRes->Get("fResATv3LowESE");
      fResCTv3LowESE = (TF1*)fFileRes->Get("fResCTv3LowESE");
      fResACv3HighESE = (TF1*)fFileRes->Get("fResACv3HighESE");
      fResATv3HighESE = (TF1*)fFileRes->Get("fResATv3HighESE");
      fResCTv3HighESE = (TF1*)fFileRes->Get("fResCTv3HighESE");
    }
    else {
      printf("Using resolution vs ESE...\n");
      for(Int_t i = 0; i < 5; ++i) {
	fResACv2vsESE[i] = (TF1*)fFileRes->Get(Form("fResACv2vsESE_%d",i));
	fResATv2vsESE[i] = (TF1*)fFileRes->Get(Form("fResATv2vsESE_%d",i));
	fResCTv2vsESE[i] = (TF1*)fFileRes->Get(Form("fResCTv2vsESE_%d",i));
	fResACv3vsESE[i] = (TF1*)fFileRes->Get(Form("fResACv3vsESE_%d",i));
	fResATv3vsESE[i] = (TF1*)fFileRes->Get(Form("fResATv3vsESE_%d",i));
	fResCTv3vsESE[i] = (TF1*)fFileRes->Get(Form("fResCTv3vsESE_%d",i));
      }
    }
  }
  
  Double_t parCent[6] = {0.0252185, 0.975157, 0.677622, 0.0290711, -0.000545769, 5.83219e-06};
    
  fHighCenCut = new TF1("fHighCenCut", "[0]+[1]*x + 6*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fHighCenCut->SetParameters(parCent);
    
  fLowCenCut = new TF1("fLowCenCut", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fLowCenCut->SetParameters(parCent);
    
    
  Double_t parV0[8] = {24.819, 0.952938, 0.0775265, 184.972, 6.9353, 0.0433712, -0.000192522, 7.67571e-07};
  fV0MultOfOnCut = new TF1("fV0MultOfOnCut", "[0]+[1]*x - 5.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0MultOfOnCut->SetParameters(parV0);

  ftrk = new TClonesArray("AliBasicParticleST",500);
  
  const Int_t nzvtx = 10;
  Double_t zvtxbins[nzvtx+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
 
  const Int_t nCentralityBins  = 14;
  Double_t centBins[] = {0.,1.,2.,3.,4.,5., 10.,20.,30.,40.,50.,60.,70.,80.,90.};

  fPoolMgr = new AliEventPoolManager(1000, 50000, nCentralityBins, centBins, nzvtx, zvtxbins);
  //fPoolMgr = new AliEventPoolManager(500, 30000, fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), nzvtx, zvtxbins);
  if (!fPoolMgr) return;
  fPoolMgr->SetTargetValues(50000, 0.1, 5);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::UserExec(Option_t *) {

  Int_t cutIndex = 0;
  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) return;  

  fRunN = fAOD->GetRunNumber();
  
  fHistEvStat->Fill(cutIndex++);
  // Trigger selection
  TString trigStr(fAOD->GetFiredTriggerClasses());
  if (trigStr.Contains("CINT7-B-NOPF")) fHistEvStat->Fill(cutIndex);
  cutIndex++;
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  if (trigStr.Contains("CINT7-B-NOPF") &&
      aodV0->GetV0ADecision()==1 &&
      aodV0->GetV0CDecision()==1) fHistEvStat->Fill(cutIndex);
  cutIndex++;
  
  if (!(IsTriggerFired())) return;
  fHistEvStat->Fill(cutIndex++);
  
  fv0mult = GetV0Multiplicity();
  fv0multonline = aodV0->GetTriggerChargeA()+aodV0->GetTriggerChargeC();
  fHistV0Multiplicity  -> Fill(fv0mult,fv0multonline);
  fitsmult = GetITSMultiplicity();
  fHistITSMultiplicity -> Fill(fitsmult);

  Double_t percentile;
  AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
  if(!multSelection) return;
  if(fCentMethod == "") fCentMethod = "V0M";
  //percentile = multSelection->GetMultiplicityPercentile(fCentMethod.Data());
  percentile = multSelection->GetMultiplicityPercentile(fCentMethod);
  fHistCentrality->Fill(percentile);

  fv0mpercentile = multSelection->GetMultiplicityPercentile("V0M");
  fv0meqpercentile = multSelection->GetMultiplicityPercentile("V0MEq");
  fv0apercentile = multSelection->GetMultiplicityPercentile("V0A");
  fv0cpercentile = multSelection->GetMultiplicityPercentile("V0C");
  fcl0percentile = multSelection->GetMultiplicityPercentile("CL0");
  fcl1percentile = multSelection->GetMultiplicityPercentile("CL1");

  Int_t centBin = GetCentBin(percentile);
  if (centBin<0) return;
  fHistEvStat->Fill(cutIndex++);

  if (fAOD->IsIncompleteDAQ()) return;
  fHistEvStat->Fill(cutIndex++);

  // Vertex selection
  const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
  if (trkVtx->GetNContributors() < 2) return;
  fHistEvStat->Fill(cutIndex++);
  Double_t covTrk[6] = {0};
  trkVtx->GetCovarianceMatrix(covTrk);
  Double_t zvtxTrk = trkVtx->GetZ();
  
  const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return;
  fHistEvStat->Fill(cutIndex++);

  fzvtx = spdVtx->GetZ();
  if (fzvtx >= fZvtxAxis->GetXmax()) return;
  if (fzvtx <= fZvtxAxis->GetXmin()) return;
  fHistEvStat->Fill(cutIndex++);

  Double_t dz = zvtxTrk - fzvtx;
  Double_t errTot = TMath::Sqrt(covTrk[5]+cov[5]);
  Double_t errTrk = TMath::Sqrt(covTrk[5]);
  Double_t nsigTot = dz/errTot;
  Double_t nsigTrk = dz/errTrk;
  if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrk)>20)
    return; // bad vertexing
  fHistEvStat->Fill(cutIndex++);
  
  fPileup1_Before->Fill(fv0multonline,fV0MultOfOnCut->Eval(fv0mult));
  if (fRemovePileup) {
    if (fv0multonline < fV0MultOfOnCut->Eval(fv0mult))
      return;
    fPileup1_After->Fill(fv0multonline,fV0MultOfOnCut->Eval(fv0mult));
  }
  fHistEvStat->Fill(cutIndex++);

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  nITSTrkls = aodTrkl->GetNumberOfTracklets();

  fPileup2_Before->Fill(fitsmult,nITSTrkls);
  if (fRemovePileup2) {
    if (fitsmult > (400.+4.*Double_t(nITSTrkls)))
      return;
    fPileup2_After->Fill(fitsmult,nITSTrkls);
  }
  fHistEvStat->Fill(cutIndex++);

  fPileup3_Before_Low->Fill(fcl0percentile, fLowCenCut->Eval(fv0mpercentile));
  fPileup3_Before_High->Fill(fcl0percentile,fHighCenCut->Eval(fv0mpercentile));
  if (fRemovePileup3) {
    if (fcl0percentile < fLowCenCut->Eval(fv0mpercentile))
      return;
    fPileup3_After_Low->Fill(fcl0percentile, fLowCenCut->Eval(fv0mpercentile));
    if (fcl0percentile > fHighCenCut->Eval(fv0mpercentile))
      return;
    fPileup3_After_High->Fill(fcl0percentile,fHighCenCut->Eval(fv0mpercentile));
  }

  fHistEvStat->Fill(cutIndex++);

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks();
  Int_t multTpc = 0;
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if (!aodTrk) continue;
    if (aodTrk->TestFilterBit(128)) multTpc++;
  }
  fHistESDvsTPC->Fill(multTpc,multEsd);
  Double_t multESDTPCDif = Double_t(multEsd) - Double_t(multTpc)*3.38;
  if (multESDTPCDif > 1000.) return;
  fHistEvStat->Fill(cutIndex++);
  
  if (!ComputeQ(fAOD,fzvtx)) return;
  
  fHistEvStat->Fill(cutIndex++);

  Int_t zvtxBin = GetZvtxBin(fzvtx);

  fHistV0CVsCent->Fill(percentile,fAOD->GetVZEROData()->GetMTotV0C());
  
  ProcessEvent(percentile, centBin, zvtxBin);

  PostData(1, fOutputList); 
  PostData(2, fOutputList1); 
}


void AliAnalysisTaskSEPbPbCorrelationsJetV2::ProcessEvent(Double_t percentile, Int_t centBin, Int_t zvtxBin)
{
  AliCodeTimerAuto("",0);

  ftrk->Clear();
  
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  fHistNtrVsCent->Fill(percentile,nTracks);

  // Calculate the q vect resolutions
  Double_t resA2=1.,resC2=1.,resT2=1.,resA2LowESE=1.,resC2LowESE=1.,resT2LowESE=1.,resA2HighESE=1.,resC2HighESE=1.,resT2HighESE=1.; 
  Double_t resA3=1.,resC3=1.,resT3=1.,resA3LowESE=1.,resC3LowESE=1.,resT3LowESE=1.,resA3HighESE=1.,resC3HighESE=1.,resT3HighESE=1.; 
  if (fUseRes) CalcResolutions(percentile,
			       resA2, resC2, resT2,
			       resA2LowESE, resC2LowESE, resT2LowESE,
			       resA2HighESE, resC2HighESE, resT2HighESE,
			       resA3, resC3, resT3,
			       resA3LowESE, resC3LowESE, resT3LowESE,
			       resA3HighESE, resC3HighESE, resT3HighESE);

  // check if mixed event pool is ready
  AliEventPool* pool = fPoolMgr->GetEventPool(percentile,fzvtx);
  Bool_t poolReady = (pool->IsReady() || pool->NTracksInPool() > 5000 || pool->GetCurrentNEvents() > 5);

  fHistCentVsZ->Fill(fzvtx,percentile);
  if (poolReady) fHistCentVsZMixed->Fill(fzvtx,percentile,pool->GetCurrentNEvents());

  
  for (Int_t iTr=0; iTr<nTracks; iTr++) {
    AliAODTrack *track = (AliAODTrack*) fAOD->GetTrack(iTr);
    if (!track) continue;
    if (track->TestFilterBit(768) && TMath::Abs(track->Eta()) < 0.8 && track->GetTPCNcls() >= 70) {
      // same event
      FillHistogramsdPhidEta(iTr,track->Charge(),track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
			     resA2, resC2, resT2,
			     resA3, resC3, resT3);
      
      // Add tracklet to the array for mixing pool
      new ((*ftrk)[ftrk->GetEntriesFast()]) AliBasicParticleST(track->Eta(),track->Phi(),track->Pt(),track->Charge());

      // mixed events
      if (poolReady) {
	for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) { // loop over mixed events
	  TObjArray *mixedTracks = pool->GetEvent(jMix);
	  FillHistogramsdPhidEtaMixed(iTr,track->Charge(),track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
				      mixedTracks);
	}
      }
      // end of mixed events
    }
  }


  // update event mixing pool
  pool->UpdatePool(ftrk);

  ////////////////////////////////////////////////////////////////////
  Int_t nhard = 0;
  Double_t phihard[nTracks], etahard[nTracks];
  Int_t indexhard[nTracks];
  for (Int_t iTr=0; iTr<nTracks; iTr++) {
    AliAODTrack *track = (AliAODTrack*) fAOD->GetTrack(iTr);
    if (!track) continue;
    if (track->TestFilterBit(768) && TMath::Abs(track->Eta()) < 0.8 && track->GetTPCNcls() >= 70 && track->Pt() >= fMinHardPt) {
      phihard[nhard] = track->Phi();
      etahard[nhard] = track->Eta();
      indexhard[nhard] = iTr;
      nhard++;
    }
  }  
  
  for (Int_t iTr=0; iTr<nTracks; iTr++) {
    AliAODTrack *track = (AliAODTrack*) fAOD->GetTrack(iTr);
    if (!track) continue;
    if (track->TestFilterBit(768) && TMath::Abs(track->Eta()) < 0.4 && track->GetTPCNcls() >= 70) {
      FillHistogramsV2(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
		       resA2, resC2, resT2,
		       resA2LowESE, resC2LowESE, resT2LowESE,
		       resA2HighESE, resC2HighESE, resT2HighESE, 0);
/*
      FillHistogramsV3(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
		       resA3, resC3, resT3,
		       resA3LowESE, resC3LowESE, resT3LowESE,
		       resA3HighESE, resC3HighESE, resT3HighESE, 0);

      Bool_t acceptV2 = kTRUE;
      Bool_t acceptV22 = kTRUE;
      Bool_t acceptV3 = kTRUE;
      Bool_t acceptV32 = kTRUE;
      for(Int_t ihard = 0; ihard < nhard; ++ihard) {
	if (indexhard[ihard] == iTr) continue;
	Double_t dphi = track->Phi()-phihard[ihard];
	if (dphi>TMath::Pi())  dphi -= TMath::TwoPi();
	if (dphi<-TMath::Pi()) dphi += TMath::TwoPi();
	if (TMath::Sqrt(dphi*dphi+(track->Eta()-etahard[ihard])*(track->Eta()-etahard[ihard]))<0.4) {
	  acceptV2 = kFALSE;
	  acceptV22 = kFALSE;
	  acceptV3 = kFALSE;
	  acceptV32 = kFALSE;
	}
	{
	  Double_t phihardv2 = phihard[ihard]+TMath::Pi()/2.;
	  if (phihardv2>TMath::TwoPi()) phihardv2 -= TMath::TwoPi();
	  Double_t dphiv2 = track->Phi()-phihardv2;
	  if (dphiv2>TMath::Pi())  dphiv2 -= TMath::TwoPi();
	  if (dphiv2<-TMath::Pi()) dphiv2 += TMath::TwoPi();
	  if (TMath::Sqrt(dphiv2*dphiv2+(track->Eta()-etahard[ihard])*(track->Eta()-etahard[ihard]))<0.4) {
	    acceptV22 = kFALSE;
	  }
	}
	{
	  Double_t phihardv3 = phihard[ihard]+TMath::Pi()/3.;
	  if (phihardv3>TMath::TwoPi()) phihardv3 -= TMath::TwoPi();
	  Double_t dphiv3 = track->Phi()-phihardv3;
	  if (dphiv3>TMath::Pi())  dphiv3 -= TMath::TwoPi();
	  if (dphiv3<-TMath::Pi()) dphiv3 += TMath::TwoPi();
	  if (TMath::Sqrt(dphiv3*dphiv3+(track->Eta()-etahard[ihard])*(track->Eta()-etahard[ihard]))<0.4) {
	    acceptV32 = kFALSE;
	  }
	}
      }
    
      if (acceptV2)
	FillHistogramsV2(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
			 resA2, resC2, resT2,
			 resA2LowESE, resC2LowESE, resT2LowESE,
			 resA2HighESE, resC2HighESE, resT2HighESE, 1);
      if (acceptV22)
	FillHistogramsV2(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
			 resA2, resC2, resT2,
			 resA2LowESE, resC2LowESE, resT2LowESE,
			 resA2HighESE, resC2HighESE, resT2HighESE, 2);
      if (acceptV3)
	FillHistogramsV3(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
			 resA3, resC3, resT3,
			 resA3LowESE, resC3LowESE, resT3LowESE,
			 resA3HighESE, resC3HighESE, resT3HighESE, 1);
      if (acceptV32)
	FillHistogramsV3(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
			 resA3, resC3, resT3,
			 resA3LowESE, resC3LowESE, resT3LowESE,
			 resA3HighESE, resC3HighESE, resT3HighESE, 2);
    */
    }
  }

  Double_t Qv0aQv0c2  = Qxa2Cor*Qxc2Cor  + Qya2Cor*Qyc2Cor;
  Double_t Qv0aQtrkl2 = Qxa2Cor*Qxtr2Cor + Qya2Cor*Qytr2Cor;
  Double_t Qv0cQtrkl2 = Qxc2Cor*Qxtr2Cor + Qyc2Cor*Qytr2Cor;
    
  fHistACv2->Fill(percentile, Qv0aQv0c2);
  fHistATv2->Fill(percentile, Qv0aQtrkl2);
  fHistCTv2->Fill(percentile, Qv0cQtrkl2);
  fHistACv2vsESE->Fill(percentile, percq2, Qv0aQv0c2);
  fHistATv2vsESE->Fill(percentile, percq2, Qv0aQtrkl2);
  fHistCTv2vsESE->Fill(percentile, percq2, Qv0cQtrkl2);

  fHistQxT2->Fill(percentile,Qxtr2Cor);
  fHistQyT2->Fill(percentile,Qytr2Cor);
  fHistQxT3->Fill(percentile,Qxtr3Cor);
  fHistQyT3->Fill(percentile,Qytr3Cor);

  if (percq2 < fLowESEcut) {
    fHistACv2LowESE->Fill(percentile, Qv0aQv0c2);
    fHistATv2LowESE->Fill(percentile, Qv0aQtrkl2);
    fHistCTv2LowESE->Fill(percentile, Qv0cQtrkl2);
  }
  if (percq2 > fHighESEcut) {
    fHistACv2HighESE->Fill(percentile, Qv0aQv0c2);
    fHistATv2HighESE->Fill(percentile, Qv0aQtrkl2);
    fHistCTv2HighESE->Fill(percentile, Qv0cQtrkl2);
  }
    
  Double_t Qv0aQv0c3  = Qxa3Cor*Qxc3Cor  + Qya3Cor*Qyc3Cor;
  Double_t Qv0aQtrkl3 = Qxa3Cor*Qxtr3Cor + Qya3Cor*Qytr3Cor;
  Double_t Qv0cQtrkl3 = Qxc3Cor*Qxtr3Cor + Qyc3Cor*Qytr3Cor;

  fHistACv3->Fill(percentile, Qv0aQv0c3);
  fHistATv3->Fill(percentile, Qv0aQtrkl3);
  fHistCTv3->Fill(percentile, Qv0cQtrkl3);
  fHistACv3vsESE->Fill(percentile, percq3, Qv0aQv0c3);
  fHistATv3vsESE->Fill(percentile, percq3, Qv0aQtrkl3);
  fHistCTv3vsESE->Fill(percentile, percq3, Qv0cQtrkl3);

  if (percq3 < fLowESEcut) {
    fHistACv3LowESE->Fill(percentile, Qv0aQv0c3);
    fHistATv3LowESE->Fill(percentile, Qv0aQtrkl3);
    fHistCTv3LowESE->Fill(percentile, Qv0cQtrkl3);
  }
  if (percq3 > fHighESEcut) {
    fHistACv3HighESE->Fill(percentile, Qv0aQv0c3);
    fHistATv3HighESE->Fill(percentile, Qv0aQtrkl3);
    fHistCTv3HighESE->Fill(percentile, Qv0cQtrkl3);
  }

  fHistv2ESEvsCent->Fill(percentile,percq2);
  fHistv3ESEvsCent->Fill(percentile,percq3);
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsV2(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
					     Double_t resA2, Double_t resC2, Double_t resT2,
					     Double_t resA2LowESE, Double_t resC2LowESE, Double_t resT2LowESE,
					     Double_t resA2HighESE, Double_t resC2HighESE, Double_t resT2HighESE,
					     Int_t index)
{
  {
    Double_t u2x = TMath::Cos(2.*phi);
    Double_t u2y = TMath::Sin(2.*phi);
    
    fHistSP2A[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2);
    fHistSP2C[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2);
    fHistSP2T[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2);
    fHistSP2AvsESE[centrality][zvtxBin][index]->Fill(pt,percq2,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2);
    fHistSP2CvsESE[centrality][zvtxBin][index]->Fill(pt,percq2,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2);
    fHistSP2TvsESE[centrality][zvtxBin][index]->Fill(pt,percq2,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2);
    if (percq2 < fLowESEcut) {
      fHistSP2ALowESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2LowESE);
      fHistSP2CLowESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2LowESE);
      fHistSP2TLowESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2LowESE);
    }
    if (percq2 > fHighESEcut) {
      fHistSP2AHighESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2HighESE);
      fHistSP2CHighESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2HighESE);
      fHistSP2THighESE[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2HighESE);
    }
  }
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsV3(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
					     Double_t resA3, Double_t resC3, Double_t resT3,
					     Double_t resA3LowESE, Double_t resC3LowESE, Double_t resT3LowESE,
					     Double_t resA3HighESE, Double_t resC3HighESE, Double_t resT3HighESE,
					     Int_t index)
{
  {
    Double_t u3x = TMath::Cos(3.*phi);
    Double_t u3y = TMath::Sin(3.*phi);

    fHistSP3A[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxa3Cor+u3y*Qya3Cor)/resA3);
    fHistSP3C[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxc3Cor+u3y*Qyc3Cor)/resC3);
    fHistSP3T[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxtr3Cor+u3y*Qytr3Cor)/resT3);    
    fHistSP3AvsESE[centrality][zvtxBin][index]->Fill(pt,percq3,(u3x*Qxa3Cor+u3y*Qya3Cor)/resA3);
    fHistSP3CvsESE[centrality][zvtxBin][index]->Fill(pt,percq3,(u3x*Qxc3Cor+u3y*Qyc3Cor)/resC3);
    fHistSP3TvsESE[centrality][zvtxBin][index]->Fill(pt,percq3,(u3x*Qxtr3Cor+u3y*Qytr3Cor)/resT3);    
    if (percq3 < fLowESEcut) {
      fHistSP3ALowESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxa3Cor+u3y*Qya3Cor)/resA3LowESE);
      fHistSP3CLowESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxc3Cor+u3y*Qyc3Cor)/resC3LowESE);
      fHistSP3TLowESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxtr3Cor+u3y*Qytr3Cor)/resT3LowESE);
    }
    if (percq3 > fHighESEcut) {
      fHistSP3AHighESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxa3Cor+u3y*Qya3Cor)/resA3HighESE);
      fHistSP3CHighESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxc3Cor+u3y*Qyc3Cor)/resC3HighESE);
      fHistSP3THighESE[centrality][zvtxBin][index]->Fill(pt,(u3x*Qxtr3Cor+u3y*Qytr3Cor)/resT3HighESE);
    }
  }
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsdPhidEta(Int_t trIndex, Short_t charge, Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
						   Double_t resA2, Double_t resC2, Double_t resT2,
						   Double_t resA3, Double_t resC3, Double_t resT3)
{
  Int_t ptBin = fPtTrigAxis->FindBin(pt);
  if (ptBin<1 || ptBin>fNbinsPtTrig) return;

  Double_t u2x = TMath::Cos(2.*phi);
  Double_t u2y = TMath::Sin(2.*phi);

  Double_t u3x = TMath::Cos(3.*phi);
  Double_t u3y = TMath::Sin(3.*phi);

  Double_t binscont_trig[3];
  binscont_trig[0] = pt; 
  binscont_trig[1] = fzvtx;
  binscont_trig[2] = percentile;
  fHistTrig->Fill(binscont_trig,0);

  Double_t binscont[6];


  const Int_t nTracks = fAOD->GetNumberOfTracks();
  for (Int_t iTr=0; iTr<nTracks; iTr++) {
    AliAODTrack *track = (AliAODTrack*) fAOD->GetTrack(iTr);
    if (!track) continue;
    if (iTr == trIndex) continue;
    if (track->TestFilterBit(768) && TMath::Abs(track->Eta()) < 0.8 && track->GetTPCNcls() >= 70) {
      // Get assoc pt bin
      Int_t assocPtBin = fPtAssocAxis->FindBin(track->Pt());
      if (assocPtBin<1 || assocPtBin>fNbinsAssocPt) continue;
      Double_t dphi = phi - track->Phi();
      if (dphi >  1.5*TMath::Pi()) dphi -= TMath::TwoPi();
      if (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
      Double_t deta = eta - track->Eta();
      fHistSP2AdPhidEta[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2);
      if (0) { // for output file size reasons
	fHistSP2CdPhidEta[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2);
	fHistSP2TdPhidEta[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2);
      }
      //fHistdPhidEta[centrality][ptBin-1][assocPtBin-1]->Fill(dphi,deta,fzvtx);
      binscont[0] = dphi;
      binscont[1] = deta;
      binscont[2] = fzvtx;
      binscont[3] = pt;
      binscont[4] = track->Pt();
      binscont[5] = percentile;
      fHistdPhidEtaPt->Fill(binscont,0);
      // fill same-sign track pair histos
      if (charge*track->Charge()>0) {
	fHistSP2AdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2);
	if (0) { // for output file size reasons
	  fHistSP2CdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2);
	  fHistSP2TdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxtr2Cor+u2y*Qytr2Cor)/resT2);
	}  
//	fHistdPhidEtaSS[centrality][ptBin-1][assocPtBin-1]->Fill(dphi,deta,fzvtx);
        fHistdPhidEtaPt_SS->Fill(binscont,0);
	fHistSP3AdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u3x*Qxa3Cor+u3y*Qya3Cor)/resA3);
      }
    }
  }
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsdPhidEtaMixed(Int_t trIndex, Short_t charge, Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin, TObjArray *mixedTracks)
{
  Double_t binscont[6];

  Int_t ptBin = fPtTrigAxis->FindBin(pt);
  if (ptBin<1 || ptBin>fNbinsPtTrig) return;
  for (Int_t iTr=0; iTr<mixedTracks->GetEntriesFast(); iTr++) {
    AliBasicParticleST *track = (AliBasicParticleST*) mixedTracks->UncheckedAt(iTr);
    // Get assoc pt bin
    Int_t assocPtBin = fPtAssocAxis->FindBin(track->Pt());
    if (assocPtBin<1 || assocPtBin>fNbinsAssocPt) continue;
    Double_t dphi = phi - track->Phi();
    if (dphi >  1.5*TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
    Double_t deta = eta - track->Eta();

//    fHistdPhidEtaMixed[centrality][ptBin-1][assocPtBin-1]->Fill(dphi,deta,fzvtx);
    binscont[0] = dphi;
    binscont[1] = deta;
    binscont[2] = fzvtx;
    binscont[3] = pt;
    binscont[4] = track->Pt();
    binscont[5] = percentile;
    fHistdPhidEtaPt_Mixed->Fill(binscont,0);  
  
    if (charge*track->Charge()>0) {
      //fHistdPhidEtaSSMixed[centrality][ptBin-1][assocPtBin-1]->Fill(dphi,deta,fzvtx);
      fHistdPhidEtaPt_Mixed_SS->Fill(binscont,0);  
    }
  }
}

//====================================================================================================================================================

Bool_t AliAnalysisTaskSEPbPbCorrelationsJetV2::IsTriggerFired() {
  
  Bool_t isSelected = kFALSE;
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7); 
  return isSelected;
}

//====================================================================================================================================================

Float_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetV0Multiplicity() {
  
  Float_t multiplicity=0;
  for (Int_t iChannel=0; iChannel<64; iChannel++) multiplicity += fAOD->GetVZEROData()->GetMultiplicity(iChannel);
  return multiplicity;

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetCentBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsCentrality) {
    AliInfo(Form("WARNING : only %d centrality bins (out of the %d proposed) will be considered",fNMaxBinsCentrality,nBins));
    nBins = fNMaxBinsCentrality;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one centrality bin must be considered");
    nBins = 1;
  }
  
  fNbinsCent = nBins;
  fCentAxis  = new TAxis(fNbinsCent, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetZvtxBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsZvtx) {
    AliInfo(Form("WARNING : only %d Zvtx bins (out of the %d proposed) will be considered",fNMaxBinsZvtx,nBins));
    nBins = fNMaxBinsZvtx;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one Zvtx bin must be considered");
    nBins = 1;
  }
  
  fNbinsZvtx = nBins;
  fZvtxAxis  = new TAxis(fNbinsZvtx, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetPtBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsPt) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsPt,nBins));
    nBins = fNMaxBinsPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fNbinsPt = nBins;
  fPtAxis  = new TAxis(fNbinsPt, limits);

}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetTrigPtBinning(Int_t nBins, Double_t *limits) {

/*
  if (nBins>fNMaxBinsPt) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsPt,nBins));
    nBins = fNMaxBinsPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
*/
  fNbinsPtTrig = nBins;
  fPtTrigAxis  = new TAxis(fNbinsPtTrig, limits);

}



//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetAssocPtBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsAssocPt) {
    AliInfo(Form("WARNING : only %d assoc pt bins (out of the %d proposed) will be considered",fNMaxBinsAssocPt,nBins));
    nBins = fNMaxBinsAssocPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one assoc pt bin must be considered");
    nBins = 1;
  }
  
  fNbinsAssocPt = nBins;
  fPtAssocAxis  = new TAxis(fNbinsAssocPt, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetEtaBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsEta) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsEta,nBins));
    nBins = fNMaxBinsEta;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fEtaAxis  = new TAxis(nBins, limits);

}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetCentBin() {

  Double_t percentile;
    AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
    percentile = multSelection->GetMultiplicityPercentile(fCentMethod.Data());

  Int_t bin = fCentAxis->FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetCentBin(Double_t percentile) {

  Int_t bin = fCentAxis->FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetZvtxBin(Double_t zvtx) {

  Int_t bin = fZvtxAxis->FindBin(zvtx) - 1;
  if (bin >= fNbinsZvtx) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Double_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetITSMultiplicity() {

  Double_t multiplicity = ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(0)+
    ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(1);

  return multiplicity;

}

Double_t AliAnalysisTaskSEPbPbCorrelationsJetV2::CalcCorrectedPhi(Double_t phi, Double_t dPhi) const{
  phi += 39./34.*dPhi;
  if (phi < 0.)  phi+=2*TMath::Pi();
  if (phi > 2.*TMath::Pi()) phi-=2*TMath::Pi();
  return phi;
}


//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::Terminate(Option_t *) {

  if (fPoolMgr)    delete fPoolMgr;

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }  

  fOutputList1 = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputList1) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//__________________________________________________________
Short_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetVertexZ(Double_t vtxZ) const
{
  // this method gives the z vtx bin used to access flow vector calibrations
    Short_t zvtx = -10;
    
    if (vtxZ >= -10. && vtxZ < -8.)
        zvtx = 0;
    if (vtxZ >= -8. && vtxZ < -6.)
        zvtx = 1;
    if (vtxZ >= -6. && vtxZ < -4.)
        zvtx = 2;
    if (vtxZ >= -4. && vtxZ < -3.)
        zvtx = 3;
    if (vtxZ >= -3. && vtxZ < -2.)
        zvtx = 4;
    if (vtxZ >= -2. && vtxZ < -1.)
        zvtx = 5;
    if (vtxZ >= -1. && vtxZ < 0)
        zvtx = 6;
    if (vtxZ >= 0 && vtxZ < 1.)
        zvtx = 7;
    if (vtxZ >= 1. && vtxZ < 2.)
        zvtx = 8;
    if (vtxZ >= 2. && vtxZ < 3.)
        zvtx = 9;
    if (vtxZ >= 3. && vtxZ < 4.)
        zvtx = 10;
    if (vtxZ >= 4. && vtxZ < 6.)
        zvtx = 11;
    if (vtxZ >= 6. && vtxZ < 8.)
        zvtx = 12;
    if (vtxZ >= 8. && vtxZ <= 10.)
        zvtx = 13;
    
    return zvtx;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEPbPbCorrelationsJetV2::ComputeQ(AliAODEvent* aod, Double_t Zvtx)
{  

    AliAODTracklets* aodTrkl = (AliAODTracklets*)aod->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    
    AliAODVZERO* aodV0 = (AliAODVZERO*)aod->GetVZEROData();

    Int_t iCentV0 = Int_t(fv0mpercentile);
    if (iCentV0 >= 90)
      return kFALSE;
    
    Int_t iCentSPD = Int_t(fcl1percentile);
    if (iCentSPD >= 90)
      return kFALSE;
    
    
    Short_t zvt = GetVertexZ(Zvtx);
    if (zvt < 0)
      return kFALSE;
    
    //Tracklets
    Qx2trkcut = Qy2trkcut = 0;
    Qx3trkcut = Qy3trkcut = 0;
    multtrkcut = 0;
    for (Int_t it = 0; it < nITSTrkls; it++){

      // We do not use these cuts anymore
      //        if (TMath::Abs(aodTrkl->GetEta(it)) > 0.5 || TMath::Abs(aodTrkl->GetDeltaPhi(it)) >= 0.005)
      //            continue;
        
        Double_t corPhi = CalcCorrectedPhi(aodTrkl->GetPhi(it), aodTrkl->GetDeltaPhi(it));
        
        Qx2trkcut += TMath::Cos(2.*corPhi);
        Qy2trkcut += TMath::Sin(2.*corPhi);
        
        Qx3trkcut += TMath::Cos(3.*corPhi);
        Qy3trkcut += TMath::Sin(3.*corPhi);
        
        multtrkcut++;
    }
    
    //V0 info
    Qxa2 = Qya2 = 0;
    Qxc2 = Qyc2 = 0;
    
    Qxa3 = Qya3 = 0;
    Qxc3 = Qyc3 = 0;
    
    sumMa = sumMc = 0;
    
    for (Int_t iV0 = 0; iV0 < 64; iV0++) {
        
        Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
        Float_t multv0 = aodV0->GetMultiplicity(iV0);
        
        if (iV0 < 32){
            
            Double_t multCorC = -10;
            
            if (iV0 < 8)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
            else if (iV0 >= 8 && iV0 < 16)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
            else if (iV0 >= 16 && iV0 < 24)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
            else if (iV0 >= 24 && iV0 < 32)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);
            
            if (multCorC < 0){
              //  cout<<"Problem with multiplicity in V0C"<<endl;
                continue;
            }
            
            Qxc2 += TMath::Cos(2.*phiV0) * multCorC;
            Qyc2 += TMath::Sin(2.*phiV0) * multCorC;
            
            Qxc3 += TMath::Cos(3.*phiV0) * multCorC;
            Qyc3 += TMath::Sin(3.*phiV0) * multCorC;
 
            sumMc = sumMc + multCorC;
            
        } else {
            
            Double_t multCorA = -10;
            
            if (iV0 >= 32 && iV0 < 40)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
            else if (iV0 >= 40 && iV0 < 48)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
            else if (iV0 >= 48 && iV0 < 56)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
            else if (iV0 >= 56 && iV0 < 64)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);
            
            if (multCorA < 0){
            //    cout<<"Problem with multiplicity in V0A"<<endl;
                continue;
            }
            
            Qxa2 += TMath::Cos(2.*phiV0) * multCorA;
            Qya2 += TMath::Sin(2.*phiV0) * multCorA;
            
            Qxa3 += TMath::Cos(3.*phiV0) * multCorA;
            Qya3 += TMath::Sin(3.*phiV0) * multCorA;

            sumMa = sumMa + multCorA;
            
        }
    }

    if (sumMa <= 0 || sumMc <= 0 || multtrkcut <= 0)
      return kFALSE;
    
    
    
    //might need to check the sigma values != 0
    if (fQy2sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQy3sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQy2sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQy3sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQy2sTrk[zvt]->GetBinContent(iCentV0+1)<1e-8 ||
	fQy3sTrk[zvt]->GetBinContent(iCentV0+1)<1e-8)
      return kFALSE;
    if (fQx2sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQx3sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQx2sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQx3sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQx2sTrk[zvt]->GetBinContent(iCentV0+1)<1e-8 ||
	fQx3sTrk[zvt]->GetBinContent(iCentV0+1)<1e-8)
      return kFALSE;
    
    Qya2Cor = (Qya2 - fQy2mV0A[zvt]->GetBinContent(iCentSPD+1))/fQy2sV0A[zvt]->GetBinContent(iCentSPD+1);
    Qxa2Cor = (Qxa2 - fQx2mV0A[zvt]->GetBinContent(iCentSPD+1))/fQx2sV0A[zvt]->GetBinContent(iCentSPD+1);
    
    Qya3Cor = (Qya3 - fQy3mV0A[zvt]->GetBinContent(iCentSPD+1))/fQy3sV0A[zvt]->GetBinContent(iCentSPD+1);
    Qxa3Cor = (Qxa3 - fQx3mV0A[zvt]->GetBinContent(iCentSPD+1))/fQx3sV0A[zvt]->GetBinContent(iCentSPD+1);
    
    
    Qyc2Cor = (Qyc2 - fQy2mV0C[zvt]->GetBinContent(iCentSPD+1))/fQy2sV0C[zvt]->GetBinContent(iCentSPD+1);
    Qxc2Cor = (Qxc2 - fQx2mV0C[zvt]->GetBinContent(iCentSPD+1))/fQx2sV0C[zvt]->GetBinContent(iCentSPD+1);
    
    Qyc3Cor = (Qyc3 - fQy3mV0C[zvt]->GetBinContent(iCentSPD+1))/fQy3sV0C[zvt]->GetBinContent(iCentSPD+1);
    Qxc3Cor = (Qxc3 - fQx3mV0C[zvt]->GetBinContent(iCentSPD+1))/fQx3sV0C[zvt]->GetBinContent(iCentSPD+1);
    
    
    
    Qytr2Cor = (Qy2trkcut - fQy2mTrk[zvt]->GetBinContent(iCentV0+1))/fQy2sTrk[zvt]->GetBinContent(iCentV0+1);
    Qxtr2Cor = (Qx2trkcut - fQx2mTrk[zvt]->GetBinContent(iCentV0+1))/fQx2sTrk[zvt]->GetBinContent(iCentV0+1);
    
    Qytr3Cor = (Qy3trkcut - fQy3mTrk[zvt]->GetBinContent(iCentV0+1))/fQy3sTrk[zvt]->GetBinContent(iCentV0+1);
    Qxtr3Cor = (Qx3trkcut - fQx3mTrk[zvt]->GetBinContent(iCentV0+1))/fQx3sTrk[zvt]->GetBinContent(iCentV0+1);

    if (fEP) { // in case we want to do EP analysis
      Double_t norm = TMath::Sqrt(Qxa2Cor*Qxa2Cor+Qya2Cor*Qya2Cor);
      if (norm > 0) {
	Qxa2Cor /= norm;
	Qya2Cor /= norm;
      }
      else {
	Qxa2Cor = Qya2Cor = 0;
      }
      norm = TMath::Sqrt(Qxa3Cor*Qxa3Cor+Qya3Cor*Qya3Cor);
      if (norm > 0) {
	Qxa3Cor /= norm;
	Qya3Cor /= norm;
      }
      else {
	Qxa3Cor = Qya3Cor = 0;
      }
      norm = TMath::Sqrt(Qxc2Cor*Qxc2Cor+Qyc2Cor*Qyc2Cor);
      if (norm > 0) {
	Qxc2Cor /= norm;
	Qyc2Cor /= norm;
      }
      else {
	Qxc2Cor = Qyc2Cor = 0;
      }
      norm = TMath::Sqrt(Qxc3Cor*Qxc3Cor+Qyc3Cor*Qyc3Cor);
      if (norm > 0) {
	Qxc3Cor /= norm;
	Qyc3Cor /= norm;
      }
      else {
	Qxc3Cor = Qyc3Cor = 0;
      }
      norm = TMath::Sqrt(Qxtr2Cor*Qxtr2Cor+Qytr2Cor*Qytr2Cor);
      if (norm > 0) {
	Qxtr2Cor /= norm;
	Qytr2Cor /= norm;
      }
      else {
	Qxtr2Cor = Qytr2Cor = 0;
      }
      norm = TMath::Sqrt(Qxtr3Cor*Qxtr3Cor+Qytr3Cor*Qytr3Cor);
      if (norm > 0) {
	Qxtr3Cor /= norm;
	Qytr3Cor /= norm;
      }
      else {
	Qxtr3Cor = Qytr3Cor = 0;
      }
    }
    
    if (fESEdet == 2) {
      // tracklets
      Double_t Qytr2CorESE = Qy2trkcut - fQy2mTrk[zvt]->GetBinContent(iCentV0+1);
      Double_t Qxtr2CorESE = Qx2trkcut - fQx2mTrk[zvt]->GetBinContent(iCentV0+1);
      
      Double_t q2Tr = TMath::Sqrt((Qxtr2CorESE*Qxtr2CorESE + Qytr2CorESE*Qytr2CorESE)/multtrkcut);
      Double_t percq2Tr = 100.*fSplQ2trk[iCentV0]->Eval(q2Tr);
      
      Double_t Qytr3CorESE = Qy3trkcut - fQy3mTrk[zvt]->GetBinContent(iCentV0+1);
      Double_t Qxtr3CorESE = Qx3trkcut - fQx3mTrk[zvt]->GetBinContent(iCentV0+1);
      
      Double_t q3Tr = TMath::Sqrt((Qxtr3CorESE*Qxtr3CorESE + Qytr3CorESE*Qytr3CorESE)/multtrkcut);
      Double_t percq3Tr = 100.*fSplQ3trk[iCentV0]->Eval(q3Tr);

      percq2 = percq2Tr;
      percq3 = percq3Tr;
    }
    else if (fESEdet == 0) {
      // V0A
      Double_t Qya2CorESE = Qya2 - fQy2mV0A[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxa2CorESE = Qxa2 - fQx2mV0A[zvt]->GetBinContent(iCentSPD+1);

      Double_t q2V0A = TMath::Sqrt((Qxa2CorESE*Qxa2CorESE + Qya2CorESE*Qya2CorESE)/sumMa);
      Double_t percq2V0A = 100.*fSplQ2V0A[iCentV0]->Eval(q2V0A);

      Double_t Qya3CorESE = Qya3 - fQy3mV0A[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxa3CorESE = Qxa3 - fQx3mV0A[zvt]->GetBinContent(iCentSPD+1);

      Double_t q3V0A = TMath::Sqrt((Qxa3CorESE*Qxa3CorESE + Qya3CorESE*Qya3CorESE)/sumMa);
      Double_t percq3V0A = 100.*fSplQ3V0A[iCentV0]->Eval(q3V0A);
      percq2 = percq2V0A;
      percq3 = percq3V0A;
    }
    else if (fESEdet == 1) {
      // V0C
      Double_t Qyc2CorESE = Qyc2 - fQy2mV0C[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxc2CorESE = Qxc2 - fQx2mV0C[zvt]->GetBinContent(iCentSPD+1);
    
      Double_t q2V0C = TMath::Sqrt((Qxc2CorESE*Qxc2CorESE + Qyc2CorESE*Qyc2CorESE)/sumMc);
      Double_t percq2V0C = 100.*fSplQ2V0C[iCentV0]->Eval(q2V0C);
    
      Double_t Qyc3CorESE = Qyc3 - fQy3mV0C[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxc3CorESE = Qxc3 - fQx3mV0C[zvt]->GetBinContent(iCentSPD+1);
    
      Double_t q3V0C = TMath::Sqrt((Qxc3CorESE*Qxc3CorESE + Qyc3CorESE*Qyc3CorESE)/sumMc);
      Double_t percq3V0C = 100.*fSplQ3V0C[iCentV0]->Eval(q3V0C);

      percq2 = percq2V0C;
      percq3 = percq3V0C;
    }
    else if (fESEdet == 3) {
      // V0A + V0C
      Double_t Qya2CorESE = Qya2 - fQy2mV0A[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxa2CorESE = Qxa2 - fQx2mV0A[zvt]->GetBinContent(iCentSPD+1);

      Double_t Qyc2CorESE = Qyc2 - fQy2mV0C[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxc2CorESE = Qxc2 - fQx2mV0C[zvt]->GetBinContent(iCentSPD+1);

      Double_t q2V0AC = TMath::Sqrt((Qxa2CorESE*Qxa2CorESE + Qya2CorESE*Qya2CorESE)/sumMa) + TMath::Sqrt((Qxc2CorESE*Qxc2CorESE + Qyc2CorESE*Qyc2CorESE)/sumMc);
      Double_t percq2V0AC = 100.*fSplQ2V0AC[iCentV0]->Eval(q2V0AC);

      Double_t Qya3CorESE = Qya3 - fQy3mV0A[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxa3CorESE = Qxa3 - fQx3mV0A[zvt]->GetBinContent(iCentSPD+1);

      Double_t Qyc3CorESE = Qyc3 - fQy3mV0C[zvt]->GetBinContent(iCentSPD+1);
      Double_t Qxc3CorESE = Qxc3 - fQx3mV0C[zvt]->GetBinContent(iCentSPD+1);
     
      Double_t q3V0AC = TMath::Sqrt((Qxa3CorESE*Qxa3CorESE + Qya3CorESE*Qya3CorESE)/sumMa) + TMath::Sqrt((Qxc3CorESE*Qxc3CorESE + Qyc3CorESE*Qyc3CorESE)/sumMc);
      Double_t percq3V0AC = 100.*fSplQ3V0AC[iCentV0]->Eval(q3V0AC);

      percq2 = percq2V0AC;
      percq3 = percq3V0AC;
    }
    else {
      printf("Wrong choice of ESE detector\n");
      return kFALSE;
    }
    
    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPbPbCorrelationsJetV2::OpenInfoCalbration(Int_t run)
{
  TFile *foadb = (TFile*)gROOT->GetListOfFiles()->FindObject("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibV0TrklNoEtaCutRun2.root");
  if (!foadb) {
    foadb = TFile::Open("./calibV0TrklNoEtaCutRun2.root");
  }
    
  if (!foadb) {
    if (!gGrid) {
        TGrid::Connect("alien://");
    }
    
    foadb = TFile::Open("alien:///alice/cern.ch/user/s/sitang/EP_Cali_run2/calibV0TrklNoEtaCutRun2.root");
  }
    
    if(!foadb){
        printf("OADB V0-Trkl calibration file cannot be opened\n");
        return;
    }
    
    if (!cont) cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
    if(!cont){
        printf("OADB object hMultV0BefCorPfpx is not available in the file\n");
        return;
    }
    if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
    }
    fMultV0 = ((TH1D*) cont->GetObject(run));
    
    
    
    for (Int_t k = 0; k < 14; k++){
        
        if (!contQx2am[k]) contQx2am[k] = (AliOADBContainer*) foadb->Get(Form("fqxa2m_%d", k));
        if(!contQx2am[k]){
            printf("OADB object fqxa2m is not available in the file\n");
            return;
        }
        if(!(contQx2am[k]->GetObject(run))){
            printf("OADB object fqxa2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0A[k]= ((TH1D*) contQx2am[k]->GetObject(run));
        
        
        if (!contQy2am[k]) contQy2am[k] = (AliOADBContainer*) foadb->Get(Form("fqya2m_%d", k));
        if(!contQy2am[k]){
            printf("OADB object fqya2m is not available in the file\n");
            return;
        }
        if(!(contQy2am[k]->GetObject(run))){
            printf("OADB object fqya2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0A[k]= ((TH1D*) contQy2am[k]->GetObject(run));
        
        
        
        if (!contQx2as[k]) contQx2as[k] = (AliOADBContainer*) foadb->Get(Form("fqxa2s_%d", k));
        if(!contQx2as[k]){
            printf("OADB object fqxa2s is not available in the file\n");
            return;
        }
        if(!(contQx2as[k]->GetObject(run))){
            printf("OADB object fqxa2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0A[k]= ((TH1D*) contQx2as[k]->GetObject(run));
        
        
        if (!contQy2as[k]) contQy2as[k] = (AliOADBContainer*) foadb->Get(Form("fqya2s_%d", k));
        if(!contQy2as[k]){
            printf("OADB object fqya2s is not available in the file\n");
            return;
        }
        if(!(contQy2as[k]->GetObject(run))){
            printf("OADB object fqya2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0A[k]= ((TH1D*) contQy2as[k]->GetObject(run));
        
        
        
        
        if (!contQx2cm[k]) contQx2cm[k] = (AliOADBContainer*) foadb->Get(Form("fqxc2m_%d", k));
        if(!contQx2cm[k]){
            printf("OADB object fqxc2m is not available in the file\n");
            return;
        }
        if(!(contQx2cm[k]->GetObject(run))){
            printf("OADB object fqxc2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0C[k]= ((TH1D*) contQx2cm[k]->GetObject(run));
        
        
        if (!contQy2cm[k]) contQy2cm[k] = (AliOADBContainer*) foadb->Get(Form("fqyc2m_%d", k));
        if(!contQy2cm[k]){
            printf("OADB object fqyc2m is not available in the file\n");
            return;
        }
        if(!(contQy2cm[k]->GetObject(run))){
            printf("OADB object fqyc2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0C[k]= ((TH1D*) contQy2cm[k]->GetObject(run));
        
        
        if (!contQx2cs[k]) contQx2cs[k] = (AliOADBContainer*) foadb->Get(Form("fqxc2s_%d", k));
        if(!contQx2cs[k]){
            printf("OADB object fqxc2s is not available in the file\n");
            return;
        }
        if(!(contQx2cs[k]->GetObject(run))){
            printf("OADB object fqxc2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0C[k]= ((TH1D*) contQx2cs[k]->GetObject(run));
        
        
        if (!contQy2cs[k]) contQy2cs[k] = (AliOADBContainer*) foadb->Get(Form("fqyc2s_%d", k));
        if(!contQy2cs[k]){
            printf("OADB object fqyc2s is not available in the file\n");
            return;
        }
        if(!(contQy2cs[k]->GetObject(run))){
            printf("OADB object fqyc2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0C[k]= ((TH1D*) contQy2cs[k]->GetObject(run));
        
        
        
        
        if (!contQx2trm[k]) contQx2trm[k] = (AliOADBContainer*) foadb->Get(Form("fqxtr2m_%d", k));
        if(!contQx2trm[k]){
            printf("OADB object fqxtr2m is not available in the file\n");
            return;
        }
        if(!(contQx2trm[k]->GetObject(run))){
            printf("OADB object fqxtr2m is not available for run %i\n", run);
            return;
        }
        fQx2mTrk[k]= ((TH1D*) contQx2trm[k]->GetObject(run));
        
        
        if (!contQy2trm[k]) contQy2trm[k] = (AliOADBContainer*) foadb->Get(Form("fqytr2m_%d", k));
        if(!contQy2trm[k]){
            printf("OADB object fqytr2m is not available in the file\n");
            return;
        }
        if(!(contQy2trm[k]->GetObject(run))){
            printf("OADB object fqytr2m is not available for run %i\n", run);
            return;
        }
        fQy2mTrk[k]= ((TH1D*) contQy2trm[k]->GetObject(run));
        
        
        if (!contQx2trs[k]) contQx2trs[k] = (AliOADBContainer*) foadb->Get(Form("fqxtr2s_%d", k));
        if(!contQx2trs[k]){
            printf("OADB object fqxtr2s is not available in the file\n");
            return;
        }
        if(!(contQx2trs[k]->GetObject(run))){
            printf("OADB object fqxtr2s is not available for run %i\n", run);
            return;
        }
        fQx2sTrk[k]= ((TH1D*) contQx2trs[k]->GetObject(run));
        
        
        if (!contQy2trs[k]) contQy2trs[k] = (AliOADBContainer*) foadb->Get(Form("fqytr2s_%d", k));
        if(!contQy2trs[k]){
            printf("OADB object fqytr2s is not available in the file\n");
            return;
        }
        if(!(contQy2trs[k]->GetObject(run))){
            printf("OADB object fqytr2s is not available for run %i\n", run);
            return;
        }
        fQy2sTrk[k]= ((TH1D*) contQy2trs[k]->GetObject(run));
        
        
        
        
        if (!contQx3am[k]) contQx3am[k] = (AliOADBContainer*) foadb->Get(Form("fqxa3m_%d", k));
        if(!contQx3am[k]){
            printf("OADB object fqxa3m is not available in the file\n");
            return;
        }
        if(!(contQx3am[k]->GetObject(run))){
            printf("OADB object fqxa3m is not available for run %i\n", run);
            return;
        }
        fQx3mV0A[k]= ((TH1D*) contQx3am[k]->GetObject(run));
        
        
        if (!contQy3am[k]) contQy3am[k] = (AliOADBContainer*) foadb->Get(Form("fqya3m_%d", k));
        if(!contQy3am[k]){
            printf("OADB object fqya3m is not available in the file\n");
            return;
        }
        if(!(contQy3am[k]->GetObject(run))){
            printf("OADB object fqya3m is not available for run %i\n", run);
            return;
        }
        fQy3mV0A[k]= ((TH1D*) contQy3am[k]->GetObject(run));
        
        
        if (!contQx3as[k]) contQx3as[k] = (AliOADBContainer*) foadb->Get(Form("fqxa3s_%d", k));
        if(!contQx3as[k]){
            printf("OADB object fqxa3s is not available in the file\n");
            return;
        }
        if(!(contQx3as[k]->GetObject(run))){
            printf("OADB object fqxa3s is not available for run %i\n", run);
            return;
        }
        fQx3sV0A[k]= ((TH1D*) contQx3as[k]->GetObject(run));
        
        
        if (!contQy3as[k]) contQy3as[k] = (AliOADBContainer*) foadb->Get(Form("fqya3s_%d", k));
        if(!contQy3as[k]){
            printf("OADB object fqya3s is not available in the file\n");
            return;
        }
        if(!(contQy3as[k]->GetObject(run))){
            printf("OADB object fqya3s is not available for run %i\n", run);
            return;
        }
        fQy3sV0A[k]= ((TH1D*) contQy3as[k]->GetObject(run));
        
        
        
        
        if (!contQx3cm[k]) contQx3cm[k] = (AliOADBContainer*) foadb->Get(Form("fqxc3m_%d", k));
        if(!contQx3cm[k]){
            printf("OADB object fqxc3m is not available in the file\n");
            return;
        }
        if(!(contQx3cm[k]->GetObject(run))){
            printf("OADB object fqxc3m is not available for run %i\n", run);
            return;
        }
        fQx3mV0C[k]= ((TH1D*) contQx3cm[k]->GetObject(run));
        
        
        if (!contQy3cm[k]) contQy3cm[k] = (AliOADBContainer*) foadb->Get(Form("fqyc3m_%d", k));
        if(!contQy3cm[k]){
            printf("OADB object fqyc3m is not available in the file\n");
            return;
        }
        if(!(contQy3cm[k]->GetObject(run))){
            printf("OADB object fqyc3m is not available for run %i\n", run);
            return;
        }
        fQy3mV0C[k]= ((TH1D*) contQy3cm[k]->GetObject(run));
        
        
        if (!contQx3cs[k]) contQx3cs[k] = (AliOADBContainer*) foadb->Get(Form("fqxc3s_%d", k));
        if(!contQx3cs[k]){
            printf("OADB object fqxc3s is not available in the file\n");
            return;
        }
        if(!(contQx3cs[k]->GetObject(run))){
            printf("OADB object fqxc3s is not available for run %i\n", run);
            return;
        }
        fQx3sV0C[k]= ((TH1D*) contQx3cs[k]->GetObject(run));
        
        
        if (!contQy3cs[k]) contQy3cs[k] = (AliOADBContainer*) foadb->Get(Form("fqyc3s_%d", k));
        if(!contQy3cs[k]){
            printf("OADB object fqyc3s is not available in the file\n");
            return;
        }
        if(!(contQy3cs[k]->GetObject(run))){
            printf("OADB object fqyc3s is not available for run %i\n", run);
            return;
        }
        fQy3sV0C[k]= ((TH1D*) contQy3cs[k]->GetObject(run));

        
        
        
        if (!contQx3trm[k]) contQx3trm[k] = (AliOADBContainer*) foadb->Get(Form("fqxtr3m_%d", k));
        if(!contQx3trm[k]){
            printf("OADB object fqxtr3m is not available in the file\n");
            return;
        }
        if(!(contQx3trm[k]->GetObject(run))){
            printf("OADB object fqxtr3m is not available for run %i\n", run);
            return;
        }
        fQx3mTrk[k]= ((TH1D*) contQx3trm[k]->GetObject(run));
        
        
        if (!contQy3trm[k]) contQy3trm[k] = (AliOADBContainer*) foadb->Get(Form("fqytr3m_%d", k));
        if(!contQy3trm[k]){
            printf("OADB object fqytr3m is not available in the file\n");
            return;
        }
        if(!(contQy3trm[k]->GetObject(run))){
            printf("OADB object fqytr3m is not available for run %i\n", run);
            return;
        }
        fQy3mTrk[k]= ((TH1D*) contQy3trm[k]->GetObject(run));
        
        
        if (!contQx3trs[k]) contQx3trs[k] = (AliOADBContainer*) foadb->Get(Form("fqxtr3s_%d", k));
        if(!contQx3trs[k]){
            printf("OADB object fqxtr3s is not available in the file\n");
            return;
        }
        if(!(contQx3trs[k]->GetObject(run))){
            printf("OADB object fqxtr3s is not available for run %i\n", run);
            return;
        }
        fQx3sTrk[k]= ((TH1D*) contQx3trs[k]->GetObject(run));
        
        
        if (!contQy3trs[k]) contQy3trs[k] = (AliOADBContainer*) foadb->Get(Form("fqytr3s_%d", k));
        if(!contQy3trs[k]){
            printf("OADB object fqytr3s is not available in the file\n");
            return;
        }
        if(!(contQy3trs[k]->GetObject(run))){
            printf("OADB object fqytr3s is not available for run %i\n", run);
            return;
        }
        fQy3sTrk[k]= ((TH1D*) contQy3trs[k]->GetObject(run));
        
    }

    
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::CalcResolutions(Double_t percentile,
					      Double_t &resA2, Double_t &resC2, Double_t &resT2,
					      Double_t &resA2LowESE, Double_t &resC2LowESE, Double_t &resT2LowESE,
					      Double_t &resA2HighESE, Double_t &resC2HighESE, Double_t &resT2HighESE,
					      Double_t &resA3, Double_t &resC3, Double_t &resT3,
					      Double_t &resA3LowESE, Double_t &resC3LowESE, Double_t &resT3LowESE,
					      Double_t &resA3HighESE, Double_t &resC3HighESE, Double_t &resT3HighESE)
{
  if (fAverageRes) {
    // average resolution
    // v2
    {
      Double_t ac2 = fResACv2->Eval(percentile);
      Double_t at2 = fResATv2->Eval(percentile);
      Double_t ct2 = fResCTv2->Eval(percentile);
      if (ac2 <= 0 || at2 <= 0 || ct2 <= 0) {
	resA2 = resC2 = resT2 = 1e6;
      }
      else {
	resA2 = TMath::Sqrt(ac2*at2/ct2);
	resC2 = TMath::Sqrt(ac2*ct2/at2);
	resT2 = TMath::Sqrt(at2*ct2/ac2);
      }
    }
    {
      Double_t ac2l = fResACv2LowESE->Eval(percentile);
      Double_t at2l = fResATv2LowESE->Eval(percentile);
      Double_t ct2l = fResCTv2LowESE->Eval(percentile);
      if (ac2l <= 0 || at2l <= 0 || ct2l <= 0) {
	resA2LowESE = resC2LowESE = resT2LowESE = 1e6;
      }
      else {
	resA2LowESE = TMath::Sqrt(ac2l*at2l/ct2l);
	resC2LowESE = TMath::Sqrt(ac2l*ct2l/at2l);
	resT2LowESE = TMath::Sqrt(at2l*ct2l/ac2l);
      }
    }
    {
      Double_t ac2h = fResACv2HighESE->Eval(percentile);
      Double_t at2h = fResATv2HighESE->Eval(percentile);
      Double_t ct2h = fResCTv2HighESE->Eval(percentile);
      if (ac2h <= 0 || at2h <= 0 || ct2h <= 0) {
	resA2HighESE = resC2HighESE = resT2HighESE = 1e6;
      }
      else {
	resA2HighESE = TMath::Sqrt(ac2h*at2h/ct2h);
	resC2HighESE = TMath::Sqrt(ac2h*ct2h/at2h);
	resT2HighESE = TMath::Sqrt(at2h*ct2h/ac2h);
      }
    }
    
    // v3  
    {
      Double_t ac3 = fResACv3->Eval(percentile);
      Double_t at3 = fResATv3->Eval(percentile);
      Double_t ct3 = fResCTv3->Eval(percentile);
      if (ac3 <= 0 || at3 <= 0 || ct3 <= 0) {
	resA3 = resC3 = resT3 = 1e6;
      }
      else {
	resA3 = TMath::Sqrt(ac3*at3/ct3);
	resC3 = TMath::Sqrt(ac3*ct3/at3);
	resT3 = TMath::Sqrt(at3*ct3/ac3);
      }
    }
    {
      Double_t ac3l = fResACv3LowESE->Eval(percentile);
      Double_t at3l = fResATv3LowESE->Eval(percentile);
      Double_t ct3l = fResCTv3LowESE->Eval(percentile);
      if (ac3l <= 0 || at3l <= 0 || ct3l <= 0) {
	resA3LowESE = resC3LowESE = resT3LowESE = 1e6;
      }
      else {
	resA3LowESE = TMath::Sqrt(ac3l*at3l/ct3l);
	resC3LowESE = TMath::Sqrt(ac3l*ct3l/at3l);
	resT3LowESE = TMath::Sqrt(at3l*ct3l/ac3l);
      }
    }
    {
      Double_t ac3h = fResACv3HighESE->Eval(percentile);
      Double_t at3h = fResATv3HighESE->Eval(percentile);
      Double_t ct3h = fResCTv3HighESE->Eval(percentile);
      if (ac3h <= 0 || at3h <= 0 || ct3h <= 0) {
	resA3HighESE = resC3HighESE = resT3HighESE = 1e6;
      }
      else {
	resA3HighESE = TMath::Sqrt(ac3h*at3h/ct3h);
	resC3HighESE = TMath::Sqrt(ac3h*ct3h/at3h);
	resT3HighESE = TMath::Sqrt(at3h*ct3h/ac3h);
      }
    }
  }
  else {
    // Resolution vs ESE
    // v2
    {
      Int_t index = Int_t(percq2/20.);
      if (index < 0) index = 0;
      if (index >= 5) index = 4;
      Double_t ac2 = fResACv2vsESE[index]->Eval(percentile);
      Double_t at2 = fResATv2vsESE[index]->Eval(percentile);
      Double_t ct2 = fResCTv2vsESE[index]->Eval(percentile);
      if (ac2 <= 0 || at2 <= 0 || ct2 <= 0) {
	resA2 = resC2 = resT2 = 1e6;
	resA2LowESE = resC2LowESE = resT2LowESE = 1e6;
	resA2HighESE = resC2HighESE = resT2HighESE = 1e6;
      }
      else {
	resA2 = resA2LowESE = resA2HighESE = TMath::Sqrt(ac2*at2/ct2);
	resC2 = resC2LowESE = resC2HighESE = TMath::Sqrt(ac2*ct2/at2);
	resT2 = resT2LowESE = resT2HighESE = TMath::Sqrt(at2*ct2/ac2);
      }
    }
    
    // v3  
    {
      Int_t index = Int_t(percq3/20.);
      if (index < 0) index = 0;
      if (index >= 5) index = 4;
      Double_t ac3 = fResACv3vsESE[index]->Eval(percentile);
      Double_t at3 = fResATv3vsESE[index]->Eval(percentile);
      Double_t ct3 = fResCTv3vsESE[index]->Eval(percentile);
      if (ac3 <= 0 || at3 <= 0 || ct3 <= 0) {
	resA3 = resC3 = resT3 = 1e6;
	resA3LowESE = resC3LowESE = resT3LowESE = 1e6;
	resA3HighESE = resC3HighESE = resT3HighESE = 1e6;
      }
      else {
	resA3 = resA3LowESE = resA3HighESE = TMath::Sqrt(ac3*at3/ct3);
	resC3 = resC3LowESE = resC3HighESE = TMath::Sqrt(ac3*ct3/at3);
	resT3 = resT3LowESE = resT3HighESE = TMath::Sqrt(at3*ct3/ac3);
      }
    }
  }
}
