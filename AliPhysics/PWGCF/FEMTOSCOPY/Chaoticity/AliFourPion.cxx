#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TF1.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"


#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"

#include "AliFourPion.h"

#define PI 3.1415927
#define G_Coeff 0.006399 // 2*pi*alpha*M_pion
#define FmToGeV 0.19733 // conversion of Fm to GeV
#define kappa3 0.15 // kappa3 Edgeworth coefficient (non-Gaussian features of C2)
#define kappa4 0.32 // kappa4 Edgeworth coefficient (non-Gaussian features of C2)
#define kappa3Fit 0.1 // kappa3 for c4QS fit
#define kappa4Fit 0.5 // kappa4 for c4QS fit

// Author: Dhevan Gangadharan

ClassImp(AliFourPion)

//________________________________________________________________________
AliFourPion::AliFourPion():
AliAnalysisTaskSE(),
  fname(0),
  fAOD(0x0), 
  fOutputList(0x0),
  fPIDResponse(0x0),
  fEC(0x0),
  fEvt(0x0),
  fTempStruct(0x0),
  fRandomNumber(0x0),
  fLEGO(kTRUE),
  fMCcase(kFALSE),
  fAODcase(kTRUE),
  fCollisionType(0),
  fGenerateSignal(kFALSE),
  fGeneratorOnly(kFALSE),
  fTabulatePairs(kFALSE),
  fOnlineCorrection(kTRUE),
  fInterpCorrection(kFALSE),
  fInterpolationType(1),
  fOneDInterpolation(0),
  fMixedChargeCut(kFALSE),
  fRMax(11),
  fRstartMC(5.0),
  ffcSq(0.7),
  ffcSqMRC(0.6),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fEAtype(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fFSIindexSmallSystem(9),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fTriggerType(0),
  fEventCounter(0),
  fEventsToMix(0),
  fEventMixingEDbins(0),
  fMultLimits(),
  fMinPt(0.16),
  fMaxPt(1.0),
  fQcut(0),
  fQLowerCut(0.005),
  fNormQcutLow(0.15),
  fNormQcutHigh(0.2),
  fKupperBound(0),
  fQupperBoundQ2(0.),
  fQupperBoundQ3(0.),
  fQupperBoundQ4(0.),
  fQbinsQ2(1),
  fQbinsQ3(1),
  fQbinsQ4(1),
  fQupperBoundWeights(0.),
  fQbinsQinv3D(0),
  fQupperBoundQinv3D(0.),
  fKstepT(),
  fKstepY(),
  fKmeanT(),
  fKmeanY(),
  fKmiddleT(),
  fKmiddleY(),
  fKmeanTOneD(),
  fQstep(0),
  fQstepWeights(0),
  fQmean(),
  fDampStart(0),
  fDampStep(0),
  fChargeSelection(kFALSE),
  fq2Binning(0),
  fLowMultBinning(0),
  fQdirectionBinning(0),
  fq2Index(0),
  fq2CutLow(0.1),
  fq2CutHigh(0.11),
  fTPCTOFboundry(0),
  fTOFboundry(0),
  fSigmaCutTPC(2.0),
  fSigmaCutTOF(2.0),
  fMinSepPairEta(0.02),
  fMinSepPairPhi(0.045),
  fShareQuality(0),
  fShareFraction(0),
  fTrueMassP(0), 
  fTrueMassPi(0), 
  fTrueMassK(0), 
  fTrueMassKs(0), 
  fTrueMassLam(0),
  fKtIndexL(0),
  fKtIndexH(0),
  fQoIndexL(0),
  fQoIndexH(0),
  fQsIndexL(0),
  fQsIndexH(0),
  fQlIndexL(0),
  fQlIndexH(0),
  fQinvIndexL(0),
  fQinvIndexH(0),
  fDummyB(0),
  fKT3transition(0.3),
  fKT4transition(0.3),
  farrP1(),
  farrP2(),
  fIC(),
  fDefaultsCharSwitch(),
  fLowQPairSwitch_E0E0(),
  fLowQPairSwitch_E0E1(),
  fLowQPairSwitch_E0E2(),
  fLowQPairSwitch_E0E3(),
  fLowQPairSwitch_E1E1(),
  fLowQPairSwitch_E1E2(),
  fLowQPairSwitch_E1E3(),
  fLowQPairSwitch_E2E3(),
  fNormQPairSwitch_E0E0(),
  fNormQPairSwitch_E0E1(),
  fNormQPairSwitch_E0E2(),
  fNormQPairSwitch_E0E3(),
  fNormQPairSwitch_E1E1(),
  fNormQPairSwitch_E1E2(),
  fNormQPairSwitch_E1E3(),
  fNormQPairSwitch_E2E3(),
  fMomResC2SC(0x0),
  fMomResC2MC(0x0),
  fWeightmuonCorrection(0x0),
  fqOutFcn(0x0),
  fqSideFcn(0x0),
  fqLongFcn(0x0)
{
  // Default constructor
  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t term=0; term<2; term++){
	    
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2=0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW = 0x0;
	    
	  }// term_2
	  
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t term=0; term<5; term++){
	      
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild = 0x0;
	      	      
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild = 0x0;
		
	      }// term_4

	    }// c4
	  }//c3
	}//c2
      }//c1
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD = 0x0;
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD = 0x0;
	}
      }
      
    }// ED
  }// Mbin
  
  // Initialze EA
  for(Int_t i=0; i<2; i++){
    fPbPbc3FitEA[i]=0x0;
    fpPbc3FitEA[i]=0x0;
    fppc3FitEA[i]=0x0;
  }

  // Initialize FSI histograms
  for(Int_t i=0; i<13; i++){
    fFSIss[i]=0x0; 
    fFSIos[i]=0x0;
  }

  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<fKbinsT; i++){// Kt iterator
    for(Int_t j=0; j<fCentBins; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
      fNormWeight2[i][j]=0x0;
      if(i==0) fNormWeightOneD[j]=0x0;
    }
  }
  
  for(Int_t FT=0; FT<2; FT++){// c3 or C3
    for(Int_t i=0; i<7; i++){// EW/LG
      for(Int_t j=0; j<50; j++){// GIndex
	ExchangeAmp[i][j][FT]=0x0;
      }
    }
  }

}
//________________________________________________________________________
AliFourPion::AliFourPion(const Char_t *name) 
  : AliAnalysisTaskSE(name), 
  fname(name),
  fAOD(0x0), 
  fOutputList(0x0),
  fPIDResponse(0x0),
  fEC(0x0),
  fEvt(0x0),
  fTempStruct(0x0),
  fRandomNumber(0x0),
  fLEGO(kTRUE),
  fMCcase(kFALSE),
  fAODcase(kTRUE),
  fCollisionType(0),
  fGenerateSignal(kFALSE),
  fGeneratorOnly(kFALSE),
  fTabulatePairs(kFALSE),
  fOnlineCorrection(kTRUE),
  fInterpCorrection(kFALSE),
  fInterpolationType(1),
  fOneDInterpolation(0),
  fMixedChargeCut(kFALSE),
  fRMax(11),
  fRstartMC(5.0),
  ffcSq(0.7),
  ffcSqMRC(0.6),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fEAtype(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fFSIindexSmallSystem(9),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fTriggerType(0),
  fEventCounter(0),
  fEventsToMix(0),
  fEventMixingEDbins(0),
  fMultLimits(),
  fMinPt(0.16),
  fMaxPt(1.0),
  fQcut(0),
  fQLowerCut(0.005),
  fNormQcutLow(0.15),
  fNormQcutHigh(0.2),
  fKupperBound(0),
  fQupperBoundQ2(0.),
  fQupperBoundQ3(0.),
  fQupperBoundQ4(0.),
  fQbinsQ2(1),
  fQbinsQ3(1),
  fQbinsQ4(1),
  fQupperBoundWeights(0.),
  fQbinsQinv3D(0),
  fQupperBoundQinv3D(0.),
  fKstepT(),
  fKstepY(),
  fKmeanT(),
  fKmeanY(),
  fKmiddleT(),
  fKmiddleY(),
  fKmeanTOneD(),
  fQstep(0),
  fQstepWeights(0),
  fQmean(),
  fDampStart(0),
  fDampStep(0),
  fChargeSelection(kFALSE),
  fq2Binning(0),
  fLowMultBinning(0),
  fQdirectionBinning(0),
  fq2Index(0),
  fq2CutLow(0.1),
  fq2CutHigh(0.11),
  fTPCTOFboundry(0),
  fTOFboundry(0),
  fSigmaCutTPC(2.0),
  fSigmaCutTOF(2.0),
  fMinSepPairEta(0.02),
  fMinSepPairPhi(0.045),
  fShareQuality(0),
  fShareFraction(0),
  fTrueMassP(0), 
  fTrueMassPi(0), 
  fTrueMassK(0), 
  fTrueMassKs(0), 
  fTrueMassLam(0),
  fKtIndexL(0),
  fKtIndexH(0),
  fQoIndexL(0),
  fQoIndexH(0),
  fQsIndexL(0),
  fQsIndexH(0),
  fQlIndexL(0),
  fQlIndexH(0),
  fQinvIndexL(0),
  fQinvIndexH(0),
  fDummyB(0),
  fKT3transition(0.3),
  fKT4transition(0.3),
  farrP1(),
  farrP2(),
  fIC(),
  fDefaultsCharSwitch(),
  fLowQPairSwitch_E0E0(),
  fLowQPairSwitch_E0E1(),
  fLowQPairSwitch_E0E2(),
  fLowQPairSwitch_E0E3(),
  fLowQPairSwitch_E1E1(),
  fLowQPairSwitch_E1E2(),
  fLowQPairSwitch_E1E3(),
  fLowQPairSwitch_E2E3(),
  fNormQPairSwitch_E0E0(),
  fNormQPairSwitch_E0E1(),
  fNormQPairSwitch_E0E2(),
  fNormQPairSwitch_E0E3(),
  fNormQPairSwitch_E1E1(),
  fNormQPairSwitch_E1E2(),
  fNormQPairSwitch_E1E3(),
  fNormQPairSwitch_E2E3(),
  fMomResC2SC(0x0),
  fMomResC2MC(0x0),
  fWeightmuonCorrection(0x0),
  fqOutFcn(0x0),
  fqSideFcn(0x0),
  fqLongFcn(0x0)
{
  // Main constructor
  fAODcase=kTRUE;
  
  

  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t term=0; term<2; term++){
	    
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2=0x0;
	    
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL = 0x0;
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW = 0x0;
	    
	  }// term_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t term=0; term<5; term++){
	      
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor = 0x0;
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild = 0x0;
	      
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild = 0x0;
		
	      }// term_4
	    }// c4
	  }//c3
	}//c2
      }//c1
      
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD = 0x0;
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD = 0x0;
	}
      }
      
    }// ED
  }// Mbin
  
  // Initialze EA
  for(Int_t i=0; i<2; i++){
    fPbPbc3FitEA[i]=0x0;
    fpPbc3FitEA[i]=0x0;
    fppc3FitEA[i]=0x0;
  }

  // Initialize FSI histograms
  for(Int_t i=0; i<13; i++){
    fFSIss[i]=0x0; 
    fFSIos[i]=0x0;
  }
  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<fKbinsT; i++){// Kt iterator
    for(Int_t j=0; j<fCentBins; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
      fNormWeight2[i][j]=0x0;
      if(i==0) fNormWeightOneD[j]=0x0;
    }
  }
  
  for(Int_t FT=0; FT<2; FT++){// c3 or C3
    for(Int_t i=0; i<7; i++){// EW/LG
      for(Int_t j=0; j<50; j++){// GIndex
	ExchangeAmp[i][j][FT]=0x0;
      }
    }
  }

  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliFourPion::AliFourPion(const AliFourPion &obj) 
  : AliAnalysisTaskSE(obj.fname),
    fname(obj.fname),
    fAOD(obj.fAOD), 
    //fESD(obj.fESD), 
    fOutputList(obj.fOutputList),
    fPIDResponse(obj.fPIDResponse),
    fEC(obj.fEC),
    fEvt(obj.fEvt),
    fTempStruct(obj.fTempStruct),
    fRandomNumber(obj.fRandomNumber),
    fLEGO(obj.fLEGO),
    fMCcase(obj.fMCcase),
    fAODcase(obj.fAODcase),
    fCollisionType(obj.fCollisionType),
    fGenerateSignal(obj.fGenerateSignal),
    fGeneratorOnly(obj.fGeneratorOnly),
    fTabulatePairs(obj.fTabulatePairs),
    fOnlineCorrection(obj.fOnlineCorrection),
    fInterpCorrection(obj.fInterpCorrection),
    fInterpolationType(obj.fInterpolationType),
    fOneDInterpolation(obj.fOneDInterpolation),
    fMixedChargeCut(obj.fMixedChargeCut),
    fRMax(obj.fRMax),
    fRstartMC(obj.fRstartMC),
    ffcSq(obj.ffcSq),
    ffcSqMRC(obj.ffcSqMRC),
    fFilterBit(obj.fFilterBit),
    fMaxChi2NDF(obj.fMaxChi2NDF),
    fMinTPCncls(obj.fMinTPCncls),
    fEAtype(obj.fEAtype),
    fBfield(obj.fBfield),
    fMbin(obj.fMbin),
    fFSIindex(obj.fFSIindex),
    fFSIindexSmallSystem(obj.fFSIindexSmallSystem),
    fEDbin(obj.fEDbin),
    fMbins(obj.fMbins),
    fMultLimit(obj.fMultLimit),
    fCentBinLowLimit(obj.fCentBinLowLimit),
    fCentBinHighLimit(obj.fCentBinHighLimit),
    fTriggerType(obj.fTriggerType),
    fEventCounter(obj.fEventCounter),
    fEventsToMix(obj.fEventsToMix),
    fEventMixingEDbins(obj.fEventMixingEDbins),
    fMultLimits(),
    fMinPt(obj.fMinPt),
    fMaxPt(obj.fMaxPt),
    fQcut(obj.fQcut),
    fQLowerCut(obj.fQLowerCut),
    fNormQcutLow(obj.fNormQcutLow),
    fNormQcutHigh(obj.fNormQcutHigh),
    fKupperBound(obj.fKupperBound),
    fQupperBoundQ2(obj.fQupperBoundQ2),
    fQupperBoundQ3(obj.fQupperBoundQ3),
    fQupperBoundQ4(obj.fQupperBoundQ4),
    fQbinsQ2(obj.fQbinsQ2),
    fQbinsQ3(obj.fQbinsQ3),
    fQbinsQ4(obj.fQbinsQ4),
    fQupperBoundWeights(obj.fQupperBoundWeights),
    fQbinsQinv3D(obj.fQbinsQinv3D),
    fQupperBoundQinv3D(obj.fQupperBoundQinv3D),
    fKstepT(),
    fKstepY(),
    fKmeanT(),
    fKmeanY(),
    fKmiddleT(),
    fKmiddleY(),
    fKmeanTOneD(),
    fQstep(obj.fQstep),
    fQstepWeights(obj.fQstepWeights),
    fQmean(),
    fDampStart(obj.fDampStart),
    fDampStep(obj.fDampStep),
    fChargeSelection(obj.fChargeSelection),
    fq2Binning(obj.fq2Binning),
    fLowMultBinning(obj.fLowMultBinning),
    fQdirectionBinning(obj.fQdirectionBinning),
    fq2Index(obj.fq2Index),
    fq2CutLow(obj.fq2CutLow),
    fq2CutHigh(obj.fq2CutHigh),
    fTPCTOFboundry(obj.fTPCTOFboundry),
    fTOFboundry(obj.fTOFboundry),
    fSigmaCutTPC(obj.fSigmaCutTPC),
    fSigmaCutTOF(obj.fSigmaCutTOF),
    fMinSepPairEta(obj.fMinSepPairEta),
    fMinSepPairPhi(obj.fMinSepPairPhi),
    fShareQuality(obj.fShareQuality),
    fShareFraction(obj.fShareFraction),
    fTrueMassP(obj.fTrueMassP), 
    fTrueMassPi(obj.fTrueMassPi), 
    fTrueMassK(obj.fTrueMassK), 
    fTrueMassKs(obj.fTrueMassKs), 
    fTrueMassLam(obj.fTrueMassLam),
    fKtIndexL(obj.fKtIndexL),
    fKtIndexH(obj.fKtIndexH),
    fQoIndexL(obj.fQoIndexL),
    fQoIndexH(obj.fQoIndexH),
    fQsIndexL(obj.fQsIndexL),
    fQsIndexH(obj.fQsIndexH),
    fQlIndexL(obj.fQlIndexL),
    fQlIndexH(obj.fQlIndexH),
    fQinvIndexL(obj.fQinvIndexL),
    fQinvIndexH(obj.fQinvIndexH),
    fDummyB(obj.fDummyB),
    fKT3transition(obj.fKT3transition),
    fKT4transition(obj.fKT4transition),
    farrP1(),
    farrP2(),
    fIC(),
    fDefaultsCharSwitch(),
    fLowQPairSwitch_E0E0(),
    fLowQPairSwitch_E0E1(),
    fLowQPairSwitch_E0E2(),
    fLowQPairSwitch_E0E3(),
    fLowQPairSwitch_E1E1(),
    fLowQPairSwitch_E1E2(),
    fLowQPairSwitch_E1E3(),
    fLowQPairSwitch_E2E3(),
    fNormQPairSwitch_E0E0(),
    fNormQPairSwitch_E0E1(),
    fNormQPairSwitch_E0E2(),
    fNormQPairSwitch_E0E3(),
    fNormQPairSwitch_E1E1(),
    fNormQPairSwitch_E1E2(),
    fNormQPairSwitch_E1E3(),
    fNormQPairSwitch_E2E3(),
    fMomResC2SC(obj.fMomResC2SC),
    fMomResC2MC(obj.fMomResC2MC),
    fWeightmuonCorrection(obj.fWeightmuonCorrection),
    fqOutFcn(obj.fqOutFcn),
    fqSideFcn(obj.fqSideFcn),
    fqLongFcn(obj.fqLongFcn)
{
  // Copy Constructor
  
  for(Int_t i=0; i<2; i++){
    fPbPbc3FitEA[i]=obj.fPbPbc3FitEA[i];
    fpPbc3FitEA[i]=obj.fpPbc3FitEA[i];
    fppc3FitEA[i]=obj.fppc3FitEA[i];
  }
  
  for(Int_t i=0; i<13; i++){
    fFSIss[i]=obj.fFSIss[i]; 
    fFSIos[i]=obj.fFSIos[i];
  }
  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<fKbinsT; i++){// Kt iterator
    for(Int_t j=0; j<fCentBins; j++){// Mbin iterator
      fNormWeight[i][j]=obj.fNormWeight[i][j];
      fNormWeight2[i][j]=obj.fNormWeight2[i][j];
      if(i==0) fNormWeightOneD[j]=obj.fNormWeightOneD[j];
    }
  }
  
  for(Int_t FT=0; FT<2; FT++){// c3 or C3
    for(Int_t i=0; i<7; i++){// EW/LG
      for(Int_t j=0; j<50; j++){// GIndex
	ExchangeAmp[i][j][FT]=obj.ExchangeAmp[i][j][FT];
      }
    }
  }
  
}
//________________________________________________________________________
AliFourPion &AliFourPion::operator=(const AliFourPion &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fname = obj.fname;
  fAOD = obj.fAOD; 
  fOutputList = obj.fOutputList;
  fPIDResponse = obj.fPIDResponse;
  fEC = obj.fEC;
  fEvt = obj.fEvt;
  fTempStruct = obj.fTempStruct;
  fRandomNumber = obj.fRandomNumber;
  fLEGO = obj.fLEGO;
  fMCcase = obj.fMCcase;
  fAODcase = obj.fAODcase;
  fCollisionType = obj.fCollisionType; 
  fGenerateSignal = obj.fGenerateSignal;
  fGeneratorOnly = obj.fGeneratorOnly;
  fTabulatePairs = obj.fTabulatePairs;
  fOnlineCorrection = obj.fOnlineCorrection;
  fInterpCorrection = obj.fInterpCorrection;
  fInterpolationType = obj.fInterpolationType;
  fOneDInterpolation = obj.fOneDInterpolation;
  fMixedChargeCut = obj.fMixedChargeCut;
  fRMax = obj.fRMax;
  fRstartMC = obj.fRstartMC;
  ffcSq = obj.ffcSq;
  ffcSqMRC = obj.ffcSqMRC;
  fFilterBit = obj.fFilterBit;
  fMaxChi2NDF = obj.fMaxChi2NDF;
  fMinTPCncls = obj.fMinTPCncls;
  fEAtype = obj.fEAtype;
  fBfield = obj.fBfield;
  fMbin = obj.fMbin;
  fFSIindex = obj.fFSIindex;
  fFSIindexSmallSystem = obj.fFSIindexSmallSystem;
  fEDbin = obj.fEDbin;
  fMbins = obj.fMbins;
  fMultLimit = obj.fMultLimit;
  fCentBinLowLimit = obj.fCentBinLowLimit;
  fCentBinHighLimit = obj.fCentBinHighLimit;
  fTriggerType = obj.fTriggerType;
  fEventCounter = obj.fEventCounter;
  fEventsToMix = obj.fEventsToMix;
  fEventMixingEDbins = obj.fEventMixingEDbins;
  fMinPt = obj.fMinPt;
  fMaxPt = obj.fMaxPt;
  fQcut = obj.fQcut;
  fQLowerCut = obj.fQLowerCut;
  fKupperBound = obj.fKupperBound;
  fQupperBoundQ2 = obj.fQupperBoundQ2;
  fQupperBoundQ3 = obj.fQupperBoundQ3;
  fQupperBoundQ4 = obj.fQupperBoundQ4;
  fQbinsQ2 = obj.fQbinsQ2;
  fQbinsQ3 = obj.fQbinsQ3;
  fQbinsQ4 = obj.fQbinsQ4;
  fQupperBoundWeights = obj.fQupperBoundWeights;
  fQbinsQinv3D = obj.fQbinsQinv3D;
  fQupperBoundQinv3D = obj.fQupperBoundQinv3D;
  fQstep = obj.fQstep;
  fQstepWeights = obj.fQstepWeights;
  fDampStart = obj.fDampStart;
  fDampStep = obj.fDampStep;
  fChargeSelection = obj.fChargeSelection;
  fq2Binning = obj.fq2Binning;
  fLowMultBinning = obj.fLowMultBinning;
  fQdirectionBinning = obj.fQdirectionBinning;
  fq2Index = obj.fq2Index;
  fq2CutLow = obj.fq2CutLow;
  fq2CutHigh = obj.fq2CutHigh;
  fTPCTOFboundry = obj.fTPCTOFboundry;
  fTOFboundry = obj.fTOFboundry;
  fSigmaCutTPC = obj.fSigmaCutTPC;
  fSigmaCutTOF = obj.fSigmaCutTOF;
  fMinSepPairEta = obj.fMinSepPairEta;
  fMinSepPairPhi = obj.fMinSepPairPhi;
  fShareQuality = obj.fShareQuality;
  fShareFraction = obj.fShareFraction;
  fTrueMassP = obj.fTrueMassP; 
  fTrueMassPi = obj.fTrueMassPi; 
  fTrueMassK = obj.fTrueMassK; 
  fTrueMassKs = obj.fTrueMassKs; 
  fTrueMassLam = obj.fTrueMassLam;
  fKtIndexL = obj.fKtIndexL;
  fKtIndexH = obj.fKtIndexH;
  fQoIndexL = obj.fQoIndexL;
  fQoIndexH = obj.fQoIndexH;
  fQsIndexL = obj.fQsIndexL;
  fQsIndexH = obj.fQsIndexH;
  fQlIndexL = obj.fQlIndexL;
  fQlIndexH = obj.fQlIndexH;
  fQinvIndexL = obj.fQinvIndexL;
  fQinvIndexH = obj.fQinvIndexH;
  fDummyB = obj.fDummyB;
  fKT3transition = obj.fKT3transition;
  fKT4transition = obj.fKT4transition;
  fMomResC2SC = obj.fMomResC2SC;
  fMomResC2MC = obj.fMomResC2MC;
  fWeightmuonCorrection = obj.fWeightmuonCorrection;
  fqOutFcn = obj.fqOutFcn;
  fqSideFcn = obj.fqSideFcn;
  fqLongFcn = obj.fqLongFcn;
  
  for(Int_t i=0; i<2; i++){
    fPbPbc3FitEA[i]=obj.fPbPbc3FitEA[i];
    fpPbc3FitEA[i]=obj.fpPbc3FitEA[i];
    fppc3FitEA[i]=obj.fppc3FitEA[i];
  }
  
  for(Int_t i=0; i<13; i++){
    fFSIss[i]=obj.fFSIss[i]; 
    fFSIos[i]=obj.fFSIos[i];
  }
  
  for(Int_t i=0; i<fKbinsT; i++){// Kt iterator
    for(Int_t j=0; j<fCentBins; j++){// Mbin iterator
      fNormWeight[i][j]=obj.fNormWeight[i][j];
      fNormWeight2[i][j]=obj.fNormWeight2[i][j];
      if(i==0) fNormWeightOneD[j]=obj.fNormWeightOneD[j];
    }
  }
  
  for(Int_t FT=0; FT<2; FT++){// c3 or C3
    for(Int_t i=0; i<7; i++){// EW/LG
      for(Int_t j=0; j<50; j++){// GIndex
	ExchangeAmp[i][j][FT]=obj.ExchangeAmp[i][j][FT];
      }
    }
  }
 
  return (*this);
}
//________________________________________________________________________
AliFourPion::~AliFourPion()
{
  // Destructor
  if(fAOD) delete fAOD; 
  //if(fESD) delete fESD; 
  if(fOutputList) delete fOutputList;
  if(fPIDResponse) delete fPIDResponse;
  if(fEC) delete fEC;
  if(fEvt) delete fEvt;
  if(fTempStruct) delete [] fTempStruct;
  if(fRandomNumber) delete fRandomNumber;
  if(fMomResC2SC) delete fMomResC2SC;
  if(fMomResC2MC) delete fMomResC2MC;
  if(fWeightmuonCorrection) delete fWeightmuonCorrection;
  
  for(Int_t i=0; i<2; i++){
    if(fPbPbc3FitEA[i]) delete fPbPbc3FitEA[i];
    if(fpPbc3FitEA[i]) delete fpPbc3FitEA[i];
    if(fppc3FitEA[i]) delete fppc3FitEA[i];
  }
  
  for(Int_t j=0; j<kMultLimitPbPb; j++){
    if(fLowQPairSwitch_E0E0[j]) delete [] fLowQPairSwitch_E0E0[j];
    if(fLowQPairSwitch_E0E1[j]) delete [] fLowQPairSwitch_E0E1[j];
    if(fLowQPairSwitch_E0E2[j]) delete [] fLowQPairSwitch_E0E2[j];
    if(fLowQPairSwitch_E0E3[j]) delete [] fLowQPairSwitch_E0E3[j];
    if(fLowQPairSwitch_E1E1[j]) delete [] fLowQPairSwitch_E1E1[j];
    if(fLowQPairSwitch_E1E2[j]) delete [] fLowQPairSwitch_E1E2[j];
    if(fLowQPairSwitch_E1E3[j]) delete [] fLowQPairSwitch_E1E3[j];
    if(fLowQPairSwitch_E2E3[j]) delete [] fLowQPairSwitch_E2E3[j];
    //
    if(fNormQPairSwitch_E0E0[j]) delete [] fNormQPairSwitch_E0E0[j];
    if(fNormQPairSwitch_E0E1[j]) delete [] fNormQPairSwitch_E0E1[j];
    if(fNormQPairSwitch_E0E2[j]) delete [] fNormQPairSwitch_E0E2[j];
    if(fNormQPairSwitch_E0E3[j]) delete [] fNormQPairSwitch_E0E3[j];
    if(fNormQPairSwitch_E1E1[j]) delete [] fNormQPairSwitch_E1E1[j];
    if(fNormQPairSwitch_E1E2[j]) delete [] fNormQPairSwitch_E1E2[j];
    if(fNormQPairSwitch_E1E3[j]) delete [] fNormQPairSwitch_E1E3[j];
    if(fNormQPairSwitch_E2E3[j]) delete [] fNormQPairSwitch_E2E3[j];
  }
  
  //
  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t term=0; term<2; term++){
	    
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fBuild) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fBuild;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW;
	    //
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv;
	    if(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW) delete Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW;
	  }// term_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t term=0; term<5; term++){
		
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3;
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3;
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor;
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted;
	      //
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild;
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms33D) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms33D;
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor3D) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor3D;
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted;
		//
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuild) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuild;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuild) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuild;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuild) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuild;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildFromFits) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildFromFits;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildFromFits) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildFromFits;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildFromFits) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildFromFits;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildFromFits) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildFromFits;
	      }// term_4

	    }//c4
	  }//c3
	}//c2
      }//c1
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  if(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD) delete KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD;
	  if(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD) delete KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD;
	}
      }
      
    }// ED
  }// Mbin
  
   
  for(Int_t i=0; i<13; i++){
    if(fFSIss[i]) delete fFSIss[i]; 
    if(fFSIos[i]) delete fFSIos[i];
  }
  
  for(Int_t i=0; i<fKbinsT; i++){// Kt iterator
    for(Int_t j=0; j<fCentBins; j++){// Mbin iterator
      if(fNormWeight[i][j]) delete fNormWeight[i][j];
      if(fNormWeight2[i][j]) delete fNormWeight2[i][j];
      if(i==0) {if(fNormWeightOneD[j]) delete fNormWeightOneD[j];}
    }
  }
  
  for(Int_t FT=0; FT<2; FT++){// c3 or C3
    for(Int_t i=0; i<7; i++){// EW/LG
      for(Int_t j=0; j<50; j++){// GIndex
	if(ExchangeAmp[i][j][FT]) delete ExchangeAmp[i][j][FT];
      }
    }
  }
  if(fqOutFcn) delete fqOutFcn;
  if(fqSideFcn) delete fqSideFcn;
  if(fqLongFcn) delete fqLongFcn;
 
}
//________________________________________________________________________
void AliFourPion::ParInit()
{
  cout<<"AliFourPion MyInit() call"<<endl;
  cout<<"lego:"<<fLEGO<<"  MCcase:"<<fMCcase<<"  CollisionType:"<<fCollisionType<<"  TabulatePairs:"<<fTabulatePairs<<"  GenSignal:"<<fGenerateSignal<<"  CentLow:"<<fCentBinLowLimit<<"  CentHigh:"<<fCentBinHighLimit<<"  RMax:"<<fRMax<<"  fc^2:"<<ffcSq<<"  FB:"<<fFilterBit<<"  MaxChi2/NDF:"<<fMaxChi2NDF<<"  MinTPCncls:"<<fMinTPCncls<<"  MinPairSepEta:"<<fMinSepPairEta<<"  MinPairSepPhi:"<<fMinSepPairPhi<<"  NsigTPC:"<<fSigmaCutTPC<<"  NsigTOF:"<<fSigmaCutTOF<<endl;
  
  fRandomNumber = new TRandom3();
  fRandomNumber->SetSeed(0);
  
  //
  fEventCounter=0;
  fEventsToMix=3;
  if(fq2Binning) fEventMixingEDbins=9;// was 81 with q1 and q2 binning
  else fEventMixingEDbins=2;// for both Z-vertex binning and LowMultBinning

  if(fMCcase) fEventMixingEDbins=2;
  
  fTPCTOFboundry = 0.6;// TPC pid used below this momentum, TOF above but below TOF_boundry
  fTOFboundry = 2.1;// TOF pid used below this momentum
  
  ////////////////////////////////////////////////
  // PadRow Pair Cuts
  fShareQuality = .5;// max
  fShareFraction = .05;// max
  ////////////////////////////////////////////////
  
  // pp and pPb mult limits
  fMultLimits[0]=0, fMultLimits[1]=5; fMultLimits[2]=10; fMultLimits[3]=15; fMultLimits[4]=20;
  fMultLimits[5]=30, fMultLimits[6]=40; fMultLimits[7]=50; fMultLimits[8]=70; fMultLimits[9]=100;
  fMultLimits[10]=150;
  
  
  
  if(fCollisionType==0) {// PbPb
    fMultLimit=kMultLimitPbPb;
    fMbins=fCentBins;
    fQcut=0.1;
    fRstartMC = 5.0;
    fQbinsQinv3D = 20;
    fQupperBoundQinv3D = 0.1;
    fQupperBoundWeights = 0.2;
  }else {// pPb & pp
    fMultLimit=kMultLimitpp; 
    fMbins=1; 
    fQcut=0.5;
    fRstartMC = 1.0;
    fQbinsQinv3D = 60;
    fQupperBoundQinv3D = 0.6;
    fQupperBoundWeights = 0.4;
  }
  
  
  fKupperBound = 1.0;
  //
  fKstepY[0] = 1.6;
  fKmeanY[0] = 0;// central y
  fKmiddleY[0] = 0;

  // 7x1 (Kt: 0-0.2, 0.2-0.25, 0.25-0.3, 0.3-0.35, 0.35-0.4, 0.4-0.45, 0.45-1.0)
  if(fKbinsT==7){
    fKstepT[0] = 0.2; fKstepT[1] = 0.05; fKstepT[2] = 0.05; fKstepT[3] = 0.05; fKstepT[4] = 0.05; fKstepT[5] = 0.05; fKstepT[6] = 0.55;
    fKmeanT[0] = 0.188; fKmeanT[1] = 0.227; fKmeanT[2] = 0.275; fKmeanT[3] = 0.324; fKmeanT[4] = 0.374; fKmeanT[5] = 0.424; fKmeanT[6] = 0.552; 
    fKmiddleT[0] = 0.1; fKmiddleT[1] = 0.225; fKmiddleT[2] = 0.275; fKmiddleT[3] = 0.325; fKmiddleT[4] = 0.375; fKmiddleT[5] = 0.425; fKmiddleT[6] = 0.725;
  }
  // 6x1 (Kt: 0-0.2, 0.2-0.24, 0.24-0.3, 0.3-0.35, 0.35-0.45, 0.45-1.0)
  if(fKbinsT==6){
    fKstepT[0] = 0.2; fKstepT[1] = 0.04; fKstepT[2] = 0.06; fKstepT[3] = 0.05; fKstepT[4] = 0.1; fKstepT[5] = 0.55;
    fKmeanT[0] = 0.188; fKmeanT[1] = 0.222; fKmeanT[2] = 0.270; fKmeanT[3] = 0.324; fKmeanT[4] = 0.395; fKmeanT[5] = 0.551; 
    fKmiddleT[0] = 0.1; fKmiddleT[1] = 0.22; fKmiddleT[2] = 0.27; fKmiddleT[3] = 0.325; fKmiddleT[4] = 0.4; fKmiddleT[5] = 0.725;
  }
  // 4x1 (Kt: 0-0.25, 0.25-0.35, 0.35-0.45, 0.45-1.0)
  if(fKbinsT==4){
    fKstepT[0] = 0.25; fKstepT[1] = 0.1; fKstepT[2] = 0.1; fKstepT[3] = 0.55;
    fKmeanT[0] = 0.212; fKmeanT[1] = 0.299; fKmeanT[2] = 0.398; fKmeanT[3] = 0.576;
    fKmiddleT[0] = 0.125; fKmiddleT[1] = 0.3; fKmiddleT[2] = 0.4; fKmiddleT[3] = 0.725;
  }
  // 3x1 (Kt: 0-0.3, 0.3-0.45, 0.45-1.0)
  if(fKbinsT==3){
    fKstepT[0] = 0.3; fKstepT[1] = 0.15; fKstepT[2] = 0.55;
    fKmeanT[0] = 0.240; fKmeanT[1] = 0.369; fKmeanT[2] = 0.576;
    fKmiddleT[0] = 0.15; fKmiddleT[1] = 0.375; fKmiddleT[2] = 0.725;
  }
  // 2x1 (Kt: 0-0.35, 0.35-1.0)
  if(fKbinsT==2){
    fKstepT[0] = 0.35; fKstepT[1] = 0.65;
    fKmeanT[0] = 0.264; fKmeanT[1] = 0.500;
    fKmiddleT[0] = 0.175; fKmiddleT[1] = 0.675;
  }
  //
  // Qinv kT bins
  for(Int_t i=0; i<fKbinsTOneD; i++){
    fKmeanTOneD[i] = 0.16 + (i+0.5)*0.03;
  }
  
 
  //
  fQupperBoundQ2 = 2.0;
  fQupperBoundQ3 = 0.6;
  fQupperBoundQ4 = 0.6;
  fQbinsQ2 = int(fQupperBoundQ2/0.005);
  fQbinsQ3 = int(fQupperBoundQ3/0.005);
  fQbinsQ4 = int(fQupperBoundQ4/0.005);
  fQstepWeights = fQupperBoundWeights/Float_t(kQbinsWeights);
  for(Int_t i=0; i<kQbinsWeights; i++) {fQmean[i]=(i+0.5)*fQstepWeights;}
  
  //
  fDampStart = 0.5;// was 0.3, then 0.5
  fDampStep = 0.02;
  
  //
  
  fEC = new AliFourPionEventCollection **[fEventMixingEDbins];
  for(UShort_t i=0; i<fEventMixingEDbins; i++){
    
    fEC[i] = new AliFourPionEventCollection *[fMbinsMixing];

    for(UShort_t j=0; j<fMbinsMixing; j++){
      
      fEC[i][j] = new AliFourPionEventCollection(fEventsToMix+1, fMultLimit, kMCarrayLimit, fMCcase);
    }
  }
  
  for(Int_t i=0; i<kMultLimitPbPb; i++) fDefaultsCharSwitch[i]='0';
  for(Int_t i=0; i<kMultLimitPbPb; i++) {
    fLowQPairSwitch_E0E0[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E1[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E2[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E1[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E2[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E2E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    //
    fNormQPairSwitch_E0E0[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E1[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E2[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E1[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E2[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E2E3[i] = new TArrayC(kMultLimitPbPb,fDefaultsCharSwitch);
  }
  
  fTempStruct = new AliFourPionTrackStruct[fMultLimit];
  
  
  fTrueMassP=0.93827, fTrueMassPi=0.13957, fTrueMassK=0.493677, fTrueMassKs=0.497614, fTrueMassLam=1.11568;
  

  // Set weights, Coulomb corrections, etc. if not in LEGO train
  if(!fLEGO) {
    SetFSICorrelations(fLEGO);// Read in 2-particle and 3-particle FSI correlations
    if(!fTabulatePairs) SetWeightArrays(fLEGO);// Set Weight Array
    if(!fMCcase && !fTabulatePairs) SetMomResCorrections(fLEGO);// Read Momentum resolution file
    if(!fMCcase && !fTabulatePairs) SetMuonCorrections(fLEGO);// Read Muon corrections
    if(!fMCcase && !fTabulatePairs) Setc3FitEAs(fLEGO);// Read EAs from c3 fits
  }
  


  // Pair-Exchange amplitudes from c3 fits
  TString *EWequation = new TString("[0]*exp(-pow(x*[1]/0.19733,2)/2.) * ( 1 + [2]/(6.*pow(2.,1.5))*(8.*pow(x*[1]/0.19733,3) - 12.*pow(x*[1]/0.19733,1)) + [3]/(24.*pow(2.,2))*(16.*pow(x*[1]/0.19733,4) -48.*pow(x*[1]/0.19733,2) + 12) + [4]/(120.*pow(2.,2.5))*(32.*pow(x*[1]/0.19733,5) - 160.*pow(x*[1]/0.19733,3) + 120.*x*[1]/0.19733) + [5]/(720.*pow(2.,3.))*(64.*pow(x*[1]/0.19733,6) - 480.*pow(x*[1]/0.19733,4) + 720.*pow(x*[1]/0.19733,2) - 120.))");
  TString *LGequation = new TString("[0]*exp(-x*[1]/0.19733/2.) * ( 1 + [2]*(-x*[1]/0.19733 + 1) + [3]/2.*(pow(x*[1]/0.19733,2) - 4.*x*[1]/0.19733 + 2) + [4]/6.*(-pow(x*[1]/0.19733,3) + 9.*pow(x*[1]/0.19733,2) - 18.*x*[1]/0.19733 + 6) + [5]/24.*(pow(x*[1]/0.19733,4) - 16.*pow(x*[1]/0.19733,3) + 72.*pow(x*[1]/0.19733,2) - 96.*x*[1]/0.19733 +24.))");
  
  if(!fMCcase && !fTabulatePairs){
    for(Int_t FT=0; FT<2; FT++){// c3 or C3
      for(Int_t i=0; i<7; i++){// Rcoh index
	for(Int_t j=0; j<50; j++){// G index
	  TString *nameEA=new TString("ExchangeAmp");
	  *nameEA += FT;
	  *nameEA += i;
	  *nameEA += j;
	  if(fEAtype==0) ExchangeAmp[i][j][FT] = new TF1(nameEA->Data(), EWequation->Data(), 0,1.0);// Edgeworth
	  else ExchangeAmp[i][j][FT] = new TF1(nameEA->Data(), LGequation->Data(), 0,1.0);// Laguerre
	  //
	  if(fCollisionType==0){
	    ExchangeAmp[i][j][FT]->FixParameter(0, fPbPbc3FitEA[FT]->GetBinContent(i+1, 1, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(1, fPbPbc3FitEA[FT]->GetBinContent(i+1, 2, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(2, fPbPbc3FitEA[FT]->GetBinContent(i+1, 3, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(3, fPbPbc3FitEA[FT]->GetBinContent(i+1, 4, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(4, fPbPbc3FitEA[FT]->GetBinContent(i+1, 5, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(5, fPbPbc3FitEA[FT]->GetBinContent(i+1, 6, j+1));
	  }else if(fCollisionType==1){
	    ExchangeAmp[i][j][FT]->FixParameter(0, fpPbc3FitEA[FT]->GetBinContent(i+1, 1, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(1, fpPbc3FitEA[FT]->GetBinContent(i+1, 2, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(2, fpPbc3FitEA[FT]->GetBinContent(i+1, 3, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(3, fpPbc3FitEA[FT]->GetBinContent(i+1, 4, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(4, fpPbc3FitEA[FT]->GetBinContent(i+1, 5, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(5, fpPbc3FitEA[FT]->GetBinContent(i+1, 6, j+1));
	  }else{
	    ExchangeAmp[i][j][FT]->FixParameter(0, fppc3FitEA[FT]->GetBinContent(i+1, 1, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(1, fppc3FitEA[FT]->GetBinContent(i+1, 2, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(2, fppc3FitEA[FT]->GetBinContent(i+1, 3, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(3, fppc3FitEA[FT]->GetBinContent(i+1, 4, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(4, fppc3FitEA[FT]->GetBinContent(i+1, 5, j+1));
	    ExchangeAmp[i][j][FT]->FixParameter(5, fppc3FitEA[FT]->GetBinContent(i+1, 6, j+1));
	  }
	}
      }
    }
  }

  // Interpolation correction factors for 4D 7 kT binning scheme with cubic interpolation
  // M0
  Float_t fIC_tempK1_Offline_M0[20]={1.15874, 1.15874, 1.02972, 0.983506, 0.978256, 0.978137, 0.99731, 0.988222, 0.976087, 1.01034, 1.00952, 1.05129, 1.33569, 1.57477, 3.00402, 7.37385, 20.8854, 79.9247, 466.958, 2362.59};
  Float_t fIC_tempK2_Offline_M0[20]={1.17498, 1.17498, 1.04839, 1.01326, 1.00137, 1.00125, 1.00243, 0.99748, 0.992767, 1.00106, 0.992621, 0.983214, 1.1318, 1.40818, 2.19906, 5.2095, 16.6987, 61.8091, 220.492, 1356.37};
  Float_t fIC_tempK3_Offline_M0[20]={1.19759, 1.19759, 1.07763, 1.0422, 1.03393, 1.02126, 1.02089, 1.02327, 1.0085, 1.00487, 0.971789, 0.993542, 1.0836, 1.44404, 2.5961, 6.18613, 23.6349, 92.2967, 423.414, 2252.86};
  Float_t fIC_tempK4_Offline_M0[20]={1.20016, 1.20016, 1.08142, 1.06597, 1.05223, 1.04274, 1.02726, 1.01956, 1.01569, 1.00886, 0.976751, 1.02217, 1.11237, 1.6687, 3.11013, 8.07875, 33.3551, 137.644, 776.742, 6874.55};
  Float_t fIC_tempK5_Offline_M0[20]={1.17724, 1.17724, 1.07073, 1.0894, 1.06995, 1.04995, 1.03323, 1.02118, 1.00556, 0.993756, 0.989652, 1.03069, 1.17341, 1.76291, 3.57525, 5.23737, 38.1559, 176.201, 837.799, 7824.81};
  Float_t fIC_tempK6_Offline_M0[20]={1.11908, 1.11908, 1.07268, 1.10707, 1.0777, 1.05077, 1.02427, 1.01298, 1.0043, 0.995263, 0.971299, 1.01614, 1.24185, 1.81242, 3.63278, 11.1267, 37.3259, 177.369, 913.754, 8618.34};
  Float_t fIC_tempK7_Offline_M0[20]={1.17818, 1.17818, 1.13705, 1.13692, 1.11537, 1.07695, 1.04205, 1.02878, 1.02658, 1.02337, 1.0331, 1.10945, 1.37393, 2.01098, 4.10517, 10.6846, 41.204, 196.267, 946.015, 10840.9};
  // M2
  Float_t fIC_tempK1_Offline_M2[20]={1.13123, 1.13123, 1.02803, 0.994159, 0.980033, 0.989634, 0.992672, 1.02432, 1.05428, 1.04889, 1.07725, 1.07841, 1.09708, 1.07032, 1.12155, 1.44804, 3.49901, 6.48294, 0.975108, 0.946954};
  Float_t fIC_tempK2_Offline_M2[20]={1.14932, 1.14932, 1.04885, 1.02059, 1.00685, 1.00726, 1.00741, 1.00643, 0.995777, 0.98466, 0.947865, 0.924781, 0.890824, 0.867119, 0.881113, 0.957089, 2.51592, 3.9161, 1.01534, 0.625349};
  Float_t fIC_tempK3_Offline_M2[20]={1.15592, 1.15592, 1.05566, 1.03865, 1.02689, 1.02093, 1.00489, 0.997721, 0.989688, 0.963098, 0.886343, 0.897781, 0.791003, 0.750446, 0.748841, 1.01915, 3.96752, 1.90596, 0.523598, 0.414215};
  Float_t fIC_tempK4_Offline_M2[20]={1.14644, 1.14644, 1.06822, 1.05527, 1.04913, 1.03348, 1.02193, 0.981811, 0.959563, 0.9284, 0.883132, 0.835031, 0.780229, 0.740518, 0.732312, 1.21148, 1.91578, 0.817937, 0.39292, 0.357011};
  Float_t fIC_tempK5_Offline_M2[20]={1.11776, 1.11776, 1.0404, 1.06973, 1.05121, 1.03497, 1.02112, 0.981788, 0.960187, 0.929131, 0.8798, 0.827386, 0.78019, 0.811157, 0.864437, 1.07885, 1.10748, 0.660085, 0.468907, 0.26552};
  Float_t fIC_tempK6_Offline_M2[20]={1.10626, 1.10626, 1.04509, 1.08403, 1.07131, 1.03575, 1.01936, 0.970014, 0.950104, 0.931275, 0.902527, 0.844336, 0.834994, 0.815202, 0.810353, 0.943186, 1.09361, 0.618069, 0.396422, 0.339816};
  Float_t fIC_tempK7_Offline_M2[20]={1.14824, 1.14824, 1.11668, 1.11007, 1.09798, 1.07503, 1.03675, 1.00177, 0.963043, 0.941321, 0.919049, 0.881138, 0.888227, 0.872429, 0.90361, 0.966787, 1.03367, 0.704508, 0.478153, 0.425762};
  // M4
  Float_t fIC_tempK1_Offline_M4[20]={1.1374, 1.1374, 1.02185, 0.997306, 0.999296, 0.993486, 1.02246, 1.01117, 1.07119, 1.08395, 1.13366, 1.17204, 1.21245, 1.17899, 1.26667, 1.0581, 1.50895, 1.95813, 3.65374, 1.15169};
  Float_t fIC_tempK2_Offline_M4[20]={1.13045, 1.13045, 1.03681, 1.0128, 1.00759, 1.01019, 0.992187, 1.00156, 1.0025, 0.989377, 0.968774, 0.955621, 0.93443, 0.862969, 0.837983, 0.818954, 0.978139, 1.18539, 2.63576, 0.789694};
  Float_t fIC_tempK3_Offline_M4[20]={1.12677, 1.12677, 1.05576, 1.02281, 1.01658, 1.01044, 0.983065, 0.974013, 0.959087, 0.932451, 0.887496, 0.816551, 0.788911, 0.798848, 0.672707, 0.596951, 0.724822, 0.985777, 1.55475, 0.512027};
  Float_t fIC_tempK4_Offline_M4[20]={1.1037, 1.1037, 1.04596, 1.03772, 1.02463, 1.00747, 0.987855, 0.967885, 0.931725, 0.892358, 0.834885, 0.801171, 0.750001, 0.715771, 0.654448, 0.573221, 0.667588, 0.853278, 0.653164, 0.360429};
  Float_t fIC_tempK5_Offline_M4[20]={1.15762, 1.15762, 1.02482, 1.05618, 1.02964, 1.01251, 0.988161, 0.966585, 0.940479, 0.864523, 0.850332, 0.755725, 0.742726, 0.688804, 0.688238, 0.621091, 0.732821, 0.696838, 0.486071, 0.295078};
  Float_t fIC_tempK6_Offline_M4[20]={1.02007, 1.02007, 1.00296, 1.05498, 1.042, 1.03672, 0.987081, 0.944023, 0.922929, 0.878624, 0.817037, 0.768768, 0.703576, 0.704773, 0.660366, 0.690419, 0.612527, 0.518631, 0.457327, 0.273755};
  Float_t fIC_tempK7_Offline_M4[20]={0.945336, 0.945336, 1.07273, 1.0781, 1.07007, 1.04254, 1.00869, 0.963936, 0.942431, 0.900622, 0.874, 0.83447, 0.81639, 0.741021, 0.703714, 0.655082, 0.68339, 0.682747, 0.535062, 0.355305};
  // M6
  Float_t fIC_tempK1_Offline_M6[20]={1.15097, 1.15097, 1.01525, 1.02325, 1.02101, 1.00372, 1.00764, 1.02211, 1.03287, 1.05653, 1.13955, 1.1666, 1.15742, 1.03046, 1.37421, 1.18267, 1.30114, 1.1483, 1.48882, 1.73755};
  Float_t fIC_tempK2_Offline_M6[20]={1.12252, 1.12252, 1.02312, 1.00526, 1.01376, 0.998435, 0.996786, 0.995867, 0.993003, 0.982585, 0.938568, 0.947641, 0.933897, 0.974274, 0.831059, 0.758825, 0.847514, 0.791063, 0.858551, 1.05155};
  Float_t fIC_tempK3_Offline_M6[20]={1.07383, 1.07383, 1.02249, 1.01156, 1.01327, 1.00279, 0.990091, 0.977096, 0.954393, 0.907036, 0.870537, 0.835714, 0.821147, 0.703437, 0.69385, 0.701033, 0.633508, 0.632785, 0.524889, 0.661071};
  Float_t fIC_tempK4_Offline_M6[20]={1.13541, 1.13541, 1.02829, 1.01309, 1.01646, 1.00418, 0.991764, 0.973578, 0.926166, 0.91091, 0.848218, 0.819249, 0.730999, 0.70245, 0.631557, 0.575834, 0.566385, 0.538136, 0.550584, 0.710427};
  Float_t fIC_tempK5_Offline_M6[20]={1.12525, 1.12525, 0.997262, 1.03232, 1.01936, 1.01814, 0.972475, 0.947976, 0.929823, 0.860047, 0.807019, 0.759762, 0.775143, 0.652177, 0.601297, 0.587284, 0.574577, 0.551835, 0.636213, 1.02383};
  Float_t fIC_tempK6_Offline_M6[20]={1.07832, 1.07832, 0.999026, 1.05192, 1.03598, 1.00859, 0.983789, 0.950884, 0.914573, 0.838778, 0.806613, 0.761976, 0.713618, 0.670622, 0.657188, 0.651494, 0.685313, 0.616319, 0.722241, 1.07535};
  Float_t fIC_tempK7_Offline_M6[20]={0.982898, 0.982898, 1.02901, 1.06818, 1.05011, 1.01388, 0.99956, 0.94955, 0.924168, 0.889595, 0.813131, 0.823018, 0.782903, 0.751848, 0.725681, 0.728617, 0.832607, 0.805799, 0.926431, 1.09914};
  // M8
  Float_t fIC_tempK1_Offline_M8[20]={-1.14881, -1.14881, 1.01761, 1.00793, 0.998043, 1.02086, 1.04191, 1.04094, 1.02651, 1.02699, 1.06824, 1.056, 1.2019, 1.08276, 1.16456, 1.18241, 1.25049, 1.22569, 1.63738, 1.05058};
  Float_t fIC_tempK2_Offline_M8[20]={-1.12199, -1.12199, 1.02699, 1.00251, 0.99294, 1.00036, 0.982521, 0.980228, 0.977508, 0.961233, 0.964523, 0.940707, 0.925457, 0.870468, 0.951827, 0.747598, 0.871441, 0.811534, 0.855158, 0.965166};
  Float_t fIC_tempK3_Offline_M8[20]={-1.13065, -1.13065, 1.02997, 1.00678, 0.99411, 0.996577, 0.969551, 0.959171, 0.940933, 0.896582, 0.886173, 0.858992, 0.855163, 0.832783, 0.742629, 0.712666, 0.675089, 0.553073, 0.597996, 0.740062};
  Float_t fIC_tempK4_Offline_M8[20]={-1.06918, -1.06918, 1.05707, 1.00292, 0.991253, 0.97448, 0.974132, 0.951746, 0.906052, 0.907403, 0.869378, 0.825885, 0.780601, 0.722887, 0.694694, 0.674391, 0.628061, 0.628201, 0.821505, 1.05173};
  Float_t fIC_tempK5_Offline_M8[20]={-0.99603, -0.99603, 0.966565, 1.00141, 0.996606, 0.971611, 0.945016, 0.939809, 0.928412, 0.900481, 0.866328, 0.771432, 0.766459, 0.742102, 0.665262, 0.642106, 0.628457, 0.637371, 0.781444, 1.2608};
  Float_t fIC_tempK6_Offline_M8[20]={-0.922137, -0.922137, 0.992591, 1.01865, 1.0108, 0.987423, 0.96401, 0.896552, 0.886862, 0.895748, 0.833005, 0.774291, 0.735689, 0.672574, 0.659043, 0.575224, 0.666163, 0.782534, 0.808509, 0.991679};
  Float_t fIC_tempK7_Offline_M8[20]={-0.93778, -0.93778, 1.0315, 1.02039, 1.0344, 1.0235, 0.980848, 0.928473, 0.903296, 0.880838, 0.848799, 0.827534, 0.799069, 0.743353, 0.761462, 0.794397, 0.855292, 0.867931, 1.0435, 1.07612};
  //
  //
  /*Float_t fIC_tempK1_Online[20]={1.08977, 1.08977, 1.03018, 1.02195, 1.02445, 1.01626, 1.02466, 1.02598, 1.01903, 1.04232, 1.04898, 1.09417, 1.292, 1.72562, 3.46667, 9.12937, 31.0891, 147.542, 825.509, 6269.37};
  Float_t fIC_tempK2_Online[20]={1.10279, 1.10279, 1.02865, 1.02268, 1.01911, 1.01757, 1.01651, 1.01223, 1.01343, 1.01847, 1.02396, 1.0222, 1.19166, 1.52826, 2.5289, 6.95869, 26.9192, 107.863, 505.356, 5256.24};
  Float_t fIC_tempK3_Online[20]={1.13141, 1.13141, 1.02827, 1.00963, 1.00813, 0.999384, 1.00072, 0.997054, 0.97849, 0.970826, 0.940245, 0.951111, 1.01939, 1.35676, 2.32194, 5.87639, 22.3804, 86.5628, 482.493, 4202.39};
  Float_t fIC_tempK4_Online[20]={1.14242, 1.14242, 1.01967, 1.00638, 0.999392, 0.996161, 0.984372, 0.97497, 0.966034, 0.929085, 0.901963, 0.928844, 0.961555, 1.30751, 2.1384, 5.31829, 19.8114, 78.2922, 396.502, 3743.63};
  Float_t fIC_tempK5_Online[20]={1.15841, 1.15841, 0.999762, 1.0077, 0.99768, 0.992972, 0.989344, 0.972073, 0.951255, 0.935728, 0.907556, 0.907257, 0.976454, 1.23242, 2.19768, 5.22021, 18.728, 78.2627, 354.761, 3122.74};
  Float_t fIC_tempK6_Online[20]={1.14786, 1.14786, 1.00226, 1.00958, 0.995811, 0.989573, 0.978303, 0.959762, 0.951294, 0.928133, 0.89688, 0.903638, 0.966431, 1.18643, 2.05014, 5.4632, 17.0296, 72.0803, 375.084, 3251.85};
  Float_t fIC_tempK7_Online[20]={1.22518, 1.22518, 1.04616, 1.02128, 1.0144, 1.00255, 0.986872, 0.97028, 0.960213, 0.940531, 0.923523, 0.940603, 1.02054, 1.31177, 2.21925, 5.16703, 17.5741, 68.9239, 315.734, 2797.83};*/
  //
  for(Int_t i=0; i<20; i++){
    for(Int_t j=0; j<5; j++){
      if(j==0){
	fIC[j][0][i] = fIC_tempK1_Offline_M0[i];
	fIC[j][1][i] = fIC_tempK2_Offline_M0[i];
	fIC[j][2][i] = fIC_tempK3_Offline_M0[i];
	fIC[j][3][i] = fIC_tempK4_Offline_M0[i];
	fIC[j][4][i] = fIC_tempK5_Offline_M0[i];
	fIC[j][5][i] = fIC_tempK6_Offline_M0[i];
	fIC[j][6][i] = fIC_tempK7_Offline_M0[i];
      }else if(j==1){
	fIC[j][0][i] = fIC_tempK1_Offline_M2[i];
	fIC[j][1][i] = fIC_tempK2_Offline_M2[i];
	fIC[j][2][i] = fIC_tempK3_Offline_M2[i];
	fIC[j][3][i] = fIC_tempK4_Offline_M2[i];
	fIC[j][4][i] = fIC_tempK5_Offline_M2[i];
	fIC[j][5][i] = fIC_tempK6_Offline_M2[i];
	fIC[j][6][i] = fIC_tempK7_Offline_M2[i];
      }else if(j==2){
	fIC[j][0][i] = fIC_tempK1_Offline_M4[i];
	fIC[j][1][i] = fIC_tempK2_Offline_M4[i];
	fIC[j][2][i] = fIC_tempK3_Offline_M4[i];
	fIC[j][3][i] = fIC_tempK4_Offline_M4[i];
	fIC[j][4][i] = fIC_tempK5_Offline_M4[i];
	fIC[j][5][i] = fIC_tempK6_Offline_M4[i];
	fIC[j][6][i] = fIC_tempK7_Offline_M4[i];
      }else if(j==3){
	fIC[j][0][i] = fIC_tempK1_Offline_M6[i];
	fIC[j][1][i] = fIC_tempK2_Offline_M6[i];
	fIC[j][2][i] = fIC_tempK3_Offline_M6[i];
	fIC[j][3][i] = fIC_tempK4_Offline_M6[i];
	fIC[j][4][i] = fIC_tempK5_Offline_M6[i];
	fIC[j][5][i] = fIC_tempK6_Offline_M6[i];
	fIC[j][6][i] = fIC_tempK7_Offline_M6[i];
      }else{
	fIC[j][0][i] = fIC_tempK1_Offline_M8[i];
	fIC[j][1][i] = fIC_tempK2_Offline_M8[i];
	fIC[j][2][i] = fIC_tempK3_Offline_M8[i];
	fIC[j][3][i] = fIC_tempK4_Offline_M8[i];
	fIC[j][4][i] = fIC_tempK5_Offline_M8[i];
	fIC[j][5][i] = fIC_tempK6_Offline_M8[i];
	fIC[j][6][i] = fIC_tempK7_Offline_M8[i];
      }
    }
  }
  
  // mean +- RMS
  fqOutFcn = new TF1("fqOutFcn","[0] + [1]*x + ([2] + [3]*x)",0,1.0);
  fqOutFcn->FixParameter(0, 5.25711e-03);// mean
  fqOutFcn->FixParameter(1, 2.07835);// mean
  fqOutFcn->FixParameter(2, -1.82312e-03);// rms
  fqOutFcn->FixParameter(3, 8.74446e-01);// rms
  fqSideFcn = new TF1("fqSideFcn","[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3) + ([4] + [5]*x)",0,1.0);
  fqSideFcn->FixParameter(0, 1.47222e-02);// mean
  fqSideFcn->FixParameter(1, 8.42044e-01);// mean
  fqSideFcn->FixParameter(2, 2.82436);// mean
  fqSideFcn->FixParameter(3, -5.47032);// mean
  fqSideFcn->FixParameter(4, 1.92920e-03);// rms
  fqSideFcn->FixParameter(5, 3.91368e-01);// rms
  fqLongFcn = new TF1("fqLongFcn","[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3) + ([4] + [5]*x)",0,1.0);
  fqLongFcn->FixParameter(0, 1.94838e-02);// mean 
  fqLongFcn->FixParameter(1, 9.22609e-01);// mean
  fqLongFcn->FixParameter(2, 1.32699);// mean
  fqLongFcn->FixParameter(3, -1.24812);// mean
  fqLongFcn->FixParameter(4, -2.23316e-04);// rms
  fqLongFcn->FixParameter(5, 3.98127e-01);// rms
  
  cout<<"End ParInit"<<endl;
  /////////////////////////////////////////////
  /////////////////////////////////////////////
 
}
//________________________________________________________________________
void AliFourPion::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  ParInit();// Initialize my settings

  cout<<"Start histogram creation"<<endl;
  fOutputList = new TList();
  fOutputList->SetOwner();
  
  TH3F *fVertexDist = new TH3F("fVertexDist","Vertex Distribution",20,-1.,1., 20,-1.,1., 600,-30.,30.);
  fVertexDist->GetXaxis()->SetTitle("X Vertex (cm)");
  fVertexDist->GetYaxis()->SetTitle("Y Vertex (cm)");
  fVertexDist->GetZaxis()->SetTitle("Z Vertex (cm)");
  fOutputList->Add(fVertexDist);
  
  
  TH2F *fDCAxyDistPlus = new TH2F("fDCAxyDistPlus","DCA distribution",300,0.,3., 50,0.,5.);
  fOutputList->Add(fDCAxyDistPlus);
  TH2F *fDCAzDistPlus = new TH2F("fDCAzDistPlus","DCA distribution",300,0.,3., 50,0.,5.);
  fOutputList->Add(fDCAzDistPlus);
  TH2F *fDCAxyDistMinus = new TH2F("fDCAxyDistMinus","DCA distribution",300,0.,3., 50,0.,5.);
  fOutputList->Add(fDCAxyDistMinus);
  TH2F *fDCAzDistMinus = new TH2F("fDCAzDistMinus","DCA distribution",300,0.,3., 50,0.,5.);
  fOutputList->Add(fDCAzDistMinus);
  
  
  TH1F *fEvents1 = new TH1F("fEvents1","Events vs. fMbin",fMbins,.5,fMbins+.5);
  fOutputList->Add(fEvents1);
  TH1F *fEvents2 = new TH1F("fEvents2","Events vs. fMbin",fMbins,.5,fMbins+.5);
  fOutputList->Add(fEvents2);
  
  TH1F *fMultDist0 = new TH1F("fMultDist0","Multiplicity Distribution",fMultLimit,-.5,fMultLimit-.5);
  fMultDist0->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist0);

  TH1F *fMultDist1 = new TH1F("fMultDist1","Multiplicity Distribution",fMultLimit,-.5,fMultLimit-.5);
  fMultDist1->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist1);
  
  TH1F *fMultDist2 = new TH1F("fMultDist2","Multiplicity Distribution",fMultLimit,-.5,fMultLimit-.5);
  fMultDist2->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist2);
  
  TH1F *fMultDist3 = new TH1F("fMultDist3","Multiplicity Distribution",fMultLimit,-.5,fMultLimit-.5);
  fMultDist3->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist3);
  
  TH3F *fChPtEtaDist = new TH3F("fChPtEtaDist","fChPtEtaDist",2,-1.1,1.1, 300,0.,3., 28,-1.4,1.4);
  fOutputList->Add(fChPtEtaDist);
  TH3F *fChPhiPtDist = new TH3F("fChPhiPtDist","fChPhiPtDist",2,-1.1,1.1, 120,0.,2*PI, 300,0.,3.);
  fOutputList->Add(fChPhiPtDist);
  
  TH2F *fCentEtaDist = new TH2F("fCentEtaDist","",10,-.5,9.5, 28,-1.4,1.4);
  fOutputList->Add(fCentEtaDist);
  TH2F *fCentPtDist = new TH2F("fCentPtDist","",10,-.5,9.5, 600,0.,3.);
  fOutputList->Add(fCentPtDist);

  TH3F *fTOFResponse = new TH3F("fTOFResponse","TOF relative time",20,0.,100., 200,0.,2., 4000,-20000.,20000.);
  fOutputList->Add(fTOFResponse);
  TH3F *fTPCResponse = new TH3F("fTPCResponse","TPCsignal",20,0.,100., 200,0.,2., 1000,0.,1000.);
  fOutputList->Add(fTPCResponse);
 
  TH1F *fRejectedPairs = new TH1F("fRejectedPairs","",400,0.,2.);
  fOutputList->Add(fRejectedPairs);
  TH1F *fRejectedPairsWeighting = new TH1F("fAcceptedPairsWeighting","",400,0.,2.);
  fOutputList->Add(fRejectedPairsWeighting);
  TH1F *fTotalPairsWeighting = new TH1F("fTotalPairsWeighting","",400,0.,2.);
  fOutputList->Add(fTotalPairsWeighting);
  //
  TH1F *fRejectedPairsMC = new TH1F("fRejectedPairsMC","",400,0.,2.);
  fOutputList->Add(fRejectedPairsMC);
  TH1F *fRejectedPairsWeightingMC = new TH1F("fAcceptedPairsWeightingMC","",400,0.,2.);
  fOutputList->Add(fRejectedPairsWeightingMC);
  TH1F *fTotalPairsWeightingMC = new TH1F("fTotalPairsWeightingMC","",400,0.,2.);
  fOutputList->Add(fTotalPairsWeightingMC);
  
  TH1I *fRejectedEvents = new TH1I("fRejectedEvents","",fMbins,0.5,fMbins+.5);
  fOutputList->Add(fRejectedEvents);
    
  TH3F *fPairsDetaDPhiNum = new TH3F("fPairsDetaDPhiNum","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiNum);
  TH3F *fPairsDetaDPhiDen = new TH3F("fPairsDetaDPhiDen","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiDen);
  TH3F *fPairsShareFracDPhiNum = new TH3F("fPairsShareFracDPhiNum","",10,-.5,9.5, 159,0.,1., 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiNum);
  TH3F *fPairsShareFracDPhiDen = new TH3F("fPairsShareFracDPhiDen","",10,-.5,9.5, 159,0.,1., 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiDen);
  TH3D* fPairsPadRowNum = new TH3D("fPairsPadRowNum","", 20,0.,1., 159,0.,1., 40,0.,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowNum);
  TH3D* fPairsPadRowDen = new TH3D("fPairsPadRowDen","", 20,0.,1., 159,0.,1., 40,0.,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowDen);



  TH2D *fResonanceOSPairs = new TH2D("fResonanceOSPairs","",fMbins,.5,fMbins+.5, 1000,0.,2.);
  if(fMCcase) fOutputList->Add(fResonanceOSPairs);
  TH2D *fAllOSPairs = new TH2D("fAllOSPairs","",fMbins,.5,fMbins+.5, 1000,0.,2.);
  if(fMCcase) fOutputList->Add(fAllOSPairs);
  
  TH3D *fPrimarySCPionPairs = new TH3D("fPrimarySCPionPairs","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  if(fMCcase) fOutputList->Add(fPrimarySCPionPairs);
  TH3D *fAllSCPionPairs = new TH3D("fAllSCPionPairs","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0,0.2);
  if(fMCcase) fOutputList->Add(fAllSCPionPairs);
  TH3D *fPrimaryMCPionPairs = new TH3D("fPrimaryMCPionPairs","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  if(fMCcase) fOutputList->Add(fPrimaryMCPionPairs);
  TH3D *fAllMCPionPairs = new TH3D("fAllMCPionPairs","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  if(fMCcase) fOutputList->Add(fAllMCPionPairs);

  //
  TH1D *fMuonParents = new TH1D("fMuonParents","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fMuonParents);
  TH1D *fSecondaryMuonParents = new TH1D("fSecondaryMuonParents","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fSecondaryMuonParents);
  TH3D *fMuonPionDeltaQinv = new TH3D("fMuonPionDeltaQinv","",2,-0.5,1.5, 20,0.,1., 100,-0.2,0.2);
  if(fMCcase) fOutputList->Add(fMuonPionDeltaQinv);
  TH1D *fPionCandidates = new TH1D("fPionCandidates","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fPionCandidates);
  //  
  

  TProfile *fAvgMult = new TProfile("fAvgMult","",fMbins,.5,fMbins+.5, 0,1500,"");
  fOutputList->Add(fAvgMult);

  TH2D *fTrackChi2NDF = new TH2D("fTrackChi2NDF","",20,0.,100., 100,0.,10.);
  fOutputList->Add(fTrackChi2NDF);
  TH2D *fTrackTPCncls = new TH2D("fTrackTPCncls","",20,0.,100., 110,50.,160.);
  fOutputList->Add(fTrackTPCncls);


  TH1D *fTPNRejects3pion1 = new TH1D("fTPNRejects3pion1","",fQbinsQ3,0.,fQupperBoundQ3);
  fOutputList->Add(fTPNRejects3pion1);
  TH1D *fTPNRejects3pion2 = new TH1D("fTPNRejects3pion2","",fQbinsQ3,0.,fQupperBoundQ3);
  fOutputList->Add(fTPNRejects3pion2);
  TH1D *fTPNRejects4pion1 = new TH1D("fTPNRejects4pion1","",fQbinsQ4,0.,fQupperBoundQ4);
  fOutputList->Add(fTPNRejects4pion1);

  /*TH3D *fKT3DistTerm1 = new TH3D("fKT3DistTerm1","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  TH3D *fKT3DistTerm5 = new TH3D("fKT3DistTerm5","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  fOutputList->Add(fKT3DistTerm1);
  fOutputList->Add(fKT3DistTerm5);
  TH3D *fKT4DistTerm1 = new TH3D("fKT4DistTerm1","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  TH3D *fKT4DistTerm13 = new TH3D("fKT4DistTerm13","",fMbins,.5,fMbins+.5, 20,0.,1., 20,0.,0.2);
  fOutputList->Add(fKT4DistTerm1);
  fOutputList->Add(fKT4DistTerm13);


  TProfile2D *fKT3AvgpT = new TProfile2D("fKT3AvgpT","",fMbins,.5,fMbins+.5, 4,-0.5,3.5, 0.,1.0,"");
  fOutputList->Add(fKT3AvgpT);
  TProfile2D *fKT4AvgpT = new TProfile2D("fKT4AvgpT","",fMbins,.5,fMbins+.5, 4,-0.5,3.5, 0.,1.0,"");
  fOutputList->Add(fKT4AvgpT);
  TH3D* fQ3AvgpTENsum0 = new TH3D("fQ3AvgpTENsum0","", 4,-0.5,3.5, fQbinsQ3,0,fQupperBoundQ3, 180,0.1,1.0);
  fOutputList->Add(fQ3AvgpTENsum0);
  TH3D* fQ3AvgpTENsum3 = new TH3D("fQ3AvgpTENsum3","", 4,-0.5,3.5, fQbinsQ3,0,fQupperBoundQ3, 180,0.1,1.0);
  fOutputList->Add(fQ3AvgpTENsum3);
  TH3D* fQ3AvgpTENsum6 = new TH3D("fQ3AvgpTENsum6","", 4,-0.5,3.5, fQbinsQ3,0,fQupperBoundQ3, 180,0.1,1.0);
  fOutputList->Add(fQ3AvgpTENsum6);
  //
  TH3D* fQ4AvgpTENsum0 = new TH3D("fQ4AvgpTENsum0","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpTENsum0);
  TH3D* fQ4AvgpTENsum1 = new TH3D("fQ4AvgpTENsum1","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpTENsum1);
  TH3D* fQ4AvgpTENsum2 = new TH3D("fQ4AvgpTENsum2","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpTENsum2);
  TH3D* fQ4AvgpTENsum3 = new TH3D("fQ4AvgpTENsum3","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpTENsum3);
  TH3D* fQ4AvgpTENsum6 = new TH3D("fQ4AvgpTENsum6","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpTENsum6);
  //
  TH3D* fQ4AvgkTENsum0 = new TH3D("fQ4AvgkTENsum0","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 100,0,1.0);
  fOutputList->Add(fQ4AvgkTENsum0);
  TH3D* fQ4AvgkTENsum6 = new TH3D("fQ4AvgkTENsum6","", 4,-0.5,3.5, fQbinsQ4,0,fQupperBoundQ4, 100,0,1.0);
  fOutputList->Add(fQ4AvgkTENsum6);*/
  //
  TH1D *fMCWeight3DTerm1SC = new TH1D("fMCWeight3DTerm1SC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm1SCden = new TH1D("fMCWeight3DTerm1SCden","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm2SC = new TH1D("fMCWeight3DTerm2SC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm2SCden = new TH1D("fMCWeight3DTerm2SCden","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm1MC = new TH1D("fMCWeight3DTerm1MC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm1MCden = new TH1D("fMCWeight3DTerm1MCden","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm2MC = new TH1D("fMCWeight3DTerm2MC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm2MCden = new TH1D("fMCWeight3DTerm2MCden","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm3MC = new TH1D("fMCWeight3DTerm3MC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm3MCden = new TH1D("fMCWeight3DTerm3MCden","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm4MC = new TH1D("fMCWeight3DTerm4MC","", 20,0.,0.2);
  TH1D *fMCWeight3DTerm4MCden = new TH1D("fMCWeight3DTerm4MCden","", 20,0.,0.2);
  fOutputList->Add(fMCWeight3DTerm1SC);
  fOutputList->Add(fMCWeight3DTerm1SCden);
  fOutputList->Add(fMCWeight3DTerm2SC);
  fOutputList->Add(fMCWeight3DTerm2SCden);
  fOutputList->Add(fMCWeight3DTerm1MC);
  fOutputList->Add(fMCWeight3DTerm1MCden);
  fOutputList->Add(fMCWeight3DTerm2MC);
  fOutputList->Add(fMCWeight3DTerm2MCden);
  fOutputList->Add(fMCWeight3DTerm3MC);
  fOutputList->Add(fMCWeight3DTerm3MCden);
  fOutputList->Add(fMCWeight3DTerm4MC);
  fOutputList->Add(fMCWeight3DTerm4MCden);


      
  for(Int_t mb=0; mb<fMbins; mb++){
    if(fCollisionType==0) {if((mb < fCentBinLowLimit) || (mb > fCentBinHighLimit)) continue;}
    
    for(Int_t edB=0; edB<fEDbins; edB++){
      
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t term=0; term<2; term++){
	    
	    TString *nameEx2 = new TString("TwoParticle_Charge1_");
	    *nameEx2 += c1;
	    nameEx2->Append("_Charge2_");
	    *nameEx2 += c2;
	    nameEx2->Append("_M_");
	    *nameEx2 += mb;
	    nameEx2->Append("_ED_");
	    *nameEx2 += edB;
	    nameEx2->Append("_Term_");
	    *nameEx2 += term+1;
	    
	    if( (c1+c2)==1 ) {if(c1!=0) continue;}// skip degenerate histogram
	    
	    TString *nameUnitMult=new TString(nameEx2->Data());
	    nameUnitMult->Append("_UnitMult");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fUnitMultBin = new TH2D(nameUnitMult->Data(),"Two Particle Distribution",21,0.5,21.5, fQbinsQ2,0.,fQupperBoundQ2);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fUnitMultBin);

	    if(edB==0){
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2 = new TH2D(nameEx2->Data(),"Two Particle Distribution",20,0.,1., fQbinsQ2,0.,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2);
	      //
	      if(c1==c2 && term==1){
		TString *nameBuild=new TString(nameEx2->Data());
		nameBuild->Append("_Build");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fBuild = new TH3D(nameBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, 20,0.,1., fQbinsQ2,0.,fQupperBoundQ2);
		fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fBuild);
	      }
	      TString *nameEx2QW=new TString(nameEx2->Data());
	      nameEx2QW->Append("_QW");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2QW = new TH2D(nameEx2QW->Data(),"Two Particle Distribution",20,0.,1., fQbinsQ2,0.,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2QW);
	      TString *nameAvgP=new TString(nameEx2->Data());
	      nameAvgP->Append("_AvgP");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fAvgP = new TProfile2D(nameAvgP->Data(),"",10,0.,1, fQbinsQ2,0.,fQupperBoundQ2, 0.,1.0,"");
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fAvgP);
	      
	    
	      if(fMCcase){
		// Momentum resolution histos
		TString *nameIdeal = new TString(nameEx2->Data());
		nameIdeal->Append("_Ideal");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal = new TH2D(nameIdeal->Data(),"Two Particle Distribution",11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal);
		TString *nameSmeared = new TString(nameEx2->Data());
		nameSmeared->Append("_Smeared");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared = new TH2D(nameSmeared->Data(),"Two Particle Distribution",11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared);
		//
		// Muon correction histos
		TString *nameMuonIdeal=new TString(nameEx2->Data());
		nameMuonIdeal->Append("_MuonIdeal");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonIdeal = new TH2D(nameMuonIdeal->Data(),"", 11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonIdeal);
		TString *nameMuonSmeared=new TString(nameEx2->Data());
		nameMuonSmeared->Append("_MuonSmeared");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonSmeared = new TH2D(nameMuonSmeared->Data(),"", 11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonSmeared);
		//
		TString *nameMuonPionK2=new TString(nameEx2->Data());
		nameMuonPionK2->Append("_MuonPionK2");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonPionK2 = new TH2D(nameMuonPionK2->Data(),"", 11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonPionK2);
		//
		TString *namePionPionK2=new TString(nameEx2->Data());
		namePionPionK2->Append("_PionPionK2");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPionPionK2 = new TH2D(namePionPionK2->Data(),"", 11,0.5,11.5, fQbinsQ2,0.,fQupperBoundQ2);
		if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPionPionK2);
		//
		//
		TString *nameEx2MC=new TString(nameEx2->Data());
		nameEx2MC->Append("_MCqinv");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv = new TH1D(nameEx2MC->Data(),"", fQbinsQ2,0.,fQupperBoundQ2);
		fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv);
		TString *nameEx2MCQW=new TString(nameEx2->Data());
		nameEx2MCQW->Append("_MCqinvQW");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = new TH1D(nameEx2MCQW->Data(),"", fQbinsQ2,0.,fQupperBoundQ2);
		fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW);
		//
		TString *nameEx2PIDpurityDen=new TString(nameEx2->Data());
		nameEx2PIDpurityDen->Append("_PIDpurityDen");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = new TH2D(nameEx2PIDpurityDen->Data(),"Two Particle Distribution",20,0.,1, fQbinsQ2,0.,fQupperBoundQ2);
		fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen);
		TString *nameEx2PIDpurityNum=new TString(nameEx2->Data());
		nameEx2PIDpurityNum->Append("_PIDpurityNum");
		Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = new TH3D(nameEx2PIDpurityNum->Data(),"Two Particle Distribution",16,0.5,16.5, 20,0.,1, fQbinsQ2,0.,fQupperBoundQ2);
		fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum);
	      }
	      TString *nameEx2OSLB1 = new TString(nameEx2->Data()); 
	      nameEx2OSLB1->Append("_osl_b1");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL);
	      nameEx2OSLB1->Append("_QW");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW);
	      //
	      TString *nameEx2OSLB2 = new TString(nameEx2->Data()); 
	      nameEx2OSLB2->Append("_osl_b2");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL);
	      nameEx2OSLB2->Append("_QW");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW);
	    }
	    
	  }// term_2
	  
	  
	  
	  // skip 3-particle if Tabulate6DPairs is true
	  if(fTabulatePairs) continue;
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t term=0; term<5; term++){
	      
	      TString *namePC3 = new TString("ThreeParticle_Charge1_");
	      *namePC3 += c1;
	      namePC3->Append("_Charge2_");
	      *namePC3 += c2;
	      namePC3->Append("_Charge3_");
	      *namePC3 += c3;
	      namePC3->Append("_M_");
	      *namePC3 += mb;
	      namePC3->Append("_ED_");
	      *namePC3 += edB;
	      namePC3->Append("_Term_");
	      *namePC3 += term+1;
	      
	      ///////////////////////////////////////
	      // skip degenerate histograms
	      if( (c1+c2+c3)==1) {if(c3!=1) continue;}
	      if( (c1+c2+c3)==2) {if(c1!=0) continue;}
	      /////////////////////////////////////////
	      
	      
	      TString *nameNorm=new TString(namePC3->Data());
	      nameNorm->Append("_Norm");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3 = new TH1D(nameNorm->Data(),"Norm",1,-0.5,0.5);
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fNorm3);
	      //
	      
	      TString *name1DQ=new TString(namePC3->Data());
	      name1DQ->Append("_1D");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3 = new TH1D(name1DQ->Data(),"", fQbinsQ3,0.,fQupperBoundQ3);
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3);
	      if(c1==0 && c2==0 && c3==0){
		TString *name3DQ=new TString(namePC3->Data());
		name3DQ->Append("_3D");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms33D = new TH3D(name3DQ->Data(),"", fQbinsQinv3D,0.,fQupperBoundQinv3D, fQbinsQinv3D,0.,fQupperBoundQinv3D, fQbinsQinv3D,0.,fQupperBoundQinv3D);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms33D);
	      }
	      //
	      TString *nameKfactor=new TString(namePC3->Data());
	      nameKfactor->Append("_Kfactor");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor = new TProfile(nameKfactor->Data(),"", fQbinsQ3,0.,fQupperBoundQ3, 0.,100., "");
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor);
	      if(c1==0 && c2==0 && c3==0){
		TString *nameKfactor3D=new TString(namePC3->Data());
		nameKfactor3D->Append("_Kfactor3D");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor3D = new TProfile3D(nameKfactor3D->Data(),"", fQbinsQinv3D,0.,fQupperBoundQinv3D, fQbinsQinv3D,0.,fQupperBoundQinv3D, fQbinsQinv3D,0.,fQupperBoundQinv3D, "");
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor3D);
	      }
	      //
	      //TString *nameKfactorW=new TString(namePC3->Data());
	      //nameKfactorW->Append("_KfactorWeighted");
	      //Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted = new TProfile(nameKfactorW->Data(),"", fQbinsQ3,0.,fQupperBoundQ3, 0.,100., "");
	      //fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted);
	      //
	      //TString *nameMeanQinv=new TString(namePC3->Data());
	      //nameMeanQinv->Append("_MeanQinv");
	      //Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMeanQinv = new TProfile(nameMeanQinv->Data(),"", fQbinsQ3,0.,fQupperBoundQ3, 0.,.2, "");
	      //fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMeanQinv);
	      
	      if(fMCcase==kTRUE){
		// Momentum resolution correction histos
		TString *nameMomResIdeal=new TString(namePC3->Data());
		nameMomResIdeal->Append("_Ideal");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fIdeal = new TH2D(nameMomResIdeal->Data(),"", 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fIdeal);
		TString *nameMomResSmeared=new TString(namePC3->Data());
		nameMomResSmeared->Append("_Smeared");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fSmeared = new TH2D(nameMomResSmeared->Data(),"", 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fSmeared);
		// Muon correction histos
		TString *nameMuonIdeal=new TString(namePC3->Data());
		nameMuonIdeal->Append("_MuonIdeal");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonIdeal = new TH3D(nameMuonIdeal->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonIdeal);
		TString *nameMuonSmeared=new TString(namePC3->Data());
		nameMuonSmeared->Append("_MuonSmeared");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonSmeared = new TH3D(nameMuonSmeared->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonSmeared);
		//
		TString *nameMuonPionK3=new TString(namePC3->Data());
		nameMuonPionK3->Append("_MuonPionK3");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonPionK3 = new TH3D(nameMuonPionK3->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonPionK3);
		//
		TString *namePionPionK3=new TString(namePC3->Data());
		namePionPionK3->Append("_PionPionK3");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fPionPionK3 = new TH3D(namePionPionK3->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0.,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fPionPionK3);
		
	      }// MCcase
	      //
	      if(c1==c2 && c1==c3 && term==4 ){
		TString *nameBuild=new TString(namePC3->Data());
		nameBuild->Append("_Build");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild = new TH2D(nameBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0.,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuild);
		//
		TString *nameCumulantBuild=new TString(namePC3->Data());
		nameCumulantBuild->Append("_CumulantBuild");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fCumulantBuild = new TH2D(nameCumulantBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0.,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fCumulantBuild);
		//
		TString *nameBuildNeg=new TString(namePC3->Data());
		nameBuildNeg->Append("_BuildNeg");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuildNeg = new TH2D(nameBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0.,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuildNeg);
		//
		TString *nameCumulantBuildNeg=new TString(namePC3->Data());
		nameCumulantBuildNeg->Append("_CumulantBuildNeg");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fCumulantBuildNeg = new TH2D(nameCumulantBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0.,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fCumulantBuildNeg);
		//
		TString *nameBuildErr=new TString(namePC3->Data());
		nameBuildErr->Append("_BuildErr");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuildErr = new TH2D(nameBuildErr->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0.,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fBuildErr);
	      }// term=4
	      
	    }// term_3
	    
	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		TString *namePC4 = new TString("FourParticle_Charge1_");
		*namePC4 += c1;
		namePC4->Append("_Charge2_");
		*namePC4 += c2;
		namePC4->Append("_Charge3_");
		*namePC4 += c3;
		namePC4->Append("_Charge4_");
		*namePC4 += c4;
		namePC4->Append("_M_");
		*namePC4 += mb;
		namePC4->Append("_ED_");
		*namePC4 += edB;
		namePC4->Append("_Term_");
		*namePC4 += term+1;
		
		///////////////////////////////////////
		// skip degenerate histograms
		if( (c1+c2+c3+c4)==1) {if(c4!=1) continue;}
		if( (c1+c2+c3+c4)==2) {if(c3+c4!=2) continue;}
		if( (c1+c2+c3+c4)==3) {if(c1!=0) continue;}
		/////////////////////////////////////////
		
		TString *nameNorm=new TString(namePC4->Data());
		nameNorm->Append("_Norm");
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4 = new TH1D(nameNorm->Data(),"Norm",1,-0.5,0.5);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4);
		//
		TString *name1DQ=new TString(namePC4->Data());
		name1DQ->Append("_1D");
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = new TH1D(name1DQ->Data(),"", fQbinsQ4,0.,fQupperBoundQ4);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4);
		//
		TString *nameKfactor=new TString(namePC4->Data());
		nameKfactor->Append("_Kfactor");
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = new TProfile(nameKfactor->Data(),"", fQbinsQ4,0.,fQupperBoundQ4, 0.,100., "");
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor);
		//
		//TString *nameKfactorW=new TString(namePC4->Data());
		//nameKfactorW->Append("_KfactorWeighted");
		//Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted = new TProfile(nameKfactorW->Data(),"", fQbinsQ4,0.,fQupperBoundQ4, 0.,100., "");
		//fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted);
		//
		if(c1==c2 && c1==c3 && c1==c4 && term==12 ){
		  TString *nameBuild=new TString(namePC4->Data());
		  nameBuild->Append("_Build");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild = new TH2D(nameBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuild);
		  //
		  TString *namePrimeBuild=new TString(namePC4->Data());
		  namePrimeBuild->Append("_PrimeBuild");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuild = new TH2D(namePrimeBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuild);
		  //
		  TString *namePrimePrimeBuild=new TString(namePC4->Data());
		  namePrimePrimeBuild->Append("_PrimePrimeBuild");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuild = new TH2D(namePrimePrimeBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuild);
		  //
		  TString *nameCumulantBuild=new TString(namePC4->Data());
		  nameCumulantBuild->Append("_CumulantBuild");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuild = new TH2D(nameCumulantBuild->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuild);
		  //
		  //
		  TString *nameBuildNeg=new TString(namePC4->Data());
		  nameBuildNeg->Append("_BuildNeg");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildNeg = new TH2D(nameBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildNeg);
		  //
		  TString *namePrimeBuildNeg=new TString(namePC4->Data());
		  namePrimeBuildNeg->Append("_PrimeBuildNeg");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildNeg = new TH2D(namePrimeBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildNeg);
		  //
		  TString *namePrimePrimeBuildNeg=new TString(namePC4->Data());
		  namePrimePrimeBuildNeg->Append("_PrimePrimeBuildNeg");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildNeg = new TH2D(namePrimePrimeBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildNeg);
		  //
		  TString *nameCumulantBuildNeg=new TString(namePC4->Data());
		  nameCumulantBuildNeg->Append("_CumulantBuildNeg");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildNeg = new TH2D(nameCumulantBuildNeg->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildNeg);
		  //
		  //
		  TString *nameBuildErr=new TString(namePC4->Data());
		  nameBuildErr->Append("_BuildErr");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildErr = new TH2D(nameBuildErr->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildErr);
		  //
		  if(c1==0 && c2==0 && c3==0 && c4==0){
		    TString *nameBuildFromFits=new TString(namePC4->Data());
		    nameBuildFromFits->Append("_BuildFromFits");
		    Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildFromFits = new TH3D(nameBuildFromFits->Data(),"", 2,0.5,2.5, kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fBuildFromFits);
		    //
		    TString *namePrimeBuildFromFits=new TString(namePC4->Data());
		    namePrimeBuildFromFits->Append("_PrimeBuildFromFits");
		    Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildFromFits = new TH3D(namePrimeBuildFromFits->Data(),"", 2,0.5,2.5, kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimeBuildFromFits);
		    //
		    TString *namePrimePrimeBuildFromFits=new TString(namePC4->Data());
		    namePrimePrimeBuildFromFits->Append("_PrimePrimeBuildFromFits");
		    Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildFromFits = new TH3D(namePrimePrimeBuildFromFits->Data(),"", 2,0.5,2.5, kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPrimePrimeBuildFromFits);
		    //
		    TString *nameCumulantBuildFromFits=new TString(namePC4->Data());
		    nameCumulantBuildFromFits->Append("_CumulantBuildFromFits");
		    Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildFromFits = new TH3D(nameCumulantBuildFromFits->Data(),"", 2,0.5,2.5, kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0.,fQupperBoundQ4);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fCumulantBuildFromFits);
		  }
		}
		
		if(fMCcase==kTRUE){
		  // Momentum resolution correction histos
		  TString *nameMomResIdeal=new TString(namePC4->Data());
		  nameMomResIdeal->Append("_Ideal");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fIdeal = new TH2D(nameMomResIdeal->Data(),"", 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fIdeal);
		  TString *nameMomResSmeared=new TString(namePC4->Data());
		  nameMomResSmeared->Append("_Smeared");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fSmeared = new TH2D(nameMomResSmeared->Data(),"", 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fSmeared);
		  // Muon correction histos
		  TString *nameMuonIdeal=new TString(namePC4->Data());
		  nameMuonIdeal->Append("_MuonIdeal");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonIdeal = new TH3D(nameMuonIdeal->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonIdeal);
		  TString *nameMuonSmeared=new TString(namePC4->Data());
		  nameMuonSmeared->Append("_MuonSmeared");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonSmeared = new TH3D(nameMuonSmeared->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonSmeared);
		  //
		  TString *nameMuonPionK4=new TString(namePC4->Data());
		  nameMuonPionK4->Append("_MuonPionK4");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonPionK4 = new TH3D(nameMuonPionK4->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonPionK4);
		  //
		  TString *namePionPionK4=new TString(namePC4->Data());
		  namePionPionK4->Append("_PionPionK4");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPionPionK4 = new TH3D(namePionPionK4->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0.,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPionPionK4);
		  
		}// MCcase
		
		
	      }
	    }
	    
	  }//c3
	}//c2
      }//c1
    }// ED
  }// mbin


  
  if(fTabulatePairs){
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	for(Int_t mb=0; mb<fMbins; mb++){
	  for(Int_t edB=0; edB<fEDbins; edB++){
	    
	    TString *nameNum = new TString("TPN_num_Kt_");
	    *nameNum += tKbin;
	    nameNum->Append("_Ky_");
	    *nameNum += yKbin;
	    nameNum->Append("_M_");
	    *nameNum += mb;
	    nameNum->Append("_ED_");
	    *nameNum += edB;
	    
	    TString *nameDen = new TString("TPN_den_Kt_");
	    *nameDen += tKbin;
	    nameDen->Append("_Ky_");
	    *nameDen += yKbin;
	    nameDen->Append("_M_");
	    *nameDen += mb;
	    nameDen->Append("_ED_");
	    *nameDen += edB;
	    
	    if(edB<=1){
	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD = new TH3D(nameNum->Data(),"", kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD);
	      
	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD = new TH3D(nameDen->Data(),"", kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights, kQbinsWeights,0.,fQupperBoundWeights);
	      fOutputList->Add(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD);
	    }	
	  
	  }
	}
      }
    }
    
  }
 
  
  TProfile *fQsmearMean = new TProfile("fQsmearMean","",2,0.5,2.5, -0.2,0.2,"");
  fOutputList->Add(fQsmearMean);
  TProfile *fQsmearSq = new TProfile("fQsmearSq","",2,0.5,2.5, -2,2,"");
  fOutputList->Add(fQsmearSq);
  TH2D *fQ2Res = new TH2D("fQ2Res","",20,0.,1, 200,-.2,.2);
  fOutputList->Add(fQ2Res);
  TH2D *fQ3Res = new TH2D("fQ3Res","",20,0.,1, 200,-.3,.3);
  fOutputList->Add(fQ3Res);
  TH2D *fQ4Res = new TH2D("fQ4Res","",20,0.,1, 200,-.4,.4);
  fOutputList->Add(fQ4Res);
  
  TH2D *DistQinv4pion = new TH2D("DistQinv4pion","",6,0.5,6.5, fQbinsQ2,0.,fQupperBoundQ2);
  fOutputList->Add(DistQinv4pion);
  TH2D *DistQinvMC4pion = new TH2D("DistQinvMC4pion","",6,0.5,6.5, fQbinsQ2,0.,fQupperBoundQ2);
  if(fMCcase) fOutputList->Add(DistQinvMC4pion);

  TH2D *fAvgQ12VersusQ3 = new TH2D("fAvgQ12VersusQ3","",40,0.,0.2, fQbinsQ3,0.,fQupperBoundQ3);
  fOutputList->Add(fAvgQ12VersusQ3);
  TH2D *fAvgQ13VersusQ3 = new TH2D("fAvgQ13VersusQ3","",40,0.,0.2, fQbinsQ3,0.,fQupperBoundQ3);
  fOutputList->Add(fAvgQ13VersusQ3);
  TH2D *fAvgQ23VersusQ3 = new TH2D("fAvgQ23VersusQ3","",40,0.,0.2, fQbinsQ3,0.,fQupperBoundQ3);
  fOutputList->Add(fAvgQ23VersusQ3);

  TH1D *fDistPionParents4 = new TH1D("fDistPionParents4","",4,0.5,4.5);
  fOutputList->Add(fDistPionParents4);

  TH2D *fDistTPCNclsFindable = new TH2D("fDistTPCNclsFindable","", 100,0.,0.5, 201,-0.5,200.5);
  fDistTPCNclsFindable->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsFindable->GetYaxis()->SetTitle("Ncls Findable");
  fOutputList->Add(fDistTPCNclsFindable);
  TProfile *fProfileTPCNclsFindable = new TProfile("fProfileTPCNclsFindable","",100,0.,0.5, 0.,200., "");
  fProfileTPCNclsFindable->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsFindable->GetYaxis()->SetTitle("<Ncls Findable>");
  fOutputList->Add(fProfileTPCNclsFindable);
  //
  TH2D *fDistTPCNclsCrossed = new TH2D("fDistTPCNclsCrossed","",100,0.,0.5, 201,-0.5,200.5);
  fDistTPCNclsCrossed->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsCrossed->GetYaxis()->SetTitle("Ncls Crossed");
  fOutputList->Add(fDistTPCNclsCrossed);
  TProfile *fProfileTPCNclsCrossed = new TProfile("fProfileTPCNclsCrossed","",100,0.,0.5, 0.,200., "");
  fProfileTPCNclsCrossed->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsCrossed->GetYaxis()->SetTitle("<Ncls Crossed>");
  fOutputList->Add(fProfileTPCNclsCrossed);
  //
  TH2D *fDistTPCNclsFindableRatio = new TH2D("fDistTPCNclsFindableRatio","",100,0.,0.5, 100,0.,1.);
  fDistTPCNclsFindableRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsFindableRatio->GetYaxis()->SetTitle("Ncls / Ncls Findable");
  fOutputList->Add(fDistTPCNclsFindableRatio);
  TProfile *fProfileTPCNclsFindableRatio = new TProfile("fProfileTPCNclsFindableRatio","",100,0.,0.5, 0.,1., "");
  fProfileTPCNclsFindableRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsFindableRatio->GetYaxis()->SetTitle("<Ncls / Ncls Findable>");
  fOutputList->Add(fProfileTPCNclsFindableRatio);
  //
  TH2D *fDistTPCNclsCrossedRatio = new TH2D("fDistTPCNclsCrossedRatio","",100,0.,0.5, 100,0.,1.);
  fDistTPCNclsCrossedRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsCrossedRatio->GetYaxis()->SetTitle("Ncls / Ncls Crossed");
  fOutputList->Add(fDistTPCNclsCrossedRatio);
  TProfile *fProfileTPCNclsCrossedRatio = new TProfile("fProfileTPCNclsCrossedRatio","",100,0.,0.5, 0.,1., "");
  fProfileTPCNclsCrossedRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsCrossedRatio->GetYaxis()->SetTitle("<Ncls / Ncls Crossed>");
  fOutputList->Add(fProfileTPCNclsCrossedRatio);

  TH2D *fc4QSFitNum = new TH2D("fc4QSFitNum","",7,0.5,7.5, fQbinsQ4,0.,fQupperBoundQ4);
  fOutputList->Add(fc4QSFitNum);
  TH2D *fc4QSFitDen = new TH2D("fc4QSFitDen","",7,0.5,7.5, fQbinsQ4,0.,fQupperBoundQ4);
  fOutputList->Add(fc4QSFitDen);

  TH3F *fq1Dist = new TH3F("fq1Dist","",fMbins,.5,fMbins+.5, 4,0.5,4.5, 200,0,10);// Mult, pT bin, q2 bin
  fOutputList->Add(fq1Dist);
  TH3F *fq2Dist = new TH3F("fq2Dist","",fMbins,.5,fMbins+.5, 4,0.5,4.5, 200,0,10);// Mult, pT bin, q2 bin
  fOutputList->Add(fq2Dist);

  TH2D *fLowPtDist = new TH2D("fLowPtDist","",fMbins,.5,fMbins+.5, 500,0.5,500.5);
  fOutputList->Add(fLowPtDist);
  /*TH3D *fPtMultCrossing1 = new TH3D("fPtMultCrossing1","",fMbins,.5,fMbins+.5, 500,0.5,500.5, 500,0.5,500.5);
  TH3D *fPtMultCrossing2 = new TH3D("fPtMultCrossing2","",fMbins,.5,fMbins+.5, 500,0.5,500.5, 500,0.5,500.5);
  TH3D *fPtMultCrossing3 = new TH3D("fPtMultCrossing3","",fMbins,.5,fMbins+.5, 500,0.5,500.5, 500,0.5,500.5);
  fOutputList->Add(fPtMultCrossing1);
  fOutputList->Add(fPtMultCrossing2);
  fOutputList->Add(fPtMultCrossing3);*/

  TH2D *fcumulant2INT = new TH2D("fcumulant2INT","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fcumulant2DIFF = new TH2D("fcumulant2DIFF","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fcumulant2EN = new TH2D("fcumulant2EN","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH1D *fcumulant2TotalINT = new TH1D("fcumulant2TotalINT","", fMbins,.5,fMbins+.5);
  TH2D *fSingleSumCosINT = new TH2D("fSingleSumCosINT","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fSingleSumSinINT = new TH2D("fSingleSumSinINT","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fSingleSumENINT = new TH2D("fSingleSumENINT","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fSingleSumCosDIFF = new TH2D("fSingleSumCosDIFF","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fSingleSumSinDIFF = new TH2D("fSingleSumSinDIFF","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH2D *fSingleSumENDIFF = new TH2D("fSingleSumENDIFF","", fMbins,.5,fMbins+.5, 40,0.5,40.5);
  TH1D *fSingleSumCosTotalINT = new TH1D("fSingleSumCosTotalINT","", fMbins,.5,fMbins+.5);
  TH1D *fSingleSumSinTotalINT = new TH1D("fSingleSumSinTotalINT","", fMbins,.5,fMbins+.5);
  TH1D *fSingleSumTotalEN = new TH1D("fSingleSumTotalEN","", fMbins,.5,fMbins+.5);
  fOutputList->Add(fcumulant2INT);
  fOutputList->Add(fcumulant2DIFF);
  fOutputList->Add(fcumulant2EN);
  fOutputList->Add(fcumulant2TotalINT);
  fOutputList->Add(fSingleSumCosINT);
  fOutputList->Add(fSingleSumSinINT);
  fOutputList->Add(fSingleSumENINT);
  fOutputList->Add(fSingleSumCosDIFF);
  fOutputList->Add(fSingleSumSinDIFF);
  fOutputList->Add(fSingleSumENDIFF);
  fOutputList->Add(fSingleSumCosTotalINT);
  fOutputList->Add(fSingleSumSinTotalINT);
  fOutputList->Add(fSingleSumTotalEN);

  TH2D *fQoutSum = new TH2D("fQoutSum","",fQbinsQ4,0.,fQupperBoundQ4, 100, 0,1.0);
  TH2D *fQsideSum = new TH2D("fQsideSum","",fQbinsQ4,0.,fQupperBoundQ4, 100, 0,1.0);
  TH2D *fQlongSum = new TH2D("fQlongSum","",fQbinsQ4,0.,fQupperBoundQ4, 100, 0,1.0);
  fOutputList->Add(fQoutSum);
  fOutputList->Add(fQsideSum);
  fOutputList->Add(fQlongSum);
  
  TH2D *fkTDist = new TH2D("fkTDist","", fMbins,.5,fMbins+.5, 100,0.,1);
  fOutputList->Add(fkTDist);

  
  cout<<"End histogram creation"<<endl;
  ////////////////////////////////////
  ///////////////////////////////////  
  
  PostData(1, fOutputList);
  
}

//________________________________________________________________________
void AliFourPion::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  //cout<<"===========  Event # "<<fEventCounter+1<<"  ==========="<<endl;
  fEventCounter++;
  if(fEventCounter%1000==0) cout<<"===========  Event # "<<fEventCounter<<"  ==========="<<endl;
  
  if(!fAODcase) {cout<<"ESDs not supported"<<endl; return;}
  
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}
  
     
  // Trigger Cut
  if(fAOD->GetRunNumber() >= 136851 && fAOD->GetRunNumber() <= 139517){// 10h data
    Bool_t isSelected1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(!isSelected1 && !fMCcase) {return;}
  }else if(fAOD->GetRunNumber() >= 167693 && fAOD->GetRunNumber() <= 170593){// 11h data
    Bool_t isSelected1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
    Bool_t isSelected2 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
    if(!isSelected1 && !isSelected2 && !fMCcase) {return;}
  }else {
    Bool_t isSelected[4]={kFALSE};
    isSelected[0] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    isSelected[1] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
    isSelected[2] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    isSelected[3] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kHighMult);
    if(!isSelected[fTriggerType] && !fMCcase) return;
  }
  
  
  AliCentrality *centrality;// for AODs and ESDs
  

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  
  TClonesArray *mcArray = 0x0;
  if(fMCcase){
    if(fAODcase){ 
      mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
      if(!mcArray || mcArray->GetEntriesFast() >= kMCarrayLimit){
	cout<<"No MC particle branch found or Array too large!!"<<endl;
	return;
      }
    }
  }
  
  
  UInt_t status=0;
  Int_t positiveTracks=0, negativeTracks=0;
  Int_t myTracks=0, pionCount=0, kaonCount=0, protonCount=0;
   
  Double_t vertex[3]={0};
  Int_t mixingEDbin=0;
  Double_t zstep=2*10/Double_t(fEventMixingEDbins), zstart=-10.;
  /////////////////////////////////////////////////
  
  Float_t cosDIFF[40]={0};
  Float_t sinDIFF[40]={0};
  Float_t countDIFF[40]={0};
  Float_t cosINT[40]={0};
  Float_t sinINT[40]={0};
  Float_t countINT[40]={0};
  Float_t cosTotalINT=0, sinTotalINT=0, countTotalINT=0;
  
  Float_t centralityPercentile=0;
  Float_t cStep=5.0, cStepMixing=5.0, cStart=0;
  Int_t MbinMixing=0;
  
  if(fAODcase){// AOD case
    if(fCollisionType==0){
      centrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
      centralityPercentile = centrality->GetCentralityPercentile("V0M");
      if(centralityPercentile == 0) {cout<<"Centrality = 0, skipping event"<<endl; return;}
      if((centralityPercentile < 5*fCentBinLowLimit) || (centralityPercentile>= 5*(fCentBinHighLimit+1))) {return;}
      cout<<"Centrality % = "<<centralityPercentile<<endl;
    }
    const AliAODVertex* primaryVertexAOD = (AliAODVertex*)fAOD->GetPrimaryVertex();
    ((TH1F*)fOutputList->FindObject("fMultDist0"))->Fill(fAOD->GetNumberOfTracks());
   
    // Pile-up rejection
    AliAnalysisUtils *AnaUtil=new AliAnalysisUtils();
    if(fCollisionType!=0) AnaUtil->SetUseMVPlpSelection(kTRUE);// use Multi-Vertex tool for pp and pPb
    else AnaUtil->SetUseMVPlpSelection(kFALSE);
    Bool_t pileUpCase=AnaUtil->IsPileUpEvent(fAOD); 
    if(pileUpCase) return;
   
    ////////////////////////////////
    // Vertexing
    ((TH1F*)fOutputList->FindObject("fMultDist1"))->Fill(fAOD->GetNumberOfTracks());
    vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
    
    if(fabs(vertex[2]) > 10) {return;} // Z-Vertex Cut 
    ((TH3F*)fOutputList->FindObject("fVertexDist"))->Fill(vertex[0], vertex[1], vertex[2]);
    
    if(!fMCcase && primaryVertexAOD->GetNContributors() < 1) {return;}
   
    ((TH1F*)fOutputList->FindObject("fMultDist2"))->Fill(fAOD->GetNumberOfTracks());
 
    fBfield = fAOD->GetMagneticField();
    
    for(Int_t i=0; i<fEventMixingEDbins; i++){
      if( (vertex[2] >= zstart+i*zstep) && (vertex[2] < zstart+(i+1)*zstep) ){
	mixingEDbin=i;
	break;
      }
    }
    
    
    /////////////////////////////
    // Create Shuffled index list
    Int_t randomIndex[fAOD->GetNumberOfTracks()];
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) randomIndex[i]=i;
    Shuffle(randomIndex,0,fAOD->GetNumberOfTracks()-1);
    /////////////////////////////
    
   
    // Track loop
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* aodtrack = (AliAODTrack*)fAOD->GetTrack(randomIndex[i]);
      if (!aodtrack) continue;
      if(myTracks >= fMultLimit) {cout<<"More tracks than Track Limit"<<endl; return;}
    
      status=aodtrack->GetStatus();
      
      if(!aodtrack->TestFilterBit(BIT(fFilterBit))) continue;// AOD filterBit cut
      ((TH2D*)fOutputList->FindObject("fTrackChi2NDF"))->Fill(centralityPercentile, aodtrack->Chi2perNDF());
      ((TH2D*)fOutputList->FindObject("fTrackTPCncls"))->Fill(centralityPercentile, aodtrack->GetTPCncls());
      if(aodtrack->GetTPCNcls() < fMinTPCncls) continue;// TPC nCluster cut
      if(aodtrack->Chi2perNDF() > fMaxChi2NDF) continue;

      if(fFilterBit != 7){
	Bool_t goodTrackOtherFB = kFALSE;
	for (Int_t j = 0; j < fAOD->GetNumberOfTracks(); j++) {
	  AliAODTrack* aodtrack2 = (AliAODTrack*)fAOD->GetTrack(randomIndex[j]);
	  if(!aodtrack2) continue;
	  if(!aodtrack2->TestFilterBit(BIT(fFilterBit))) continue;
	  
	  if(-(aodtrack->GetID()+1)==aodtrack2->GetID()) {goodTrackOtherFB=kTRUE; break;}
	  
	}
	if(!goodTrackOtherFB) continue;
      }
      
      
      if(aodtrack->Pt() < 0.16) continue;
      if(fabs(aodtrack->Eta()) > 0.8) continue;
     
      Bool_t goodMomentum = aodtrack->GetPxPyPz( fTempStruct[myTracks].fP);
      if(!goodMomentum) continue; 
      aodtrack->GetXYZ( fTempStruct[myTracks].fX);
      
      
      Double_t dca2[2]={0};
      dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - vertex[0],2) + pow(fTempStruct[myTracks].fX[1] - vertex[1],2));
      dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - vertex[2],2));
      Double_t dca3d = sqrt( pow(dca2[0],2) + pow(dca2[1],2));
             
      fTempStruct[myTracks].fStatus = status;
      fTempStruct[myTracks].fFiltermap = aodtrack->GetFilterMap();
      fTempStruct[myTracks].fId = aodtrack->GetID();
      //
      fTempStruct[myTracks].fLabel = aodtrack->GetLabel();
      fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
      if(fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2*PI;
      fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0],2) + pow(fTempStruct[myTracks].fP[1],2));
      fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt,2) + pow(fTempStruct[myTracks].fP[2],2) );
      fTempStruct[myTracks].fEta = aodtrack->Eta();
      fTempStruct[myTracks].fCharge = aodtrack->Charge();
      fTempStruct[myTracks].fDCAXY = dca2[0];
      fTempStruct[myTracks].fDCAZ = dca2[1];
      fTempStruct[myTracks].fDCA = dca3d;
      fTempStruct[myTracks].fClusterMap = aodtrack->GetTPCClusterMap();
      fTempStruct[myTracks].fSharedMap = aodtrack->GetTPCSharedMap();
      

      ////////////////////////////////////////////////////    
      // c2 and d2 calculations
      Int_t BOI= int(fTempStruct[myTracks].fPt / 0.1);
      if(BOI<40){
	cosDIFF[BOI] += cos(2*fTempStruct[myTracks].fPhi);
	sinDIFF[BOI] += sin(2*fTempStruct[myTracks].fPhi);
	countDIFF[BOI]++;
	for(Int_t boi=0; boi<40; boi++){
	  if(boi==BOI) continue;
	  cosINT[boi] += cos(2*fTempStruct[myTracks].fPhi);
	  sinINT[boi] += sin(2*fTempStruct[myTracks].fPhi);
	  countINT[boi]++;
	}
      }
      if(fTempStruct[myTracks].fPt > 0.2){
	cosTotalINT += cos(2*fTempStruct[myTracks].fPhi);
	sinTotalINT += sin(2*fTempStruct[myTracks].fPhi);
	countTotalINT++;
      }
      ///////////////////////////////////////////////////

      if(fTempStruct[myTracks].fMom > 0.9999) continue;// upper P bound
            
     
     
      // PID section
      fTempStruct[myTracks].fElectron = kFALSE;
      fTempStruct[myTracks].fPion = kFALSE;
      fTempStruct[myTracks].fKaon = kFALSE;
      fTempStruct[myTracks].fProton = kFALSE;
      
      Float_t nSigmaTPC[5];
      Float_t nSigmaTOF[5];
      nSigmaTPC[0]=10; nSigmaTPC[1]=10; nSigmaTPC[2]=10; nSigmaTPC[3]=10; nSigmaTPC[4]=10;
      nSigmaTOF[0]=10; nSigmaTOF[1]=10; nSigmaTOF[2]=10; nSigmaTOF[3]=10; nSigmaTOF[4]=10;
      fTempStruct[myTracks].fTOFhit = kFALSE;// default
      Float_t signalTPC=0, signalTOF=0;
      Double_t integratedTimesTOF[10]={0};

    
      Bool_t DoPIDWorkAround=kTRUE;
      //if(fFilterBit == 7) DoPIDWorkAround=kTRUE;
      if(fMCcase && fCollisionType!=0) DoPIDWorkAround=kFALSE;
      if(DoPIDWorkAround==kFALSE && fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kPion)) < 900) {
	nSigmaTPC[0]=fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kElectron);
	nSigmaTPC[1]=fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kMuon);
	nSigmaTPC[2]=fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kPion);
	nSigmaTPC[3]=fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kKaon);
	nSigmaTPC[4]=fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kProton);
	//
	nSigmaTOF[0]=fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kElectron);
	nSigmaTOF[1]=fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kMuon);
	nSigmaTOF[2]=fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kPion);
	nSigmaTOF[3]=fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kKaon);
	nSigmaTOF[4]=fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kProton);
	signalTPC = aodtrack->GetTPCsignal();
	if( (status&AliESDtrack::kTOFpid)!=0 && (status&AliESDtrack::kTIME)!=0 && (status&AliESDtrack::kTOFout)!=0 && (status&AliESDtrack::kTOFmismatch)<=0){// good tof hit
	  fTempStruct[myTracks].fTOFhit = kTRUE;
	  signalTOF = aodtrack->GetTOFsignal();
	  aodtrack->GetIntegratedTimes(integratedTimesTOF);
	}else fTempStruct[myTracks].fTOFhit = kFALSE;
	
      }else {// FilterBit 7 PID workaround
	
	for(Int_t j = 0; j < fAOD->GetNumberOfTracks(); j++) {
	  AliAODTrack* aodTrack2 = (AliAODTrack*)fAOD->GetTrack(j);
	  if (!aodTrack2) continue;
	  if(aodtrack->GetID() != (-aodTrack2->GetID() - 1)) continue;// (-aodTrack2->GetID() - 1)
	  
	  UInt_t status2=aodTrack2->GetStatus();
	  
	  nSigmaTPC[0]=fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kElectron);
	  nSigmaTPC[1]=fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kMuon);
	  nSigmaTPC[2]=fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kPion);
	  nSigmaTPC[3]=fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kKaon);
	  nSigmaTPC[4]=fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kProton);
	  //
	  nSigmaTOF[0]=fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kElectron);
	  nSigmaTOF[1]=fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kMuon);
	  nSigmaTOF[2]=fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kPion);
	  nSigmaTOF[3]=fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kKaon);
	  nSigmaTOF[4]=fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kProton);
	  signalTPC = aodTrack2->GetTPCsignal();
	  
	  if( (status2&AliESDtrack::kTOFpid)!=0 && (status2&AliESDtrack::kTIME)!=0 && (status2&AliESDtrack::kTOFout)!=0 && (status2&AliESDtrack::kTOFmismatch)<=0){// good tof hit
	    fTempStruct[myTracks].fTOFhit = kTRUE;
	    signalTOF = aodTrack2->GetTOFsignal();
	    aodTrack2->GetIntegratedTimes(integratedTimesTOF);
	  }else fTempStruct[myTracks].fTOFhit = kFALSE;
	  
	  
	}// aodTrack2
      }// FilterBit 7 PID workaround
      
     
      ///////////////////
      ((TH3F*)fOutputList->FindObject("fTPCResponse"))->Fill(centralityPercentile, fTempStruct[myTracks].fMom, signalTPC);
      if(fTempStruct[myTracks].fTOFhit) {
	((TH3F*)fOutputList->FindObject("fTOFResponse"))->Fill(centralityPercentile, fTempStruct[myTracks].fMom, signalTOF - integratedTimesTOF[3]);
      }
      ///////////////////
      
      // Use TOF if good hit and above threshold
      if(fTempStruct[myTracks].fTOFhit && fTempStruct[myTracks].fMom > fTPCTOFboundry){
	if(fabs(nSigmaTOF[0])<fSigmaCutTOF) fTempStruct[myTracks].fElectron = kTRUE;// Electron candidate
	if(fabs(nSigmaTOF[2])<fSigmaCutTOF) fTempStruct[myTracks].fPion = kTRUE;// Pion candidate
	if(fabs(nSigmaTOF[3])<fSigmaCutTOF) fTempStruct[myTracks].fKaon = kTRUE;// Kaon candidate
	if(fabs(nSigmaTOF[4])<fSigmaCutTOF) fTempStruct[myTracks].fProton = kTRUE;// Proton candidate
      }else {// TPC info instead
	if(fabs(nSigmaTPC[0])<fSigmaCutTPC) fTempStruct[myTracks].fElectron = kTRUE;// Electron candidate
	if(fabs(nSigmaTPC[2])<fSigmaCutTPC) fTempStruct[myTracks].fPion = kTRUE;// Pion candidate
	if(fabs(nSigmaTPC[3])<fSigmaCutTPC) fTempStruct[myTracks].fKaon = kTRUE;// Kaon candidate
	if(fabs(nSigmaTPC[4])<fSigmaCutTPC) fTempStruct[myTracks].fProton = kTRUE;// Proton candidate
      }
      
    
      // Ensure there is only 1 candidate per track
      if(fTempStruct[myTracks].fElectron && fTempStruct[myTracks].fMom < 0.45) continue;// Remove electron band
      if(!fTempStruct[myTracks].fPion && !fTempStruct[myTracks].fKaon && !fTempStruct[myTracks].fProton) continue;
      if(!fTempStruct[myTracks].fPion) continue;// only pions
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;
     


   
      if(fTempStruct[myTracks].fCharge==+1) {
	((TH2F*)fOutputList->FindObject("fDCAxyDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }else {
	((TH2F*)fOutputList->FindObject("fDCAxyDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }
     
      ((TH3F*)fOutputList->FindObject("fChPhiPtDist"))->Fill(aodtrack->Charge(), aodtrack->Phi(), aodtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fChPtEtaDist"))->Fill(aodtrack->Charge(), aodtrack->Pt(), aodtrack->Eta());
      ((TH2F*)fOutputList->FindObject("fCentPtDist"))->Fill(int(centralityPercentile/5.), aodtrack->Pt());
      ((TH2F*)fOutputList->FindObject("fCentEtaDist"))->Fill(int(centralityPercentile/5.), aodtrack->Eta());

      ((TH2D*)fOutputList->FindObject("fDistTPCNclsFindable"))->Fill(aodtrack->Pt(), aodtrack->GetTPCNclsF());
      ((TProfile*)fOutputList->FindObject("fProfileTPCNclsFindable"))->Fill(aodtrack->Pt(), aodtrack->GetTPCNclsF());
      //
      ((TH2D*)fOutputList->FindObject("fDistTPCNclsCrossed"))->Fill(aodtrack->Pt(), aodtrack->GetTPCNCrossedRows());
      ((TProfile*)fOutputList->FindObject("fProfileTPCNclsCrossed"))->Fill(aodtrack->Pt(), aodtrack->GetTPCNCrossedRows());
      //
      if(aodtrack->GetTPCNclsF() > 0){
	((TH2D*)fOutputList->FindObject("fDistTPCNclsFindableRatio"))->Fill(aodtrack->Pt(), double(aodtrack->GetTPCNcls())/double(aodtrack->GetTPCNclsF()));
	((TProfile*)fOutputList->FindObject("fProfileTPCNclsFindableRatio"))->Fill(aodtrack->Pt(), double(aodtrack->GetTPCNcls())/double(aodtrack->GetTPCNclsF()));
      }      
      // 
      ((TH2D*)fOutputList->FindObject("fDistTPCNclsCrossedRatio"))->Fill(aodtrack->Pt(), aodtrack->GetTPCFoundFraction());
      ((TProfile*)fOutputList->FindObject("fProfileTPCNclsCrossedRatio"))->Fill(aodtrack->Pt(), aodtrack->GetTPCFoundFraction());
      

      if(fTempStruct[myTracks].fPion) {// pions
	fTempStruct[myTracks].fEaccepted = sqrt(pow(fTempStruct[myTracks].fMom,2) + pow(fTrueMassPi,2)); 
	fTempStruct[myTracks].fKey = 1;
      }else if(fTempStruct[myTracks].fKaon){// kaons
	fTempStruct[myTracks].fEaccepted = sqrt(pow(fTempStruct[myTracks].fMom,2) + pow(fTrueMassK,2));;
	fTempStruct[myTracks].fKey = 10;
      }else{// protons
	fTempStruct[myTracks].fEaccepted = sqrt(pow(fTempStruct[myTracks].fMom,2) + pow(fTrueMassP,2));;
	fTempStruct[myTracks].fKey = 100;
      }
      
           

      if(aodtrack->Charge() > 0) positiveTracks++;
      else negativeTracks++;
      
      if(fTempStruct[myTracks].fPion) pionCount++;
      if(fTempStruct[myTracks].fKaon) kaonCount++;
      if(fTempStruct[myTracks].fProton) protonCount++;

      myTracks++;
      
      if(fMCcase){// muon mothers
	AliAODMCParticle *tempMCTrack=(AliAODMCParticle*)mcArray->At(abs(aodtrack->GetLabel()));
	if(abs(tempMCTrack->GetPdgCode())==13 && tempMCTrack->GetMother()>0){// muons
	  AliAODMCParticle *parent=(AliAODMCParticle*)mcArray->At(tempMCTrack->GetMother());
	  if(parent->IsPhysicalPrimary()){
	    ((TH1D*)fOutputList->FindObject("fMuonParents"))->Fill(abs(parent->GetPdgCode()));
	  }else ((TH1D*)fOutputList->FindObject("fSecondaryMuonParents"))->Fill(abs(parent->GetPdgCode()));
	}
	((TH1D*)fOutputList->FindObject("fPionCandidates"))->Fill(abs(tempMCTrack->GetPdgCode()));
      }
    }
    //cout<<"kinkcount = "<<kinkcount<<"   pionkinks = "<<pionkinks<<"   primarypionkinks = "<<primarypionkinks<<endl;
  }else {// ESD tracks
    cout<<"ESDs not supported currently"<<endl;
    return;
  }

  // Generator info only
  if(fMCcase && fGeneratorOnly){
    myTracks=0; pionCount=0; kaonCount=0; protonCount=0;// reset track counters
    for(Int_t mctrackN=0; mctrackN<mcArray->GetEntriesFast(); mctrackN++){
      if(myTracks >= fMultLimit) {cout<<"More tracks than Track Limit"<<endl; return;}
      if(myTracks >= 1300) continue;// additional cut to limit high mult events which exceed pair # limits
      
      AliAODMCParticle *mcParticle = (AliAODMCParticle*)mcArray->At(mctrackN);
      if(!mcParticle) continue;
      if(fabs(mcParticle->Eta())>0.8) continue;
      if(mcParticle->Charge()!=-3 && mcParticle->Charge()!=+3) continue;// x3 by convention
      if(mcParticle->Pt() < 0.16 || mcParticle->Pt() > 1.0) continue;
      if(!mcParticle->IsPrimary()) continue;
      if(!mcParticle->IsPhysicalPrimary()) continue;
      if(abs(mcParticle->GetPdgCode())!=211) continue;
      
      fTempStruct[myTracks].fP[0] = mcParticle->Px();
      fTempStruct[myTracks].fP[1] = mcParticle->Py();
      fTempStruct[myTracks].fP[2] = mcParticle->Pz();
      fTempStruct[myTracks].fX[0] = 0.; fTempStruct[myTracks].fX[1] = 0.; fTempStruct[myTracks].fX[2] = 0.;
      
      fTempStruct[myTracks].fId = myTracks;// use my track counter 
      fTempStruct[myTracks].fLabel = mctrackN;
      fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
      if(fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2*PI;
      fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0],2) + pow(fTempStruct[myTracks].fP[1],2));
      fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt,2) + pow(fTempStruct[myTracks].fP[2],2) );
      fTempStruct[myTracks].fEta = mcParticle->Eta();
      fTempStruct[myTracks].fCharge = int(mcParticle->Charge()/3.);
      fTempStruct[myTracks].fDCAXY = 0.;
      fTempStruct[myTracks].fDCAZ = 0.;
      fTempStruct[myTracks].fDCA = 0.;
      fTempStruct[myTracks].fPion = kTRUE;
      fTempStruct[myTracks].fEaccepted = sqrt(pow(fTempStruct[myTracks].fMom,2) + pow(fTrueMassPi,2)); 
      fTempStruct[myTracks].fKey = 1;
      
      myTracks++;
      pionCount++;
    }
  }
  
  if(myTracks >= 1) {
    ((TH1F*)fOutputList->FindObject("fMultDist3"))->Fill(myTracks);
  }
 
 
  //cout<<"There are "<<myTracks<<"  myTracks"<<endl;
  //cout<<"pionCount = "<<pionCount<<"   kaonCount = "<<kaonCount<<"   protonCount = "<<protonCount<<endl;
  //return;

  /////////////////////////////////////////
  // Pion Multiplicity Cut (To ensure all Correlation orders are present in each event)
  if(myTracks < 4) {return;}
  /////////////////////////////////////////
 

  ////////////////////////////////
  ///////////////////////////////
  // Mbin determination
  //
  // Mbin set to Pion Count Only for pp!!!!!!!
  fMbin=-1;
  if(fCollisionType!=0){
    
    if(fCollisionType==1){// p-Pb
      if(pionCount >= fMultLimits[3] && pionCount < fMultLimits[10]) fMbin=0;// only 1 bin
    }
    if(fCollisionType==2){// pp
      if(pionCount >= fMultLimits[2] && pionCount < fMultLimits[10]) fMbin=0;// only 1 bin
    }

    for(Int_t i=0; i<fMbinsMixing; i++){// event-mixing M bin
      if( ( pionCount >= fMultLimits[i]) && ( pionCount < fMultLimits[i+1]) ){
	MbinMixing=i;// 0 = lowest mult
	break;
      }
    }

  }else{
    for(Int_t i=0; i<fCentBins; i++){// correlation analysis M bin
      if( (centralityPercentile >= cStart+i*cStep) && (centralityPercentile < cStart+(i+1)*cStep) ){
	fMbin=i;// 0 = most central
	break;
      }
    }
    for(Int_t i=0; i<fMbinsMixing; i++){// event-mixing M bin
      if( (centralityPercentile >= cStart+i*cStepMixing) && (centralityPercentile < cStart+(i+1)*cStepMixing) ){
	MbinMixing=i;// 0 = most central
	break;
      }
    }
  }
  
  if(fMbin==-1) {return;}
  
  
  //////////////////////////////////////////////////////////////////////////
  // c2 and d2 calculation
  for(Int_t boi=0; boi<40; boi++){
    if(countINT[boi]>0) {
      Float_t cumulant2 = cosINT[boi]*cosINT[boi] + sinINT[boi]*sinINT[boi] - countINT[boi];
      cumulant2 /= countINT[boi]*(countINT[boi]-1);
      Float_t Ncomb = countINT[boi]*(countINT[boi]-1);
      ((TH2D*)fOutputList->FindObject("fcumulant2INT"))->Fill(fMbin+1, boi+1, cumulant2 * Ncomb);
      ((TH2D*)fOutputList->FindObject("fSingleSumCosINT"))->Fill(fMbin+1, boi+1, cosINT[boi] / countINT[boi] * Ncomb);
      ((TH2D*)fOutputList->FindObject("fSingleSumSinINT"))->Fill(fMbin+1, boi+1, sinINT[boi] / countINT[boi] * Ncomb);
      ((TH2D*)fOutputList->FindObject("fSingleSumENINT"))->Fill(fMbin+1, boi+1, Ncomb);
    }
    if(countDIFF[boi]>0) {
      Float_t Ncomb = countDIFF[boi]*countINT[boi];
      if(countINT[boi]>0) {
	Float_t cumulant2 = cosDIFF[boi]*cosINT[boi] + sinDIFF[boi]*sinINT[boi];
	cumulant2 /= countDIFF[boi]*countINT[boi];
	((TH2D*)fOutputList->FindObject("fcumulant2DIFF"))->Fill(fMbin+1, boi+1, cumulant2 * Ncomb);
	((TH2D*)fOutputList->FindObject("fcumulant2EN"))->Fill(fMbin+1, boi+1, Ncomb);
      }
      ((TH2D*)fOutputList->FindObject("fSingleSumCosDIFF"))->Fill(fMbin+1, boi+1, cosDIFF[boi] / countDIFF[boi] * Ncomb);
      ((TH2D*)fOutputList->FindObject("fSingleSumSinDIFF"))->Fill(fMbin+1, boi+1, sinDIFF[boi] / countDIFF[boi] * Ncomb);
      ((TH2D*)fOutputList->FindObject("fSingleSumENDIFF"))->Fill(fMbin+1, boi+1, Ncomb);
    }
  }
  //
  if(countTotalINT > 0){
    Float_t cumulant2 = cosTotalINT*cosTotalINT + sinTotalINT*sinTotalINT - countTotalINT;
    cumulant2 /= countTotalINT*(countTotalINT-1);
    Float_t Ncomb = countTotalINT*(countTotalINT-1);
    ((TH1D*)fOutputList->FindObject("fcumulant2TotalINT"))->Fill(fMbin+1, cumulant2 * Ncomb);
    ((TH1D*)fOutputList->FindObject("fSingleSumCosTotalINT"))->Fill(fMbin+1, cosTotalINT / countTotalINT * Ncomb);
    ((TH1D*)fOutputList->FindObject("fSingleSumSinTotalINT"))->Fill(fMbin+1, sinTotalINT / countTotalINT * Ncomb);
    ((TH1D*)fOutputList->FindObject("fSingleSumTotalEN"))->Fill(fMbin+1, Ncomb);
  }

  // q1 vector
  Float_t Q1x[4]={0};
  Float_t Q1y[4]={0};
  Float_t Mq[4]={0};
  Float_t qVect1[4]={0};
  Float_t Psi1[4]={0};
  Int_t qindex=0;
  // q2 vector
  Float_t Q2x[4]={0};
  Float_t Q2y[4]={0};
  Float_t qVect2[4]={0};
  Float_t Psi2[4]={0};
  for(Int_t i=0; i<myTracks; i++){
    if(fChargeSelection && fTempStruct[i].fCharge !=-1) continue;

    if(fTempStruct[i].fPt < 0.25) qindex=0;// was 0.28
    else if(fTempStruct[i].fPt < 0.35) qindex=1;// was 0.4
    else if(fTempStruct[i].fPt < 0.5) qindex=2;
    else qindex=3;
    
    Q1x[qindex] += cos(fTempStruct[i].fPhi);
    Q1y[qindex] += sin(fTempStruct[i].fPhi);
    Q2x[qindex] += cos(2*fTempStruct[i].fPhi);
    Q2y[qindex] += sin(2*fTempStruct[i].fPhi);
    Mq[qindex]++;
  }
  for(Int_t i=0; i<4; i++){ 
    qVect1[i] = sqrt(pow(Q1x[i],2)+pow(Q1y[i],2)); 
    qVect2[i] = sqrt(pow(Q2x[i],2)+pow(Q2y[i],2));
    if(Mq[i] > 0) {qVect1[i] /= sqrt(Mq[i]); qVect2[i] /= sqrt(Mq[i]);}
    ((TH3F*)fOutputList->FindObject("fq1Dist"))->Fill(fMbin+1, i+1, qVect1[i]);
    ((TH3F*)fOutputList->FindObject("fq2Dist"))->Fill(fMbin+1, i+1, qVect2[i]);
    Psi1[i] = atan2(Q1y[i],Q1x[i]) + PI;// 0 to +2PI
    Psi2[i] = atan2(Q2y[i],Q2x[i]) / 2.;// -PI/2 to +PI/2
    Psi2[i] = fabs(Psi2[i]);// 0 to +PI/2
  }
  ((TH2D*)fOutputList->FindObject("fLowPtDist"))->Fill(fMbin+1, Mq[fq2Index]);
  /*((TH3D*)fOutputList->FindObject("fPtMultCrossing1"))->Fill(fMbin+1, Mq[0], Mq[2]);
  ((TH3D*)fOutputList->FindObject("fPtMultCrossing2"))->Fill(fMbin+1, Mq[1], Mq[2]);
  ((TH3D*)fOutputList->FindObject("fPtMultCrossing3"))->Fill(fMbin+1, Mq[2], Mq[3]);*/
  //
  if(fq2Binning){// bin in q2
    if(qVect2[fq2Index] < fq2CutLow) fEDbin = 0;
    else fEDbin = 1;
  
    Int_t Inq1=0;
    if(qVect1[fq2Index] > 0.5 && qVect1[fq2Index] < 1.5) Inq1=1;
    else Inq1=2;
    //
    Int_t Inq2=0;
    if(qVect2[fq2Index] > 0.5 && qVect2[fq2Index] < 1.5) Inq2=1;
    else Inq2=2;
    //
    //if(!fMCcase) mixingEDbin = (Inq1*3*3*3) + int(Psi1[fq2Index]/(2*PI/3.) - 0.000001)*3*3 + (Inq2*3) + int(Psi2[fq2Index]/(PI/6.) - 0.000001);// q1 and q2
    if(!fMCcase) mixingEDbin = (Inq2*3) + int(Psi2[fq2Index]/(PI/6.) - 0.000001);// q2 only
  }
  //
  if(fLowMultBinning){
    Float_t HighPtMult_L = 300 - 55*fMbin;// approximate binning for 0.35-0.5 counting interval
    Float_t HighPtMult_H = 330 - 55*fMbin;
    ////////////////////////////////////////////////////////
    if(Mq[2] < HighPtMult_L || Mq[2] > HighPtMult_H) return;// remove event completely 
    ////////////////////////////////////////////////////////
    Float_t meanLowPtMult[10] ={175, 147, 125, 103, 86, 71, 58, 48, 38, 30};// was 170. - 25.*fMbin;// approximate mean values vs. fMbin
    //Float_t sigmaLowPtMult = 0.1*meanLowPtMult;// approximate sigma values
    if(Mq[fq2Index] < meanLowPtMult[fMbin]) fEDbin=0;
    else fEDbin = 1;
    //
    mixingEDbin = fEDbin;
  }
  
  //////////////////////////////////////////////////////////////////////////
  

  
  ///////////////////
  // can only be called after fMbin has been set
  // Radius parameter only matters for Monte-Carlo data
  SetFSIindex(fRMax);
  ///////////////////  
  
  Int_t rBinForTPNMomRes = 10;
  if(fMbin==0) {rBinForTPNMomRes=10;}// 10 fm with EW (fRMax should be 11 for normal running)
  else if(fMbin==1) {rBinForTPNMomRes=9;}
  else if(fMbin<=3) {rBinForTPNMomRes=8;}
  else if(fMbin<=5) {rBinForTPNMomRes=7;}
  else {rBinForTPNMomRes=6;}

  //////////////////////////////////////////////////
  if(!fq2Binning && !fLowMultBinning && fQdirectionBinning==0) fEDbin=0;// Extra Dimension bin (Kt3, q2,....)
  //////////////////////////////////////////////////
  
  
  
  ((TH1F*)fOutputList->FindObject("fEvents1"))->Fill(fMbin+1);
  ((TProfile*)fOutputList->FindObject("fAvgMult"))->Fill(fMbin+1., pionCount);

  ////////////////////////////////////
  // Add event to buffer if > 0 tracks
  if(myTracks > 0){
    fEC[mixingEDbin][MbinMixing]->FIFOShift();
    (fEvt) = fEC[mixingEDbin][MbinMixing]->fEvtStr;
    (fEvt)->fNtracks = myTracks;
    (fEvt)->fFillStatus = 1;
    for(Int_t i=0; i<myTracks; i++) (fEvt)->fTracks[i] = fTempStruct[i];
    if(fMCcase){
      (fEvt)->fMCarraySize = mcArray->GetEntriesFast();
      for(Int_t i=0; i<mcArray->GetEntriesFast(); i++) {
	AliAODMCParticle *tempMCTrack = (AliAODMCParticle*)mcArray->At(i);
	(fEvt)->fMCtracks[i].fPx = tempMCTrack->Px();
	(fEvt)->fMCtracks[i].fPy = tempMCTrack->Py();
	(fEvt)->fMCtracks[i].fPz = tempMCTrack->Pz();
	(fEvt)->fMCtracks[i].fPtot = sqrt(pow(tempMCTrack->Px(),2)+pow(tempMCTrack->Py(),2)+pow(tempMCTrack->Pz(),2));
	(fEvt)->fMCtracks[i].fPdgCode = tempMCTrack->GetPdgCode();
	(fEvt)->fMCtracks[i].fMotherLabel = tempMCTrack->GetMother();
      }	
    }
  }
  
  
  
  Float_t qinv12=0, qinv13=0, qinv14=0, qinv23=0, qinv24=0, qinv34=0;
  Float_t qout=0, qside=0, qlong=0;
  Float_t qout12=0, qside12=0, qlong12=0;
  Float_t qout13=0, qside13=0, qlong13=0;
  Float_t qout14=0, qside14=0, qlong14=0;
  Float_t qout23=0, qside23=0, qlong23=0;
  Float_t qout24=0, qside24=0, qlong24=0;
  Float_t qout34=0, qside34=0, qlong34=0;
  Float_t kT12=0, kT13=0, kT14=0, kT23=0, kT24=0, kT34=0;
  Float_t q3=0, q3MC=0;
  Float_t q4=0, q4MC=0;
  Int_t ch1=0, ch2=0, ch3=0, ch4=0;
  Int_t bin1=0, bin2=0, bin3=0, bin4=0;
  Float_t pVect1[4]={0}; 
  Float_t pVect2[4]={0};
  Float_t pVect3[4]={0};
  Float_t pVect4[4]={0};
  Float_t pVect1MC[4]={0}; 
  Float_t pVect2MC[4]={0};
  Float_t pVect3MC[4]={0};
  Float_t pVect4MC[4]={0};
  Float_t Pparent1[4]={0};
  Float_t Pparent2[4]={0};
  Float_t Pparent3[4]={0};
  Float_t Pparent4[4]={0};
  Float_t weight12=0, weight13=0, weight14=0, weight23=0, weight24=0, weight34=0;
  Float_t weight12Err=0, weight13Err=0, weight14Err=0, weight23Err=0, weight24Err=0, weight34Err=0;
  Float_t weight12CC[3]={0};
  Float_t weight13CC[3]={0};
  Float_t weight14CC[3]={0};
  Float_t weight23CC[3]={0};
  Float_t weight24CC[3]={0};
  Float_t weight34CC[3]={0};
  //Float_t weight12CC_e=0, weight13CC_e=0, weight14CC_e=0, weight23CC_e=0, weight24CC_e=0, weight34CC_e=0;
  Float_t weightTotal=0;//, weightTotalErr=0;
  Float_t weightPrime=0, weightPrimePrime=0, weightCumulant=0;
  Float_t qinv12MC=0, qinv13MC=0, qinv14MC=0, qinv23MC=0, qinv24MC=0, qinv34MC=0; 
  Float_t parentQinv12=0, parentQinv13=0, parentQinv14=0, parentQinv23=0, parentQinv24=0, parentQinv34=0;
  Float_t parentQ3=0;
  Float_t FSICorr12=0, FSICorr13=0, FSICorr14=0, FSICorr23=0, FSICorr24=0, FSICorr34=0;
  Bool_t pionParent1=kFALSE, pionParent2=kFALSE, pionParent3=kFALSE, pionParent4=kFALSE;
  Bool_t FilledMCpair12=kFALSE, FilledMCtriplet123=kFALSE;
  Bool_t Positive1stTripletWeights=kTRUE, Positive2ndTripletWeights=kTRUE;
  Float_t T12=0, T13=0, T14=0, T23=0, T24=0, T34=0;
  Float_t t12=0, t13=0, t14=0, t23=0, t24=0, t34=0;
  Int_t momBin12=1, momBin13=1, momBin14=1, momBin23=1, momBin24=1, momBin34=1;
  Float_t MomResCorr12=1.0, MomResCorr13=1.0, MomResCorr14=1.0, MomResCorr23=1.0, MomResCorr24=1.0, MomResCorr34=1.0;
  //
  AliAODMCParticle *mcParticle1=0x0;
  AliAODMCParticle *mcParticle2=0x0;
  

  ////////////////////
  Int_t EDindex3=0, EDindex4=0;

  // reset to defaults
  for(Int_t i=0; i<fMultLimit; i++) {
    fLowQPairSwitch_E0E0[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E1[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E2[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E0E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E1[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E2[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E1E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fLowQPairSwitch_E2E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    //
    fNormQPairSwitch_E0E0[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E1[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E2[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E0E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E1[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E2[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E1E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
    fNormQPairSwitch_E2E3[i]->Set(kMultLimitPbPb,fDefaultsCharSwitch);
  }
 
  
  //////////////////////////////////////////
  // make low-q pair storage and normalization-pair storage
  // 
  for(Int_t en1=0; en1<=2; en1++){// 1st event number (en1=0 is the same event as current event)
    for(Int_t en2=en1; en2<=3; en2++){// 2nd event number (en2=0 is the same event as current event)
      if(en1>1 && en1==en2) continue;
      
      for (Int_t i=0; i<(fEvt+en1)->fNtracks; i++) {// 1st particle
	for (Int_t j=i+1; j<(fEvt+en2)->fNtracks; j++) {// 2nd particle
	  
	  
	  pVect1[0]=(fEvt+en1)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	  pVect1[1]=(fEvt+en1)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	  pVect1[2]=(fEvt+en1)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	  pVect1[3]=(fEvt+en1)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	  ch1 = Int_t(((fEvt+en1)->fTracks[i].fCharge + 1)/2.);
	  ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	  
	  qinv12 = GetQinv(pVect1, pVect2);
	  kT12 = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;
	  SetFillBins2(ch1, ch2, bin1, bin2);
	  
	  if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting/merging cut)
	  if(ch1 == ch2 && !fGeneratorOnly){
	    Int_t tempChGroup[2]={0,0};
	    if(en1==0 && en2==1) ((TH1F*)fOutputList->FindObject("fTotalPairsWeighting"))->Fill(qinv12, MCWeight(tempChGroup, 10, ffcSqMRC, qinv12, 0.));
	    if(!AcceptPair((fEvt+en1)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	      if(en1==0 && en2==0) ((TH1F*)fOutputList->FindObject("fRejectedPairs"))->Fill(qinv12);
	      continue;
	    }
	    if(en1==0 && en2==1) ((TH1F*)fOutputList->FindObject("fAcceptedPairsWeighting"))->Fill(qinv12, MCWeight(tempChGroup, 10, ffcSqMRC, qinv12, 0.));
	  }
	  if(fMixedChargeCut && ch1 != ch2 && !fGeneratorOnly && !fMCcase){// remove +- low-q pairs to keep balance between ++ and +- contributions to multi-particle Q3,Q4 projections
	    Int_t tempChGroup[2]={0,1};
	    if(en1==0 && en2==1) ((TH1F*)fOutputList->FindObject("fTotalPairsWeightingMC"))->Fill(qinv12, MCWeight(tempChGroup, 10, ffcSqMRC, qinv12, 0.));
	    if(!AcceptPairPM((fEvt+en1)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	      if(en1==0 && en2==0) ((TH1F*)fOutputList->FindObject("fRejectedPairsMC"))->Fill(qinv12);
	      continue;
	    }
	    if(en1==0 && en2==1) ((TH1F*)fOutputList->FindObject("fAcceptedPairsWeightingMC"))->Fill(qinv12, MCWeight(tempChGroup, 10, ffcSqMRC, qinv12, 0.));
	  }
	  
	  GetQosl(pVect1, pVect2, qout, qside, qlong);
	  if( (en1+en2==0)) {
	    if(!fGenerateSignal) Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[0].fTerms2->Fill(kT12, qinv12);
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[0].fTerms2QW->Fill(kT12, qinv12, qinv12);
	    if(qinv12<0.1) ((TH2D*)fOutputList->FindObject("fkTDist"))->Fill(fMbin+1, kT12);
	    // osl frame
	    if((kT12 > 0.2) && (kT12 < 0.3)){
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[0].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[0].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    if((kT12 > 0.6) && (kT12 < 0.7)){  
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[1].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[1].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    // unit mult bins
	    /*if( (fEvt+en1)->fNtracks%100==0){
	      Int_t kTindex=0;
	      if(kT12>0.3) kTindex=1;
	      Int_t UnitMultBin = int((fEvt+en1)->fNtracks / 100.) + 1;
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[kTindex].TwoPT[0].fUnitMultBin->Fill(UnitMultBin, qinv12);
	      }*/
	  
	  }
	  if( (en1+en2==1)) {
	    if(!fGenerateSignal) {
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fTerms2->Fill(kT12, qinv12);
	      // Build 2-pion correlations from previous 2-pion tabulations
	      if(!fTabulatePairs && bin1==bin2 && qinv12<fQcut){
		GetWeight(pVect1, pVect2, weight12, weight12Err);
		//
		if(fOnlineCorrection){
		  momBin12 = fMomResC2SC->GetYaxis()->FindBin(qinv12);
		  Float_t MuonCorr12 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin12);
		  MomResCorr12 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin12);
		  FSICorr12 = FSICorrelation(ch1,ch2, qinv12);
		  weight12CC[2] = ((weight12+1)*MomResCorr12 - ffcSq*FSICorr12 - (1-ffcSq));
		  weight12CC[2] /= FSICorr12*ffcSq;
		  weight12CC[2] *= MuonCorr12;
		}else weight12CC[2] = weight12;

		if(weight12CC[2]<0) weight12CC[2]=0;
		
		Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fBuild->Fill(4, kT12, qinv12, 1);
		Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fBuild->Fill(5, kT12, qinv12, weight12CC[2]);
	      }
	    }
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fTerms2QW->Fill(kT12, qinv12, qinv12);
	    // osl frame
	    if((kT12 > 0.2) && (kT12 < 0.3)){  
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[0].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[0].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    if((kT12 > 0.6) && (kT12 < 0.7)){  
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[1].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[en2].OSL_ktbin[1].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    // unit mult bins
	    /*if( (fEvt+en1)->fNtracks%100==0){
	      Int_t kTindex=0;
	      if(kT12>0.3) kTindex=1;
	      Int_t UnitMultBin = int((fEvt+en1)->fNtracks / 100.) + 1;
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[kTindex].TwoPT[1].fUnitMultBin->Fill(UnitMultBin, qinv12);
	      }*/
	  }

	  /////////////////////////////////////////////////////
	  if(fTabulatePairs && en1==0 && en2<=1 && bin1==bin2){
	    if(fChargeSelection && (fEvt+en1)->fTracks[i].fCharge !=-1) continue;
	    Float_t kY = 0;
	    Int_t kTbin=-1, kYbin=-1;
	    Bool_t PairToReject=kFALSE;
	    if((fEvt+en1)->fTracks[i].fPt < fMinPt || (fEvt+en1)->fTracks[i].fPt > fMaxPt) PairToReject=kTRUE;
	    if((fEvt+en2)->fTracks[j].fPt < fMinPt || (fEvt+en2)->fTracks[j].fPt > fMaxPt) PairToReject=kTRUE;
	    if(!PairToReject){
	      for(Int_t kIt=0; kIt<fKbinsT; kIt++) {if(kT12 < (fKmiddleT[kIt] + fKstepT[kIt]/2.)) {kTbin = kIt; break;}} 
	      for(Int_t kIt=0; kIt<fKbinsY; kIt++) {if(kY < (fKmiddleY[kIt] + fKstepY[kIt]/2.)) {kYbin = kIt; break;}}
	      if((kTbin<0) || (kYbin<0)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	      if((kTbin>=fKbinsT) || (kYbin>=fKbinsY)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	      if(fGenerateSignal && en2==0) {
		Int_t chGroup2[2]={ch1,ch2};
		Float_t WInput = MCWeight(chGroup2, fRMax, ffcSqMRC, qinv12, kT12);
		KT[kTbin].KY[kYbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fTerms2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong), WInput);
	      }else KT[kTbin].KY[kYbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fTerms2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    }
	  }
	  
	  //////////////////////////////////////////////////////////////////////////////
	 
	  if(qinv12 <= fQcut) {
	    if(en1==0 && en2==0) {fLowQPairSwitch_E0E0[i]->AddAt('1',j);}
	    if(en1==0 && en2==1) {fLowQPairSwitch_E0E1[i]->AddAt('1',j);}
	    if(en1==0 && en2==2) {fLowQPairSwitch_E0E2[i]->AddAt('1',j);}
	    if(en1==0 && en2==3) {fLowQPairSwitch_E0E3[i]->AddAt('1',j);}
	    if(en1==1 && en2==1) {fLowQPairSwitch_E1E1[i]->AddAt('1',j);}
	    if(en1==1 && en2==2) {fLowQPairSwitch_E1E2[i]->AddAt('1',j);}
	    if(en1==1 && en2==3) {fLowQPairSwitch_E1E3[i]->AddAt('1',j);}
	    if(en1==2 && en2==3) {fLowQPairSwitch_E2E3[i]->AddAt('1',j);}
	  }
	  if((qinv12 >= fNormQcutLow) && (qinv12 < fNormQcutHigh)) {
	    if(en1==0 && en2==0) {fNormQPairSwitch_E0E0[i]->AddAt('1',j);}
	    if(en1==0 && en2==1) {fNormQPairSwitch_E0E1[i]->AddAt('1',j);}
	    if(en1==0 && en2==2) {fNormQPairSwitch_E0E2[i]->AddAt('1',j);}
	    if(en1==0 && en2==3) {fNormQPairSwitch_E0E3[i]->AddAt('1',j);}
	    if(en1==1 && en2==1) {fNormQPairSwitch_E1E1[i]->AddAt('1',j);}
	    if(en1==1 && en2==2) {fNormQPairSwitch_E1E2[i]->AddAt('1',j);}
	    if(en1==1 && en2==3) {fNormQPairSwitch_E1E3[i]->AddAt('1',j);}
	    if(en1==2 && en2==3) {fNormQPairSwitch_E2E3[i]->AddAt('1',j);}
	  }
	  
	}
      }
    }
  }
    
 
  ///////////////////////////////////////////////////  
  // Do not use pairs from events with too many pairs
  
  ((TH1F*)fOutputList->FindObject("fEvents2"))->Fill(fMbin+1);
  
  ///////////////////////////////////////////////////
  
 

  if(fTabulatePairs) return;

  //TF1 *SCpairWeight = new TF1("SCpairWeight","[0] + [1]*x + [2]*exp(-[3]*x)",0,0.2);// same-charge pair weight for monte-carlo data without two-track cuts.
  //SCpairWeight->FixParameter(0, 0.959);
  //SCpairWeight->FixParameter(1, 0.278);
  //SCpairWeight->FixParameter(2, -1.759);
  //SCpairWeight->FixParameter(3, 115.107);

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Normalization counting of 3- and 4-particle terms
  for(Int_t en2=0; en2<=1; en2++){// 2nd event number (en2=0 is the same event as current event)
    for(Int_t en3=en2; en3<=2; en3++){// 3rd event number
      if(en2==0 && en3>2) continue;// not needed config
      if(en2==1 && en3==en2) continue;// not needed config
      for(Int_t en4=en3; en4<=3; en4++){// 4th event number
	if(en3==0 && en4>1) continue;// not needed config
	if(en3==1 && en4==3) continue;// not needed configs
	if(en3==2 && (en2+en3+en4)!=6) continue;// not needed configs
	
	for (Int_t i=0; i<myTracks; i++) {// 1st particle
	  pVect1[1]=(fEvt)->fTracks[i].fP[0];
	  pVect1[2]=(fEvt)->fTracks[i].fP[1];
	  pVect1[3]=(fEvt)->fTracks[i].fP[2];
	  ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	  
	  for (Int_t j=i+1; j<(fEvt+en2)->fNtracks; j++) {// 2nd particle
	    if(en2==0) {if(fNormQPairSwitch_E0E0[i]->At(j)=='0') continue;}
	    else {if(fNormQPairSwitch_E0E1[i]->At(j)=='0') continue;}
	    
	    pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	    pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	    pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	    ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	   
	    for (Int_t k=j+1; k<(fEvt+en3)->fNtracks; k++) {// 3rd particle
	      if(en3==0) {
		if(fNormQPairSwitch_E0E0[i]->At(k)=='0') continue;
		if(fNormQPairSwitch_E0E0[j]->At(k)=='0') continue;
	      }else if(en3==1){
		if(fNormQPairSwitch_E0E1[i]->At(k)=='0') continue;
		if(fNormQPairSwitch_E0E1[j]->At(k)=='0') continue;
	      }else{
		if(fNormQPairSwitch_E0E2[i]->At(k)=='0') continue;
		if(fNormQPairSwitch_E1E2[j]->At(k)=='0') continue;
	      }
	      
	      pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
	      pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
	      pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
	      ch3 = Int_t(((fEvt+en3)->fTracks[k].fCharge + 1)/2.);
	      Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
	      SetFillBins3(ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);
	      
	      Float_t KT3 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
	      if(!fq2Binning && !fLowMultBinning){
		if(KT3<=fKT3transition) EDindex3=0;
		else EDindex3=1;
	      }else{
		EDindex3 = fEDbin;
		if(KT3>fKT3transition) {
		  EDindex3=2+fEDbin;
		}
	      }
	      
	      if(en2==0 && en3==0 && en4==0) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fNorm3->Fill(0);
	      if(en2==1 && en3==2 && en4==3) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fNorm3->Fill(0);
	      if(en2==0 && en3==1 && en4==2) {
		if(fill2) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fNorm3->Fill(0);
		if(fill3) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fNorm3->Fill(0);
		if(fill4) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fNorm3->Fill(0);
	      }
	      
	      
	      for (Int_t l=k+1; l<(fEvt+en4)->fNtracks; l++) {// 4th particle
		if(en4==0){
		  if(fNormQPairSwitch_E0E0[i]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E0E0[j]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E0E0[k]->At(l)=='0') continue;
		}else if(en4==1){
		  if(en3==0){
		    if(fNormQPairSwitch_E0E1[i]->At(l)=='0') continue;
		    if(fNormQPairSwitch_E0E1[j]->At(l)=='0') continue;
		    if(fNormQPairSwitch_E0E1[k]->At(l)=='0') continue;
		  }else{
		    if(fNormQPairSwitch_E0E1[i]->At(l)=='0') continue;
		    if(fNormQPairSwitch_E0E1[j]->At(l)=='0') continue;
		    if(fNormQPairSwitch_E1E1[k]->At(l)=='0') continue;
		  }
		}else if(en4==2){
		  if(fNormQPairSwitch_E0E2[i]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E0E2[j]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E1E2[k]->At(l)=='0') continue;
		}else{
		  if(fNormQPairSwitch_E0E3[i]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E1E3[j]->At(l)=='0') continue;
		  if(fNormQPairSwitch_E2E3[k]->At(l)=='0') continue;
		}
		
		pVect4[1]=(fEvt+en4)->fTracks[l].fP[0];
		pVect4[2]=(fEvt+en4)->fTracks[l].fP[1];
		pVect4[3]=(fEvt+en4)->fTracks[l].fP[2];
		ch4 = Int_t(((fEvt+en4)->fTracks[l].fCharge + 1)/2.);
		Float_t KT4 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1]+pVect4[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2]+pVect4[2],2))/4.;
		
		if(!fq2Binning && !fLowMultBinning){
		  if(KT4<=fKT4transition) EDindex4=0;
		  else EDindex4=1;
		}else{
		  EDindex4 = fEDbin;
		  if(KT4>fKT4transition) {
		    EDindex4=2+fEDbin;
		  }
		}

		Bool_t FillTerms[13]={kFALSE};
		SetFillBins4(ch1, ch2, ch3, ch4, bin1, bin2, bin3, bin4, en2+en3+en4, FillTerms);
		//
		for(int ft=0; ft<13; ft++) {
		  if(FillTerms[ft]) Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[ft].fNorm4->Fill(0.); 
		}
		
		
	      }
	    }
	  }
	}  
	
      }
    }
  }
    


  

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //
    //
    // Start the Main Correlation Analysis
    //
    //
    ///////////////////////////////////////////////////////////////////////
  


    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    for(Int_t en2=0; en2<=1; en2++){// 2nd event number (en2=0 is the same event as current event)
      for(Int_t en3=en2; en3<=2; en3++){// 3rd event number
	if(en2==0 && en3>2) continue;// not needed config
	if(en2==1 && en3==en2) continue;// not needed config
	for(Int_t en4=en3; en4<=3; en4++){// 4th event number
	  if(en3==0 && en4>1) continue;// not needed config
	  if(en3==1 && en4==3) continue;// not needed configs
	  if(en3==2 && (en2+en3+en4)!=6) continue;// not needed configs
	  
	  Int_t ENsum=en2+en3+en4;// 0 or 1 or 3 or 6
	 
	  /////////////////////////////////////////////////////////////
	  for (Int_t i=0; i<myTracks; i++) {// 1st particle
	    pVect1[0]=(fEvt)->fTracks[i].fEaccepted;
	    pVect1[1]=(fEvt)->fTracks[i].fP[0];
	    pVect1[2]=(fEvt)->fTracks[i].fP[1];
	    pVect1[3]=(fEvt)->fTracks[i].fP[2];
	    ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	    if((fEvt)->fTracks[i].fPt < fMinPt) continue; 
	    if((fEvt)->fTracks[i].fPt > fMaxPt) continue;

	    /////////////////////////////////////////////////////////////
	    for (Int_t j=i+1; j<(fEvt+en2)->fNtracks; j++) {// 2nd particle
	      if(en2==0) {if(fLowQPairSwitch_E0E0[i]->At(j)=='0') continue;}
	      else {if(fLowQPairSwitch_E0E1[i]->At(j)=='0') continue;}
	      if((fEvt+en2)->fTracks[j].fPt < fMinPt) continue; 
	      if((fEvt+en2)->fTracks[j].fPt > fMaxPt) continue;
	      
	      pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	      ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	      qinv12 = GetQinv(pVect1, pVect2);
	      kT12 = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;
	      GetQosl(pVect1, pVect2, qout12, qside12, qlong12);

	      SetFillBins2(ch1, ch2, bin1, bin2);
	      Int_t kTindex=0;
	      if(kT12<=0.3) kTindex=0;
	      else kTindex=1;
	      
	      FSICorr12 = FSICorrelation(ch1,ch2, qinv12);
	      
	      // two particle terms filled during tabulation of low-q pairs
	      
	      
	      if(fMCcase){
		FilledMCpair12=kFALSE;

		if(ch1==ch2 && fMbin==0 && qinv12<0.2 && ENsum!=2 && ENsum!=3 && ENsum!=6){
		  for(Int_t rstep=0; rstep<10; rstep++){
		    Float_t coeff = (rstep)*0.2*(0.18/1.2);
		    Float_t phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
		    if(phi1 > 2*PI) phi1 -= 2*PI;
		    if(phi1 < 0) phi1 += 2*PI;
		    Float_t phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
		    if(phi2 > 2*PI) phi2 -= 2*PI;
		    if(phi2 < 0) phi2 += 2*PI;
		    Float_t deltaphi = phi1 - phi2;
		    if(deltaphi > PI) deltaphi -= PI;
		    if(deltaphi < -PI) deltaphi += PI;
		    
		    if(ENsum==0) ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiNum"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
		    else ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiDen"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
		  }
		  
		}// pair selection

		// Check that label does not exceed stack size
		if((fEvt)->fTracks[i].fLabel < (fEvt)->fMCarraySize && (fEvt+en2)->fTracks[j].fLabel < (fEvt+en2)->fMCarraySize){
		  if(ENsum==0 && abs((fEvt+en2)->fTracks[j].fLabel) == abs((fEvt)->fTracks[i].fLabel)) continue;
		  pVect1MC[0]=sqrt(pow((fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		  pVect2MC[0]=sqrt(pow((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPtot,2)+pow(fTrueMassPi,2));
		  pVect1MC[1]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPx; pVect2MC[1]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPx;
		  pVect1MC[2]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPy; pVect2MC[2]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPy;
		  pVect1MC[3]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPz; pVect2MC[3]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPz;
		  qinv12MC = GetQinv(pVect1MC, pVect2MC);
		  Int_t chGroup2[2]={ch1,ch2};

		  if(fGenerateSignal && (ENsum==0 || ENsum==6)){
		    if(ENsum==0) {
		      Float_t WInput = MCWeight(chGroup2, fRMax, ffcSqMRC, qinv12MC, 0.);
		      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[0].fTerms2->Fill(kT12, qinv12, WInput);
		    }else{
		      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fTerms2->Fill(kT12, qinv12);
		    }		  
		  }
		  
		  if(qinv12<0.1 && ch1==ch2 && ENsum==0) {
		    ((TProfile*)fOutputList->FindObject("fQsmearMean"))->Fill(1.,qinv12-qinv12MC); 
		    ((TProfile*)fOutputList->FindObject("fQsmearSq"))->Fill(1.,1000.*pow(qinv12-qinv12MC,2));
		    ((TH2D*)fOutputList->FindObject("fQ2Res"))->Fill(kT12, qinv12-qinv12MC);
		  }
		  		  
		  // secondary contamination
		  if(ENsum==0){
		    mcParticle1 = (AliAODMCParticle*)mcArray->At(abs((fEvt)->fTracks[i].fLabel));
		    mcParticle2 = (AliAODMCParticle*)mcArray->At(abs((fEvt+en2)->fTracks[j].fLabel));
		    if(!mcParticle1 || !mcParticle2) continue;
		    if(abs(mcParticle1->GetPdgCode())==211 && abs(mcParticle2->GetPdgCode())==211){
		      if(ch1==ch2) {
			((TH3D*)fOutputList->FindObject("fAllSCPionPairs"))->Fill(fMbin+1, kT12, qinv12);
			if(!mcParticle1->IsSecondaryFromWeakDecay() && !mcParticle2->IsSecondaryFromWeakDecay()) {
			  ((TH3D*)fOutputList->FindObject("fPrimarySCPionPairs"))->Fill(fMbin+1, kT12, qinv12);
			}	      
		      }else{
			((TH3D*)fOutputList->FindObject("fAllMCPionPairs"))->Fill(fMbin+1, kT12, qinv12);
			if(!mcParticle1->IsSecondaryFromWeakDecay() && !mcParticle2->IsSecondaryFromWeakDecay()) {
			  ((TH3D*)fOutputList->FindObject("fPrimaryMCPionPairs"))->Fill(fMbin+1, kT12, qinv12);
			}
		      }
		    }
		  }
		  
		  if(ENsum==6){// all mixed events
		  
		    Float_t rForQW=5.0;
		    if(fFSIindex<=1) rForQW=10;
		    else if(fFSIindex==2) rForQW=9;
		    else if(fFSIindex==3) rForQW=8;
		    else if(fFSIindex==4) rForQW=7;
		    else if(fFSIindex==5) rForQW=6;
		    else if(fFSIindex==6) rForQW=5;
		    else if(fFSIindex==7) rForQW=4;
		    else if(fFSIindex==8) rForQW=3;
		    else rForQW=2;
		    
		    
		    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fMCqinv->Fill(qinv12MC, MCWeight(chGroup2, rForQW, ffcSqMRC, qinv12MC, 0.));// was 4,5
		    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fMCqinvQW->Fill(qinv12MC, qinv12MC*MCWeight(chGroup2, rForQW, ffcSqMRC, qinv12MC, 0.));// was 4,5
		    // pion purity
		    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fPIDpurityDen->Fill(kT12, qinv12);
		    Int_t SCNumber = 1;
		    Int_t PdgCodeSum = abs((fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPdgCode) + abs((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPdgCode);
		    if(PdgCodeSum==22) SCNumber=1;// e-e
		    else if(PdgCodeSum==24) SCNumber=2;// e-mu
		    else if(PdgCodeSum==222) SCNumber=3;// e-pi
		    else if(PdgCodeSum==332) SCNumber=4;// e-k
		    else if(PdgCodeSum==2223) SCNumber=5;// e-p
		    else if(PdgCodeSum==26) SCNumber=6;// mu-mu
		    else if(PdgCodeSum==224) SCNumber=7;// mu-pi
		    else if(PdgCodeSum==334) SCNumber=8;// mu-k
		    else if(PdgCodeSum==2225) SCNumber=9;// mu-p
		    else if(PdgCodeSum==422) SCNumber=10;// pi-pi
		    else if(PdgCodeSum==532) SCNumber=11;// pi-k
		    else if(PdgCodeSum==2423) SCNumber=12;// pi-p
		    else if(PdgCodeSum==642) SCNumber=13;// k-k
		    else if(PdgCodeSum==2533) SCNumber=14;// k-p
		    else if(PdgCodeSum==4424) SCNumber=15;// p-p
		    else {SCNumber=16;}
		    
		    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fPIDpurityNum->Fill(SCNumber, kT12, qinv12);
		    
		    ///////////////////////
		    // muon contamination
		    Pparent1[0]=pVect1MC[0]; Pparent1[1]=pVect1MC[1]; Pparent1[2]=pVect1MC[2]; Pparent1[3]=pVect1MC[3];
		    Pparent2[0]=pVect2MC[0]; Pparent2[1]=pVect2MC[1]; Pparent2[2]=pVect2MC[2]; Pparent2[3]=pVect2MC[3];
		    pionParent1=kFALSE; pionParent2=kFALSE;
		    FilledMCpair12=kTRUE;
		    //
		    if(abs((fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPdgCode)==13 || abs((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPdgCode)==13){// muon check
		      Int_t MotherLabel1 = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fMotherLabel;
		      if(abs((fEvt)->fMCtracks[MotherLabel1].fPdgCode)==211) {
			pionParent1=kTRUE;
			Pparent1[1] = (fEvt)->fMCtracks[MotherLabel1].fPx; Pparent1[2] = (fEvt)->fMCtracks[MotherLabel1].fPy; Pparent1[3] = (fEvt)->fMCtracks[MotherLabel1].fPz;
			Pparent1[0] = sqrt(pow(Pparent1[1],2)+pow(Pparent1[2],2)+pow(Pparent1[3],2)+pow(fTrueMassPi,2));
		      }
		      // 
		      if(abs((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPdgCode)==13) {
			Int_t MotherLabel2 = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fMotherLabel;
			if(abs((fEvt+en2)->fMCtracks[MotherLabel2].fPdgCode)==211) {
			  pionParent2=kTRUE;
			  Pparent2[1] = (fEvt+en2)->fMCtracks[MotherLabel2].fPx; Pparent2[2] = (fEvt+en2)->fMCtracks[MotherLabel2].fPy; Pparent2[3] = (fEvt+en2)->fMCtracks[MotherLabel2].fPz;
			  Pparent2[0] = sqrt(pow(Pparent2[1],2)+pow(Pparent2[2],2)+pow(Pparent2[3],2)+pow(fTrueMassPi,2));
			}
		      }
		      
		      parentQinv12 = GetQinv(Pparent1, Pparent2);
		      
		      if(pionParent1 || pionParent2){
			if(parentQinv12 > 0.001 && parentQinv12 < 0.3){
			  Float_t muonPionK12 = FSICorrelation(ch1, ch2, qinv12MC);
			  Float_t pionPionK12 = FSICorrelation(ch1, ch2, parentQinv12);
			  for(Int_t term=1; term<=2; term++){
			    for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			      Float_t Rvalue = fRstartMC+Riter;
			      Float_t WInput = 1.0;
			      if(term==1) {
				WInput = MCWeight(chGroup2, Rvalue, 1.0, parentQinv12, 0.);
			      }else{
				muonPionK12 = 1.0; pionPionK12=1.0;
			      }
			      
			      Charge1[bin1].Charge2[bin2].MB[0].EDB[0].TwoPT[term-1].fMuonSmeared->Fill(Rvalue, qinv12MC, WInput);
			      Charge1[bin1].Charge2[bin2].MB[0].EDB[0].TwoPT[term-1].fMuonIdeal->Fill(Rvalue, parentQinv12, WInput);
			      Charge1[bin1].Charge2[bin2].MB[0].EDB[0].TwoPT[term-1].fMuonPionK2->Fill(Rvalue, qinv12MC, muonPionK12);
			      Charge1[bin1].Charge2[bin2].MB[0].EDB[0].TwoPT[term-1].fPionPionK2->Fill(Rvalue, parentQinv12, pionPionK12);
			    }// Riter
			  }// term loop
			  
			  if(ch1==ch2) ((TH3D*)fOutputList->FindObject("fMuonPionDeltaQinv"))->Fill(0., kT12, qinv12MC-parentQinv12);
			  else ((TH3D*)fOutputList->FindObject("fMuonPionDeltaQinv"))->Fill(1., kT12, qinv12MC-parentQinv12);
			}// parentQ check
		      }// pion parent check
		    }// muon check
		  
		    
		    Int_t indexq2 = qinv12 / 0.005;
		    if(indexq2 >=200) indexq2=199; 
		    Float_t WSpectrum = 1.0;
		    	    
		    // momentum resolution
		    for(Int_t Riter=0; Riter<fRVALUES; Riter++){
		      Float_t Rvalue = fRstartMC+Riter;
		      Float_t WInput = MCWeight(chGroup2, Rvalue, ffcSqMRC, qinv12MC, 0.);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[0].fIdeal->Fill(Rvalue, qinv12MC, WInput * WSpectrum);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[1].fIdeal->Fill(Rvalue, qinv12MC, WSpectrum);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[0].fSmeared->Fill(Rvalue, qinv12, WInput * WSpectrum);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[1].fSmeared->Fill(Rvalue, qinv12, WSpectrum);
		    }
		    
		  }// ENsum check
		}// MC array check
	      }// MC case
	      
	     
	     
	      /////////////////////////////////////////////////////////////
	      for (Int_t k=j+1; k<(fEvt+en3)->fNtracks; k++) {// 3rd particle
		if(en3==0) {
		  if(fLowQPairSwitch_E0E0[i]->At(k)=='0') continue;
		  if(fLowQPairSwitch_E0E0[j]->At(k)=='0') continue;
		}else if(en3==1){
		  if(fLowQPairSwitch_E0E1[i]->At(k)=='0') continue;
		  if(fLowQPairSwitch_E0E1[j]->At(k)=='0') continue;
		}else{
		  if(fLowQPairSwitch_E0E2[i]->At(k)=='0') continue;
		  if(fLowQPairSwitch_E1E2[j]->At(k)=='0') continue;
		}
		if((fEvt+en3)->fTracks[k].fPt < fMinPt) continue; 
		if((fEvt+en3)->fTracks[k].fPt > fMaxPt) continue;

		pVect3[0]=(fEvt+en3)->fTracks[k].fEaccepted;
		pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
		pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
		pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
		ch3 = Int_t(((fEvt+en3)->fTracks[k].fCharge + 1)/2.);
		qinv13 = GetQinv(pVect1, pVect3);
		qinv23 = GetQinv(pVect2, pVect3);
		q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));
		GetQosl(pVect1, pVect3, qout13, qside13, qlong13);
		GetQosl(pVect2, pVect3, qout23, qside23, qlong23);
		Int_t chGroup3[3]={ch1,ch2,ch3};
		Float_t QinvMCGroup3[3]={0};
		Float_t kTGroup3[3]={0};
		FilledMCtriplet123 = kFALSE;
		if(fMCcase){
		  if((fEvt+en3)->fTracks[k].fLabel == (fEvt+en2)->fTracks[j].fLabel) continue;
		  if((fEvt+en3)->fTracks[k].fLabel == (fEvt)->fTracks[i].fLabel) continue;
		  
		  pVect3MC[0]=sqrt(pow((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		  pVect3MC[1]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPx;
		  pVect3MC[2]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPy;
		  pVect3MC[3]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPz;
		  qinv13MC = GetQinv(pVect1MC, pVect3MC);
		  qinv23MC = GetQinv(pVect2MC, pVect3MC);
		  QinvMCGroup3[0] = qinv12MC; QinvMCGroup3[1] = qinv13MC; QinvMCGroup3[2] = qinv23MC;
		}
		
		
		Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
		SetFillBins3(ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);
		
		Float_t KT3 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
		if(!fq2Binning && !fLowMultBinning){
		  if(KT3<=fKT3transition) EDindex3=0;
		  else EDindex3=1;
		}else{
		  EDindex3 = fEDbin;
		  if(KT3>fKT3transition) {
		    EDindex3=2+fEDbin;
		  }
		}

		FSICorr13 = FSICorrelation(ch1,ch3, qinv13);
		FSICorr23 = FSICorrelation(ch2,ch3, qinv23);
		
		if(!fGenerateSignal && !fMCcase) {
		  momBin12 = fMomResC2SC->GetYaxis()->FindBin(qinv12);
		  momBin13 = fMomResC2SC->GetYaxis()->FindBin(qinv13);
		  momBin23 = fMomResC2SC->GetYaxis()->FindBin(qinv23);
		  if(qinv12 > 0.1) momBin12 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		  if(qinv13 > 0.1) momBin13 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		  if(qinv23 > 0.1) momBin23 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		  //
		  if(ch1==ch2) MomResCorr12 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin12);
		  else MomResCorr12 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin12);
		  if(ch1==ch3) MomResCorr13 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin13);
		  else MomResCorr13 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin13);
		  if(ch2==ch3) MomResCorr23 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin23);
		  else MomResCorr23 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin23);
		}
		if(ENsum==0) {
		  Float_t Winput=1.0;
		  if(fMCcase && fGenerateSignal) Winput = MCWeight3(1, fRMax, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fTerms3->Fill(q3, Winput); 
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fKfactor->Fill(q3, 1/(FSICorr12*FSICorr13*FSICorr23), Winput);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fKfactorWeighted->Fill(q3, 1/(FSICorr12*FSICorr13*FSICorr23), MomResCorr12*MomResCorr13*MomResCorr23 * Winput);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fMeanQinv->Fill(q3, qinv12);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fMeanQinv->Fill(q3, qinv13);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[0].fMeanQinv->Fill(q3, qinv23);
		  if(bin1==bin2 && bin1==bin3){
		    Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[0].fTerms33D->Fill(qinv12, qinv13, qinv23);
		    Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[0].fKfactor3D->Fill(qinv12, qinv13, qinv23, 1/(FSICorr12*FSICorr13*FSICorr23));
		  }
		}
		if(ENsum==6) {
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fTerms3->Fill(q3);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fMeanQinv->Fill(q3, qinv12);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fMeanQinv->Fill(q3, qinv13);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fMeanQinv->Fill(q3, qinv23);
		  if(bin1==bin2 && bin1==bin3) Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[4].fTerms33D->Fill(qinv12, qinv13, qinv23);
		}
		if(ENsum==3){
		  Float_t Winput=1.0;
		  if(fill2) {
		    if(fMCcase && fGenerateSignal) Winput = MCWeight3(2, fRMax, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fTerms3->Fill(q3, Winput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fKfactor->Fill(q3, 1/(FSICorr12), Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12 * Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fMeanQinv->Fill(q3, qinv12);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fMeanQinv->Fill(q3, qinv13);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[1].fMeanQinv->Fill(q3, qinv23);
		    if(bin1==bin2 && bin1==bin3){
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[1].fTerms33D->Fill(qinv12, qinv13, qinv23);
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[1].fKfactor3D->Fill(qinv12, qinv13, qinv23, 1/(FSICorr12));
		    }
		  }if(fill3) {
		    if(fMCcase && fGenerateSignal) Winput = MCWeight3(3, fRMax, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fTerms3->Fill(q3, Winput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fKfactor->Fill(q3, 1/(FSICorr12), Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12 * Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fMeanQinv->Fill(q3, qinv12);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fMeanQinv->Fill(q3, qinv13);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[2].fMeanQinv->Fill(q3, qinv23);
		    if(bin1==bin2 && bin1==bin3){
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[2].fTerms33D->Fill(qinv13, qinv12, qinv23);
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[2].fKfactor3D->Fill(qinv13, qinv12, qinv23, 1/(FSICorr12));
		    }
		  }if(fill4) {
		    if(fMCcase && fGenerateSignal) Winput = MCWeight3(4, fRMax, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fTerms3->Fill(q3, Winput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fKfactor->Fill(q3, 1/(FSICorr12), Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12 * Winput);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fMeanQinv->Fill(q3, qinv12);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fMeanQinv->Fill(q3, qinv13);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[3].fMeanQinv->Fill(q3, qinv23);
		    if(bin1==bin2 && bin1==bin3){
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[3].fTerms33D->Fill(qinv13, qinv23, qinv12);
		      Charge1[0].Charge2[0].Charge3[0].MB[fMbin].EDB[EDindex3].ThreePT[3].fKfactor3D->Fill(qinv13, qinv23, qinv12, 1/(FSICorr12));
		    }
		  }
		}
		
		// C3 Building
		if(ENsum==6 && ch1==ch2 && ch1==ch3){
		  Positive1stTripletWeights = kTRUE;
		  //
		  GetWeight(pVect1, pVect2, weight12, weight12Err);
		  GetWeight(pVect1, pVect3, weight13, weight13Err);
		  GetWeight(pVect2, pVect3, weight23, weight23Err);
		  
		  
		  if(sqrt(fabs(weight12*weight13*weight23)) > 1.0) {// weight should never be larger than 1
		    if(fMbin==0 && bin1==0 && EDindex3==0) {
		      ((TH1D*)fOutputList->FindObject("fTPNRejects3pion1"))->Fill(q3, sqrt(fabs(weight12*weight13*weight23)));
		    }
		  }
		  
		  if(fOnlineCorrection){
		    Float_t MuonCorr12=1.0, MuonCorr13=1.0, MuonCorr23=1.0;
		    if(!fGenerateSignal && !fMCcase) {
		      MuonCorr12 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin12);
		      MuonCorr13 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin13);
		      MuonCorr23 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin23);
		    }
		    
		    // no MRC, no Muon Correction
		    weight12CC[0] = ((weight12+1) - ffcSq*FSICorr12 - (1-ffcSq));
		    weight12CC[0] /= FSICorr12*ffcSq;
		    weight13CC[0] = ((weight13+1) - ffcSq*FSICorr13 - (1-ffcSq));
		    weight13CC[0] /= FSICorr13*ffcSq;
		    weight23CC[0] = ((weight23+1) - ffcSq*FSICorr23 - (1-ffcSq));
		    weight23CC[0] /= FSICorr23*ffcSq;
		    if(weight12CC[0] > 0 && weight13CC[0] > 0 && weight23CC[0] > 0){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuild->Fill(1, q3, sqrt(weight12CC[0]*weight13CC[0]*weight23CC[0]));
		    }
		    // no Muon Correction
		    weight12CC[1] = ((weight12+1)*MomResCorr12 - ffcSq*FSICorr12 - (1-ffcSq));
		    weight12CC[1] /= FSICorr12*ffcSq;
		    weight13CC[1] = ((weight13+1)*MomResCorr13 - ffcSq*FSICorr13 - (1-ffcSq));
		    weight13CC[1] /= FSICorr13*ffcSq;
		    weight23CC[1] = ((weight23+1)*MomResCorr23 - ffcSq*FSICorr23 - (1-ffcSq));
		    weight23CC[1] /= FSICorr23*ffcSq;
		    if(weight12CC[1] > 0 && weight13CC[1] > 0 && weight23CC[1] > 0){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuild->Fill(2, q3, sqrt(weight12CC[1]*weight13CC[1]*weight23CC[1]));
		    }
		    // both Corrections
		    weight12CC[2] = ((weight12+1)*MomResCorr12 - ffcSq*FSICorr12 - (1-ffcSq));
		    weight12CC[2] /= FSICorr12*ffcSq;
		    weight12CC[2] *= MuonCorr12;
		    weight13CC[2] = ((weight13+1)*MomResCorr13 - ffcSq*FSICorr13 - (1-ffcSq));
		    weight13CC[2] /= FSICorr13*ffcSq;
		    weight13CC[2] *= MuonCorr13;
		    weight23CC[2] = ((weight23+1)*MomResCorr23 - ffcSq*FSICorr23 - (1-ffcSq));
		    weight23CC[2] /= FSICorr23*ffcSq;
		    weight23CC[2] *= MuonCorr23;
		    //
		  }else{
		    weight12CC[2] = weight12;
		    weight13CC[2] = weight13;
		    weight23CC[2] = weight23;
		  }
		  
		  
		  if(weight12CC[2] < 0 || weight13CC[2] < 0 || weight23CC[2] < 0) {// C2^QS can never be less than unity
		    if(fMbin==0 && bin1==0 && EDindex3==0) {
		      ((TH1D*)fOutputList->FindObject("fTPNRejects3pion2"))->Fill(q3, sqrt(fabs(weight12CC[2]*weight13CC[2]*weight23CC[2])));
		    }
		    if(weight12CC[2] < 0) weight12CC[2]=0;
		    if(weight13CC[2] < 0) weight13CC[2]=0;
		    if(weight23CC[2] < 0) weight23CC[2]=0;
		    Positive1stTripletWeights = kFALSE;
		  }
		  /////////////////////////////////////////////////////
		  weightTotal = sqrt(weight12CC[2]*weight13CC[2]*weight23CC[2]);
		  /////////////////////////////////////////////////////
		  if(Positive1stTripletWeights){
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuild->Fill(3, q3, weightTotal);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuild->Fill(4, q3, 1);
		  }else{
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuildNeg->Fill(4, q3, 1);
		  }
		
		  //
		  // Full Weight reconstruction
		  
		  for(Int_t RcohIndex=0; RcohIndex<kRcohsteps; RcohIndex++){// Rcoh 
		    t12 = exp(-pow(RcohIndex/FmToGeV * qinv12,2)/2.);
		    t23 = exp(-pow(RcohIndex/FmToGeV * qinv23,2)/2.);
		    t13 = exp(-pow(RcohIndex/FmToGeV * qinv13,2)/2.);
		    for(Int_t GIndex=0; GIndex<kGsteps; GIndex++){
		      Int_t FillBin = 5 + RcohIndex*kGsteps + GIndex;
		      Float_t G = 0.02*GIndex;
		      if(RcohIndex!=6){
			T12 = (-2*G*(1-G)*t12 + sqrt(pow(2*G*(1-G)*t12,2) + 4*pow(1-G,2)*weight12CC[2])) / (2*pow(1-G,2));
			T13 = (-2*G*(1-G)*t13 + sqrt(pow(2*G*(1-G)*t13,2) + 4*pow(1-G,2)*weight13CC[2])) / (2*pow(1-G,2));
			T23 = (-2*G*(1-G)*t23 + sqrt(pow(2*G*(1-G)*t23,2) + 4*pow(1-G,2)*weight23CC[2])) / (2*pow(1-G,2));
		      }else{// Full Size
			T12 = sqrt(weight12CC[2] / (1-G*G));
			T13 = sqrt(weight13CC[2] / (1-G*G));
			T23 = sqrt(weight23CC[2] / (1-G*G));
			t12 = T12;
			t13 = T13;
			t23 = T23;
		      }
		      weightTotal = 2*G*(1-G)*(T12*t12 + T13*t13 + T23*t23) + pow(1-G,2)*(T12*T12 + T13*T13 + T23*T23);
		      weightTotal += 2*G*pow(1-G,2)*(T12*T13*t23 + T12*T23*t13 + T13*T23*t12) + 2*pow(1-G,3)*T12*T13*T23;
		      weightCumulant = 2*G*pow(1-G,2)*(T12*T13*t23 + T12*T23*t13 + T13*T23*t12) + 2*pow(1-G,3)*T12*T13*T23;
		      
		      if(Positive1stTripletWeights){
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuild->Fill(FillBin, q3, weightTotal);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fCumulantBuild->Fill(FillBin, q3, weightCumulant);
		      }else{
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuildNeg->Fill(FillBin, q3, weightTotal);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fCumulantBuildNeg->Fill(FillBin, q3, weightCumulant);
		      }
		    }
		  }
		  //
		  //weight12CC_e = weight12Err*MomResCorr12 / FSICorr12 / ffcSq * MuonCorr12;
		  //weight13CC_e = weight13Err*MomResCorr13 / FSICorr13 / ffcSq * MuonCorr13;
		  //weight23CC_e = weight23Err*MomResCorr23 / FSICorr23 / ffcSq * MuonCorr23;
		  //if(weight12CC[2]*weight13CC[2]*weight23CC[2] > 0){
		  //weightTotalErr = pow(2 * sqrt(3) * weight12CC_e*weight13CC[2]*weight23CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight23CC[2]),2);
		  //}
		  //weightTotalErr += pow(weight12CC_e,2) + pow(weight13CC_e,2) + pow(weight23CC_e,2);
		  //Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[EDindex3].ThreePT[4].fBuildErr->Fill(4, q3, weightTotalErr);
		  
		  //}// 1st r3 den check
		  
		}// C3 Building section end
		
		
		/*if(ch1==ch2 && ch1==ch3){
		   Float_t pt1=sqrt(pow(pVect1[1],2)+pow(pVect1[2],2));
		   Float_t pt2=sqrt(pow(pVect2[1],2)+pow(pVect2[2],2));
		   Float_t pt3=sqrt(pow(pVect3[1],2)+pow(pVect3[2],2));
		  if(ENsum==0){
		    ((TH3D*)fOutputList->FindObject("fKT3DistTerm1"))->Fill(fMbin+1, KT3, q3);
		    if(q3<0.1){
		      ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, EDindex3, pt1);
		      ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, EDindex3, pt2);
		      ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, EDindex3, pt3);
		    }
		  }
		  if(fMbin==0){
		    if(ENsum==0){
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum0"))->Fill(EDindex3, q3, pt1);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum0"))->Fill(EDindex3, q3, pt2);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum0"))->Fill(EDindex3, q3, pt3);
		    }
		    if(ENsum==3){
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum3"))->Fill(EDindex3, q3, pt1);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum3"))->Fill(EDindex3, q3, pt2);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum3"))->Fill(EDindex3, q3, pt3);
		    }
		    if(ENsum==6){
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum6"))->Fill(EDindex3, q3, pt1);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum6"))->Fill(EDindex3, q3, pt2);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpTENsum6"))->Fill(EDindex3, q3, pt3);
		    }
		  }
		  
		  }
		if(ch1==ch2 && ch1==ch3 && ENsum==6) ((TH3D*)fOutputList->FindObject("fKT3DistTerm5"))->Fill(fMbin+1, KT3, q3);
	      	*/
		

		
		if(fMCcase && ENsum==6 && FilledMCpair12){// for momentum resolution and muon correction
		  if((fEvt+en3)->fTracks[k].fLabel < (fEvt+en3)->fMCarraySize){
		    
		    pVect3MC[0]=sqrt(pow((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		    pVect3MC[1]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPx;
		    pVect3MC[2]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPy;
		    pVect3MC[3]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPz;
		    qinv13MC = GetQinv(pVect1MC, pVect3MC);
		    qinv23MC = GetQinv(pVect2MC, pVect3MC);
		    
		    q3MC = sqrt(pow(qinv12MC,2)+pow(qinv13MC,2)+pow(qinv23MC,2));
		    if(q3<0.1 && ch1==ch2 && ch1==ch3) ((TH2D*)fOutputList->FindObject("fQ3Res"))->Fill(KT3, q3-q3MC);
		    
		    //Float_t TripletWeightTTC=1.0;// same-charge weights to mimic two-track depletion of same-charge pairs
		    //if(ch1==ch2 && qinv12>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv12);
		    //if(ch1==ch3 && qinv13>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv13);
		    //if(ch2==ch3 && qinv23>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv23);
		    
		    
		    Pparent3[0]=pVect3MC[0]; Pparent3[1]=pVect3MC[1]; Pparent3[2]=pVect3MC[2]; Pparent3[3]=pVect3MC[3];
		    pionParent3=kFALSE;
		    
		    if(abs((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPdgCode)==13){// muon check
		      Int_t MotherLabel3 = (fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fMotherLabel;
		      if(abs((fEvt+en3)->fMCtracks[MotherLabel3].fPdgCode)==211) {
			pionParent3=kTRUE;
			Pparent3[1] = (fEvt+en3)->fMCtracks[MotherLabel3].fPx; Pparent3[2] = (fEvt+en3)->fMCtracks[MotherLabel3].fPy; Pparent3[3] = (fEvt+en3)->fMCtracks[MotherLabel3].fPz;
			Pparent3[0] = sqrt(pow(Pparent3[1],2)+pow(Pparent3[2],2)+pow(Pparent3[3],2)+pow(fTrueMassPi,2));
		      }
		    }
		    
		    parentQinv13 = GetQinv(Pparent1, Pparent3);
		    parentQinv23 = GetQinv(Pparent2, Pparent3);
		    parentQ3 = sqrt(pow(parentQinv12,2) + pow(parentQinv13,2) + pow(parentQinv23,2));
		    
		    if(parentQinv12 > 0.001 && parentQinv13 > 0.001 && parentQinv23 > 0.001 && parentQ3 < 0.5){
		      FilledMCtriplet123=kTRUE;
		      if(pionParent1 || pionParent2 || pionParent3) {// want at least one pion-->muon
			
			Float_t parentQinvGroup3[3]={parentQinv12, parentQinv13, parentQinv23};
			//Float_t parentkTGroup3[3]={float(sqrt(pow(Pparent1[1]+Pparent2[1],2) + pow(Pparent1[2]+Pparent2[2],2))/2.),
			//float(sqrt(pow(Pparent1[1]+Pparent3[1],2) + pow(Pparent1[2]+Pparent3[2],2))/2.),
			//float(sqrt(pow(Pparent2[1]+Pparent3[1],2) + pow(Pparent2[2]+Pparent3[2],2))/2.)};
			Float_t parentkTGroup3[3]={0};
			
			((TH2D*)fOutputList->FindObject("fAvgQ12VersusQ3"))->Fill(parentQ3, parentQinv12);
			((TH2D*)fOutputList->FindObject("fAvgQ13VersusQ3"))->Fill(parentQ3, parentQinv13);
			((TH2D*)fOutputList->FindObject("fAvgQ23VersusQ3"))->Fill(parentQ3, parentQinv23);
			
			if(q3MC>=sqrt(3.)*fQLowerCut){
			  for(Int_t term=1; term<=4; term++){
			    if(term==1) {}
			    else if(term==2) {if(!pionParent1 && !pionParent2) continue;}
			    else if(term==3) {if(!pionParent1 && !pionParent3) continue;}
			    else {if(!pionParent2 && !pionParent3) continue;}
			    for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			      Float_t Rvalue = fRstartMC+Riter;
			      Float_t WInput = MCWeight3(term, Rvalue, 1.0, chGroup3, parentQinvGroup3, parentkTGroup3);
			      Float_t WInputParentFSI = MCWeightFSI3(term, Rvalue, 1.0, chGroup3, parentQinvGroup3);
			      Float_t WInputFSI = MCWeightFSI3(term, Rvalue, 1.0, chGroup3, QinvMCGroup3);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonSmeared->Fill(1, Rvalue, q3MC, WInput);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonIdeal->Fill(1, Rvalue, parentQ3, WInput);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonPionK3->Fill(1, Rvalue, q3MC, WInputFSI);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fPionPionK3->Fill(1, Rvalue, parentQ3, WInputParentFSI);
			      //
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonSmeared->Fill(2, Rvalue, q3MC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonIdeal->Fill(2, Rvalue, parentQ3);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonPionK3->Fill(2, Rvalue, q3MC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fPionPionK3->Fill(2, Rvalue, parentQ3);
			    }// Riter
			  }// term loop
			}// q3MC min
		      }// pion parent check
		    }// parentQ check (muon correction)
		      
		    
		    Int_t indexq3 = q3 / 0.005;
		    if(indexq3 >=35) indexq3=34; 
		    Float_t WSpectrum = 1;

		    // 3-pion momentum resolution
		    for(Int_t term=1; term<=5; term++){
		      for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			Float_t Rvalue = fRstartMC+Riter;
			Float_t WInput = MCWeight3(term, Rvalue, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
			if(q3MC>=sqrt(3.)*fQLowerCut) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[EDindex3].ThreePT[term-1].fIdeal->Fill(Rvalue, q3MC, WInput*WSpectrum);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[EDindex3].ThreePT[term-1].fSmeared->Fill(Rvalue, q3, WInput*WSpectrum);
		      }
		    }
		  }// 3rd particle label check
		}// MCcase and ENsum==6
		
		
		
		
		/////////////////////////////////////////////////////////////
		for (Int_t l=k+1; l<(fEvt+en4)->fNtracks; l++) {// 4th particle
		  if(en4==0){
		    if(fLowQPairSwitch_E0E0[i]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E0E0[j]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E0E0[k]->At(l)=='0') continue;
		  }else if(en4==1){
		    if(en3==0){
		      if(fLowQPairSwitch_E0E1[i]->At(l)=='0') continue;
		      if(fLowQPairSwitch_E0E1[j]->At(l)=='0') continue;
		      if(fLowQPairSwitch_E0E1[k]->At(l)=='0') continue;
		    }else{ 
		      if(fLowQPairSwitch_E0E1[i]->At(l)=='0') continue;
		      if(fLowQPairSwitch_E0E1[j]->At(l)=='0') continue;
		      if(fLowQPairSwitch_E1E1[k]->At(l)=='0') continue;
		    }
		  }else if(en4==2){
		    if(fLowQPairSwitch_E0E2[i]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E0E2[j]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E1E2[k]->At(l)=='0') continue;
		  }else{
		    if(fLowQPairSwitch_E0E3[i]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E1E3[j]->At(l)=='0') continue;
		    if(fLowQPairSwitch_E2E3[k]->At(l)=='0') continue;
		  }
		  if((fEvt+en4)->fTracks[l].fPt < fMinPt) continue; 
		  if((fEvt+en4)->fTracks[l].fPt > fMaxPt) continue;
		  
		  pVect4[0]=(fEvt+en4)->fTracks[l].fEaccepted;
		  pVect4[1]=(fEvt+en4)->fTracks[l].fP[0];
		  pVect4[2]=(fEvt+en4)->fTracks[l].fP[1];
		  pVect4[3]=(fEvt+en4)->fTracks[l].fP[2];
		  ch4 = Int_t(((fEvt+en4)->fTracks[l].fCharge + 1)/2.);
		  qinv14 = GetQinv(pVect1, pVect4);
		  qinv24 = GetQinv(pVect2, pVect4);
		  qinv34 = GetQinv(pVect3, pVect4);
		  q4 = sqrt(pow(q3,2) + pow(qinv14,2) + pow(qinv24,2) + pow(qinv34,2));
		  GetQosl(pVect1, pVect4, qout14, qside14, qlong14);
		  GetQosl(pVect2, pVect4, qout24, qside24, qlong24);
		  GetQosl(pVect3, pVect4, qout34, qside34, qlong34);
		  Int_t chGroup4[4]={ch1,ch2,ch3,ch4};
		  Float_t QinvMCGroup4[6]={0};
		  Float_t kTGroup4[6]={0};
		  
		  
		  if(fMCcase){// for momentum resolution and muon correction
		    if((fEvt+en4)->fTracks[l].fLabel == (fEvt+en3)->fTracks[k].fLabel) continue;
		    if((fEvt+en4)->fTracks[l].fLabel == (fEvt+en2)->fTracks[j].fLabel) continue;
		    if((fEvt+en4)->fTracks[l].fLabel == (fEvt)->fTracks[i].fLabel) continue;
		    
		    pVect4MC[0]=sqrt(pow((fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		    pVect4MC[1]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPx;
		    pVect4MC[2]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPy;
		    pVect4MC[3]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPz;
		    qinv14MC = GetQinv(pVect1MC, pVect4MC);
		    qinv24MC = GetQinv(pVect2MC, pVect4MC);
		    qinv34MC = GetQinv(pVect3MC, pVect4MC);
		    
		    QinvMCGroup4[0] = qinv12MC; QinvMCGroup4[1] = qinv13MC; QinvMCGroup4[2] = qinv14MC;
		    QinvMCGroup4[3] = qinv23MC; QinvMCGroup4[4] = qinv24MC; QinvMCGroup4[5] = qinv34MC;
		    
		  }
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6){
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(1, qinv12); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(2, qinv13);
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(3, qinv14); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(4, qinv23);
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(5, qinv24); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(6, qinv34);
		  }
		  
		  Float_t KT4 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1]+pVect4[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2]+pVect4[2],2))/4.;
		  Float_t QoutSum = fabs(qout12) + fabs(qout13) + fabs(qout14) + fabs(qout23) + fabs(qout24) + fabs(qout34);
		  Float_t QsideSum = fabs(qside12) + fabs(qside13) + fabs(qside14) + fabs(qside23) + fabs(qside24) + fabs(qside34);
		  Float_t QlongSum = fabs(qlong12) + fabs(qlong13) + fabs(qlong14) + fabs(qlong23) + fabs(qlong24) + fabs(qlong34);
		  
		  if(!fq2Binning && !fLowMultBinning && fQdirectionBinning==0){
		    if(KT4<=fKT4transition) EDindex4=0;
		    else EDindex4=1;
		  }else{
		    if(KT4<=fKT4transition){
		      if(ch1==ch2 && ch1==ch3 && ch1==ch4){
			((TH2D*)fOutputList->FindObject("fQoutSum"))->Fill(q4, QoutSum);
			((TH2D*)fOutputList->FindObject("fQsideSum"))->Fill(q4, QsideSum);
			((TH2D*)fOutputList->FindObject("fQlongSum"))->Fill(q4, QlongSum);
		      }
		    }
		    if(fQdirectionBinning==1){
		      if(QoutSum < fqOutFcn->Eval(q4)) fEDbin=0;
		      else fEDbin=1;
		    }else if(fQdirectionBinning==2){
		      if(QsideSum < fqSideFcn->Eval(q4)) fEDbin=0;
		      else fEDbin=1;
		    }else if(fQdirectionBinning==3){
		      if(QlongSum < fqLongFcn->Eval(q4)) fEDbin=0;
		      else fEDbin=1;
		    }else {}
		    //
		    EDindex4 = fEDbin;
		    if(KT4>fKT4transition) EDindex4=2+fEDbin;
		  }
		  
		  FSICorr14 = FSICorrelation(ch1,ch4, qinv14);
		  FSICorr24 = FSICorrelation(ch2,ch4, qinv24);
		  FSICorr34 = FSICorrelation(ch3,ch4, qinv34);
		  
		  if(!fGenerateSignal && !fMCcase) {
		    momBin14 = fMomResC2SC->GetYaxis()->FindBin(qinv14);
		    momBin24 = fMomResC2SC->GetYaxis()->FindBin(qinv24);
		    momBin34 = fMomResC2SC->GetYaxis()->FindBin(qinv34);		  
		    if(qinv14 > 0.1) momBin14 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		    if(qinv24 > 0.1) momBin24 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		    if(qinv34 > 0.1) momBin34 = fMomResC2SC->GetYaxis()->FindBin(0.095);
		    //
		    if(ch1==ch4) MomResCorr14 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin14);
		    else MomResCorr14 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin14);
		    if(ch2==ch4) MomResCorr24 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin24);
		    else MomResCorr24 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin24);
		    if(ch3==ch4) MomResCorr34 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin34);
		    else MomResCorr34 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin34);
		  }
		  
		  Bool_t FillTerms[13]={kFALSE};
		  SetFillBins4(ch1, ch2, ch3, ch4, bin1, bin2, bin3, bin4, ENsum, FillTerms);
		  //
		  for(int ft=0; ft<13; ft++) {
		    Float_t FSIfactor = 1.0;
		    Float_t MomResWeight = 1.0;
		    Float_t WInput = 1.0;
		    if(fMCcase && fGenerateSignal) WInput = MCWeight4(ft+1, fRMax, ffcSqMRC, chGroup4, QinvMCGroup4, kTGroup4);
		    if(ft==0) {
		      FSIfactor = 1/(FSICorr12 * FSICorr13 * FSICorr14 * FSICorr23 * FSICorr24 * FSICorr34);
		      MomResWeight = MomResCorr12 * MomResCorr13 * MomResCorr14 * MomResCorr23 * MomResCorr24 * MomResCorr34;
		    }else if(ft<=4) {
		      FSIfactor = 1/(FSICorr12 * FSICorr13 * FSICorr23);
		      MomResWeight = MomResCorr12 * MomResCorr13 * MomResCorr23;
		    }else if(ft<=10) {
		      FSIfactor = 1/(FSICorr12);
		      MomResWeight = MomResCorr12;
		    }else if(ft==11) {
		      FSIfactor = 1/(FSICorr12 * FSICorr34);
		      MomResWeight = MomResCorr12 * MomResCorr34;
		    }else {FSIfactor = 1.0; MomResWeight = 1.0;}
		    
		    if(FillTerms[ft]) {
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[ft].fTerms4->Fill(q4, WInput);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[ft].fKfactor->Fill(q4, FSIfactor, WInput);
		      //Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[ft].fKfactorWeighted->Fill(q4, FSIfactor, MomResWeight*WInput);
		    }
		  }
		  
		  /////////////////////////////////////////////////////////////
		  // C4 building
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6 && !fMCcase){
		    Positive2ndTripletWeights=kTRUE;
		    //
		    GetWeight(pVect1, pVect4, weight14, weight14Err);
		    GetWeight(pVect2, pVect4, weight24, weight24Err);
		    GetWeight(pVect3, pVect4, weight34, weight34Err);
		    
		    if(fOnlineCorrection){
		      Float_t MuonCorr14=1.0, MuonCorr24=1.0, MuonCorr34=1.0;
		      if(!fGenerateSignal && !fMCcase) {
			MuonCorr14 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin14);
			MuonCorr24 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin24);
			MuonCorr34 = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin34);
		      }
		      
		      // no MRC, no Muon Correction
		      weight14CC[0] = ((weight14+1) - ffcSq*FSICorr14 - (1-ffcSq));
		      weight14CC[0] /= FSICorr14*ffcSq;
		      weight24CC[0] = ((weight24+1) - ffcSq*FSICorr24 - (1-ffcSq));
		      weight24CC[0] /= FSICorr24*ffcSq;
		      weight34CC[0] = ((weight34+1) - ffcSq*FSICorr34 - (1-ffcSq));
		      weight34CC[0] /= FSICorr34*ffcSq;
		      if(weight14CC[0] > 0 && weight24CC[0] > 0 && weight34CC[0] > 0 && weight12CC[0] > 0 && weight13CC[0] > 0 && weight23CC[0] > 0){
			weightTotal  = sqrt(weight12CC[0]*weight13CC[0]*weight24CC[0]*weight34CC[0]);
			weightTotal += sqrt(weight12CC[0]*weight14CC[0]*weight23CC[0]*weight34CC[0]);
			weightTotal += sqrt(weight13CC[0]*weight14CC[0]*weight23CC[0]*weight24CC[0]);
			weightTotal /= 3.;
			Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuild->Fill(1, q4, weightTotal);
		      }
		      // no Muon Correction
		      weight14CC[1] = ((weight14+1)*MomResCorr14 - ffcSq*FSICorr14 - (1-ffcSq));
		      weight14CC[1] /= FSICorr14*ffcSq;
		      weight24CC[1] = ((weight24+1)*MomResCorr24 - ffcSq*FSICorr24 - (1-ffcSq));
		      weight24CC[1] /= FSICorr24*ffcSq;
		      weight34CC[1] = ((weight34+1)*MomResCorr34 - ffcSq*FSICorr34 - (1-ffcSq));
		      weight34CC[1] /= FSICorr34*ffcSq;
		      if(weight14CC[1] > 0 && weight24CC[1] > 0 && weight34CC[1] > 0 && weight12CC[1] > 0 && weight13CC[1] > 0 && weight23CC[1] > 0){
			weightTotal  = sqrt(weight12CC[1]*weight13CC[1]*weight24CC[1]*weight34CC[1]);
			weightTotal += sqrt(weight12CC[1]*weight14CC[1]*weight23CC[1]*weight34CC[1]);
			weightTotal += sqrt(weight13CC[1]*weight14CC[1]*weight23CC[1]*weight24CC[1]);
			weightTotal /= 3.;
			Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuild->Fill(2, q4, weightTotal);
		      }
		      // both corrections
		      weight14CC[2] = ((weight14+1)*MomResCorr14 - ffcSq*FSICorr14 - (1-ffcSq));
		      weight14CC[2] /= FSICorr14*ffcSq;
		      weight14CC[2] *= MuonCorr14;
		      weight24CC[2] = ((weight24+1)*MomResCorr24 - ffcSq*FSICorr24 - (1-ffcSq));
		      weight24CC[2] /= FSICorr24*ffcSq;
		      weight24CC[2] *= MuonCorr24;
		      weight34CC[2] = ((weight34+1)*MomResCorr34 - ffcSq*FSICorr34 - (1-ffcSq));
		      weight34CC[2] /= FSICorr34*ffcSq;
		      weight34CC[2] *= MuonCorr34;
		    }else{
		      weight14CC[2] = weight14;
		      weight24CC[2] = weight24;
		      weight34CC[2] = weight34;
		    }
		    
		    
		    if(weight14CC[2] < 0 || weight24CC[2] < 0 || weight34CC[2] < 0) {// C2^QS can never be less than unity
		      if(fMbin==0 && bin1==0 && EDindex4==0) {
			((TH1D*)fOutputList->FindObject("fTPNRejects4pion1"))->Fill(q4, sqrt(fabs(weight12CC[2]*weight23CC[2]*weight34CC[2]*weight14CC[2])));
		      }
		      if(weight14CC[2] < 0) weight14CC[2]=0;
		      if(weight24CC[2] < 0) weight24CC[2]=0;
		      if(weight34CC[2] < 0) weight34CC[2]=0;
		      Positive2ndTripletWeights=kFALSE;
		    }
		    /////////////////////////////////////////////////////
		    weightTotal  = sqrt(weight12CC[2]*weight13CC[2]*weight24CC[2]*weight34CC[2]);
		    weightTotal += sqrt(weight12CC[2]*weight14CC[2]*weight23CC[2]*weight34CC[2]);
		    weightTotal += sqrt(weight13CC[2]*weight14CC[2]*weight23CC[2]*weight24CC[2]);
		    weightTotal /= 3.;
		    if(Positive1stTripletWeights && Positive2ndTripletWeights){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuild->Fill(3, q4, weightTotal);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuild->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuild->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuild->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuild->Fill(4, q4, 1);
		    }else{
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildNeg->Fill(3, q4, weightTotal);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildNeg->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuildNeg->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuildNeg->Fill(4, q4, 1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuildNeg->Fill(4, q4, 1);
		    }
		   
		    
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildFromFits->Fill(1, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuildFromFits->Fill(1, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuildFromFits->Fill(1, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuildFromFits->Fill(1, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildFromFits->Fill(2, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuildFromFits->Fill(2, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuildFromFits->Fill(2, 4, q4, 1);
		    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuildFromFits->Fill(2, 4, q4, 1);
		    
		    // Full Weight reconstruction
		    for(Int_t type=0; type<3; type++){// C2 interpolation, c3 fit, C3 fit
		      for(Int_t RcohIndex=0; RcohIndex<kRcohsteps; RcohIndex++){// Rcoh=0,1,2,3,4,5 fm, then Rcoh=Rch
			if(RcohIndex>1 && RcohIndex<6) continue;// save cpu time
			if(fCollisionType!=0 && RcohIndex!=6) continue;// save cpu time
			t12 = exp(-pow(RcohIndex/FmToGeV * qinv12,2)/2.);
			t13 = exp(-pow(RcohIndex/FmToGeV * qinv13,2)/2.);
			t14 = exp(-pow(RcohIndex/FmToGeV * qinv14,2)/2.);
			t23 = exp(-pow(RcohIndex/FmToGeV * qinv23,2)/2.);
			t24 = exp(-pow(RcohIndex/FmToGeV * qinv24,2)/2.);
			t34 = exp(-pow(RcohIndex/FmToGeV * qinv34,2)/2.);
						
			for(Int_t GIndex=0; GIndex<kGsteps; GIndex++){
			  Int_t FillBin = 5 + RcohIndex*kGsteps + GIndex;
			  Float_t G = 0.02*GIndex;
			  
			  if(type==0){// From C2
			    Float_t a = pow(1-G,2);
			    Float_t b = 2*G*(1-G);
			    if(RcohIndex!=6){
			      T12 = (-b*t12 + sqrt(pow(b*t12,2) + 4*a*weight12CC[2])) / (2*a);
			      T13 = (-b*t13 + sqrt(pow(b*t13,2) + 4*a*weight13CC[2])) / (2*a);
			      T14 = (-b*t14 + sqrt(pow(b*t14,2) + 4*a*weight14CC[2])) / (2*a);
			      T23 = (-b*t23 + sqrt(pow(b*t23,2) + 4*a*weight23CC[2])) / (2*a);
			      T24 = (-b*t24 + sqrt(pow(b*t24,2) + 4*a*weight24CC[2])) / (2*a);
			      T34 = (-b*t34 + sqrt(pow(b*t34,2) + 4*a*weight34CC[2])) / (2*a);
			    }else{
			      T12 = sqrt(weight12CC[2] / (1-G*G));
			      T13 = sqrt(weight13CC[2] / (1-G*G));
			      T14 = sqrt(weight14CC[2] / (1-G*G));
			      T23 = sqrt(weight23CC[2] / (1-G*G));
			      T24 = sqrt(weight24CC[2] / (1-G*G));
			      T34 = sqrt(weight34CC[2] / (1-G*G));
			      t12 = T12;
			      t13 = T13;
			      t14 = T14;
			      t23 = T23;
			      t24 = T24;
			      t34 = T34;
			    }
			  }else{// from c3 or C3 fit
			    T12 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv12);
			    T13 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv13);
			    T14 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv14);
			    T23 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv23);
			    T24 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv24);
			    T34 = ExchangeAmp[RcohIndex][GIndex][type-1]->Eval(qinv34);
			    if(RcohIndex==6){
			      t12 = T12;
			      t13 = T13;
			      t14 = T14;
			      t23 = T23;
			      t24 = T24;
			      t34 = T34;
			    }
			  }
			  			  			  
			  // Build the correlation functions
			  weightTotal = 2*G*(1-G)*(T12*t12 + T13*t13 + T14*t14 + T23*t23 + T24*t24 + T34*t34);// 2-pion
			  weightTotal += pow(1-G,2)*(T12*T12 + T13*T13 + T14*T14 + T23*T23 + T24*T24 + T34*T34);// 2-pion fully chaotic
			  weightTotal += 2*G*pow(1-G,3)*(T12*t12*T34*T34 + T12*T12*T34*t34 + T13*t13*T24*T24 + T13*T13*T24*t24 + T14*t14*T23*T23 + T14*T14*T23*t23);// 2-pair
			  weightTotal += pow(1-G,4)*(pow(T12,2)*pow(T34,2) + pow(T13,2)*pow(T24,2) + pow(T14,2)*pow(T23,2));// 2-pair fully chaotic
			  weightTotal += 2*G*pow(1-G,2)*(T12*T13*t23 + T12*T23*t13 + T13*T23*t12  + T12*T14*t24 + T12*T24*t14 + T14*T24*t12);// 3-pion
			  weightTotal += 2*G*pow(1-G,2)*(T13*T14*t34 + T13*T34*t14 + T14*T34*t13  + T23*T24*t34 + T23*T34*t24 + T24*T34*t23);// 3-pion
			  weightTotal += 2*pow(1-G,3)*(T12*T13*T23 + T12*T14*T24 + T13*T14*T34 + T23*T24*T34);// 3-pion fully chaotic
			  weightTotal += 2*G*pow(1-G,3)*(T12*t23*T34*T14 + T12*T23*t34*T14 + T12*T23*T34*t14 + t12*T23*T34*T14);// 4-pion
			  weightTotal += 2*G*pow(1-G,3)*(T12*t24*T34*T13 + T12*T24*T34*t13 + T12*T24*t34*T13 + t12*T24*T34*T13);// 4-pion
			  weightTotal += 2*G*pow(1-G,3)*(T13*T23*t24*T14 + T13*t23*T24*T14 + T13*T23*T24*t14 + t13*T23*T24*T14);// 4-pion
			  weightTotal += 2*pow(1-G,4)*(T12*T23*T34*T14 + T12*T24*T34*T13 + T13*T23*T24*T14);// 4-pion fully chaotic
			  //
			  weightPrime = weightTotal - 2*G*(1-G)*(T12*t12 + T13*t13 + T14*t14 + T23*t23 + T24*t24 + T34*t34);
			  weightPrime -= pow(1-G,2)*(T12*T12 + T13*T13 + T14*T14 + T23*T23 + T24*T24 + T34*T34);
			  weightPrimePrime = weightPrime - 2*G*pow(1-G,3)*(T12*t12*T34*T34 + T12*T12*T34*t34 + T13*t13*T24*T24 + T13*T13*T24*t24 + T14*t14*T23*T23 + T14*T14*T23*t23);
			  weightPrimePrime -= pow(1-G,4)*(pow(T12,2)*pow(T34,2) + pow(T13,2)*pow(T24,2) + pow(T14,2)*pow(T23,2));
			  weightCumulant = 2*G*pow(1-G,3)*(T12*t23*T34*T14 + T12*T23*t34*T14 + T12*T23*T34*t14 + t12*T23*T34*T14);
			  weightCumulant += 2*G*pow(1-G,3)*(T12*t24*T34*T13 + T12*T24*T34*t13 + T12*T24*t34*T13 + t12*T24*T34*T13);
			  weightCumulant += 2*G*pow(1-G,3)*(T13*T23*t24*T14 + T13*t23*T24*T14 + T13*T23*T24*t14 + t13*T23*T24*T14);
			  weightCumulant += 2*pow(1-G,4)*(T12*T23*T34*T14 + T12*T24*T34*T13 + T13*T23*T24*T14);
			  
			  
			  if(type==0){
			    if(Positive1stTripletWeights && Positive2ndTripletWeights){
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuild->Fill(FillBin, q4, weightTotal);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuild->Fill(FillBin, q4, weightPrime);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuild->Fill(FillBin, q4, weightPrimePrime);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuild->Fill(FillBin, q4, weightCumulant);
			    }else{
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildNeg->Fill(FillBin, q4, weightTotal);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuildNeg->Fill(FillBin, q4, weightPrime);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuildNeg->Fill(FillBin, q4, weightPrimePrime);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuildNeg->Fill(FillBin, q4, weightCumulant);
			    }
			  }else{
			    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildFromFits->Fill(type, FillBin, q4, weightTotal);
			    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimeBuildFromFits->Fill(type, FillBin, q4, weightPrime);
			    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fPrimePrimeBuildFromFits->Fill(type, FillBin, q4, weightPrimePrime);
			    Charge1[0].Charge2[0].Charge3[0].Charge4[0].MB[fMbin].EDB[EDindex4].FourPT[12].fCumulantBuildFromFits->Fill(type, FillBin, q4, weightCumulant);
			  }
			  
			}// GIndex 
		      }// RcohIndex
		    }// type
		    // stat errors
		    //weight14CC_e = weight14Err*MomResCorr14 / FSICorr14 / ffcSq * MuonCorr14;
		    //weight24CC_e = weight24Err*MomResCorr24 / FSICorr24 / ffcSq * MuonCorr24;
		    //weight34CC_e = weight34Err*MomResCorr34 / FSICorr34 / ffcSq * MuonCorr34;
		    //if(weight12CC[2]*weight13CC[2]*weight24CC[2]*weight34CC[2] > 0){
		    //weightTotalErr = pow( 6 * 2 * weight12CC_e*weight13CC[2]*weight24CC[2]*weight34CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight24CC[2]*weight34CC[2]),2);
		    //}
		    //if(weight12CC[2]*weight13CC[2]*weight23CC[2] > 0){
		    //weightTotalErr += pow( 8 * sqrt(3) * weight12CC_e*weight13CC[2]*weight23CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight23CC[2]),2);
		    //}
		    //weightTotalErr += 2*(pow(weight12CC_e*weight34CC[2],2) + pow(weight13CC_e*weight24CC[2],2) + pow(weight14CC_e*weight23CC[2],2));
		    //weightTotalErr += pow(weight12CC_e,2) + pow(weight13CC_e,2) + pow(weight14CC_e,2) + pow(weight23CC_e,2) + pow(weight24CC_e,2) + pow(weight34CC_e,2);
		    //Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[EDindex4].FourPT[12].fBuildErr->Fill(4, q4, weightTotalErr);
		  
		    // Radius estimations for c4
		    if(fMbin==0 && EDindex4==0 && fCollisionType==0){
		      for(Int_t Rindex=0; Rindex<7; Rindex++){
			Float_t R = (6. + Rindex)/FmToGeV;
			Float_t arg12=qinv12*R;
			Float_t arg13=qinv13*R;
			Float_t arg14=qinv14*R;
			Float_t arg23=qinv23*R;
			Float_t arg24=qinv24*R;
			Float_t arg34=qinv34*R;
			// Exchange Amplitudes
			Float_t EA12 = exp(-pow(arg12,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12));
			Float_t EA13 = exp(-pow(arg13,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12));
			Float_t EA14 = exp(-pow(arg14,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg14,3) - 12.*arg14) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg14,4) -48.*pow(arg14,2) + 12));
			Float_t EA23 = exp(-pow(arg23,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12));
			Float_t EA24 = exp(-pow(arg24,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg24,3) - 12.*arg24) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg24,4) -48.*pow(arg24,2) + 12));
			Float_t EA34 = exp(-pow(arg34,2)/2.)*(1 + kappa3Fit/(6.*pow(2.,1.5))*(8.*pow(arg34,3) - 12.*arg34) + kappa4Fit/(24.*pow(2.,2))*(16.*pow(arg34,4) -48.*pow(arg34,2) + 12));
			//
			Float_t TotalCorrelation = 1 + 2*(EA12*EA13*EA24*EA34 + EA12*EA14*EA23*EA34 + EA13*EA14*EA23*EA24);
			((TH2D*)fOutputList->FindObject("fc4QSFitNum"))->Fill(Rindex+1, q4, TotalCorrelation);
			((TH2D*)fOutputList->FindObject("fc4QSFitDen"))->Fill(Rindex+1, q4);
		      }
		    }
		  }// SC and ENsum=6
		  /////////////////////////////////////////////////////////////
		  
		  /*if(ch1==ch2 && ch1==ch3 && ch1==ch4){
		    Float_t pt1=sqrt(pow(pVect1[1],2)+pow(pVect1[2],2));
		    Float_t pt2=sqrt(pow(pVect2[1],2)+pow(pVect2[2],2));
		    Float_t pt3=sqrt(pow(pVect3[1],2)+pow(pVect3[2],2));
		    Float_t pt4=sqrt(pow(pVect4[1],2)+pow(pVect4[2],2));
		    if(ENsum==0){
		      ((TH3D*)fOutputList->FindObject("fKT4DistTerm1"))->Fill(fMbin+1, KT4, q4);
		      if(q4<0.105){
			((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, EDindex4, pt1);
			((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, EDindex4, pt2);
			((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, EDindex4, pt3);
			((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, EDindex4, pt4);
		      }
		    }
		    if(fMbin==0){
		      kT13 = sqrt(pow(pVect1[1]+pVect3[1],2) + pow(pVect1[2]+pVect3[2],2))/2.;
		      kT14 = sqrt(pow(pVect1[1]+pVect4[1],2) + pow(pVect1[2]+pVect4[2],2))/2.;
		      kT23 = sqrt(pow(pVect2[1]+pVect3[1],2) + pow(pVect2[2]+pVect3[2],2))/2.;
		      kT24 = sqrt(pow(pVect2[1]+pVect4[1],2) + pow(pVect2[2]+pVect4[2],2))/2.;
		      kT34 = sqrt(pow(pVect3[1]+pVect4[1],2) + pow(pVect3[2]+pVect4[2],2))/2.;
		      if(ENsum==0){
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum0"))->Fill(EDindex4, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum0"))->Fill(EDindex4, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum0"))->Fill(EDindex4, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum0"))->Fill(EDindex4, q4, pt4);
			//
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT12);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT13);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT14);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT23);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT24);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum0"))->Fill(EDindex4, q4, kT34);
		      }else if(ENsum==1){
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum1"))->Fill(EDindex4, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum1"))->Fill(EDindex4, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum1"))->Fill(EDindex4, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum1"))->Fill(EDindex4, q4, pt4);
			}else if(ENsum==2){
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum2"))->Fill(EDindex4, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum2"))->Fill(EDindex4, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum2"))->Fill(EDindex4, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum2"))->Fill(EDindex4, q4, pt4);
		      }else if(ENsum==3){
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum3"))->Fill(EDindex4, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum3"))->Fill(EDindex4, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum3"))->Fill(EDindex4, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum3"))->Fill(EDindex4, q4, pt4);
		      }else{// 6
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum6"))->Fill(EDindex4, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum6"))->Fill(EDindex4, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum6"))->Fill(EDindex4, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpTENsum6"))->Fill(EDindex4, q4, pt4);
			//
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT12);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT13);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT14);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT23);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT24);
			((TH3D*)fOutputList->FindObject("fQ4AvgkTENsum6"))->Fill(EDindex4, q4, kT34);
		      }
		      
		    }
		    }
		  
		    if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6) ((TH3D*)fOutputList->FindObject("fKT4DistTerm13"))->Fill(fMbin+1, KT4, q4);
		  */
		  
		  // momenumtum resolution and muon corrections
		  if(fMCcase && ENsum==6 && FilledMCtriplet123){// for momentum resolution and muon correction
		    if((fEvt+en4)->fTracks[l].fLabel < (fEvt+en4)->fMCarraySize){
		      
		      pVect4MC[0]=sqrt(pow((fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		      pVect4MC[1]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPx;
		      pVect4MC[2]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPy;
		      pVect4MC[3]=(fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPz;
		      qinv14MC = GetQinv(pVect1MC, pVect4MC);
		      qinv24MC = GetQinv(pVect2MC, pVect4MC);
		      qinv34MC = GetQinv(pVect3MC, pVect4MC);
		      
		      q4MC = sqrt(pow(q3MC,2) + pow(qinv14MC,2) +  pow(qinv24MC,2) +  pow(qinv34MC,2));
		      if(q4<0.1 && ch1==ch2 && ch1==ch3 && ch1==ch4) ((TH2D*)fOutputList->FindObject("fQ4Res"))->Fill(KT4, q4-q4MC);
		      if(ch1==ch2 && ch1==ch3 && ch1==ch4) {
			((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(1, qinv12MC); ((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(2, qinv13MC);
			((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(3, qinv14MC); ((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(4, qinv23MC);
			((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(5, qinv24MC); ((TH2D*)fOutputList->FindObject("DistQinvMC4pion"))->Fill(6, qinv34MC);
		      }

		      //Float_t QuadWeightTTC=1.0;// same-charge weights to mimic two-track depletion of same-charge pairs
		      //if(ch1==ch2 && qinv12>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv12);
		      //if(ch1==ch3 && qinv13>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv13);
		      //if(ch1==ch4 && qinv14>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv14);
		      //if(ch2==ch3 && qinv23>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv23);
		      //if(ch2==ch4 && qinv24>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv24);
		      //if(ch3==ch4 && qinv34>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv34);
		      

		      		      
		      Pparent4[0]=pVect4MC[0]; Pparent4[1]=pVect4MC[1]; Pparent4[2]=pVect4MC[2]; Pparent4[3]=pVect4MC[3];
		      pionParent4=kFALSE;
		      if(abs((fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fPdgCode)==13){// muon check
			Int_t MotherLabel4 = (fEvt+en4)->fMCtracks[abs((fEvt+en4)->fTracks[l].fLabel)].fMotherLabel;
			if(abs((fEvt+en4)->fMCtracks[MotherLabel4].fPdgCode)==211) {
			  pionParent4=kTRUE;
			  Pparent4[1] = (fEvt+en4)->fMCtracks[MotherLabel4].fPx; Pparent4[2] = (fEvt+en4)->fMCtracks[MotherLabel4].fPy; Pparent4[3] = (fEvt+en4)->fMCtracks[MotherLabel4].fPz;
			  Pparent4[0] = sqrt(pow(Pparent4[1],2)+pow(Pparent4[2],2)+pow(Pparent4[3],2)+pow(fTrueMassPi,2));
			}
		      }

		      parentQinv14 = GetQinv(Pparent1, Pparent4);
		      parentQinv24 = GetQinv(Pparent2, Pparent4);
		      parentQinv34 = GetQinv(Pparent3, Pparent4);
		      Float_t parentQ4 = sqrt(pow(parentQ3,2) + pow(parentQinv14,2) + pow(parentQinv24,2) + pow(parentQinv34,2));
		      
		      if(parentQinv14 > 0.001 && parentQinv24 > 0.001 && parentQinv34 > 0.001 && parentQ4 < 0.5){
			if(pionParent1 || pionParent2 || pionParent3 || pionParent4) {// want at least one pion-->muon
			 
			  if(pionParent1) ((TH1D*)fOutputList->FindObject("fDistPionParents4"))->Fill(1);
			  if(pionParent2) ((TH1D*)fOutputList->FindObject("fDistPionParents4"))->Fill(2);
			  if(pionParent3) ((TH1D*)fOutputList->FindObject("fDistPionParents4"))->Fill(3);
			  if(pionParent4) ((TH1D*)fOutputList->FindObject("fDistPionParents4"))->Fill(4);
			  Float_t parentQinvGroup4[6]={parentQinv12, parentQinv13, parentQinv14, parentQinv23, parentQinv24, parentQinv34};
			  Float_t parentkTGroup4[6]={0};
			  
			  for(Int_t term=1; term<=12; term++){
			    if(term==1) {}
			    else if(term==2) {if(!pionParent1 && !pionParent2 && !pionParent3) continue;}
			    else if(term==3) {if(!pionParent1 && !pionParent2 && !pionParent4) continue;}
			    else if(term==4) {if(!pionParent1 && !pionParent3 && !pionParent4) continue;}
			    else if(term==5) {if(!pionParent2 && !pionParent3 && !pionParent4) continue;}
			    else if(term==6) {if(!pionParent1 && !pionParent2) continue;}
			    else if(term==7) {if(!pionParent1 && !pionParent3) continue;}
			    else if(term==8) {if(!pionParent1 && !pionParent4) continue;}
			    else if(term==9) {if(!pionParent2 && !pionParent3) continue;}
			    else if(term==10) {if(!pionParent2 && !pionParent4) continue;}
			    else if(term==11) {if(!pionParent3 && !pionParent4) continue;}
			    else {} 
			    for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			      Float_t Rvalue = fRstartMC+Riter;
			      Float_t WInput = MCWeight4(term, Rvalue, 1.0, chGroup4, parentQinvGroup4, parentkTGroup4);
			      Float_t WInputParentFSI = MCWeightFSI4(term, Rvalue, 1.0, chGroup4, parentQinvGroup4);
			      Float_t WInputFSI = MCWeightFSI4(term, Rvalue, 1.0, chGroup4, QinvMCGroup4);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonSmeared->Fill(1, Rvalue, q4MC, WInput);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonIdeal->Fill(1, Rvalue, parentQ4, WInput);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonPionK4->Fill(1, Rvalue, q4MC, WInputFSI);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fPionPionK4->Fill(1, Rvalue, parentQ4, WInputParentFSI);
			      //
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonSmeared->Fill(2, Rvalue, q4MC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonIdeal->Fill(2, Rvalue, parentQ4);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonPionK4->Fill(2, Rvalue, q4MC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fPionPionK4->Fill(2, Rvalue, parentQ4);
			    }// Riter
			  }// term loop
			  
			}// pion parent check
		      }// parentQ check (muon correction)
		    
		      Int_t indexq4 = q4 / 0.005;
		      if(indexq4 >=50) indexq4=49; 
		      Float_t WSpectrum = 1.0;

		      // 4-pion momentum resolution
		      for(Int_t term=1; term<=13; term++){
			for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			  Float_t Rvalue = fRstartMC+Riter;
			  Float_t WInput = MCWeight4(term, Rvalue, ffcSqMRC, chGroup4, QinvMCGroup4, kTGroup4);
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[EDindex4].FourPT[term-1].fIdeal->Fill(Rvalue, q4MC, WInput*WSpectrum);
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[EDindex4].FourPT[term-1].fSmeared->Fill(Rvalue, q4, WInput*WSpectrum);
			}
		      }
		    
		    }// label check particle 4
		  }// MCcase
		  
		}// 4th particle
	      }// 3rd particle
	    }// 2nd particle
	  }// 1st particle
	  
	}// en4
      }// en3
    }// en2
    
    

  
  
  
  
  // Post output data.
  PostData(1, fOutputList);
  
}
//________________________________________________________________________
void AliFourPion::Terminate(Option_t *) 
{
  // Called once at the end of the query
 
  cout<<"Done"<<endl;

}
//________________________________________________________________________
Bool_t AliFourPion::AcceptPair(AliFourPionTrackStruct first, AliFourPionTrackStruct second)
{
  
  if(fabs(first.fEta-second.fEta) > fMinSepPairEta) return kTRUE;
  
  // propagate through B field to r=1m
  Float_t phi1 = first.fPhi - asin(first.fCharge*(0.1*fBfield)*0.15/first.fPt);// 0.15 for D=1m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  Float_t phi2 = second.fPhi - asin(second.fCharge*(0.1*fBfield)*0.15/second.fPt);// 0.15 for D=1m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  Float_t deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
    
  
  // propagate through B field to r=1.6m
  phi1 = first.fPhi - asin(first.fCharge*(0.1*fBfield)*0.24/first.fPt);// mine. 0.24 for D=1.6m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  phi2 = second.fPhi - asin(second.fCharge*(0.1*fBfield)*0.24/second.fPt);// mine. 0.24 for D=1.6m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
  
  
   
  //
  
  //Int_t ncl1 = first.fClusterMap.GetNbits();
  //Int_t ncl2 = second.fClusterMap.GetNbits();
  //Int_t sumCls = 0; Int_t sumSha = 0; Int_t sumQ = 0;
  //Double_t shfrac = 0; Double_t qfactor = 0;
  //for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
  //if (first.fClusterMap.TestBitNumber(imap) && second.fClusterMap.TestBitNumber(imap)) {// Both clusters
  //  if (first.fSharedMap.TestBitNumber(imap) && second.fSharedMap.TestBitNumber(imap)) { // Shared
  //	sumQ++;
  //	sumCls+=2;
  //	sumSha+=2;}
  //  else {sumQ--; sumCls+=2;}
  //}
  //else if (first.fClusterMap.TestBitNumber(imap) || second.fClusterMap.TestBitNumber(imap)) {// Non shared
  //  sumQ++;
  //  sumCls++;}
  //}
  //if (sumCls>0) {
  //qfactor = sumQ*1.0/sumCls;
  //shfrac = sumSha*1.0/sumCls;
  //}
  
  //if(qfactor > fShareQuality || shfrac > fShareFraction) return kFALSE;
  
  
  return kTRUE;
  

}
//________________________________________________________________________
Bool_t AliFourPion::AcceptPairPM(AliFourPionTrackStruct first, AliFourPionTrackStruct second)
{// optional pair cuts for +- pairs
  
  if(fabs(first.fEta-second.fEta) > fMinSepPairEta) return kTRUE;
  
  // propagate through B field to r=1m
  Float_t phi1 = first.fPhi - asin(1.*(0.1*fBfield)*0.15/first.fPt);// 0.15 for D=1m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  Float_t phi2 = second.fPhi - asin(1.*(0.1*fBfield)*0.15/second.fPt);// 0.15 for D=1m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  Float_t deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
    
  
  // propagate through B field to r=1.6m
  phi1 = first.fPhi - asin(1.*(0.1*fBfield)*0.24/first.fPt);// mine. 0.24 for D=1.6m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  phi2 = second.fPhi - asin(1.*(0.1*fBfield)*0.24/second.fPt);// mine. 0.24 for D=1.6m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
  
  return kTRUE;
  
}
//________________________________________________________________________
Float_t AliFourPion::Gamov(Int_t chargeBin1, Int_t chargeBin2, Float_t qinv)
{
  Float_t arg = G_Coeff/qinv;
  
  if(chargeBin1==chargeBin2) return (exp(arg)-1)/(arg);
  else {return (exp(-arg)-1)/(-arg);}
  
}
//________________________________________________________________________
void AliFourPion::Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
{
  Int_t j, k;
  Int_t a = i2 - i1;
  for (Int_t i = i1; i < i2+1; i++) {
    j = (Int_t) (gRandom->Rndm() * a);
    k = iarr[j];
    iarr[j] = iarr[i];
    iarr[i] = k;
  }
}


//________________________________________________________________________
Float_t AliFourPion::GetQinv(Float_t track1[], Float_t track2[]){
  
  Float_t qinv = sqrt( fabs(pow(track1[1]-track2[1],2) + pow(track1[2]-track2[2],2) + pow(track1[3]-track2[3],2) - pow(track1[0]-track2[0],2)) );
  return qinv;
  
}
//________________________________________________________________________
void AliFourPion::GetQosl(Float_t track1[], Float_t track2[], Float_t& qout, Float_t& qside, Float_t& qlong){
 
  Float_t p0 = track1[0] + track2[0];
  Float_t px = track1[1] + track2[1];
  Float_t py = track1[2] + track2[2];
  Float_t pz = track1[3] + track2[3];
  
  Float_t mt = sqrt(p0*p0 - pz*pz);
  Float_t pt = sqrt(px*px + py*py);
  
  Float_t v0 = track1[0] - track2[0];
  Float_t vx = track1[1] - track2[1];
  Float_t vy = track1[2] - track2[2];
  Float_t vz = track1[3] - track2[3];
  
  qout = (px*vx + py*vy)/pt;
  qside = (px*vy - py*vx)/pt;
  qlong = (p0*vz - pz*v0)/mt;
}
//________________________________________________________________________
void AliFourPion::SetWeightArrays(Bool_t legoCase, TH3F *histos[AliFourPion::fKbinsT][AliFourPion::fCentBins], TH3F *histos2[AliFourPion::fKbinsT][AliFourPion::fCentBins], TH2F *histos1D[AliFourPion::fCentBins]){
  
  if(legoCase){
    cout<<"LEGO call to SetWeightArrays"<<endl;
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t mb=0; mb<fMbins; mb++){
	fNormWeight[tKbin][mb] = (TH3F*)histos[tKbin][mb]->Clone();
	fNormWeight[tKbin][mb]->SetDirectory(0);
	fNormWeight2[tKbin][mb] = (TH3F*)histos2[tKbin][mb]->Clone();
	fNormWeight2[tKbin][mb]->SetDirectory(0);
	if(tKbin==0){	
	  fNormWeightOneD[mb] = (TH2F*)histos1D[mb]->Clone();
	  fNormWeightOneD[mb]->SetDirectory(0);
	}
      }
    }
    
  }else{
    TString *wFileName=new TString("WeightFile_");
    if(fCollisionType==0) wFileName->Append("PbPb.root");
    else if(fCollisionType==1) wFileName->Append("pPb.root");
    else wFileName->Append("pp.root");
    TFile *wFile = new TFile(wFileName->Data(),"READ");
    if(!wFile->IsOpen()) {cout<<"No Weight File!!!!!!!!!!"<<endl; return;}
    else cout<<"Good Weight File Found!"<<endl;
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t mb=0; mb<fMbins; mb++){
	for(Int_t q2bin=0; q2bin<2; q2bin++){
	  TString *name = new TString("Weight_Kt_");
	  *name += tKbin;
	  name->Append("_Ky_0");
	  name->Append("_M_");
	  *name += mb;
	  name->Append("_ED_");
	  *name += q2bin;
	  
	  if(q2bin==0) {
	    fNormWeight[tKbin][mb] = (TH3F*)wFile->Get(name->Data());
	    fNormWeight[tKbin][mb]->SetDirectory(0);
	  }else{
	    fNormWeight2[tKbin][mb] = (TH3F*)wFile->Get(name->Data());
	    fNormWeight2[tKbin][mb]->SetDirectory(0);
	  }
	  if(tKbin==0 && q2bin==0){
	    TString *name1D = new TString("Weight_M_");
	    *name1D += mb;
	    name1D->Append("_1D");
	    fNormWeightOneD[mb] = (TH2F*)wFile->Get(name1D->Data());
	    fNormWeightOneD[mb]->SetDirectory(0);
	  }
	}//q2bin
      }//mb
    }//kt
    
    wFile->Close();
  }
  
  cout<<"Done reading weight file"<<endl;
  
}
//________________________________________________________________________
void AliFourPion::GetWeight(Float_t track1[], Float_t track2[], Float_t& wgt, Float_t& wgtErr){
  
  Float_t kt=sqrt( pow(track1[1]+track2[1],2) + pow(track1[2]+track2[2],2))/2.;
  //
  Float_t qOut=0,qSide=0,qLong=0;
  GetQosl(track1, track2, qOut, qSide, qLong);
  qOut = fabs(qOut);
  qSide = fabs(qSide);
  qLong = fabs(qLong);
  Float_t wd=0, xd=0, yd=0, zd=0;
  //Float_t qInvtemp=GetQinv(track1, track2);
  
  //
  Int_t q2bin=0;
  if(fq2Binning || fLowMultBinning) q2bin = fEDbin;
  //
  
  if(kt < fKmeanT[0]) {fKtIndexL=0; fKtIndexH=1;}
  else if(kt >= fKmeanT[fKbinsT-1]) {fKtIndexL=fKbinsT-2; fKtIndexH=fKbinsT-1;}
  else {
    for(Int_t i=0; i<fKbinsT-1; i++){
      if((kt >= fKmeanT[i]) && (kt < fKmeanT[i+1])) {fKtIndexL=i; fKtIndexH=i+1; break;}
    }
  }
  wd = (kt-fKmeanT[fKtIndexL])/(fKmeanT[fKtIndexH]-fKmeanT[fKtIndexL]);
  if(fMaxPt<=0.251) {fKtIndexL=0; fKtIndexH=0; wd=0;}
  if(fMinPt>0.249 && fKtIndexL==0) {fKtIndexL=1; wd=0;}
  
  
  //
  if(fOneDInterpolation){
    // different kT binning for Qinv weights
    if(kt < fKmeanTOneD[0]) {fKtIndexL=0; fKtIndexH=1;}
    else if(kt >= fKmeanTOneD[fKbinsTOneD-1]) {fKtIndexL=fKbinsTOneD-2; fKtIndexH=fKbinsTOneD-1;}
    else {
      for(Int_t i=0; i<fKbinsTOneD-1; i++){
	if((kt >= fKmeanTOneD[i]) && (kt < fKmeanTOneD[i+1])) {fKtIndexL=i; fKtIndexH=i+1; break;}
      }
    }
    wd = (kt-fKmeanTOneD[fKtIndexL])/(fKmeanTOneD[fKtIndexH]-fKmeanTOneD[fKtIndexL]);
    //
    Float_t qInv=GetQinv(track1, track2);
    Int_t MaxBinsQinv=20;
    if(fCollisionType!=0) MaxBinsQinv=100;
    if(qInv < 0.0075) {fQinvIndexL=1; fQinvIndexH=1; xd=0;}
    else if(qInv >= 0.0975 && fCollisionType==0) {fQinvIndexL=19; fQinvIndexH=19; xd=1;}
    else if(qInv >= 0.4975 && fCollisionType!=0) {fQinvIndexL=99; fQinvIndexH=99; xd=1;}
    else {
      for(Int_t i=0; i<MaxBinsQinv-1; i++){
	if((qInv >= (i+0.5)*0.005) && (qInv < (i+1.5)*0.005)) {fQinvIndexL=i; fQinvIndexH=i+1; break;}
      }
      xd = (qInv-(fQinvIndexL+0.5)*0.005)/(0.005);
    }
    if(fInterpolationType==0){// Linear Interpolation of qinv,kt
      Float_t c0 = fNormWeightOneD[fMbin]->GetBinContent(fKtIndexL+1, fQinvIndexL+1)*(1-xd) + fNormWeightOneD[fMbin]->GetBinContent(fKtIndexL+1, fQinvIndexH+1)*xd;
      Float_t c1 = fNormWeightOneD[fMbin]->GetBinContent(fKtIndexH+1, fQinvIndexL+1)*(1-xd) + fNormWeightOneD[fMbin]->GetBinContent(fKtIndexH+1, fQinvIndexH+1)*xd;
      // kT interpolation
      wgt = (c0*(1-wd) + c1*wd);
    }else{// cubic Interpolation
      Float_t fOneDarrP1[4]={0};
      Float_t fOneDarrP2[4]={0};
      for(Int_t x=0; x<4; x++){
	Int_t binInv = fQinvIndexL + x;
	if(binInv<=0) binInv = 1;
	if(binInv>MaxBinsQinv) binInv = MaxBinsQinv;
	fOneDarrP1[x] = fNormWeightOneD[fMbin]->GetBinContent(fKtIndexL+1, binInv);
	fOneDarrP2[x] = fNormWeightOneD[fMbin]->GetBinContent(fKtIndexH+1, binInv);
      }
      Float_t coord[3]={xd, yd, zd}; 
      Float_t c0 = cubicInterpolate (fOneDarrP1, xd);
      Float_t c1 = cubicInterpolate (fOneDarrP2, xd);
      // kT interpolation
      wgt = c0*(1-wd) + c1*wd;
    }
    
    // simplified stat error 
    Float_t avgErr = fNormWeightOneD[fMbin]->GetBinError(fKtIndexL+1, fQinvIndexL+1);
    avgErr += fNormWeightOneD[fMbin]->GetBinError(fKtIndexH+1, fQinvIndexH+1);
    avgErr /= 2.;
    //
    wgtErr = avgErr;
  }else{
    //
    if(qOut < fQmean[0]) {fQoIndexL=0; fQoIndexH=0; xd=0;}
    else if(qOut >= fQmean[kQbinsWeights-1]) {fQoIndexL=kQbinsWeights-1; fQoIndexH=kQbinsWeights-1; xd=1;}
    else {
      for(Int_t i=0; i<kQbinsWeights-1; i++){
	if((qOut >= fQmean[i]) && (qOut < fQmean[i+1])) {fQoIndexL=i; fQoIndexH=i+1; break;}
      }
      xd = (qOut-fQmean[fQoIndexL])/(fQmean[fQoIndexH]-fQmean[fQoIndexL]);
    }
    //
    if(qSide < fQmean[0]) {fQsIndexL=0; fQsIndexH=0; yd=0;}
    else if(qSide >= fQmean[kQbinsWeights-1]) {fQsIndexL=kQbinsWeights-1; fQsIndexH=kQbinsWeights-1; yd=1;}
    else {
      for(Int_t i=0; i<kQbinsWeights-1; i++){
	if((qSide >= fQmean[i]) && (qSide < fQmean[i+1])) {fQsIndexL=i; fQsIndexH=i+1; break;}
      }
      yd = (qSide-fQmean[fQsIndexL])/(fQmean[fQsIndexH]-fQmean[fQsIndexL]);
    }
    //
    if(qLong < fQmean[0]) {fQlIndexL=0; fQlIndexH=0; zd=0;}
    else if(qLong >= fQmean[kQbinsWeights-1]) {fQlIndexL=kQbinsWeights-1; fQlIndexH=kQbinsWeights-1; zd=1;}
    else {
      for(Int_t i=0; i<kQbinsWeights-1; i++){
	if((qLong >= fQmean[i]) && (qLong < fQmean[i+1])) {fQlIndexL=i; fQlIndexH=i+1; break;}
      }
      zd = (qLong-fQmean[fQlIndexL])/(fQmean[fQlIndexH]-fQmean[fQlIndexL]);
    }
    
    //
    if(fInterpolationType==0){// Linear Interpolation of osl
      // w interpolation (kt)
      Float_t c000, c100, c010, c001, c110, c101, c011, c111;
      if(q2bin==0){
	c000 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*wd;
	c100 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*wd;
	c010 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*wd;
	c001 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*wd;
	c110 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*wd;
	c101 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*wd;
	c011 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*wd;
	c111 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*wd;
      }else{
	c000 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*wd;
	c100 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*wd;
	c010 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*wd;
	c001 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*wd;
	c110 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*wd;
	c101 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*wd;
	c011 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*wd;
	c111 = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight2[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*wd;
      }
      //
      
      // x interpolation (qOut)
      Float_t c00 = c000*(1-xd) + c100*xd;
      Float_t c10 = c010*(1-xd) + c110*xd;
      Float_t c01 = c001*(1-xd) + c101*xd;
      Float_t c11 = c011*(1-xd) + c111*xd;
      // y interpolation (qSide)
      Float_t c0 = c00*(1-yd) + c10*yd;
      Float_t c1 = c01*(1-yd) + c11*yd;
      // z interpolation (qLong)
      wgt = (c0*(1-zd) + c1*zd);
      // testing
      //Float_t c0_2 = fNormWeight[fKtIndexL][fMbin]->Interpolate(qOut,qSide,qLong);
      //Float_t c1_2 = fNormWeight[fKtIndexH][fMbin]->Interpolate(qOut,qSide,qLong);
      //Float_t wgt_2 = c0_2*(1-wd) + c1_2*wd;
      //if(qinvtemp<0.02) cout<<qinvtemp<<"  "<<wgt<<"  "<<wgt_2<<"  "<<ExceedsMeans<<"  "<<qOut<<" "<<qSide<<" "<<qLong<<endl;
    }else{// cubic interpolation of osl
      
      for(Int_t x=0; x<4; x++){
	for(Int_t y=0; y<4; y++){
	  for(Int_t z=0; z<4; z++){
	    Int_t binO = fQoIndexL + x;
	    Int_t binS = fQsIndexL + y;
	    Int_t binL = fQlIndexL + z;
	    if(binO<=0) binO = 1;
	    if(binS<=0) binS = 1;
	    if(binL<=0) binL = 1;
	    if(binO>kQbinsWeights) binO = kQbinsWeights;
	    if(binS>kQbinsWeights) binS = kQbinsWeights;
	    if(binL>kQbinsWeights) binL = kQbinsWeights;
	    if(q2bin==0) {
	      farrP1[x][y][z] = fNormWeight[fKtIndexL][fMbin]->GetBinContent(binO,binS,binL);
	      farrP2[x][y][z] = fNormWeight[fKtIndexH][fMbin]->GetBinContent(binO,binS,binL);
	    }else{
	      farrP1[x][y][z] = fNormWeight2[fKtIndexL][fMbin]->GetBinContent(binO,binS,binL);
	      farrP2[x][y][z] = fNormWeight2[fKtIndexH][fMbin]->GetBinContent(binO,binS,binL);
	    }
	    
	  }
	}
      }
      Float_t coord[3]={xd, yd, zd}; 
      Float_t c0 = nCubicInterpolate(3, (Float_t*) farrP1, coord);
      Float_t c1 = nCubicInterpolate(3, (Float_t*) farrP2, coord);
      // kT interpolation
      wgt = c0*(1-wd) + c1*wd;
      //
      if(fCollisionType==0 && fInterpCorrection){
	Float_t qInvtemp=GetQinv(track1, track2);
	Int_t qInvIndex=int(qInvtemp/0.005);
	if(qInvIndex>19) qInvIndex=19;
	Int_t KtIndex = fKtIndexL;
	if(kt > (fKmiddleT[fKtIndexL]+fKstepT[fKtIndexL]/2.)) KtIndex++;
	Int_t MIndex = int(fMbin/2.);
	if(MIndex>4) MIndex=4;
	//
	wgt *= fIC[MIndex][KtIndex][qInvIndex];
      }
    }
    // simplified stat error 
    Float_t avgErr = fNormWeight[fKtIndexL][fMbin]->GetBinError(fQoIndexH+1,fQsIndexH+1,fQlIndexH+1);
    avgErr += fNormWeight[fKtIndexH][fMbin]->GetBinError(fQoIndexL+1,fQsIndexL+1,fQlIndexL+1);
    avgErr /= 2.;
    //
    wgtErr = avgErr;
  }
  
}
//________________________________________________________________________
Float_t AliFourPion::MCWeight(Int_t c[2], Float_t R, Float_t fcSq, Float_t qinv, Float_t k12){
  
  Float_t radius = R/0.19733;// convert to GeV (starts at 5 fm, was 3 fm)
  Float_t r12=radius*(1-k12/2.0);
  SetFSIindex(R);
  Float_t coulCorr12 = FSICorrelation(c[0], c[1], qinv);
  if(c[0]==c[1]){
    Float_t arg=qinv*r12;
    Float_t EW = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg,3) - 12.*arg);
    EW += kappa4/(24.*pow(2.,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
    return ((1-fcSq) + fcSq*(1 + exp(-pow(qinv*r12,2))*pow(EW,2))*coulCorr12);
  }else {
    return ((1-fcSq) + fcSq*coulCorr12);
  }
    
}
//________________________________________________________________________
Float_t AliFourPion::MCWeightOSL(Int_t charge1, Int_t charge2, Int_t r, Int_t dampIndex, Float_t qinv, Float_t qo, Float_t qs, Float_t ql){
  
  Float_t radiusOut = Float_t(r)/0.19733;// convert to GeV (starts at 5 fm, was 3 fm)
  Float_t radiusSide = radiusOut;
  Float_t radiusLong = radiusOut;
  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t coulCorr12 = FSICorrelation(charge1, charge2, qinv);
  if(charge1==charge2){
    return ((1-myDamp) + myDamp*(1 + exp(-pow(qo*radiusOut,2)) * exp(-pow(qs*radiusSide,2)) * exp(-pow(ql*radiusLong,2)))*coulCorr12);
  }else {
    return ((1-myDamp) + myDamp*coulCorr12);
  }
    
}

//________________________________________________________________________
Float_t AliFourPion::MCWeight3(Int_t term, Float_t R, Float_t fcSq, Int_t c[3], Float_t qinv[3], Float_t kT[3]){
  // FSI + QS correlations
  if(term==5) return 1.0;
  
  Float_t radius=R/0.19733;
  Float_t r12=radius*(1-kT[0]/2.0);
  Float_t r13=radius*(1-kT[1]/2.0);
  Float_t r23=radius*(1-kT[2]/2.0);
 
  Float_t fc = sqrt(fcSq);
  
  SetFSIindex(R);
  Float_t Kfactor12 = FSICorrelation(c[0],c[1], qinv[0]);// K2
  Float_t Kfactor13 = FSICorrelation(c[0],c[2], qinv[1]);// K2
  Float_t Kfactor23 = FSICorrelation(c[1],c[2], qinv[2]);// K2
  
  if(c[0]==c[1] && c[0]==c[2]){// all three of the same charge
    Float_t arg12=qinv[0]*r12;
    Float_t arg13=qinv[1]*r13;
    Float_t arg23=qinv[2]*r23;
    Float_t EW12 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
    EW12 += kappa4/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
    Float_t EW13 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
    EW13 += kappa4/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
    Float_t EW23 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
    EW23 += kappa4/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
    if(term==1){
      Float_t C3QS = 1 + exp(-pow(qinv[0]*r12,2))*pow(EW12,2) + exp(-pow(qinv[1]*r13,2))*pow(EW13,2) + exp(-pow(qinv[2]*r23,2))*pow(EW23,2);
      C3QS += 2*exp(-(pow(r12,2)*pow(qinv[0],2) + pow(r13,2)*pow(qinv[1],2) + pow(r23,2)*pow(qinv[2],2))/2.)*EW12*EW13*EW23;
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*(1+exp(-pow(qinv[0]*r12,2))*pow(EW12,2))*Kfactor12;
      C3 += pow(fc,2)*(1-fc)*(1+exp(-pow(qinv[1]*r13,2))*pow(EW13,2))*Kfactor13;
      C3 += pow(fc,2)*(1-fc)*(1+exp(-pow(qinv[2]*r23,2))*pow(EW23,2))*Kfactor23;
      C3 += pow(fc,3)*C3QS*Kfactor12*Kfactor13*Kfactor23;
      return C3;
    }else if(term==2){
      return ((1-fcSq) + fcSq*(1 + exp(-pow(qinv[0]*r12,2))*pow(EW12,2))*Kfactor12);
    }else if(term==3){
      return ((1-fcSq) + fcSq*(1 + exp(-pow(qinv[1]*r13,2))*pow(EW13,2))*Kfactor13);
    }else if(term==4){
      return ((1-fcSq) + fcSq*(1 + exp(-pow(qinv[2]*r23,2))*pow(EW23,2))*Kfactor23);
    }else return 1.0;
    
  }else{// mixed charge case
    Float_t arg=qinv[0]*r12;
    Float_t KfactorSC = Kfactor12;
    Float_t KfactorMC1 = Kfactor13;
    Float_t KfactorMC2 = Kfactor23;
    if(c[0]==c[2]) {arg=qinv[1]*r13; KfactorSC = Kfactor13; KfactorMC1 = Kfactor12; KfactorMC2 = Kfactor23;} 
    if(c[1]==c[2]) {arg=qinv[2]*r23; KfactorSC = Kfactor23; KfactorMC1 = Kfactor12; KfactorMC2 = Kfactor13;}
    Float_t EW = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg,3) - 12.*arg);
    EW += kappa4/(24.*pow(2.,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
    if(term==1){
      Float_t C3QS = 1 + exp(-pow(arg,2))*pow(EW,2);
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*(1+exp(-pow(arg,2))*pow(EW,2))*KfactorSC;
      C3 += pow(fc,2)*(1-fc)*KfactorMC1;
      C3 += pow(fc,2)*(1-fc)*KfactorMC2;
      C3 += pow(fc,3)*C3QS*KfactorSC*KfactorMC1*KfactorMC2;
      return C3;
    }else if(term==2){
      if( (c[0]+c[1]+c[2])==1) return ((1-fcSq) + fcSq*(1 + exp(-pow(arg,2))*pow(EW,2))*KfactorSC);
      else return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else if(term==3){
      return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else if(term==4){
      if( (c[0]+c[1]+c[2])==2) return ((1-fcSq) + fcSq*(1 + exp(-pow(arg,2))*pow(EW,2))*KfactorSC);
      else return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else return 1.0;
  }
  
}
//________________________________________________________________________
Float_t AliFourPion::MCWeightFSI3(Int_t term, Float_t R, Float_t fcSq, Int_t c[3], Float_t qinv[3]){
  // FSI only (no QS correlations)
  if(term==5) return 1.0;
  
  Float_t fc = sqrt(fcSq);
  SetFSIindex(R);
  Float_t Kfactor12 = FSICorrelation(c[0],c[1], qinv[0]);// K2
  Float_t Kfactor13 = FSICorrelation(c[0],c[2], qinv[1]);// K2
  Float_t Kfactor23 = FSICorrelation(c[1],c[2], qinv[2]);// K2
  
  if(c[0]==c[1] && c[0]==c[2]){// all three of the same charge
    if(term==1){
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*Kfactor12;
      C3 += pow(fc,2)*(1-fc)*Kfactor13;
      C3 += pow(fc,2)*(1-fc)*Kfactor23;
      C3 += pow(fc,3)*Kfactor12*Kfactor13*Kfactor23;
      return C3;
    }else if(term==2){
      return ((1-fcSq) + fcSq*Kfactor12);
    }else if(term==3){
      return ((1-fcSq) + fcSq*Kfactor13);
    }else if(term==4){
      return ((1-fcSq) + fcSq*Kfactor23);
    }else return 1.0;
    
  }else{// mixed charge case
    Float_t KfactorSC = Kfactor12;
    Float_t KfactorMC1 = Kfactor13;
    Float_t KfactorMC2 = Kfactor23;
    if(c[0]==c[2]) {KfactorSC = Kfactor13; KfactorMC1 = Kfactor12; KfactorMC2 = Kfactor23;} 
    if(c[1]==c[2]) {KfactorSC = Kfactor23; KfactorMC1 = Kfactor12; KfactorMC2 = Kfactor13;}
    if(term==1){
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*KfactorSC;
      C3 += pow(fc,2)*(1-fc)*KfactorMC1;
      C3 += pow(fc,2)*(1-fc)*KfactorMC2;
      C3 += pow(fc,3)*KfactorSC*KfactorMC1*KfactorMC2;
      return C3;
    }else if(term==2){
      if( (c[0]+c[1]+c[2])==1) return ((1-fcSq) + fcSq*KfactorSC);
      else return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else if(term==3){
      return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else if(term==4){
      if( (c[0]+c[1]+c[2])==2) return ((1-fcSq) + fcSq*KfactorSC);
      else return ((1-fcSq) + fcSq*KfactorMC1);// doesn't matter MC1 or MC2
    }else return 1.0;
  }
  
}
//________________________________________________________________________
Float_t AliFourPion::MCWeight4(Int_t term, Float_t R, Float_t fcSq, Int_t c[4], Float_t qinv[6], Float_t kT[6]){
  if(term==13) return 1.0;

  // Charge ordering:
  // ----, ---+, --++, -+++, ++++
  //
  // term ordering:
  // Term 1: 1-2-3-4  (all same-event)
  // Term 2: 1-2-3 4 (particle 4 from different event)
  // Term 3: 1-2-4 3 (particle 3 from different event)
  // Term 4: 1-3-4 2 (particle 2 from different event)
  // Term 5: 2-3-4 1 (particle 1 from different event)
  // Term 6: 1-2 3 4 (particle 1 and 2 from same event)
  // Term 7: 1-3 2 4
  // Term 8: 1-4 2 3
  // Term 9: 2-3 1 4
  // Term 10: 2-4 1 3
  // Term 11: 3-4 1 2
  // Term 12: 1 2 3 4 (all from different events)

  Float_t radius = R/0.19733;
  Float_t r[6]={0};
  r[0]=radius*(1-kT[0]/2.0);
  r[1]=radius*(1-kT[1]/2.0);
  r[2]=radius*(1-kT[2]/2.0);
  r[3]=radius*(1-kT[3]/2.0);
  r[4]=radius*(1-kT[4]/2.0);
  r[5]=radius*(1-kT[5]/2.0);
    
  Int_t ChargeSum=c[0]+c[1]+c[2]+c[3];
 
  Float_t fc = sqrt(fcSq);
  SetFSIindex(R);
  Float_t Kfactor12 = FSICorrelation(c[0],c[1], qinv[0]);// K2
  Float_t Kfactor13 = FSICorrelation(c[0],c[2], qinv[1]);// K2
  Float_t Kfactor14 = FSICorrelation(c[0],c[3], qinv[2]);// K2
  Float_t Kfactor23 = FSICorrelation(c[1],c[2], qinv[3]);// K2
  Float_t Kfactor24 = FSICorrelation(c[1],c[3], qinv[4]);// K2
  Float_t Kfactor34 = FSICorrelation(c[2],c[3], qinv[5]);// K2
  Float_t arg12=qinv[0]*r[0];
  Float_t arg13=qinv[1]*r[1];
  Float_t arg14=qinv[2]*r[2];
  Float_t arg23=qinv[3]*r[3];
  Float_t arg24=qinv[4]*r[4];
  Float_t arg34=qinv[5]*r[5];
  // Exchange Amplitudes
  Float_t EA12 = exp(-pow(arg12,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12) + kappa4/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12));
  Float_t EA13 = exp(-pow(arg13,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13) + kappa4/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12));
  Float_t EA14 = exp(-pow(arg14,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg14,3) - 12.*arg14) + kappa4/(24.*pow(2.,2))*(16.*pow(arg14,4) -48.*pow(arg14,2) + 12));
  Float_t EA23 = exp(-pow(arg23,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23) + kappa4/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12));
  Float_t EA24 = exp(-pow(arg24,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg24,3) - 12.*arg24) + kappa4/(24.*pow(2.,2))*(16.*pow(arg24,4) -48.*pow(arg24,2) + 12));
  Float_t EA34 = exp(-pow(arg34,2)/2.)*(1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg34,3) - 12.*arg34) + kappa4/(24.*pow(2.,2))*(16.*pow(arg34,4) -48.*pow(arg34,2) + 12));
  
  if(c[0]==c[1] && c[0]==c[2] && c[0]==c[3]){// ---- and ++++ configuration
    
    if(term==1){
      Float_t C4QS = 1 + pow(EA12,2) + pow(EA13,2) + pow(EA14,2) + pow(EA23,2) + pow(EA24,2) + pow(EA34,2);// baseline + single pairs
      C4QS += pow(EA12,2) * pow(EA34,2);// 2-pairs
      C4QS += pow(EA13,2) * pow(EA24,2);// 2-pairs
      C4QS += pow(EA14,2) * pow(EA23,2);// 2-pairs
      C4QS += 2*EA12*EA13*EA23 + 2*EA12*EA14*EA24 + 2*EA13*EA14*EA34 + 2*EA23*EA24*EA34;// 3-particle exhange
      C4QS += 3*EA12*EA23*EA34*EA14 + 3*EA12*EA13*EA34*EA24;// 4-particle exchange
      Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
      C4 += pow(fc,2)*pow(1-fc,2)*( (1 + pow(EA12,2))*Kfactor12 + (1 + pow(EA13,2))*Kfactor13 + (1 + pow(EA14,2))*Kfactor14 );
      C4 += pow(fc,2)*pow(1-fc,2)*( (1 + pow(EA23,2))*Kfactor23 + (1 + pow(EA24,2))*Kfactor24 + (1 + pow(EA34,2))*Kfactor34);
      C4 += pow(fc,3)*(1-fc)*(1 + pow(EA12,2) + pow(EA13,2) + pow(EA23,2) + 2*EA12*EA13*EA23) * Kfactor12*Kfactor13*Kfactor23;
      C4 += pow(fc,3)*(1-fc)*(1 + pow(EA12,2) + pow(EA14,2) + pow(EA24,2) + 2*EA12*EA14*EA24) * Kfactor12*Kfactor14*Kfactor24;
      C4 += pow(fc,3)*(1-fc)*(1 + pow(EA13,2) + pow(EA14,2) + pow(EA34,2) + 2*EA13*EA14*EA34) * Kfactor13*Kfactor14*Kfactor34;
      C4 += pow(fc,3)*(1-fc)*(1 + pow(EA23,2) + pow(EA24,2) + pow(EA34,2) + 2*EA23*EA24*EA34) * Kfactor23*Kfactor24*Kfactor34;
      C4 += pow(fc,4)*C4QS*Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
      return C4;
    }else if(term<=5){
      Float_t EA1=0, EA2=0, EA3=0, Kpair1=0, Kpair2=0, Kpair3=0;
      if(term==2)      {EA1=EA12; EA2=EA13; EA3=EA23; Kpair1=Kfactor12; Kpair2=Kfactor13; Kpair3=Kfactor23;}
      else if(term==3) {EA1=EA12; EA2=EA14; EA3=EA24; Kpair1=Kfactor12; Kpair2=Kfactor14; Kpair3=Kfactor24;}
      else if(term==4) {EA1=EA13; EA2=EA14; EA3=EA34; Kpair1=Kfactor13; Kpair2=Kfactor14; Kpair3=Kfactor34;}
      else             {EA1=EA23; EA2=EA24; EA3=EA34; Kpair1=Kfactor23; Kpair2=Kfactor24; Kpair3=Kfactor34;}
      Float_t C3QS = 1 + pow(EA1,2) + pow(EA2,2) + pow(EA3,2);
      C3QS += 2*EA1*EA2*EA3;
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*( (1+pow(EA1,2))*Kpair1 + (1+pow(EA2,2))*Kpair2 + (1+pow(EA3,2))*Kpair3 );
      C3 += pow(fc,3)*C3QS*Kpair1*Kpair2*Kpair3;
      return C3;
    }else if(term<=11){
      if(term==6)       return ((1-fcSq) + fcSq*(1 + pow(EA12,2))*Kfactor12);
      else if(term==7)  return ((1-fcSq) + fcSq*(1 + pow(EA13,2))*Kfactor13);
      else if(term==8)  return ((1-fcSq) + fcSq*(1 + pow(EA14,2))*Kfactor14);
      else if(term==9)  return ((1-fcSq) + fcSq*(1 + pow(EA23,2))*Kfactor23);
      else if(term==10) return ((1-fcSq) + fcSq*(1 + pow(EA24,2))*Kfactor24);
      else              return ((1-fcSq) + fcSq*(1 + pow(EA34,2))*Kfactor34);
    }else if(term==12){
      Float_t C22 = (1-fcSq) + fcSq*(1 + pow(EA12,2))*Kfactor12;
      C22 *= (1-fcSq) + fcSq*(1 + pow(EA34,2))*Kfactor34;
      return C22;
    }else return 1.0;
    
  }else{// mixed charge case
    if( ChargeSum==1 || ChargeSum==3){// ---+ and -+++ configuration
      Float_t EA1=0, EA2=0, EA3=0, Kpair1=0, Kpair2=0, Kpair3=0, Kpair4=0, Kpair5=0, Kpair6=0;
      Int_t c_OddOneOut = 1;
      if(ChargeSum==3) c_OddOneOut = 0;
      //
      if(c[0]==c_OddOneOut) {EA1=EA23; EA2=EA24; EA3=EA34; Kpair1=Kfactor23; Kpair2=Kfactor24; Kpair3=Kfactor34;   Kpair4=Kfactor12; Kpair5=Kfactor13; Kpair6=Kfactor14;}
      else if(c[1]==c_OddOneOut) {EA1=EA13; EA2=EA14; EA3=EA34; Kpair1=Kfactor13; Kpair2=Kfactor14; Kpair3=Kfactor34;   Kpair4=Kfactor12; Kpair5=Kfactor23; Kpair6=Kfactor24;}
      else if(c[2]==c_OddOneOut) {EA1=EA12; EA2=EA14; EA3=EA24; Kpair1=Kfactor12; Kpair2=Kfactor14; Kpair3=Kfactor24;   Kpair4=Kfactor13; Kpair5=Kfactor23; Kpair6=Kfactor34;}
      else {EA1=EA12; EA2=EA13; EA3=EA23; Kpair1=Kfactor12; Kpair2=Kfactor13; Kpair3=Kfactor23;   Kpair4=Kfactor14; Kpair5=Kfactor24; Kpair6=Kfactor34;}
      
      if(term==1){
	Float_t C3QS = 1 + pow(EA1,2) + pow(EA2,2) + pow(EA3,2) + 2*EA1*EA2*EA3;
	Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
	C4 += pow(fc,2)*pow(1-fc,2)*( (1 + pow(EA1,2))*Kpair1 + (1 + pow(EA2,2))*Kpair2 + (1 + pow(EA3,2))*Kpair3 );
	C4 += pow(fc,2)*pow(1-fc,2)*( Kpair4 + Kpair5 + Kpair6 );
	C4 += pow(fc,3)*(1-fc)*(1 + pow(EA1,2) + pow(EA2,2) +  pow(EA3,2) + 2*EA1*EA2*EA3) * Kpair1*Kpair2*Kpair3;
	C4 += pow(fc,3)*(1-fc)*( (1 + pow(EA1,2))*Kpair1*Kpair4*Kpair5 + (1+pow(EA2,2))*Kpair2*Kpair4*Kpair6 + (1+pow(EA3,2))*Kpair3*Kpair5*Kpair6);// doesn't matter which MC K's
	C4 += pow(fc,4)*C3QS*Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
	return C4;
      }else if( (term==2 && ChargeSum==1) || (term==5 && ChargeSum==3)){
	Float_t C3QS = 1 + pow(EA1,2) + pow(EA2,2) + pow(EA3,2) + 2*EA1*EA2*EA3;
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*(1+pow(EA1,2))*Kpair1;
	C3 += pow(fc,2)*(1-fc)*(1+pow(EA2,2))*Kpair2;
	C3 += pow(fc,2)*(1-fc)*(1+pow(EA3,2))*Kpair3;
	C3 += pow(fc,3)*C3QS*Kpair1*Kpair2*Kpair3;
	return C3;
      }else if(term<=5){// one SC pair, two MC pairs
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*(1+pow(EA1,2))*Kpair1;// any SC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair4;// any MC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair5;// any MC pair will do
	C3 += pow(fc,3)*(1+pow(EA1,2))*Kpair1*Kpair4*Kpair5;
	return C3;
      }else if(term==6 || term==7){
	if(ChargeSum==1) return ((1-fcSq) + fcSq*(1 + pow(EA1,2))*Kpair1);// any SC pair will do
	else return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else if(term==8){
	return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else if(term==9){
	return ((1-fcSq) + fcSq*(1 + pow(EA1,2))*Kpair1);// any SC pair will do
      }else if(term==10 || term==11){
	if(ChargeSum==3) return ((1-fcSq) + fcSq*(1 + pow(EA1,2))*Kpair1);// any SC pair will do
	else return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else return 1.0;// for 12 and 13
    }else{// --++ configuration
      Float_t EA1=0, EA2=0, Kpair1=0, Kpair2=0, Kpair3=0, Kpair4=0, Kpair5=0, Kpair6=0;
      if(c[0]==c[1]) {EA1=EA12; EA2=EA34; Kpair1=Kfactor12; Kpair2=Kfactor34;    Kpair3=Kfactor13; Kpair4=Kfactor14; Kpair5=Kfactor23; Kpair6=Kfactor24;}
      else if(c[0]==c[2]) {EA1=EA13; EA2=EA24; Kpair1=Kfactor13; Kpair2=Kfactor24;    Kpair3=Kfactor12; Kpair4=Kfactor14; Kpair5=Kfactor23; Kpair6=Kfactor34;}
      else {EA1=EA14; EA2=EA23; Kpair1=Kfactor14; Kpair2=Kfactor23;    Kpair3=Kfactor12; Kpair4=Kfactor13; Kpair5=Kfactor24; Kpair6=Kfactor34;}
      //
      if(term==1){
	Float_t C2QS = 1 + pow(EA1,2)*pow(EA2,2);
	Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
	C4 += pow(fc,2)*pow(1-fc,2)*( (1 + pow(EA1,2))*Kpair1 + (1 + pow(EA2,2))*Kpair2 );
	C4 += pow(fc,2)*pow(1-fc,2)*( Kpair3 + Kpair4 + Kpair5 + Kpair6 );
	C4 += pow(fc,3)*(1-fc)*( (1 + pow(EA1,2))*Kpair1*Kpair3*Kpair4 + (1 + pow(EA2,2))*Kpair2*Kpair3*Kpair4);
	C4 += pow(fc,3)*(1-fc)*( (1 + pow(EA1,2))*Kpair1*Kpair5*Kpair6 + (1 + pow(EA2,2))*Kpair2*Kpair5*Kpair6);// doesn't matter which two MC K's used
	C4 += pow(fc,4)*C2QS*Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
	return C4;
      }else if(term<=5){
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*(1+pow(EA1,2))*Kpair1;// any SC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair4;// any MC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair6;// any MC pair will do
	C3 += pow(fc,3)*(1+pow(EA1,2))*Kpair1*Kpair4*Kpair6;
	return C3;
      }else if(term==6 || term==11){
	return ((1-fcSq) + fcSq*(1 + pow(EA1,2))*Kpair1);// any SC pair will do
      }else if(term!=12 && term !=13){
	return ((1-fcSq) + fcSq*Kpair3);// any MC pair will do
      }else if(term==12){
	Float_t C22 = (1-fcSq) + fcSq*(1 + pow(EA1,2))*Kpair1;
	C22 *= (1-fcSq) + fcSq*(1 + pow(EA2,2))*Kpair2;
	return C22;
      }else return 1.0;
    }
  }
  
}
//________________________________________________________________________
Float_t AliFourPion::MCWeightFSI4(Int_t term, Float_t R, Float_t fcSq, Int_t c[4], Float_t qinv[6]){
  if(term==13) return 1.0;
    
  Int_t ChargeSum=c[0]+c[1]+c[2]+c[3];
 
  Float_t fc = sqrt(fcSq);
  SetFSIindex(R);
  Float_t Kfactor12 = FSICorrelation(c[0],c[1], qinv[0]);// K2
  Float_t Kfactor13 = FSICorrelation(c[0],c[2], qinv[1]);// K2
  Float_t Kfactor14 = FSICorrelation(c[0],c[3], qinv[2]);// K2
  Float_t Kfactor23 = FSICorrelation(c[1],c[2], qinv[3]);// K2
  Float_t Kfactor24 = FSICorrelation(c[1],c[3], qinv[4]);// K2
  Float_t Kfactor34 = FSICorrelation(c[2],c[3], qinv[5]);// K2
  
  if(c[0]==c[1] && c[0]==c[2] && c[0]==c[3]){// ---- and ++++ configuration
    
    if(term==1){
      Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
      C4 += pow(fc,2)*pow(1-fc,2)*( Kfactor12 + Kfactor13 + Kfactor14 );
      C4 += pow(fc,2)*pow(1-fc,2)*( Kfactor23 + Kfactor24 + Kfactor34 );
      C4 += pow(fc,3)*(1-fc) * Kfactor12*Kfactor13*Kfactor23;
      C4 += pow(fc,3)*(1-fc) * Kfactor12*Kfactor14*Kfactor24;
      C4 += pow(fc,3)*(1-fc) * Kfactor13*Kfactor14*Kfactor34;
      C4 += pow(fc,3)*(1-fc) * Kfactor23*Kfactor24*Kfactor34;
      C4 += pow(fc,4) * Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
      return C4;
    }else if(term<=5){
      Float_t Kpair1=0, Kpair2=0, Kpair3=0;
      if(term==2) {Kpair1=Kfactor12; Kpair2=Kfactor13; Kpair3=Kfactor23;}
      else if(term==3) {Kpair1=Kfactor12; Kpair2=Kfactor14; Kpair3=Kfactor24;}
      else if(term==4) {Kpair1=Kfactor13; Kpair2=Kfactor14; Kpair3=Kfactor34;}
      else {Kpair1=Kfactor23; Kpair2=Kfactor24; Kpair3=Kfactor34;}
      Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      C3 += pow(fc,2)*(1-fc)*( Kpair1 + Kpair2 + Kpair3 );
      C3 += pow(fc,3)*Kpair1*Kpair2*Kpair3;
      return C3;
    }else if(term<=11){
      if(term==6) return ((1-fcSq) + fcSq*Kfactor12);
      else if(term==7) return ((1-fcSq) + fcSq*Kfactor13);
      else if(term==8) return ((1-fcSq) + fcSq*Kfactor14);
      else if(term==9) return ((1-fcSq) + fcSq*Kfactor23);
      else if(term==10) return ((1-fcSq) + fcSq*Kfactor24);
      else return ((1-fcSq) + fcSq*Kfactor34);
    }else if(term==12){
      Float_t C22 = (1-fcSq) + fcSq*Kfactor12;
      C22 *= (1-fcSq) + fcSq*Kfactor34;
      return C22;
    }else return 1.0;
    
  }else{// mixed charge case
    if( ChargeSum==1 || ChargeSum==3){// ---+ and -+++ configuration
      Float_t Kpair1=0, Kpair2=0, Kpair3=0, Kpair4=0, Kpair5=0, Kpair6=0;
      Int_t c_OddOneOut = 1;
      if(ChargeSum==3) c_OddOneOut = 0;
      //
      if(c[0]==c_OddOneOut) {Kpair1=Kfactor23; Kpair2=Kfactor24; Kpair3=Kfactor34;   Kpair4=Kfactor12; Kpair5=Kfactor13; Kpair6=Kfactor14;}
      else if(c[1]==c_OddOneOut) {Kpair1=Kfactor13; Kpair2=Kfactor14; Kpair3=Kfactor34;   Kpair4=Kfactor12; Kpair5=Kfactor23; Kpair6=Kfactor24;}
      else if(c[2]==c_OddOneOut) {Kpair1=Kfactor12; Kpair2=Kfactor14; Kpair3=Kfactor24;   Kpair4=Kfactor13; Kpair5=Kfactor23; Kpair6=Kfactor34;}
      else {Kpair1=Kfactor12; Kpair2=Kfactor13; Kpair3=Kfactor23;   Kpair4=Kfactor14; Kpair5=Kfactor24; Kpair6=Kfactor34;}
      
      if(term==1){
	Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
	C4 += pow(fc,2)*pow(1-fc,2)*( Kpair1 + Kpair2 + Kpair3 );
	C4 += pow(fc,2)*pow(1-fc,2)*( Kpair4 + Kpair5 + Kpair6 );
	C4 += pow(fc,3)*(1-fc)*Kpair1*Kpair2*Kpair3;
	C4 += pow(fc,3)*(1-fc)*(Kpair1*Kpair4*Kpair5 + Kpair2*Kpair4*Kpair6 + Kpair3*Kpair5*Kpair6);// doesn't matter which two MC K's used
	C4 += pow(fc,4)*Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
	return C4;
      }else if( (term==2 && ChargeSum==1) || (term==5 && ChargeSum==3)){
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*Kpair1;
	C3 += pow(fc,2)*(1-fc)*Kpair2;
	C3 += pow(fc,2)*(1-fc)*Kpair3;
	C3 += pow(fc,3)*Kpair1*Kpair2*Kpair3;
	return C3;
      }else if(term<=5){// one SC pair, two MC pairs
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*Kpair1;// any SC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair4;// any MC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair5;// any MC pair will do
	C3 += pow(fc,3)*Kpair1*Kpair4*Kpair5;
	return C3;
      }else if(term==6 || term==7){
	if(ChargeSum==1) return ((1-fcSq) + fcSq*Kpair1);// any SC pair will do
	else return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else if(term==8){
	return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else if(term==9){
	return ((1-fcSq) + fcSq*Kpair1);// any SC pair will do
      }else if(term==10 || term==11){
	if(ChargeSum==3) return ((1-fcSq) + fcSq*Kpair1);// any SC pair will do
	else return ((1-fcSq) + fcSq*Kpair4);// any MC pair will do
      }else return 1.0;// 12 and 13
    }else{// --++ configuration
      Float_t Kpair1=0, Kpair2=0, Kpair3=0, Kpair4=0, Kpair5=0, Kpair6=0;
      if(c[0]==c[1]) {Kpair1=Kfactor12; Kpair2=Kfactor34;    Kpair3=Kfactor13; Kpair4=Kfactor14; Kpair5=Kfactor23; Kpair6=Kfactor24;}
      else if(c[0]==c[2]) {Kpair1=Kfactor13; Kpair2=Kfactor24;    Kpair3=Kfactor12; Kpair4=Kfactor14; Kpair5=Kfactor23; Kpair6=Kfactor34;}
      else {Kpair1=Kfactor14; Kpair2=Kfactor23;    Kpair3=Kfactor12; Kpair4=Kfactor13; Kpair5=Kfactor24; Kpair6=Kfactor34;}
      //
      if(term==1){
	Float_t C4 = pow(1-fc,4) + 4*fc*pow(1-fc,3);
	C4 += pow(fc,2)*pow(1-fc,2)*( Kpair1 + Kpair2 + Kpair3 + Kpair4 + Kpair5 + Kpair6);
	C4 += pow(fc,3)*(1-fc)*( Kpair1*Kpair3*Kpair4 + Kpair2*Kpair3*Kpair4 + Kpair1*Kpair5*Kpair6 + Kpair2*Kpair5*Kpair6);
	C4 += pow(fc,4)*Kfactor12*Kfactor13*Kfactor14*Kfactor23*Kfactor24*Kfactor34;
	return C4;
      }else if(term<=5){
	Float_t C3 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
	C3 += pow(fc,2)*(1-fc)*Kpair1;// any SC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair4;// any MC pair will do
	C3 += pow(fc,2)*(1-fc)*Kpair6;// any MC pair will do
	C3 += pow(fc,3)*Kpair1*Kpair4*Kpair6;
	return C3;
      }else if(term==6 || term==11){
	return ((1-fcSq) + fcSq*Kpair1);// any SC pair will do
      }else if(term !=12 && term !=13){
	return ((1-fcSq) + fcSq*Kpair3);// any MC pair will do
      }else if(term==12){
	Float_t C22 = (1-fcSq) + fcSq*Kpair1;
	C22 *= (1-fcSq) + fcSq*Kpair2;
	return C22;
      }else return 1.0;
    }
  }
  
}
//________________________________________________________________________
void AliFourPion::SetMomResCorrections(Bool_t legoCase, TH2D *temp2DSC, TH2D *temp2DMC){
  
 
  if(legoCase){
    cout<<"LEGO call to SetMomResCorrections"<<endl;
    fMomResC2SC = (TH2D*)temp2DSC->Clone();
    fMomResC2SC->SetDirectory(0);
    fMomResC2MC = (TH2D*)temp2DMC->Clone();
    fMomResC2MC->SetDirectory(0);
  }else {
    TFile *momResFile = new TFile("MomResFile.root","READ");
    if(!momResFile->IsOpen()) {
      cout<<"No momentum resolution file found"<<endl;
      AliFatal("No momentum resolution file found.  Kill process.");
    }else {cout<<"Good Momentum Resolution File Found!"<<endl;}
    
    TH2D *temp2DSC2 = (TH2D*)momResFile->Get("MRC_C2_SC");
    fMomResC2SC = (TH2D*)temp2DSC2->Clone();
    fMomResC2SC->SetDirectory(0);
    //
    TH2D *temp2DMC2 = (TH2D*)momResFile->Get("MRC_C2_MC");
    fMomResC2MC = (TH2D*)temp2DMC2->Clone();
    fMomResC2MC->SetDirectory(0);
    //
    momResFile->Close();
  }

  
  for(Int_t bx=1; bx<=fMomResC2SC->GetNbinsX(); bx++){
    for(Int_t by=1; by<=fMomResC2SC->GetNbinsY(); by++){
      if(fMomResC2SC->GetBinContent(bx,by) > 1.5) fMomResC2SC->SetBinContent(bx,by, 1.0);// Maximum is ~1.02 
      if(fMomResC2SC->GetBinContent(bx,by) < 0.8) fMomResC2SC->SetBinContent(bx,by, 1.0);// Minimum is ~0.8
      if(fMomResC2MC->GetBinContent(bx,by) > 1.5) fMomResC2MC->SetBinContent(bx,by, 1.0);// Maximum is ~1.02 
      if(fMomResC2MC->GetBinContent(bx,by) < 0.8) fMomResC2MC->SetBinContent(bx,by, 1.0);// Minimum is ~0.8
    }
  }
  
  cout<<"Done reading momentum resolution file"<<endl;
}
//________________________________________________________________________
void AliFourPion::SetFSICorrelations(Bool_t legoCase, TH1D *tempss[13], TH1D *tempos[13]){
  // read in 2-particle and 3-particle FSI correlations = K2 & K3
  // 2-particle input histo from file is binned in qinv.  3-particle in qinv of each pair
  if(legoCase){
    cout<<"LEGO call to SetFSICorrelations"<<endl;
    for(Int_t MB=0; MB<13; MB++) {
      fFSIss[MB] = (TH1D*)tempss[MB]->Clone();
      fFSIos[MB] = (TH1D*)tempos[MB]->Clone();
      //
      fFSIss[MB]->SetDirectory(0);
      fFSIos[MB]->SetDirectory(0);
    }
  }else {
    cout<<"non LEGO call to SetFSICorrelations"<<endl;
    TFile *fsifile = new TFile("KFile.root","READ");
    if(!fsifile->IsOpen()) {
      cout<<"No FSI file found"<<endl;
      AliFatal("No FSI file found.  Kill process.");
    }else {cout<<"Good FSI File Found!"<<endl;}
    
    TH1D *temphistoSS[13];
    TH1D *temphistoOS[13];
    for(Int_t MB=0; MB<13; MB++) {
      TString *nameK2SS = new TString("K2ss_");
      *nameK2SS += MB;
      temphistoSS[MB] = (TH1D*)fsifile->Get(nameK2SS->Data());
      //
      TString *nameK2OS = new TString("K2os_");
      *nameK2OS += MB;
      temphistoOS[MB] = (TH1D*)fsifile->Get(nameK2OS->Data());
      //
      fFSIss[MB] = (TH1D*)temphistoSS[MB]->Clone();
      fFSIos[MB] = (TH1D*)temphistoOS[MB]->Clone();
      fFSIss[MB]->SetDirectory(0);
      fFSIos[MB]->SetDirectory(0);
    }
    //
    
    fsifile->Close();
  }
  
  cout<<"Done reading FSI file"<<endl;
}
//________________________________________________________________________
Float_t AliFourPion::FSICorrelation(Int_t charge1, Int_t charge2, Float_t qinv){
  // returns 2-particle Coulomb correlations = K2
  Int_t qbinL = fFSIss[fFSIindex]->GetXaxis()->FindBin(qinv-fFSIss[fFSIindex]->GetXaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSIss[fFSIindex]->GetNbinsX()) {
    if(charge1!=charge2) {
      Float_t ScaleFac = (fFSIos[fFSIindex]->GetBinContent(fFSIos[fFSIindex]->GetNbinsX()-1) - 1);
      ScaleFac /= (Gamov(charge1, charge2, fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(fFSIos[fFSIindex]->GetNbinsX()-1)) - 1);
      return ( (Gamov(charge1, charge2, qinv)-1)*ScaleFac + 1); 
    }else{
      Float_t ScaleFac = (fFSIss[fFSIindex]->GetBinContent(fFSIss[fFSIindex]->GetNbinsX()-1) - 1);
      ScaleFac /= (Gamov(charge1, charge2, fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(fFSIss[fFSIindex]->GetNbinsX()-1)) - 1);
      return ( (Gamov(charge1, charge2, qinv)-1)*ScaleFac + 1);
    }
  }
  
  Float_t slope=0;
  if(charge1==charge2){
    slope = fFSIss[fFSIindex]->GetBinContent(qbinL) - fFSIss[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSIss[fFSIindex]->GetBinContent(qbinL));
  }else {
    slope = fFSIos[fFSIindex]->GetBinContent(qbinL) - fFSIos[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSIos[fFSIindex]->GetBinContent(qbinL));
  }
}
//________________________________________________________________________
void AliFourPion::SetFillBins2(Int_t c1, Int_t c2, Int_t &b1, Int_t &b2){
  if((c1+c2)==1) {b1=0; b2=1;}// Re-assign to merge degenerate histos
  else {b1=c1; b2=c2;}
}
//________________________________________________________________________
void AliFourPion::SetFillBins3(Int_t c1, Int_t c2, Int_t c3, Short_t part, Int_t &b1, Int_t &b2, Int_t &b3, Bool_t &fill2, Bool_t &fill3, Bool_t &fill4){
    
  // "part" specifies which pair is from the same event.  Only relevant for terms 2-4 
  Bool_t seSS=kFALSE;
  if(part==1) {// default case (irrelevant for term 1 and term 5)
    if(c1==c2) seSS=kTRUE;
  }
  if(part==2){
    if(c1==c3) seSS=kTRUE;
  }
  
  
  // fill2, fill3, fill4 are only used for Cumulant Terms 2,3,4
  if( (c1+c2+c3)==1) {
    b1=0; b2=0; b3=1;// Re-assign to merge degenerate histos
    //
    if(seSS) fill2=kTRUE;
    else {fill3=kTRUE; fill4=kTRUE;}
    //
  }else if( (c1+c2+c3)==2) {
    b1=0; b2=1; b3=1;
    //
    if(!seSS) {fill2=kTRUE; fill3=kTRUE;}
    else fill4=kTRUE;
    //
  }else {
    b1=c1; b2=c2; b3=c3;
    fill2=kTRUE; fill3=kTRUE; fill4=kTRUE;
  }
  
}
//________________________________________________________________________
void AliFourPion::SetFillBins4(Int_t c1, Int_t c2, Int_t c3, Int_t c4, Int_t &b1, Int_t &b2, Int_t &b3, Int_t &b4, Int_t ENsum, Bool_t fillTerm[13]){
  
  // fill2, fill3, fill4 are only used for Cumulant Terms 2,3,4
  if( (c1+c2+c3+c4)==0 || (c1+c2+c3+c4)==4) {// all of the same charge: ---- or ++++
    
    b1=c1; b2=c2; b3=c3; b4=c4;
    if(ENsum==0) fillTerm[0]=kTRUE;
    else if(ENsum==1) {fillTerm[1]=kTRUE; fillTerm[2]=kTRUE; fillTerm[3]=kTRUE; fillTerm[4]=kTRUE;}
    else if(ENsum==2) {fillTerm[11]=kTRUE;}
    else if(ENsum==3) {fillTerm[5]=kTRUE; fillTerm[6]=kTRUE; fillTerm[7]=kTRUE; fillTerm[8]=kTRUE; fillTerm[9]=kTRUE; fillTerm[10]=kTRUE;}
    else fillTerm[12]=kTRUE;
  
  }else if( (c1+c2+c3+c4)==1) {// one positive charge: ---+
  
    b1=0; b2=0; b3=0; b4=1;// Re-assign to merge degenerate histos
    if(ENsum==0) fillTerm[0]=kTRUE;
    else if(ENsum==1){
      if(c4==1) fillTerm[1]=kTRUE;
      else {fillTerm[2]=kTRUE; fillTerm[3]=kTRUE; fillTerm[4]=kTRUE;}
    }else if(ENsum==2){
      fillTerm[11]=kTRUE;
    }else if(ENsum==3){
      if(c3==1 || c4==1) {fillTerm[5]=kTRUE; fillTerm[6]=kTRUE; fillTerm[8]=kTRUE;} 
      else {fillTerm[7]=kTRUE; fillTerm[9]=kTRUE; fillTerm[10]=kTRUE;}
    }else fillTerm[12]=kTRUE;
  
  }else if( (c1+c2+c3+c4)==2) {// two positive charges: --++
    
    b1=0; b2=0; b3=1; b4=1;// Re-assign to merge degenerate histos
    if(ENsum==0) fillTerm[0]=kTRUE;
    else if(ENsum==1){
      if(c4==1) {fillTerm[1]=kTRUE; fillTerm[2]=kTRUE;}
      else {fillTerm[3]=kTRUE; fillTerm[4]=kTRUE;}
    }else if(ENsum==2){
      if( (c1+c2)==0) fillTerm[11]=kTRUE;
    }else if(ENsum==3){
      if( (c1+c2)==0) fillTerm[5]=kTRUE;
      else if( (c1+c2)==1) {fillTerm[6]=kTRUE; fillTerm[7]=kTRUE; fillTerm[8]=kTRUE; fillTerm[9]=kTRUE;}
      else fillTerm[10]=kTRUE;
    }else fillTerm[12]=kTRUE;

  }else{// three positive charges
    
    b1=0; b2=1; b3=1; b4=1;// Re-assign to merge degenerate histos
    if(ENsum==0) fillTerm[0]=kTRUE;
    else if(ENsum==1){
      if(c4==0) fillTerm[4]=kTRUE;
      else {fillTerm[1]=kTRUE; fillTerm[2]=kTRUE; fillTerm[3]=kTRUE;}
    }else if(ENsum==2){
      fillTerm[11]=kTRUE;
    }else if(ENsum==3){
      if(c3==0 || c4==0) {fillTerm[8]=kTRUE; fillTerm[9]=kTRUE; fillTerm[10]=kTRUE;}
      else {fillTerm[5]=kTRUE; fillTerm[6]=kTRUE; fillTerm[7]=kTRUE;}
    }else fillTerm[12]=kTRUE;
  
  }
  
}
//________________________________________________________________________
void AliFourPion::SetFSIindex(Float_t R){
  if(!fMCcase){
    if(fCollisionType==0){
      if(fMbin==0) fFSIindex = 0;//0-5%
      else if(fMbin==1) fFSIindex = 1;//5-10%
      else if(fMbin<=3) fFSIindex = 2;//10-20%
      else if(fMbin<=5) fFSIindex = 3;//20-30%
      else if(fMbin<=7) fFSIindex = 4;//30-40%
      else if(fMbin<=9) fFSIindex = 5;//40-50%
      else if(fMbin<=12) fFSIindex = 6;//40-50%
      else if(fMbin<=15) fFSIindex = 7;//40-50%
      else if(fMbin<=18) fFSIindex = 8;//40-50%
      else fFSIindex = 8;//90-100%
    }else fFSIindex = fFSIindexSmallSystem;// pPb and pp
  }else{// FSI binning for MC 
    if(R>=10.) fFSIindex = 0;
    else if(R>=9.) fFSIindex = 1;
    else if(R>=8.) fFSIindex = 2;
    else if(R>=7.) fFSIindex = 3;
    else if(R>=6.) fFSIindex = 4;
    else if(R>=5.) fFSIindex = 5;
    else if(R>=4.) fFSIindex = 6;
    else if(R>=3.) fFSIindex = 7;
    else if(R>=2.) fFSIindex = 8;
    else fFSIindex = 9;
  }
}
//________________________________________________________________________
void AliFourPion::SetMuonCorrections(Bool_t legoCase, TH2D *tempMuon){
  if(legoCase){
    cout<<"LEGO call to SetMuonCorrections"<<endl;
    fWeightmuonCorrection = (TH2D*)tempMuon->Clone();
    fWeightmuonCorrection->SetDirectory(0);
  }else {
    cout<<"non LEGO call to SetMuonCorrections"<<endl;
    TFile *MuonFile=new TFile("MuonCorrection.root","READ");
    if(!MuonFile->IsOpen()) {
      cout<<"No Muon file found"<<endl;
      AliFatal("No Muon file found.  Kill process.");
    }else {cout<<"Good Muon File Found!"<<endl;}
    
    fWeightmuonCorrection = (TH2D*)MuonFile->Get("WeightmuonCorrection");
    fWeightmuonCorrection->SetDirectory(0);
    //
    MuonFile->Close();
  }
  cout<<"Done reading Muon file"<<endl;
}
//________________________________________________________________________
void AliFourPion::Setc3FitEAs(Bool_t legoCase, TH3D *histoPbPb[2], TH3D *histopPb[2], TH3D *histopp[2]){
  
  if(legoCase){
    cout<<"LEGO call to Setc3FitEAs"<<endl;
    for(Int_t FT=0; FT<2; FT++){
      fPbPbc3FitEA[FT] = (TH3D*)histoPbPb[FT]->Clone();
      fpPbc3FitEA[FT] = (TH3D*)histopPb[FT]->Clone();
      fppc3FitEA[FT] = (TH3D*)histopp[FT]->Clone();
      fPbPbc3FitEA[FT]->SetDirectory(0);
      fpPbc3FitEA[FT]->SetDirectory(0);
      fppc3FitEA[FT]->SetDirectory(0);
    }
  }else{
    cout<<"non LEGO call to Setc3FitEAs"<<endl;
    TFile *EAfile = new TFile("c3EAfile.root","READ");
    if(!EAfile->IsOpen()) {
      cout<<"No EA file found"<<endl;
      AliFatal("No EA file found.  Kill process.");
    }else {cout<<"Good EA File Found!"<<endl;}
    for(Int_t FT=0; FT<2; FT++){
      if(FT==0){
	fPbPbc3FitEA[FT] = (TH3D*)EAfile->Get("PbPbEA_c3");
	fpPbc3FitEA[FT] = (TH3D*)EAfile->Get("pPbEA_c3");
	fppc3FitEA[FT] = (TH3D*)EAfile->Get("ppEA_c3");
      }else{
	fPbPbc3FitEA[FT] = (TH3D*)EAfile->Get("PbPbEA_C3");
	fpPbc3FitEA[FT] = (TH3D*)EAfile->Get("pPbEA_C3");
	fppc3FitEA[FT] = (TH3D*)EAfile->Get("ppEA_C3");
      }
      fPbPbc3FitEA[FT]->SetDirectory(0);
      fpPbc3FitEA[FT]->SetDirectory(0);
      fppc3FitEA[FT]->SetDirectory(0);
    }
    EAfile->Close();
  }
  cout<<"Done reading EA file"<<endl;
}
//________________________________________________________________________
Float_t AliFourPion::cubicInterpolate (Float_t p[4], Float_t x) {
  return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));// Paulinternet
}
//________________________________________________________________________
Float_t AliFourPion::nCubicInterpolate (Int_t n, Float_t* p, Float_t coordinates[]) {
 
  if (n == 1) {
    return cubicInterpolate(p, *coordinates);
  }
  else {
    Float_t arr[4];
    Int_t skip = 1 << (n - 1) * 2;
    arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
    arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
    arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
    arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
    return cubicInterpolate(arr, *coordinates);
  }
}
