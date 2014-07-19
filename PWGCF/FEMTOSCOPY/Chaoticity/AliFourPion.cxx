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
  fPbPbcase(kTRUE),
  fGenerateSignal(kFALSE),
  fGeneratorOnly(kFALSE),
  fTabulatePairs(kFALSE),
  fLinearInterpolation(kTRUE),
  fMixedChargeCut(kFALSE),
  fRMax(11),
  ffcSq(0.7),
  ffcSqMRC(0.6),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fMinPt(0.16),
  fMaxPt(1.0),
  fQcut(0),
  fQLowerCut(0),
  fNormQcutLow(0),
  fNormQcutHigh(0),
  fKupperBound(0),
  fQupperBoundQ2(0),
  fQupperBoundQ3(0),
  fQupperBoundQ4(0),
  fQbinsQ2(1),
  fQbinsQ3(1),
  fQbinsQ4(1),
  fQupperBoundWeights(0),
  fKstepT(),
  fKstepY(),
  fKmeanT(),
  fKmeanY(),
  fKmiddleT(),
  fKmiddleY(),
  fQstep(0),
  fQstepWeights(0),
  fQmean(),
  fDampStart(0),
  fDampStep(0),
  fTPCTOFboundry(0),
  fTOFboundry(0),
  fSigmaCutTPC(2.0),
  fSigmaCutTOF(2.0),
  fMinSepPairEta(0.03),
  fMinSepPairPhi(0.04),
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
  fDummyB(0),
  fKT3transition(0.3),
  fKT4transition(0.3),
  farrP1(),
  farrP2(),
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
  fWeightmuonCorrection(0x0)
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
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm = 0x0;
	      	      
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm = 0x0;
		
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
  
  // Initialize FSI histograms
  for(Int_t i=0; i<12; i++){
    fFSIss[i]=0x0; 
    fFSIos[i]=0x0;
  }


  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
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
  fPbPbcase(kTRUE),
  fGenerateSignal(kFALSE),
  fGeneratorOnly(kFALSE),
  fTabulatePairs(kFALSE),
  fLinearInterpolation(kTRUE),
  fMixedChargeCut(kFALSE),
  fRMax(11),
  ffcSq(0.7),
  ffcSqMRC(0.6),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fMinPt(0.16),
  fMaxPt(1.0),
  fQcut(0),
  fQLowerCut(0),
  fNormQcutLow(0),
  fNormQcutHigh(0),
  fKupperBound(0),
  fQupperBoundQ2(0),
  fQupperBoundQ3(0),
  fQupperBoundQ4(0),
  fQbinsQ2(1),
  fQbinsQ3(1),
  fQbinsQ4(1),
  fQupperBoundWeights(0),
  fKstepT(),
  fKstepY(),
  fKmeanT(),
  fKmeanY(),
  fKmiddleT(),
  fKmiddleY(),
  fQstep(0),
  fQstepWeights(0),
  fQmean(),
  fDampStart(0),
  fDampStep(0),
  fTPCTOFboundry(0),
  fTOFboundry(0),
  fSigmaCutTPC(2.0),
  fSigmaCutTOF(2.0),
  fMinSepPairEta(0.03),
  fMinSepPairPhi(0.04),
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
  fDummyB(0),
  fKT3transition(0.3),
  fKT4transition(0.3),
  farrP1(),
  farrP2(),
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
  fWeightmuonCorrection(0x0)
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
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm = 0x0;
	      
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm = 0x0;
		
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
  
  // Initialize FSI histograms
  for(Int_t i=0; i<12; i++){
    fFSIss[i]=0x0; 
    fFSIos[i]=0x0;
  }
  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
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
    fPbPbcase(obj.fPbPbcase),
    fGenerateSignal(obj.fGenerateSignal),
    fGeneratorOnly(obj.fGeneratorOnly),
    fTabulatePairs(obj.fTabulatePairs),
    fLinearInterpolation(obj.fLinearInterpolation),
    fMixedChargeCut(obj.fMixedChargeCut),
    fRMax(obj.fRMax),
    ffcSq(obj.ffcSq),
    ffcSqMRC(obj.ffcSqMRC),
    fFilterBit(obj.fFilterBit),
    fMaxChi2NDF(obj.fMaxChi2NDF),
    fMinTPCncls(obj.fMinTPCncls),
    fBfield(obj.fBfield),
    fMbin(obj.fMbin),
    fFSIindex(obj.fFSIindex),
    fEDbin(obj.fEDbin),
    fMbins(obj.fMbins),
    fMultLimit(obj.fMultLimit),
    fCentBinLowLimit(obj.fCentBinLowLimit),
    fCentBinHighLimit(obj.fCentBinHighLimit),
    fEventCounter(obj.fEventCounter),
    fEventsToMix(obj.fEventsToMix),
    fZvertexBins(obj.fZvertexBins),
    fMultLimits(),
    fMinPt(obj.fMinPt),
    fMaxPt(obj.fMaxPt),
    fQcut(obj.fQcut),
    fQLowerCut(obj.fQLowerCut),
    fNormQcutLow(0),
    fNormQcutHigh(0),
    fKupperBound(obj.fKupperBound),
    fQupperBoundQ2(obj.fQupperBoundQ2),
    fQupperBoundQ3(obj.fQupperBoundQ3),
    fQupperBoundQ4(obj.fQupperBoundQ4),
    fQbinsQ2(obj.fQbinsQ2),
    fQbinsQ3(obj.fQbinsQ3),
    fQbinsQ4(obj.fQbinsQ4),
    fQupperBoundWeights(obj.fQupperBoundWeights),
    fKstepT(),
    fKstepY(),
    fKmeanT(),
    fKmeanY(),
    fKmiddleT(),
    fKmiddleY(),
    fQstep(obj.fQstep),
    fQstepWeights(obj.fQstepWeights),
    fQmean(),
    fDampStart(obj.fDampStart),
    fDampStep(obj.fDampStep),
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
    fDummyB(obj.fDummyB),
    fKT3transition(obj.fKT3transition),
    fKT4transition(obj.fKT4transition),
    farrP1(),
    farrP2(),
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
    fWeightmuonCorrection(obj.fWeightmuonCorrection)
{
  // Copy Constructor
  
  for(Int_t i=0; i<12; i++){
    fFSIss[i]=obj.fFSIss[i]; 
    fFSIos[i]=obj.fFSIos[i];
  }
  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
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
  fLEGO = fLEGO;
  fMCcase = obj.fMCcase;
  fAODcase = obj.fAODcase;
  fPbPbcase = obj.fPbPbcase; 
  fGenerateSignal = obj.fGenerateSignal;
  fGeneratorOnly = obj.fGeneratorOnly;
  fTabulatePairs = obj.fTabulatePairs;
  fLinearInterpolation = obj.fLinearInterpolation;
  fMixedChargeCut = obj.fMixedChargeCut;
  fRMax = obj.fRMax;
  ffcSq = obj.ffcSq;
  ffcSqMRC = obj.ffcSqMRC;
  fFilterBit = obj.fFilterBit;
  fMaxChi2NDF = obj.fMaxChi2NDF;
  fMinTPCncls = obj.fMinTPCncls;
  fBfield = obj.fBfield;
  fMbin = obj.fMbin;
  fFSIindex = obj.fFSIindex;
  fEDbin = obj.fEDbin;
  fMbins = obj.fMbins;
  fMultLimit = obj.fMultLimit;
  fCentBinLowLimit = obj.fCentBinLowLimit;
  fCentBinHighLimit = obj.fCentBinHighLimit;
  fEventCounter = obj.fEventCounter;
  fEventsToMix = obj.fEventsToMix;
  fZvertexBins = obj.fZvertexBins;
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
  fQstep = obj.fQstep;
  fQstepWeights = obj.fQstepWeights;
  fDampStart = obj.fDampStart;
  fDampStep = obj.fDampStep;
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
  fDummyB = obj.fDummyB;
  fKT3transition = obj.fKT3transition;
  fKT4transition = obj.fKT4transition;
  fMomResC2SC = obj.fMomResC2SC;
  fMomResC2MC = obj.fMomResC2MC;
  fWeightmuonCorrection = obj.fWeightmuonCorrection;
  
  for(Int_t i=0; i<12; i++){
    fFSIss[i]=obj.fFSIss[i]; 
    fFSIos[i]=obj.fFSIos[i];
  }
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=obj.fNormWeight[i][j];
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
	      if(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm) delete Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm;
		
	    }// term_3

	    for(Int_t c4=0; c4<2; c4++){
	      for(Int_t term=0; term<13; term++){
		
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fNorm4;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor;
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted;
		//
		if(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm) delete Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm;
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
  
   
  for(Int_t i=0; i<12; i++){
    if(fFSIss[i]) delete fFSIss[i]; 
    if(fFSIos[i]) delete fFSIos[i];
  }
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      if(fNormWeight[i][j]) delete fNormWeight[i][j];
    }
  }
  
}
//________________________________________________________________________
void AliFourPion::ParInit()
{
  cout<<"AliFourPion MyInit() call"<<endl;
  cout<<"lego:"<<fLEGO<<"  MCcase:"<<fMCcase<<"  PbPbcase:"<<fPbPbcase<<"  TabulatePairs:"<<fTabulatePairs<<"  GenSignal:"<<fGenerateSignal<<"  CentLow:"<<fCentBinLowLimit<<"  CentHigh:"<<fCentBinHighLimit<<"  RMax:"<<fRMax<<"  fc^2:"<<ffcSq<<"  FB:"<<fFilterBit<<"  MaxChi2/NDF:"<<fMaxChi2NDF<<"  MinTPCncls:"<<fMinTPCncls<<"  MinPairSepEta:"<<fMinSepPairEta<<"  MinPairSepPhi:"<<fMinSepPairPhi<<"  NsigTPC:"<<fSigmaCutTPC<<"  NsigTOF:"<<fSigmaCutTOF<<endl;

  fRandomNumber = new TRandom3();
  fRandomNumber->SetSeed(0);
    
  //
  fEventCounter=0;
  fEventsToMix=3;
  fZvertexBins=2;//2
  
  fTPCTOFboundry = 0.6;// TPC pid used below this momentum, TOF above but below TOF_boundry
  fTOFboundry = 2.1;// TOF pid used below this momentum
  
  ////////////////////////////////////////////////
  // PadRow Pair Cuts
  fShareQuality = .5;// max
  fShareFraction = .05;// max
  ////////////////////////////////////////////////
  
  
  fMultLimits[0]=0, fMultLimits[1]=2, fMultLimits[2]=4, fMultLimits[3]=6, fMultLimits[4]=8, fMultLimits[5]=10;
  fMultLimits[6]=12, fMultLimits[7]=14, fMultLimits[8]=16, fMultLimits[9]=18, fMultLimits[10]=20, fMultLimits[11]=150;
  
    
  
  if(fPbPbcase) {// PbPb
    fMultLimit=kMultLimitPbPb;
    fMbins=fCentBins;
    fQcut=0.1;
    fNormQcutLow = 0.15;// 0.15
    fNormQcutHigh = 0.2;// 0.175
  }else {// pp
    fMultLimit=kMultLimitpp; 
    fMbins=kMultBinspp; 
    fQcut=0.6;
    fNormQcutLow = 1.0;
    fNormQcutHigh = 1.5;
  }
  
  fQLowerCut = 0.005;// was 0.005
  fKupperBound = 1.0;
  //
  fKstepY[0] = 1.6;
  fKmeanY[0] = 0;// central y
  fKmiddleY[0] = 0;

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
  fQupperBoundWeights = 0.2;
  fQupperBoundQ2 = 2.0;
  fQupperBoundQ3 = 0.6;
  fQupperBoundQ4 = 0.6;
  fQbinsQ2 = fQupperBoundQ2/0.005;
  fQbinsQ3 = fQupperBoundQ3/0.005;
  fQbinsQ4 = fQupperBoundQ4/0.005;
  fQstepWeights = fQupperBoundWeights/Float_t(kQbinsWeights);
  for(Int_t i=0; i<kQbinsWeights; i++) {fQmean[i]=(i+0.5)*fQstepWeights;}
  //
  fDampStart = 0.5;// was 0.3, then 0.5
  fDampStep = 0.02;
  
  //
  
  fEC = new AliFourPionEventCollection **[fZvertexBins];
  for(UShort_t i=0; i<fZvertexBins; i++){
    
    fEC[i] = new AliFourPionEventCollection *[fMbins];

    for(UShort_t j=0; j<fMbins; j++){
      
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

  

  // Set weights, Coulomb corrections, and Momentum resolution corrections manually if not on LEGO
  if(!fLEGO) {
    SetFSICorrelations(fLEGO);// Read in 2-particle and 3-particle FSI correlations
    if(!fTabulatePairs) SetWeightArrays(fLEGO);// Set Weight Array
    if(!fMCcase && !fTabulatePairs) SetMomResCorrections(fLEGO);// Read Momentum resolution file
    if(!fMCcase && !fTabulatePairs) SetMuonCorrections(fLEGO);// Read Muon corrections
  }
  
  /////////////////////////////////////////////
  /////////////////////////////////////////////
  
}
//________________________________________________________________________
void AliFourPion::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  ParInit();// Initialize my settings


  fOutputList = new TList();
  fOutputList->SetOwner();
  
  TH3F *fVertexDist = new TH3F("fVertexDist","Vertex Distribution",20,-1,1, 20,-1,1, 600,-30,30);
  fVertexDist->GetXaxis()->SetTitle("X Vertex (cm)");
  fVertexDist->GetYaxis()->SetTitle("Y Vertex (cm)");
  fVertexDist->GetZaxis()->SetTitle("Z Vertex (cm)");
  fOutputList->Add(fVertexDist);
  
  
  TH2F *fDCAxyDistPlus = new TH2F("fDCAxyDistPlus","DCA distribution",300,0,3., 50,0,5);
  fOutputList->Add(fDCAxyDistPlus);
  TH2F *fDCAzDistPlus = new TH2F("fDCAzDistPlus","DCA distribution",300,0,3., 50,0,5);
  fOutputList->Add(fDCAzDistPlus);
  TH2F *fDCAxyDistMinus = new TH2F("fDCAxyDistMinus","DCA distribution",300,0,3., 50,0,5);
  fOutputList->Add(fDCAxyDistMinus);
  TH2F *fDCAzDistMinus = new TH2F("fDCAzDistMinus","DCA distribution",300,0,3., 50,0,5);
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
  
  TH3F *fPtEtaDist = new TH3F("fPtEtaDist","fPtEtaDist",2,-1.1,1.1, 300,0,3., 28,-1.4,1.4);
  fOutputList->Add(fPtEtaDist);

  TH3F *fPhiPtDist = new TH3F("fPhiPtDist","fPhiPtDist",2,-1.1,1.1, 120,0,2*PI, 300,0,3.);
  fOutputList->Add(fPhiPtDist);
  
  TH3F *fTOFResponse = new TH3F("fTOFResponse","TOF relative time",20,0,100, 200,0,2, 4000,-20000,20000);
  fOutputList->Add(fTOFResponse);
  TH3F *fTPCResponse = new TH3F("fTPCResponse","TPCsignal",20,0,100, 200,0,2, 1000,0,1000);
  fOutputList->Add(fTPCResponse);
 
  TH1F *fRejectedPairs = new TH1F("fRejectedPairs","",400,0,2);
  fOutputList->Add(fRejectedPairs);
  TH1F *fRejectedPairsWeighting = new TH1F("fAcceptedPairsWeighting","",400,0,2);
  fOutputList->Add(fRejectedPairsWeighting);
  TH1F *fTotalPairsWeighting = new TH1F("fTotalPairsWeighting","",400,0,2);
  fOutputList->Add(fTotalPairsWeighting);
  //
  TH1F *fRejectedPairsMC = new TH1F("fRejectedPairsMC","",400,0,2);
  fOutputList->Add(fRejectedPairsMC);
  TH1F *fRejectedPairsWeightingMC = new TH1F("fAcceptedPairsWeightingMC","",400,0,2);
  fOutputList->Add(fRejectedPairsWeightingMC);
  TH1F *fTotalPairsWeightingMC = new TH1F("fTotalPairsWeightingMC","",400,0,2);
  fOutputList->Add(fTotalPairsWeightingMC);
  
  TH1I *fRejectedEvents = new TH1I("fRejectedEvents","",fMbins,0.5,fMbins+.5);
  fOutputList->Add(fRejectedEvents);
    
  TH3F *fPairsDetaDPhiNum = new TH3F("fPairsDetaDPhiNum","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiNum);
  TH3F *fPairsDetaDPhiDen = new TH3F("fPairsDetaDPhiDen","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiDen);
  TH3F *fPairsShareFracDPhiNum = new TH3F("fPairsShareFracDPhiNum","",10,-.5,9.5, 159,0,1, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiNum);
  TH3F *fPairsShareFracDPhiDen = new TH3F("fPairsShareFracDPhiDen","",10,-.5,9.5, 159,0,1, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiDen);
  TH3D* fPairsPadRowNum = new TH3D("fPairsPadRowNum","", 20,0,1, 159,0,1, 40,0,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowNum);
  TH3D* fPairsPadRowDen = new TH3D("fPairsPadRowDen","", 20,0,1, 159,0,1, 40,0,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowDen);



  TH2D *fResonanceOSPairs = new TH2D("fResonanceOSPairs","",fMbins,.5,fMbins+.5, 1000,0,2);
  if(fMCcase) fOutputList->Add(fResonanceOSPairs);
  TH2D *fAllOSPairs = new TH2D("fAllOSPairs","",fMbins,.5,fMbins+.5, 1000,0,2);
  if(fMCcase) fOutputList->Add(fAllOSPairs);
  
  TH3D *fPrimarySCPionPairs = new TH3D("fPrimarySCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fPrimarySCPionPairs);
  TH3D *fAllSCPionPairs = new TH3D("fAllSCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fAllSCPionPairs);
  TH3D *fPrimaryMCPionPairs = new TH3D("fPrimaryMCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fPrimaryMCPionPairs);
  TH3D *fAllMCPionPairs = new TH3D("fAllMCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fAllMCPionPairs);

  //
  TH1D *fMuonParents = new TH1D("fMuonParents","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fMuonParents);
  TH1D *fSecondaryMuonParents = new TH1D("fSecondaryMuonParents","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fSecondaryMuonParents);
  TH3D *fMuonPionDeltaQinv = new TH3D("fMuonPionDeltaQinv","",2,-0.5,1.5, 20,0,1, 100,-0.2,0.2);
  if(fMCcase) fOutputList->Add(fMuonPionDeltaQinv);
  TH1D *fPionCandidates = new TH1D("fPionCandidates","",500,0.5,500.5);
  if(fMCcase) fOutputList->Add(fPionCandidates);
  //  
  

  TProfile *fAvgMult = new TProfile("fAvgMult","",fMbins,.5,fMbins+.5, 0,1500,"");
  fOutputList->Add(fAvgMult);

  TH2D *fTrackChi2NDF = new TH2D("fTrackChi2NDF","",20,0,100, 100,0,10);
  fOutputList->Add(fTrackChi2NDF);
  TH2D *fTrackTPCncls = new TH2D("fTrackTPCncls","",20,0,100, 110,50,160);
  fOutputList->Add(fTrackTPCncls);


  TH1D *fTPNRejects3pion1 = new TH1D("fTPNRejects3pion1","",fQbinsQ3,0,fQupperBoundQ3);
  fOutputList->Add(fTPNRejects3pion1);
  TH1D *fTPNRejects3pion2 = new TH1D("fTPNRejects3pion2","",fQbinsQ3,0,fQupperBoundQ3);
  fOutputList->Add(fTPNRejects3pion2);
  TH1D *fTPNRejects4pion1 = new TH1D("fTPNRejects4pion1","",fQbinsQ4,0,fQupperBoundQ4);
  fOutputList->Add(fTPNRejects4pion1);

  TH3D *fKT3DistTerm1 = new TH3D("fKT3DistTerm1","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  TH3D *fKT3DistTerm5 = new TH3D("fKT3DistTerm5","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  fOutputList->Add(fKT3DistTerm1);
  fOutputList->Add(fKT3DistTerm5);
  TH3D *fKT4DistTerm1 = new TH3D("fKT4DistTerm1","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  TH3D *fKT4DistTerm13 = new TH3D("fKT4DistTerm13","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  fOutputList->Add(fKT4DistTerm1);
  fOutputList->Add(fKT4DistTerm13);


  TProfile2D *fKT3AvgpT = new TProfile2D("fKT3AvgpT","",fMbins,.5,fMbins+.5, 2,-0.5,1.5, 0.,1.0,"");
  fOutputList->Add(fKT3AvgpT);
  TProfile2D *fKT4AvgpT = new TProfile2D("fKT4AvgpT","",fMbins,.5,fMbins+.5, 2,-0.5,1.5, 0.,1.0,"");
  fOutputList->Add(fKT4AvgpT);
  TH3D* fQ3AvgpT = new TH3D("fQ3AvgpT","", 2,-0.5,1.5, fQbinsQ3,0,fQupperBoundQ3, 180,0.1,1.0);
  fOutputList->Add(fQ3AvgpT);
  TH3D* fQ4AvgpT = new TH3D("fQ4AvgpT","", 2,-0.5,1.5, fQbinsQ4,0,fQupperBoundQ4, 180,0.1,1.0);
  fOutputList->Add(fQ4AvgpT);


  TH1D *fMCWeight3DTerm1SC = new TH1D("fMCWeight3DTerm1SC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm1SCden = new TH1D("fMCWeight3DTerm1SCden","", 20,0,0.2);
  TH1D *fMCWeight3DTerm2SC = new TH1D("fMCWeight3DTerm2SC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm2SCden = new TH1D("fMCWeight3DTerm2SCden","", 20,0,0.2);
  TH1D *fMCWeight3DTerm1MC = new TH1D("fMCWeight3DTerm1MC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm1MCden = new TH1D("fMCWeight3DTerm1MCden","", 20,0,0.2);
  TH1D *fMCWeight3DTerm2MC = new TH1D("fMCWeight3DTerm2MC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm2MCden = new TH1D("fMCWeight3DTerm2MCden","", 20,0,0.2);
  TH1D *fMCWeight3DTerm3MC = new TH1D("fMCWeight3DTerm3MC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm3MCden = new TH1D("fMCWeight3DTerm3MCden","", 20,0,0.2);
  TH1D *fMCWeight3DTerm4MC = new TH1D("fMCWeight3DTerm4MC","", 20,0,0.2);
  TH1D *fMCWeight3DTerm4MCden = new TH1D("fMCWeight3DTerm4MCden","", 20,0,0.2);
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
    if((mb < fCentBinLowLimit) || (mb > fCentBinHighLimit)) continue;
    
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
	    
	    
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2 = new TH2D(nameEx2->Data(),"Two Particle Distribution",20,0,1, fQbinsQ2,0,fQupperBoundQ2);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2);
	    TString *nameEx2QW=new TString(nameEx2->Data());
	    nameEx2QW->Append("_QW");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2QW = new TH2D(nameEx2QW->Data(),"Two Particle Distribution",20,0,1, fQbinsQ2,0,fQupperBoundQ2);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fTerms2QW);
	    TString *nameAvgP=new TString(nameEx2->Data());
	    nameAvgP->Append("_AvgP");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fAvgP = new TProfile2D(nameAvgP->Data(),"",10,0,1, fQbinsQ2,0,fQupperBoundQ2, 0.,1.0,"");
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fAvgP);
	    
	    TString *nameUnitMult=new TString(nameEx2->Data());
	    nameUnitMult->Append("_UnitMult");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fUnitMultBin = new TH2D(nameUnitMult->Data(),"Two Particle Distribution",21,0.5,21.5, fQbinsQ2,0,fQupperBoundQ2);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fUnitMultBin);
	    
	    if(fMCcase){
	      // Momentum resolution histos
	      TString *nameIdeal = new TString(nameEx2->Data());
	      nameIdeal->Append("_Ideal");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal = new TH2D(nameIdeal->Data(),"Two Particle Distribution",11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fIdeal);
	      TString *nameSmeared = new TString(nameEx2->Data());
	      nameSmeared->Append("_Smeared");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared = new TH2D(nameSmeared->Data(),"Two Particle Distribution",11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fSmeared);
	      //
	      // Muon correction histos
	      TString *nameMuonIdeal=new TString(nameEx2->Data());
	      nameMuonIdeal->Append("_MuonIdeal");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonIdeal = new TH2D(nameMuonIdeal->Data(),"", 11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonIdeal);
	      TString *nameMuonSmeared=new TString(nameEx2->Data());
	      nameMuonSmeared->Append("_MuonSmeared");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonSmeared = new TH2D(nameMuonSmeared->Data(),"", 11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonSmeared);
	      //
	      TString *nameMuonPionK2=new TString(nameEx2->Data());
	      nameMuonPionK2->Append("_MuonPionK2");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonPionK2 = new TH2D(nameMuonPionK2->Data(),"", 11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMuonPionK2);
	      //
	      TString *namePionPionK2=new TString(nameEx2->Data());
	      namePionPionK2->Append("_PionPionK2");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPionPionK2 = new TH2D(namePionPionK2->Data(),"", 11,0.5,11.5, fQbinsQ2,0,fQupperBoundQ2);
	      if(mb==0 && edB==0) fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPionPionK2);
	      //
	      //
	      TString *nameEx2MC=new TString(nameEx2->Data());
	      nameEx2MC->Append("_MCqinv");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv = new TH1D(nameEx2MC->Data(),"", fQbinsQ2,0,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinv);
	      TString *nameEx2MCQW=new TString(nameEx2->Data());
	      nameEx2MCQW->Append("_MCqinvQW");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = new TH1D(nameEx2MCQW->Data(),"", fQbinsQ2,0,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW);
	      //
	      TString *nameEx2PIDpurityDen=new TString(nameEx2->Data());
	      nameEx2PIDpurityDen->Append("_PIDpurityDen");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = new TH2D(nameEx2PIDpurityDen->Data(),"Two Particle Distribution",20,0,1, fQbinsQ2,0,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen);
	      TString *nameEx2PIDpurityNum=new TString(nameEx2->Data());
	      nameEx2PIDpurityNum->Append("_PIDpurityNum");
	      Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = new TH3D(nameEx2PIDpurityNum->Data(),"Two Particle Distribution",16,0.5,16.5, 20,0,1, fQbinsQ2,0,fQupperBoundQ2);
	      fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum);
	    }
	    TString *nameEx2OSLB1 = new TString(nameEx2->Data()); 
	    nameEx2OSLB1->Append("_osl_b1");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSL);
	    nameEx2OSLB1->Append("_QW");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fTerms2OSLQW);
	    //
	    TString *nameEx2OSLB2 = new TString(nameEx2->Data()); 
	    nameEx2OSLB2->Append("_osl_b2");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSL);
	    nameEx2OSLB2->Append("_QW");
	    Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	    fOutputList->Add(Charge1[c1].Charge2[c2].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fTerms2OSLQW);
	    
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
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3 = new TH1D(name1DQ->Data(),"", fQbinsQ3,0,fQupperBoundQ3);
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTerms3);
	      //
	      TString *nameKfactor=new TString(namePC3->Data());
	      nameKfactor->Append("_Kfactor");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor = new TProfile(nameKfactor->Data(),"", fQbinsQ3,0,fQupperBoundQ3, 0,100, "");
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactor);
	      //
	      TString *nameKfactorW=new TString(namePC3->Data());
	      nameKfactorW->Append("_KfactorWeighted");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted = new TProfile(nameKfactorW->Data(),"", fQbinsQ3,0,fQupperBoundQ3, 0,100, "");
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fKfactorWeighted);
	      //
	      TString *nameMeanQinv=new TString(namePC3->Data());
	      nameMeanQinv->Append("_MeanQinv");
	      Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMeanQinv = new TProfile(nameMeanQinv->Data(),"", fQbinsQ3,0,fQupperBoundQ3, 0,.2, "");
	      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMeanQinv);
	      
	      if(fMCcase==kTRUE){
		// Momentum resolution correction histos
		TString *nameMomResIdeal=new TString(namePC3->Data());
		nameMomResIdeal->Append("_Ideal");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fIdeal = new TH2D(nameMomResIdeal->Data(),"", 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fIdeal);
		TString *nameMomResSmeared=new TString(namePC3->Data());
		nameMomResSmeared->Append("_Smeared");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fSmeared = new TH2D(nameMomResSmeared->Data(),"", 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fSmeared);
		// Muon correction histos
		TString *nameMuonIdeal=new TString(namePC3->Data());
		nameMuonIdeal->Append("_MuonIdeal");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonIdeal = new TH3D(nameMuonIdeal->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonIdeal);
		TString *nameMuonSmeared=new TString(namePC3->Data());
		nameMuonSmeared->Append("_MuonSmeared");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonSmeared = new TH3D(nameMuonSmeared->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonSmeared);
		//
		TString *nameMuonPionK3=new TString(namePC3->Data());
		nameMuonPionK3->Append("_MuonPionK3");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonPionK3 = new TH3D(nameMuonPionK3->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fMuonPionK3);
		//
		TString *namePionPionK3=new TString(namePC3->Data());
		namePionPionK3->Append("_PionPionK3");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fPionPionK3 = new TH3D(namePionPionK3->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ3,0,fQupperBoundQ3);
		if(mb==0 && edB==0 && term<4) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fPionPionK3);
		
	      }// MCcase
	      //
	      if(c1==c2 && c1==c3 && term==4 ){
		TString *nameTwoPartNorm=new TString(namePC3->Data());
		nameTwoPartNorm->Append("_TwoPartNorm");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm = new TH2D(nameTwoPartNorm->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNorm);
		//
		TString *nameTwoPartNegNorm=new TString(namePC3->Data());
		nameTwoPartNegNorm->Append("_TwoPartNegNorm");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNegNorm = new TH2D(nameTwoPartNegNorm->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNegNorm);
		//
		TString *nameTwoPartNormErr=new TString(namePC3->Data());
		nameTwoPartNormErr->Append("_TwoPartNormErr");
		Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNormErr = new TH2D(nameTwoPartNormErr->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ3,0,fQupperBoundQ3);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].MB[mb].EDB[edB].ThreePT[term].fTwoPartNormErr);
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
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4 = new TH1D(name1DQ->Data(),"", fQbinsQ4,0,fQupperBoundQ4);
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTerms4);
		//
		TString *nameKfactor=new TString(namePC4->Data());
		nameKfactor->Append("_Kfactor");
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor = new TProfile(nameKfactor->Data(),"", fQbinsQ4,0,fQupperBoundQ4, 0,100, "");
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactor);
		//
		TString *nameKfactorW=new TString(namePC4->Data());
		nameKfactorW->Append("_KfactorWeighted");
		Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted = new TProfile(nameKfactorW->Data(),"", fQbinsQ4,0,fQupperBoundQ4, 0,100, "");
		fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fKfactorWeighted);
		//
		if(c1==c2 && c1==c3 && c1==c4 && term==12 ){
		  TString *nameTwoPartNorm=new TString(namePC4->Data());
		  nameTwoPartNorm->Append("_TwoPartNorm");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm = new TH2D(nameTwoPartNorm->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNorm);
		  //
		  TString *nameTwoPartNegNorm=new TString(namePC4->Data());
		  nameTwoPartNegNorm->Append("_TwoPartNegNorm");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNegNorm = new TH2D(nameTwoPartNegNorm->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNegNorm);
		  //
		  TString *nameTwoPartNormErr=new TString(namePC4->Data());
		  nameTwoPartNormErr->Append("_TwoPartNormErr");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNormErr = new TH2D(nameTwoPartNormErr->Data(),"", kDENtypes,0.5,kDENtypes+0.5, fQbinsQ4,0,fQupperBoundQ4);
		  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fTwoPartNormErr);
		}
		
		if(fMCcase==kTRUE){
		  // Momentum resolution correction histos
		  TString *nameMomResIdeal=new TString(namePC4->Data());
		  nameMomResIdeal->Append("_Ideal");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fIdeal = new TH2D(nameMomResIdeal->Data(),"", 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
		  if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fIdeal);
		  TString *nameMomResSmeared=new TString(namePC4->Data());
		  nameMomResSmeared->Append("_Smeared");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fSmeared = new TH2D(nameMomResSmeared->Data(),"", 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
		  if(mb==0) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fSmeared);
		  // Muon correction histos
		  TString *nameMuonIdeal=new TString(namePC4->Data());
		  nameMuonIdeal->Append("_MuonIdeal");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonIdeal = new TH3D(nameMuonIdeal->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonIdeal);
		  TString *nameMuonSmeared=new TString(namePC4->Data());
		  nameMuonSmeared->Append("_MuonSmeared");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonSmeared = new TH3D(nameMuonSmeared->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonSmeared);
		  //
		  TString *nameMuonPionK4=new TString(namePC4->Data());
		  nameMuonPionK4->Append("_MuonPionK4");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonPionK4 = new TH3D(nameMuonPionK4->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
		  if(mb==0 && edB==0 && term<12) fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fMuonPionK4);
		  //
		  TString *namePionPionK4=new TString(namePC4->Data());
		  namePionPionK4->Append("_PionPionK4");
		  Charge1[c1].Charge2[c2].Charge3[c3].Charge4[c4].MB[mb].EDB[edB].FourPT[term].fPionPionK4 = new TH3D(namePionPionK4->Data(),"", 2,0.5,2.5, 11,0.5,11.5, fQbinsQ4,0,fQupperBoundQ4);
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
	    
	    if(edB==0){
	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD = new TH3D(nameNum->Data(),"", kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	      fOutputList->Add(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fTerms2ThreeD);
	      
	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fTerms2ThreeD = new TH3D(nameDen->Data(),"", kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
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
  TH2D *fQ2Res = new TH2D("fQ2Res","",20,0,1, 200,-.2,.2);
  fOutputList->Add(fQ2Res);
  TH2D *fQ3Res = new TH2D("fQ3Res","",20,0,1, 200,-.3,.3);
  fOutputList->Add(fQ3Res);
  TH2D *fQ4Res = new TH2D("fQ4Res","",20,0,1, 200,-.4,.4);
  fOutputList->Add(fQ4Res);
  
  TH2D *DistQinv4pion = new TH2D("DistQinv4pion","",6,0.5,6.5, 20,0,0.1);
  fOutputList->Add(DistQinv4pion);
  TH2D *DistQinvMC4pion = new TH2D("DistQinvMC4pion","",6,0.5,6.5, 20,0,0.1);
  if(fMCcase) fOutputList->Add(DistQinvMC4pion);

  TH2D *fAvgQ12VersusQ3 = new TH2D("fAvgQ12VersusQ3","",10,0,0.1, 20,0,0.1);
  fOutputList->Add(fAvgQ12VersusQ3);
  TH2D *fAvgQ13VersusQ3 = new TH2D("fAvgQ13VersusQ3","",10,0,0.1, 20,0,0.1);
  fOutputList->Add(fAvgQ13VersusQ3);
  TH2D *fAvgQ23VersusQ3 = new TH2D("fAvgQ23VersusQ3","",10,0,0.1, 20,0,0.1);
  fOutputList->Add(fAvgQ23VersusQ3);

  TH1D *fDistPionParents4 = new TH1D("fDistPionParents4","",4,0.5,4.5);
  fOutputList->Add(fDistPionParents4);

  TH2D *fDistTPCNclsFindable = new TH2D("fDistTPCNclsFindable","", 100,0,0.5, 201,-0.5,200.5);
  fDistTPCNclsFindable->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsFindable->GetYaxis()->SetTitle("Ncls Findable");
  fOutputList->Add(fDistTPCNclsFindable);
  TProfile *fProfileTPCNclsFindable = new TProfile("fProfileTPCNclsFindable","",100,0,0.5, 0,200, "");
  fProfileTPCNclsFindable->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsFindable->GetYaxis()->SetTitle("<Ncls Findable>");
  fOutputList->Add(fProfileTPCNclsFindable);
  //
  TH2D *fDistTPCNclsCrossed = new TH2D("fDistTPCNclsCrossed","",100,0,0.5, 201,-0.5,200.5);
  fDistTPCNclsCrossed->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsCrossed->GetYaxis()->SetTitle("Ncls Crossed");
  fOutputList->Add(fDistTPCNclsCrossed);
  TProfile *fProfileTPCNclsCrossed = new TProfile("fProfileTPCNclsCrossed","",100,0,0.5, 0,200, "");
  fProfileTPCNclsCrossed->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsCrossed->GetYaxis()->SetTitle("<Ncls Crossed>");
  fOutputList->Add(fProfileTPCNclsCrossed);
  //
  TH2D *fDistTPCNclsFindableRatio = new TH2D("fDistTPCNclsFindableRatio","",100,0,0.5, 100,0,1);
  fDistTPCNclsFindableRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsFindableRatio->GetYaxis()->SetTitle("Ncls / Ncls Findable");
  fOutputList->Add(fDistTPCNclsFindableRatio);
  TProfile *fProfileTPCNclsFindableRatio = new TProfile("fProfileTPCNclsFindableRatio","",100,0,0.5, 0,1, "");
  fProfileTPCNclsFindableRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsFindableRatio->GetYaxis()->SetTitle("<Ncls / Ncls Findable>");
  fOutputList->Add(fProfileTPCNclsFindableRatio);
  //
  TH2D *fDistTPCNclsCrossedRatio = new TH2D("fDistTPCNclsCrossedRatio","",100,0,0.5, 100,0,1);
  fDistTPCNclsCrossedRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fDistTPCNclsCrossedRatio->GetYaxis()->SetTitle("Ncls / Ncls Crossed");
  fOutputList->Add(fDistTPCNclsCrossedRatio);
  TProfile *fProfileTPCNclsCrossedRatio = new TProfile("fProfileTPCNclsCrossedRatio","",100,0,0.5, 0,1, "");
  fProfileTPCNclsCrossedRatio->GetXaxis()->SetTitle("pT (GeV/c)"); fProfileTPCNclsCrossedRatio->GetYaxis()->SetTitle("<Ncls / Ncls Crossed>");
  fOutputList->Add(fProfileTPCNclsCrossedRatio);

  TH2D *fc4QSFitNum = new TH2D("fc4QSFitNum","",7,0.5,7.5, fQbinsQ4,0,fQupperBoundQ4);
  fOutputList->Add(fc4QSFitNum);
  TH2D *fc4QSFitDen = new TH2D("fc4QSFitDen","",7,0.5,7.5, fQbinsQ4,0,fQupperBoundQ4);
  fOutputList->Add(fc4QSFitDen);
  
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
  }else {return;}

  ///////////////////////////////////////////////////////////
  const AliAODVertex *primaryVertexAOD;
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
  Int_t zbin=0;
  Double_t zstep=2*10/Double_t(fZvertexBins), zstart=-10.;
  /////////////////////////////////////////////////

  
  Float_t centralityPercentile=0;
  Float_t cStep=5.0, cStart=0;
  
 
  if(fAODcase){// AOD case
    
    if(fPbPbcase){
      centrality = fAOD->GetCentrality();
      centralityPercentile = centrality->GetCentralityPercentile("V0M");
      if(centralityPercentile == 0) {cout<<"Centrality = 0, skipping event"<<endl; return;}
      if((centralityPercentile < 5*fCentBinLowLimit) || (centralityPercentile>= 5*(fCentBinHighLimit+1))) {/*cout<<"Centrality out of Range.  Skipping Event"<<endl;*/ return;}
      cout<<"Centrality % = "<<centralityPercentile<<endl;
    }
    
    
    ((TH1F*)fOutputList->FindObject("fMultDist0"))->Fill(fAOD->GetNumberOfTracks());

    // Pile-up rejection
    AliAnalysisUtils *AnaUtil=new AliAnalysisUtils();
    if(!fPbPbcase) AnaUtil->SetUseMVPlpSelection(kTRUE);// use Multi-Vertex tool for pp and pPb
    else AnaUtil->SetUseMVPlpSelection(kFALSE);
    Bool_t pileUpCase=AnaUtil->IsPileUpEvent(fAOD); 
    if(pileUpCase) return;
    
    ////////////////////////////////
    // Vertexing
    ((TH1F*)fOutputList->FindObject("fMultDist1"))->Fill(fAOD->GetNumberOfTracks());
    primaryVertexAOD = fAOD->GetPrimaryVertex();
    vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
    
    if(fabs(vertex[2]) > 10) {cout<<"Zvertex Out of Range. Skip Event"<<endl; return;} // Z-Vertex Cut 
    ((TH3F*)fOutputList->FindObject("fVertexDist"))->Fill(vertex[0], vertex[1], vertex[2]);
    
    if(!fMCcase && primaryVertexAOD->GetNContributors() < 1) {cout<<"Bad Vertex. Skip Event"<<endl; return;}
   
    ((TH1F*)fOutputList->FindObject("fMultDist2"))->Fill(fAOD->GetNumberOfTracks());
 
    fBfield = fAOD->GetMagneticField();
    
    for(Int_t i=0; i<fZvertexBins; i++){
      if( (vertex[2] >= zstart+i*zstep) && (vertex[2] < zstart+(i+1)*zstep) ){
	zbin=i;
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
      AliAODTrack* aodtrack = fAOD->GetTrack(randomIndex[i]);
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
	  AliAODTrack* aodtrack2 = fAOD->GetTrack(randomIndex[j]);
	  if(!aodtrack2) continue;
	  if(!aodtrack2->TestFilterBit(BIT(fFilterBit))) continue;
	  
	  if(-(aodtrack->GetID()+1)==aodtrack2->GetID()) {goodTrackOtherFB=kTRUE; break;}
	  
	}
	if(!goodTrackOtherFB) continue;
      }
      
      
      if(aodtrack->Pt() < fMinPt) continue;
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
      
    
      
      if(fTempStruct[myTracks].fMom > fMaxPt) continue;// upper P bound
      //if(fTempStruct[myTracks].fPt > 0.9999) continue;// upper P bound
      //if(fTempStruct[myTracks].fP[2] > 0.9999) continue;// upper P bound
     
     
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
      if(fMCcase && !fPbPbcase) DoPIDWorkAround=kFALSE;
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
	  AliAODTrack* aodTrack2 = fAOD->GetTrack(j);
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
	  
	  //if(aodTrack2->Pt()<0.2) cout<<aodTrack2->GetTPCNclsF()<<"  "<<aodTrack2->GetTPCNCrossedRows()<<"  "<<aodTrack2->GetTPCNcls()<<"  "<<aodTrack2->GetTPCFoundFraction()<<endl;


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
      //if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;// superfluous
      ////////////////////////
      //if(fTempStruct[myTracks].fProton && fTempStruct[myTracks].fMom < 0.25) continue;//extra cut for protons// superfluous


   
      if(fTempStruct[myTracks].fCharge==+1) {
	((TH2F*)fOutputList->FindObject("fDCAxyDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }else {
	((TH2F*)fOutputList->FindObject("fDCAxyDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }
     
      ((TH3F*)fOutputList->FindObject("fPhiPtDist"))->Fill(aodtrack->Charge(), aodtrack->Phi(), aodtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fPtEtaDist"))->Fill(aodtrack->Charge(), aodtrack->Pt(), aodtrack->Eta());
      
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
      if(mcParticle->Pt() < fMinPt || mcParticle->Pt() > fMaxPt) continue;
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
  if(myTracks < 4) {cout<<"Less than 4 tracks. Skipping Event."<<endl; return;}
  /////////////////////////////////////////
 

  ////////////////////////////////
  ///////////////////////////////
  // Mbin determination
  //
  // Mbin set to Pion Count Only for pp!!!!!!!
  fMbin=-1;
  if(!fPbPbcase){
    for(Int_t i=0; i<kMultBinspp; i++){
      if( ( pionCount > fMultLimits[i]) && ( pionCount <= fMultLimits[i+1]) ) { fMbin=i; break;}
      // Mbin 0 has 1 pion
    }
  }else{
    for(Int_t i=0; i<fCentBins; i++){
      if( (centralityPercentile >= cStart+i*cStep) && (centralityPercentile < cStart+(i+1)*cStep) ){
	fMbin=i;// 0 = most central
	break;
      }
    }
  }
  
  if(fMbin==-1) {cout<<"Bad Mbin+++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl; return;}
  
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
  fEDbin=0;// Extra Dimension bin (Kt, (Kt-Psi),....)
  //////////////////////////////////////////////////
  

  
  ((TH1F*)fOutputList->FindObject("fEvents1"))->Fill(fMbin+1);
  ((TProfile*)fOutputList->FindObject("fAvgMult"))->Fill(fMbin+1., pionCount);

  ////////////////////////////////////
  // Add event to buffer if > 0 tracks
  if(myTracks > 0){
    fEC[zbin][fMbin]->FIFOShift();
    (fEvt) = fEC[zbin][fMbin]->fEvtStr;
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
  Float_t kT12=0;
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
  Float_t qinv12MC=0, qinv13MC=0, qinv14MC=0, qinv23MC=0, qinv24MC=0, qinv34MC=0; 
  Float_t parentQinv12=0, parentQinv13=0, parentQinv14=0, parentQinv23=0, parentQinv24=0, parentQinv34=0;
  Float_t parentQ3=0;
  Float_t FSICorr12=0, FSICorr13=0, FSICorr14=0, FSICorr23=0, FSICorr24=0, FSICorr34=0;
  Bool_t pionParent1=kFALSE, pionParent2=kFALSE, pionParent3=kFALSE, pionParent4=kFALSE;
  Bool_t FilledMCpair12=kFALSE, FilledMCtriplet123=kFALSE;
  Bool_t Positive1stTripletWeights=kTRUE, Positive2ndTripletWeights=kTRUE;
  Float_t T12=0, T13=0, T14=0, T23=0, T24=0, T34=0;
  Int_t momBin12=1, momBin13=1, momBin14=1, momBin23=1, momBin24=1, momBin34=1;
  Float_t MomResCorr12=1.0, MomResCorr13=1.0, MomResCorr14=1.0, MomResCorr23=1.0, MomResCorr24=1.0, MomResCorr34=1.0;
  //
  AliAODMCParticle *mcParticle1=0x0;
  AliAODMCParticle *mcParticle2=0x0;
  

  ////////////////////
  //Int_t PairCount[7]={0};
  //Int_t NormPairCount[7]={0};
  Int_t KT3index=0, KT4index=0;

  // reset to defaults
  for(Int_t i=0; i<kMultLimitPbPb; i++) {
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
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[0].fTerms2->Fill(kT12, qinv12);
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[0].fTerms2QW->Fill(kT12, qinv12, qinv12);
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
	    if( (fEvt+en1)->fNtracks%100==0){
	      Int_t kTindex=0;
	      if(kT12>0.3) kTindex=1;
	      Int_t UnitMultBin = int((fEvt+en1)->fNtracks / 100.) + 1;
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[kTindex].TwoPT[0].fUnitMultBin->Fill(UnitMultBin, qinv12);
	    }
	    
	  }
	  if( (en1+en2==1)) {
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fTerms2->Fill(kT12, qinv12);
	    Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[0].TwoPT[1].fTerms2QW->Fill(kT12, qinv12, qinv12);
	    // osl frame
	    if((kT12 > 0.2) && (kT12 < 0.3)){  
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    if((kT12 > 0.6) && (kT12 < 0.7)){  
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fTerms2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fTerms2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	    }
	    // unit mult bins
	    if( (fEvt+en1)->fNtracks%100==0){
	      Int_t kTindex=0;
	      if(kT12>0.3) kTindex=1;
	      Int_t UnitMultBin = int((fEvt+en1)->fNtracks / 100.) + 1;
	      Charge1[bin1].Charge2[bin2].MB[fMbin].EDB[kTindex].TwoPT[1].fUnitMultBin->Fill(UnitMultBin, qinv12);
	    }
	  }
	  //////////////////////////////////////////
	  if(fTabulatePairs && en1==0 && en2<=1 && bin1==bin2){
	    Float_t kY = 0;
	    Int_t kTbin=-1, kYbin=-1;
	    
	    for(Int_t kIt=0; kIt<fKbinsT; kIt++) {if(kT12 < (fKmiddleT[kIt] + fKstepT[kIt]/2.)) {kTbin = kIt; break;}} 
	    for(Int_t kIt=0; kIt<fKbinsY; kIt++) {if(kY < (fKmiddleY[kIt] + fKstepY[kIt]/2.)) {kYbin = kIt; break;}}
	    if((kTbin<0) || (kYbin<0)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    if((kTbin>=fKbinsT) || (kYbin>=fKbinsY)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    if(fGenerateSignal && en2==0) {
	      Int_t chGroup2[2]={ch1,ch2};
	      Float_t WInput = MCWeight(chGroup2, fRMax, ffcSqMRC, qinv12, kT12);
	      KT[kTbin].KY[kYbin].MB[fMbin].EDB[0].TwoPT[en2].fTerms2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong), WInput);
	    }else KT[kTbin].KY[kYbin].MB[fMbin].EDB[0].TwoPT[en2].fTerms2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong));
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
    
  //cout<<PairCount[0]<<"  "<<PairCount[1]<<"  "<<PairCount[2]<<"  "<<PairCount[3]<<"  "<<PairCount[4]<<"  "<<PairCount[5]<<"  "<<PairCount[6]<<endl;
  //cout<<NormPairCount[0]<<"  "<<NormPairCount[1]<<"  "<<NormPairCount[2]<<"  "<<NormPairCount[3]<<"  "<<NormPairCount[4]<<"  "<<NormPairCount[5]<<"  "<<NormPairCount[6]<<endl;
  ///////////////////////////////////////////////////  
  // Do not use pairs from events with too many pairs
  
  ((TH1F*)fOutputList->FindObject("fEvents2"))->Fill(fMbin+1);
  
  ///////////////////////////////////////////////////
  
  
  if(fTabulatePairs) return;

  /*TF1 *SCpairWeight = new TF1("SCpairWeight","[0] + [1]*x + [2]*exp(-[3]*x)",0,0.2);// same-charge pair weight for monte-carlo data without two-track cuts.
  SCpairWeight->FixParameter(0, 0.959);
  SCpairWeight->FixParameter(1, 0.278);
  SCpairWeight->FixParameter(2, -1.759);
  SCpairWeight->FixParameter(3, 115.107);*/

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
	      if(KT3<=fKT3transition) KT3index=0;
	      else KT3index=1;
	      
	      if(en2==0 && en3==0 && en4==0) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fNorm3->Fill(0);
	      if(en2==1 && en3==2 && en4==3) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fNorm3->Fill(0);
	      if(en2==0 && en3==1 && en4==2) {
		if(fill2) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fNorm3->Fill(0);
		if(fill3) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fNorm3->Fill(0);
		if(fill4) Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fNorm3->Fill(0);
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
		if(KT4<=fKT4transition) KT4index=0;
		else KT4index=1;
		
		Bool_t FillTerms[13]={kFALSE};
		SetFillBins4(ch1, ch2, ch3, ch4, bin1, bin2, bin3, bin4, en2+en3+en4, FillTerms);
		//
		for(int ft=0; ft<13; ft++) {
		  if(FillTerms[ft]) Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[ft].fNorm4->Fill(0.); 
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

	    /////////////////////////////////////////////////////////////
	    for (Int_t j=i+1; j<(fEvt+en2)->fNtracks; j++) {// 2nd particle
	      if(en2==0) {if(fLowQPairSwitch_E0E0[i]->At(j)=='0') continue;}
	      else {if(fLowQPairSwitch_E0E1[i]->At(j)=='0') continue;}

	      pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	      ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	      qinv12 = GetQinv(pVect1, pVect2);
	      kT12 = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;
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
		    Int_t chGroup2[2]={ch1,ch2};

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
			      Float_t Rvalue = 5+Riter;
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
		  
		    
		    // momentum resolution
		    for(Int_t Riter=0; Riter<fRVALUES; Riter++){
		      Float_t Rvalue = 5+Riter;
		      Float_t WInput = MCWeight(chGroup2, Rvalue, ffcSqMRC, qinv12MC, 0.);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[0].fIdeal->Fill(Rvalue, qinv12MC, WInput);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[1].fIdeal->Fill(Rvalue, qinv12MC);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[0].fSmeared->Fill(Rvalue, qinv12, WInput);
		      Charge1[bin1].Charge2[bin2].MB[0].EDB[kTindex].TwoPT[1].fSmeared->Fill(Rvalue, qinv12);
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
		
		pVect3[0]=(fEvt+en3)->fTracks[k].fEaccepted;
		pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
		pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
		pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
		ch3 = Int_t(((fEvt+en3)->fTracks[k].fCharge + 1)/2.);
		qinv13 = GetQinv(pVect1, pVect3);
		qinv23 = GetQinv(pVect2, pVect3);
		q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));
		
		FilledMCtriplet123 = kFALSE;
		
		Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
		SetFillBins3(ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);
		
		Float_t KT3 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
		if(KT3<=fKT3transition) KT3index=0;
		else KT3index=1;
		
		FSICorr13 = FSICorrelation(ch1,ch3, qinv13);
		FSICorr23 = FSICorrelation(ch2,ch3, qinv23);
		if(!fGenerateSignal && !fMCcase) {
		  momBin12 = fMomResC2SC->GetYaxis()->FindBin(qinv12);
		  momBin13 = fMomResC2SC->GetYaxis()->FindBin(qinv13);
		  momBin23 = fMomResC2SC->GetYaxis()->FindBin(qinv23);		  
		  if(momBin12 >= 20) momBin12 = 19;
		  if(momBin13 >= 20) momBin13 = 19;
		  if(momBin23 >= 20) momBin23 = 19;
		  //
		  if(ch1==ch2) MomResCorr12 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin12);
		  else MomResCorr12 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin12);
		  if(ch1==ch3) MomResCorr13 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin13);
		  else MomResCorr13 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin13);
		  if(ch2==ch3) MomResCorr23 = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin23);
		  else MomResCorr23 = fMomResC2MC->GetBinContent(rBinForTPNMomRes, momBin23);
		}
		
		if(ENsum==0) {
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fTerms3->Fill(q3); 
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fKfactor->Fill(q3, 1/(FSICorr12*FSICorr13*FSICorr23));
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fKfactorWeighted->Fill(q3, 1/(FSICorr12*FSICorr13*FSICorr23), MomResCorr12*MomResCorr13*MomResCorr23);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fMeanQinv->Fill(q3, qinv12);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fMeanQinv->Fill(q3, qinv13);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[0].fMeanQinv->Fill(q3, qinv23);
		}
		if(ENsum==6) {
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTerms3->Fill(q3);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fMeanQinv->Fill(q3, qinv12);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fMeanQinv->Fill(q3, qinv13);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fMeanQinv->Fill(q3, qinv23);
		}
		if(ENsum==3){
		  if(fill2) {
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fTerms3->Fill(q3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fKfactor->Fill(q3, 1/(FSICorr12));
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fMeanQinv->Fill(q3, qinv12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fMeanQinv->Fill(q3, qinv13);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[1].fMeanQinv->Fill(q3, qinv23);
		  }if(fill3) {
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fTerms3->Fill(q3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fKfactor->Fill(q3, 1/(FSICorr12));
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fMeanQinv->Fill(q3, qinv12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fMeanQinv->Fill(q3, qinv13);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[2].fMeanQinv->Fill(q3, qinv23);
		  }if(fill4) {
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fTerms3->Fill(q3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fKfactor->Fill(q3, 1/(FSICorr12));
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fKfactorWeighted->Fill(q3, 1/(FSICorr12), MomResCorr12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fMeanQinv->Fill(q3, qinv12);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fMeanQinv->Fill(q3, qinv13);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[3].fMeanQinv->Fill(q3, qinv23);
		  }
		}
		
		// r3 denominator
		if(ENsum==6 && ch1==ch2 && ch1==ch3){
		  Positive1stTripletWeights = kTRUE;
		  //
		  GetWeight(pVect1, pVect2, weight12, weight12Err);
		  GetWeight(pVect1, pVect3, weight13, weight13Err);
		  GetWeight(pVect2, pVect3, weight23, weight23Err);
		  
		  if(sqrt(fabs(weight12*weight13*weight23)) > 1.0) {// weight should never be larger than 1
		    if(fMbin==0 && bin1==0) {
		      ((TH1D*)fOutputList->FindObject("fTPNRejects3pion1"))->Fill(q3, sqrt(fabs(weight12*weight13*weight23)));
		    }
		  }else{
		    
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNorm->Fill(1, q3, sqrt(weight12CC[0]*weight13CC[0]*weight23CC[0]));
		    }
		    // no Muon Correction
		    weight12CC[1] = ((weight12+1)*MomResCorr12 - ffcSq*FSICorr12 - (1-ffcSq));
		    weight12CC[1] /= FSICorr12*ffcSq;
		    weight13CC[1] = ((weight13+1)*MomResCorr13 - ffcSq*FSICorr13 - (1-ffcSq));
		    weight13CC[1] /= FSICorr13*ffcSq;
		    weight23CC[1] = ((weight23+1)*MomResCorr23 - ffcSq*FSICorr23 - (1-ffcSq));
		    weight23CC[1] /= FSICorr23*ffcSq;
		    if(weight12CC[1] > 0 && weight13CC[1] > 0 && weight23CC[1] > 0){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNorm->Fill(2, q3, sqrt(weight12CC[1]*weight13CC[1]*weight23CC[1]));
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
		    
		    if(weight12CC[2] < 0 || weight13CC[2] < 0 || weight23CC[2] < 0) {// C2^QS can never be less than unity
		      if(fMbin==0 && bin1==0) {
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNorm->Fill(3, q3, weightTotal);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNorm->Fill(4, q3, 1);
		    }else{
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNegNorm->Fill(4, q3, 1);
		    }
		    //
		    // Full Weight reconstruction
		    
		    for(Int_t RcohIndex=0; RcohIndex<2; RcohIndex++){// Rcoh=0, then Rcoh=Rch
		      for(Int_t GIndex=0; GIndex<50; GIndex++){
			Int_t FillBin = 5 + RcohIndex*50 + GIndex;
			Float_t G = 0.02*GIndex;
			if(RcohIndex==0){
			  T12 = (-2*G*(1-G) + sqrt(pow(2*G*(1-G),2) + 4*pow(1-G,2)*weight12CC[2])) / (2*pow(1-G,2));
			  T13 = (-2*G*(1-G) + sqrt(pow(2*G*(1-G),2) + 4*pow(1-G,2)*weight13CC[2])) / (2*pow(1-G,2));
			  T23 = (-2*G*(1-G) + sqrt(pow(2*G*(1-G),2) + 4*pow(1-G,2)*weight23CC[2])) / (2*pow(1-G,2));
			  weightTotal = 2*G*(1-G)*(T12 + T13 + T23) + pow(1-G,2)*(T12*T12 + T13*T13 + T23*T23);
			  weightTotal += 2*G*pow(1-G,2)*(T12*T13 + T12*T23 + T13*T23) + 2*pow(1-G,3)*T12*T13*T23;
			}else{
			  T12 = sqrt(weight12CC[2] / (1-G*G));
			  T13 = sqrt(weight13CC[2] / (1-G*G));
			  T23 = sqrt(weight23CC[2] / (1-G*G));
			  weightTotal = (1-G*G)*(T12*T12 + T13*T13 + T23*T23);
			  weightTotal += (6*G*pow(1-G,2) + 2*pow(1-G,3)) * T12*T13*T23;
			}
			if(Positive1stTripletWeights){
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNorm->Fill(FillBin, q3, weightTotal);
			}else{
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNegNorm->Fill(FillBin, q3, weightTotal);
			}
		      }
		    }
		    //
		    /*weight12CC_e = weight12Err*MomResCorr12 / FSICorr12 / ffcSq * MuonCorr12;
		      weight13CC_e = weight13Err*MomResCorr13 / FSICorr13 / ffcSq * MuonCorr13;
		      weight23CC_e = weight23Err*MomResCorr23 / FSICorr23 / ffcSq * MuonCorr23;
		      if(weight12CC[2]*weight13CC[2]*weight23CC[2] > 0){
		      weightTotalErr = pow(2 * sqrt(3) * weight12CC_e*weight13CC[2]*weight23CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight23CC[2]),2);
		      }
		      weightTotalErr += pow(weight12CC_e,2) + pow(weight13CC_e,2) + pow(weight23CC_e,2);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[fMbin].EDB[KT3index].ThreePT[4].fTwoPartNormErr->Fill(4, q3, weightTotalErr);*/
		    
		  }// 1st r3 den check
		  
		}// r3 den
		

		if(ch1==ch2 && ch1==ch3 && ENsum==0){
		  ((TH3D*)fOutputList->FindObject("fKT3DistTerm1"))->Fill(fMbin+1, KT3, q3);
		  if(q3<0.1){
		    Float_t pt1=sqrt(pow(pVect1[1],2)+pow(pVect1[2],2));
		    Float_t pt2=sqrt(pow(pVect2[1],2)+pow(pVect2[2],2));
		    Float_t pt3=sqrt(pow(pVect3[1],2)+pow(pVect3[2],2));
		    ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, KT3index, pt1);
		    ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, KT3index, pt2);
		    ((TProfile2D*)fOutputList->FindObject("fKT3AvgpT"))->Fill(fMbin+1, KT3index, pt3);
		    if(fMbin==0){
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpT"))->Fill(KT3index, q3, pt1);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpT"))->Fill(KT3index, q3, pt2);
		      ((TH3D*)fOutputList->FindObject("fQ3AvgpT"))->Fill(KT3index, q3, pt3);
		    }
		  }
		  
		}
		if(ch1==ch2 && ch1==ch3 && ENsum==6) ((TH3D*)fOutputList->FindObject("fKT3DistTerm5"))->Fill(fMbin+1, KT3, q3);
	      	
		

		
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
		    
		    Float_t TripletWeightTTC=1.0;// same-charge weights to mimic two-track depletion of same-charge pairs
		    //if(ch1==ch2 && qinv12>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv12);
		    //if(ch1==ch3 && qinv13>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv13);
		    //if(ch2==ch3 && qinv23>0.006) TripletWeightTTC *= SCpairWeight->Eval(qinv23);
		    
		    Int_t chGroup3[3]={ch1,ch2,ch3};
		    Float_t QinvMCGroup3[3]={qinv12MC,qinv13MC,qinv23MC};
		    //Float_t kTGroup3[3]={float(sqrt(pow(pVect1MC[1]+pVect2MC[1],2) + pow(pVect1MC[2]+pVect2MC[2],2))/2.),
		    //float(sqrt(pow(pVect1MC[1]+pVect3MC[1],2) + pow(pVect1MC[2]+pVect3MC[2],2))/2.),
		    //float(sqrt(pow(pVect2MC[1]+pVect3MC[1],2) + pow(pVect2MC[2]+pVect3MC[2],2))/2.)};
		    Float_t kTGroup3[3]={0};
		    
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

			for(Int_t term=1; term<=4; term++){
			  if(term==1) {}
			  else if(term==2) {if(!pionParent1 && !pionParent2) continue;}
			  else if(term==3) {if(!pionParent1 && !pionParent3) continue;}
			  else {if(!pionParent2 && !pionParent3) continue;}
			  for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			    Float_t Rvalue = 5+Riter;
			    Float_t WInput = MCWeight3(term, Rvalue, 1.0, chGroup3, parentQinvGroup3, parentkTGroup3);
			    Float_t WInputParentFSI = MCWeightFSI3(term, Rvalue, 1.0, chGroup3, parentQinvGroup3);
			    Float_t WInputFSI = MCWeightFSI3(term, Rvalue, 1.0, chGroup3, QinvMCGroup3);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonSmeared->Fill(1, Rvalue, q3MC, WInput*TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonIdeal->Fill(1, Rvalue, parentQ3, WInput*TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonPionK3->Fill(1, Rvalue, q3MC, WInputFSI*TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fPionPionK3->Fill(1, Rvalue, parentQ3, WInputParentFSI*TripletWeightTTC);
			    //
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonSmeared->Fill(2, Rvalue, q3MC, TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonIdeal->Fill(2, Rvalue, parentQ3, TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fMuonPionK3->Fill(2, Rvalue, q3MC, TripletWeightTTC);
			    Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[0].ThreePT[term-1].fPionPionK3->Fill(2, Rvalue, parentQ3, TripletWeightTTC);
			  }// Riter
			}// term loop
		    
		      }// pion parent check
		    }// parentQ check (muon correction)
		    
		    // 3-pion momentum resolution
		    for(Int_t term=1; term<=5; term++){
		      for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			Float_t Rvalue = 5+Riter;
			Float_t WInput = MCWeight3(term, Rvalue, ffcSqMRC, chGroup3, QinvMCGroup3, kTGroup3);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[KT3index].ThreePT[term-1].fIdeal->Fill(Rvalue, q3MC, WInput*TripletWeightTTC);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].MB[0].EDB[KT3index].ThreePT[term-1].fSmeared->Fill(Rvalue, q3, WInput*TripletWeightTTC);
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

		  pVect4[0]=(fEvt+en4)->fTracks[l].fEaccepted;
		  pVect4[1]=(fEvt+en4)->fTracks[l].fP[0];
		  pVect4[2]=(fEvt+en4)->fTracks[l].fP[1];
		  pVect4[3]=(fEvt+en4)->fTracks[l].fP[2];
		  ch4 = Int_t(((fEvt+en4)->fTracks[l].fCharge + 1)/2.);
		  qinv14 = GetQinv(pVect1, pVect4);
		  qinv24 = GetQinv(pVect2, pVect4);
		  qinv34 = GetQinv(pVect3, pVect4);
		  q4 = sqrt(pow(q3,2) + pow(qinv14,2) + pow(qinv24,2) + pow(qinv34,2));
		  
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6){
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(1, qinv12); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(2, qinv13); 
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(3, qinv14); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(4, qinv23); 
		    ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(5, qinv24); ((TH2D*)fOutputList->FindObject("DistQinv4pion"))->Fill(6, qinv34);
		  }
		  
		  Float_t KT4 = sqrt(pow(pVect1[1]+pVect2[1]+pVect3[1]+pVect4[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2]+pVect4[2],2))/4.;
		  if(KT4<=fKT4transition) KT4index=0;
		  else KT4index=1;
		  
		  FSICorr14 = FSICorrelation(ch1,ch4, qinv14);
		  FSICorr24 = FSICorrelation(ch2,ch4, qinv24);
		  FSICorr34 = FSICorrelation(ch3,ch4, qinv34);
		  
		  if(!fGenerateSignal && !fMCcase) {
		    momBin14 = fMomResC2SC->GetYaxis()->FindBin(qinv14);
		    momBin24 = fMomResC2SC->GetYaxis()->FindBin(qinv24);
		    momBin34 = fMomResC2SC->GetYaxis()->FindBin(qinv34);		  
		    if(momBin14 >= 20) momBin14 = 19;
		    if(momBin24 >= 20) momBin24 = 19;
		    if(momBin34 >= 20) momBin34 = 19;
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[ft].fTerms4->Fill(q4);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[ft].fKfactor->Fill(q4, FSIfactor);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[ft].fKfactorWeighted->Fill(q4, FSIfactor, MomResWeight);
		    }
		  }
		  
		  /////////////////////////////////////////////////////////////
		  // r4{2}
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6){
		    Positive2ndTripletWeights=kTRUE;
		    //
		    GetWeight(pVect1, pVect4, weight14, weight14Err);
		    GetWeight(pVect2, pVect4, weight24, weight24Err);
		    GetWeight(pVect3, pVect4, weight34, weight34Err);
		    
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNorm->Fill(1, q4, weightTotal);
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNorm->Fill(2, q4, weightTotal);
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
 
		    if(weight14CC[2] < 0 || weight24CC[2] < 0 || weight34CC[2] < 0) {// C2^QS can never be less than unity
		      if(fMbin==0 && bin1==0 && KT4index==0) {
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
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNorm->Fill(3, q4, weightTotal);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNorm->Fill(4, q4, 1);
		    }else{
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNegNorm->Fill(4, q4, 1);
		    }
		    // Full Weight reconstruction
		    for(Int_t RcohIndex=0; RcohIndex<2; RcohIndex++){// Rcoh=0, then Rcoh=Rch
		      for(Int_t GIndex=0; GIndex<50; GIndex++){
			Int_t FillBin = 5 + RcohIndex*50 + GIndex;
			Float_t G = 0.02*GIndex;
			if(RcohIndex==0){// Rcoh=0
			  Float_t a = pow(1-G,2);
			  Float_t b = 2*G*(1-G);
			  T12 = (-b + sqrt(pow(b,2) + 4*a*weight12CC[2])) / (2*a);
			  T13 = (-b + sqrt(pow(b,2) + 4*a*weight13CC[2])) / (2*a);
			  T14 = (-b + sqrt(pow(b,2) + 4*a*weight14CC[2])) / (2*a);
			  T23 = (-b + sqrt(pow(b,2) + 4*a*weight23CC[2])) / (2*a);
			  T24 = (-b + sqrt(pow(b,2) + 4*a*weight24CC[2])) / (2*a);
			  T34 = (-b + sqrt(pow(b,2) + 4*a*weight34CC[2])) / (2*a);
			  weightTotal = 2*G*(1-G)*(T12 + T13 + T14 + T23 + T24 + T34) + pow(1-G,2)*(T12*T12 + T13*T13 + T14*T14 + T23*T23 + T24*T24 + T34*T34);// 2-pion
			  weightTotal += 2*G*pow(1-G,3)*(T12*T34*T34 + T12*T12*T34 + T13*T24*T24 + T13*T13*T24 + T14*T23*T23 + T14*T14*T23);// 2-pair
			  weightTotal += pow(1-G,4)*(pow(T12,2)*pow(T34,2) + pow(T13,2)*pow(T24,2) + pow(T14,2)*pow(T23,2));// 2-pair fully chaotic
			  weightTotal += 2*G*pow(1-G,2)*(T12*T13 + T12*T23 + T13*T23  + T12*T14 + T12*T24 + T14*T24);// 3-pion
			  weightTotal += 2*G*pow(1-G,2)*(T13*T14 + T13*T34 + T14*T34  + T23*T24 + T23*T34 + T24*T34);// 3-pion
			  weightTotal += 2*pow(1-G,3)*(T12*T13*T23 + T12*T14*T24 + T13*T14*T34 + T23*T24*T34);// 3-pion fully chaotic
			  weightTotal += 2*G*pow(1-G,3)*(T12*T14*T34 + T12*T14*T23 + T12*T23*T34 + T14*T23*T34);// 4-pion
			  weightTotal += 2*G*pow(1-G,3)*(T12*T13*T34 + T12*T34*T24 + T12*T24*T13 + T13*T24*T34);// 4-pion
			  weightTotal += 2*G*pow(1-G,3)*(T14*T13*T23 + T14*T13*T24 + T13*T23*T24 + T14*T24*T23);// 4-pion
			  weightTotal += 2*pow(1-G,4)*(T12*T13*T24*T34 + T12*T14*T23*T34 + T13*T14*T23*T24);// 4-pion fully chaotic
			}else{// Rcoh=Rch
			  T12 = sqrt(weight12CC[2] / (1-G*G));
			  T13 = sqrt(weight13CC[2] / (1-G*G));
			  T14 = sqrt(weight14CC[2] / (1-G*G));
			  T23 = sqrt(weight23CC[2] / (1-G*G));
			  T24 = sqrt(weight24CC[2] / (1-G*G));
			  T34 = sqrt(weight34CC[2] / (1-G*G));
			  weightTotal = (1-G*G)*(T12*T12 + T13*T13 + T14*T14 + T23*T23 + T24*T24 + T34*T34);// 2-pion
			  weightTotal += (4*G*pow(1-G,3)+pow(1-G,4))*(pow(T12,2)*pow(T34,2) + pow(T13,2)*pow(T24,2) + pow(T14,2)*pow(T23,2));// 2-pair
			  weightTotal += (6*G*pow(1-G,2) + 2*pow(1-G,3))*(T12*T13*T23);// 3-pion
			  weightTotal += (6*G*pow(1-G,2) + 2*pow(1-G,3))*(T12*T14*T24);// 3-pion
			  weightTotal += (6*G*pow(1-G,2) + 2*pow(1-G,3))*(T13*T14*T34);// 3-pion
			  weightTotal += (6*G*pow(1-G,2) + 2*pow(1-G,3))*(T23*T24*T34);// 3-pion
			  weightTotal += (8*G*pow(1-G,3) + 2*pow(1-G,4))*(T12*T13*T24*T34);// 4-pion
			  weightTotal += (8*G*pow(1-G,3) + 2*pow(1-G,4))*(T12*T14*T23*T34);// 4-pion
			  weightTotal += (8*G*pow(1-G,3) + 2*pow(1-G,4))*(T13*T14*T23*T24);// 4-pion
			}
			if(Positive1stTripletWeights && Positive2ndTripletWeights){
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNorm->Fill(FillBin, q4, weightTotal);
			}else{
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNegNorm->Fill(FillBin, q4, weightTotal);
			}
		      }
		    }
		    // stat errors
		    /*weight14CC_e = weight14Err*MomResCorr14 / FSICorr14 / ffcSq * MuonCorr14;
		      weight24CC_e = weight24Err*MomResCorr24 / FSICorr24 / ffcSq * MuonCorr24;
		      weight34CC_e = weight34Err*MomResCorr34 / FSICorr34 / ffcSq * MuonCorr34;
		      if(weight12CC[2]*weight13CC[2]*weight24CC[2]*weight34CC[2] > 0){
		      weightTotalErr = pow( 6 * 2 * weight12CC_e*weight13CC[2]*weight24CC[2]*weight34CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight24CC[2]*weight34CC[2]),2);
		      }
		      if(weight12CC[2]*weight13CC[2]*weight23CC[2] > 0){
		      weightTotalErr += pow( 8 * sqrt(3) * weight12CC_e*weight13CC[2]*weight23CC[2] / sqrt(weight12CC[2]*weight13CC[2]*weight23CC[2]),2);
		      }
		      weightTotalErr += 2*(pow(weight12CC_e*weight34CC[2],2) + pow(weight13CC_e*weight24CC[2],2) + pow(weight14CC_e*weight23CC[2],2));
		      weightTotalErr += pow(weight12CC_e,2) + pow(weight13CC_e,2) + pow(weight14CC_e,2) + pow(weight23CC_e,2) + pow(weight24CC_e,2) + pow(weight34CC_e,2);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[fMbin].EDB[KT4index].FourPT[12].fTwoPartNormErr->Fill(4, q4, weightTotalErr);
		    */
		    if(fMbin==0 && KT4index==0){
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
		  
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==0){
		    ((TH3D*)fOutputList->FindObject("fKT4DistTerm1"))->Fill(fMbin+1, KT4, q4);
		    if(q4<0.105){
		      Float_t pt1=sqrt(pow(pVect1[1],2)+pow(pVect1[2],2));
		      Float_t pt2=sqrt(pow(pVect2[1],2)+pow(pVect2[2],2));
		      Float_t pt3=sqrt(pow(pVect3[1],2)+pow(pVect3[2],2));
		      Float_t pt4=sqrt(pow(pVect4[1],2)+pow(pVect4[2],2));
		      ((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, KT4index, pt1);
		      ((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, KT4index, pt2);
		      ((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, KT4index, pt3);
		      ((TProfile2D*)fOutputList->FindObject("fKT4AvgpT"))->Fill(fMbin+1, KT4index, pt4);
		      if(fMbin==0){
			((TH3D*)fOutputList->FindObject("fQ4AvgpT"))->Fill(KT4index, q4, pt1);
			((TH3D*)fOutputList->FindObject("fQ4AvgpT"))->Fill(KT4index, q4, pt2);
			((TH3D*)fOutputList->FindObject("fQ4AvgpT"))->Fill(KT4index, q4, pt3);
			((TH3D*)fOutputList->FindObject("fQ4AvgpT"))->Fill(KT4index, q4, pt4);
		      }
		    }
		  }
		  if(ch1==ch2 && ch1==ch3 && ch1==ch4 && ENsum==6) ((TH3D*)fOutputList->FindObject("fKT4DistTerm13"))->Fill(fMbin+1, KT4, q4);


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

		      Float_t QuadWeightTTC=1.0;// same-charge weights to mimic two-track depletion of same-charge pairs
		      //if(ch1==ch2 && qinv12>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv12);
		      //if(ch1==ch3 && qinv13>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv13);
		      //if(ch1==ch4 && qinv14>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv14);
		      //if(ch2==ch3 && qinv23>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv23);
		      //if(ch2==ch4 && qinv24>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv24);
		      //if(ch3==ch4 && qinv34>0.006) QuadWeightTTC *= SCpairWeight->Eval(qinv34);
		      

		      Int_t chGroup4[4]={ch1,ch2,ch3,ch4};
		      Float_t QinvMCGroup4[6]={qinv12MC, qinv13MC, qinv14MC, qinv23MC, qinv24MC, qinv34MC};
		      Float_t kTGroup4[6]={0};
		      
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
			      Float_t Rvalue = 5+Riter;
			      Float_t WInput = MCWeight4(term, Rvalue, 1.0, chGroup4, parentQinvGroup4, parentkTGroup4);
			      Float_t WInputParentFSI = MCWeightFSI4(term, Rvalue, 1.0, chGroup4, parentQinvGroup4);
			      Float_t WInputFSI = MCWeightFSI4(term, Rvalue, 1.0, chGroup4, QinvMCGroup4);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonSmeared->Fill(1, Rvalue, q4MC, WInput*QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonIdeal->Fill(1, Rvalue, parentQ4, WInput*QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonPionK4->Fill(1, Rvalue, q4MC, WInputFSI*QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fPionPionK4->Fill(1, Rvalue, parentQ4, WInputParentFSI*QuadWeightTTC);
			      //
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonSmeared->Fill(2, Rvalue, q4MC, QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonIdeal->Fill(2, Rvalue, parentQ4, QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fMuonPionK4->Fill(2, Rvalue, q4MC, QuadWeightTTC);
			      Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[0].FourPT[term-1].fPionPionK4->Fill(2, Rvalue, parentQ4, QuadWeightTTC);
			    }// Riter
			  }// term loop
			  
			}// pion parent check
		      }// parentQ check (muon correction)
		      
		      // 4-pion momentum resolution
		      for(Int_t term=1; term<=13; term++){
			for(Int_t Riter=0; Riter<fRVALUES; Riter++){
			  Float_t Rvalue = 5+Riter;
			  Float_t WInput = MCWeight4(term, Rvalue, ffcSqMRC, chGroup4, QinvMCGroup4, kTGroup4);
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[KT4index].FourPT[term-1].fIdeal->Fill(Rvalue, q4MC, WInput*QuadWeightTTC);
			  Charge1[bin1].Charge2[bin2].Charge3[bin3].Charge4[bin4].MB[0].EDB[KT4index].FourPT[term-1].fSmeared->Fill(Rvalue, q4, WInput*QuadWeightTTC);
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
  
  /* Int_t ncl1 = first.fClusterMap.GetNbits();
  Int_t ncl2 = second.fClusterMap.GetNbits();
  Int_t sumCls = 0; Int_t sumSha = 0; Int_t sumQ = 0;
  Double_t shfrac = 0; Double_t qfactor = 0;
  for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
    if (first.fClusterMap.TestBitNumber(imap) && second.fClusterMap.TestBitNumber(imap)) {// Both clusters
      if (first.fSharedMap.TestBitNumber(imap) && second.fSharedMap.TestBitNumber(imap)) { // Shared
	sumQ++;
	sumCls+=2;
	sumSha+=2;}
      else {sumQ--; sumCls+=2;}
    }
    else if (first.fClusterMap.TestBitNumber(imap) || second.fClusterMap.TestBitNumber(imap)) {// Non shared
      sumQ++;
      sumCls++;}
  }
  if (sumCls>0) {
    qfactor = sumQ*1.0/sumCls;
    shfrac = sumSha*1.0/sumCls;
  }
  
  if(qfactor > fShareQuality || shfrac > fShareFraction) return kFALSE;
  */
  
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
void AliFourPion::SetWeightArrays(Bool_t legoCase, TH3F *histos[AliFourPion::fKbinsT][AliFourPion::fCentBins]){

  if(legoCase){
    cout<<"LEGO call to SetWeightArrays"<<endl;
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t mb=0; mb<fCentBins; mb++){
	fNormWeight[tKbin][mb] = (TH3F*)histos[tKbin][mb]->Clone();
	fNormWeight[tKbin][mb]->SetDirectory(0);
      }
    }
    
  }else{
    
    TFile *wFile = new TFile("WeightFile.root","READ");
    if(!wFile->IsOpen()) {cout<<"No Weight File!!!!!!!!!!"<<endl; return;}
    else cout<<"Good Weight File Found!"<<endl;
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t mb=0; mb<fCentBins; mb++){
		    
	TString *name = new TString("Weight_Kt_");
	*name += tKbin;
	name->Append("_Ky_0");
	name->Append("_M_");
	*name += mb;
	name->Append("_ED_0");
	
	
	fNormWeight[tKbin][mb] = (TH3F*)wFile->Get(name->Data());
	fNormWeight[tKbin][mb]->SetDirectory(0);
	
	
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
  //Float_t qinvtemp=GetQinv(0,track1, track2);
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
  if(fLinearInterpolation){// Linear Interpolation of osl
    // w interpolation (kt)
    Float_t c000 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexL+1)*wd;
    Float_t c100 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexL+1)*wd;
    Float_t c010 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexL+1)*wd;
    Float_t c001 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexL+1, fQlIndexH+1)*wd;
    Float_t c110 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1)*wd;
    Float_t c101 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1)*wd;
    Float_t c011 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1)*wd;
    Float_t c111 = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*(1-wd) + fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1)*wd;
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
	  farrP1[x][y][z] = fNormWeight[fKtIndexL][fMbin]->GetBinContent(binO,binS,binL);
	  farrP2[x][y][z] = fNormWeight[fKtIndexH][fMbin]->GetBinContent(binO,binS,binL);
	}
      }
    }
    Float_t coord[3]={xd, yd, zd}; 
    Float_t c0 = nCubicInterpolate(3, (Float_t*) farrP1, coord);
    Float_t c1 = nCubicInterpolate(3, (Float_t*) farrP2, coord);
    // kT interpolation
    wgt = c0*(1-wd) + c1*wd;
  }
  ////
  
  // simplified stat error 
  Float_t avgErr = fNormWeight[fKtIndexL][fMbin]->GetBinError(fQoIndexH+1,fQsIndexH+1,fQlIndexH+1);
  avgErr += fNormWeight[fKtIndexH][fMbin]->GetBinError(fQoIndexL+1,fQsIndexL+1,fQlIndexL+1);
  avgErr /= 2.;
  //
  wgtErr = avgErr;
  
 
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
void AliFourPion::SetFSICorrelations(Bool_t legoCase, TH1D *tempss[12], TH1D *tempos[12]){
  // read in 2-particle and 3-particle FSI correlations = K2 & K3
  // 2-particle input histo from file is binned in qinv.  3-particle in qinv of each pair
  if(legoCase){
    cout<<"LEGO call to SetFSICorrelations"<<endl;
    for(Int_t MB=0; MB<12; MB++) {
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
    
    TH1D *temphistoSS[12];
    TH1D *temphistoOS[12];
    for(Int_t MB=0; MB<12; MB++) {
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
    if(fPbPbcase){
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
    }else fFSIindex = 9;// pp and pPb
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
