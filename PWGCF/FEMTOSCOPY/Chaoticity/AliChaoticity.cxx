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

#include "AliChaoticity.h"

#define PI 3.1415927
#define G_Coeff 0.006399 // 2*pi*alpha*M_pion
#define kappa3 0.24 // kappa3 Edgeworth coefficient (non-Gaussian features of C2)
#define kappa4 0.16 // kappa4 Edgeworth coefficient (non-Gaussian features of C2)


// Author: Dhevan Gangadharan

ClassImp(AliChaoticity)

//________________________________________________________________________
AliChaoticity::AliChaoticity():
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
  fPdensityExplicitLoop(kFALSE),
  fPdensityPairCut(kTRUE),
  fTabulatePairs(kFALSE),
  fRMax(11),
  fFixedLambdaBinMomRes(9),
  fFixedLambdaBinr3(10),
  fFilterBit(7),
  fBfield(0),
  fMbin(0),
  fFSIbin(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fQcut(),
  fQLowerCut(0),
  fNormQcutLow(),
  fNormQcutHigh(),
  fKupperBound(0),
  fQupperBound(0),
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
  fMinSepPair(0.035),
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
  fDefaultsCharMult(),
  fDefaultsCharSE(),
  fDefaultsCharME(),
  fDefaultsInt(),
  fPairLocationSE(),
  fPairLocationME(),
  fTripletSkip1(),
  fTripletSkip2(),
  fOtherPairLocation1(),
  fOtherPairLocation2(),
  fNormPairSwitch(),
  fPairSplitCut(),
  fNormPairs(),
  fMomResC2(0x0)
  
{
  // Default constructor
  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2=0x0;
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW = 0x0;
	      
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
		for(Int_t dt=0; dt<kDENtypes; dt++){
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormSmeared = 0x0;
		}//dt
		
	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD = 0x0;
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD = 0x0;
	}
      }
      
    }// ED
  }// Mbin
  
  // Initialize FSI histograms
  for(Int_t i=0; i<2; i++){
    fFSI2SS[i]=0x0; 
    fFSI2OS[i]=0x0;
  }
  for(Int_t i=0; i<6; i++){
    fFSIOmega0SS[i]=0x0; 
    fFSIOmega0OS[i]=0x0;
  }


  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
    }
  }
  

}
//________________________________________________________________________
AliChaoticity::AliChaoticity(const Char_t *name) 
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
  fPdensityExplicitLoop(kFALSE),
  fPdensityPairCut(kTRUE),
  fTabulatePairs(kFALSE),
  fRMax(11),
  fFixedLambdaBinMomRes(9),
  fFixedLambdaBinr3(10),
  fFilterBit(7),
  fBfield(0),
  fMbin(0),
  fFSIbin(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fQcut(),
  fQLowerCut(0),
  fNormQcutLow(),
  fNormQcutHigh(),
  fKupperBound(0),
  fQupperBound(0),
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
  fMinSepPair(0.035),
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
  fDefaultsCharMult(),
  fDefaultsCharSE(),
  fDefaultsCharME(),
  fDefaultsInt(),
  fPairLocationSE(),
  fPairLocationME(),
  fTripletSkip1(),
  fTripletSkip2(),
  fOtherPairLocation1(),
  fOtherPairLocation2(),
  fNormPairSwitch(),
  fPairSplitCut(),
  fNormPairs(),
  fMomResC2(0x0)

{
  // Main constructor
  fAODcase=kTRUE;
  fPdensityExplicitLoop = kFALSE;
  fPdensityPairCut = kTRUE;
  

  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2=0x0;
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW = 0x0;
	      
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
		for(Int_t dt=0; dt<kDENtypes; dt++){
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared = 0x0;
		  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormSmeared = 0x0;
		}//dt
		
	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD = 0x0;
	  KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD = 0x0;
	}
      }
      
    }// ED
  }// Mbin
  
  // Initialize FSI histograms
  for(Int_t i=0; i<2; i++){
    fFSI2SS[i]=0x0; 
    fFSI2OS[i]=0x0;
  }
  for(Int_t i=0; i<6; i++){
    fFSIOmega0SS[i]=0x0; 
    fFSIOmega0OS[i]=0x0;
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
AliChaoticity::AliChaoticity(const AliChaoticity &obj) 
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
    fPdensityExplicitLoop(obj.fPdensityExplicitLoop),
    fPdensityPairCut(obj.fPdensityPairCut),
    fTabulatePairs(obj.fTabulatePairs),
    fRMax(obj.fRMax),
    fFixedLambdaBinMomRes(obj.fFixedLambdaBinMomRes),
    fFixedLambdaBinr3(obj.fFixedLambdaBinr3),
    fFilterBit(obj.fFilterBit),
    fBfield(obj.fBfield),
    fMbin(obj.fMbin),
    fFSIbin(obj.fFSIbin),
    fEDbin(obj.fEDbin),
    fMbins(obj.fMbins),
    fMultLimit(obj.fMultLimit),
    fCentBinLowLimit(obj.fCentBinLowLimit),
    fCentBinHighLimit(obj.fCentBinHighLimit),
    fEventCounter(obj.fEventCounter),
    fEventsToMix(obj.fEventsToMix),
    fZvertexBins(obj.fZvertexBins),
    fMultLimits(),
    fQcut(),
    fQLowerCut(obj.fQLowerCut),
    fNormQcutLow(),
    fNormQcutHigh(),
    fKupperBound(obj.fKupperBound),
    fQupperBound(obj.fQupperBound),
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
    fMinSepPair(obj.fMinSepPair),
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
    fDefaultsCharMult(),
    fDefaultsCharSE(),
    fDefaultsCharME(),
    fDefaultsInt(),
    fPairLocationSE(),
    fPairLocationME(),
    fTripletSkip1(),
    fTripletSkip2(),
    fOtherPairLocation1(),
    fOtherPairLocation2(),
    fNormPairSwitch(),
    fPairSplitCut(),
    fNormPairs(),
    fMomResC2(obj.fMomResC2)
{
  // Copy constructor  
  for(Int_t i=0; i<2; i++){
    fFSI2SS[i]=obj.fFSI2SS[i]; 
    fFSI2OS[i]=obj.fFSI2OS[i];
  }
  for(Int_t i=0; i<6; i++){
    fFSIOmega0SS[i]=obj.fFSIOmega0SS[i]; 
    fFSIOmega0OS[i]=obj.fFSIOmega0OS[i];
  }
  
  // Initialize fNormWeight and fNormWeightErr to 0
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=0x0;
    }
  }
  

}
//________________________________________________________________________
AliChaoticity &AliChaoticity::operator=(const AliChaoticity &obj) 
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
  fPdensityExplicitLoop = obj.fPdensityExplicitLoop;
  fPdensityPairCut = obj.fPdensityPairCut;
  fTabulatePairs = obj.fTabulatePairs;
  fRMax = obj.fRMax;
  fFixedLambdaBinMomRes = obj.fFixedLambdaBinMomRes;
  fFixedLambdaBinr3 = obj.fFixedLambdaBinr3;
  fFilterBit = obj.fFilterBit;
  fBfield = obj.fBfield;
  fMbin = obj.fMbin;
  fFSIbin = obj.fFSIbin;
  fEDbin = obj.fEDbin;
  fMbins = obj.fMbins;
  fMultLimit = obj.fMultLimit;
  fCentBinLowLimit = obj.fCentBinLowLimit;
  fCentBinHighLimit = obj.fCentBinHighLimit;
  fEventCounter = obj.fEventCounter;
  fEventsToMix = obj.fEventsToMix;
  fZvertexBins = obj.fZvertexBins;
  fQLowerCut = obj.fQLowerCut;
  fKupperBound = obj.fKupperBound;
  fQupperBound = obj.fQupperBound;
  fQupperBoundWeights = obj.fQupperBoundWeights;
  fQstep = obj.fQstep;
  fQstepWeights = obj.fQstepWeights;
  fDampStart = obj.fDampStart;
  fDampStep = obj.fDampStep;
  fTPCTOFboundry = obj.fTPCTOFboundry;
  fTOFboundry = obj.fTOFboundry;
  fSigmaCutTPC = obj.fSigmaCutTPC;
  fSigmaCutTOF = obj.fSigmaCutTOF;
  fMinSepPair = obj.fMinSepPair;
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
  fMomResC2 = obj.fMomResC2;

  for(Int_t i=0; i<2; i++){
    fFSI2SS[i]=obj.fFSI2SS[i]; 
    fFSI2OS[i]=obj.fFSI2OS[i];
  }
  for(Int_t i=0; i<6; i++){
    fFSIOmega0SS[i]=obj.fFSIOmega0SS[i]; 
    fFSIOmega0OS[i]=obj.fFSIOmega0OS[i];
  }
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      fNormWeight[i][j]=obj.fNormWeight[i][j];
    }
  }
  
  return (*this);
}
//________________________________________________________________________
AliChaoticity::~AliChaoticity()
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
  if(fMomResC2) delete fMomResC2;
  

  for(Int_t i=0; i<fMultLimit; i++){
    if(fPairLocationSE[i]) delete [] fPairLocationSE[i];
    if(fPairLocationME[i]) delete [] fPairLocationME[i];
    for(Int_t j=0; j<2; j++){
      if(fOtherPairLocation1[j][i]) delete [] fOtherPairLocation1[j][i];
      if(fOtherPairLocation2[j][i]) delete [] fOtherPairLocation2[j][i];
    }
    for(Int_t j=0; j<3; j++) if(fNormPairSwitch[j][i]) delete [] fNormPairSwitch[j][i];
    for(Int_t j=0; j<4; j++) if(fPairSplitCut[j][i]) delete [] fPairSplitCut[j][i];
  }
  for(Int_t i=0; i<kPairLimit; i++) if(fTripletSkip1[i]) delete [] fTripletSkip1[i];
  for(Int_t i=0; i<2*kPairLimit; i++) if(fTripletSkip2[i]) delete [] fTripletSkip2[i];
  for(Int_t i=0; i<3; i++) if(fNormPairs[i]) delete [] fNormPairs[i];
  //
  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2;
	      
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW;
	      //
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW;
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Terms) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Terms;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Terms) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Terms;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsIdeal) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsIdeal;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsIdeal) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsIdeal;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSmeared) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSmeared;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSmeared) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSmeared;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Q3W) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Q3W;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Q3W) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Q3W;
		//
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK3;
		//
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK2) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK2;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK2) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK2;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK2) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK2;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK2) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK2;

		//
		for(Int_t dt=0; dt<kDENtypes; dt++){
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared;
		  if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormSmeared;

		}//dt
		
	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
      for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
	for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	  if(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD) delete KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD;
	  if(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD) delete KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD;
	}
      }
      
    }// ED
  }// Mbin
  
  if(fMomResC2) delete fMomResC2;
 
  for(Int_t i=0; i<2; i++){
    if(fFSI2SS[i]) delete fFSI2SS[i]; 
    if(fFSI2OS[i]) delete fFSI2OS[i];
  }
  for(Int_t i=0; i<6; i++){
    if(fFSIOmega0SS[i]) delete fFSIOmega0SS[i]; 
    if(fFSIOmega0OS[i]) delete fFSIOmega0OS[i];
  }
  for(Int_t i=0; i<3; i++){// Kt iterator
    for(Int_t j=0; j<10; j++){// Mbin iterator
      if(fNormWeight[i][j]) delete fNormWeight[i][j];
    }
  }
  
}
//________________________________________________________________________
void AliChaoticity::ParInit()
{
  cout<<"AliChaoticity MyInit() call"<<endl;
  cout<<"lego:"<<fLEGO<<"  MCcase:"<<fMCcase<<"  PbPbcase:"<<fPbPbcase<<"  TabulatePairs:"<<fTabulatePairs<<"  GenSignal:"<<fGenerateSignal<<"  CentLow:"<<fCentBinLowLimit<<"  CentHigh:"<<fCentBinHighLimit<<"  RMax:"<<fRMax<<"  LambdaBinMomRes:"<<fFixedLambdaBinMomRes<<"  LambdaBinr3:"<<fFixedLambdaBinr3<<"  FB:"<<fFilterBit<<"  MinPairSep:"<<fMinSepPair<<"  NsigTPC:"<<fSigmaCutTPC<<"  NsigTOF:"<<fSigmaCutTOF<<endl;

  fRandomNumber = new TRandom3();
  fRandomNumber->SetSeed(0);
    
  //
  fEventCounter=0;
  if(fPdensityExplicitLoop) fEventsToMix=3;
  else if(fPdensityPairCut && !fPdensityExplicitLoop) fEventsToMix=2;
  else fEventsToMix=0;
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
    fQcut[0]=0.1;//pi-pi, pi-k, pi-p
    fQcut[1]=0.1;//k-k
    fQcut[2]=0.6;//the rest
    fNormQcutLow[0] = 0.15;//0.15
    fNormQcutHigh[0] = 0.175;//0.175
    fNormQcutLow[1] = 1.34;//1.34
    fNormQcutHigh[1] = 1.4;//1.4
    fNormQcutLow[2] = 1.1;//1.1
    fNormQcutHigh[2] = 1.4;//1.4
  }
  else {// pp
    fMultLimit=kMultLimitpp; 
    fMbins=kMultBinspp; 
    fQcut[0]=0.6;
    fQcut[1]=0.6;
    fQcut[2]=0.6;
    fNormQcutLow[0] = 1.0;
    fNormQcutHigh[0] = 1.5;
    fNormQcutLow[1] = 1.0;
    fNormQcutHigh[1] = 1.5;
    fNormQcutLow[2] = 1.0;
    fNormQcutHigh[2] = 1.5;
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
  fQupperBound = 0.1;
  fQstep = fQupperBound/Float_t(kQbins);
  fQstepWeights = fQupperBoundWeights/Float_t(kQbinsWeights);
  for(Int_t i=0; i<kQbinsWeights; i++) {fQmean[i]=(i+0.5)*fQstepWeights;}
  //
  fDampStart = 0.5;// was 0.3
  fDampStep = 0.02;
  
  //

  
  fEC = new AliChaoticityEventCollection **[fZvertexBins];
  for(UShort_t i=0; i<fZvertexBins; i++){
    
    fEC[i] = new AliChaoticityEventCollection *[fMbins];

    for(UShort_t j=0; j<fMbins; j++){
      
      fEC[i][j] = new AliChaoticityEventCollection(fEventsToMix+1, fMultLimit, kPairLimit, kMCarrayLimit, fMCcase);
    }
  }
  
    
  for(Int_t i=0; i<fMultLimit; i++) fDefaultsCharMult[i]='0';
  for(Int_t i=0; i<kPairLimit; i++) fDefaultsCharSE[i]='0';
  for(Int_t i=0; i<2*kPairLimit; i++) fDefaultsCharME[i]='0';
  for(Int_t i=0; i<fMultLimit; i++) fDefaultsInt[i]=-1;
  for(Int_t i=0; i<fMultLimit; i++) fPairLocationSE[i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<fMultLimit; i++) fPairLocationME[i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<kPairLimit; i++) fTripletSkip1[i] = new TArrayC(fMultLimit,fDefaultsCharSE);
  for(Int_t i=0; i<2*kPairLimit; i++) fTripletSkip2[i] = new TArrayC(fMultLimit,fDefaultsCharME);
 

  // Normalization utilities
  for(Int_t i=0; i<fMultLimit; i++) fOtherPairLocation1[0][i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<fMultLimit; i++) fOtherPairLocation1[1][i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<fMultLimit; i++) fOtherPairLocation2[0][i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<fMultLimit; i++) fOtherPairLocation2[1][i] = new TArrayI(fMultLimit,fDefaultsInt);
  for(Int_t i=0; i<fMultLimit; i++) fNormPairSwitch[0][i] = new TArrayC(fMultLimit,fDefaultsCharMult);
  for(Int_t i=0; i<fMultLimit; i++) fNormPairSwitch[1][i] = new TArrayC(fMultLimit,fDefaultsCharMult);
  for(Int_t i=0; i<fMultLimit; i++) fNormPairSwitch[2][i] = new TArrayC(fMultLimit,fDefaultsCharMult);
 
  // Track Merging/Splitting utilities
  for(Int_t i=0; i<fMultLimit; i++) fPairSplitCut[0][i] = new TArrayC(fMultLimit,fDefaultsCharMult);// P11
  for(Int_t i=0; i<fMultLimit; i++) fPairSplitCut[1][i] = new TArrayC(fMultLimit,fDefaultsCharMult);// P12
  for(Int_t i=0; i<fMultLimit; i++) fPairSplitCut[2][i] = new TArrayC(fMultLimit,fDefaultsCharMult);// P13
  for(Int_t i=0; i<fMultLimit; i++) fPairSplitCut[3][i] = new TArrayC(fMultLimit,fDefaultsCharMult);// P23

  
  fNormPairs[0] = new AliChaoticityNormPairStruct[kNormPairLimit];
  fNormPairs[1] = new AliChaoticityNormPairStruct[kNormPairLimit];
  

  fTempStruct = new AliChaoticityTrackStruct[fMultLimit];
   
   
  fTrueMassP=0.93827, fTrueMassPi=0.13957, fTrueMassK=0.493677, fTrueMassKs=0.497614, fTrueMassLam=1.11568;

  

  // Set weights, Coulomb corrections, and Momentum resolution corrections manually if not on LEGO
  if(!fLEGO) {
    SetFSICorrelations(fLEGO);// Read in 2-particle and 3-particle FSI correlations
    if(!fTabulatePairs) SetWeightArrays(fLEGO);// Set Weight Array
    if(!fMCcase && !fTabulatePairs) SetMomResCorrections(fLEGO);// Read Momentum resolution file
    //if(!fTabulatePairs) SetMomResCorrections(fLEGO);// Read Momentum resolution file
  }
  
  /////////////////////////////////////////////
  /////////////////////////////////////////////
  
}
//________________________________________________________________________
void AliChaoticity::UserCreateOutputObjects()
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
 
  TH1F *fRejectedPairs = new TH1F("fRejectedPairs","",200,0,2);
  fOutputList->Add(fRejectedPairs);
  TH1I *fRejectedEvents = new TH1I("fRejectedEvents","",fMbins,0.5,fMbins+.5);
  fOutputList->Add(fRejectedEvents);
  
  TH3F *fPairsDetaDPhiNum = new TH3F("fPairsDetaDPhiNum","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiNum);
  TH3F *fPairsDetaDPhiDen = new TH3F("fPairsDetaDPhiDen","",10,-.5,9.5, 200,-0.2,0.2, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsDetaDPhiDen);
 
  TH2D *fResonanceOSPairs = new TH2D("fResonanceOSPairs","",fMbins,.5,fMbins+.5, 1000,0,2);
  if(fMCcase) fOutputList->Add(fResonanceOSPairs);
  TH2D *fAllOSPairs = new TH2D("fAllOSPairs","",fMbins,.5,fMbins+.5, 1000,0,2);
  if(fMCcase) fOutputList->Add(fAllOSPairs);
  
  TProfile *fAvgMult = new TProfile("fAvgMult","",fMbins,.5,fMbins+.5, 0,1500,"");
  fOutputList->Add(fAvgMult);

  TH3D *fTPNRejects1 = new TH3D("fTPNRejects1","",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
  fOutputList->Add(fTPNRejects1);
  TH3D *fTPNRejects2 = new TH3D("fTPNRejects2","",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
  fOutputList->Add(fTPNRejects2);
  TH3D *fTPNRejects3 = new TH3D("fTPNRejects3","",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
  fOutputList->Add(fTPNRejects3);
  TH3D *fTPNRejects4 = new TH3D("fTPNRejects4","",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
  fOutputList->Add(fTPNRejects4);
  TH3D *fTPNRejects5 = new TH3D("fTPNRejects5","",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
  fOutputList->Add(fTPNRejects5);


  TH3D *fKt3DistTerm1 = new TH3D("fKt3DistTerm1","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  TH3D *fKt3DistTerm5 = new TH3D("fKt3DistTerm5","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  fOutputList->Add(fKt3DistTerm1);
  fOutputList->Add(fKt3DistTerm5);

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


  if(fPdensityExplicitLoop || fPdensityPairCut){
    
    for(Int_t mb=0; mb<fMbins; mb++){
      if((mb < fCentBinLowLimit) || (mb > fCentBinHighLimit)) continue;

      for(Int_t edB=0; edB<fEDbins; edB++){
	for(Int_t c1=0; c1<2; c1++){
	  for(Int_t c2=0; c2<2; c2++){
	    for(Int_t sc=0; sc<kSCLimit2; sc++){
	      for(Int_t term=0; term<2; term++){
		
		TString *nameEx2 = new TString("Explicit2_Charge1_");
		*nameEx2 += c1;
		nameEx2->Append("_Charge2_");
		*nameEx2 += c2;
		nameEx2->Append("_SC_");
		*nameEx2 += sc;
		nameEx2->Append("_M_");
		*nameEx2 += mb;
		nameEx2->Append("_ED_");
		*nameEx2 += edB;
		nameEx2->Append("_Term_");
		*nameEx2 += term+1;
		
		if(sc==0 || sc==3 || sc==5){
		  if( (c1+c2)==1 ) {if(c1!=0) continue;}// skip degenerate histogram
		}
		
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2 = new TH2D(nameEx2->Data(),"Two Particle Distribution",20,0,1, 400,0,2);
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2);
		TString *nameEx2QW=new TString(nameEx2->Data());
		nameEx2QW->Append("_QW");
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW = new TH2D(nameEx2QW->Data(),"Two Particle Distribution",20,0,1, 400,0,2);
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW);
		TString *nameAvgP=new TString(nameEx2->Data());
		nameAvgP->Append("_AvgP");
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP = new TProfile2D(nameAvgP->Data(),"",10,0,1, 400,0,2, 0.,1.0,"");
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP);
		
		// Momentum resolution histos
		if(fMCcase && sc==0){
		  TString *nameIdeal = new TString(nameEx2->Data());
		  nameIdeal->Append("_Ideal");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = new TH2D(nameIdeal->Data(),"Two Particle Distribution",fRVALUES*kNDampValues,-0.5,fRVALUES*kNDampValues-0.5, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal);
		  TString *nameSmeared = new TString(nameEx2->Data());
		  nameSmeared->Append("_Smeared");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = new TH2D(nameSmeared->Data(),"Two Particle Distribution",fRVALUES*kNDampValues,-0.5,fRVALUES*kNDampValues-0.5, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared);
		  //
		  TString *nameEx2MC=new TString(nameEx2->Data());
		  nameEx2MC->Append("_MCqinv");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv = new TH1D(nameEx2MC->Data(),"",400,0,2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv);
		  TString *nameEx2MCQW=new TString(nameEx2->Data());
		  nameEx2MCQW->Append("_MCqinvQW");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = new TH1D(nameEx2MCQW->Data(),"",400,0,2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW);
		  //
		  TString *nameEx2PIDpurityDen=new TString(nameEx2->Data());
		  nameEx2PIDpurityDen->Append("_PIDpurityDen");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = new TH2D(nameEx2PIDpurityDen->Data(),"Two Particle Distribution",20,0,1, 400,0,2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen);
		  TString *nameEx2PIDpurityNum=new TString(nameEx2->Data());
		  nameEx2PIDpurityNum->Append("_PIDpurityNum");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = new TH2D(nameEx2PIDpurityNum->Data(),"Two Particle Distribution",20,0,1, 400,0,2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum);
		}
		if(sc==0){
		  
		  TString *nameEx2OSLB1 = new TString(nameEx2->Data()); 
		  nameEx2OSLB1->Append("_osl_b1");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSL);
		  nameEx2OSLB1->Append("_QW");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW = new TH3D(nameEx2OSLB1->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[0].fExplicit2OSLQW);
		  //
		  TString *nameEx2OSLB2 = new TString(nameEx2->Data()); 
		  nameEx2OSLB2->Append("_osl_b2");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSL);
		  nameEx2OSLB2->Append("_QW");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW = new TH3D(nameEx2OSLB2->Data(),"Two Particle Distribution",kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].OSL_ktbin[1].fExplicit2OSLQW);
		  
		}

	      }// term_2
	    }// SC_2
	    
	    // skip 3-particle if Tabulate6DPairs is true
	    if(fTabulatePairs) continue;
 
	    for(Int_t c3=0; c3<2; c3++){
	      for(Int_t sc=0; sc<kSCLimit3; sc++){
		for(Int_t term=0; term<5; term++){
		  TString *nameEx3 = new TString("Explicit3_Charge1_");
		  *nameEx3 += c1;
		  nameEx3->Append("_Charge2_");
		  *nameEx3 += c2;
		  nameEx3->Append("_Charge3_");
		  *nameEx3 += c3;
		  nameEx3->Append("_SC_");
		  *nameEx3 += sc;
		  nameEx3->Append("_M_");
		  *nameEx3 += mb;
		  nameEx3->Append("_ED_");
		  *nameEx3 += edB;
		  nameEx3->Append("_Term_");
		  *nameEx3 += term+1;
		  
		  TString *namePC3 = new TString("PairCut3_Charge1_");
		  *namePC3 += c1;
		  namePC3->Append("_Charge2_");
		  *namePC3 += c2;
		  namePC3->Append("_Charge3_");
		  *namePC3 += c3;
		  namePC3->Append("_SC_");
		  *namePC3 += sc;
		  namePC3->Append("_M_");
		  *namePC3 += mb;
		  namePC3->Append("_ED_");
		  *namePC3 += edB;
		  namePC3->Append("_Term_");
		  *namePC3 += term+1;
	      
		  ///////////////////////////////////////
		  // skip degenerate histograms
		  if(sc==0 || sc==6 || sc==9){// Identical species
		    if( (c1+c2+c3)==1) {if(c3!=1) continue;}
		    if( (c1+c2+c3)==2) {if(c1!=0) continue;}
		  }else if(sc!=5){
		    if( (c1+c2)==1) {if(c1!=0) continue;}
		  }else {}// do nothing for pi-k-p case
		  
		  /////////////////////////////////////////
	      
		  

		  if(fPdensityExplicitLoop){
		    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3 = new TH1D(nameEx3->Data(),"Three Particle Distribution",200,0,2);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fExplicit3);
		    //
		    nameEx3->Append("_Norm");
		    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3 = new TH1D(nameEx3->Data(),"Explicit_3 Norm",1,-0.5,0.5);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNormEx3);
		  }
		  if(fPdensityPairCut){
		    TString *nameNorm=new TString(namePC3->Data());
		    nameNorm->Append("_Norm");
		    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = new TH1D(nameNorm->Data(),"Norm",1,-0.5,0.5);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3);
		    //
		    if(sc<=2){
		      TString *name3DQ=new TString(namePC3->Data());
		      name3DQ->Append("_3D");
		      Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = new TH3D(name3DQ->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
		      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3);
		      //
		      
		      const int NEdgesPos=16;
		      double lowEdges4vectPos[NEdgesPos]={0};
		      lowEdges4vectPos[0]=0.0;
		      lowEdges4vectPos[1]=0.00005;// best resolution at low Q^2
		      for(int edge=2; edge<NEdgesPos; edge++){
			lowEdges4vectPos[edge] = lowEdges4vectPos[edge-1] + lowEdges4vectPos[1]*(edge);
		      }
		      const int NEdges=2*NEdgesPos-1;
		      double lowEdges4vect[NEdges]={0};
		      for(int edge=0; edge<NEdges; edge++){
			if(edge<NEdgesPos-1) lowEdges4vect[edge] = -lowEdges4vectPos[NEdgesPos-1-edge];
			else if(edge==NEdgesPos-1) lowEdges4vect[edge] = 0;
			else lowEdges4vect[edge] = lowEdges4vectPos[edge-NEdgesPos+1];
			//if(c1==c2 && c1==c3) cout<<lowEdges4vect[edge]<<endl;
		      }
		      
		      /*
		      const int NEdgesPos=16;
		      double lowEdges4vectPos[NEdgesPos]={0};
		      lowEdges4vectPos[0]=0.0;
		      lowEdges4vectPos[1]=0.0002;// was 0.0005, then 0.0002
		      for(int edge=2; edge<NEdgesPos; edge++){
			lowEdges4vectPos[edge] = lowEdges4vectPos[edge-1] + lowEdges4vectPos[1];
		      }
		      const int NEdges=2*NEdgesPos-1;
		      double lowEdges4vect[NEdges]={0};
		      for(int edge=0; edge<NEdges; edge++){
			if(edge<NEdgesPos-1) lowEdges4vect[edge] = -lowEdges4vectPos[NEdgesPos-1-edge];
			else if(edge==NEdgesPos-1) lowEdges4vect[edge] = 0;
			else lowEdges4vect[edge] = lowEdges4vectPos[edge-NEdgesPos+1];
		      }
		      */
		      if(c1==c2 && c1==c3 && sc==0 && fMCcase==kFALSE){
			TString *name4vect1=new TString(namePC3->Data());
			TString *name4vect2=new TString(namePC3->Data());
			name4vect1->Append("_4VectProd1");
			name4vect2->Append("_4VectProd2");
			// use 3.75e6 MeV^4 as the resolution on QprodSum
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Terms = new TH3D(name4vect1->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Terms);
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Terms = new TH3D(name4vect2->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Terms);
		      }
		      if(sc==0 && fMCcase==kTRUE){
			TString *name3DMomResIdeal=new TString(namePC3->Data());
			name3DMomResIdeal->Append("_Ideal");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal = new TH3D(name3DMomResIdeal->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal);
			TString *name3DMomResSmeared=new TString(namePC3->Data());
			name3DMomResSmeared->Append("_Smeared");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared = new TH3D(name3DMomResSmeared->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared);
			//
			TString *name3DMomResQW12=new TString(namePC3->Data());
			name3DMomResQW12->Append("_QW12");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fQW12 = new TH3D(name3DMomResQW12->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fQW12);
			TString *name3DMomResQW13=new TString(namePC3->Data());
			name3DMomResQW13->Append("_QW13");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fQW13 = new TH3D(name3DMomResQW13->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fQW13);
			//
			if(term==0){
			  TString *name3DSumK3=new TString(namePC3->Data());
			  name3DSumK3->Append("_SumK3");
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSumK3 = new TH3D(name3DSumK3->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSumK3);
			  TString *name3DEnK3=new TString(namePC3->Data());
			  name3DEnK3->Append("_EnK3");
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fEnK3 = new TH3D(name3DEnK3->Data(),"", kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fEnK3);
			}

			if(c1==c2 && c1==c3){
			  TString *name4vect1Ideal=new TString(namePC3->Data());
			  TString *name4vect1Smeared=new TString(namePC3->Data());
			  TString *name4vect2Ideal=new TString(namePC3->Data());
			  TString *name4vect2Smeared=new TString(namePC3->Data());
			  TString *name4vect1Q3W=new TString(namePC3->Data());
			  TString *name4vect2Q3W=new TString(namePC3->Data());
			  name4vect1Ideal->Append("_4VectProd1Ideal");
			  name4vect1Smeared->Append("_4VectProd1Smeared");
			  name4vect2Ideal->Append("_4VectProd2Ideal");
			  name4vect2Smeared->Append("_4VectProd2Smeared");
			  name4vect1Q3W->Append("_4VectProd1Q3W");
			  name4vect2Q3W->Append("_4VectProd2Q3W");
			  // 
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsIdeal = new TH3D(name4vect1Ideal->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsIdeal);
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSmeared = new TH3D(name4vect1Smeared->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSmeared);
			  //
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsIdeal = new TH3D(name4vect2Ideal->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsIdeal);
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSmeared = new TH3D(name4vect2Smeared->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSmeared);
			  //
			  if(term==0){// average Q3 in each FVP cell
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Q3W = new TH3D(name4vect1Q3W->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1Q3W);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Q3W = new TH3D(name4vect2Q3W->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2Q3W);
			  }
			  //
			  if(term==0){
			    TString *name4vect1SumK3=new TString(namePC3->Data());
			    TString *name4vect2SumK3=new TString(namePC3->Data());
			    TString *name4vect1EnK3=new TString(namePC3->Data());
			    TString *name4vect2EnK3=new TString(namePC3->Data());
			    name4vect1SumK3->Append("_4VectProd1SumK3");
			    name4vect2SumK3->Append("_4VectProd2SumK3");
			    name4vect1EnK3->Append("_4VectProd1EnK3");
			    name4vect2EnK3->Append("_4VectProd2EnK3");
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK3 = new TH3D(name4vect1SumK3->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK3);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK3 = new TH3D(name4vect2SumK3->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK3);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK3 = new TH3D(name4vect1EnK3->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK3);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK3 = new TH3D(name4vect2EnK3->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK3);
			  }// term 0
			  if(term > 0 && term < 4){
			    TString *name4vect1SumK2=new TString(namePC3->Data());
			    TString *name4vect2SumK2=new TString(namePC3->Data());
			    TString *name4vect1EnK2=new TString(namePC3->Data());
			    TString *name4vect2EnK2=new TString(namePC3->Data());
			    name4vect1SumK2->Append("_4VectProd1SumK2");
			    name4vect2SumK2->Append("_4VectProd2SumK2");
			    name4vect1EnK2->Append("_4VectProd1EnK2");
			    name4vect2EnK2->Append("_4VectProd2EnK2");
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK2 = new TH3D(name4vect1SumK2->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsSumK2);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK2 = new TH3D(name4vect2SumK2->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsSumK2);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK2 = new TH3D(name4vect1EnK2->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd1TermsEnK2);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK2 = new TH3D(name4vect2EnK2->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].f4VectProd2TermsEnK2);
			  }// terms 1,2,3
			}
		      }// MCcase
		      //
		      if(c1==c2 && c1==c3 && term==4 && sc==0){
			for(Int_t dt=0; dt<kDENtypes; dt++){
			  TString *nameDenType=new TString("PairCut3_Charge1_");
			  *nameDenType += c1;
			  nameDenType->Append("_Charge2_");
			  *nameDenType += c2;
			  nameDenType->Append("_Charge3_");
			  *nameDenType += c3;
			  nameDenType->Append("_SC_");
			  *nameDenType += sc;
			  nameDenType->Append("_M_");
			  *nameDenType += mb;
			  nameDenType->Append("_ED_");
			  *nameDenType += edB;
			  nameDenType->Append("_TPN_");
			  *nameDenType += dt;
			  
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm = new TH3D(nameDenType->Data(),"",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNorm);
			  // neglect errors for TPN
			  //nameDenType->Append("_Err");
			  //Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNormErr = new TH3D(nameDenType->Data(),"",kQbins,0,fQupperBound, kQbins,0,fQupperBound, kQbins,0,fQupperBound);
			  //fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].fTwoPartNormErr);
			  //
			  TString *name4vect1TPN=new TString(nameDenType->Data());
			  TString *name4vect2TPN=new TString(nameDenType->Data());
			  name4vect1TPN->Append("_4VectProd1");
			  name4vect2TPN->Append("_4VectProd2");
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm = new TH3D(name4vect1TPN->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNorm);
			  Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm = new TH3D(name4vect2TPN->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			  fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNorm);
			  //
			  if(fMCcase){
			    TString *name4vect1TPNIdeal=new TString(nameDenType->Data());
			    TString *name4vect2TPNIdeal=new TString(nameDenType->Data());
			    TString *name4vect1TPNSmeared=new TString(nameDenType->Data());
			    TString *name4vect2TPNSmeared=new TString(nameDenType->Data());
			    name4vect1TPNIdeal->Append("_4VectProd1Ideal");
			    name4vect2TPNIdeal->Append("_4VectProd2Ideal");
			    name4vect1TPNSmeared->Append("_4VectProd1Smeared");
			    name4vect2TPNSmeared->Append("_4VectProd2Smeared");
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal = new TH3D(name4vect1TPNIdeal->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormIdeal);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal = new TH3D(name4vect2TPNIdeal->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormIdeal);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared = new TH3D(name4vect1TPNSmeared->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd1TwoPartNormSmeared);
			    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormSmeared = new TH3D(name4vect2TPNSmeared->Data(),"",NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect, NEdges-1,lowEdges4vect);
			    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].DT[dt].f4VectProd2TwoPartNormSmeared);
			  }

			}
					
		      }// term=4
		    }// c and sc exclusion
		  }// PdensityPairCut
		}// term_3
	      }// SC_3
	    }//c3
	  }//c2
	}//c1
      }// ED
    }// mbin
  }// Pdensity Method

  
  if(fTabulatePairs){
    
    for(Int_t tKbin=0; tKbin<fKbinsT; tKbin++){
      for(Int_t yKbin=0; yKbin<fKbinsY; yKbin++){
	for(Int_t mb=0; mb<fMbins; mb++){
	  for(Int_t edB=0; edB<fEDbins; edB++){
	      
	      TString *nameNum = new TString("TwoPart_num_Kt_");
	      *nameNum += tKbin;
	      nameNum->Append("_Ky_");
	      *nameNum += yKbin;
	      nameNum->Append("_M_");
	      *nameNum += mb;
	      nameNum->Append("_ED_");
	      *nameNum += edB;

	      TString *nameDen = new TString("TwoPart_den_Kt_");
	      *nameDen += tKbin;
	      nameDen->Append("_Ky_");
	      *nameDen += yKbin;
	      nameDen->Append("_M_");
	      *nameDen += mb;
	      nameDen->Append("_ED_");
	      *nameDen += edB;


	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD = new TH3D(nameNum->Data(),"", kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	      fOutputList->Add(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[0].fExplicit2ThreeD);
	      
	      KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD = new TH3D(nameDen->Data(),"", kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights, kQbinsWeights,0,fQupperBoundWeights);
	      fOutputList->Add(KT[tKbin].KY[yKbin].MB[mb].EDB[edB].TwoPT[1].fExplicit2ThreeD);
	    }
	  }
	}
      }
 
  }
  
  
  TProfile *fQsmearMean = new TProfile("fQsmearMean","",2,0.5,2.5, -0.2,0.2,"");
  fOutputList->Add(fQsmearMean);
  TProfile *fQsmearSq = new TProfile("fQsmearSq","",2,0.5,2.5, -2,2,"");
  fOutputList->Add(fQsmearSq);
  TH1D *fQDist = new TH1D("fQDist","",200,-.2,.2);
  fOutputList->Add(fQDist);
  
  

  ////////////////////////////////////
  ///////////////////////////////////  
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliChaoticity::Exec(Option_t *) 
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
    

    
    
    ////////////////////////////////
    // Vertexing
    ((TH1F*)fOutputList->FindObject("fMultDist1"))->Fill(fAOD->GetNumberOfTracks());
    primaryVertexAOD = fAOD->GetPrimaryVertex();
    vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
    
    if(fabs(vertex[2]) > 10) {cout<<"Zvertex Out of Range. Skip Event"<<endl; return;} // Z-Vertex Cut 
    ((TH3F*)fOutputList->FindObject("fVertexDist"))->Fill(vertex[0], vertex[1], vertex[2]);
    
    if(fAOD->IsPileupFromSPD()) {cout<<"PileUpEvent. Skip Event"<<endl; return;} // Reject Pile-up events
    if(primaryVertexAOD->GetNContributors() < 1) {cout<<"Bad Vertex. Skip Event"<<endl; return;}
   
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
      
      if(aodtrack->Pt() < 0.16) continue;
      if(fabs(aodtrack->Eta()) > 0.8) continue;
           
     
      Bool_t goodMomentum = aodtrack->GetPxPyPz( fTempStruct[myTracks].fP);
      if(!goodMomentum) continue; 
      aodtrack->GetXYZ( fTempStruct[myTracks].fX);
            
      Float_t dca2[2];
      Float_t dca3d;

      dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - vertex[0],2) + pow(fTempStruct[myTracks].fX[1] - vertex[1],2));
      dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - vertex[2],2));
      dca3d = sqrt( pow(dca2[0],2) + pow(dca2[1],2));
             
      fTempStruct[myTracks].fStatus = status;
      fTempStruct[myTracks].fFiltermap = aodtrack->GetFilterMap();
      fTempStruct[myTracks].fId = aodtrack->GetID();
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
      
    
      
      if(fTempStruct[myTracks].fMom > 0.9999) continue;// upper P bound
      if(fTempStruct[myTracks].fPt > 0.9999) continue;// upper P bound
      if(fTempStruct[myTracks].fP[2] > 0.9999) continue;// upper P bound

      if(fTempStruct[myTracks].fCharge==+1) {
	((TH2F*)fOutputList->FindObject("fDCAxyDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }else {
	((TH2F*)fOutputList->FindObject("fDCAxyDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }
     
      ((TH3F*)fOutputList->FindObject("fPhiPtDist"))->Fill(aodtrack->Charge(), aodtrack->Phi(), aodtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fPtEtaDist"))->Fill(aodtrack->Charge(), aodtrack->Pt(), aodtrack->Eta());

            
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

      if(fFilterBit != 7) {
	nSigmaTPC[0]=fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kElectron));
	nSigmaTPC[1]=fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kMuon));
	nSigmaTPC[2]=fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kPion));
	nSigmaTPC[3]=fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kKaon));
	nSigmaTPC[4]=fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kProton));
	//
	nSigmaTOF[0]=fabs(fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kElectron));
	nSigmaTOF[1]=fabs(fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kMuon));
	nSigmaTOF[2]=fabs(fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kPion));
	nSigmaTOF[3]=fabs(fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kKaon));
	nSigmaTOF[4]=fabs(fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kProton));
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
	  
	  nSigmaTPC[0]=fabs(fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kElectron));
	  nSigmaTPC[1]=fabs(fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kMuon));
	  nSigmaTPC[2]=fabs(fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kPion));
	  nSigmaTPC[3]=fabs(fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kKaon));
	  nSigmaTPC[4]=fabs(fPIDResponse->NumberOfSigmasTPC(aodTrack2,AliPID::kProton));
	  //
	  nSigmaTOF[0]=fabs(fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kElectron));
	  nSigmaTOF[1]=fabs(fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kMuon));
	  nSigmaTOF[2]=fabs(fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kPion));
	  nSigmaTOF[3]=fabs(fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kKaon));
	  nSigmaTOF[4]=fabs(fPIDResponse->NumberOfSigmasTOF(aodTrack2,AliPID::kProton));
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
	if(nSigmaTOF[0]<fSigmaCutTOF) fTempStruct[myTracks].fElectron = kTRUE;// Electron candidate
	if(nSigmaTOF[2]<fSigmaCutTOF) fTempStruct[myTracks].fPion = kTRUE;// Pion candidate
	if(nSigmaTOF[3]<fSigmaCutTOF) fTempStruct[myTracks].fKaon = kTRUE;// Kaon candidate
	if(nSigmaTOF[4]<fSigmaCutTOF) fTempStruct[myTracks].fProton = kTRUE;// Proton candidate
      }else {// TPC info instead
	if(nSigmaTPC[0]<fSigmaCutTPC) fTempStruct[myTracks].fElectron = kTRUE;// Electron candidate
	if(nSigmaTPC[2]<fSigmaCutTPC) fTempStruct[myTracks].fPion = kTRUE;// Pion candidate
	if(nSigmaTPC[3]<fSigmaCutTPC) fTempStruct[myTracks].fKaon = kTRUE;// Kaon candidate
	if(nSigmaTPC[4]<fSigmaCutTPC) fTempStruct[myTracks].fProton = kTRUE;// Proton candidate
      }
               
      
      // Ensure there is only 1 candidate per track
      if(fTempStruct[myTracks].fElectron && fTempStruct[myTracks].fMom < 0.45) continue;// Remove electron band
      if(!fTempStruct[myTracks].fPion && !fTempStruct[myTracks].fKaon && !fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;
      ////////////////////////
      if(fTempStruct[myTracks].fProton && fTempStruct[myTracks].fMom < 0.25) continue;//extra cut for protons

      if(!fTempStruct[myTracks].fPion) continue;// only pions
      
          
           
    
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
      
    }
  }else {// ESD tracks
    cout<<"ESDs not supported currently"<<endl;
    return;
  }
  
  
  if(myTracks >= 1) {
    ((TH1F*)fOutputList->FindObject("fMultDist3"))->Fill(myTracks);
  }
 
 
  //cout<<"There are "<<myTracks<<"  myTracks"<<endl;
  //cout<<"pionCount = "<<pionCount<<"   kaonCount = "<<kaonCount<<"   protonCount = "<<protonCount<<endl;

  /////////////////////////////////////////
  // Pion Multiplicity Cut (To ensure all Correlation orders are present in each event)
  if(myTracks < 3) {cout<<"Less than 3 tracks. Skipping Event."<<endl; return;}
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
  
  fFSIbin=0;
  if(fMbin==0) fFSIbin = 0;//0-5%
  else if(fMbin==1) fFSIbin = 1;//5-10%
  else if(fMbin<=3) fFSIbin = 2;//10-20%
  else if(fMbin<=5) fFSIbin = 3;//20-30%
  else if(fMbin<=7) fFSIbin = 4;//30-40%
  else fFSIbin = 5;//40-50%

  Int_t rIndexForTPNMomRes = fRMax-6;
  //if(fMbin<=1) {rIndexForTPNMomRes=fRMax;}
  //else if(fMbin<=3) {rIndexForTPNMomRes=fRMax-1;}
  //else if(fMbin<=5) {rIndexForTPNMomRes=fRMax-2;}
  //else {rIndexForTPNMomRes=fRMax-3;}
  if(fMbin==0) {rIndexForTPNMomRes=fRMax-6;}// 10 fm with EW (fRMax should be 11 for normal running)
  else if(fMbin==1) {rIndexForTPNMomRes=fRMax-7;}
  else if(fMbin<=3) {rIndexForTPNMomRes=fRMax-8;}
  else if(fMbin<=5) {rIndexForTPNMomRes=fRMax-9;}
  else {rIndexForTPNMomRes=fRMax-10;}

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
      }	
    }
  }
    
  
  
  Float_t qinv12=0, qinv13=0, qinv23=0;
  Float_t qinv12Flat=0;
  Float_t qout=0, qside=0, qlong=0;
  Float_t qoutFlat=0, qsideFlat=0, qlongFlat=0;
  Float_t qoutMC=0, qsideMC=0, qlongMC=0;
  Float_t firstQ=0, secondQ=0, thirdQ=0;
  Float_t firstQMC=0, secondQMC=0, thirdQMC=0;
  Float_t transK12=0, rapK12=0, transK3=0;
  Int_t transKbin=0, rapKbin=0;
  Float_t q3=0, q3MC=0;
  Int_t ch1=0, ch2=0, ch3=0;
  Short_t key1=0, key2=0, key3=0;
  Int_t bin1=0, bin2=0, bin3=0;
  Float_t pVect1[4]={0}; 
  Float_t pVect2[4]={0};
  Float_t pVect3[4]={0}; 
  Float_t pVect1MC[4]={0}; 
  Float_t pVect2MC[4]={0};
  Float_t pVect3MC[4]={0};
  Float_t pVect2Flat[4]={0};
  Float_t pVect3Flat[4]={0};
  Int_t index1=0, index2=0, index3=0;
  Float_t weight12=0, weight13=0, weight23=0;
  Float_t weight12Err=0, weight13Err=0, weight23Err=0;
  Float_t weight12CC=0, weight13CC=0, weight23CC=0;
  Float_t weightTotal=0;//, weightTotalErr=0;
  Float_t qinv12MC=0, qinv13MC=0, qinv23MC=0; 
  Float_t Qsum1v1=0, Qsum2=0, Qsum3v1=0, Qsum1v2=0, Qsum3v2=0;
  Float_t Qsum1v1MC=0, Qsum2MC=0, Qsum3v1MC=0, Qsum1v2MC=0, Qsum3v2MC=0;
  //
  AliAODMCParticle *mcParticle1=0x0;
  AliAODMCParticle *mcParticle2=0x0;
  

  if(fPdensityPairCut){
    ////////////////////
    Int_t pairCountSE=0, pairCountME=0;
    Int_t normPairCount[2]={0};
    Int_t numOtherPairs1[2][fMultLimit];
    Int_t numOtherPairs2[2][fMultLimit];
    Bool_t exitCode=kFALSE;
    Int_t tempNormFillCount[2][2][2][10][5];


    // reset to defaults
    for(Int_t i=0; i<fMultLimit; i++) {
      fPairLocationSE[i]->Set(fMultLimit,fDefaultsInt);
      fPairLocationME[i]->Set(fMultLimit,fDefaultsInt);
           
      // Normalization Utilities
      fOtherPairLocation1[0][i]->Set(fMultLimit,fDefaultsInt);
      fOtherPairLocation1[1][i]->Set(fMultLimit,fDefaultsInt);
      fOtherPairLocation2[0][i]->Set(fMultLimit,fDefaultsInt);
      fOtherPairLocation2[1][i]->Set(fMultLimit,fDefaultsInt);
      fNormPairSwitch[0][i]->Set(fMultLimit,fDefaultsCharMult);
      fNormPairSwitch[1][i]->Set(fMultLimit,fDefaultsCharMult);
      fNormPairSwitch[2][i]->Set(fMultLimit,fDefaultsCharMult);
      numOtherPairs1[0][i]=0;
      numOtherPairs1[1][i]=0;
      numOtherPairs2[0][i]=0;
      numOtherPairs2[1][i]=0;
      
      // Track Merging/Splitting Utilities
      fPairSplitCut[0][i]->Set(fMultLimit,fDefaultsCharMult);// P11
      fPairSplitCut[1][i]->Set(fMultLimit,fDefaultsCharMult);// P12
      fPairSplitCut[2][i]->Set(fMultLimit,fDefaultsCharMult);// P13
      fPairSplitCut[3][i]->Set(fMultLimit,fDefaultsCharMult);// P23
    }

    // Reset the temp Normalization counters
    for(Int_t i=0; i<2; i++){// Charge1
      for(Int_t j=0; j<2; j++){// Charge2
	for(Int_t k=0; k<2; k++){// Charge3
	  for(Int_t l=0; l<10; l++){// FillIndex (species Combination)
	    for(Int_t m=0; m<5; m++){// Term (Cumulant term)
	      tempNormFillCount[i][j][k][l][m] = 0;
	    }
	  }
	}
      }
    }
  	      
 
    ///////////////////////////////////////////////////////
    // Start the pairing process
    // P11 pairing
    // 1st Particle
    
    for (Int_t i=0; i<myTracks; i++) {
         
      Int_t en2=0;
   
      // 2nd particle
      for (Int_t j=i+1; j<(fEvt+en2)->fNtracks; j++) {
	
	key1 = (fEvt)->fTracks[i].fKey;
	key2 = (fEvt+en2)->fTracks[j].fKey;
	Short_t fillIndex2 = FillIndex2part(key1+key2);
	Short_t qCutBin = SetQcutBin(fillIndex2);
	Short_t normBin = SetNormBin(fillIndex2);
	pVect1[0]=(fEvt)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	pVect1[1]=(fEvt)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	pVect1[2]=(fEvt)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	pVect1[3]=(fEvt)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	
	//
	
	qinv12 = GetQinv(fillIndex2, pVect1, pVect2);
	GetQosl(pVect1, pVect2, qout, qside, qlong);
	transK12 = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;


	if(fGenerateSignal){// Flatten the Q-dist to increase pair population at low-q (testing purposes only)
	  /*Float_t Qflattened = 0.005 + 0.2*gRandom->Rndm();
	  Float_t theta12 = PI*gRandom->Rndm();
	  Float_t phi12 = 2*PI*gRandom->Rndm();
	  pVect2Flat[1] = pVect1[1] + Qflattened*sin(theta12)*cos(phi12);
	  pVect2Flat[2] = pVect1[2] + Qflattened*sin(theta12)*sin(phi12);
	  pVect2Flat[3] = pVect1[3] + Qflattened*cos(theta12);
	  pVect2Flat[0] = sqrt(pow(pVect2Flat[1],2)+pow(pVect2Flat[2],2)+pow(pVect2Flat[3],2)+pow(fTrueMassPi,2));*/
	  //
	  pVect2Flat[0]=pVect2[0]; pVect2Flat[1]=pVect2[1]; pVect2Flat[2]=pVect2[2]; pVect2Flat[3]=pVect2[3]; 
	  //
	  qinv12Flat = GetQinv(fillIndex2, pVect1, pVect2Flat);
	  GetQosl(pVect1, pVect2Flat, qoutFlat, qsideFlat, qlongFlat);
	}
	
	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting/merging cut)
	
	
	//

	///////////////////////////////
	ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	SetFillBins2(fillIndex2, key1, key2, ch1, ch2, bin1, bin2);
	
	if(fMCcase && ch1==ch2 && fMbin==0){
	  for(Int_t rstep=0; rstep<10; rstep++){
	    Float_t coeff = (rstep)*0.2*(0.18/1.2);
	    // propagate through B field to r=1.2m
	    Float_t phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	    if(phi1 > 2*PI) phi1 -= 2*PI;
	    if(phi1 < 0) phi1 += 2*PI;
	    Float_t phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	    if(phi2 > 2*PI) phi2 -= 2*PI;
	    if(phi2 < 0) phi2 += 2*PI;
	    Float_t deltaphi = phi1 - phi2;
	    if(deltaphi > PI) deltaphi -= PI;
	    if(deltaphi < -PI) deltaphi += PI;
	    ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiNum"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
	  }
	}

	// Pair Splitting/Merging cut
	if(ch1 == ch2){
	  if(!AcceptPair((fEvt)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	    fPairSplitCut[0][i]->AddAt('1',j);
	    ((TH1F*)fOutputList->FindObject("fRejectedPairs"))->Fill(qinv12);
	    continue;
	  }
	}

	// HIJING tests 
	if(fMCcase && fillIndex2==0){
	  
	  // Check that label does not exceed stack size
	  if((fEvt)->fTracks[i].fLabel < (fEvt)->fMCarraySize && (fEvt+en2)->fTracks[j].fLabel < (fEvt+en2)->fMCarraySize){
	    pVect1MC[0]=sqrt(pow((fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
	    pVect2MC[0]=sqrt(pow((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPtot,2)+pow(fTrueMassPi,2));
	    pVect1MC[1]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPx; pVect2MC[1]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPx;
	    pVect1MC[2]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPy; pVect2MC[2]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPy;
	    pVect1MC[3]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPz; pVect2MC[3]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPz;
	    qinv12MC = GetQinv(fillIndex2, pVect1MC, pVect2MC);
	    GetQosl(pVect1MC, pVect2MC, qoutMC, qsideMC, qlongMC);
	    if(qinv12<0.1 && ch1==ch2) {
	      ((TProfile*)fOutputList->FindObject("fQsmearMean"))->Fill(1.,qinv12-qinv12MC); 
	      ((TProfile*)fOutputList->FindObject("fQsmearSq"))->Fill(1.,1000.*pow(qinv12-qinv12MC,2));
	      ((TH1D*)fOutputList->FindObject("fQDist"))->Fill(qinv12-qinv12MC);
	    }
	    
	    //if(transK12 <= 0.35) fEDbin=0;
	    //else fEDbin=1;

	    /*for(Int_t rIter=0; rIter<fRVALUES; rIter++){// 3fm to 8fm + 1 Therminator setting
	      for(Int_t myDampIt=0; myDampIt<kNDampValues; myDampIt++){
		Int_t denIndex = rIter*kNDampValues + myDampIt;
		Float_t WInput = MCWeight(ch1,ch2, rIter+kRmin, myDampIt, qinv12MC);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[0].fIdeal->Fill(denIndex, qinv12MC, WInput);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[1].fIdeal->Fill(denIndex, qinv12MC);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[0].fSmeared->Fill(denIndex, qinv12, WInput);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[1].fSmeared->Fill(denIndex, qinv12);
	      }
	      }*/
	    //fEDbin=0;

	    mcParticle1 = (AliAODMCParticle*)mcArray->At(abs((fEvt)->fTracks[i].fLabel));
	    mcParticle2 = (AliAODMCParticle*)mcArray->At(abs((fEvt+en2)->fTracks[j].fLabel));
	    
	    //HIJING resonance test
	    if(ch1 != ch2){
	      ((TH1F*)fOutputList->FindObject("fAllOSPairs"))->Fill(fMbin+1, qinv12);
	      if(abs(mcParticle1->GetPdgCode())==211 && abs(mcParticle2->GetPdgCode())==211){// Pions
		if(mcParticle1->GetMother() == mcParticle2->GetMother() && mcParticle1->GetMother() >=0){
		  ((TH1F*)fOutputList->FindObject("fResonanceOSPairs"))->Fill(fMbin+1, qinv12);
		}
	      }
	    }
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMCqinv->Fill(qinv12MC, MCWeight(ch1,ch2,10,10,qinv12MC));// was 4,5
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMCqinvQW->Fill(qinv12MC, qinv12MC*MCWeight(ch1,ch2,10,10,qinv12MC));// was 4,5
	    // pion purity	    
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fPIDpurityDen->Fill(transK12, qinv12);
	    if(abs(mcParticle1->GetPdgCode())==211 && abs(mcParticle2->GetPdgCode())==211){// Pions
	      Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fPIDpurityNum->Fill(transK12, qinv12);
	    }

	  }// label check 2
	}// MC case
	
	//////////////////////////////////////////
	// 2-particle term
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2->Fill(transK12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2QW->Fill(transK12, qinv12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fAvgP->Fill(transK12, qinv12, (fEvt)->fTracks[i].fMom);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fAvgP->Fill(transK12, qinv12, (fEvt+en2)->fTracks[j].fMom);
	
	// osl frame
	if(fillIndex2==0){
	  if((transK12 > 0.2) && (transK12 < 0.3)){  
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fExplicit2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fExplicit2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	  }
	  if((transK12 > 0.6) && (transK12 < 0.7)){  
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fExplicit2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fExplicit2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	  }
	}
	
	//////////////////////////////////////////
	if(fTabulatePairs){
	  if(fillIndex2==0 && bin1==bin2){
	    rapK12 = 0;
	    transKbin=-1; rapKbin=-1;
	    
	    for(Int_t kIt=0; kIt<fKbinsT; kIt++) {if(transK12 < (fKmiddleT[kIt] + fKstepT[kIt]/2.)) {transKbin = kIt; break;}} 
	    for(Int_t kIt=0; kIt<fKbinsY; kIt++) {if(rapK12 < (fKmiddleY[kIt] + fKstepY[kIt]/2.)) {rapKbin = kIt; break;}}
	    if((transKbin<0) || (rapKbin<0)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    if((transKbin>=fKbinsT) || (rapKbin>=fKbinsY)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    Float_t WInput = 1.0;
	    if(fGenerateSignal) {
	      WInput = MCWeight(ch1,ch2, fRMax, fFixedLambdaBinMomRes, qinv12Flat);
	      KT[transKbin].KY[rapKbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2ThreeD->Fill(fabs(qoutFlat), fabs(qsideFlat), fabs(qlongFlat), WInput);
	    }else KT[transKbin].KY[rapKbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    
	    continue;
	  }
	}
	

	// exit out of loop if there are too many pairs  
	if(pairCountSE >= kPairLimit) {exitCode=kTRUE; continue;}// Too many SE pairs
	if(exitCode) continue;

	//////////////////////////
	// Enforce the Qcut
	if(qinv12 <= fQcut[qCutBin]) {
	 
	  ///////////////////////////
	  // particle 1
	  (fEvt)->fPairsSE[pairCountSE].fP1[0] = (fEvt)->fTracks[i].fP[0];
	  (fEvt)->fPairsSE[pairCountSE].fP1[1] = (fEvt)->fTracks[i].fP[1];
	  (fEvt)->fPairsSE[pairCountSE].fP1[2] = (fEvt)->fTracks[i].fP[2];
	  (fEvt)->fPairsSE[pairCountSE].fE1 = (fEvt)->fTracks[i].fEaccepted;
	  (fEvt)->fPairsSE[pairCountSE].fCharge1 = (fEvt)->fTracks[i].fCharge;
	  (fEvt)->fPairsSE[pairCountSE].fIndex1 = i;
	  (fEvt)->fPairsSE[pairCountSE].fKey1 = key1;
	  (fEvt)->fPairsSE[pairCountSE].fLabel1 = (fEvt)->fTracks[i].fLabel;
	  if(fMCcase && ((fEvt)->fTracks[i].fLabel < (fEvt)->fMCarraySize)){
	    (fEvt)->fPairsSE[pairCountSE].fP1MC[0] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPx;
	    (fEvt)->fPairsSE[pairCountSE].fP1MC[1] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPy;
	    (fEvt)->fPairsSE[pairCountSE].fP1MC[2] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPz;
	  }
	  // particle 2
	  (fEvt)->fPairsSE[pairCountSE].fP2[0] = (fEvt+en2)->fTracks[j].fP[0];
	  (fEvt)->fPairsSE[pairCountSE].fP2[1] = (fEvt+en2)->fTracks[j].fP[1];
	  (fEvt)->fPairsSE[pairCountSE].fP2[2] = (fEvt+en2)->fTracks[j].fP[2];
	  (fEvt)->fPairsSE[pairCountSE].fE2 = (fEvt+en2)->fTracks[j].fEaccepted;
	  (fEvt)->fPairsSE[pairCountSE].fCharge2 = (fEvt+en2)->fTracks[j].fCharge;
	  (fEvt)->fPairsSE[pairCountSE].fIndex2 = j;
	  (fEvt)->fPairsSE[pairCountSE].fKey2 = key2;
	  (fEvt)->fPairsSE[pairCountSE].fLabel2 = (fEvt+en2)->fTracks[j].fLabel;
	  if(fMCcase && ((fEvt+en2)->fTracks[j].fLabel < (fEvt+en2)->fMCarraySize)){
	    (fEvt)->fPairsSE[pairCountSE].fP2MC[0] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPx;
	    (fEvt)->fPairsSE[pairCountSE].fP2MC[1] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPy;
	    (fEvt)->fPairsSE[pairCountSE].fP2MC[2] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPz;
	  }
	
	  (fEvt)->fPairsSE[pairCountSE].fQinv = qinv12;
	  
	  fPairLocationSE[i]->AddAt(pairCountSE,j);
	  
	  pairCountSE++;
	  
	}
	
	
	/////////////////////////////////////////////////////////
	// Normalization Region
	
	if((qinv12 >= fNormQcutLow[normBin]) && (qinv12 < fNormQcutHigh[normBin])){
	  // particle 1
	  fNormPairs[en2][normPairCount[en2]].fCharge1 = (fEvt)->fTracks[i].fCharge;
	  fNormPairs[en2][normPairCount[en2]].fIndex1 = i;
	  fNormPairs[en2][normPairCount[en2]].fKey1 = (fEvt)->fTracks[i].fKey;
	  // particle 2
	  fNormPairs[en2][normPairCount[en2]].fCharge2 = (fEvt+en2)->fTracks[j].fCharge;
	  fNormPairs[en2][normPairCount[en2]].fIndex2 = j;
	  fNormPairs[en2][normPairCount[en2]].fKey2 = (fEvt+en2)->fTracks[j].fKey;
	  
	  
	  //other past pairs with particle j
	  for(Int_t pastpair=0; pastpair<numOtherPairs2[0][j]; pastpair++){
	    Int_t locationOtherPair = fOtherPairLocation2[0][j]->At(pastpair);
	    if(locationOtherPair < 0) continue;// no pair there
	    Int_t indexOther1 = i;
	    Int_t indexOther2 = fNormPairs[0][ locationOtherPair ].fIndex1;
	    
	    // Both possible orderings of other indexes
	    if( (fNormPairSwitch[0][indexOther1]->At(indexOther2)=='1') || (fNormPairSwitch[0][indexOther2]->At(indexOther1)=='1')) {
	      
	      // 1 and 2 are from SE
	      ch3 = Int_t((fNormPairs[0][ locationOtherPair ].fCharge1 + 1)/2.);
	      key3 = fNormPairs[0][ locationOtherPair ].fKey1;
	      Short_t fillIndex3 = FillIndex3part(key1+key2+key3);
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 0, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	      
	      tempNormFillCount[bin1][bin2][bin3][fillIndex3][0]++;
	    }
	    
	  }// pastpair P11 loop
	  
	  
	  fNormPairSwitch[en2][i]->AddAt('1',j);	    
	  fOtherPairLocation1[en2][i]->AddAt(normPairCount[en2], numOtherPairs1[en2][i]);// location of otherpair with i as 1st particle
	  fOtherPairLocation2[en2][j]->AddAt(normPairCount[en2], numOtherPairs2[en2][j]);// location of otherpair with j as 2nd particle
	  
	  numOtherPairs1[en2][i]++;
	  numOtherPairs2[en2][j]++;
	  
	  
	  normPairCount[en2]++;
	  if(normPairCount[en2] >= kNormPairLimit) exitCode=kTRUE;
	  
	}// Norm Region
	
      }// j particle
    }// i particle
    
    
    
    //////////////////////////////////////////////
    // P12 pairing
    // 1st Particle
    for (Int_t i=0; i<myTracks; i++) {
         
      Int_t en2=1;

      // 2nd particle
      for (Int_t j=0; j<(fEvt+en2)->fNtracks; j++) {
	  	
	key1 = (fEvt)->fTracks[i].fKey;
	key2 = (fEvt+en2)->fTracks[j].fKey;
	Short_t fillIndex2 = FillIndex2part(key1+key2);
	Short_t qCutBin = SetQcutBin(fillIndex2);
	Short_t normBin = SetNormBin(fillIndex2);
	pVect1[0]=(fEvt)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	pVect1[1]=(fEvt)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	pVect1[2]=(fEvt)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	pVect1[3]=(fEvt)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	
	qinv12 = GetQinv(fillIndex2, pVect1, pVect2);
	GetQosl(pVect1, pVect2, qout, qside, qlong);
	transK12 = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;
	//if(transK12 <= 0.35) fEDbin=0;
	//else fEDbin=1;

	if(fGenerateSignal){// Flatten the Q-dist to increase pair population at low-q (testing purposes only)
	  /*Float_t Qflattened = 0.005 + 0.2*gRandom->Rndm();
	  Float_t theta12 = PI*gRandom->Rndm();
	  Float_t phi12 = 2*PI*gRandom->Rndm();
	  pVect2Flat[1] = pVect1[1] + Qflattened*sin(theta12)*cos(phi12);
	  pVect2Flat[2] = pVect1[2] + Qflattened*sin(theta12)*sin(phi12);
	  pVect2Flat[3] = pVect1[3] + Qflattened*cos(theta12);
	  pVect2Flat[0] = sqrt(pow(pVect2Flat[1],2)+pow(pVect2Flat[2],2)+pow(pVect2Flat[3],2)+pow(fTrueMassPi,2));*/
	  //
	  pVect2Flat[0]=pVect2[0]; pVect2Flat[1]=pVect2[1]; pVect2Flat[2]=pVect2[2]; pVect2Flat[3]=pVect2[3]; 
	  //
	  qinv12Flat = GetQinv(fillIndex2, pVect1, pVect2Flat);
	  GetQosl(pVect1, pVect2Flat, qoutFlat, qsideFlat, qlongFlat);
	}

	
	///////////////////////////////
	ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	SetFillBins2(fillIndex2, key1, key2, ch1, ch2, bin1, bin2);
	
	if(fMCcase){
	  if((fEvt)->fTracks[i].fLabel < (fEvt)->fMCarraySize && (fEvt+en2)->fTracks[j].fLabel < (fEvt+en2)->fMCarraySize){
	    pVect1MC[0]=sqrt(pow((fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
	    pVect2MC[0]=sqrt(pow((fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPtot,2)+pow(fTrueMassPi,2));
	    pVect1MC[1]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPx; pVect2MC[1]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPx;
	    pVect1MC[2]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPy; pVect2MC[2]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPy;
	    pVect1MC[3]=(fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPz; pVect2MC[3]=(fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPz;
	    qinv12MC = GetQinv(fillIndex2, pVect1MC, pVect2MC);
	    //
	    for(Int_t rIter=0; rIter<fRVALUES; rIter++){
	      for(Int_t myDampIt=0; myDampIt<kNDampValues; myDampIt++){
		Int_t denIndex = rIter*kNDampValues + myDampIt;
		Float_t WInput = MCWeight(ch1,ch2, rIter+kRmin, myDampIt, qinv12MC);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[0].fIdeal->Fill(denIndex, qinv12MC, WInput);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[1].fIdeal->Fill(denIndex, qinv12MC);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[0].fSmeared->Fill(denIndex, qinv12, WInput);
		Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[1].fSmeared->Fill(denIndex, qinv12);
	      }
	    }
	    if(qinv12 > fQcut[qCutBin]) continue;
	    
	    /////////////////////////////////////////////////////
	    if(!fTabulatePairs) {
	      // 3-particle MRC
	      Short_t fillIndex3 = 0;
	      key1=1; key2=1; key3=1;
	      Int_t en3 = 2;
	      
	      for (Int_t k=0; k<(fEvt+en3)->fNtracks; k++) {
		if((fEvt+en3)->fTracks[k].fLabel < (fEvt+en3)->fMCarraySize){
		  ch3 = Int_t(((fEvt+en3)->fTracks[k].fCharge + 1)/2.);
		  pVect3[0]=(fEvt+en3)->fTracks[k].fEaccepted;
		  pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
		  pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
		  pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
		  qinv13 = GetQinv(0, pVect1, pVect3);
		  qinv23 = GetQinv(0, pVect2, pVect3);
		  
		if(qinv13 > fQcut[qCutBin] || qinv23 > fQcut[qCutBin]) continue;

		
		pVect3MC[0]=sqrt(pow((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		pVect3MC[1]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPx;
		pVect3MC[2]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPy;
		pVect3MC[3]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPz;
		qinv13MC = GetQinv(0, pVect1MC, pVect3MC);
		qinv23MC = GetQinv(0, pVect2MC, pVect3MC);
		
		
		q3MC = sqrt(pow(qinv12MC,2)+pow(qinv13MC,2)+pow(qinv23MC,2));

		//
		// The below call to SetFillBins3 will work for all 3-particle terms since all are for fully mixed events. part is set to 1, but only matters for terms 2-4.
		Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
		SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);
	  
	
		for(Int_t jj=1; jj<=5; jj++){// term loop
		  
		  if(jj==2) {if(!fill2) continue;}//12
		  if(jj==3) {if(!fill3) continue;}//13
		  if(jj==4) {if(!fill4) continue;}//23
		
		  Float_t WInput=1.0;
		  Double_t K3=1.0;
		  ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, 1, jj, firstQ, secondQ, thirdQ);
		  ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, 1, jj, firstQMC, secondQMC, thirdQMC);

		  if(ch1==ch2 && ch1==ch3){// same charge
		    WInput = MCWeight3D(kTRUE, jj, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		    if(jj==1) {
		      K3 = FSICorrelationTherm2(+1,+1, firstQMC)*FSICorrelationTherm2(+1,+1, secondQMC)*FSICorrelationTherm2(+1,+1, thirdQMC);// GRS
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SC"))->Fill(q3MC, WInput);
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SCden"))->Fill(q3MC);
		    }else if(jj!=5){
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SC"))->Fill(q3MC, WInput);
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SCden"))->Fill(q3MC);
		    }
		  }else {// mixed charge
		    if(bin1==bin2) {
		      WInput = MCWeight3D(kFALSE, jj, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		      if(jj==1) K3 = FSICorrelationTherm2(+1,+1, firstQMC)*FSICorrelationTherm2(+1,-1, secondQMC)*FSICorrelationTherm2(+1,-1, thirdQMC);// GRS
		    }else {
		      if(jj==1 || jj==5) WInput = MCWeight3D(kFALSE, jj, fFixedLambdaBinMomRes, thirdQMC, secondQMC, firstQMC);// thirdQMC is ss
		      else WInput = MCWeight3D(kFALSE, 6-jj, fFixedLambdaBinMomRes, thirdQMC, secondQMC, firstQMC);
		    
		      if(jj==1) K3 = FSICorrelationTherm2(+1,+1, thirdQMC)*FSICorrelationTherm2(+1,-1, secondQMC)*FSICorrelationTherm2(+1,-1, firstQMC);// GRS
		    }
		    if(jj==1){
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1MC"))->Fill(q3MC, WInput);
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1MCden"))->Fill(q3MC);
		    }else{
		      if(bin1==bin2){
			if(jj==2){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MCden"))->Fill(q3MC);
			}else if(jj==3){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MCden"))->Fill(q3MC);
			}else if(jj==4){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MCden"))->Fill(q3MC);
			}else{}
		      }else{
			if(jj==2){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MCden"))->Fill(q3MC);
			}else if(jj==3){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MCden"))->Fill(q3MC);
			}else if(jj==4){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MCden"))->Fill(q3MC);
			}else{}
		      }
		      
		    }
		  }
		  
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fIdeal->Fill(firstQMC, secondQMC, thirdQMC, WInput);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fSmeared->Fill(firstQ, secondQ, thirdQ, WInput);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fQW12->Fill(firstQMC, secondQMC, thirdQMC, WInput*firstQMC);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fQW13->Fill(firstQMC, secondQMC, thirdQMC, WInput*secondQMC);
		  if(jj==1){
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fSumK3->Fill(firstQMC, secondQMC, thirdQMC, WInput/K3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fEnK3->Fill(firstQMC, secondQMC, thirdQMC, WInput);
		  }
		  
		  if(ch1==ch2 && ch1==ch3){
		    if(jj==1){
		      FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		      FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    }else if(jj==2) {
		      FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		      FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    }else if(jj==3){ 
		      FourVectProdTerms(pVect1, pVect3, pVect2, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		      FourVectProdTerms(pVect1MC, pVect3MC, pVect2MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    }else if(jj==4) {
		      FourVectProdTerms(pVect3, pVect1, pVect2, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		      FourVectProdTerms(pVect3MC, pVect1MC, pVect2MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    }else {
		      FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		      FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    }
		    
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsIdeal->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsSmeared->Fill(Qsum1v1, Qsum2, Qsum3v1, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsIdeal->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsSmeared->Fill(Qsum1v2, Qsum2, Qsum3v2, WInput);
		    if(jj==1){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1Q3W->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput*q3);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2Q3W->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput*q3);
		    }
		    //
		    if(qinv12MC > fQLowerCut && qinv13MC > fQLowerCut && qinv23MC > fQLowerCut){
		      // does not really matter if MC or real data triplets are used to average 1/K3...but better to use umsmeared values
		      if(jj==1){
			WInput = MCWeight3D(kTRUE, 1, 25, firstQMC, secondQMC, thirdQMC);// pure 3-pion (lambda=1)
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsSumK3->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput/K3);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsSumK3->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput/K3);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsEnK3->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsEnK3->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		      }if(jj>1 && jj<=4){
			Float_t InteractingQ=qinv12MC;
			Double_t K2 = FSICorrelationTherm2(+1,+1, InteractingQ);// K2 from Therminator source
			WInput = MCWeight3D(kTRUE, jj, 25, firstQMC, secondQMC, thirdQMC);// pure 2-pion (lambda=1)
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsSumK2->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput/K2);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsSumK2->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput/K2);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsEnK2->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsEnK2->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		      }
		    }

		  }// same charges
		
		}// jj
	      }// MCarray check
	      }// 3rd particle
	    }// TabulatePairs Check
	  }// MCarray check
	  
	  if(ch1==ch2 && fMbin==0){
	    for(Int_t rstep=0; rstep<10; rstep++){
	      Float_t coeff = (rstep)*0.2*(0.18/1.2);
	      // propagate through B field to r=1.2m
	      Float_t phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	      if(phi1 > 2*PI) phi1 -= 2*PI;
	      if(phi1 < 0) phi1 += 2*PI;
	      Float_t phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	      if(phi2 > 2*PI) phi2 -= 2*PI;
	      if(phi2 < 0) phi2 += 2*PI;
	      Float_t deltaphi = phi1 - phi2;
	      if(deltaphi > PI) deltaphi -= PI;
	      if(deltaphi < -PI) deltaphi += PI;
	      ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiDen"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
	    }
	  }
	  
	}// fMCcase
	
	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting cut)
	if(ch1 == ch2){
	  if(!AcceptPair((fEvt)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	    fPairSplitCut[1][i]->AddAt('1',j);
	    continue;
	  }
	}
	
	//////////////////////////////////////////
	// 2-particle term
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2->Fill(transK12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2QW->Fill(transK12, qinv12, qinv12);
	
	// osl frame
	if(fillIndex2==0){
	  if((transK12 > 0.2) && (transK12 < 0.3)){
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fExplicit2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[0].fExplicit2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	  }
	  if((transK12 > 0.6) && (transK12 < 0.7)){  
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fExplicit2OSL->Fill(fabs(qout), fabs(qside), fabs(qlong));
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].OSL_ktbin[1].fExplicit2OSLQW->Fill(fabs(qout), fabs(qside), fabs(qlong), qinv12);
	  }
	}
	//////////////////////////////////////////
	if(fTabulatePairs){
	  if(fillIndex2==0 && bin1==bin2){
	    rapK12 = 0;
	    transKbin=-1; rapKbin=-1;
	    
	    for(Int_t kIt=0; kIt<fKbinsT; kIt++) {if(transK12 < (fKmiddleT[kIt] + fKstepT[kIt]/2.)) {transKbin = kIt; break;}} 
	    for(Int_t kIt=0; kIt<fKbinsY; kIt++) {if(rapK12 < (fKmiddleY[kIt] + fKstepY[kIt]/2.)) {rapKbin = kIt; break;}}
	    if((transKbin<0) || (rapKbin<0)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    if((transKbin>=fKbinsT) || (rapKbin>=fKbinsY)) {cout<<"problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl; continue;}
	    
	    if(fGenerateSignal) KT[transKbin].KY[rapKbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2ThreeD->Fill(fabs(qoutFlat), fabs(qsideFlat), fabs(qlongFlat));
	    else KT[transKbin].KY[rapKbin].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2ThreeD->Fill(fabs(qout), fabs(qside), fabs(qlong));
	   
	    continue;
	  }
	}
	
	
	if(pairCountME >= 2*kPairLimit) {exitCode=kTRUE; continue;}// Too many SE pairs
	if(exitCode) continue;

	if(qinv12 <= fQcut[qCutBin]) {
	  ///////////////////////////
	  
	  // particle 1
	  (fEvt)->fPairsME[pairCountME].fP1[0] = (fEvt)->fTracks[i].fP[0];
	  (fEvt)->fPairsME[pairCountME].fP1[1] = (fEvt)->fTracks[i].fP[1];
	  (fEvt)->fPairsME[pairCountME].fP1[2] = (fEvt)->fTracks[i].fP[2];
	  (fEvt)->fPairsME[pairCountME].fE1 = (fEvt)->fTracks[i].fEaccepted;
	  (fEvt)->fPairsME[pairCountME].fCharge1 = (fEvt)->fTracks[i].fCharge;
	  (fEvt)->fPairsME[pairCountME].fIndex1 = i;
	  (fEvt)->fPairsME[pairCountME].fKey1 = key1;
	  (fEvt)->fPairsME[pairCountME].fLabel1 = (fEvt)->fTracks[i].fLabel;
	  if(fMCcase && ((fEvt)->fTracks[i].fLabel < (fEvt)->fMCarraySize)){
	    (fEvt)->fPairsME[pairCountME].fP1MC[0] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPx;
	    (fEvt)->fPairsME[pairCountME].fP1MC[1] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPy;
	    (fEvt)->fPairsME[pairCountME].fP1MC[2] = (fEvt)->fMCtracks[abs((fEvt)->fTracks[i].fLabel)].fPz;
	  }
	  // particle 2
	  (fEvt)->fPairsME[pairCountME].fP2[0] = (fEvt+en2)->fTracks[j].fP[0];
	  (fEvt)->fPairsME[pairCountME].fP2[1] = (fEvt+en2)->fTracks[j].fP[1];
	  (fEvt)->fPairsME[pairCountME].fP2[2] = (fEvt+en2)->fTracks[j].fP[2];
	  (fEvt)->fPairsME[pairCountME].fE2 = (fEvt+en2)->fTracks[j].fEaccepted;
	  (fEvt)->fPairsME[pairCountME].fCharge2 = (fEvt+en2)->fTracks[j].fCharge;
	  (fEvt)->fPairsME[pairCountME].fIndex2 = j;
	  (fEvt)->fPairsME[pairCountME].fKey2 = key2;
	  (fEvt)->fPairsME[pairCountME].fLabel2 = (fEvt+en2)->fTracks[j].fLabel;
	  if(fMCcase && ((fEvt+en2)->fTracks[j].fLabel < (fEvt+en2)->fMCarraySize)){
	    (fEvt)->fPairsME[pairCountME].fP2MC[0] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPx;
	    (fEvt)->fPairsME[pairCountME].fP2MC[1] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPy;
	    (fEvt)->fPairsME[pairCountME].fP2MC[2] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[j].fLabel)].fPz;
	  }
	  
	  (fEvt)->fPairsME[pairCountME].fQinv = qinv12;
	  
	  fPairLocationME[i]->AddAt(Int_t(pairCountME),j);
	  
	  pairCountME++;
	  
	}
	
	if((qinv12 >= fNormQcutLow[normBin]) && (qinv12 < fNormQcutHigh[normBin])){
	  // particle 1
	  fNormPairs[en2][normPairCount[en2]].fCharge1 = (fEvt)->fTracks[i].fCharge;
	  fNormPairs[en2][normPairCount[en2]].fIndex1 = i;
	  fNormPairs[en2][normPairCount[en2]].fKey1 = (fEvt)->fTracks[i].fKey;
	  // particle 2
	  fNormPairs[en2][normPairCount[en2]].fCharge2 = (fEvt+en2)->fTracks[j].fCharge;
	  fNormPairs[en2][normPairCount[en2]].fIndex2 = j;
	  fNormPairs[en2][normPairCount[en2]].fKey2 = (fEvt+en2)->fTracks[j].fKey;
	  
	  //other past pairs in P11 with particle i
	  for(Int_t pastpairP11=0; pastpairP11<numOtherPairs2[0][i]; pastpairP11++){// past pair in P11 with i as 1st and 2nd particle
	    Int_t locationOtherPairP11 = fOtherPairLocation2[0][i]->At(pastpairP11);// i is 2nd particle
	    if(locationOtherPairP11 < 0) continue;// no pair there
	    Int_t indexOther1P11 = fNormPairs[0][ locationOtherPairP11 ].fIndex1; 
	    	    
	    //Check other past pairs in P12
	    if( (fNormPairSwitch[1][indexOther1P11]->At(j)=='0')) continue;
	    
	    // 1 and 3 are from SE
	    ch3 = Int_t((fNormPairs[0][ locationOtherPairP11 ].fCharge1 + 1)/2.);// charge of second particle in P11
	    key3 = fNormPairs[0][ locationOtherPairP11 ].fKey1;
	    Short_t fillIndex3 = FillIndex3part(key1+key2+key3);
	    Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
	    SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 2, bin1, bin2, bin3, fill2, fill3, fill4);
	    
	 	    
	    if(fill2) tempNormFillCount[bin1][bin2][bin3][fillIndex3][1]++;
	    if(fill3) tempNormFillCount[bin1][bin2][bin3][fillIndex3][2]++;
	    if(fill4) tempNormFillCount[bin1][bin2][bin3][fillIndex3][3]++;
	    

	  }// P11 loop
	  
	  
	  fNormPairSwitch[en2][i]->AddAt('1',j);	    
	  fOtherPairLocation1[en2][i]->AddAt(normPairCount[en2], numOtherPairs1[en2][i]);// location of otherpair with i as 1st particle
	  fOtherPairLocation2[en2][j]->AddAt(normPairCount[en2], numOtherPairs2[en2][j]);// location of otherpair with j as 2nd particle
	  
	  numOtherPairs1[en2][i]++;
	  numOtherPairs2[en2][j]++;
	  
	  normPairCount[en2]++;
	  if(normPairCount[en2] >= kNormPairLimit) exitCode=kTRUE;

	}// Norm Region
	

      }
    }
    
 
    ///////////////////////////////////////
    // P13 pairing (just for Norm counting of term5)
    for (Int_t i=0; i<myTracks; i++) {
      
      // exit out of loop if there are too many pairs
      // dont bother with this loop if exitCode is set.
      if(exitCode) break;
      
      // 2nd particle
      Int_t en2=2;
      
      for (Int_t j=0; j<(fEvt+en2)->fNtracks; j++) {
	
	key1 = (fEvt)->fTracks[i].fKey;
	key2 = (fEvt+en2)->fTracks[j].fKey;
	Short_t fillIndex2 = FillIndex2part(key1+key2);
	Short_t normBin = SetNormBin(fillIndex2);
	pVect1[0]=(fEvt)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	pVect1[1]=(fEvt)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	pVect1[2]=(fEvt)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	pVect1[3]=(fEvt)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];

	qinv12 = GetQinv(fillIndex2, pVect1, pVect2);
	
	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting cut)
	
	ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	
	if(ch1 == ch2){
	  if(!AcceptPair((fEvt)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	    fPairSplitCut[2][i]->AddAt('1',j);
	    continue;
	  }
	}
	
	/////////////////////////////////////////////////////////
	// Normalization Region
	
	if((qinv12 >= fNormQcutLow[normBin]) && (qinv12 < fNormQcutHigh[normBin])){
	
	  fNormPairSwitch[en2][i]->AddAt('1',j);	    
	
	}// Norm Region
      }
    }


   
    ///////////////////////////////////////
    // P23 pairing (just for Norm counting of term5)
    Int_t en1=1;
    for (Int_t i=0; i<(fEvt+en1)->fNtracks; i++) {
      
      // exit out of loop if there are too many pairs
      // dont bother with this loop if exitCode is set.
      if(exitCode) break;
      
      // 2nd event
      Int_t en2=2;
      // 2nd particle
      for (Int_t j=0; j<(fEvt+en2)->fNtracks; j++) {
	
	if(exitCode) break;

	key1 = (fEvt+en1)->fTracks[i].fKey;
	key2 = (fEvt+en2)->fTracks[j].fKey;
	Short_t fillIndex2 = FillIndex2part(key1+key2);
	Short_t normBin = SetNormBin(fillIndex2);
	pVect1[0]=(fEvt+en1)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	pVect1[1]=(fEvt+en1)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	pVect1[2]=(fEvt+en1)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	pVect1[3]=(fEvt+en1)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];

	qinv12 = GetQinv(fillIndex2, pVect1, pVect2);

	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting cut)
	
	///////////////////////////////
	ch1 = Int_t(((fEvt+en1)->fTracks[i].fCharge + 1)/2.);
	ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	
	if(ch1 == ch2){
	  if(!AcceptPair((fEvt+en1)->fTracks[i], (fEvt+en2)->fTracks[j])) {
	    fPairSplitCut[3][i]->AddAt('1',j);
	    continue;
	  }
	}

	if((qinv12 < fNormQcutLow[normBin]) || (qinv12 >= fNormQcutHigh[normBin])) continue;
	
	Int_t index1P23 = i;
	Int_t index2P23 = j;
	
	for(Int_t pastpairP12=0; pastpairP12<numOtherPairs2[1][index1P23]; pastpairP12++){// loop in P12 with i as 2nd particle
	  Int_t locationOtherPairP12 = fOtherPairLocation2[1][index1P23]->At(pastpairP12);
	  if(locationOtherPairP12 < 0) continue; // no pair there
	  Int_t index1P12 = fNormPairs[1][ locationOtherPairP12 ].fIndex1;
	  
	 	  
	  //Check other past pair status in P13
	  if( (fNormPairSwitch[2][index1P12]->At(index2P23)=='0')) continue;
	  
	  // all from different event
	  ch3 = Int_t((fNormPairs[1][ locationOtherPairP12 ].fCharge1 + 1)/2.);// charge of first particle in P12
	  key3 = fNormPairs[1][ locationOtherPairP12 ].fKey1;
	  Short_t fillIndex3 = FillIndex3part(key1+key2+key3);
	  SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 3, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	  
	  tempNormFillCount[bin1][bin2][bin3][fillIndex3][4]++;
	}
      }
    }
    
    
  
    
    ///////////////////////////////////////////////////  
    // Do not use pairs from events with too many pairs
    if(exitCode) {
      cout<<"SE or ME or Norm PairCount too large.  Discarding all pairs and skipping event"<<endl;
      (fEvt)->fNpairsSE = 0;
      (fEvt)->fNpairsME = 0;
      ((TH1F*)fOutputList->FindObject("fRejectedEvents"))->Fill(fMbin+1);
      return;// Skip event
    }else{
      (fEvt)->fNpairsSE = pairCountSE;
      (fEvt)->fNpairsME = pairCountME;  
      ((TH1F*)fOutputList->FindObject("fEvents2"))->Fill(fMbin+1);
    }
    ///////////////////////////////////////////////////


    //cout<<"pairCountSE = "<<pairCountSE<<"   pairCountME = "<<pairCountME<<endl;
    //cout<<"Start Main analysis"<<endl;
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //
    //
    // Start the Main Correlation Analysis
    //
    //
    ///////////////////////////////////////////////////////////////////////
    
    
    
    /////////////////////////////////////////////////////////    
    // Skip 3-particle part if Tabulate6DPairs is set to true
    if(fTabulatePairs) return;
    /////////////////////////////////////////////////////////

    // Set the Normalization counters
    for(Int_t termN=0; termN<5; termN++){
      
      if(termN==0){
	if((fEvt)->fNtracks ==0) continue;
      }else if(termN<4){
	if((fEvt)->fNtracks ==0) continue;
	if((fEvt+1)->fNtracks ==0) continue;
      }else {
	if((fEvt)->fNtracks ==0) continue;
	if((fEvt+1)->fNtracks ==0) continue;
	if((fEvt+2)->fNtracks ==0) continue;
      }
     
      for(Int_t sc=0; sc<kSCLimit3; sc++){
	
	for(Int_t c1=0; c1<2; c1++){
	  for(Int_t c2=0; c2<2; c2++){
	    for(Int_t c3=0; c3<2; c3++){
	      
	      if(sc==0 || sc==6 || sc==9){// Identical species
		if( (c1+c2+c3)==1) {if(c1!=0 || c2!=0 || c3!=1) continue;}
		if( (c1+c2+c3)==2) {if(c1!=0) continue;}
	      }else if(sc!=5){
		if( (c1+c2)==1) {if(c1!=0) continue;}
	      }else {}// do nothing for pi-k-p case
	      Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[fMbin].EDB[fEDbin].ThreePT[termN].fNorm3->Fill(0.,tempNormFillCount[c1][c2][c3][sc][termN]);
	    }
	  }
	}
      }
    }
    
    
    
    /////////////////////////////////////////////
    // Calculate Pair-Cut Correlations
    for(Int_t en1case=0; en1case<2; en1case++){// limit at 2 (normal)
      
      Int_t nump1=0;
      if(en1case==0) nump1 = (fEvt)->fNpairsSE;
      if(en1case==1) nump1 = (fEvt)->fNpairsME;
     
      // 1st pair
      for(Int_t p1=0; p1<nump1; p1++){
	
	if(en1case==0){
	  ch1 = Int_t(((fEvt)->fPairsSE[p1].fCharge1 + 1)/2.);
	  ch2 = Int_t(((fEvt)->fPairsSE[p1].fCharge2 + 1)/2.);
	  pVect1[0] = (fEvt)->fPairsSE[p1].fE1; pVect2[0] = (fEvt)->fPairsSE[p1].fE2;
	  pVect1[1] = (fEvt)->fPairsSE[p1].fP1[0]; pVect2[1] = (fEvt)->fPairsSE[p1].fP2[0]; 
	  pVect1[2] = (fEvt)->fPairsSE[p1].fP1[1]; pVect2[2] = (fEvt)->fPairsSE[p1].fP2[1];
	  pVect1[3] = (fEvt)->fPairsSE[p1].fP1[2]; pVect2[3] = (fEvt)->fPairsSE[p1].fP2[2];
	  index1 = (fEvt)->fPairsSE[p1].fIndex1; index2 = (fEvt)->fPairsSE[p1].fIndex2;
	  key1 = (fEvt)->fPairsSE[p1].fKey1; key2 = (fEvt)->fPairsSE[p1].fKey2;
	  /*pVect1MC[1] = (fEvt)->fPairsSE[p1].fP1MC[0]; pVect2MC[1] = (fEvt)->fPairsSE[p1].fP2MC[0]; 
	  pVect1MC[2] = (fEvt)->fPairsSE[p1].fP1MC[1]; pVect2MC[2] = (fEvt)->fPairsSE[p1].fP2MC[1];
	  pVect1MC[3] = (fEvt)->fPairsSE[p1].fP1MC[2]; pVect2MC[3] = (fEvt)->fPairsSE[p1].fP2MC[2];
	  pVect1MC[0] = sqrt(pow(pVect1MC[1],2)+pow(pVect1MC[2],2)+pow(pVect1MC[3],2)+pow(fTrueMassPi,2));
	  pVect2MC[0] = sqrt(pow(pVect2MC[1],2)+pow(pVect2MC[2],2)+pow(pVect2MC[3],2)+pow(fTrueMassPi,2));
	  */
	  qinv12 = (fEvt)->fPairsSE[p1].fQinv;
	}
	if(en1case==1){
	  ch1 = Int_t(((fEvt)->fPairsME[p1].fCharge1 + 1)/2.);
	  ch2 = Int_t(((fEvt)->fPairsME[p1].fCharge2 + 1)/2.);
	  pVect1[0] = (fEvt)->fPairsME[p1].fE1; pVect2[0] = (fEvt)->fPairsME[p1].fE2; 
	  pVect1[1] = (fEvt)->fPairsME[p1].fP1[0]; pVect2[1] = (fEvt)->fPairsME[p1].fP2[0]; 
	  pVect1[2] = (fEvt)->fPairsME[p1].fP1[1]; pVect2[2] = (fEvt)->fPairsME[p1].fP2[1];
	  pVect1[3] = (fEvt)->fPairsME[p1].fP1[2]; pVect2[3] = (fEvt)->fPairsME[p1].fP2[2];
	  index1 = (fEvt)->fPairsME[p1].fIndex1; index2 = (fEvt)->fPairsME[p1].fIndex2;
	  key1 = (fEvt)->fPairsME[p1].fKey1; key2 = (fEvt)->fPairsME[p1].fKey2;
	  /*pVect1MC[1] = (fEvt)->fPairsME[p1].fP1MC[0]; pVect2MC[1] = (fEvt)->fPairsME[p1].fP2MC[0]; 
	  pVect1MC[2] = (fEvt)->fPairsME[p1].fP1MC[1]; pVect2MC[2] = (fEvt)->fPairsME[p1].fP2MC[1];
	  pVect1MC[3] = (fEvt)->fPairsME[p1].fP1MC[2]; pVect2MC[3] = (fEvt)->fPairsME[p1].fP2MC[2];
	  pVect1MC[0] = sqrt(pow(pVect1MC[1],2)+pow(pVect1MC[2],2)+pow(pVect1MC[3],2)+pow(fTrueMassPi,2));
	  pVect2MC[0] = sqrt(pow(pVect2MC[1],2)+pow(pVect2MC[2],2)+pow(pVect2MC[3],2)+pow(fTrueMassPi,2));
	  */
	  qinv12 = (fEvt)->fPairsME[p1].fQinv;
	}

	/*if(fGenerateSignal){
	  Bool_t goodFlattenedPair=kFALSE;
	  while(!goodFlattenedPair){
	    Float_t Qflattened = fQLowerCut + (fQcut[0]-fQLowerCut)*gRandom->Rndm();
	    Float_t theta12 = PI*gRandom->Rndm();
	    Float_t phi12 = 2*PI*gRandom->Rndm();
	    pVect2Flat[1] = pVect1[1] + Qflattened*sin(theta12)*cos(phi12);
	    pVect2Flat[2] = pVect1[2] + Qflattened*sin(theta12)*sin(phi12);
	    pVect2Flat[3] = pVect1[3] + Qflattened*cos(theta12);
	    pVect2Flat[0] = sqrt(pow(pVect2Flat[1],2)+pow(pVect2Flat[2],2)+pow(pVect2Flat[3],2)+pow(fTrueMassPi,2));
	    //
	    //pVect2Flat[0]=pVect2[0]; pVect2Flat[1]=pVect2[1]; pVect2Flat[2]=pVect2[2]; pVect2Flat[3]=pVect2[3];
	    //
	    qinv12 = GetQinv(0, pVect1, pVect2Flat);
	    if(qinv12 < fQcut[0] && qinv12>fQLowerCut) goodFlattenedPair=kTRUE;
	  }
	  }*/
	
	// en2 buffer
	for(Int_t en2=0; en2<3; en2++){
	  //////////////////////////////////////

	  Bool_t skipcase=kTRUE;
	  Short_t config=-1, part=-1;
	  if(en1case==0 && en2==0) {skipcase=kFALSE; config=1; part=0;}// P11T1
	  if(en1case==0 && en2==1) {skipcase=kFALSE; config=2; part=1;}// P11T2
	  if(en1case==1 && en2==0) {skipcase=kFALSE; config=2; part=2;}// P12T1
	  if(en1case==1 && en2==2) {skipcase=kFALSE; config=3; part=3;}// P12T3
	  	 
	  if(skipcase) continue;
	
	  
	  // 3-particle terms
	  // 3rd particle
	  for(Int_t k=0; k<(fEvt+en2)->fNtracks; k++){
	    index3 = k;
	    

	    // remove auto-correlations and duplicate triplets
	    if(config==1){
	      if( index1 == index3) continue;
	      if( index2 == index3) continue;
	      if(fPairSplitCut[0][index1]->At(index2)=='1') continue;// Track splitting/merging
	     
	      // skip the switched off triplets
	      if(fTripletSkip1[fPairLocationSE[index1]->At(index2)]->At(index3)=='1') {
		fTripletSkip1[fPairLocationSE[index1]->At(index2)]->AddAt('0',index3);// Reset
		continue;
	      }
	      ///////////////////////////////
	      // Turn off 1st possible degenerate triplet
	      if(index1 < index3){// verify correct id ordering ( index1 < k )
		if(fPairLocationSE[index1]->At(index3) >= 0){
		  fTripletSkip1[fPairLocationSE[index1]->At(index3)]->AddAt('1',index2);
		}
		if(fPairSplitCut[0][index1]->At(index3)=='1') continue;// Track splitting/merging
	      }else {// or k < index1
		if(fPairLocationSE[index3]->At(index1) >= 0){
		  fTripletSkip1[fPairLocationSE[index3]->At(index1)]->AddAt('1',index2);
		}
		if(fPairSplitCut[0][index3]->At(index1)=='1') continue;// Track splitting/merging
	      }
	      // turn off 2nd possible degenerate triplet
	      if(index2 < index3){// verify correct id ordering (index2 < k)
		if(fPairLocationSE[index2]->At(index3) >= 0){
		  fTripletSkip1[fPairLocationSE[index2]->At(index3)]->AddAt('1',index1);
		}
		if(fPairSplitCut[0][index2]->At(index3)=='1') continue;// Track splitting/merging
	      }else {// or k < index2
		if(fPairLocationSE[index3]->At(index2) >= 0){
		  fTripletSkip1[fPairLocationSE[index3]->At(index2)]->AddAt('1',index1);
		}
		if(fPairSplitCut[0][index3]->At(index2)=='1') continue;// Track splitting/merging
	      }

	    }// end config 1
	    
	    if(config==2 && part==1){// SE pair and third particle from next event. P11T2
	      ///////////////////////////////
	      // Turn off 1st possible degenerate triplet
	      if(fPairLocationME[index1]->At(index3) >= 0){
		fTripletSkip2[fPairLocationME[index1]->At(index3)]->AddAt('1',index2);
	      }
	      
	      // turn off 2nd possible degenerate triplet
	      if(fPairLocationME[index2]->At(index3) >= 0){
		fTripletSkip2[fPairLocationME[index2]->At(index3)]->AddAt('1',index1);
	      }
	      
	      if(fPairSplitCut[0][index1]->At(index2)=='1') continue;// Track splitting/merging
	      if(fPairSplitCut[1][index1]->At(index3)=='1') continue;// Track splitting/merging
	      if(fPairSplitCut[1][index2]->At(index3)=='1') continue;// Track splitting/merging
	    }// end config 2 part 1

	    if(config==2 && part==2){// P12T1
	      if( index1 == index3) continue;
	      
	      // skip the switched off triplets
	      if(fTripletSkip2[fPairLocationME[index1]->At(index2)]->At(index3)=='1') {
		fTripletSkip2[fPairLocationME[index1]->At(index2)]->AddAt('0',index3);// Reset
		continue;
	      }
	      // turn off another possible degenerate
	      if(fPairLocationME[index3]->At(index2) >= 0){
		fTripletSkip2[fPairLocationME[index3]->At(index2)]->AddAt('1',index1);
	      }// end config 2 part 2

	      if(fPairSplitCut[1][index1]->At(index2)=='1') continue;// Track splitting/merging
	      if(index1 < index3) {if(fPairSplitCut[0][index1]->At(index3)=='1') continue;}// Track splitting/merging
	      else {if(fPairSplitCut[0][index3]->At(index1)=='1') continue;}// Track splitting/merging
	      if(fPairSplitCut[1][index3]->At(index2)=='1') continue;// Track splitting/merging
	    }
	    if(config==3){// P12T3
	      if(fPairSplitCut[1][index1]->At(index2)=='1') continue;// Track splitting/merging
	      if(fPairSplitCut[2][index1]->At(index3)=='1') continue;// Track splitting/merging
	      if(fPairSplitCut[3][index2]->At(index3)=='1') continue;// Track splitting/merging
	    }// end config 3
	    
	    

	    ch3 = Int_t(((fEvt+en2)->fTracks[k].fCharge + 1)/2.);
	    key3 = (fEvt+en2)->fTracks[k].fKey;
	    Short_t fillIndex3 = FillIndex3part(key1+key2+key3);
	    Short_t fillIndex13 = FillIndex2part(key1+key3);
	    Short_t fillIndex23 = FillIndex2part(key2+key3);
	    Short_t qCutBin13 = SetQcutBin(fillIndex13);
	    Short_t qCutBin23 = SetQcutBin(fillIndex23);
	    pVect3[0] = (fEvt+en2)->fTracks[k].fEaccepted;
	    pVect3[1] = (fEvt+en2)->fTracks[k].fP[0];
	    pVect3[2] = (fEvt+en2)->fTracks[k].fP[1];
	    pVect3[3] = (fEvt+en2)->fTracks[k].fP[2];
	    qinv13 = GetQinv(fillIndex13, pVect1, pVect3);
	    qinv23 = GetQinv(fillIndex23, pVect2, pVect3);

	    if(qinv13 < fQLowerCut) continue;
	    if(qinv23 < fQLowerCut) continue;
	    if(qinv13 > fQcut[qCutBin13]) continue;
	    if(qinv23 > fQcut[qCutBin23]) continue;

	    /*if(fGenerateSignal){
	      Bool_t goodFlattenedTriplet=kFALSE;
	      while(!goodFlattenedTriplet){
		Float_t Qflattened = fQLowerCut + (fQcut[0]-fQLowerCut)*gRandom->Rndm();
		Float_t theta13 = PI*gRandom->Rndm();
		Float_t phi13 = 2*PI*gRandom->Rndm();
		pVect3Flat[1] = pVect1[1] + Qflattened*sin(theta13)*cos(phi13);
		pVect3Flat[2] = pVect1[2] + Qflattened*sin(theta13)*sin(phi13);
		pVect3Flat[3] = pVect1[3] + Qflattened*cos(theta13);
		pVect3Flat[0] = sqrt(pow(pVect3Flat[1],2)+pow(pVect3Flat[2],2)+pow(pVect3Flat[3],2)+pow(fTrueMassPi,2));
		//
		pVect3Flat[0]=pVect3[0]; pVect3Flat[1]=pVect3[1]; pVect3Flat[2]=pVect3[2]; pVect3Flat[3]=pVect3[3];
		//
		qinv13 = GetQinv(0, pVect1, pVect3Flat);
		qinv23 = GetQinv(0, pVect2Flat, pVect3Flat);
		if(qinv13 < fQcut[qCutBin13] && qinv23 < fQcut[qCutBin23]) {
		  if(qinv13>fQLowerCut && qinv23>fQLowerCut) goodFlattenedTriplet=kTRUE;
		}
	      }
	      }*/
	    
	    
	    /*if(fMCcase){
	      pVect3MC[1] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPx;
	      pVect3MC[2] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPy;
	      pVect3MC[3] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPz;
	      pVect3MC[0] = sqrt(pow(pVect3MC[1],2)+pow(pVect3MC[2],2)+pow(pVect3MC[3],2)+pow(fTrueMassPi,2));
	      qinv12MC = GetQinv(0, pVect1MC, pVect2MC);
	      qinv13MC = GetQinv(0, pVect1MC, pVect3MC);
	      qinv23MC = GetQinv(0, pVect2MC, pVect3MC);
	      q3MC = sqrt(pow(qinv12MC,2) + pow(qinv13MC,2) + pow(qinv23MC,2));
	      }*/
	    
	  
	    // if all three pair cuts are the same then the case (config=2 && term=2) never reaches here.
	    
	    q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));
	    transK3 = sqrt( pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
	    //if(transK3<0.35) fEDbin=0;
	    //else fEDbin=1;
	    firstQ=0; secondQ=0; thirdQ=0;
	    //firstQMC=0; secondQMC=0; thirdQMC=0;
	    
	    //
	    
	    //	    
	    if(config==1) {// 123
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 0, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	      
	      if(fillIndex3 <= 2){
		ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, 0, 1, firstQ, secondQ, thirdQ);
		//if(fillIndex3==0 && fMCcase) ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, 0, 1, firstQMC, secondQMC, thirdQMC);
		Float_t WInput = 1.0;
		if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, 1, fFixedLambdaBinMomRes, firstQ, secondQ, thirdQ);
		////
		
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fTerms3->Fill(firstQ, secondQ, thirdQ, WInput);
		////
		//
		if(fillIndex3==0 && ch1==ch2 && ch1==ch3 && fMCcase==kFALSE){
		  FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1Terms->Fill(Qsum1v1, Qsum2, Qsum3v1);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2Terms->Fill(Qsum1v2, Qsum2, Qsum3v2);
		  ((TH3D*)fOutputList->FindObject("fKt3DistTerm1"))->Fill(fMbin+1, transK3, q3);
		}		
		//////////////////////////////////////
		// Momentum resolution and <K3> calculation
		/*if(fillIndex3==0 && fMCcase){
		 
		  WInput = 1.0;
		  Double_t K3=1.0;
		  if(ch1==ch2 && ch1==ch3){// same charge
		    WInput = MCWeight3D(kTRUE, 1, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		    //K3 = FSICorrelationOmega0(kTRUE, firstQMC, secondQMC, thirdQMC);// Omega0 method
		    K3 = FSICorrelationTherm2(+1,+1, firstQMC)*FSICorrelationTherm2(+1,+1, secondQMC)*FSICorrelationTherm2(+1,+1, thirdQMC);// GRS
		    ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SC"))->Fill(q3MC, WInput);
		    ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SCden"))->Fill(q3MC);
		  }else {// mixed charge
		    if(bin1==bin2) {
		      WInput = MCWeight3D(kFALSE, 1, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		      //K3 = FSICorrelationOmega0(kFALSE, firstQMC, secondQMC, thirdQMC);// Omega0 method
		      K3 = FSICorrelationTherm2(+1,+1, firstQMC)*FSICorrelationTherm2(+1,-1, secondQMC)*FSICorrelationTherm2(+1,-1, thirdQMC);// GRS
		    }else {
		      WInput = MCWeight3D(kFALSE, 1, fFixedLambdaBinMomRes, thirdQMC, secondQMC, firstQMC);// thirdQMC is ss 
		      //K3 = FSICorrelationOmega0(kFALSE, thirdQMC, secondQMC, firstQMC);// Omega0 method
		      K3 = FSICorrelationTherm2(+1,+1, thirdQMC)*FSICorrelationTherm2(+1,-1, secondQMC)*FSICorrelationTherm2(+1,-1, firstQMC);// GRS
		    }
		    ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1MC"))->Fill(q3MC, WInput);
		    ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1MCden"))->Fill(q3MC);
		  }
		  
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fIdeal->Fill(firstQMC, secondQMC, thirdQMC, WInput);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fSmeared->Fill(firstQ, secondQ, thirdQ, WInput);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fQW12->Fill(firstQMC, secondQMC, thirdQMC, WInput*firstQMC);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fQW13->Fill(firstQMC, secondQMC, thirdQMC, WInput*secondQMC);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fSumK3->Fill(firstQMC, secondQMC, thirdQMC, WInput/K3);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fEnK3->Fill(firstQMC, secondQMC, thirdQMC, WInput);
		  if(ch1==ch2 && ch1==ch3){
		    FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		    FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1TermsIdeal->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1TermsSmeared->Fill(Qsum1v1, Qsum2, Qsum3v1, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2TermsIdeal->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2TermsSmeared->Fill(Qsum1v2, Qsum2, Qsum3v2, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1Q3W->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput*q3);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2Q3W->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput*q3);
		    //
		    if(qinv12MC > fQLowerCut && qinv13MC > fQLowerCut && qinv23MC > fQLowerCut){
		      // does not really matter if MC or real data triplets are used to average 1/K3...but better to use umsmeared values
		      WInput = MCWeight3D(kTRUE, 1, 25, firstQMC, secondQMC, thirdQMC);// pure 3-pion (lambda=1)
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1TermsSumK3->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput/K3);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2TermsSumK3->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput/K3);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd1TermsEnK3->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].f4VectProd2TermsEnK3->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		    }
		  }
		}// fMCcase
		*/
	      }
	      
	    }else if(config==2){// 12, 13, 23
	      
	      Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, part, bin1, bin2, bin3, fill2, fill3, fill4);
	  
	      // loop over terms 2-4
	      for(Int_t jj=2; jj<5; jj++){
		if(jj==2) {if(!fill2) continue;}//12
		if(jj==3) {if(!fill3) continue;}//13
		if(jj==4) {if(!fill4) continue;}//23
	
		if(fillIndex3 <= 2){
		  ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, part, jj, firstQ, secondQ, thirdQ);
		  //if(fillIndex3==0 && fMCcase) ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, part, jj, firstQMC, secondQMC, thirdQMC);
		  Float_t WInput = 1.0;
		  if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, jj, fFixedLambdaBinMomRes, firstQ, secondQ, thirdQ);
		  ////
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fTerms3->Fill(firstQ, secondQ, thirdQ, WInput);
		  ////
		  if(fillIndex3==0 && ch1==ch2 && ch1==ch3){
		    if(part==1){// P11T2
		      if(jj==2) {
			FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }if(jj==3){ 
			FourVectProdTerms(pVect1, pVect3, pVect2, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect1MC, pVect3MC, pVect2MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }if(jj==4) {
			FourVectProdTerms(pVect3, pVect1, pVect2, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect3MC, pVect1MC, pVect2MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }		    
		    }else{// P12T1
		      if(jj==2) {
			FourVectProdTerms(pVect1, pVect3, pVect2, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect1MC, pVect3MC, pVect2MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }if(jj==3) {
			FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }if(jj==4) {
			FourVectProdTerms(pVect2, pVect1, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
			//if(fMCcase) FourVectProdTerms(pVect2MC, pVect1MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);// 4-vector product sums
		      }		    
		    }
		    if(!fMCcase){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1Terms->Fill(Qsum1v1, Qsum2, Qsum3v1);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2Terms->Fill(Qsum1v2, Qsum2, Qsum3v2);
		    }
		  }
		  //////////////////////////////////////
		  // Momentum resolution calculation
		  /*if(fillIndex3==0 && fMCcase){
		    WInput = 1.0;
		    if(ch1==ch2 && ch1==ch3){// same charge
		      WInput = MCWeight3D(kTRUE, jj, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SC"))->Fill(q3MC, WInput);
		      ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SCden"))->Fill(q3MC);
		    }else {// mixed charge
		      if(bin1==bin2) WInput = MCWeight3D(kFALSE, jj, fFixedLambdaBinMomRes, firstQMC, secondQMC, thirdQMC);
		      else WInput = MCWeight3D(kFALSE, 6-jj, fFixedLambdaBinMomRes, thirdQMC, secondQMC, firstQMC);// thirdQMC is ss
		      
		      if(bin1==bin2){
			if(jj==2){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2MCden"))->Fill(q3MC);
			}else if(jj==3){
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm3MCden"))->Fill(q3MC);
			}else {
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MC"))->Fill(q3MC, WInput);
			  ((TH1D*)fOutputList->FindObject("fMCWeight3DTerm4MCden"))->Fill(q3MC);
			}
		      }
		    }
		    //
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fIdeal->Fill(firstQMC, secondQMC, thirdQMC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fSmeared->Fill(firstQ, secondQ, thirdQ, WInput);
		    //
		    if(ch1==ch2 && ch1==ch3){
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsIdeal->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsSmeared->Fill(Qsum1v1, Qsum2, Qsum3v1, WInput);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsIdeal->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsSmeared->Fill(Qsum1v2, Qsum2, Qsum3v2, WInput);
		      //
		      if(qinv12MC > fQLowerCut && qinv13MC > fQLowerCut && qinv23MC > fQLowerCut){
			// does not really matter if MC or real data triplets are used to average 1/K2...but better to use umsmeared values
			Float_t InteractingQ=0;
			if(part==1) {InteractingQ=qinv12MC;}// 12 from SE
			else {InteractingQ=qinv13MC;}// 13 from SE
			Double_t K2 = FSICorrelationTherm2(+1,+1, InteractingQ);// K2 from Therminator source
			WInput = MCWeight3D(kTRUE, jj, 25, firstQMC, secondQMC, thirdQMC);// pure 2-pion (lambda=1)
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsSumK2->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput/K2);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsSumK2->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput/K2);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd1TermsEnK2->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, WInput);
			Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].f4VectProd2TermsEnK2->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, WInput);
		      }
		    }
		  }// fMCcase
		  */
		}
	      }
	      
	    }else {// config 3: All particles from different events
	      
	      // "enhancement" differs from 1.0 only when Qinv goes over fQcut
	      //Float_t enhancement=1.0;
	      //Int_t nUnderCut=0;
	      //if(qinv13<fQcut[qCutBin13]) nUnderCut++;
	      //if(qinv23<fQcut[qCutBin23]) nUnderCut++;
	      //if(nUnderCut==0) enhancement = (1+1+1)/1.;// 1 LowQ pair
	      //if(nUnderCut==1) enhancement = (1+2)/2.;// 2 LowQ pair
	      //if(nUnderCut==2) enhancement = 1.;// 3 LowQ pair
	      
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 3, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	      
	      if(ch1==ch2 && ch1==ch3 && fillIndex3==0) {
		FourVectProdTerms(pVect1, pVect2, pVect3, Qsum1v1, Qsum2, Qsum3v1, Qsum1v2, Qsum3v2);// 4-vector product sums
		if(!fMCcase) ((TH3D*)fOutputList->FindObject("fKt3DistTerm5"))->Fill(fMbin+1, transK3, q3);
	      }	      
	      //if(fMCcase && ch1==ch2 && ch1==ch3 && fillIndex3==0) FourVectProdTerms(pVect1MC, pVect2MC, pVect3MC, Qsum1v1MC, Qsum2MC, Qsum3v1MC, Qsum1v2MC, Qsum3v2MC);
	      
	      if(fillIndex3 <= 2){
		ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, part, 5, firstQ, secondQ, thirdQ);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fTerms3->Fill(firstQ, secondQ, thirdQ);
		if(fillIndex3==0 && ch1==ch2 && ch1==ch3 && fMCcase==kFALSE){
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd1Terms->Fill(Qsum1v1, Qsum2, Qsum3v1);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd2Terms->Fill(Qsum1v2, Qsum2, Qsum3v2);
		}
		//////////////////////////////////////
		// Momentum resolution calculation
		/*if(fillIndex3==0 && fMCcase){
		  ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, part, 5, firstQMC, secondQMC, thirdQMC);
		  
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fIdeal->Fill(firstQMC, secondQMC, thirdQMC);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fSmeared->Fill(firstQ, secondQ, thirdQ);
		  if(ch1==ch2 && ch1==ch3){
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd1TermsIdeal->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd1TermsSmeared->Fill(Qsum1v1, Qsum2, Qsum3v1);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd2TermsIdeal->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].f4VectProd2TermsSmeared->Fill(Qsum1v2, Qsum2, Qsum3v2);
		    
		  }
		}// fMCcase
		*/
	      }
	     
	      if(fillIndex3 !=0) continue;// only calculate TPN for pi-pi-pi
	      if(ch1!=ch2 || ch1!=ch3) continue;// only calcualte TPN for ss
	      
	      
	      //if(fMCcase) continue;// only calcualte TPN for real data
	      if(!fGenerateSignal){
		GetWeight(pVect1, pVect2, pVect1, pVect2, weight12, weight12Err);
		GetWeight(pVect1, pVect3, pVect1, pVect3, weight13, weight13Err);
		GetWeight(pVect2, pVect3, pVect2, pVect3, weight23, weight23Err);
	      }else {
		GetWeight(pVect1, pVect2Flat, pVect1, pVect2, weight12, weight12Err);
		GetWeight(pVect1, pVect3Flat, pVect1, pVect3, weight13, weight13Err);
		GetWeight(pVect2Flat, pVect3Flat, pVect2, pVect3, weight23, weight23Err);
	      }
	      if(sqrt(fabs(weight12*weight13*weight23)) > 1.0) {
		if(fMbin==0 && bin1==0) {
		  ((TH3F*)fOutputList->FindObject("fTPNRejects1"))->Fill(qinv12, qinv13, qinv23, sqrt(fabs(weight12*weight13*weight23)));
		}
		continue;// weight should never be larger than 1
	      }
	      	 
	      
	      Float_t myDamp = fDampStart + (fDampStep)*fFixedLambdaBinr3;// lambdabin=0.52 for v1 draft, 0.7 is more realistic
	      Int_t denIndex = 0;
	      Int_t momResIndex = rIndexForTPNMomRes*kNDampValues + fFixedLambdaBinMomRes;// lambdabin=0.52 for v1 draft, 0.4 is more realistic

	      Float_t coulCorr12 = FSICorrelationTherm2(+1,+1, qinv12);
	      Float_t coulCorr13 = FSICorrelationTherm2(+1,+1, qinv13);
	      Float_t coulCorr23 = FSICorrelationTherm2(+1,+1, qinv23);
	      if(coulCorr12 < 0.1 || coulCorr13 < 0.1 || coulCorr23 < 0.1) {// Safety check
		if(fMbin==0 && bin1==0) {
		  ((TH3F*)fOutputList->FindObject("fTPNRejects2"))->Fill(qinv12, qinv13, qinv23, sqrt(fabs(weight12*weight13*weight23)));
		}
		continue;
	      }
	      Float_t MomResCorr12=1.0, MomResCorr13=1.0, MomResCorr23=1.0;
	      if(!fGenerateSignal && !fMCcase) {
		Int_t momBin12 = fMomResC2->GetYaxis()->FindBin(qinv12);
		Int_t momBin13 = fMomResC2->GetYaxis()->FindBin(qinv13);
		Int_t momBin23 = fMomResC2->GetYaxis()->FindBin(qinv23);		  
		if(momBin12 >= kQbins) momBin12 = kQbins-1;
		if(momBin13 >= kQbins) momBin13 = kQbins-1;
		if(momBin23 >= kQbins) momBin23 = kQbins-1;
		MomResCorr12 = fMomResC2->GetBinContent(momResIndex+1, momBin12);
		MomResCorr13 = fMomResC2->GetBinContent(momResIndex+1, momBin13);
		MomResCorr23 = fMomResC2->GetBinContent(momResIndex+1, momBin23);
		if(MomResCorr12 > 1.2 || MomResCorr13 > 1.2 || MomResCorr23 > 1.2) {// Safety check
		  if(fMbin==0 && bin1==0) {
		    ((TH3F*)fOutputList->FindObject("fTPNRejects3"))->Fill(qinv12, qinv13, qinv23, sqrt(fabs(weight12*weight13*weight23)));
		  }
		  continue;
		}
	      }
	      weight12CC = ((weight12+1)*MomResCorr12 - myDamp*coulCorr12 - (1-myDamp));
	      weight12CC /= coulCorr12*myDamp;
	      weight13CC = ((weight13+1)*MomResCorr13 - myDamp*coulCorr13 - (1-myDamp));
	      weight13CC /= coulCorr13*myDamp;
	      weight23CC = ((weight23+1)*MomResCorr23 - myDamp*coulCorr23 - (1-myDamp));
	      weight23CC /= coulCorr23*myDamp;
	      
	      if(weight12CC < 0 || weight13CC < 0 || weight23CC < 0) {
		if(fMbin==0 && bin1==0) {
		  weightTotal = sqrt(fabs(weight12CC*weight13CC*weight23CC));
		  ((TH3F*)fOutputList->FindObject("fTPNRejects4"))->Fill(qinv12, qinv13, qinv23, weightTotal);
		}
		continue;// C2^QS can never be less than unity
	      }
	      
	      /////////////////////////////////////////////////////
	      weightTotal = sqrt(weight12CC*weight13CC*weight23CC);
	      /////////////////////////////////////////////////////

	      if(weightTotal > 1.5) {
		if(fMbin==0 && bin1==0) {
		  ((TH3F*)fOutputList->FindObject("fTPNRejects5"))->Fill(qinv12, qinv13, qinv23, weightTotal);
		}
		continue;// C2^QS never be greater than 1.0 in theory. Can be slightly larger than 1.0 with fluctuations
	      }

	    
	      
	      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].fTwoPartNorm->Fill(qinv12, qinv13, qinv23, weightTotal);
	      
	      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd1TwoPartNorm->Fill(Qsum1v1, Qsum2, Qsum3v1, weightTotal);
	      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd2TwoPartNorm->Fill(Qsum1v2, Qsum2, Qsum3v2, weightTotal);
	      /*if(fMCcase){
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd1TwoPartNormIdeal->Fill(Qsum1v1MC, Qsum2MC, Qsum3v1MC, weightTotal);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd1TwoPartNormSmeared->Fill(Qsum1v1, Qsum2, Qsum3v1, weightTotal);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd2TwoPartNormIdeal->Fill(Qsum1v2MC, Qsum2MC, Qsum3v2MC, weightTotal);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].DT[denIndex].f4VectProd2TwoPartNormSmeared->Fill(Qsum1v2, Qsum2, Qsum3v2, weightTotal);
		}*/
		  
		  
	      // Save cpu time and memory by skipping r3 denominator calculation below.  den errors are negligible compared to num errors.
	      /*
		if(weightTotal > 0.0001){// tiny numbers cause a Float_ting point exception below
		weightTotalErr = pow((weight12Err*coulCorr12)*weight13CC*weight23CC,2);
		weightTotalErr += pow(weight12CC*(weight13Err*coulCorr13)*weight23CC,2);
		weightTotalErr += pow(weight12CC*weight13CC*(weight23Err*coulCorr23),2);
		weightTotalErr /= pow(2*weightTotal,2);
		
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].TwoPartNormErr->Fill(denIndex, q3, weightTotalErr);
		}
	      */
	      
	      //}// damp iter
	      //}// R iter
	      
	      
	      
	    }
	  }// end 3rd particle
	}// en2
	
	
      }// p1
    }//en1
    
    ///////////////////
  }// end of PdensityPairs
  
  
 
    
   
    
  
  ////////////////////////////////////////////////////////
  // Pdensity Method with Explicit Loops
  if(fPdensityExplicitLoop){
    
    ////////////////////////////////////
    // 2nd, 3rd, and 4th order Correlations
    
    // First Particle
    for (Int_t i=0; i<myTracks; i++) {
      ch1 = Int_t( ((fEvt)->fTracks[i].fCharge + 1)/2. );
      pVect1[0] = (fEvt)->fTracks[i].fEaccepted;
      pVect1[1] = (fEvt)->fTracks[i].fP[0];
      pVect1[2] = (fEvt)->fTracks[i].fP[1];
      pVect1[3] = (fEvt)->fTracks[i].fP[2];
      key1 = (fEvt)->fTracks[i].fKey;

      // Second Event
      for(Int_t en2=0; en2<fEventsToMix+1; en2++){
	Int_t startbin2=0;
	if(en2==0) startbin2=i+1;
	
	// Second Particle
	for (Int_t j=startbin2; j<(fEvt+en2)->fNtracks; j++) {
	  ch2 = Int_t( ((fEvt+en2)->fTracks[j].fCharge + 1)/2. );
	  pVect2[0] = (fEvt+en2)->fTracks[j].fEaccepted;
	  pVect2[1] = (fEvt+en2)->fTracks[j].fP[0];
	  pVect2[2] = (fEvt+en2)->fTracks[j].fP[1];
	  pVect2[3] = (fEvt+en2)->fTracks[j].fP[2];
	  key2 = (fEvt+en2)->fTracks[j].fKey;

	  Short_t fillIndex12 = FillIndex2part(key1+key2);
	  qinv12 = GetQinv(fillIndex12, pVect1, pVect2);
	  
	  if(qinv12 < fQLowerCut) continue;

	  
	  // 2-particle part is filled always during pair creator
	  
	  // Third Event
	  for(Int_t en3=en2; en3<fEventsToMix+1; en3++){
	    Int_t startbin3=0;
	    if(en3==en2) startbin3=j+1;
	    else startbin3=0;
	    
	    
	    // Third Particle
	    for (Int_t k=startbin3; k<(fEvt+en3)->fNtracks; k++) {
	      ch3 = Int_t( ((fEvt+en3)->fTracks[k].fCharge + 1)/2. );
	      pVect3[0] = (fEvt+en3)->fTracks[k].fEaccepted;
	      pVect3[1] = (fEvt+en3)->fTracks[k].fP[0];
	      pVect3[2] = (fEvt+en3)->fTracks[k].fP[1];
	      pVect3[3] = (fEvt+en3)->fTracks[k].fP[2];
	      key3 = (fEvt+en3)->fTracks[k].fKey;
	      
	      Short_t fillIndex3 = FillIndex3part(key1+key2+key3);
	      Short_t fillIndex13 = FillIndex2part(key1+key3);
	      qinv13 = GetQinv(fillIndex13, pVect1, pVect3);
	      Short_t fillIndex23 = FillIndex2part(key2+key3);
	      qinv23 = GetQinv(fillIndex23, pVect2, pVect3);

	      
	      if(qinv13 < fQLowerCut) continue;
	      if(qinv23 < fQLowerCut) continue;
	      
	   	      
	      q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));
	      
	      Short_t normBin12 = SetNormBin(fillIndex12);
	      Short_t normBin13 = SetNormBin(fillIndex13);
	      Short_t normBin23 = SetNormBin(fillIndex23);

	      
	      if(en3==0 && en2==0) {// 123
		SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 0, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
		
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fExplicit3->Fill(q3);// 123
		
		if((qinv12>=fNormQcutLow[normBin12]) && (qinv13>=fNormQcutLow[normBin13]) && (qinv23>=fNormQcutLow[normBin23])) {
		  if((qinv12<fNormQcutHigh[normBin12]) && (qinv13<fNormQcutHigh[normBin13]) && (qinv23<fNormQcutHigh[normBin23])) {
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fNormEx3->Fill(0.);
		  }
		}
		
	      }else if((en2==0 && en3==1) ) {// 12-3, 13-2, 23-1
		Float_t gFact=1;
		
		Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
		SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);

			
		if(fill2){
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[1].fExplicit3->Fill(q3, gFact);// 12
		  if((qinv12>=fNormQcutLow[normBin12]) && (qinv13>=fNormQcutLow[normBin13]) && (qinv23>=fNormQcutLow[normBin23])) {
		    if((qinv12<fNormQcutHigh[normBin12]) && (qinv13<fNormQcutHigh[normBin13]) && (qinv23<fNormQcutHigh[normBin23])) {
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[1].fNormEx3->Fill(0.);
		    }
		  }
		}
		if(fill3){
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[2].fExplicit3->Fill(q3, gFact);// 12
		  if((qinv12>=fNormQcutLow[normBin12]) && (qinv13>=fNormQcutLow[normBin13]) && (qinv23>=fNormQcutLow[normBin23])) {
		    if((qinv12<fNormQcutHigh[normBin12]) && (qinv13<fNormQcutHigh[normBin13]) && (qinv23<fNormQcutHigh[normBin23])) {
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[2].fNormEx3->Fill(0.);
		    }
		  }
		}
		if(fill4){
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[3].fExplicit3->Fill(q3, gFact);// 12
		  if((qinv12>=fNormQcutLow[normBin12]) && (qinv13>=fNormQcutLow[normBin13]) && (qinv23>=fNormQcutLow[normBin23])) {
		    if((qinv12<fNormQcutHigh[normBin12]) && (qinv13<fNormQcutHigh[normBin13]) && (qinv23<fNormQcutHigh[normBin23])) {
		      Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[3].fNormEx3->Fill(0.);
		    }
		  }
		}
		
	      }else if(en2==1 && en3==2){// all uncorrelated events
		SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 3, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
		
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fExplicit3->Fill(q3);
		if((qinv12>=fNormQcutLow[normBin12]) && (qinv13>=fNormQcutLow[normBin13]) && (qinv23>=fNormQcutLow[normBin23])) {
		  if((qinv12<fNormQcutHigh[normBin12]) && (qinv13<fNormQcutHigh[normBin13]) && (qinv23<fNormQcutHigh[normBin23])) {
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fNormEx3->Fill(0.);
		  }
		}	
		Short_t qCutBin12 = SetQcutBin(fillIndex12);
		Short_t qCutBin13 = SetQcutBin(fillIndex13);
		Short_t qCutBin23 = SetQcutBin(fillIndex23);
		
		if( (qinv12 < fQcut[qCutBin12]) || (qinv13 < fQcut[qCutBin13]) || (qinv23 < fQcut[qCutBin23])){
		  
		  Int_t nUnderCut=0;
		  if(qinv12<fQcut[qCutBin12]) nUnderCut++;
		  if(qinv13<fQcut[qCutBin13]) nUnderCut++;
		  if(qinv23<fQcut[qCutBin23]) nUnderCut++;
		  
		}
		
	      }else {}
  	    
	      
  	    }// 3rd particle
	  }// 3rd event
	  
	}// 2nd particle
      }// 2nd event

    }// 1st particle
		 
    
  
       
  }// End of PdensityExplicit
  

 
  
  // Post output data.
  PostData(1, fOutputList);
  
}
//________________________________________________________________________
void AliChaoticity::Terminate(Option_t *) 
{
  // Called once at the end of the query
 
  cout<<"Done"<<endl;

}
//________________________________________________________________________
Bool_t AliChaoticity::AcceptPair(AliChaoticityTrackStruct first, AliChaoticityTrackStruct second)
{
 
  if(fabs(first.fEta-second.fEta) > fMinSepPair) return kTRUE;
  
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

  //cout<<deltaphi<<"  "<<fMinSepPair<<"  "<<fMinSepTPCEntranceEta<<endl;
  if(deltaphi < fMinSepPair) return kFALSE;// Min Separation
  
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

  if(deltaphi < fMinSepPair) return kFALSE;// Min Separation

   
  //
  
  Int_t ncl1 = first.fClusterMap.GetNbits();
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
  
  
  return kTRUE;
  

}
//________________________________________________________________________
Float_t AliChaoticity::GamovFactor(Int_t chargeBin1, Int_t chargeBin2, Float_t qinv)
{
  Float_t arg = G_Coeff/qinv;
  
  if(chargeBin1==chargeBin2) return (exp(arg)-1)/(arg);
  else {return (exp(-arg)-1)/(-arg);}
  
}
//________________________________________________________________________
void AliChaoticity::Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
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
Short_t AliChaoticity::FillIndex2part(Short_t key){

  if(key==2) return 0;// pi-pi
  else if(key==11) return 1;// pi-k
  else if(key==101) return 2;// pi-p
  else if(key==20) return 3;// k-k
  else if(key==110) return 4;// k-p
  else return 5;// p-p
}
//________________________________________________________________________
Short_t AliChaoticity::FillIndex3part(Short_t key){
  
  if(key==3) return 0;// pi-pi-pi
  else if(key==12) return 1;// pi-pi-k
  else if(key==21) return 2;// k-k-pi
  else if(key==102) return 3;// pi-pi-p
  else if(key==201) return 4;// p-p-pi
  else if(key==111) return 5;// pi-k-p
  else if(key==30) return 6;// k-k-k
  else if(key==120) return 7;// k-k-p
  else if(key==210) return 8;// p-p-k
  else return 9;// p-p-p
  
}
//________________________________________________________________________
Short_t AliChaoticity::SetQcutBin(Short_t fi){// fi=FillIndex
  if(fi <= 2) return 0;
  else if(fi==3) return 1;
  else return 2;
}
//________________________________________________________________________
Short_t AliChaoticity::SetNormBin(Short_t fi){// fi=FillIndex
  if(fi==0) return 0;
  else if(fi <= 2) return 1;
  else return 2;
}
//________________________________________________________________________
void AliChaoticity::SetFillBins2(Short_t fi, Short_t key1, Short_t key2, Int_t c1, Int_t c2, Int_t &b1, Int_t &b2){
  
  if(fi==0 || fi==3 || fi==5){// Identical species
    if((c1+c2)==1) {b1=0; b2=1;}// Re-assign to merge degenerate histos
    else {b1=c1; b2=c2;}
  }else {// Mixed species
    if(key1 < key2) { b1=c1; b2=c2;}
    else {b1=c2; b2=c1;}
  }
  
}
//________________________________________________________________________
void AliChaoticity::SetFillBins3(Short_t fi, Short_t key1, Short_t key2, Short_t key3, Int_t c1, Int_t c2, Int_t c3, Short_t part, Int_t &b1, Int_t &b2, Int_t &b3, Bool_t &fill2, Bool_t &fill3, Bool_t &fill4){
  
  
  // seSS, seSK, SE_keysum only used to determine which terms to fill (only used for terms 2-4)
  // part only matters for terms 2-4
  Bool_t seSS=kFALSE;
  Bool_t seSK=kFALSE;
  Short_t seKeySum=0;// only used for pi-k-p case
  if(part==1) {// default case (irrelevant for term 1 and term 5)
    if(c1==c2) seSS=kTRUE;
    if(key1==key2) seSK=kTRUE;
    seKeySum = key1+key2;
  }
  if(part==2){
    if(c1==c3) seSS=kTRUE;
    if(key1==key3) seSK=kTRUE;
    seKeySum = key1+key3;
  }
  
  
  // fill2, fill3, fill4 are only used for Cumulant Terms 2,3,4
  
  if(fi==0 || fi==6 || fi==9){// Identical species
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
  }else if(fi != 5){// all the rest except pi-k-p
    if(key1==key2){
      b3=c3;
      if( (c1+c2)==1) {b1=0; b2=1;}
      else {b1=c1; b2=c2;}
    }else if(key1==key3){
      b3=c2;
      if( (c1+c3)==1) {b1=0; b2=1;}
      else {b1=c1; b2=c3;}
    }else {// Key2==Key3
      b3=c1;
      if( (c2+c3)==1) {b1=0; b2=1;}
      else {b1=c2; b2=c3;}
    }
    //////////////////////////////
    if(seSK) fill2=kTRUE;// Same keys from Same Event
    else {// Different keys from Same Event
      if( (c1+c2+c3)==1) {
	if(b3==0) {
	  if(seSS) fill3=kTRUE;
	  else fill4=kTRUE;
	}else{fill3=kTRUE; fill4=kTRUE;}// b3=1 so fill both
      }else if( (c1+c2+c3)==2) {
	if(b3==1) {
	  if(seSS) fill4=kTRUE;
	  else fill3=kTRUE;
	}else{fill3=kTRUE; fill4=kTRUE;}// b3=0 so fill both
      }else{fill3=kTRUE; fill4=kTRUE;}// all same charge so fill both
    }
    /////////////////////////////
  }else {// pi-k-p  (no charge ordering applies since all are unique)
    if(key1==1){
      if(key2==10) {b1=c1; b2=c2; b3=c3;}// pi-k-p
      else {b1=c1; b2=c3; b3=c2;}// pi-p-k
    }else if(key1==10){
      if(key2==1) {b1=c2; b2=c1; b3=c3;}// k-pi-p
      else {b1=c3; b2=c1; b3=c2;}// k-p-pi
    }else {// key1==100
      if(key2==1) {b1=c2; b2=c3; b3=c1;}// p-pi-k
      else {b1=c3; b2=c2; b3=c1;}// p-k-pi
    }
    ////////////////////////////////////
    if(seKeySum==11) fill2=kTRUE;
    else if(seKeySum==101) fill3=kTRUE;
    else fill4=kTRUE;
    ////////////////////////////////////
  }
  
}
//________________________________________________________________________
void AliChaoticity::ArrangeQs(Short_t fi, Short_t key1, Short_t key2, Short_t key3, Int_t c1, Int_t c2, Int_t c3, Float_t q12, Float_t q13, Float_t q23, Short_t part, Short_t term, Float_t &fQ, Float_t &sQ, Float_t &tQ){
 
  // for terms 2-4: start by setting q12(part 1) or q13(part 2)
  if(fi==0 || fi==6 || fi==9){// Identical species
    if( (c1+c2+c3)==1) {// fQ=ss, sQ=os, tQ=os
      if(term==1 || term==5){
	if(c1==c2) {fQ=q12; sQ=q13; tQ=q23;}
	else if(c1==c3) {fQ=q13; sQ=q12; tQ=q23;}
	else {fQ=q23; sQ=q12; tQ=q13;}
      }else if(term==2 && part==1){
	fQ=q12; sQ=q13; tQ=q23;
      }else if(term==2 && part==2){
	fQ=q13; sQ=q12; tQ=q23;
      }else if(term==3 && part==1){
	sQ=q12; 
	if(c1==c3) {fQ=q13; tQ=q23;}
	else {fQ=q23; tQ=q13;}
      }else if(term==3 && part==2){
	sQ=q13;
	if(c1==c2) {fQ=q12; tQ=q23;}
	else {fQ=q23; tQ=q12;}
      }else if(term==4 && part==1){
	tQ=q12;
	if(c1==c3) {fQ=q13; sQ=q23;}
	else {fQ=q23; sQ=q13;}
      }else if(term==4 && part==2){
	tQ=q13;
	if(c1==c2) {fQ=q12; sQ=q23;}
	else {fQ=q23; sQ=q12;}
      }else cout<<"problem!!!!!!!!!!!!!"<<endl;
    }else if( (c1+c2+c3)==2) {// fQ=os, sQ=os, tQ=ss
      if(term==1 || term==5){
	if(c1==c2) {tQ=q12; sQ=q13; fQ=q23;}
	else if(c1==c3) {tQ=q13; sQ=q12; fQ=q23;}
	else {tQ=q23; sQ=q12; fQ=q13;}
      }else if(term==2 && part==1){
	fQ=q12; 
	if(c1==c3) {tQ=q13; sQ=q23;}
	else {tQ=q23; sQ=q13;}
      }else if(term==2 && part==2){
	fQ=q13; 
	if(c1==c2) {tQ=q12; sQ=q23;}
	else {tQ=q23; sQ=q12;}
      }else if(term==3 && part==1){
	sQ=q12; 
	if(c1==c3) {tQ=q13; fQ=q23;}
	else {tQ=q23; fQ=q13;}
      }else if(term==3 && part==2){
	sQ=q13; 
	if(c1==c2) {tQ=q12; fQ=q23;}
	else {tQ=q23; fQ=q12;}
      }else if(term==4 && part==1){
	tQ=q12; sQ=q13; fQ=q23;
      }else if(term==4 && part==2){
	tQ=q13; sQ=q12; fQ=q23;
      }else cout<<"problem!!!!!!!!!!!!!"<<endl;
    }else {// fQ=ss, sQ=ss, tQ=ss
      if(term==1 || term==5) {fQ=q12; sQ=q13; tQ=q23;}
      else if(term==2 && part==1) {fQ=q12; sQ=q13; tQ=q23;}
      else if(term==2 && part==2) {fQ=q13; sQ=q12; tQ=q23;}
      else if(term==3 && part==1) {sQ=q12; fQ=q13; tQ=q23;}
      else if(term==3 && part==2) {sQ=q13; fQ=q12; tQ=q23;}
      else if(term==4 && part==1) {tQ=q12; fQ=q13; sQ=q23;}
      else if(term==4 && part==2) {tQ=q13; fQ=q12; sQ=q23;}
    }
  }else if(fi != 5){// all the rest except pi-k-p	
    if(key1==key2){
      fQ=q12;
      if(c1==c2){
	// cases not explicity shown below are not possible
	if(term==1 || term==5) {sQ=q13; tQ=q23;}
	else if(term==2 && part==1) {sQ=q13; tQ=q23;}
	else if(term==3 && part==2) {sQ=q13; tQ=q23;}
	else if(term==4 && part==2) {tQ=q13; sQ=q23;}
	else cout<<"problem!!!!!!!!!!!!!"<<endl;
      }else if(c3==0){
	if(c1==c3) {sQ=q13; tQ=q23;}
	else {sQ=q23; tQ=q13;}
      }else {//c3==1
	if(c1==c3) {tQ=q13; sQ=q23;}
	else {tQ=q23; sQ=q13;}
      }
    }else if(key1==key3){
      fQ=q13;
      if(c1==c3){
	// cases not explicity shown below are not possible
	if(term==1 || term==5) {sQ=q12; tQ=q23;}
	else if(term==2 && part==2) {sQ=q12; tQ=q23;}
	else if(term==3 && part==1) {sQ=q12; tQ=q23;}
	else if(term==4 && part==1) {tQ=q12; sQ=q23;}
	else cout<<"problem!!!!!!!!!!!!!!!!!!!!!!"<<endl;
      }else if(c2==0){
	if(c1==c2) {sQ=q12; tQ=q23;}
	else {sQ=q23; tQ=q12;}
      }else {//c2==1
	if(c1==c2) {tQ=q12; sQ=q23;}
	else {tQ=q23; sQ=q12;}
      }
    }else {// key2==key3
      fQ=q23;
      if(c2==c3){
	// cases not explicity shown below are not possible
	if(term==1 || term==5) {sQ=q12; tQ=q13;}
	else if(term==3 && part==1) {sQ=q12; tQ=q13;}
	else if(term==3 && part==2) {sQ=q13; tQ=q12;}
	else if(term==4 && part==1) {tQ=q12; sQ=q13;}
	else if(term==4 && part==2) {tQ=q13; sQ=q12;}
	else cout<<"problem!!!!!!!!!!!!!!!!!!!!!!"<<endl;
      }else if(c1==0){
	if(c1==c2) {sQ=q12; tQ=q13;}
	else {sQ=q13; tQ=q12;}
      }else {//c1==1
	if(c1==c2) {tQ=q12; sQ=q13;}
	else {tQ=q13; sQ=q12;}
      }
    }
  }else {// pi-k-p
    if(key1==1){
      if(key2==10) {fQ=q12; sQ=q13; tQ=q23;}// pi-k-p
      else {fQ=q13; sQ=q12; tQ=q23;}// pi-p-k
    }else if(key1==10){
      if(key2==1) {fQ=q12; sQ=q23; tQ=q13;}// k-pi-p
      else {fQ=q13; sQ=q23; tQ=q12;}// k-p-pi
    }else {// key1==100
      if(key2==1) {fQ=q23; sQ=q12; tQ=q13;}// p-pi-k
      else {fQ=q23; sQ=q13; tQ=q12;}// p-k-pi
    }
    
  }


}
//________________________________________________________________________
Float_t AliChaoticity::GetQinv(Short_t fi, Float_t track1[], Float_t track2[]){
  
  Float_t qinv=1.0;
  
  if(fi==0 || fi==3 || fi==5){// identical masses
    qinv = sqrt( pow(track1[1]-track2[1],2) + pow(track1[2]-track2[2],2) + pow(track1[3]-track2[3],2) - pow(track1[0]-track2[0],2));
  }else{// different masses
    Float_t px = track1[1] + track2[1]; 
    Float_t py = track1[2] + track2[2]; 
    Float_t pz = track1[3] + track2[3];    
    Float_t pSquared = pow(track1[0]+track2[0],2) - px*px - py*py - pz*pz;
    Float_t deltaDOTsum = (track1[0]-track2[0])*(track1[0]+track2[0]);
    deltaDOTsum -= (track1[1]-track2[1])*px + (track1[2]-track2[2])*py + (track1[3]-track2[3])*pz;
    
    qinv =  pow( (track1[1]-track2[1]) - deltaDOTsum*px/(pSquared),2);
    qinv += pow( (track1[2]-track2[2]) - deltaDOTsum*py/(pSquared),2);
    qinv += pow( (track1[3]-track2[3]) - deltaDOTsum*pz/(pSquared),2);
    qinv -= pow( (track1[0]-track2[0]) - deltaDOTsum*(track1[0]+track2[0])/(pSquared),2);
    qinv = sqrt(qinv);
  }
  
  return qinv;
  
}
//________________________________________________________________________
void AliChaoticity::GetQosl(Float_t track1[], Float_t track2[], Float_t& qout, Float_t& qside, Float_t& qlong){
 
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
void AliChaoticity::SetWeightArrays(Bool_t legoCase, TH3F *histos[AliChaoticity::fKbinsT][AliChaoticity::fCentBins]){

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
void AliChaoticity::GetWeight(Float_t track1[], Float_t track2[], Float_t track3[], Float_t track4[], Float_t& wgt, Float_t& wgtErr){
  
  Float_t kt=sqrt( pow(track3[1]+track4[1],2) + pow(track3[2]+track4[2],2))/2.;
  //
  Float_t qOut=0,qSide=0,qLong=0;
  GetQosl(track1, track2, qOut, qSide, qLong);
  qOut = fabs(qOut);
  qSide = fabs(qSide);
  qLong = fabs(qLong);
  Float_t wd=0, xd=0, yd=0, zd=0;
  //
  
  if(kt < fKmeanT[0]) {fKtIndexL=0; fKtIndexH=1;}// fKtIndexL=0; fKtIndexH=0; no extrapolation
  else if(kt >= fKmeanT[fKbinsT-1]) {fKtIndexL=fKbinsT-2; fKtIndexH=fKbinsT-1;}// fKtIndexL=fKbinsT-1; fKtIndexH=fKbinsT-1; no extrapolation
  else {
    for(Int_t i=0; i<fKbinsT-1; i++){
      if((kt >= fKmeanT[i]) && (kt < fKmeanT[i+1])) {fKtIndexL=i; fKtIndexH=i+1; break;}
    }
  }
  //
  /////////
  if(qOut < fQmean[0]) {fQoIndexL=0; fQoIndexH=0; xd=0;}
  else if(qOut >= fQmean[kQbinsWeights-1]) {fQoIndexL=kQbinsWeights-1; fQoIndexH=kQbinsWeights-1; xd=0;}
  else {
    for(Int_t i=0; i<kQbinsWeights-1; i++){
      if((qOut >= fQmean[i]) && (qOut < fQmean[i+1])) {fQoIndexL=i; fQoIndexH=i+1; break;}
    }
    xd = (qOut-fQmean[fQoIndexL])/(fQmean[fQoIndexH]-fQmean[fQoIndexL]);
  }
  //
  if(qSide < fQmean[0]) {fQsIndexL=0; fQsIndexH=0; yd=0;}
  else if(qSide >= fQmean[kQbinsWeights-1]) {fQsIndexL=kQbinsWeights-1; fQsIndexH=kQbinsWeights-1; yd=0;}
  else {
    for(Int_t i=0; i<kQbinsWeights-1; i++){
      if((qSide >= fQmean[i]) && (qSide < fQmean[i+1])) {fQsIndexL=i; fQsIndexH=i+1; break;}
    }
    yd = (qSide-fQmean[fQsIndexL])/(fQmean[fQsIndexH]-fQmean[fQsIndexL]);
  }
  //
  if(qLong < fQmean[0]) {fQlIndexL=0; fQlIndexH=0; zd=0;}
  else if(qLong >= fQmean[kQbinsWeights-1]) {fQlIndexL=kQbinsWeights-1; fQlIndexH=kQbinsWeights-1; zd=0;}
  else {
    for(Int_t i=0; i<kQbinsWeights-1; i++){
      if((qLong >= fQmean[i]) && (qLong < fQmean[i+1])) {fQlIndexL=i; fQlIndexH=i+1; break;}
    }
    zd = (qLong-fQmean[fQlIndexL])/(fQmean[fQlIndexH]-fQmean[fQlIndexL]);
  }
  //

  
  //Float_t min = fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1,fQsIndexH+1,fQlIndexH+1);
  Float_t minErr = fNormWeight[fKtIndexL][fMbin]->GetBinError(fQoIndexH+1,fQsIndexH+1,fQlIndexH+1);
  /*
  Float_t deltaW=0;
  // kt
  deltaW += (fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1) - min)*(kt-fKmeanT[fKtIndexL])/((fKstepT[fKtIndexL]+fKstepT[fKtIndexH])/2.);
  // Qo 
  deltaW += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1) - min)*(qOut-fQmean[fQoIndexL])/fQstepWeights;
  // Qs
  deltaW += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1) - min)*(qSide-fQmean[fQsIndexL])/fQstepWeights;
  // Ql
  deltaW += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1) - min)*(qLong-fQmean[fQlIndexL])/fQstepWeights;
  //
  wgt = min + deltaW;
  */
  
 
  //
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
  

  ////
  
  // Denominator errors negligible compared to numerator so do not waste cpu time below.  
  //Float_t deltaWErr=0;
  // Kt
  /*
  deltaWErr += (fNormWeight[fKtIndexH][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexH+1) - minErr)*(kt-fKmeanT[fKtIndexL])/((fKstepT[fKtIndexL]+fKstepT[fKtIndexH])/2.);
  // Qo 
  deltaWErr += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexL+1, fQsIndexH+1, fQlIndexH+1) - minErr)*(qOut-fQmean[fQoIndexL])/fQstepWeights;
  // Qs
  deltaWErr += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexL+1, fQlIndexH+1) - minErr)*(qSide-fQmean[fQsIndexL])/fQstepWeights;
  // Ql
  deltaWErr += (fNormWeight[fKtIndexL][fMbin]->GetBinContent(fQoIndexH+1, fQsIndexH+1, fQlIndexL+1) - minErr)*(qLong-fQmean[fQlIndexL])/fQstepWeights;
  */
  wgtErr = minErr;
  
 
}
//________________________________________________________________________
Float_t AliChaoticity::MCWeight(Int_t charge1, Int_t charge2, Int_t r, Int_t dampIndex, Float_t qinv){
  
  Float_t radius = Float_t(r)/0.19733;// convert to GeV (starts at 5 fm, was 3 fm)
  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t coulCorr12 = FSICorrelationTherm2(charge1, charge2, qinv);
  if(charge1==charge2){
    Float_t arg=qinv*radius;
    Float_t EW = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg,3) - 12.*arg);
    EW += kappa4/(24.*pow(2.,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
    return ((1-myDamp) + myDamp*(1 + exp(-pow(qinv*radius,2))*pow(EW,2))*coulCorr12);
  }else {
    return ((1-myDamp) + myDamp*coulCorr12);
  }
    
}
//________________________________________________________________________
Float_t AliChaoticity::MCWeightOSL(Int_t charge1, Int_t charge2, Int_t r, Int_t dampIndex, Float_t qinv, Float_t qo, Float_t qs, Float_t ql){
  
  Float_t radiusOut = Float_t(r)/0.19733;// convert to GeV (starts at 5 fm, was 3 fm)
  Float_t radiusSide = radiusOut;
  Float_t radiusLong = radiusOut;
  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t coulCorr12 = FSICorrelationTherm2(charge1, charge2, qinv);
  if(charge1==charge2){
    return ((1-myDamp) + myDamp*(1 + exp(-pow(qo*radiusOut,2)) * exp(-pow(qs*radiusSide,2)) * exp(-pow(ql*radiusLong,2)))*coulCorr12);
  }else {
    return ((1-myDamp) + myDamp*coulCorr12);
  }
    
}

//________________________________________________________________________
Float_t AliChaoticity::MCWeight3D(Bool_t SameCharge, Int_t term, Int_t dampIndex, Float_t q12, Float_t q13, Float_t q23){
  if(term==5) return 1.0;
  
  Float_t radius=fRMax;// was in terms of bins starting at 3 fm Gaussian source
  //if(fMbin<=1) {}
  //else if(fMbin<=3) {radius = radius-1;}
  //else if(fMbin<=5) {radius = radius-2;}
  //else {radius = radius-3;}
  
  radius /= 0.19733;

  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t fc = sqrt(myDamp);
  
  if(SameCharge){// all three of the same charge
    Float_t coulCorr12 = FSICorrelationTherm2(+1,+1, q12);// K2
    Float_t coulCorr13 = FSICorrelationTherm2(+1,+1, q13);// K2
    Float_t coulCorr23 = FSICorrelationTherm2(+1,+1, q23);// K2
    Float_t arg12=q12*radius;
    Float_t arg13=q13*radius;
    Float_t arg23=q23*radius;
    Float_t EW12 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
    EW12 += kappa4/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
    Float_t EW13 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
    EW13 += kappa4/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
    Float_t EW23 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
    EW23 += kappa4/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
    if(term==1){
      Float_t c3QS = 1 + exp(-pow(q12*radius,2))*pow(EW12,2) + exp(-pow(q13*radius,2))*pow(EW13,2) + exp(-pow(q23*radius,2))*pow(EW23,2);
      c3QS += 2*exp(-pow(radius,2)*(pow(q12,2) + pow(q13,2) + pow(q23,2))/2.)*EW12*EW13*EW23;
      Float_t w123 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      w123 += pow(fc,2)*(1-fc)*(1+exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12;
      w123 += pow(fc,2)*(1-fc)*(1+exp(-pow(q13*radius,2))*pow(EW13,2))*coulCorr13;
      w123 += pow(fc,2)*(1-fc)*(1+exp(-pow(q23*radius,2))*pow(EW23,2))*coulCorr23;
      w123 += pow(fc,3)*c3QS*coulCorr12*coulCorr13*coulCorr23;// was pow(fc,3)*c3QS*FSICorrelationOmega0(kTRUE, q12, q13, q23)
      return w123;
    }else if(term==2){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12);
    }else if(term==3){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q13*radius,2))*pow(EW13,2))*coulCorr13);
    }else if(term==4){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q23*radius,2))*pow(EW23,2))*coulCorr23);
    }else return 1.0;
  
  }else{// mixed charge case pair 12 always treated as ss
    Float_t coulCorr12 = FSICorrelationTherm2(+1,+1, q12);// K2 ss
    Float_t coulCorr13 = FSICorrelationTherm2(+1,-1, q13);// K2 os
    Float_t coulCorr23 = FSICorrelationTherm2(+1,-1, q23);// K2 os
    Float_t arg12=q12*radius;
    Float_t EW12 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
    EW12 += kappa4/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
    if(term==1){
      Float_t c3QS = 1 + exp(-pow(q12*radius,2))*pow(EW12,2);
      Float_t w123 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      w123 += pow(fc,2)*(1-fc)*(1+exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12;
      w123 += pow(fc,2)*(1-fc)*coulCorr13;
      w123 += pow(fc,2)*(1-fc)*coulCorr23;
      w123 += pow(fc,3)*c3QS*coulCorr12*coulCorr13*coulCorr23;// was pow(fc,3)*c3QS*FSICorrelationOmega0(kFALSE, q12, q13, q23)
      return w123;
    }else if(term==2){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12);
    }else if(term==3){
      return ((1-myDamp) + myDamp*coulCorr13);
    }else if(term==4){
      return ((1-myDamp) + myDamp*coulCorr23);
    }else return 1.0;
  }
  
}
//________________________________________________________________________
void AliChaoticity::SetMomResCorrections(Bool_t legoCase, TH2D *temp2D){
  
 
  if(legoCase){
    cout<<"LEGO call to SetMomResCorrections"<<endl;
    fMomResC2 = (TH2D*)temp2D->Clone();
    fMomResC2->SetDirectory(0);
  }else {
    TFile *momResFile = new TFile("MomResFile.root","READ");
    if(!momResFile->IsOpen()) {
      cout<<"No momentum resolution file found"<<endl;
      AliFatal("No momentum resolution file found.  Kill process.");
    }else {cout<<"Good Momentum Resolution File Found!"<<endl;}
    
    TH2D *temp2D2 = (TH2D*)momResFile->Get("MomResHisto_pp");
    fMomResC2 = (TH2D*)temp2D2->Clone();
    fMomResC2->SetDirectory(0);
            
    momResFile->Close();
  }

  // fMomResC2->GetBinContent(1,5) should be ~1.007
  if(fMomResC2->GetBinContent(1,5) > 1.2) AliFatal("AliChaoticity: SetMomResCorrections Problem");// Additional Safety check
  if(fMomResC2->GetBinContent(1,5) < 0.95) AliFatal("AliChaoticity: SetMomResCorrections Problem");// Additional Safety check

  for(Int_t bx=1; bx<=fMomResC2->GetNbinsX(); bx++){
    for(Int_t by=1; by<=fMomResC2->GetNbinsX(); by++){
      if(fMomResC2->GetBinContent(bx,by) > 1.5) fMomResC2->SetBinContent(bx,by, 1.5);// Maximum is ~1.02 
      if(fMomResC2->GetBinContent(bx,by) < 0.95) fMomResC2->SetBinContent(bx,by, 0.95);// Minimum is ~0.98
    }
  }

  cout<<"Done reading momentum resolution file"<<endl;
}
//________________________________________________________________________
void AliChaoticity::SetFSICorrelations(Bool_t legoCase, TH2D *temp2DGaus[2], TH2D *temp2DTherm[2], TH3D *temp3Dos[6], TH3D *temp3Dss[6]){
  // read in 2-particle and 3-particle FSI correlations = K2 & K3
  // 2-particle input histo from file is binned in qinv.  3-particle in qinv of each pair
  if(legoCase){
    cout<<"LEGO call to SetFSICorrelations"<<endl;
    fFSI2SS[0] = (TH2D*)temp2DGaus[0]->Clone();
    fFSI2OS[0] = (TH2D*)temp2DGaus[1]->Clone();
    fFSI2SS[1] = (TH2D*)temp2DTherm[0]->Clone();
    fFSI2OS[1] = (TH2D*)temp2DTherm[1]->Clone();
    //
    fFSI2SS[0]->SetDirectory(0);
    fFSI2OS[0]->SetDirectory(0);
    fFSI2SS[1]->SetDirectory(0);
    fFSI2OS[1]->SetDirectory(0);

    for(Int_t CB=0; CB<6; CB++) {
      fFSIOmega0OS[CB] = (TH3D*)temp3Dos[CB]->Clone();
      fFSIOmega0SS[CB] = (TH3D*)temp3Dss[CB]->Clone();
      //
      fFSIOmega0OS[CB]->SetDirectory(0);
      fFSIOmega0SS[CB]->SetDirectory(0);
    }
  }else {
    cout<<"non LEGO call to SetFSICorrelations"<<endl;
    TFile *fsifile = new TFile("KFile.root","READ");
    if(!fsifile->IsOpen()) {
      cout<<"No FSI file found"<<endl;
      AliFatal("No FSI file found.  Kill process.");
    }else {cout<<"Good FSI File Found!"<<endl;}
    
    TH2D *temphisto2GausSS = (TH2D*)fsifile->Get("K2ssG");
    TH2D *temphisto2GausOS = (TH2D*)fsifile->Get("K2osG");
    TH2D *temphisto2ThermSS = (TH2D*)fsifile->Get("K2ssT");
    TH2D *temphisto2ThermOS = (TH2D*)fsifile->Get("K2osT");
    TH3D *temphisto3OS[6];
    TH3D *temphisto3SS[6];
    for(Int_t CB=0; CB<6; CB++) {
      TString *nameK3SS = new TString("K3ss_");
      *nameK3SS += CB;
      temphisto3SS[CB] = (TH3D*)fsifile->Get(nameK3SS->Data());
      //
      TString *nameK3OS = new TString("K3os_");
      *nameK3OS += CB;
      temphisto3OS[CB] = (TH3D*)fsifile->Get(nameK3OS->Data());
    }

    fFSI2SS[0] = (TH2D*)temphisto2GausSS->Clone();
    fFSI2OS[0] = (TH2D*)temphisto2GausOS->Clone();
    fFSI2SS[1] = (TH2D*)temphisto2ThermSS->Clone();
    fFSI2OS[1] = (TH2D*)temphisto2ThermOS->Clone();
    fFSI2SS[0]->SetDirectory(0);
    fFSI2OS[0]->SetDirectory(0);
    fFSI2SS[1]->SetDirectory(0);
    fFSI2OS[1]->SetDirectory(0);

    for(Int_t CB=0; CB<6; CB++) {
      fFSIOmega0SS[CB] = (TH3D*)temphisto3SS[CB]->Clone();
      fFSIOmega0OS[CB] = (TH3D*)temphisto3OS[CB]->Clone();
      fFSIOmega0SS[CB]->SetDirectory(0);
      fFSIOmega0OS[CB]->SetDirectory(0);
    }   
    //
    
    fsifile->Close();
  }
  /*
  // condition FSI histogram for edge effects
  for(Int_t CB=0; CB<6; CB++){
    for(Int_t ii=1; ii<=fFSIOmega0SS[CB]->GetNbinsX(); ii++){
      for(Int_t jj=1; jj<=fFSIOmega0SS[CB]->GetNbinsY(); jj++){
	for(Int_t kk=1; kk<=fFSIOmega0SS[CB]->GetNbinsZ(); kk++){
	  
	  if(fFSIOmega0SS[CB]->GetBinContent(ii,jj,kk) <=0){
	    Double_t Q12 = fFSIOmega0SS[CB]->GetXaxis()->GetBinCenter(ii);
	    Double_t Q23 = fFSIOmega0SS[CB]->GetYaxis()->GetBinCenter(jj);
	    Double_t Q13 = fFSIOmega0SS[CB]->GetZaxis()->GetBinCenter(kk);
	    //
	    Int_t Q12bin=ii;
	    Int_t Q23bin=jj;
	    Int_t Q13bin=kk;
	    Int_t AC=0;//Adjust Counter
	    Int_t AClimit=10;// maximum bin shift
	    if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q12bin++; AC++;}}
	    if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q12bin--; AC++;}}
	    //
	    if(Q13 < sqrt(pow(Q12,2)+pow(Q23,2) - 2*Q12*Q23)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q13bin++; AC++;}}
	    if(Q13 > sqrt(pow(Q12,2)+pow(Q23,2) + 2*Q12*Q23)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q13bin--; AC++;}}
	    //
	    if(Q23 < sqrt(pow(Q12,2)+pow(Q13,2) - 2*Q12*Q13)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q23bin++; AC++;}}
	    if(Q23 > sqrt(pow(Q12,2)+pow(Q13,2) + 2*Q12*Q13)) {while(fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0 && AC<AClimit) {Q23bin--; AC++;}}
	    
	    // Save cpu time by setting empty cell contents (edge effects) to nearest non-zero cell (these cells are not used very often anyway.)
	    if(AC==AClimit) {
	      fFSIOmega0SS[CB]->SetBinContent(ii,jj,kk, 1.0);
	      fFSIOmega0OS[CB]->SetBinContent(ii,jj,kk, 1.0);
	    }else {
	      fFSIOmega0SS[CB]->SetBinContent(ii,jj,kk, fFSIOmega0SS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin));
	      fFSIOmega0OS[CB]->SetBinContent(ii,jj,kk, fFSIOmega0OS[CB]->GetBinContent(Q12bin, Q23bin, Q13bin));
	    }
	  }
	  
	}
      }
    }
  }
  */
  // fFSI2SS[1]->GetBinContent(1,2) should be ~0.32
  if(fFSI2SS[1]->GetBinContent(1,2) > 1.0) AliFatal("AliChaoticity: SetFSICorrelations Problem");// Additional Safety check
  if(fFSI2SS[1]->GetBinContent(1,2) < 0.1) AliFatal("AliChaoticity: SetFSICorrelations Problem");// Additional Safety check

  for(Int_t ii=1; ii<=fFSI2SS[0]->GetNbinsX(); ii++){
      for(Int_t jj=1; jj<=fFSI2SS[0]->GetNbinsY(); jj++){
	if(fFSI2SS[0]->GetBinContent(ii,jj) > 1.0) fFSI2SS[0]->SetBinContent(ii,jj, 1.0);
	if(fFSI2SS[1]->GetBinContent(ii,jj) > 1.0) fFSI2SS[1]->SetBinContent(ii,jj, 1.0);
	if(fFSI2OS[0]->GetBinContent(ii,jj) > 10.0) fFSI2OS[0]->SetBinContent(ii,jj, 10.0);
	if(fFSI2OS[1]->GetBinContent(ii,jj) > 10.0) fFSI2OS[1]->SetBinContent(ii,jj, 10.0);
	//
	if(fFSI2SS[0]->GetBinContent(ii,jj) < 0.05) fFSI2SS[0]->SetBinContent(ii,jj, 0.05);
	if(fFSI2SS[1]->GetBinContent(ii,jj) < 0.05) fFSI2SS[1]->SetBinContent(ii,jj, 0.05);
	if(fFSI2OS[0]->GetBinContent(ii,jj) < 0.9) fFSI2OS[0]->SetBinContent(ii,jj, 0.9);
	if(fFSI2OS[1]->GetBinContent(ii,jj) < 0.9) fFSI2OS[1]->SetBinContent(ii,jj, 0.9);
      }
  }

  cout<<"Done reading FSI file"<<endl;
}
//________________________________________________________________________
Float_t AliChaoticity::FSICorrelationGaus2(Int_t charge1, Int_t charge2, Int_t rIndex, Float_t qinv){
  // returns 2-particle Coulomb correlations = K2
  if(rIndex >= fRVALUES) return 1.0;
  Int_t qbinL = fFSI2SS[0]->GetYaxis()->FindBin(qinv-fFSI2SS[0]->GetYaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSI2SS[0]->GetNbinsY()) return 1.0;
  
  Float_t slope=0;
  if(charge1==charge2){
    slope = fFSI2SS[0]->GetBinContent(rIndex+1, qbinL) - fFSI2SS[0]->GetBinContent(rIndex+1, qbinH);
    slope /= fFSI2SS[0]->GetYaxis()->GetBinCenter(qbinL) - fFSI2SS[0]->GetYaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2SS[0]->GetYaxis()->GetBinCenter(qbinL)) + fFSI2SS[0]->GetBinContent(rIndex+1, qbinL));
  }else {
    slope = fFSI2OS[0]->GetBinContent(rIndex+1, qbinL) - fFSI2OS[0]->GetBinContent(rIndex+1, qbinH);
    slope /= fFSI2OS[0]->GetYaxis()->GetBinCenter(qbinL) - fFSI2OS[0]->GetYaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2OS[0]->GetYaxis()->GetBinCenter(qbinL)) + fFSI2OS[0]->GetBinContent(rIndex+1, qbinL));
  }
}
//________________________________________________________________________
Float_t AliChaoticity::FSICorrelationTherm2(Int_t charge1, Int_t charge2, Float_t qinv){
  // returns 2-particle Coulomb correlations = K2
  Int_t qbinL = fFSI2SS[1]->GetYaxis()->FindBin(qinv-fFSI2SS[1]->GetYaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSI2SS[1]->GetNbinsY()) return 1.0;
  
  Float_t slope=0;
  if(charge1==charge2){
    slope = fFSI2SS[1]->GetBinContent(fFSIbin+1, qbinL) - fFSI2SS[1]->GetBinContent(fFSIbin+1, qbinH);
    slope /= fFSI2SS[1]->GetYaxis()->GetBinCenter(qbinL) - fFSI2SS[1]->GetYaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2SS[1]->GetYaxis()->GetBinCenter(qbinL)) + fFSI2SS[1]->GetBinContent(fFSIbin+1, qbinL));
  }else {
    slope = fFSI2OS[1]->GetBinContent(fFSIbin+1, qbinL) - fFSI2OS[1]->GetBinContent(fFSIbin+1, qbinH);
    slope /= fFSI2OS[1]->GetYaxis()->GetBinCenter(qbinL) - fFSI2OS[1]->GetYaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2OS[1]->GetYaxis()->GetBinCenter(qbinL)) + fFSI2OS[1]->GetBinContent(fFSIbin+1, qbinL));
  }
}
//________________________________________________________________________
Double_t AliChaoticity::FSICorrelationOmega0(Bool_t SameCharge, Double_t Q12, Double_t Q13, Double_t Q23){
  // returns 3d 3-particle Coulomb Correlation = K3
  Int_t Q12bin = fFSIOmega0SS[fFSIbin]->GetXaxis()->FindBin(Q12);
  Int_t Q13bin = fFSIOmega0SS[fFSIbin]->GetZaxis()->FindBin(Q13);
  Int_t Q23bin = fFSIOmega0SS[fFSIbin]->GetYaxis()->FindBin(Q23);
  Int_t index12L = int(fabs(Q12 - fFSI2SS[1]->GetYaxis()->GetBinWidth(1)/2.)/(fFSI2SS[1]->GetYaxis()->GetBinWidth(1)));
  Int_t index12H = index12L+1;
  Int_t index13L = int(fabs(Q13 - fFSI2SS[1]->GetYaxis()->GetBinWidth(1)/2.)/(fFSI2SS[1]->GetYaxis()->GetBinWidth(1)));
  Int_t index13H = index13L+1;
  Int_t index23L = int(fabs(Q23 - fFSI2SS[1]->GetYaxis()->GetBinWidth(1)/2.)/(fFSI2SS[1]->GetYaxis()->GetBinWidth(1)));
  Int_t index23H = index23L+1;

  if(SameCharge){
    if(fFSIOmega0SS[fFSIbin]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0) return 1.0;
    Double_t base = fFSIOmega0SS[fFSIbin]->GetBinContent(index12L+1, index23L+1, index13L+1);
    Double_t InterPolated = 0;
    Double_t slope12 = fFSIOmega0SS[fFSIbin]->GetBinContent(index12H+1, index23L+1, index13L+1);
    slope12 -= base;
    slope12 /= fFSIOmega0SS[fFSIbin]->GetXaxis()->GetBinWidth(1);
    InterPolated += slope12*fabs(Q12 - fFSIOmega0SS[fFSIbin]->GetXaxis()->GetBinCenter(index12L+1));
    Double_t slope23 = fFSIOmega0SS[fFSIbin]->GetBinContent(index12L+1, index23H+1, index13L+1);
    slope23 -= base;
    slope23 /= fFSIOmega0SS[fFSIbin]->GetYaxis()->GetBinWidth(1);
    InterPolated += slope23*fabs(Q23 - fFSIOmega0SS[fFSIbin]->GetYaxis()->GetBinCenter(index23L+1));
    Double_t slope13 = fFSIOmega0SS[fFSIbin]->GetBinContent(index12L+1, index23L+1, index13H+1);
    slope13 -= base;
    slope13 /= fFSIOmega0SS[fFSIbin]->GetZaxis()->GetBinWidth(1);
    InterPolated += slope13*fabs(Q13 - fFSIOmega0SS[fFSIbin]->GetZaxis()->GetBinCenter(index13L+1));
    if( (base+InterPolated) <= 0) return 1.0;
    return (base+InterPolated);
  
  }else{// mixed charge. Q12 is always designated as the same-charge pair
    if(fFSIOmega0OS[fFSIbin]->GetBinContent(Q12bin, Q23bin, Q13bin) <=0) return 1.0;
    Double_t base = fFSIOmega0OS[fFSIbin]->GetBinContent(index12L+1, index23H+1, index13H+1);
    Double_t InterPolated = 0;
    Double_t slope12 = fFSIOmega0OS[fFSIbin]->GetBinContent(index12H+1, index23H+1, index13H+1);
    slope12 -= base;
    slope12 /= fFSIOmega0OS[fFSIbin]->GetXaxis()->GetBinWidth(1);
    InterPolated += slope12*fabs(Q12 - fFSIOmega0OS[fFSIbin]->GetXaxis()->GetBinCenter(index12L+1));
    Double_t slope23 = fFSIOmega0OS[fFSIbin]->GetBinContent(index12L+1, index23L+1, index13H+1);
    slope23 -= base;
    slope23 /= fFSIOmega0OS[fFSIbin]->GetYaxis()->GetBinWidth(1);
    InterPolated += slope23*fabs(Q23 - fFSIOmega0OS[fFSIbin]->GetYaxis()->GetBinCenter(index23L+1));
    Double_t slope13 = fFSIOmega0OS[fFSIbin]->GetBinContent(index12L+1, index23H+1, index13L+1);
    slope13 -= base;
    slope13 /= fFSIOmega0OS[fFSIbin]->GetZaxis()->GetBinWidth(1);
    InterPolated += slope13*fabs(Q13 - fFSIOmega0OS[fFSIbin]->GetZaxis()->GetBinCenter(index13L+1));
    if( (base+InterPolated) <= 0) return 1.0;
    return (base+InterPolated);
    
  }
}
//________________________________________________________________________
void AliChaoticity::FourVectProdTerms(Float_t pV1[], Float_t pV2[], Float_t pV3[], Float_t& QS1v1, Float_t& QS2, Float_t& QS3v1, Float_t& QS1v2, Float_t& QS3v2){
  QS1v1 = (pV1[0]-pV2[0])*(pV2[1]-pV3[1]) - (pV1[1]-pV2[1])*(pV2[0]-pV3[0]);
  QS1v1 += (pV1[0]-pV2[0])*(pV2[2]-pV3[2]) - (pV1[2]-pV2[2])*(pV2[0]-pV3[0]);
  QS1v1 += (pV1[0]-pV2[0])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[0]-pV3[0]);
  QS2 = (pV1[1]-pV2[1])*(pV2[2]-pV3[2]) - (pV1[2]-pV2[2])*(pV2[1]-pV3[1]);
  QS3v1 = (pV1[1]-pV2[1])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[1]-pV3[1]);
  QS3v1 += (pV1[2]-pV2[2])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[2]-pV3[2]);
  //
  QS1v2 = (pV1[0]-pV2[0])*(pV2[1]-pV3[1]) - (pV1[1]-pV2[1])*(pV2[0]-pV3[0]);
  QS1v2 += (pV1[0]-pV2[0])*(pV2[2]-pV3[2]) - (pV1[2]-pV2[2])*(pV2[0]-pV3[0]);
  QS3v2 = (pV1[1]-pV2[1])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[1]-pV3[1]);
  QS3v2 += (pV1[0]-pV2[0])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[0]-pV3[0]);
  QS3v2 += (pV1[2]-pV2[2])*(pV2[3]-pV3[3]) - (pV1[3]-pV2[3])*(pV2[2]-pV3[2]);  
}
//________________________________________________________________________
