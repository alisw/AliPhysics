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
#include "AliAODTracklets.h"
#include "AliAnalysisUtils.h"

#include "AliThreePionRadii.h"

#define PI 3.1415927
#define G_Coeff 0.006399 // 2*pi*alpha*M_pion
#define kappa3 0.2 // kappa3 Edgeworth coefficient (non-Gaussian features of C2)
#define kappa4 0.45 // kappa4 Edgeworth coefficient (non-Gaussian features of C2)


// Author: Dhevan Gangadharan

ClassImp(AliThreePionRadii)

//________________________________________________________________________
AliThreePionRadii::AliThreePionRadii():
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
  fPdensityPairCut(kTRUE),
  fRMax(11),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fKt3bins(1),
  fV0Mbinning(kFALSE),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fTriggerType(0),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fQcut(),
  fQLowerCut(0),
  fQlimitC2(2.0),
  fQbinsC2(400),
  fNormQcutLow(),
  fNormQcutHigh(),
  fKupperBound(0),
  fQupperBound(0),
  fQbins(0),
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
  fNormPairs()
{
  // Default constructor
  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2 = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	      //
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = 0x0;
	      
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTermsQ3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared = 0x0;

	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
      
    }// ED
  }// Mbin
  
  // Initialize 3-pion FSI histograms
  for(Int_t i=0; i<10; i++){
    fFSI2SS[i]=0x0; 
    fFSI2OS[i]=0x0;
  }

}
//________________________________________________________________________
AliThreePionRadii::AliThreePionRadii(const Char_t *name) 
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
  fPdensityPairCut(kTRUE),
  fRMax(11),
  fFilterBit(7),
  fMaxChi2NDF(10),
  fMinTPCncls(0),
  fBfield(0),
  fMbin(0),
  fFSIindex(0),
  fEDbin(0),
  fMbins(fCentBins),
  fMultLimit(0),
  fKt3bins(1),
  fV0Mbinning(kFALSE),
  fCentBinLowLimit(0),
  fCentBinHighLimit(1),
  fTriggerType(0),
  fEventCounter(0),
  fEventsToMix(0),
  fZvertexBins(0),
  fMultLimits(),
  fQcut(),
  fQLowerCut(0),
  fQlimitC2(2.0),
  fQbinsC2(400),
  fNormQcutLow(),
  fNormQcutHigh(),
  fKupperBound(0),
  fQupperBound(0),
  fQbins(0),
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
  fNormPairs()
{
  // Main constructor
  fAODcase=kTRUE;
  fPdensityPairCut = kTRUE;
  

  for(Int_t mb=0; mb<fMbins; mb++){
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2 = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = 0x0;
	      //
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = 0x0;
	      Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = 0x0;
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTermsQ3 = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal = 0x0;
		Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared = 0x0;
	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
            
    }// ED
  }// Mbin
  
  // Initialize 3-pion FSI histograms
  for(Int_t i=0; i<10; i++){
    fFSI2SS[i]=0x0; 
    fFSI2OS[i]=0x0;
  }


  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliThreePionRadii::AliThreePionRadii(const AliThreePionRadii &obj) 
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
    fPdensityPairCut(obj.fPdensityPairCut),
    fRMax(obj.fRMax),
    fFilterBit(obj.fFilterBit),
    fMaxChi2NDF(obj.fMaxChi2NDF),
    fMinTPCncls(obj.fMinTPCncls),
    fBfield(obj.fBfield),
    fMbin(obj.fMbin),
    fFSIindex(obj.fFSIindex),
    fEDbin(obj.fEDbin),
    fMbins(obj.fMbins),
    fMultLimit(obj.fMultLimit),
    fKt3bins(obj.fKt3bins),
    fV0Mbinning(obj.fV0Mbinning),
    fCentBinLowLimit(obj.fCentBinLowLimit),
    fCentBinHighLimit(obj.fCentBinHighLimit),
    fTriggerType(obj.fTriggerType),
    fEventCounter(obj.fEventCounter),
    fEventsToMix(obj.fEventsToMix),
    fZvertexBins(obj.fZvertexBins),
    fMultLimits(),
    fQcut(),
    fQLowerCut(obj.fQLowerCut),
    fQlimitC2(obj.fQlimitC2),
    fQbinsC2(obj.fQbinsC2),
    fNormQcutLow(),
    fNormQcutHigh(),
    fKupperBound(obj.fKupperBound),
    fQupperBound(obj.fQupperBound),
    fQbins(obj.fQbins),
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
    fNormPairs()
{
  // Copy Constructor
  
  for(Int_t i=0; i<10; i++){
    fFSI2SS[i]=obj.fFSI2SS[i]; 
    fFSI2OS[i]=obj.fFSI2OS[i];
  }

}
//________________________________________________________________________
AliThreePionRadii &AliThreePionRadii::operator=(const AliThreePionRadii &obj) 
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
  fPbPbcase = obj.fPbPbcase; 
  fGenerateSignal = obj.fGenerateSignal;
  fGeneratorOnly = obj.fGeneratorOnly;
  fPdensityPairCut = obj.fPdensityPairCut;
  fRMax = obj.fRMax;
  fFilterBit = obj.fFilterBit;
  fMaxChi2NDF = obj.fMaxChi2NDF;
  fMinTPCncls = obj.fMinTPCncls;
  fBfield = obj.fBfield;
  fMbin = obj.fMbin;
  fFSIindex = obj.fFSIindex;
  fEDbin = obj.fEDbin;
  fMbins = obj.fMbins;
  fMultLimit = obj.fMultLimit;
  fKt3bins = obj.fKt3bins;
  fV0Mbinning = obj.fV0Mbinning;
  fCentBinLowLimit = obj.fCentBinLowLimit;
  fCentBinHighLimit = obj.fCentBinHighLimit;
  fTriggerType = obj.fTriggerType;
  fEventCounter = obj.fEventCounter;
  fEventsToMix = obj.fEventsToMix;
  fZvertexBins = obj.fZvertexBins;
  fQLowerCut = obj.fQLowerCut;
  fQlimitC2 = obj.fQlimitC2;
  fQbinsC2 = obj.fQbinsC2;
  fKupperBound = obj.fKupperBound;
  fQupperBound = obj.fQupperBound;
  fQbins = obj.fQbins;
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
  fDummyB = obj.fDummyB;
 
 
  for(Int_t i=0; i<10; i++){
    fFSI2SS[i]=obj.fFSI2SS[i]; 
    fFSI2OS[i]=obj.fFSI2OS[i];
  }
  
  return (*this);
}
//________________________________________________________________________
AliThreePionRadii::~AliThreePionRadii()
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
    if(fPbPbcase && ((mb < fCentBinLowLimit) || (mb > fCentBinHighLimit))) continue;
    for(Int_t edB=0; edB<fEDbins; edB++){
      for(Int_t c1=0; c1<2; c1++){
	for(Int_t c2=0; c2<2; c2++){
	  for(Int_t sc=0; sc<kSCLimit2; sc++){
	    for(Int_t term=0; term<2; term++){
	      
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW;

	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared;
	      //
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv;
	      if(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW) delete Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW;
	    }// term_2
	  }// SC_2
	  
	  for(Int_t c3=0; c3<2; c3++){
	    for(Int_t sc=0; sc<kSCLimit3; sc++){
	      for(Int_t term=0; term<5; term++){
		
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3;
		if(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3) delete Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3;
				
	      }// term_3
	    }// SC_3
	  }//c3
	}//c2
      }//c1
      
    }// ED
  }// Mbin
  
   
  for(Int_t i=0; i<10; i++){
    if(fFSI2SS[i]) delete fFSI2SS[i]; 
    if(fFSI2OS[i]) delete fFSI2OS[i];
  }
  
}
//________________________________________________________________________
void AliThreePionRadii::ParInit()
{
  cout<<"AliThreePionRadii MyInit() call"<<endl;
  cout<<"lego:"<<fLEGO<<"  MCcase:"<<fMCcase<<"  PbPbcase:"<<fPbPbcase<<"  GenSignal:"<<fGenerateSignal<<"  CentLow:"<<fCentBinLowLimit<<"  CentHigh:"<<fCentBinHighLimit<<"  RMax:"<<fRMax<<"  FB:"<<fFilterBit<<"  MaxChi2/NDF:"<<fMaxChi2NDF<<"  MinTPCncls:"<<fMinTPCncls<<"  MinPairSepEta:"<<fMinSepPairEta<<"  MinPairSepPhi:"<<fMinSepPairPhi<<"  NsigTPC:"<<fSigmaCutTPC<<"  NsigTOF:"<<fSigmaCutTOF<<endl;

  fRandomNumber = new TRandom3();
  fRandomNumber->SetSeed(0);
    
  //
  fEventCounter=0;
  if(fPdensityPairCut) fEventsToMix=2;
  else fEventsToMix=0;
  fZvertexBins=2;//2
  
  fTPCTOFboundry = 0.6;// TPC pid used below this momentum, TOF above but below TOF_boundry
  fTOFboundry = 2.1;// TOF pid used below this momentum
  
  ////////////////////////////////////////////////
  // PadRow Pair Cuts
  fShareQuality = .5;// max
  fShareFraction = .05;// max
  ////////////////////////////////////////////////
  
  
  //fMultLimits[0]=0, fMultLimits[1]=5, fMultLimits[2]=10, fMultLimits[3]=15, fMultLimits[4]=20, fMultLimits[5]=25;
  //fMultLimits[6]=30, fMultLimits[7]=35, fMultLimits[8]=40, fMultLimits[9]=45, fMultLimits[10]=kMultLimitPP;
  
  fMultLimits[0]=0, fMultLimits[1]=5; fMultLimits[2]=10; fMultLimits[3]=15; fMultLimits[4]=20;
  fMultLimits[5]=30, fMultLimits[6]=40; fMultLimits[7]=50; fMultLimits[8]=70; fMultLimits[9]=100;
  fMultLimits[10]=150, fMultLimits[11]=200; fMultLimits[12]=260; fMultLimits[13]=320; fMultLimits[14]=400;
  fMultLimits[15]=500, fMultLimits[16]=600; fMultLimits[17]=700; fMultLimits[18]=850; fMultLimits[19]=1050;
  fMultLimits[20]=2000;
  
  
  if(fPbPbcase && fCentBinLowLimit < 6) {// PbPb 0-30%, was 0-50%
    fMultLimit=kMultLimitPbPb; 
    fMbins=fCentBins; 
    fQcut[0]=0.1;//pi-pi, pi-k, pi-p
    fQcut[1]=0.1;//k-k
    fQcut[2]=0.6;//the rest
    fNormQcutLow[0] = 0.15;// was 0.15
    fNormQcutHigh[0] = 0.175;// was 0.175
    fNormQcutLow[1] = 1.34;//1.34
    fNormQcutHigh[1] = 1.4;//1.4
    fNormQcutLow[2] = 1.1;//1.1
    fNormQcutHigh[2] = 1.4;//1.4
    //
    fQlimitC2 = 2.0;
    fQbinsC2 = 400;
    fQupperBound = fQcut[0];
    fQbins = kQbins;
    //
    fDampStart = 0.5;
    fDampStep = 0.02;
  }else if(fPbPbcase && fCentBinLowLimit >= 6) {// PbPb 30-100%, was 50-100%
    fMultLimit=kMultLimitPbPb;
    fMbins=fCentBins;
    fQcut[0]=0.2;//pi-pi, pi-k, pi-p
    fQcut[1]=0.2;//k-k
    fQcut[2]=1.2;//the rest
    fNormQcutLow[0] = 0.3;// was 0.3
    fNormQcutHigh[0] = 0.35;// was 0.35
    fNormQcutLow[1] = 1.34;//1.34
    fNormQcutHigh[1] = 1.4;//1.4
    fNormQcutLow[2] = 1.1;//1.1
    fNormQcutHigh[2] = 1.4;//1.4
    //
    fQlimitC2 = 2.0;
    fQbinsC2 = 400;
    fQupperBound = fQcut[0];
    fQbins = 2*kQbins;
    //
    fDampStart = 0.5;
    fDampStep = 0.02;
  }else {// pp or pPb
    fMultLimit=kMultLimitPP;
    fMbins=fCentBins;
    fQcut[0]=2.0;// 0.4
    fQcut[1]=2.0;
    fQcut[2]=2.0;
    fNormQcutLow[0] = 1.0;// was 1.0
    fNormQcutHigh[0] = 1.2;// was 1.2
    fNormQcutLow[1] = 1.0;
    fNormQcutHigh[1] = 1.2;
    fNormQcutLow[2] = 1.0;
    fNormQcutHigh[2] = 1.2;
    //
    fQlimitC2 = 2.0;
    fQbinsC2 = 200;
    fQupperBound = 0.5;// was 0.4
    fQbins = kQbinsPP;
    //
    fDampStart = 0.5;
    fDampStep = 0.02;
  }

  fQLowerCut = 0.005;
  fKupperBound = 1.0;
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
  }
  
  /////////////////////////////////////////////
  /////////////////////////////////////////////
  
}
//________________________________________________________________________
void AliThreePionRadii::UserCreateOutputObjects()
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
  
  TH1F *fMultDist4 = new TH1F("fMultDist4","Multiplicity Distribution",fMultLimit,-.5,fMultLimit-.5);
  fMultDist4->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist4);
  
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
  TH3F *fPairsShareFracDPhiNum = new TH3F("fPairsShareFracDPhiNum","",10,-.5,9.5, 159,0,1, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiNum);
  TH3F *fPairsShareFracDPhiDen = new TH3F("fPairsShareFracDPhiDen","",10,-.5,9.5, 159,0,1, 600,-0.3,0.3);
  if(fMCcase) fOutputList->Add(fPairsShareFracDPhiDen);
  TH3D* fPairsPadRowNum = new TH3D("fPairsPadRowNum","", 20,0,1, 159,0,1, 40,0,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowNum);
  TH3D* fPairsPadRowDen = new TH3D("fPairsPadRowDen","", 20,0,1, 159,0,1, 40,0,0.2);
  if(fMCcase) fOutputList->Add(fPairsPadRowDen);


  
  TH3D *fPrimarySCPionPairs = new TH3D("fPrimarySCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fPrimarySCPionPairs);
  TH3D *fAllSCPionPairs = new TH3D("fAllSCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fAllSCPionPairs);
  TH3D *fPrimaryMCPionPairs = new TH3D("fPrimaryMCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fPrimaryMCPionPairs);
  TH3D *fAllMCPionPairs = new TH3D("fAllMCPionPairs","",fMbins,.5,fMbins+.5, 20,0,1, 20,0,0.2);
  if(fMCcase) fOutputList->Add(fAllMCPionPairs);
  //
  TH3D *fMuonContamSmearedNum2 = new TH3D("fMuonContamSmearedNum2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamSmearedNum2);
  TH3D *fMuonContamSmearedDen2 = new TH3D("fMuonContamSmearedDen2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamSmearedDen2);
  TH3D *fMuonContamIdealNum2 = new TH3D("fMuonContamIdealNum2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamIdealNum2);
  TH3D *fMuonContamIdealDen2 = new TH3D("fMuonContamIdealDen2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamIdealDen2);
  //
  TH3D *fMuonContamSmearedNum3 = new TH3D("fMuonContamSmearedNum3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamSmearedNum3);
  TH3D *fMuonContamSmearedDen3 = new TH3D("fMuonContamSmearedDen3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamSmearedDen3);
  TH3D *fMuonContamIdealNum3 = new TH3D("fMuonContamIdealNum3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamIdealNum3);
  TH3D *fMuonContamIdealDen3 = new TH3D("fMuonContamIdealDen3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonContamIdealDen3);
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
  TH3D *fMuonPionK2 = new TH3D("fMuonPionK2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  TH3D *fPionPionK2 = new TH3D("fPionPionK2","",2,-0.5,1.5, 20,0,1, 100,0,0.5);
  TH3D *fMuonPionK3 = new TH3D("fMuonPionK3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  TH3D *fPionPionK3 = new TH3D("fPionPionK3","",2,-0.5,1.5, 2,-0.5,1.5, 50,0,0.5);
  if(fMCcase) fOutputList->Add(fMuonPionK2);
  if(fMCcase) fOutputList->Add(fPionPionK2);
  if(fMCcase) fOutputList->Add(fMuonPionK3);
  if(fMCcase) fOutputList->Add(fPionPionK3);
  //
  TProfile *fAvgMult = new TProfile("fAvgMult","",fMbins,.5,fMbins+.5, 0,1500,"");
  fOutputList->Add(fAvgMult);
  TH2D *fAvgMultHisto2D = new TH2D("fAvgMultHisto2D","",fMbins,.5,fMbins+.5, 1000,0.5,2000.5);
  fOutputList->Add(fAvgMultHisto2D);
  TH2D *fAvgMultHisto2DV0C = new TH2D("fAvgMultHisto2DV0C","",fMbins,.5,fMbins+.5, 1000,0.5,2000.5);
  fOutputList->Add(fAvgMultHisto2DV0C);
  TH2D *fAvgMultHisto2DV0AplusC = new TH2D("fAvgMultHisto2DV0AplusC","",fMbins,.5,fMbins+.5, 1000,0.5,2000.5);
  fOutputList->Add(fAvgMultHisto2DV0AplusC);

  TH2D *fTrackChi2NDF = new TH2D("fTrackChi2NDF","",20,0,100, 100,0,10);
  fOutputList->Add(fTrackChi2NDF);
  TH2D *fTrackTPCncls = new TH2D("fTrackTPCncls","",20,0,100, 110,50,160);
  fOutputList->Add(fTrackTPCncls);



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

  TH2D *fdNchdEtaResponse = new TH2D("fdNchdEtaResponse","",15,0,15, 15,0,15);
  TH2D *fNpionTrueDist = new TH2D("fNpionTrueDist","",fMbins,.5,fMbins+.5, 3000,0.5,3000.5);
  TH2D *fNchTrueDist = new TH2D("fNchTrueDist","",fMbins,.5,fMbins+.5, 3000,0.5,3000.5);// default Nch mapping
  TH2D *fNchTrueDistCMS = new TH2D("fNchTrueDistCMS","",fMbins,.5,fMbins+.5, 3000,0.5,3000.5);// default Nch mapping
  TH2D *fNchTrueDistFullPt = new TH2D("fNchTrueDistFullPt","",fMbins,.5,fMbins+.5, 3000,0.5,3000.5);// full Pt Nch range mapping
  TH2D *fNchTrueDistPubMethod = new TH2D("fNchTrueDistPubMethod","",fMbins,.5,fMbins+.5, 3000,0.5,3000.5);// Published pp Nch mapping
  Float_t PubBins[9]={1.,12.,17.,23.,29.,35.,42.,52.,152.};
  TProfile *fAvgNchTrueDistvsPubMethodBin = new TProfile("fAvgNchTrueDistvsPubMethodBin","",8,PubBins,"");
  TProfile *fAvgRecRate = new TProfile("fAvgRecRate","",3000,0.5,3000.5, 0,3000, "");
  if(fMCcase) fOutputList->Add(fdNchdEtaResponse);
  if(fMCcase) fOutputList->Add(fNpionTrueDist);
  if(fMCcase) fOutputList->Add(fNchTrueDist);
  if(fMCcase) fOutputList->Add(fNchTrueDistCMS);
  if(fMCcase) fOutputList->Add(fNchTrueDistFullPt);
  if(fMCcase) fOutputList->Add(fNchTrueDistPubMethod);
  if(fMCcase) fOutputList->Add(fAvgRecRate);
  if(fMCcase) fOutputList->Add(fAvgNchTrueDistvsPubMethodBin);
  
  TH2D *fdCentVsNchdEta = new TH2D("fdCentVsNchdEta","",fMbins,.5,fMbins+.5, 15,0,15);
  if(fPbPbcase) fOutputList->Add(fdCentVsNchdEta);
  
  TH1D *fV0TotSignal = new TH1D("fV0TotSignal","",3000, 0,30000); 
  if(fV0Mbinning) fOutputList->Add(fV0TotSignal);
  
  TH2D *fMultBinVsCent = new TH2D("fMultBinVsCent","",fMbins,.5,fMbins+.5, 100,0,100);
  fOutputList->Add(fMultBinVsCent);

  TH1D *fExtendedQ3Histo_term1 = new TH1D("fExtendedQ3Histo_term1","",50,0,0.5);
  TH1D *fExtendedQ3Histo_term2 = new TH1D("fExtendedQ3Histo_term2","",50,0,0.5);
  TH1D *fExtendedQ3Histo_term5 = new TH1D("fExtendedQ3Histo_term5","",50,0,0.5);
  fOutputList->Add(fExtendedQ3Histo_term1);
  fOutputList->Add(fExtendedQ3Histo_term2);
  fOutputList->Add(fExtendedQ3Histo_term5);

  if(fPdensityPairCut){
    
    for(Int_t mb=0; mb<fMbins; mb++){
      if((mb < fCentBinLowLimit) || (mb > fCentBinHighLimit)) continue;
      
      for(Int_t edB=0; edB<fEDbins; edB++){
	if(edB >= fKt3bins) continue;
	
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
		
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2 = new TH2D(nameEx2->Data(),"Two Particle Distribution",20,0,1, fQbinsC2,0,fQlimitC2);
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2);
		//
		TString *nameMeanKt=new TString(nameEx2->Data());
		nameMeanKt->Append("_MeanKt");
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMeanKt = new TH1D(nameMeanKt->Data(),"Two Particle Distribution",200,0,1);
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMeanKt);
		//
		TString *nameEx2QW=new TString(nameEx2->Data());
		nameEx2QW->Append("_QW");
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW = new TH2D(nameEx2QW->Data(),"Two Particle Distribution",20,0,1, fQbinsC2,0,fQlimitC2);
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fExplicit2QW);
		TString *nameAvgP=new TString(nameEx2->Data());
		nameAvgP->Append("_AvgP");
		Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP = new TProfile2D(nameAvgP->Data(),"",10,0,1, fQbinsC2,0,fQlimitC2, 0.,1.0,"");
		fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fAvgP);
		
		// Momentum resolution histos
		if(fMCcase && sc==0){
		  TString *nameIdeal = new TString(nameEx2->Data());
		  nameIdeal->Append("_Ideal");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal = new TH2D(nameIdeal->Data(),"Two Particle Distribution",fRVALUES*kNDampValues,-0.5,fRVALUES*kNDampValues-0.5, fQbins,0,fQupperBound);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fIdeal);
		  TString *nameSmeared = new TString(nameEx2->Data());
		  nameSmeared->Append("_Smeared");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared = new TH2D(nameSmeared->Data(),"Two Particle Distribution",fRVALUES*kNDampValues,-0.5,fRVALUES*kNDampValues-0.5, fQbins,0,fQupperBound);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fSmeared);
		  //
		  TString *nameEx2MC=new TString(nameEx2->Data());
		  nameEx2MC->Append("_MCqinv");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv = new TH1D(nameEx2MC->Data(),"",fQbinsC2,0,fQlimitC2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinv);
		  TString *nameEx2MCQW=new TString(nameEx2->Data());
		  nameEx2MCQW->Append("_MCqinvQW");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW = new TH1D(nameEx2MCQW->Data(),"",fQbinsC2,0,fQlimitC2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fMCqinvQW);
		  //
		  TString *nameEx2PIDpurityDen=new TString(nameEx2->Data());
		  nameEx2PIDpurityDen->Append("_PIDpurityDen");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen = new TH2D(nameEx2PIDpurityDen->Data(),"Two Particle Distribution",20,0,1, fQbinsC2,0,fQlimitC2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityDen);
		  TString *nameEx2PIDpurityNum=new TString(nameEx2->Data());
		  nameEx2PIDpurityNum->Append("_PIDpurityNum");
		  Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum = new TH3D(nameEx2PIDpurityNum->Data(),"Two Particle Distribution",16,0.5,16.5, 20,0,1, fQbinsC2,0,fQlimitC2);
		  fOutputList->Add(Charge1[c1].Charge2[c2].SC[sc].MB[mb].EDB[edB].TwoPT[term].fPIDpurityNum);
		}
		
	      }// term_2
	    }// SC_2
	    
	    
 
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
	      
		  

		  
		  if(fPdensityPairCut){
		    TString *nameNorm=new TString(namePC3->Data());
		    nameNorm->Append("_Norm");
		    Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3 = new TH1D(nameNorm->Data(),"Norm",1,-0.5,0.5);
		    fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fNorm3);
		    //
		    if(sc<=2){
		      TString *nameQ3=new TString(namePC3->Data());
		      nameQ3->Append("_Q3");
		      Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTermsQ3 = new TH1D(nameQ3->Data(),"", 200,0,2.0);
		      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTermsQ3);
		      //
		      TString *name3DQ=new TString(namePC3->Data());
		      name3DQ->Append("_3D");
		      Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3 = new TH3D(name3DQ->Data(),"", fQbins,0,fQupperBound, fQbins,0,fQupperBound, fQbins,0,fQupperBound);
		      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fTerms3);
		      //
		      TString *nameMeanKt=new TString(namePC3->Data());
		      nameMeanKt->Append("_MeanKt");
		      Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fMeanKt = new TH1D(nameMeanKt->Data(),"", 200,0,1);
		      fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fMeanKt);
		      
		      if(sc==0 && fMCcase==kTRUE){
			TString *name3DMomResIdeal=new TString(namePC3->Data());
			name3DMomResIdeal->Append("_Ideal");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal = new TH1D(name3DMomResIdeal->Data(),"", 200,0,2.0);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fIdeal);
			TString *name3DMomResSmeared=new TString(namePC3->Data());
			name3DMomResSmeared->Append("_Smeared");
			Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared = new TH1D(name3DMomResSmeared->Data(),"", 200,0,2.0);
			fOutputList->Add(Charge1[c1].Charge2[c2].Charge3[c3].SC[sc].MB[mb].EDB[edB].ThreePT[term].fSmeared);
		      }// MCcase
		      
		      
		    }// sc exclusion
		  }// PdensityPairCut
		}// term_3
	      }// SC_3
	    }//c3
	  }//c2
	}//c1
      }// ED
    }// mbin
  }// Pdensity Method

  
    
  TH1D *frstar4VectDist = new TH1D("frstar4VectDist","",10000,0,100);
  fOutputList->Add(frstar4VectDist);
  
  
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
void AliThreePionRadii::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  fEventCounter++;
  if(fEventCounter%1000==0) cout<<"===========  Event # "<<fEventCounter<<"  ==========="<<endl;

  if(!fAODcase) {cout<<"ESDs not supported"<<endl; return;}
  
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}
  
  
  // Trigger Cut
  if(fPbPbcase){
    if(fAOD->GetRunNumber() >= 136851 && fAOD->GetRunNumber() <= 139517){// 10h data
      Bool_t isSelected1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      if(!isSelected1 && !fMCcase) {return;}
    }
    if(fAOD->GetRunNumber() >= 167693 && fAOD->GetRunNumber() <= 170593){// 11h data
      Bool_t isSelected1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      Bool_t isSelected2 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
      Bool_t isSelected3 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
      if(!isSelected1 && !isSelected2 && !isSelected3 && !fMCcase) {return;}
    }
  }else{// pp and pPb
    Bool_t isSelected[4]={kFALSE};
    isSelected[0] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    isSelected[1] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
    isSelected[2] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    isSelected[3] = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kHighMult);
    if(!isSelected[fTriggerType] && !fMCcase) return;
  }
  
  
  ///////////////////////////////////////////////////////////
  const AliAODVertex *primaryVertexAOD;
  AliCentrality *centrality;// for AODs and ESDs

 
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  
  TClonesArray *mcArray = 0x0;
  Int_t mcNch=0, mcNchCMS=0, mcNchFullPt=0, mcNchPubMethod=0;
  Int_t mcNpion=0;
  if(fMCcase){
    if(fAODcase){ 
      mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
      if(!mcArray || mcArray->GetEntriesFast() >= kMCarrayLimit){
	cout<<"No MC particle branch found or Array too large!!"<<endl;
	return;
      }
      
      // Count true Nch at mid-rapidity
      for(Int_t mctrackN=0; mctrackN<mcArray->GetEntriesFast(); mctrackN++){
	AliAODMCParticle *mcParticle = (AliAODMCParticle*)mcArray->At(mctrackN);
	if(!mcParticle) continue;
	if(mcParticle->Charge()!=-3 && mcParticle->Charge()!=+3) continue;// x3 by convention
	//if(!mcParticle->IsPrimary()) continue;// superfluous when IsPhysicalPrimary() is used
	if(!mcParticle->IsPhysicalPrimary()) continue;
	//
	if(fabs(mcParticle->Eta())<=0.5) mcNchPubMethod++;// Published pp binning
	if(fabs(mcParticle->Eta())<=0.8) mcNchFullPt++;// My binning in full Pt range
	
	if(mcParticle->Pt() < 0.16 || mcParticle->Pt() > 1.0) continue;
	
	//
	if(mcParticle->P() < 1.0) {
	  if(fabs(mcParticle->Eta())<=0.8) {
	    mcNch++;// My binning in my pt range
	    if(abs(mcParticle->GetPdgCode())==211) mcNpion++;
	  }
	}
	// p-Pb CMS boost counting
	Double_t newPz = mcParticle->Pz()*cosh(0.465) - mcParticle->E()*sinh(0.465);
	Double_t newP = sqrt(pow(mcParticle->Pt(),2) + pow(newPz,2));
	if(newP < 1.0){
	  Double_t newEta = 0.5 * log( (newP+newPz) / (newP-newPz));
	  if(TMath::Abs(newEta)<=0.8) {
	    mcNchCMS++;
	  }
	}
      }
    }
  }// fMCcase
  
  UInt_t status=0;
  Int_t positiveTracks=0, negativeTracks=0;
  Int_t myTracks=0, pionCount=0, kaonCount=0, protonCount=0;
  Int_t FBTracks=0, AODTracks=0;

  Double_t vertex[3]={0};
  Int_t zbin=0;
  Double_t zstep=2*10/Double_t(fZvertexBins), zstart=-10.;
  /////////////////////////////////////////////////
  
  Float_t centralityPercentile=0;
  //Float_t cStep=5.0, cStart=0;
  Int_t trackletMult = 0;

  if(fAODcase){// AOD case
    
    if(fPbPbcase){
      centrality = fAOD->GetCentrality();
      centralityPercentile = centrality->GetCentralityPercentile("V0M");
      if(centralityPercentile == 0) {/*cout<<"Centrality = 0, skipping event"<<endl;*/ return;}
      //if((centralityPercentile < 5*fCentBinLowLimit) || (centralityPercentile>= 5*(fCentBinHighLimit+1))) {/*cout<<"Centrality out of Range.  Skipping Event"<<endl;*/ return;}
      cout<<"Centrality % = "<<centralityPercentile<<endl;
    }else{
      //cout<<"AOD multiplicity = "<<fAOD->GetNumberOfTracks()<<endl;
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
    
    if(fabs(vertex[2]) > 10) {/*cout<<"Zvertex Out of Range. Skip Event"<<endl;*/ return;} // Z-Vertex Cut 
    ((TH3F*)fOutputList->FindObject("fVertexDist"))->Fill(vertex[0], vertex[1], vertex[2]);
    
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!aodtrack) AliFatal("Not a standard AOD");
      if (!aodtrack) continue;
      AODTracks++;
      if(!aodtrack->TestFilterBit(BIT(fFilterBit))) continue;// AOD filterBit cut
      FBTracks++;
    }
    ((TH1F*)fOutputList->FindObject("fMultDist2"))->Fill(FBTracks);

    //if(fAOD->IsPileupFromSPD()) {/*cout<<"PileUpEvent. Skip Event"<<endl;*/ return;} // Old Pile-up cut
    if(primaryVertexAOD->GetNContributors() < 1) {/*cout<<"Bad Vertex. Skip Event"<<endl;*/ return;}
   
    ((TH1F*)fOutputList->FindObject("fMultDist3"))->Fill(FBTracks);
 
    fBfield = fAOD->GetMagneticField();
    
    for(Int_t i=0; i<fZvertexBins; i++){
      if( (vertex[2] >= zstart+i*zstep) && (vertex[2] < zstart+(i+1)*zstep) ){
	zbin=i;
	break;
      }
    }
    
    AliAODTracklets *tracklets = (AliAODTracklets*)fAOD->GetTracklets();
    for(Int_t trackletN=0; trackletN<tracklets->GetNumberOfTracklets(); trackletN++){
      if(tracklets->GetTheta(trackletN) > 1.0904 && tracklets->GetTheta(trackletN) < 2.0512) trackletMult++;// |eta|<0.5 tracklets
    }
   
    /////////////////////////////
    // Create Shuffled index list
    Int_t randomIndex[fAOD->GetNumberOfTracks()];
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) randomIndex[i]=i;
    Shuffle(randomIndex,0,fAOD->GetNumberOfTracks()-1);
    /////////////////////////////
  
    // Track loop
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(randomIndex[i]));
      if(!aodtrack) AliFatal("Not a standard AOD");
      if (!aodtrack) continue;
      if(myTracks >= fMultLimit) {cout<<"More tracks than Track Limit"<<endl; return;}
      
      status=aodtrack->GetStatus();
      
      if(!aodtrack->TestFilterBit(BIT(7))) continue;// AOD filterBit cut
      if(aodtrack->GetTPCNcls() < 70) continue;// TPC nCluster cut
      
      // FilterBit Overlap Check
      if(fFilterBit != 7){
	Bool_t goodTrackOtherFB = kFALSE;
	if(fMCcase && fAOD->GetRunNumber()<=126437) goodTrackOtherFB=kTRUE;// FB7 to FB5 mapping in 10f6a MC does not work
	
	for (Int_t j = 0; j < fAOD->GetNumberOfTracks(); j++) {
	  AliAODTrack* aodtrack2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(randomIndex[j]));
	  if(!aodtrack2) AliFatal("Not a standard AOD");
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
	  AliAODTrack* aodTrack2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));
	  if(!aodTrack2) AliFatal("Not a standard AOD");
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

      //cout<<nSigmaTPC[2]<<endl;
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
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;
      if(fTempStruct[myTracks].fPion && fTempStruct[myTracks].fKaon && fTempStruct[myTracks].fProton) continue;
      ////////////////////////
      if(fTempStruct[myTracks].fProton && fTempStruct[myTracks].fMom < 0.25) continue;//extra cut for protons

      if(!fTempStruct[myTracks].fPion) continue;// only pions
          
     


      if(fTempStruct[myTracks].fCharge==+1) {
	((TH2F*)fOutputList->FindObject("fDCAxyDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistPlus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }else {
	((TH2F*)fOutputList->FindObject("fDCAxyDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[0]);
	((TH2F*)fOutputList->FindObject("fDCAzDistMinus"))->Fill(fTempStruct[myTracks].fPt, dca2[1]);
      }
      
      ((TH3F*)fOutputList->FindObject("fPhiPtDist"))->Fill(aodtrack->Charge(), aodtrack->Phi(), aodtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fPtEtaDist"))->Fill(aodtrack->Charge(), aodtrack->Pt(), aodtrack->Eta());

           
    
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
      
     
      ((TH2D*)fOutputList->FindObject("fTrackChi2NDF"))->Fill(centralityPercentile, aodtrack->Chi2perNDF());
      ((TH2D*)fOutputList->FindObject("fTrackTPCncls"))->Fill(centralityPercentile, aodtrack->GetTPCncls());
      if(aodtrack->Chi2perNDF() > fMaxChi2NDF) continue;
      if(aodtrack->GetTPCncls() < fMinTPCncls) continue;
      

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
    ((TH1F*)fOutputList->FindObject("fMultDist4"))->Fill(myTracks);
  }
 
 
  //cout<<"There are "<<myTracks<<"  myTracks"<<endl;
  //cout<<"pionCount = "<<pionCount<<"   kaonCount = "<<kaonCount<<"   protonCount = "<<protonCount<<endl;

  /////////////////////////////////////////
  // Pion Multiplicity Cut (To ensure all Correlation orders are present in each event)
  if(myTracks < 3) {/*cout<<"Less than 3 tracks. Skipping Event."<<endl;*/ return;}
  /////////////////////////////////////////
 

  ////////////////////////////////
  ///////////////////////////////
  // Mbin determination
  //
  fMbin=-1;
  for(Int_t i=0; i<fCentBins; i++){
    if( pionCount >= fMultLimits[i] && pionCount < fMultLimits[i+1]) {fMbin = fCentBins-i-1; break;}
  }
  
  
  fFSIindex=0;
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
  
  if(fMCcase){// FSI binning for MC 
    if(fRMax>=10) fFSIindex = 0;
    else if(fRMax>=9) fFSIindex = 1;
    else if(fRMax>=8) fFSIindex = 2;
    else if(fRMax>=7) fFSIindex = 3;
    else if(fRMax>=6) fFSIindex = 4;
    else if(fRMax>=5) fFSIindex = 5;
    else if(fRMax>=4) fFSIindex = 6;
    else if(fRMax>=3) fFSIindex = 7;
    else if(fRMax>=2) fFSIindex = 8;
    else fFSIindex = 9;
  }

  if(fV0Mbinning){
    Bool_t useV0=kFALSE;
    if(fPbPbcase) useV0=kTRUE;
    if(!fPbPbcase && fAOD->GetRunNumber() >= 195344 && fAOD->GetRunNumber() <= 195677) useV0=kTRUE;
    if(useV0){
      AliAODVZERO *vZero = fAOD->GetVZEROData();
      Float_t vZeroAmp = vZero->GetMTotV0A();
      centrality = fAOD->GetCentrality();
      centralityPercentile = centrality->GetCentralityPercentile("V0M");
      for(Int_t i=0; i<fCentBins; i++){
	if(vZeroAmp/4.4 >= fMultLimits[i] && vZeroAmp/4.4 < fMultLimits[i+1]) {fMbin = fCentBins-i-1; break;}
      }
      ((TH1D*)fOutputList->FindObject("fV0TotSignal"))->Fill(vZeroAmp);
      //cout<<centralityPercentile<<"  "<<vZeroAmp<<"  "<<fMbin<<endl;
      //
      Int_t fMbinV0C=-1;
      vZeroAmp = vZero->GetMTotV0C();
      for(Int_t i=0; i<fCentBins; i++){
	if(vZeroAmp/4.4 >= fMultLimits[i] && vZeroAmp/4.4 < fMultLimits[i+1]) {fMbinV0C = fCentBins-i-1; break;}
      }
      //
      Int_t fMbinV0AplusC=-1;
      vZeroAmp = vZero->GetMTotV0A() + vZero->GetMTotV0C();
      for(Int_t i=0; i<fCentBins; i++){
	if(vZeroAmp/4.4 >= fMultLimits[i] && vZeroAmp/4.4 < fMultLimits[i+1]) {fMbinV0AplusC = fCentBins-i-1; break;}
      }
      ((TH2D*)fOutputList->FindObject("fAvgMultHisto2DV0C"))->Fill(fMbinV0C+1., pionCount);
      ((TH2D*)fOutputList->FindObject("fAvgMultHisto2DV0AplusC"))->Fill(fMbinV0AplusC+1., pionCount);
    }
  }
  
  if(fMbin==-1) {cout<<pionCount<<"  Bad Mbin+++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl; return;}

  ((TH1F*)fOutputList->FindObject("fEvents1"))->Fill(fMbin+1);
  ((TProfile*)fOutputList->FindObject("fAvgMult"))->Fill(fMbin+1., pionCount);
  ((TH2D*)fOutputList->FindObject("fAvgMultHisto2D"))->Fill(fMbin+1., pionCount);
  if(fMCcase){
    ((TH2D*)fOutputList->FindObject("fdNchdEtaResponse"))->Fill(pow(trackletMult,1/3.), pow(mcNch,1/3.));
    ((TH2D*)fOutputList->FindObject("fNpionTrueDist"))->Fill(fMbin+1., mcNpion);
    ((TH2D*)fOutputList->FindObject("fNchTrueDist"))->Fill(fMbin+1., mcNch);// default Nch mapping
    ((TH2D*)fOutputList->FindObject("fNchTrueDistCMS"))->Fill(fMbin+1., mcNchCMS);// p-Pb CMS counting
    ((TH2D*)fOutputList->FindObject("fNchTrueDistFullPt"))->Fill(fMbin+1., mcNchFullPt);// full Pt Nch range mapping
    ((TH2D*)fOutputList->FindObject("fNchTrueDistPubMethod"))->Fill(fMbin+1., mcNchPubMethod);// Published pp Method Nch mapping
    ((TProfile*)fOutputList->FindObject("fAvgNchTrueDistvsPubMethodBin"))->Fill(mcNchPubMethod, mcNch);// Mapping of Published bins to default Nch bins
    ((TProfile*)fOutputList->FindObject("fAvgRecRate"))->Fill(mcNpion, pionCount);
  }
  if(fPbPbcase){
    ((TH2D*)fOutputList->FindObject("fdCentVsNchdEta"))->Fill(fMbin+1, pow(trackletMult,1/3.));
    centrality = fAOD->GetCentrality();
    centralityPercentile = centrality->GetCentralityPercentile("V0M");
    ((TH2D*)fOutputList->FindObject("fMultBinVsCent"))->Fill(fMbin+1, centralityPercentile);
  }

  // Mult cut
  if(fMbin < fCentBinLowLimit || fMbin > fCentBinHighLimit) {cout<<"Mult out of range"<<endl; return;}
  
  //////////////////////////////////////////////////
  fEDbin=0;// Extra Dimension bin (Kt, (Kt-Psi),....)
  //////////////////////////////////////////////////


  
  //return;// un-comment for a run to calculate Nrec to Nch Mapping 
  // to test the eta dependence of radii
  /*Int_t firstTrackCount=myTracks;
  Int_t newTrackCount=0;
  myTracks=0; pionCount=0; kaonCount=0; protonCount=0;// reset track counters
  for(Int_t newTracks=0; newTracks<firstTrackCount; newTracks++){
    
    if(fTempStruct[newTracks].fEta > -0.4) continue;
    
    fTempStruct[newTrackCount]=fTempStruct[newTracks];
  
    newTrackCount++;
    pionCount++;
  }
  myTracks=newTrackCount;// re-assign main counter
  */
  
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
  Float_t qout=0, qside=0, qlong=0;
  Float_t qoutMC=0, qsideMC=0, qlongMC=0;
  Float_t firstQ=0, secondQ=0, thirdQ=0;
  Float_t firstQMC=0, secondQMC=0, thirdQMC=0;
  Float_t transK12=0, transK3=0;
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
  Int_t index1=0, index2=0, index3=0;
  Float_t qinv12MC=0, qinv13MC=0, qinv23MC=0;
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
  	      
 

    /////////////////////////////////////////////////////////////////
    // extended range Q3 baseline
    /*for(Int_t iter=0; iter<3; iter++){
      for (Int_t i=0; i<myTracks; i++) {
	
	Int_t en2=0;
	if(iter==2) en2=1;
	Int_t start2=i+1;
	if(en2!=0) start2=0;
	// 2nd particle
	for (Int_t j=start2; j<(fEvt+en2)->fNtracks; j++) {
	  if((fEvt)->fTracks[i].fCharge != (fEvt+en2)->fTracks[j].fCharge) continue;
	  key1 = (fEvt)->fTracks[i].fKey;
	  key2 = (fEvt+en2)->fTracks[j].fKey;
	  Short_t fillIndex2 = FillIndex2part(key1+key2);
	  pVect1[0]=(fEvt)->fTracks[i].fEaccepted; pVect2[0]=(fEvt+en2)->fTracks[j].fEaccepted;
	  pVect1[1]=(fEvt)->fTracks[i].fP[0];      pVect2[1]=(fEvt+en2)->fTracks[j].fP[0];
	  pVect1[2]=(fEvt)->fTracks[i].fP[1];      pVect2[2]=(fEvt+en2)->fTracks[j].fP[1];
	  pVect1[3]=(fEvt)->fTracks[i].fP[2];      pVect2[3]=(fEvt+en2)->fTracks[j].fP[2];
	  qinv12 = GetQinv(fillIndex2, pVect1, pVect2);
	  
	  if(qinv12>0.5) continue;
	  if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) continue;
	  
	  Int_t en3=0;
	  if(iter==1) en3=1;
	  if(iter==2) en3=2;
	  Int_t start3=j+1;
	  if(iter>0) start3=0;
	  // 3nd particle
	  for (Int_t k=start3; k<(fEvt+en3)->fNtracks; k++) {
	    if((fEvt)->fTracks[i].fCharge != (fEvt+en3)->fTracks[k].fCharge) continue;
	    pVect3[0]=(fEvt+en3)->fTracks[k].fEaccepted;
	    pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
	    pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
	    pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
	    qinv13 = GetQinv(fillIndex2, pVect1, pVect3);
	    if(qinv13>0.5) continue;
	    qinv23 = GetQinv(fillIndex2, pVect2, pVect3);
	    if(qinv23>0.5) continue;
	    if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en3)->fTracks[k]))) continue;
	    if(!AcceptPair(&((fEvt+en2)->fTracks[j]), &((fEvt+en3)->fTracks[k]))) continue;

	    q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));

	    if(iter==0) ((TH1D*)fOutputList->FindObject("fExtendedQ3Histo_term1"))->Fill(q3);
	    if(iter==1) ((TH1D*)fOutputList->FindObject("fExtendedQ3Histo_term2"))->Fill(q3);
	    if(iter==2) ((TH1D*)fOutputList->FindObject("fExtendedQ3Histo_term5"))->Fill(q3);
	    
	  }
	}
      }
    }
    */
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


	//

	///////////////////////////////
	ch1 = Int_t(((fEvt)->fTracks[i].fCharge + 1)/2.);
	ch2 = Int_t(((fEvt+en2)->fTracks[j].fCharge + 1)/2.);
	SetFillBins2(fillIndex2, key1, key2, ch1, ch2, bin1, bin2);
	
	if(fMCcase && ch1==ch2 && fMbin==0 && qinv12<0.2){
	  //////////////////////////
	  // pad-row method testing
	  Float_t coeff = (5)*0.2*(0.18/1.2);// 5 to evaluate at 1.0m in TPC
	  Float_t phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	  if(phi1 > 2*PI) phi1 -= 2*PI;
	  if(phi1 < 0) phi1 += 2*PI;
	  Float_t phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	  if(phi2 > 2*PI) phi2 -= 2*PI;
	  if(phi2 < 0) phi2 += 2*PI;
	  Float_t deltaphi = phi1 - phi2;
	  if(deltaphi > PI) deltaphi -= PI;
	  if(deltaphi < -PI) deltaphi += PI;
	  
	  Int_t ncl1 = (fEvt)->fTracks[i].fClusterMap.GetNbits();
	  Int_t ncl2 = (fEvt+en2)->fTracks[j].fClusterMap.GetNbits();
	  Float_t sumCls = 0; Float_t sumSha = 0; Float_t sumQ = 0;
	  Double_t shfrac = 0; //Double_t qfactor = 0;
	  for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
	    if ((fEvt)->fTracks[i].fClusterMap.TestBitNumber(imap) && (fEvt+en2)->fTracks[j].fClusterMap.TestBitNumber(imap)) {// Both clusters
	      if ((fEvt)->fTracks[i].fSharedMap.TestBitNumber(imap) && (fEvt+en2)->fTracks[j].fSharedMap.TestBitNumber(imap)) { // Shared
		sumQ++;
		sumCls+=2;
		sumSha+=2;}
	      else {sumQ--; sumCls+=2;}
	    }
	    else if ((fEvt)->fTracks[i].fClusterMap.TestBitNumber(imap) || (fEvt+en2)->fTracks[j].fClusterMap.TestBitNumber(imap)) {// Non shared
	      sumQ++;
	      sumCls++;}
	  }
	  if (sumCls>0) {
	    //qfactor = sumQ*1.0/sumCls;
	    shfrac = sumSha*1.0/sumCls;
	  }
	  if(fabs(deltaphi)<0.07 && fabs((fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta) < 0.03){
	    ((TH3D*)fOutputList->FindObject("fPairsPadRowNum"))->Fill(transK12, shfrac, qinv12);
	  }
	  
	  for(Int_t rstep=0; rstep<10; rstep++){
	    coeff = (rstep)*0.2*(0.18/1.2);
	    phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	    if(phi1 > 2*PI) phi1 -= 2*PI;
	    if(phi1 < 0) phi1 += 2*PI;
	    phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	    if(phi2 > 2*PI) phi2 -= 2*PI;
	    if(phi2 < 0) phi2 += 2*PI;
	    deltaphi = phi1 - phi2;
	    if(deltaphi > PI) deltaphi -= PI;
	    if(deltaphi < -PI) deltaphi += PI;

	    if(fabs((fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta) < 0.03){
	      ((TH3F*)fOutputList->FindObject("fPairsShareFracDPhiNum"))->Fill(rstep, shfrac, deltaphi);
	    }
	    //if(shfrac < 0.05){
	    ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiNum"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
	    //}
	  }
	  
	  
	}// MCcase and pair selection
	
	// Pair Splitting/Merging cut
	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting/merging cut)
	if(ch1 == ch2 && !fGeneratorOnly){
	  if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) {
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
	    
	    
	   
	    mcParticle1 = (AliAODMCParticle*)mcArray->At(abs((fEvt)->fTracks[i].fLabel));
	    mcParticle2 = (AliAODMCParticle*)mcArray->At(abs((fEvt+en2)->fTracks[j].fLabel));
	    
	   
	    // secondary contamination
	    if(abs(mcParticle1->GetPdgCode())==211 && abs(mcParticle2->GetPdgCode())==211){
	      if(ch1==ch2) {
		((TH3D*)fOutputList->FindObject("fAllSCPionPairs"))->Fill(fMbin+1, transK12, qinv12);
		if(!mcParticle1->IsSecondaryFromWeakDecay() && !mcParticle2->IsSecondaryFromWeakDecay()) {
		  ((TH3D*)fOutputList->FindObject("fPrimarySCPionPairs"))->Fill(fMbin+1, transK12, qinv12);
		}	      
	      }else{
		((TH3D*)fOutputList->FindObject("fAllMCPionPairs"))->Fill(fMbin+1, transK12, qinv12);
		if(!mcParticle1->IsSecondaryFromWeakDecay() && !mcParticle2->IsSecondaryFromWeakDecay()) {
		  ((TH3D*)fOutputList->FindObject("fPrimaryMCPionPairs"))->Fill(fMbin+1, transK12, qinv12);
		}
	      }
	    }

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
	    
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMCqinv->Fill(qinv12MC, MCWeight(ch1,ch2,rForQW,10,qinv12MC));// was 4,5
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMCqinvQW->Fill(qinv12MC, qinv12MC*MCWeight(ch1,ch2,rForQW,10,qinv12MC));// was 4,5
	    // pion purity
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fPIDpurityDen->Fill(transK12, qinv12);
	    Int_t SCNumber = 1;
	    
	    if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==22) SCNumber=1;// e-e
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==24) SCNumber=2;// e-mu
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==222) SCNumber=3;// e-pi
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==332) SCNumber=4;// e-k
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==2223) SCNumber=5;// e-p
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==26) SCNumber=6;// mu-mu
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==224) SCNumber=7;// mu-pi
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==334) SCNumber=8;// mu-k
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==2225) SCNumber=9;// mu-p
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==422) SCNumber=10;// pi-pi
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==532) SCNumber=11;// pi-k
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==2423) SCNumber=12;// pi-p
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==642) SCNumber=13;// k-k
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==2533) SCNumber=14;// k-p
	    else if((abs(mcParticle1->GetPdgCode()) + abs(mcParticle2->GetPdgCode()))==4424) SCNumber=15;// p-p
	    else SCNumber=16;
	    
	    Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fPIDpurityNum->Fill(SCNumber, transK12, qinv12);
	    
	   
	    ///////////////////////
	    // muon contamination
	    if(qinv12 < fQcut[0] && ((fEvt)->fTracks[i].fLabel != (fEvt+en2)->fTracks[j].fLabel)){
	      if(abs(mcParticle1->GetPdgCode())==13 || abs(mcParticle2->GetPdgCode())==13){// muon check
		Float_t Pparent1[4]={pVect1MC[0],pVect1MC[1],pVect1MC[2],pVect1MC[3]}; 
		Float_t Pparent2[4]={pVect2MC[0],pVect2MC[1],pVect2MC[2],pVect2MC[3]};
		Bool_t pionParent1=kFALSE, pionParent2=kFALSE;
		if(abs(mcParticle1->GetPdgCode())==13) {
		  AliAODMCParticle *parent1=(AliAODMCParticle*)mcArray->At(mcParticle1->GetMother());
		  if(abs(parent1->GetPdgCode())==211) {
		    pionParent1=kTRUE;
		    Pparent1[1] = parent1->Px(); Pparent1[2] = parent1->Py(); Pparent1[3] = parent1->Pz();
		    Pparent1[0] = sqrt(pow(Pparent1[1],2)+pow(Pparent1[2],2)+pow(Pparent1[3],2)+pow(fTrueMassPi,2));
		  }
		}
		// 
		if(abs(mcParticle2->GetPdgCode())==13) {
		  AliAODMCParticle *parent2=(AliAODMCParticle*)mcArray->At(mcParticle2->GetMother());
		  if(abs(parent2->GetPdgCode())==211) {
		    pionParent2=kTRUE;
		    Pparent2[1] = parent2->Px(); Pparent2[2] = parent2->Py(); Pparent2[3] = parent2->Pz();
		    Pparent2[0] = sqrt(pow(Pparent2[1],2)+pow(Pparent2[2],2)+pow(Pparent2[3],2)+pow(fTrueMassPi,2));
		  }
		}
		
		Float_t parentQinv12 = GetQinv(0, Pparent1, Pparent2);
		Float_t WInput = 1.0, muonPionK12=1.0, pionPionK12=1.0;
		if(parentQinv12 > 0.001 && parentQinv12 < 0.3) {
		  WInput = MCWeight(ch1,ch2, fRMax, 25, parentQinv12);
		  muonPionK12 = FSICorrelation2(ch1, ch2, qinv12MC);
		  pionPionK12 = FSICorrelation2(ch1, ch2, parentQinv12);
		}
		Int_t ChComb=0;
		if(ch1 != ch2) ChComb=1;
		if(pionParent1 || pionParent2){
		  ((TH3D*)fOutputList->FindObject("fMuonContamSmearedNum2"))->Fill(ChComb, transK12, qinv12MC, WInput);
		  ((TH3D*)fOutputList->FindObject("fMuonContamSmearedDen2"))->Fill(ChComb, transK12, qinv12MC);
		  ((TH3D*)fOutputList->FindObject("fMuonContamIdealNum2"))->Fill(ChComb, transK12, parentQinv12, WInput);
		  ((TH3D*)fOutputList->FindObject("fMuonContamIdealDen2"))->Fill(ChComb, transK12, parentQinv12);
		  ((TH3D*)fOutputList->FindObject("fMuonPionDeltaQinv"))->Fill(ChComb, transK12, qinv12MC-parentQinv12);
		  ((TH3D*)fOutputList->FindObject("fMuonPionK2"))->Fill(ChComb, transK12, qinv12MC, muonPionK12);
		  ((TH3D*)fOutputList->FindObject("fPionPionK2"))->Fill(ChComb, transK12, parentQinv12, pionPionK12);
		}
		////////////////////////////////////
		// 3rd particle
		Int_t en3=0;
		for (Int_t k=j+1; k<(fEvt+en3)->fNtracks; k++) {
		  pVect3[0]=(fEvt+en3)->fTracks[k].fEaccepted;
		  pVect3[1]=(fEvt+en3)->fTracks[k].fP[0];
		  pVect3[2]=(fEvt+en3)->fTracks[k].fP[1];
		  pVect3[3]=(fEvt+en3)->fTracks[k].fP[2];
		  //
		  qinv13 = GetQinv(0, pVect1, pVect3);
		  qinv23 = GetQinv(0, pVect2, pVect3);
		  if(qinv13 > fQcut[0] || qinv23 > fQcut[0]) continue;
		  
		  if(qinv13 < fQLowerCut || qinv23 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting/merging cut)
		  if(ch1 == ch3 && !fGeneratorOnly){
		    if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en3)->fTracks[k]))) {
		      continue;
		    }
		  }
		  if(ch2 == ch3 && !fGeneratorOnly){
		    if(!AcceptPair(&((fEvt+en2)->fTracks[j]), &((fEvt+en3)->fTracks[k]))) {
		      continue;
		    }
		  }
		  
		  if((fEvt+en3)->fTracks[k].fLabel == (fEvt+en2)->fTracks[j].fLabel) continue;
		  if((fEvt+en3)->fTracks[k].fLabel == (fEvt)->fTracks[i].fLabel) continue;
		  
		  if((fEvt+en3)->fTracks[k].fLabel < (fEvt+en3)->fMCarraySize){
		    AliAODMCParticle *mcParticle3 = (AliAODMCParticle*)mcArray->At(abs((fEvt+en3)->fTracks[k].fLabel));
		    
		    ch3 = Int_t(((fEvt+en3)->fTracks[k].fCharge + 1)/2.);
		    pVect3MC[0]=sqrt(pow((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		    pVect3MC[1]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPx;
		    pVect3MC[2]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPy;
		    pVect3MC[3]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPz;
		    qinv13MC = GetQinv(0, pVect1MC, pVect3MC);
		    qinv23MC = GetQinv(0, pVect2MC, pVect3MC);
		    
		    q3MC = sqrt(pow(qinv12MC,2)+pow(qinv13MC,2)+pow(qinv23MC,2));
		    transK3 = sqrt( pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
		    Int_t K3index=0;
		    if(transK3>0.3) K3index=1;

		    Float_t Pparent3[4]={pVect3MC[0],pVect3MC[1],pVect3MC[2],pVect3MC[3]}; 
		    Bool_t pionParent3=kFALSE;
		    if(abs(mcParticle3->GetPdgCode())==13){// muon check
		      AliAODMCParticle *parent3=(AliAODMCParticle*)mcArray->At(mcParticle3->GetMother());
		      if(abs(parent3->GetPdgCode())==211) {
			pionParent3=kTRUE;
			Pparent3[1] = parent3->Px(); Pparent3[2] = parent3->Py(); Pparent3[3] = parent3->Pz();
			Pparent3[0] = sqrt(pow(Pparent3[1],2)+pow(Pparent3[2],2)+pow(Pparent3[3],2)+pow(fTrueMassPi,2));
		      }
		    }
		    
		    Float_t parentQinv13 = GetQinv(0, Pparent1, Pparent3);
		    Float_t parentQinv23 = GetQinv(0, Pparent2, Pparent3);
		    Float_t parentQ3 = sqrt(pow(parentQinv12,2) + pow(parentQinv13,2) + pow(parentQinv23,2));
		    if(parentQ3 >= 0.5) continue;
		    if(parentQinv12 < 0.001) continue;
		    if(parentQinv13 < 0.001) continue;
		    if(parentQinv23 < 0.001) continue;
		    
		    if(!pionParent1 && !pionParent2 && !pionParent3) continue;// want at least one pion-->muon
		    
		    Int_t ChCombtriplet=0;
		    if(ch1!=ch2 || ch1!=ch3 || ch2!=ch3) ChCombtriplet=1;
		    Float_t WInput3=1.0;
		    Float_t muonPionK13=1.0, muonPionK23=1.0;
		    Float_t pionPionK13=1.0, pionPionK23=1.0;
		    muonPionK13 = FSICorrelation2(ch1, ch3, qinv13MC);
		    pionPionK13 = FSICorrelation2(ch1, ch3, parentQinv13);
		    muonPionK23 = FSICorrelation2(ch2, ch3, qinv23MC);
		    pionPionK23 = FSICorrelation2(ch2, ch3, parentQinv23);
		    if(ChCombtriplet==0) WInput3 = MCWeight3D(kTRUE, 1, 25, parentQinv12, parentQinv13, parentQinv23);
		    else{
		      if(ch1==ch2) WInput3 = MCWeight3D(kFALSE, 1, 25, parentQinv12, parentQinv13, parentQinv23);
		      else if(ch1==ch3) WInput3 = MCWeight3D(kFALSE, 1, 25, parentQinv13, parentQinv12, parentQinv23);
		      else WInput3 = MCWeight3D(kFALSE, 1, 25, parentQinv23, parentQinv12, parentQinv13);
		    }
		    if(WInput3>0 && WInput3<10.) {
		      ((TH3D*)fOutputList->FindObject("fMuonContamSmearedNum3"))->Fill(ChCombtriplet, K3index, q3MC, WInput3);
		      ((TH3D*)fOutputList->FindObject("fMuonContamSmearedDen3"))->Fill(ChCombtriplet, K3index, q3MC);
		      ((TH3D*)fOutputList->FindObject("fMuonContamIdealNum3"))->Fill(ChCombtriplet, K3index, parentQ3, WInput3);
		      ((TH3D*)fOutputList->FindObject("fMuonContamIdealDen3"))->Fill(ChCombtriplet, K3index, parentQ3);
		      ((TH3D*)fOutputList->FindObject("fMuonPionK3"))->Fill(ChCombtriplet, K3index, q3MC, muonPionK12*muonPionK13*muonPionK23);
		      ((TH3D*)fOutputList->FindObject("fPionPionK3"))->Fill(ChCombtriplet, K3index, parentQ3, pionPionK12*pionPionK13*pionPionK23);
		    }
		  }//label check of 3
		}// 3rd particle
	      }// muon code check of 1 and 2
	    }// qinv12 cut
	  }// label check 2
	}// MC case
	
	//////////////////////////////////////////
	// 2-particle term
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2->Fill(transK12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2QW->Fill(transK12, qinv12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fAvgP->Fill(transK12, qinv12, (fEvt)->fTracks[i].fMom);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fAvgP->Fill(transK12, qinv12, (fEvt+en2)->fTracks[j].fMom);
	//
	if(qinv12<fQupperBound) Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMeanKt->Fill(transK12);
	
	//////////////////////////////////////////
	
	

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
	   
	    /////////////////////////////////////////////////////
	    if(qinv12 <= fQcut[qCutBin]) {// 3-particle MRC
	      
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
		  q3 = sqrt(pow(qinv12,2)+pow(qinv13,2)+pow(qinv23,2));
		  
		  if(qinv13 > fQcut[qCutBin] || qinv23 > fQcut[qCutBin]) continue;
		  
		  
		  pVect3MC[0]=sqrt(pow((fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPtot,2)+pow(fTrueMassPi,2)); 
		  pVect3MC[1]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPx;
		  pVect3MC[2]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPy;
		  pVect3MC[3]=(fEvt+en3)->fMCtracks[abs((fEvt+en3)->fTracks[k].fLabel)].fPz;
		  qinv13MC = GetQinv(0, pVect1MC, pVect3MC);
		  qinv23MC = GetQinv(0, pVect2MC, pVect3MC);
		  
		  
		  q3MC = sqrt(pow(qinv12MC,2)+pow(qinv13MC,2)+pow(qinv23MC,2));
		  transK3 = sqrt( pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
		  
		  if(qinv12 < fQLowerCut) continue;
		  if(qinv13 < fQLowerCut) continue;
		  if(qinv23 < fQLowerCut) continue;
		  if(ch1 == ch2){
		    if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) continue;
		  }
		  if(ch1 == ch3){
		    if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en3)->fTracks[k]))) continue;
		  }
		  if(ch2 == ch3){
		    if(!AcceptPair(&((fEvt+en2)->fTracks[j]), &((fEvt+en3)->fTracks[k]))) continue;
		  }

		  //
		  // The below call to SetFillBins3 will work for all 3-particle terms since all are for fully mixed events. part is set to 1, but only matters for terms 2-4.
		  Bool_t fill2=kFALSE, fill3=kFALSE, fill4=kFALSE;
		  SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 1, bin1, bin2, bin3, fill2, fill3, fill4);
		  
		  
		  for(Int_t jj=1; jj<=5; jj++){// term loop
		    
		    if(jj==2) {if(!fill2) continue;}//12
		    if(jj==3) {if(!fill3) continue;}//13
		    if(jj==4) {if(!fill4) continue;}//23
		    
		    Float_t WInput=1.0;
		    ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, 1, jj, firstQ, secondQ, thirdQ);
		    ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, 1, jj, firstQMC, secondQMC, thirdQMC);
		    
		    if(ch1==ch2 && ch1==ch3){// same charge
		      WInput = MCWeight3D(kTRUE, jj, 10, firstQMC, secondQMC, thirdQMC);
		      if(jj==1) {
			((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SC"))->Fill(q3MC, WInput);
			((TH1D*)fOutputList->FindObject("fMCWeight3DTerm1SCden"))->Fill(q3MC);
		      }else if(jj!=5){
			((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SC"))->Fill(q3MC, WInput);
			((TH1D*)fOutputList->FindObject("fMCWeight3DTerm2SCden"))->Fill(q3MC);
		      }
		    }else {// mixed charge
		      if(bin1==bin2) {
			WInput = MCWeight3D(kFALSE, jj, 10, firstQMC, secondQMC, thirdQMC);
		      }else {
			if(jj==1 || jj==5) WInput = MCWeight3D(kFALSE, jj, 10, thirdQMC, secondQMC, firstQMC);// thirdQMC is ss
			else WInput = MCWeight3D(kFALSE, 6-jj, 10, thirdQMC, secondQMC, firstQMC);
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
		    
		    
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fIdeal->Fill(q3MC, WInput);
		    Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fSmeared->Fill(q3, WInput);
		    
		    		    
		  }// jj
		}// MCarray check, 3rd particle
	      }// 3rd particle
	      
	    }// TabulatePairs Check
	    
	  }// MCarray check, 1st and 2nd particle
	  
	  // reset key's and fill bins (they were altered for 3 particle MRC calculation)
	  key1 = (fEvt)->fTracks[i].fKey;
	  key2 = (fEvt+en2)->fTracks[j].fKey;
	  SetFillBins2(fillIndex2, key1, key2, ch1, ch2, bin1, bin2);

	  
	  if(ch1==ch2 && fMbin==0 && qinv12<0.2){
	    //////////////////////////
	    // pad-row method testing
	    Float_t coeff = (5)*0.2*(0.18/1.2);// 5 to evaluate at 1.0m in TPC
	    Float_t phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	    if(phi1 > 2*PI) phi1 -= 2*PI;
	    if(phi1 < 0) phi1 += 2*PI;
	    Float_t phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	    if(phi2 > 2*PI) phi2 -= 2*PI;
	    if(phi2 < 0) phi2 += 2*PI;
	    Float_t deltaphi = phi1 - phi2;
	    if(deltaphi > PI) deltaphi -= PI;
	    if(deltaphi < -PI) deltaphi += PI;
	    
	    Int_t ncl1 = (fEvt)->fTracks[i].fClusterMap.GetNbits();
	    Int_t ncl2 = (fEvt+en2)->fTracks[j].fClusterMap.GetNbits();
	    Float_t sumCls = 0; Float_t sumSha = 0; Float_t sumQ = 0;
	    Double_t shfrac = 0; //Double_t qfactor = 0;
	    for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
	      if ((fEvt)->fTracks[i].fClusterMap.TestBitNumber(imap) && (fEvt+en2)->fTracks[j].fClusterMap.TestBitNumber(imap)) {// Both clusters
		if ((fEvt)->fTracks[i].fSharedMap.TestBitNumber(imap) && (fEvt+en2)->fTracks[j].fSharedMap.TestBitNumber(imap)) { // Shared
		  sumQ++;
		  sumCls+=2;
		  sumSha+=2;}
		else {sumQ--; sumCls+=2;}
	      }
	      else if ((fEvt)->fTracks[i].fClusterMap.TestBitNumber(imap) || (fEvt+en2)->fTracks[j].fClusterMap.TestBitNumber(imap)) {// Non shared
		sumQ++;
		sumCls++;}
	    }
	    if (sumCls>0) {
	      //qfactor = sumQ*1.0/sumCls;
	      shfrac = sumSha*1.0/sumCls;
	    }
	    if(fabs(deltaphi)<0.07 && fabs((fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta) < 0.03){
	      ((TH3D*)fOutputList->FindObject("fPairsPadRowDen"))->Fill(transK12, shfrac, qinv12);
	    }
	    
	    for(Int_t rstep=0; rstep<10; rstep++){
	      coeff = (rstep)*0.2*(0.18/1.2);
	      // propagate through B field to r=1.2m
	      phi1 = (fEvt)->fTracks[i].fPhi - asin((fEvt)->fTracks[i].fCharge*(0.1*fBfield)*coeff/(fEvt)->fTracks[i].fPt);
	      if(phi1 > 2*PI) phi1 -= 2*PI;
	      if(phi1 < 0) phi1 += 2*PI;
	      phi2 = (fEvt+en2)->fTracks[j].fPhi - asin((fEvt+en2)->fTracks[j].fCharge*(0.1*fBfield)*coeff/(fEvt+en2)->fTracks[j].fPt);
	      if(phi2 > 2*PI) phi2 -= 2*PI;
	      if(phi2 < 0) phi2 += 2*PI;
	      deltaphi = phi1 - phi2;
	      if(deltaphi > PI) deltaphi -= PI;
	      if(deltaphi < -PI) deltaphi += PI;
	      
	      if(fabs((fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta) < 0.03){
		((TH3F*)fOutputList->FindObject("fPairsShareFracDPhiDen"))->Fill(rstep, shfrac, deltaphi);
	      }
	      //if(shfrac < 0.05){
	      ((TH3F*)fOutputList->FindObject("fPairsDetaDPhiDen"))->Fill(rstep, (fEvt)->fTracks[i].fEta-(fEvt+en2)->fTracks[j].fEta, deltaphi);
	      //}
	    }
	    
	   
	    
	    
	  }// desired pair selection
	  
	
  
	}// fMCcase
	
	

	if(qinv12 < fQLowerCut) continue;// remove unwanted low-q pairs (also a type of track splitting cut)
	if(ch1 == ch2 && !fGeneratorOnly){
	  if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) {
	    fPairSplitCut[1][i]->AddAt('1',j);
	    continue;
	  }
	}
	
	//////////////////////////////////////////
	// 2-particle term
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2->Fill(transK12, qinv12);
	Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fExplicit2QW->Fill(transK12, qinv12, qinv12);
	if(qinv12<fQupperBound) Charge1[bin1].Charge2[bin2].SC[fillIndex2].MB[fMbin].EDB[fEDbin].TwoPT[en2].fMeanKt->Fill(transK12);
	
	//////////////////////////////////////////

	
	
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
	
	if(ch1 == ch2 && !fGeneratorOnly){
	  if(!AcceptPair(&((fEvt)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) {
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
	
	if(ch1 == ch2 && !fGeneratorOnly){
	  if(!AcceptPair(&((fEvt+en1)->fTracks[i]), &((fEvt+en2)->fTracks[j]))) {
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
	  qinv12 = (fEvt)->fPairsSE[p1].fQinv;
	  //
	  pVect1MC[1] = (fEvt)->fPairsSE[p1].fP1MC[0]; pVect2MC[1] = (fEvt)->fPairsSE[p1].fP2MC[0];
          pVect1MC[2] = (fEvt)->fPairsSE[p1].fP1MC[1]; pVect2MC[2] = (fEvt)->fPairsSE[p1].fP2MC[1];
          pVect1MC[3] = (fEvt)->fPairsSE[p1].fP1MC[2]; pVect2MC[3] = (fEvt)->fPairsSE[p1].fP2MC[2];
          pVect1MC[0] = sqrt(pow(pVect1MC[1],2)+pow(pVect1MC[2],2)+pow(pVect1MC[3],2)+pow(fTrueMassPi,2));
          pVect2MC[0] = sqrt(pow(pVect2MC[1],2)+pow(pVect2MC[2],2)+pow(pVect2MC[3],2)+pow(fTrueMassPi,2));
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
	  qinv12 = (fEvt)->fPairsME[p1].fQinv;
	  //
	  pVect1MC[1] = (fEvt)->fPairsME[p1].fP1MC[0]; pVect2MC[1] = (fEvt)->fPairsME[p1].fP2MC[0];
          pVect1MC[2] = (fEvt)->fPairsME[p1].fP1MC[1]; pVect2MC[2] = (fEvt)->fPairsME[p1].fP2MC[1];
          pVect1MC[3] = (fEvt)->fPairsME[p1].fP1MC[2]; pVect2MC[3] = (fEvt)->fPairsME[p1].fP2MC[2];
          pVect1MC[0] = sqrt(pow(pVect1MC[1],2)+pow(pVect1MC[2],2)+pow(pVect1MC[3],2)+pow(fTrueMassPi,2));
          pVect2MC[0] = sqrt(pow(pVect2MC[1],2)+pow(pVect2MC[2],2)+pow(pVect2MC[3],2)+pow(fTrueMassPi,2));
	}
	
	
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

	   
	    
	    if(fMCcase){
              pVect3MC[1] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPx;
              pVect3MC[2] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPy;
              pVect3MC[3] = (fEvt+en2)->fMCtracks[abs((fEvt+en2)->fTracks[k].fLabel)].fPz;
              pVect3MC[0] = sqrt(pow(pVect3MC[1],2)+pow(pVect3MC[2],2)+pow(pVect3MC[3],2)+pow(fTrueMassPi,2));
              qinv12MC = GetQinv(0, pVect1MC, pVect2MC);
              qinv13MC = GetQinv(0, pVect1MC, pVect3MC);
              qinv23MC = GetQinv(0, pVect2MC, pVect3MC);
            }

	    
	    
	    // if all three pair cuts are the same then the case (config=2 && term=2) never reaches here.
	    Float_t Kt12 = sqrt( pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.;
	    Float_t Kt13 = sqrt( pow(pVect1[1]+pVect3[1],2) + pow(pVect1[2]+pVect3[2],2))/2.;
	    Float_t Kt23 = sqrt( pow(pVect2[1]+pVect3[1],2) + pow(pVect2[2]+pVect3[2],2))/2.;
	    q3 = sqrt(pow(qinv12,2) + pow(qinv13,2) + pow(qinv23,2));
	    transK3 = sqrt( pow(pVect1[1]+pVect2[1]+pVect3[1],2) + pow(pVect1[2]+pVect2[2]+pVect3[2],2))/3.;
	    if(fKt3bins==1) fEDbin=0;
	    else if(fKt3bins==2){
	      if(transK3<0.3) fEDbin=0;
	      else fEDbin=1;
	    }else{// fKt3bins==3 is the limit set by fEDbins
	      if(transK3<0.25) fEDbin=0;
	      else if(transK3<0.35) fEDbin=1;
	      else fEDbin=2;
	    }
	    
	    firstQ=0; secondQ=0; thirdQ=0;
	    
	    
	    //
	    
	    if(config==1) {// 123
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 0, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	      
	      if(fillIndex3 <= 2){
		ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, 0, 1, firstQ, secondQ, thirdQ);
		if(fillIndex3==0 && fMCcase) ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, 0, 1, firstQMC, secondQMC, thirdQMC);
		Float_t WInput = 1.0;
		if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, 1, 10, firstQ, secondQ, thirdQ);
		//if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, 1, fFixedLambdaBinr3, firstQMC, secondQMC, thirdQMC);
		////
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fTermsQ3->Fill(q3, WInput);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fTerms3->Fill(firstQ, secondQ, thirdQ, WInput);
		//
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fMeanKt->Fill(Kt12);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fMeanKt->Fill(Kt13);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[0].fMeanKt->Fill(Kt23);
		////
		//
		if(fillIndex3==0 && ch1==ch2 && ch1==ch3 && fMCcase==kFALSE){
		  ((TH3D*)fOutputList->FindObject("fKt3DistTerm1"))->Fill(fMbin+1, transK3, q3);
		}		
		
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
		  if(fillIndex3==0 && fMCcase) ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12MC, qinv13MC, qinv23MC, part, jj, firstQMC, secondQMC, thirdQMC);
		  Float_t WInput = 1.0;
		  if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, jj, 10, firstQ, secondQ, thirdQ);
		  //if(fGenerateSignal && ch1==ch2 && ch1==ch3) WInput = MCWeight3D(kTRUE, jj, fFixedLambdaBinr3, firstQMC, secondQMC, thirdQMC);
		  ////
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fTermsQ3->Fill(q3, WInput);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fTerms3->Fill(firstQ, secondQ, thirdQ, WInput);
		  //
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fMeanKt->Fill(Kt12);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fMeanKt->Fill(Kt13);
		  Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[jj-1].fMeanKt->Fill(Kt23);
		}
	      }
	      
	    }else {// config 3: All particles from different events
	      
	      SetFillBins3(fillIndex3, key1, key2, key3, ch1, ch2, ch3, 3, bin1, bin2, bin3, fDummyB, fDummyB, fDummyB);
	      
	      if(ch1==ch2 && ch1==ch3 && fillIndex3==0) {
		if(!fMCcase) ((TH3D*)fOutputList->FindObject("fKt3DistTerm5"))->Fill(fMbin+1, transK3, q3);
	      }	      
	      
	      if(fillIndex3 <= 2){
		ArrangeQs(fillIndex3, key1, key2, key3, ch1, ch2, ch3, qinv12, qinv13, qinv23, part, 5, firstQ, secondQ, thirdQ);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fTermsQ3->Fill(q3);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fTerms3->Fill(firstQ, secondQ, thirdQ);
		//
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fMeanKt->Fill(Kt12);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fMeanKt->Fill(Kt13);
		Charge1[bin1].Charge2[bin2].Charge3[bin3].SC[fillIndex3].MB[fMbin].EDB[fEDbin].ThreePT[4].fMeanKt->Fill(Kt23);
	      }
	      
	      
	    }// config 3
	  }// end 3rd particle
	}// en2
	
	
      }// p1
    }//en1
    
    ///////////////////
  }// end of PdensityPairs

 
  
  // Post output data.
  PostData(1, fOutputList);
  
}
//________________________________________________________________________
void AliThreePionRadii::Terminate(Option_t *) 
{
  // Called once at the end of the query
 
  cout<<"Done"<<endl;

}
//________________________________________________________________________
Bool_t AliThreePionRadii::AcceptPair(AliChaoticityTrackStruct *first, AliChaoticityTrackStruct *second)
{
  
  if(fabs(first->fEta-second->fEta) > fMinSepPairEta) return kTRUE;
  
  // propagate through B field to r=1m
  Float_t phi1 = first->fPhi - asin(first->fCharge*(0.1*fBfield)*0.15/first->fPt);// 0.15 for D=1m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  Float_t phi2 = second->fPhi - asin(second->fCharge*(0.1*fBfield)*0.15/second->fPt);// 0.15 for D=1m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  Float_t deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
    
  
  // propagate through B field to r=1.6m
  phi1 = first->fPhi - asin(first->fCharge*(0.1*fBfield)*0.24/first->fPt);// mine. 0.24 for D=1.6m
  if(phi1 > 2*PI) phi1 -= 2*PI;
  if(phi1 < 0) phi1 += 2*PI;
  phi2 = second->fPhi - asin(second->fCharge*(0.1*fBfield)*0.24/second->fPt);// mine. 0.24 for D=1.6m 
  if(phi2 > 2*PI) phi2 -= 2*PI;
  if(phi2 < 0) phi2 += 2*PI;
  
  deltaphi = phi1 - phi2;
  if(deltaphi > PI) deltaphi -= 2*PI;
  if(deltaphi < -PI) deltaphi += 2*PI;
  deltaphi = fabs(deltaphi);

  if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
  
  
   
  //
  
  Int_t ncl1 = first->fClusterMap.GetNbits();
  Int_t ncl2 = second->fClusterMap.GetNbits();
  Int_t sumCls = 0; Int_t sumSha = 0; Int_t sumQ = 0;
  Double_t shfrac = 0; Double_t qfactor = 0;
  for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
    if (first->fClusterMap.TestBitNumber(imap) && second->fClusterMap.TestBitNumber(imap)) {// Both clusters
      if (first->fSharedMap.TestBitNumber(imap) && second->fSharedMap.TestBitNumber(imap)) { // Shared
	sumQ++;
	sumCls+=2;
	sumSha+=2;}
      else {sumQ--; sumCls+=2;}
    }
    else if (first->fClusterMap.TestBitNumber(imap) || second->fClusterMap.TestBitNumber(imap)) {// Non shared
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
Float_t AliThreePionRadii::GamovFactor(Int_t chargeBin1, Int_t chargeBin2, Float_t qinv)
{
  Float_t arg = G_Coeff/qinv;
  
  if(chargeBin1==chargeBin2) return (exp(arg)-1)/(arg);
  else {return (exp(-arg)-1)/(-arg);}
  
}
//________________________________________________________________________
void AliThreePionRadii::Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
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
Short_t AliThreePionRadii::FillIndex2part(Short_t key){

  if(key==2) return 0;// pi-pi
  else if(key==11) return 1;// pi-k
  else if(key==101) return 2;// pi-p
  else if(key==20) return 3;// k-k
  else if(key==110) return 4;// k-p
  else return 5;// p-p
}
//________________________________________________________________________
Short_t AliThreePionRadii::FillIndex3part(Short_t key){
  
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
Short_t AliThreePionRadii::SetQcutBin(Short_t fi){// fi=FillIndex
  if(fi <= 2) return 0;
  else if(fi==3) return 1;
  else return 2;
}
//________________________________________________________________________
Short_t AliThreePionRadii::SetNormBin(Short_t fi){// fi=FillIndex
  if(fi==0) return 0;
  else if(fi <= 2) return 1;
  else return 2;
}
//________________________________________________________________________
void AliThreePionRadii::SetFillBins2(Short_t fi, Short_t key1, Short_t key2, Int_t c1, Int_t c2, Int_t &b1, Int_t &b2){
  
  if(fi==0 || fi==3 || fi==5){// Identical species
    if((c1+c2)==1) {b1=0; b2=1;}// Re-assign to merge degenerate histos
    else {b1=c1; b2=c2;}
  }else {// Mixed species
    if(key1 < key2) { b1=c1; b2=c2;}
    else {b1=c2; b2=c1;}
  }
  
}
//________________________________________________________________________
void AliThreePionRadii::SetFillBins3(Short_t fi, Short_t key1, Short_t key2, Short_t key3, Int_t c1, Int_t c2, Int_t c3, Short_t part, Int_t &b1, Int_t &b2, Int_t &b3, Bool_t &fill2, Bool_t &fill3, Bool_t &fill4){
  
  
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
void AliThreePionRadii::ArrangeQs(Short_t fi, Short_t key1, Short_t key2, Short_t key3, Int_t c1, Int_t c2, Int_t c3, Float_t q12, Float_t q13, Float_t q23, Short_t part, Short_t term, Float_t &fQ, Float_t &sQ, Float_t &tQ){
 
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
Float_t AliThreePionRadii::GetQinv(Short_t fi, Float_t track1[], Float_t track2[]){
  
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
void AliThreePionRadii::GetQosl(Float_t track1[], Float_t track2[], Float_t& qout, Float_t& qside, Float_t& qlong){
 
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
Float_t AliThreePionRadii::MCWeight(Int_t charge1, Int_t charge2, Float_t r, Int_t dampIndex, Float_t qinv){
  
  Float_t radius = r/0.19733;// convert to GeV (starts at 5 fm, was 3 fm)
  
  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t coulCorr12 = FSICorrelation2(charge1, charge2, qinv);
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
Float_t AliThreePionRadii::MCWeight3D(Bool_t SameCharge, Int_t term, Int_t dampIndex, Float_t q12, Float_t q13, Float_t q23){
  if(term==5) return 1.0;
  
  Float_t radius=fRMax;
  radius /= 0.19733;

  Float_t myDamp = fDampStart + (fDampStep)*dampIndex;
  Float_t fc = sqrt(myDamp);
  
  if(SameCharge){// all three of the same charge
    Float_t coulCorr12 = FSICorrelation2(+1,+1, q12);// K2
    Float_t coulCorr13 = FSICorrelation2(+1,+1, q13);// K2
    Float_t coulCorr23 = FSICorrelation2(+1,+1, q23);// K2
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
      w123 += pow(fc,3)*c3QS*coulCorr12*coulCorr13*coulCorr23;
      return w123;
    }else if(term==2){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12);
    }else if(term==3){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q13*radius,2))*pow(EW13,2))*coulCorr13);
    }else if(term==4){
      return ((1-myDamp) + myDamp*(1 + exp(-pow(q23*radius,2))*pow(EW23,2))*coulCorr23);
    }else return 1.0;
  
  }else{// mixed charge case pair 12 always treated as ss
    Float_t coulCorr12 = FSICorrelation2(+1,+1, q12);// K2 ss
    Float_t coulCorr13 = FSICorrelation2(+1,-1, q13);// K2 os
    Float_t coulCorr23 = FSICorrelation2(+1,-1, q23);// K2 os
    Float_t arg12=q12*radius;
    Float_t EW12 = 1 + kappa3/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
    EW12 += kappa4/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
    if(term==1){
      Float_t c3QS = 1 + exp(-pow(q12*radius,2))*pow(EW12,2);
      Float_t w123 = pow(1-fc,3) + 3*fc*pow(1-fc,2);
      w123 += pow(fc,2)*(1-fc)*(1+exp(-pow(q12*radius,2))*pow(EW12,2))*coulCorr12;
      w123 += pow(fc,2)*(1-fc)*coulCorr13;
      w123 += pow(fc,2)*(1-fc)*coulCorr23;
      w123 += pow(fc,3)*c3QS*coulCorr12*coulCorr13*coulCorr23;
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
void AliThreePionRadii::SetFSICorrelations(Bool_t legoCase, TH1D *temp1DSS[10], TH1D *temp1DOS[10]){
  // read in 2-particle FSI correlations = K2 
  // 2-particle input histo from file is binned in qinv.
  if(legoCase){
    cout<<"LEGO call to SetFSICorrelations"<<endl;
    for(int mb=0; mb<10; mb++){
      fFSI2SS[mb] = (TH1D*)temp1DSS[mb]->Clone();
      fFSI2OS[mb] = (TH1D*)temp1DOS[mb]->Clone();
      //
      fFSI2SS[mb]->SetDirectory(0);
      fFSI2OS[mb]->SetDirectory(0);
    }
  }else {
    cout<<"non LEGO call to SetFSICorrelations"<<endl;
    TFile *fsifile = new TFile("KFile.root","READ");
    if(!fsifile->IsOpen()) {
      cout<<"No FSI file found"<<endl;
      AliFatal("No FSI file found.  Kill process.");
    }else {cout<<"Good FSI File Found!"<<endl;}
    for(int mb=0; mb<10; mb++){
      TString *nameK2ss = new TString("K2ss_");
      TString *nameK2os = new TString("K2os_");
      *nameK2ss += mb;
      *nameK2os += mb;
      TH1D *temphistoSS = (TH1D*)fsifile->Get(nameK2ss->Data());
      TH1D *temphistoOS = (TH1D*)fsifile->Get(nameK2os->Data());
      //
      fFSI2SS[mb] = (TH1D*)temphistoSS->Clone();
      fFSI2OS[mb] = (TH1D*)temphistoOS->Clone();
      fFSI2SS[mb]->SetDirectory(0);
      fFSI2OS[mb]->SetDirectory(0);
    }
    
    fsifile->Close();
  }

  for(int mb=0; mb<10; mb++){
    for(Int_t ii=1; ii<=fFSI2SS[mb]->GetNbinsX(); ii++){
      if(fFSI2SS[mb]->GetBinContent(ii) > 1.0) fFSI2SS[mb]->SetBinContent(ii, 1.0);
      if(fFSI2OS[mb]->GetBinContent(ii) > 10.0) fFSI2OS[mb]->SetBinContent(ii, 10.0);
      //
      if(fFSI2SS[mb]->GetBinContent(ii) < 0.05) fFSI2SS[mb]->SetBinContent(ii, 0.05);
      if(fFSI2OS[mb]->GetBinContent(ii) < 0.9) fFSI2OS[mb]->SetBinContent(ii, 0.9);
    }
  }
  
  cout<<"Done reading FSI file"<<endl;
}
//________________________________________________________________________
Float_t AliThreePionRadii::FSICorrelation2(Int_t charge1, Int_t charge2, Float_t qinv){
  // returns 2-particle Coulomb correlations = K2
  Int_t qbinL = fFSI2SS[fFSIindex]->GetXaxis()->FindBin(qinv-fFSI2SS[fFSIindex]->GetXaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSI2SS[fFSIindex]->GetNbinsX()) return 1.0;
  
  Float_t slope=0;
  if(charge1==charge2){
    slope = fFSI2SS[fFSIindex]->GetBinContent(qbinL) - fFSI2SS[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSI2SS[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSI2SS[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2SS[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSI2SS[fFSIindex]->GetBinContent(qbinL));
  }else {
    slope = fFSI2OS[fFSIindex]->GetBinContent(qbinL) - fFSI2OS[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSI2OS[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSI2OS[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSI2OS[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSI2OS[fFSIindex]->GetBinContent(qbinL));
  }
}
//________________________________________________________________________
