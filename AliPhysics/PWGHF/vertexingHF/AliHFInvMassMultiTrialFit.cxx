/**************************************************************************
 * Copyright(c) 2008-2019, ALICE Experiment at CERN, All rights reserved. *
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

#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFile.h>
#include "AliHFInvMassFitter.h"
#include "AliHFInvMassMultiTrialFit.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFInvMassMultiTrialFit);
/// \endcond


//_________________________________________________________________________
AliHFInvMassMultiTrialFit::AliHFInvMassMultiTrialFit() : 
  TNamed(),
  fNumOfRebinSteps(4),
  fRebinSteps(0x0),
  fNumOfFirstBinSteps(1),
  fNumOfLowLimFitSteps(6),
  fLowLimFitSteps(0x0),
  fNumOfUpLimFitSteps(6),
  fUpLimFitSteps(0x0),
  fNumOfnSigmaBinCSteps(11),
  fnSigmaBinCSteps(0x0),
  fnSigmaForBkgEval(3),
  fSigmaGausMC(0.010),
  fSigmaMCVariationUp(0.1),
  fSigmaMCVariationDw(0.1),
  fMassD(1.86484),
  fUpperMassToFix(1.868),
  fLowerMassToFix(1.862),
  fSuffix(""),
  fFitOption(0),
  fUseExpoBkg(kTRUE),
  fUseLinBkg(kTRUE),  
  fUsePol2Bkg(kTRUE),
  fUsePol3Bkg(kTRUE),
  fUsePol4Bkg(kTRUE),
  fUsePol5Bkg(kFALSE),
  fUsePowLawBkg(kFALSE),
  fUsePowLawTimesExpoBkg(kFALSE),
  fUse2GausSignal(kFALSE),
  fUse2GausSigmaRatioSignal(kFALSE),
  fFixSecondGausSig(-1.),
  fFixSecondGausFrac(-1.),
  fFixSecondGausSigRat(-1.),
  fUseFixSigUpFreeMean(kTRUE),
  fUseFixSigDownFreeMean(kTRUE),
  fUseFreeS(kTRUE),
  fUseFixedMeanFreeS(kTRUE),
  fUseFixSigFreeMean(kTRUE),
  fUseFixSigFixMean(kTRUE),
  fUseFixSigFixMeanUp(kTRUE),
  fUseFixSigFixMeanDown(kTRUE),
  fUseFreeSigFixMeanUp(kTRUE),
  fUseFreeSigFixMeanDown(kTRUE),
  fUseFixSigVarWithFixMean(kTRUE),
  fUseSecondPeak(kFALSE),
  fMassSecondPeak(1.86958),
  fSigmaSecondPeak(0.01),
  fFixMassSecondPeak(kFALSE),
  fFixSigmaSecondPeak(kFALSE),
  fSaveBkgVal(kFALSE),
  fDrawIndividualFits(kFALSE),
  fHistoRawYieldDistAll(0x0),
  fHistoRawYieldTrialAll(0x0),
  fHistoSigmaTrialAll(0x0),
  fHistoMeanTrialAll(0x0),
  fHistoChi2TrialAll(0x0),
  fHistoSignifTrialAll(0x0),
  fHistoBkgTrialAll(0x0),
  fHistoBkgInBinEdgesTrialAll(0x0),
  fHistoRawYieldDistBinC0All(0x0),
  fHistoRawYieldTrialBinC0All(0x0),
  fHistoRawYieldDistBinC1All(0x0),
  fHistoRawYieldTrialBinC1All(0x0),
  fHistoRawYieldDist(0x0),
  fHistoRawYieldTrial(0x0),
  fHistoSigmaTrial(0x0),
  fHistoMeanTrial(0x0),
  fHistoChi2Trial(0x0),
  fHistoSignifTrial(0x0),
  fHistoBkgTrial(0x0),
  fHistoBkgInBinEdgesTrial(0x0),
  fHistoRawYieldDistBinC0(0x0),
  fHistoRawYieldTrialBinC0(0x0),
  fHistoRawYieldDistBinC1(0x0),
  fHistoRawYieldTrialBinC1(0x0),
  fhTemplRefl(0x0),
  fhTemplSign(0x0),
  fFixRefloS(1.),
  fNtupleMultiTrials(0x0),
  fNtupleBinCount(0x0),
  fMinYieldGlob(0),
  fMaxYieldGlob(0),
  fMassFitters()
{
  // constructor
  Int_t rebinStep[4]={3,4,5,6};
  Double_t minMassStep[6]={1.68,1.70,1.72,1.74,1.76,1.78};
  Double_t maxMassStep[6]={2.06,2.04,2.02,2.00,1.98,1.96};
  Double_t nSigmasBC[11]={2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0};
  ConfigureRebinSteps(4,rebinStep);
  ConfigureLowLimFitSteps(6,minMassStep);
  ConfigureUpLimFitSteps(6,maxMassStep);
  ConfigurenSigmaBinCSteps(11,nSigmasBC);
}

//________________________________________________________________________
AliHFInvMassMultiTrialFit::~AliHFInvMassMultiTrialFit(){
  // destructor
  delete [] fRebinSteps;
  delete [] fLowLimFitSteps;
  delete [] fUpLimFitSteps;
  if(fhTemplRefl) delete fhTemplRefl;
  if(fhTemplSign) delete fhTemplSign;
  for (auto fitter : fMassFitters) delete fitter;
  delete fHistoRawYieldDistAll;
  delete fHistoRawYieldTrialAll;
  delete fHistoSigmaTrialAll;
  delete fHistoMeanTrialAll;
  delete fHistoChi2TrialAll;
  delete fHistoSignifTrialAll;
  delete fHistoBkgTrialAll;
  delete fHistoBkgInBinEdgesTrialAll;
  delete fHistoRawYieldDistBinC0All;
  delete fHistoRawYieldTrialBinC0All;
  delete fHistoRawYieldDistBinC1All;
  delete fHistoRawYieldTrialBinC1All;
  for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
    for(Int_t types=0; types<kNSigFuncCases; types++){
      for(Int_t igs=0; igs<kNGausSigCases; igs++){
	for(Int_t igm=0; igm<kNGausMeanCases; igm++){
	  Int_t theCase=igm*kNGausSigCases*kNBkgFuncCases*kNSigFuncCases+igs*kNBkgFuncCases*kNSigFuncCases+types*kNBkgFuncCases+typeb;
	  delete fHistoRawYieldDist[theCase];
	  delete fHistoRawYieldDistBinC0[theCase];
	  delete fHistoRawYieldDistBinC1[theCase];
	  delete fHistoRawYieldTrial[theCase];
	  delete fHistoRawYieldTrialBinC0[theCase];
	  delete fHistoRawYieldTrialBinC1[theCase];
	  delete fHistoSigmaTrial[theCase];
	  delete fHistoMeanTrial[theCase];
	  delete fHistoChi2Trial[theCase];
	  delete fHistoSignifTrial[theCase];
	  if(fHistoBkgTrial) delete fHistoBkgTrial[theCase];
	  if(fHistoBkgInBinEdgesTrial) delete fHistoBkgInBinEdgesTrial[theCase];
	}
      }
    }
  }
  delete [] fHistoRawYieldDist;
  delete [] fHistoRawYieldDistBinC0;
  delete [] fHistoRawYieldDistBinC1;
  delete [] fHistoRawYieldTrial;
  delete [] fHistoRawYieldTrialBinC0;
  delete [] fHistoRawYieldTrialBinC1;
  delete [] fHistoSigmaTrial;
  delete [] fHistoMeanTrial;
  delete [] fHistoChi2Trial;
  delete [] fHistoSignifTrial;
  delete [] fHistoBkgTrial;
  delete [] fHistoBkgInBinEdgesTrial;

  delete fNtupleMultiTrials;
  delete fNtupleBinCount;
}

//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit::CreateHistos(){
  // creates output histograms

  const Int_t nCases=kNBkgFuncCases*kNGausSigCases*kNGausMeanCases*kNSigFuncCases;

  TString funcBkg[kNBkgFuncCases]={"Expo","Lin","Pol2","Pol3","Pol4","Pol5","PowLaw","PowLawExpo"};
  TString funcSig[kNSigFuncCases]={"","2Gaus","2GausSigmaRatio"};
  TString gausSig[kNGausSigCases]={"FixedSig","FixedSigUp","FixedSigDw","FreeSig"};
  TString gausMean[kNGausMeanCases]={"FixedMean","FixedMeanUp","FixedMeanDw","FreeMean"};

  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;

  fHistoRawYieldDistAll = new TH1F(Form("hRawYieldDistAll%s",fSuffix.Data()),"  ; Raw Yield",5000,0.,50000.);
  fHistoRawYieldTrialAll = new TH1F(Form("hRawYieldTrialAll%s",fSuffix.Data())," ; Trial # ; Raw Yield",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  fHistoSigmaTrialAll = new TH1F(Form("hSigmaTrialAll%s",fSuffix.Data())," ; Trial # ; Sigma (GeV/c^{2})",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  fHistoMeanTrialAll = new TH1F(Form("hMeanTrialAll%s",fSuffix.Data())," ; Trial # ; Mean (GeV/c^{2})",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  fHistoChi2TrialAll = new TH1F(Form("hChi2TrialAll%s",fSuffix.Data()),"  ; Trial # ; #chi^{2}",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  fHistoSignifTrialAll = new TH1F(Form("hSignifTrialAll%s",fSuffix.Data()),"  ; Trial # ; Significance",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  if(fSaveBkgVal) {
    fHistoBkgTrialAll = new TH1F(Form("hBkgTrialAll%s",fSuffix.Data()),"  ; Background",nCases*totTrials,-0.5,nCases*totTrials-0.5);
    fHistoBkgInBinEdgesTrialAll = new TH1F(Form("hBkgInBinEdgesTrialAll%s",fSuffix.Data()),"  ; Background in bin edges",nCases*totTrials,-0.5,nCases*totTrials-0.5);
  }


  fHistoRawYieldDistBinC0All = new TH1F(Form("hRawYieldDistBinC0All%s",fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
  fHistoRawYieldTrialBinC0All = new TH2F(Form("hRawYieldTrialBinC0All%s",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
  fHistoRawYieldDistBinC1All = new TH1F(Form("hRawYieldDistBinC1All%s",fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
  fHistoRawYieldTrialBinC1All = new TH2F(Form("hRawYieldTrialBinC1All%s",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);

  fHistoRawYieldDist = new TH1F*[nCases];
  fHistoRawYieldTrial = new TH1F*[nCases];
  fHistoSigmaTrial = new TH1F*[nCases];
  fHistoMeanTrial = new TH1F*[nCases];
  fHistoChi2Trial = new TH1F*[nCases];
  fHistoSignifTrial = new TH1F*[nCases];
  if(fSaveBkgVal) {
    fHistoBkgTrial = new TH1F*[nCases];
    fHistoBkgInBinEdgesTrial = new TH1F*[nCases];
  }

  fHistoRawYieldDistBinC0 = new TH1F*[nCases];
  fHistoRawYieldTrialBinC0 = new TH2F*[nCases];
  fHistoRawYieldDistBinC1 = new TH1F*[nCases];
  fHistoRawYieldTrialBinC1 = new TH2F*[nCases];

  for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
    for(Int_t types=0; types<kNSigFuncCases; types++){
      for(Int_t igs=0; igs<kNGausSigCases; igs++){
	for(Int_t igm=0; igm<kNGausMeanCases; igm++){
	  Int_t theCase=igm*kNGausSigCases*kNBkgFuncCases*kNSigFuncCases+igs*kNBkgFuncCases*kNSigFuncCases+types*kNBkgFuncCases+typeb;
	  printf("typeb=  %d   types = %d   igs = %d igm=%d  theCase=%d out of %d\n",typeb,types,igs,igm,theCase,nCases);
	  fHistoRawYieldDist[theCase]=new TH1F(Form("hRawYieldDist%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data()),"  ; Raw Yield",5000,0.,50000.);
	  fHistoRawYieldDistBinC0[theCase]=new TH1F(Form("hRawYieldDistBinC0%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
	  fHistoRawYieldDistBinC1[theCase]=new TH1F(Form("hRawYieldDistBinC1%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
	  fHistoRawYieldTrial[theCase]=new TH1F(Form("hRawYieldTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Raw Yield",totTrials,-0.5,totTrials-0.5);
	  fHistoRawYieldTrialBinC0[theCase]=new TH2F(Form("hRawYieldTrialBinC0%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
	  fHistoRawYieldTrialBinC1[theCase]=new TH2F(Form("hRawYieldTrialBinC1%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
	  fHistoSigmaTrial[theCase]=new TH1F(Form("hSigmaTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Sigma (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
	  fHistoMeanTrial[theCase]=new TH1F(Form("hMeanTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Mean (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
	  fHistoChi2Trial[theCase]=new TH1F(Form("hChi2Trial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; #chi^{2}",totTrials,-0.5,totTrials-0.5);
	  fHistoSignifTrial[theCase]=new TH1F(Form("hSignifTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data())," ; Trial # ; Significance",totTrials,-0.5,totTrials-0.5);
	  if(fSaveBkgVal) {
	    fHistoBkgTrial[theCase] = new TH1F(Form("hBkgTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data()),"  ; Background",totTrials,-0.5,totTrials-0.5);
	    fHistoBkgInBinEdgesTrial[theCase] = new TH1F(Form("hBkgInBinEdgesTrial%s%s%s%s%s",funcBkg[typeb].Data(),funcSig[types].Data(),gausSig[igs].Data(),gausMean[igm].Data(),fSuffix.Data()),"  ; Background in bin edges",totTrials,-0.5,totTrials-0.5);
	  }      
	
	  fHistoChi2Trial[theCase]->SetMarkerStyle(7);
	  fHistoSignifTrial[theCase]->SetMarkerStyle(7);
	  fHistoSigmaTrial[theCase]->SetMarkerStyle(7);
	  fHistoMeanTrial[theCase]->SetMarkerStyle(7);
	  if(fSaveBkgVal) {
	    fHistoBkgTrial[theCase]->SetMarkerStyle(7);
	    fHistoBkgInBinEdgesTrial[theCase]->SetMarkerStyle(7);
	  }
	}
      }
    }
  }
  fNtupleMultiTrials = new TNtuple(Form("ntuMultiTrial%s",fSuffix.Data()),Form("ntuMultiTrial%s",fSuffix.Data()),"rebin:firstb:minfit:maxfit:bkgfunc:sigfunc:confsig:confmean:chi2:signif:mean:emean:sigma:esigma:rawy:erawy",128000);
  fNtupleMultiTrials->SetDirectory(nullptr);
  fNtupleBinCount = new TNtuple(Form("ntuBinCount%s",fSuffix.Data()),Form("ntuBinCount%s",fSuffix.Data()),"rebin:firstb:minfit:maxfit:bkgfunc:sigfunc:confsig:confmean:chi2:nSigmaBC:rawyBC0:erawyBC0:rawyBC1:erawyBC1",128000);
  fNtupleBinCount->SetDirectory(nullptr);
  return kTRUE;

}

//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit::DoMultiTrials(TH1D* hInvMassHisto, TPad* thePad){
  // perform the multiple fits

  Bool_t hOK=CreateHistos();
  if(!hOK) return kFALSE;

  Int_t itrial=0;
  Int_t itrialBC=0;
  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;

  fMinYieldGlob=999999.;
  fMaxYieldGlob=0.;
  Float_t xnt[16];
  Float_t xntBC[14];

  for(Int_t ir=0; ir<fNumOfRebinSteps; ir++){
    Int_t rebin=fRebinSteps[ir];
    for(Int_t iFirstBin=1; iFirstBin<=fNumOfFirstBinSteps; iFirstBin++) {
      TH1F* hRebinned=0x0;
      if(fNumOfFirstBinSteps==1) hRebinned=(TH1F*)AliVertexingHFUtils::RebinHisto(hInvMassHisto,rebin,-1);
      else hRebinned=(TH1F*)AliVertexingHFUtils::RebinHisto(hInvMassHisto,rebin,iFirstBin);
      for(Int_t iMinMass=0; iMinMass<fNumOfLowLimFitSteps; iMinMass++){
        Double_t minMassForFit=fLowLimFitSteps[iMinMass];
        Double_t hmin=TMath::Max(minMassForFit,hRebinned->GetBinLowEdge(2));
        for(Int_t iMaxMass=0; iMaxMass<fNumOfUpLimFitSteps; iMaxMass++){
          Double_t maxMassForFit=fUpLimFitSteps[iMaxMass];
          Double_t hmax=TMath::Min(maxMassForFit,hRebinned->GetBinLowEdge(hRebinned->GetNbinsX()));
          ++itrial;
          for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
            if(typeb==kExpoBkg && !fUseExpoBkg) continue;
            if(typeb==kLinBkg && !fUseLinBkg) continue;
            if(typeb==kPol2Bkg && !fUsePol2Bkg) continue;
            if(typeb==kPol3Bkg && !fUsePol3Bkg) continue;
            if(typeb==kPol4Bkg && !fUsePol4Bkg) continue;
            if(typeb==kPol5Bkg && !fUsePol5Bkg) continue;
            if(typeb==kPowBkg && !fUsePowLawBkg) continue;
            if(typeb==kPowTimesExpoBkg && !fUsePowLawTimesExpoBkg) continue;  
	    for(Int_t types=0; types<kNSigFuncCases; types++){
	      if(types==k2Gaus && !fUse2GausSignal) continue;
	      if(types==k2GausSigmaRatioPar && !fUse2GausSigmaRatioSignal) continue;
	      for(Int_t igs=0; igs<kNGausSigCases; igs++){
		if (igs==kFreeSig && !fUseFreeS) continue;
		for(Int_t igm=0; igm<kNGausMeanCases; igm++){
		  if (igs==kFixSig){
		    if (igm==kFreeMean && !fUseFixSigFreeMean) continue;
		    if (igm==kFixMean && !fUseFixSigFixMean) continue;
		    if (igm==kFixMeanUp && !fUseFixSigFixMeanUp) continue;
		    if (igm==kFixMeanDown && !fUseFixSigFixMeanDown) continue;
		  }
		  if (igs==kFreeSig){
		    if (igm==kFixMean  && !fUseFixedMeanFreeS) continue;
		    if (igm==kFixMeanUp && !fUseFreeSigFixMeanUp) continue;
		    if (igm==kFixMeanDown && !fUseFreeSigFixMeanDown) continue;
		  }
		  if (igs==kFixSigUp){
		    if (igm==kFreeMean && !fUseFixSigUpFreeMean) continue;
		    if (igm==kFixMean && !fUseFixSigVarWithFixMean) continue;
		    if (igm==kFixMeanUp && !fUseFixSigVarWithFixMean) continue;
		    if (igm==kFixMeanDown && !fUseFixSigVarWithFixMean) continue;
		  }
		  if (igs==kFixSigDown){
		    if (igm==kFreeMean && !fUseFixSigDownFreeMean) continue;
		    if (igm==kFixMean && !fUseFixSigVarWithFixMean) continue;
		    if (igm==kFixMeanUp && !fUseFixSigVarWithFixMean) continue;
		    if (igm==kFixMeanDown && !fUseFixSigVarWithFixMean) continue;
		  }
		  
		  Int_t theCase=igm*kNGausSigCases*kNBkgFuncCases*kNSigFuncCases+igs*kNBkgFuncCases*kNSigFuncCases+types*kNBkgFuncCases+typeb;
		  Int_t globBin=itrial+theCase*totTrials;
		  for(Int_t j=0; j<16; j++) xnt[j]=0.;
		
		  Bool_t mustDeleteFitter = kTRUE;
		  AliHFInvMassFitter*  fitter=0x0;
		  if(typeb==kExpoBkg){
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, AliHFInvMassFitter::kExpo, types);
		  }else if(typeb==kLinBkg){
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, AliHFInvMassFitter::kLin, types);
		  }else if(typeb==kPol2Bkg){
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, AliHFInvMassFitter::kPol2, types);
		  }else if(typeb==kPowBkg){
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, AliHFInvMassFitter::kPow, types);
		  }else if(typeb==kPowTimesExpoBkg){
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, AliHFInvMassFitter::kPowEx, types);
		  }else{
		    fitter=new AliHFInvMassFitter(hRebinned, hmin, hmax, 6, types);
		    if(typeb==kPol3Bkg) fitter->SetPolDegreeForBackgroundFit(3);
		    if(typeb==kPol4Bkg) fitter->SetPolDegreeForBackgroundFit(4);
		    if(typeb==kPol5Bkg) fitter->SetPolDegreeForBackgroundFit(5);
		  }
		  if(types==k2Gaus){
		    if(fFixSecondGausSig>=0.) fitter->SetFixSecondGaussianSigma(fFixSecondGausSig);
		    if(fFixSecondGausFrac>=0.) fitter->SetFixFrac2Gaus(fFixSecondGausFrac);
		  }else if(types==k2GausSigmaRatioPar){
		    if(fFixSecondGausSigRat>=0.) fitter->SetFixRatio2GausSigma(fFixSecondGausSigRat);
		    if(fFixSecondGausFrac>=0.) fitter->SetFixFrac2Gaus(fFixSecondGausFrac);
		  }
		  // D0 Reflection
		  if(fhTemplRefl && fhTemplSign){
		    TH1F *hReflModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(fhTemplRefl,hRebinned,minMassForFit,maxMassForFit);
		    TH1F *hSigModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(fhTemplSign,hRebinned,minMassForFit,maxMassForFit);
		    TH1F* hrfl=fitter->SetTemplateReflections(hReflModif,"2gaus",minMassForFit,maxMassForFit);
		    if(!hrfl) printf("ERROR in SetTemplateReflections\n");
		    if(fFixRefloS>0){
		      Double_t fixSoverRefAt=fFixRefloS*(hReflModif->Integral(hReflModif->FindBin(minMassForFit*1.0001),hReflModif->FindBin(maxMassForFit*0.999))/hSigModif->Integral(hSigModif->FindBin(minMassForFit*1.0001),hSigModif->FindBin(maxMassForFit*0.999)));
		      fitter->SetFixReflOverS(fixSoverRefAt);
		    }
		    delete hReflModif;
		    delete hSigModif;
		  }
		  if(fUseSecondPeak){
		    fitter->IncludeSecondGausPeak(fMassSecondPeak, fFixMassSecondPeak, fSigmaSecondPeak, fFixSigmaSecondPeak);
		  }
		  if(fFitOption==1) fitter->SetUseChi2Fit();
		  fitter->SetInitialGaussianMean(fMassD);
		  fitter->SetInitialGaussianSigma(fSigmaGausMC);
		  xnt[0]=rebin;
		  xnt[1]=iFirstBin;
		  xnt[2]=minMassForFit;
		  xnt[3]=maxMassForFit;
		  xnt[4]=typeb;
		  xnt[5]=types;
		  if(igs==kFixSig){
		    fitter->SetFixGaussianSigma(fSigmaGausMC);
		    xnt[6]=1;
		  }else if(igs==kFixSigUp){
		    fitter->SetFixGaussianSigma(fSigmaGausMC*(1.+fSigmaMCVariationUp));
		    xnt[6]=2;
		  }else if(igs==kFixSigDown){
		    fitter->SetFixGaussianSigma(fSigmaGausMC*(1.-fSigmaMCVariationDw));
		    xnt[6]=3;
		  }else{ // igs==kFreeSig
		    xnt[6]=0;
		  }
		  if(igm==kFixMean){
		    fitter->SetFixGaussianMean(fMassD);
		    xnt[7]=1;
		  }else if(igm==kFixMeanUp){
		    fitter->SetFixGaussianMean(fUpperMassToFix);
		    xnt[7]=2;
		  }else if(igm==kFixMeanDown){
		    fitter->SetFixGaussianMean(fLowerMassToFix);
		    xnt[7]=3;
		  }else{ // igm==kFreeMean){
		    xnt[7]=0;
		  }
		  
		  Bool_t out=kFALSE;
		  Double_t chisq=-1.;
		  Double_t sigma=0.;
		  Double_t esigma=0.;
		  Double_t pos=.0;
		  Double_t epos=.0;
		  Double_t ry=.0;
		  Double_t ery=.0;
		  Double_t significance=0.;
		  Double_t erSignif=0.;
		  Double_t bkg=0.;
		  Double_t erbkg=0.;
		  Double_t bkgBEdge=0;
		  Double_t erbkgBEdge=0;
		  if(typeb<kNBkgFuncCases){
		    printf("****** START FIT OF HISTO %s WITH REBIN %d FIRST BIN %d MASS RANGE %f-%f BACKGROUND FIT FUNCTION=%d CONFIG SIGMA/MEAN=%d %d\n",hInvMassHisto->GetName(),rebin,iFirstBin,minMassForFit,maxMassForFit,typeb,igs,igm);
		    out=fitter->MassFitter(0);
		    chisq=fitter->GetReducedChiSquare();
		    fitter->Significance(fnSigmaForBkgEval,significance,erSignif);
		    sigma=fitter->GetSigma();
		    pos=fitter->GetMean();
		    esigma=fitter->GetSigmaUncertainty();
		    if(esigma<0.00001) esigma=0.0001;
		    epos=fitter->GetMeanUncertainty();
		    if(epos<0.00001) epos=0.0001;
		    ry=fitter->GetRawYield();
		    ery=fitter->GetRawYieldError();
		    fitter->Background(fnSigmaForBkgEval,bkg,erbkg);
		    Double_t minval = hInvMassHisto->GetXaxis()->GetBinLowEdge(hInvMassHisto->FindBin(pos-fnSigmaForBkgEval*sigma));
		    Double_t maxval = hInvMassHisto->GetXaxis()->GetBinUpEdge(hInvMassHisto->FindBin(pos+fnSigmaForBkgEval*sigma));
		    fitter->Background(minval,maxval,bkgBEdge,erbkgBEdge);
		    if(out && fDrawIndividualFits && thePad){
		      thePad->Clear();
		      fitter->DrawHere(thePad, fnSigmaForBkgEval);
		      fMassFitters.push_back(fitter);
		      mustDeleteFitter = kFALSE;
		      for (auto format : fInvMassFitSaveAsFormats) {
			thePad->SaveAs(Form("FitOutput_%s_Trial%d.%s",hInvMassHisto->GetName(),globBin, format.c_str()));
		      }
		    }
		  }
		  xnt[8]=chisq;
		  if(out && chisq>0. && sigma>0.5*fSigmaGausMC && sigma<2.0*fSigmaGausMC){
		    xnt[9]=significance;
		    xnt[10]=pos;
		    xnt[11]=epos;
		    xnt[12]=sigma;
		    xnt[13]=esigma;
		    xnt[14]=ry;
		    xnt[15]=ery;
		    fHistoRawYieldDistAll->Fill(ry);
		    fHistoRawYieldTrialAll->SetBinContent(globBin,ry);
		    fHistoRawYieldTrialAll->SetBinError(globBin,ery);
		    fHistoSigmaTrialAll->SetBinContent(globBin,sigma);
		    fHistoSigmaTrialAll->SetBinError(globBin,esigma);
		    fHistoMeanTrialAll->SetBinContent(globBin,pos);
		    fHistoMeanTrialAll->SetBinError(globBin,epos);
		    fHistoChi2TrialAll->SetBinContent(globBin,chisq);
		    fHistoChi2TrialAll->SetBinError(globBin,0.00001);
		    fHistoSignifTrialAll->SetBinContent(globBin,significance);
		    fHistoSignifTrialAll->SetBinError(globBin,erSignif);
		    if(fSaveBkgVal) {
		      fHistoBkgTrialAll->SetBinContent(globBin,bkg);
		      fHistoBkgTrialAll->SetBinError(globBin,erbkg);
		      fHistoBkgInBinEdgesTrialAll->SetBinContent(globBin,bkgBEdge);
		      fHistoBkgInBinEdgesTrialAll->SetBinError(globBin,erbkgBEdge);
		    }
		    
		    if(ry<fMinYieldGlob) fMinYieldGlob=ry;
		    if(ry>fMaxYieldGlob) fMaxYieldGlob=ry;
		    fHistoRawYieldDist[theCase]->Fill(ry);
		    fHistoRawYieldTrial[theCase]->SetBinContent(itrial,ry);
		    fHistoRawYieldTrial[theCase]->SetBinError(itrial,ery);
		    fHistoSigmaTrial[theCase]->SetBinContent(itrial,sigma);
		    fHistoSigmaTrial[theCase]->SetBinError(itrial,esigma);
		    fHistoMeanTrial[theCase]->SetBinContent(itrial,pos);
		    fHistoMeanTrial[theCase]->SetBinError(itrial,epos);
		    fHistoChi2Trial[theCase]->SetBinContent(itrial,chisq);
		    fHistoChi2Trial[theCase]->SetBinError(itrial,0.00001);
		    fHistoSignifTrial[theCase]->SetBinContent(itrial,significance);
		    fHistoSignifTrial[theCase]->SetBinError(itrial,erSignif);
		    if(fSaveBkgVal) {
		      fHistoBkgTrial[theCase]->SetBinContent(itrial,bkg);
		      fHistoBkgTrial[theCase]->SetBinError(itrial,erbkg);
		      fHistoBkgInBinEdgesTrial[theCase]->SetBinContent(itrial,bkgBEdge);
		      fHistoBkgInBinEdgesTrial[theCase]->SetBinError(itrial,erbkgBEdge);
		    }
		    fNtupleMultiTrials->Fill(xnt);
		    if(types==0){
		      // bin counting done only for 1 case of signal line shape
		      for(Int_t j=0; j<9; j++) xntBC[j]=xnt[j];
		      
		      
		      for(Int_t iStepBC=0; iStepBC<fNumOfnSigmaBinCSteps; iStepBC++){
			Double_t minMassBC=fMassD-fnSigmaBinCSteps[iStepBC]*sigma;
			Double_t maxMassBC=fMassD+fnSigmaBinCSteps[iStepBC]*sigma;
			if(minMassBC>minMassForFit &&
			   maxMassBC<maxMassForFit &&
			   minMassBC>(hRebinned->GetXaxis()->GetXmin()) &&
			   maxMassBC<(hRebinned->GetXaxis()->GetXmax())){
			  Double_t cnts0,ecnts0;
			  Double_t cnts1,ecnts1;
			  cnts0=fitter->GetRawYieldBinCounting(ecnts0,fnSigmaBinCSteps[iStepBC],0,0);
			  cnts1=fitter->GetRawYieldBinCounting(ecnts1,fnSigmaBinCSteps[iStepBC],1,0);
			  xntBC[9]=fnSigmaBinCSteps[iStepBC];
			  xntBC[10]=cnts0;
			  xntBC[11]=ecnts0;
			  xntBC[12]=cnts1;
			  xntBC[13]=ecnts1;
			  ++itrialBC;
			  fHistoRawYieldDistBinC0All->Fill(cnts0);
			  fHistoRawYieldTrialBinC0All->SetBinContent(globBin,iStepBC+1,cnts0);
			  fHistoRawYieldTrialBinC0All->SetBinError(globBin,iStepBC+1,ecnts0);
			  fHistoRawYieldTrialBinC0[theCase]->SetBinContent(itrial,iStepBC+1,cnts0);
			  fHistoRawYieldTrialBinC0[theCase]->SetBinError(itrial,iStepBC+1,ecnts0);
			  fHistoRawYieldDistBinC0[theCase]->Fill(cnts0);
			  fHistoRawYieldDistBinC1All->Fill(cnts1);
			  fHistoRawYieldTrialBinC1All->SetBinContent(globBin,iStepBC+1,cnts1);
			  fHistoRawYieldTrialBinC1All->SetBinError(globBin,iStepBC+1,ecnts1);
			  fHistoRawYieldTrialBinC1[theCase]->SetBinContent(itrial,iStepBC+1,cnts1);
			  fHistoRawYieldTrialBinC1[theCase]->SetBinError(itrial,iStepBC+1,ecnts1);
			  fHistoRawYieldDistBinC1[theCase]->Fill(cnts1);
			  fNtupleBinCount->Fill(xntBC);
			}
		      }
		    }
		  }
		  if (mustDeleteFitter) delete fitter;
		}
	      }
	    }
	  }
	}
      }
      delete hRebinned;
    }
  }
  return kTRUE;
}

//________________________________________________________________________
void AliHFInvMassMultiTrialFit::SaveToRoot(TString fileName, TString option) const{
  // save histos in a root file for further analysis
  const Int_t nCases=kNBkgFuncCases*kNGausSigCases*kNGausMeanCases*kNSigFuncCases;
  TFile outHistos(fileName.Data(),option.Data());
  if (outHistos.IsZombie()) {
    Printf("Could not open file '%s'!", fileName.Data());
    return;
  }
  outHistos.cd();
  fHistoRawYieldTrialAll->Write();
  fHistoSigmaTrialAll->Write();
  fHistoMeanTrialAll->Write();
  fHistoChi2TrialAll->Write();
  fHistoSignifTrialAll->Write();
  if(fSaveBkgVal) {
    fHistoBkgTrialAll->Write();
    fHistoBkgInBinEdgesTrialAll->Write();
  }
  fHistoRawYieldDistBinC0All->Write();
  fHistoRawYieldTrialBinC0All->Write();
  fHistoRawYieldDistBinC1All->Write();
  fHistoRawYieldTrialBinC1All->Write();
  for(Int_t ic=0; ic<nCases; ic++){
    fHistoRawYieldTrial[ic]->Write();
    fHistoSigmaTrial[ic]->Write();
    fHistoMeanTrial[ic]->Write();
    fHistoChi2Trial[ic]->Write();
    fHistoSignifTrial[ic]->Write();
    if(fSaveBkgVal) {
      fHistoBkgTrial[ic]->Write();
      fHistoBkgInBinEdgesTrial[ic]->Write();
    }
    fHistoRawYieldTrialBinC0[ic]->Write();
    fHistoRawYieldDistBinC0[ic]->Write();
    fHistoRawYieldTrialBinC1[ic]->Write();
    fHistoRawYieldDistBinC1[ic]->Write();
  }
  fNtupleMultiTrials->SetDirectory(&outHistos);
  fNtupleMultiTrials->Write();
  fNtupleBinCount->SetDirectory(&outHistos);
  fNtupleBinCount->Write();
  outHistos.Close();
}

//________________________________________________________________________
void AliHFInvMassMultiTrialFit::DrawHistos(TCanvas* cry) const{
  // draw histos
  cry->Divide(2,2);
  cry->cd(1);
  fHistoSigmaTrialAll->Draw();
  cry->cd(3);
  fHistoChi2TrialAll->Draw();
  cry->cd(2);
  fHistoRawYieldTrialAll->Draw();
  cry->cd(4);
  fHistoRawYieldDistAll->Draw();
  TLatex* tmean=new TLatex(0.15,0.8,Form("mean=%.1f",fHistoRawYieldDistAll->GetMean()));
  tmean->SetNDC();
  tmean->Draw();
  TLatex* thrms=new TLatex(0.15,0.72,Form("rms=%.1f",fHistoRawYieldDistAll->GetRMS()));
  thrms->SetNDC();
  thrms->Draw();
  TLatex* tmax=new TLatex(0.6,0.8,Form("max=%.1f",fMaxYieldGlob));
  tmax->SetNDC();
  tmax->Draw();
  TLatex* tmin=new TLatex(0.6,0.72,Form("min=%.1f",fMinYieldGlob));
  tmin->SetNDC();
  tmin->Draw();
  TLatex* trms=new TLatex(0.6,0.64,Form("(max-min)/sqrt(12)=%.1f",(fMaxYieldGlob-fMinYieldGlob)/sqrt(12)));
  trms->SetNDC();
  trms->Draw();

}
//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit::DoFitWithPol3Bkg(TH1F* histoToFit, Double_t  hmin, Double_t  hmax, 
    Int_t iCase){
  //

  TH1F *hCutTmp=(TH1F*)histoToFit->Clone("hCutTmp");
  for(Int_t ib=1; ib<=hCutTmp->GetNbinsX(); ib++){
    Double_t xc=hCutTmp->GetBinCenter(ib);
    if(xc>(fMassD-5.*fSigmaGausMC) && xc<(fMassD+5.*fSigmaGausMC)){
      hCutTmp->SetBinContent(ib,0.);
      hCutTmp->SetBinError(ib,0.);
    }
  }

  hCutTmp->Fit("pol2","E0","",hmin,hmax);
  TF1* f2=(TF1*)hCutTmp->GetListOfFunctions()->FindObject("pol2");
  TF1* f3=new TF1("myPol3","pol3");
  for(Int_t i=0; i<3;i++) f3->SetParameter(i,f2->GetParameter(i));
  hCutTmp->Fit(f3,"E0","",hmin,hmax);
  Double_t quickCount=0.;
  for(Int_t ib=1; ib<=histoToFit->GetNbinsX(); ib++){
    Double_t xc=hCutTmp->GetBinCenter(ib);
    if(xc>(fMassD-3.*fSigmaGausMC) && xc<(fMassD+3.*fSigmaGausMC)){
      quickCount+=(histoToFit->GetBinContent(ib)-f3->Eval(xc));
    }
  }
  TF1* fSB=new TF1("fSB","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x");
  fSB->SetParameter(0,quickCount);
  fSB->SetParameter(1,fMassD);
  fSB->SetParameter(2,fSigmaGausMC);
  for(Int_t j=0; j<4; j++) fSB->SetParameter(j+3,f3->GetParameter(j));
  if(iCase==0) fSB->FixParameter(2,fSigmaGausMC);
  else if(iCase==1) fSB->FixParameter(2,fSigmaGausMC*(1.+fSigmaMCVariationUp));
  else if(iCase==2) fSB->FixParameter(2,fSigmaGausMC*(1.-fSigmaMCVariationDw));
  else if(iCase==4){ 
    fSB->FixParameter(1,fMassD);
    fSB->FixParameter(2,fSigmaGausMC);
  } else if(iCase==5){ 
    fSB->FixParameter(1,fMassD);
  }
  histoToFit->Fit(fSB,"ME0","",hmin,hmax);
  // quality cuts
  if(fSB->GetParError(0)<0.01*fSB->GetParameter(0)) return kFALSE;
  if(fSB->GetParError(0)>0.6*fSB->GetParameter(0)) return kFALSE;

  delete hCutTmp;
  return kTRUE;
}
//__________________________________________________________________________________
void AliHFInvMassMultiTrialFit::SetTemplatesForReflections(const TH1F *hr, const TH1F *hs) {
  /// signal and reflection templates
  if(fhTemplSign) delete fhTemplSign;
  if(fhTemplRefl) delete fhTemplRefl;
  fhTemplRefl=(TH1F*)hr->Clone("hTemplRefl");
  fhTemplSign=(TH1F*)hs->Clone("hTemplSign");
  return;
}

