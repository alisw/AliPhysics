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
#include "AliHFMassFitter.h"
#include "AliHFMassFitterVAR.h"
#include "AliHFMultiTrials.h"

/// \cond CLASSIMP
ClassImp(AliHFMultiTrials);
/// \endcond


//_________________________________________________________________________
AliHFMultiTrials::AliHFMultiTrials() : 
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
  fSigmaMCVariation(0.15),
  fMassD(1.86484),
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
  fUseFixSigUpFreeMean(kTRUE),
  fUseFixSigDownFreeMean(kTRUE),
  fUseFreeS(kTRUE),
  fUseFixedMeanFreeS(kTRUE),
  fUseFixSigFreeMean(kTRUE),
  fUseFixSigFixMean(kTRUE),
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
  fHistoRawYieldDistBinCAll(0x0),
  fHistoRawYieldTrialBinCAll(0x0),
  fHistoRawYieldDist(0x0),
  fHistoRawYieldTrial(0x0),
  fHistoSigmaTrial(0x0),
  fHistoMeanTrial(0x0),
  fHistoChi2Trial(0x0),
  fHistoSignifTrial(0x0),
  fHistoBkgTrial(0x0),
  fHistoBkgInBinEdgesTrial(0x0),
  fHistoRawYieldDistBinC(0x0),
  fHistoRawYieldTrialBinC(0x0),
  fhTemplRefl(0x0),
  fFixRefloS(0),
  fNtupleMultiTrials(0x0),
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
AliHFMultiTrials::~AliHFMultiTrials(){
  // destructor
  delete [] fRebinSteps;
  delete [] fLowLimFitSteps;
  delete [] fUpLimFitSteps;
  if(fhTemplRefl) delete fhTemplRefl;
  for (auto fitter : fMassFitters) delete fitter;
}

//________________________________________________________________________
Bool_t AliHFMultiTrials::CreateHistos(){
  // creates output histograms

  const Int_t nCases=kNBkgFuncCases*kNFitConfCases;

  TString funcBkg[kNBkgFuncCases]={"Expo","Lin","Pol2","Pol3","Pol4","Pol5","PowLaw","PowLawExpo"};
  TString gausSig[kNFitConfCases]={"FixedS","FixedSp20","FixedSm20","FreeS","FixedMeanFixedS","FixedMeanFreeS"};

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


  fHistoRawYieldDistBinCAll = new TH1F(Form("hRawYieldDistBinCAll%s",fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
  fHistoRawYieldTrialBinCAll = new TH2F(Form("hRawYieldTrialBinCAll%s",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);

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

  fHistoRawYieldDistBinC = new TH1F*[nCases];
  fHistoRawYieldTrialBinC = new TH2F*[nCases];

  for(Int_t ib=0; ib<kNBkgFuncCases; ib++){
    for(Int_t igs=0; igs<kNFitConfCases; igs++){
      Int_t theCase=igs*kNBkgFuncCases+ib;
      fHistoRawYieldDist[theCase]=new TH1F(Form("hRawYieldDist%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield",5000,0.,50000.);
      fHistoRawYieldDistBinC[theCase]=new TH1F(Form("hRawYieldDistBinC%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",5000,0.,50000.);
      fHistoRawYieldTrial[theCase]=new TH1F(Form("hRawYieldTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Raw Yield",totTrials,-0.5,totTrials-0.5);
      fHistoRawYieldTrialBinC[theCase]=new TH2F(Form("hRawYieldTrialBinC%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
      fHistoSigmaTrial[theCase]=new TH1F(Form("hSigmaTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Sigma (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
      fHistoMeanTrial[theCase]=new TH1F(Form("hMeanTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Mean (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
      fHistoChi2Trial[theCase]=new TH1F(Form("hChi2Trial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; #chi^{2}",totTrials,-0.5,totTrials-0.5);
      fHistoSignifTrial[theCase]=new TH1F(Form("hSignifTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Significance",totTrials,-0.5,totTrials-0.5);
      if(fSaveBkgVal) {
        fHistoBkgTrial[theCase] = new TH1F(Form("hBkgTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Background",totTrials,-0.5,totTrials-0.5);
        fHistoBkgInBinEdgesTrial[theCase] = new TH1F(Form("hBkgInBinEdgesTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Background in bin edges",totTrials,-0.5,totTrials-0.5);
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
  fNtupleMultiTrials = new TNtuple(Form("ntuMultiTrial%s",fSuffix.Data()),Form("ntuMultiTrial%s",fSuffix.Data()),"rebin:firstb:minfit:maxfit:bkgfunc:confsig:confmean:chi2:signif:mean:emean:sigma:esigma:rawy:erawy",128000);
  fNtupleMultiTrials->SetDirectory(nullptr);
  return kTRUE;

}

//________________________________________________________________________
Bool_t AliHFMultiTrials::DoMultiTrials(TH1D* hInvMassHisto, TPad* thePad){
  // perform the multiple fits

  Bool_t hOK=CreateHistos();
  if(!hOK) return kFALSE;

  Int_t itrial=0;
  Int_t types=0;
  Int_t itrialBC=0;
  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;

  fMinYieldGlob=999999.;
  fMaxYieldGlob=0.;
  Float_t xnt[15];

  for(Int_t ir=0; ir<fNumOfRebinSteps; ir++){
    Int_t rebin=fRebinSteps[ir];
    for(Int_t iFirstBin=1; iFirstBin<=fNumOfFirstBinSteps; iFirstBin++) {
      TH1F* hRebinned=0x0;
      if(fNumOfFirstBinSteps==1) hRebinned=RebinHisto(hInvMassHisto,rebin,-1);
      else hRebinned=RebinHisto(hInvMassHisto,rebin,iFirstBin);
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
            for(Int_t igs=0; igs<kNFitConfCases; igs++){
              if (igs==kFixSigUpFreeMean && !fUseFixSigUpFreeMean) continue;
              if (igs==kFixSigDownFreeMean && !fUseFixSigDownFreeMean) continue;
              if (igs==kFreeSigFixMean  && !fUseFixedMeanFreeS) continue;
              if (igs==kFreeSigFreeMean  && !fUseFreeS) continue;
              if (igs==kFixSigFreeMean  && !fUseFixSigFreeMean) continue;
              if (igs==kFixSigFixMean   && !fUseFixSigFixMean) continue;
              Int_t theCase=igs*kNBkgFuncCases+typeb;
              Int_t globBin=itrial+theCase*totTrials;
              for(Int_t j=0; j<15; j++) xnt[j]=0.;

              Bool_t mustDeleteFitter = kTRUE;
              AliHFMassFitterVAR*  fitter=0x0;
              //if D0 Reflection
              if(fhTemplRefl){
                fitter=new AliHFMassFitterVAR(hRebinned,hmin,hmax,1,typeb,2);
                fitter->SetTemplateReflections(fhTemplRefl);
                fitter->SetFixReflOverS(fFixRefloS,kTRUE);
              }
              else {
                if(typeb<=kPol2Bkg){
                  fitter=new AliHFMassFitterVAR(hRebinned,hmin, hmax,1,typeb,types);
                }else if(typeb==kPowBkg){
                  fitter=new AliHFMassFitterVAR(hRebinned,hmin, hmax,1,4,types);
                }else if(typeb==kPowTimesExpoBkg){
                  fitter=new AliHFMassFitterVAR(hRebinned,hmin, hmax,1,5,types);
                }else{
                  fitter=new AliHFMassFitterVAR(hRebinned,hmin, hmax,1,6,types);
                  if(typeb==kPol3Bkg) fitter->SetBackHighPolDegree(3);
                  if(typeb==kPol4Bkg) fitter->SetBackHighPolDegree(4);
                  if(typeb==kPol5Bkg) fitter->SetBackHighPolDegree(5);
                }
                fitter->SetReflectionSigmaFactor(0);
              }
              if(fFitOption==1) fitter->SetUseChi2Fit();
              fitter->SetInitialGaussianMean(fMassD);
              fitter->SetInitialGaussianSigma(fSigmaGausMC);
              xnt[0]=rebin;
              xnt[1]=iFirstBin;
              xnt[2]=minMassForFit;
              xnt[3]=maxMassForFit;
              xnt[4]=typeb;
              xnt[6]=0;
              if(igs==kFixSigFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC,kTRUE);
                xnt[5]=1;
              }else if(igs==kFixSigUpFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC*(1.+fSigmaMCVariation),kTRUE);
                xnt[5]=2;
              }else if(igs==kFixSigDownFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC*(1.-fSigmaMCVariation),kTRUE);
                xnt[5]=3;
              }else if(igs==kFreeSigFreeMean){
                xnt[5]=0;
              }else if(igs==kFixSigFixMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC,kTRUE);
                fitter->SetFixGaussianMean(fMassD,kTRUE);
                xnt[5]=1;
                xnt[6]=1;
              }else if(igs==kFreeSigFixMean){
                fitter->SetFixGaussianMean(fMassD,kTRUE);
                xnt[5]=0;
                xnt[6]=1;
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
              TF1* fB1=0x0;
              if(typeb<kNBkgFuncCases){
                printf("****** START FIT OF HISTO %s WITH REBIN %d FIRST BIN %d MASS RANGE %f-%f BACKGROUND FIT FUNCTION=%d CONFIG SIGMA/MEAN=%d\n",hInvMassHisto->GetName(),rebin,iFirstBin,minMassForFit,maxMassForFit,typeb,igs);
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
                fB1=fitter->GetBackgroundFullRangeFunc();
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
              // else{
              //   out=DoFitWithPol3Bkg(hRebinned,hmin,hmax,igs);
              //   if(out && thePad){
              // 	thePad->Clear();
              // 	hRebinned->Draw();
              // 	TF1* fSB=(TF1*)hRebinned->GetListOfFunctions()->FindObject("fSB");
              // 	fB1=new TF1("fB1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",hmin,hmax);
              // 	for(Int_t j=0; j<4; j++) fB1->SetParameter(j,fSB->GetParameter(3+j));
              // 	fB1->SetLineColor(2);
              // 	fB1->Draw("same");
              // 	fSB->SetLineColor(4);
              // 	fSB->Draw("same");
              // 	thePad->Update();
              // 	chisq=fSB->GetChisquare()/fSB->GetNDF();;
              // 	sigma=fSB->GetParameter(2);
              // 	esigma=fSB->GetParError(2);
              // 	if(esigma<0.00001) esigma=0.0001;
              // 	pos=fSB->GetParameter(1);
              // 	epos=fSB->GetParError(1);
              // 	if(epos<0.00001) epos=0.0001;
              // 	ry=fSB->GetParameter(0)/hRebinned->GetBinWidth(1);
              // 	ery=fSB->GetParError(0)/hRebinned->GetBinWidth(1);
              //   }
              // }
              xnt[7]=chisq;
              if(out && chisq>0. && sigma>0.5*fSigmaGausMC && sigma<2.0*fSigmaGausMC){
                xnt[8]=significance;
                xnt[9]=pos;
                xnt[10]=epos;
                xnt[11]=sigma;
                xnt[12]=esigma;
                xnt[13]=ry;
                xnt[14]=ery;
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

                for(Int_t iStepBC=0; iStepBC<fNumOfnSigmaBinCSteps; iStepBC++){
                  Double_t minMassBC=fMassD-fnSigmaBinCSteps[iStepBC]*sigma;
                  Double_t maxMassBC=fMassD+fnSigmaBinCSteps[iStepBC]*sigma;
                  if(minMassBC>minMassForFit &&
                      maxMassBC<maxMassForFit &&
                      minMassBC>(hRebinned->GetXaxis()->GetXmin()) &&
                      maxMassBC<(hRebinned->GetXaxis()->GetXmax())){
                    Double_t cnts,ecnts;
                    BinCount(hRebinned,fB1,1,minMassBC,maxMassBC,cnts,ecnts);
                    ++itrialBC;
                    fHistoRawYieldDistBinCAll->Fill(cnts);
                    fHistoRawYieldTrialBinCAll->SetBinContent(globBin,iStepBC+1,cnts);
                    fHistoRawYieldTrialBinCAll->SetBinError(globBin,iStepBC+1,ecnts);
                    fHistoRawYieldTrialBinC[theCase]->SetBinContent(itrial,iStepBC+1,cnts);
                    fHistoRawYieldTrialBinC[theCase]->SetBinError(itrial,iStepBC+1,ecnts);
                    fHistoRawYieldDistBinC[theCase]->Fill(cnts);
                  }
                }
              }
              if (mustDeleteFitter) delete fitter;
              fNtupleMultiTrials->Fill(xnt);
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
void AliHFMultiTrials::SaveToRoot(TString fileName, TString option) const{
  // save histos in a root file for further analysis
  const Int_t nCases=kNBkgFuncCases*kNFitConfCases;
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
  fHistoRawYieldDistBinCAll->Write(); 
  fHistoRawYieldTrialBinCAll->Write(); 
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
    fHistoRawYieldTrialBinC[ic]->Write();
    fHistoRawYieldDistBinC[ic]->Write();
  }
  fNtupleMultiTrials->SetDirectory(&outHistos);
  fNtupleMultiTrials->Write();
  outHistos.Close();
}

//________________________________________________________________________
void AliHFMultiTrials::DrawHistos(TCanvas* cry) const{
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
TH1F* AliHFMultiTrials::RebinHisto(TH1D* hOrig, Int_t reb, Int_t firstUse) const{
  // Rebin histogram, from bin firstUse to lastUse
  // Use all bins if firstUse=-1

  Int_t nBinOrig=hOrig->GetNbinsX();
  Int_t firstBinOrig=1;
  Int_t lastBinOrig=nBinOrig;
  Int_t nBinOrigUsed=nBinOrig;
  Int_t nBinFinal=nBinOrig/reb;
  if(firstUse>=1){ 
    firstBinOrig=firstUse;
    nBinFinal=(nBinOrig-firstUse+1)/reb;
    nBinOrigUsed=nBinFinal*reb;
    lastBinOrig=firstBinOrig+nBinOrigUsed-1;
  }else{
    Int_t exc=nBinOrigUsed%reb;
    if(exc!=0){
      nBinOrigUsed-=exc;
      firstBinOrig+=exc/2;
      lastBinOrig=firstBinOrig+nBinOrigUsed-1;
    }
  }

  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n",nBinOrig,nBinFinal,nBinOrigUsed,firstBinOrig,lastBinOrig);
  Float_t lowLim=hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim=hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1F* hRebin=new TH1F(Form("%s-rebin%d_%d",hOrig->GetName(),reb,firstUse),hOrig->GetTitle(),nBinFinal,lowLim,hiLim);
  Int_t lastSummed=firstBinOrig-1;
  for(Int_t iBin=1;iBin<=nBinFinal; iBin++){
    Float_t sum=0.;
    Float_t sum2=0.;
    for(Int_t iOrigBin=0;iOrigBin<reb;iOrigBin++){
      sum+=hOrig->GetBinContent(lastSummed+1);
      sum2+=hOrig->GetBinError(lastSummed+1)*hOrig->GetBinError(lastSummed+1);
      lastSummed++;
    }
    hRebin->SetBinContent(iBin,sum);
    hRebin->SetBinError(iBin,TMath::Sqrt(sum2));
  }
  return hRebin;
}

//________________________________________________________________________
void AliHFMultiTrials::BinCount(TH1F* h, TF1* fB, Int_t rebin, Double_t minMass, Double_t maxMass, Double_t& count, Double_t& ecount) const{
  // compute yield with bin couting
  Int_t minBinSum=h->FindBin(minMass);
  Int_t maxBinSum=h->FindBin(maxMass);
  Double_t cntSig=0.;
  Double_t cntErr=0.;
  for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
    Double_t bkg=fB ? fB->Eval(h->GetBinCenter(iMB))/(Double_t)rebin : 0;
    cntSig+=(h->GetBinContent(iMB)-bkg);
    cntErr+=(h->GetBinError(iMB)*h->GetBinError(iMB));
  }
  count=cntSig;
  ecount=TMath::Sqrt(cntErr);
}
//________________________________________________________________________
Bool_t AliHFMultiTrials::DoFitWithPol3Bkg(TH1F* histoToFit, Double_t  hmin, Double_t  hmax, 
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
  else if(iCase==1) fSB->FixParameter(2,fSigmaGausMC*(1.+fSigmaMCVariation));
  else if(iCase==2) fSB->FixParameter(2,fSigmaGausMC*(1.-fSigmaMCVariation));
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
TH1F* AliHFMultiTrials::SetTemplateRefl(const TH1F *h) {
  fhTemplRefl=(TH1F*)h->Clone("hTemplRefl");
  return fhTemplRefl;
}

