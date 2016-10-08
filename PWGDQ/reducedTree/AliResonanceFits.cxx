/*
***********************************************************
  Implementation of the AliResonanceFits
  Contact: i.c.arsene@cern.ch
  2014/10/09
  *********************************************************
*/

#ifndef ALIRESONANCEFITS_H
#include "AliResonanceFits.h"
#endif

#include <iostream>
using std::cout;
using std::endl;

#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TMinuit.h>
#include <THn.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>

ClassImp(AliResonanceFits)

const Float_t gkDEFAULT_FLOAT_INIT = -999;
TH1* AliResonanceFits::fTempHistSignal = 0x0;
TH1* AliResonanceFits::fTempHistBkgnd = 0x0;
Int_t AliResonanceFits::fNdf = 0;
Float_t AliResonanceFits::fFitRange[2] = {0.0, 20.0};
Float_t AliResonanceFits::fExclusionRange[2] = {0.0,0.0};
Float_t AliResonanceFits::fPtRange[2] = {gkDEFAULT_FLOAT_INIT, gkDEFAULT_FLOAT_INIT};
Bool_t AliResonanceFits::fUse2DMatching = kFALSE;

//_______________________________________________________________________________
AliResonanceFits::AliResonanceFits() :
  fSEOS(0x0),
  fSELSleg1(0x0),
  fSELSleg2(0x0),
  fMEOS(0x0),
  fMELSleg1(0x0),
  fMELSleg2(0x0),
  fEffVsPtCent(0x0),
  fWeightVsCent(0x0),
  fEventVsCent(0x0),
  fSignalMCshape(0x0),
  fCentralitySelection(kFALSE),
  fVertexSelection(kFALSE),
  fEPSelection(kFALSE),
  fPtSelection(kFALSE),
  fPlottingOption(0),
  fBkgMethod(1),
  fMatchingOption(1),
  fWeightedAveragePower(2.0),
  fMinuitFitOption(1.0),
  fFixScale(kFALSE),
  fMEMatchOption(2),
  fLSmethod(1),
  fIsProcessed(kFALSE),
  //fNdf(0),
//  fTempHistSignal(0x0),
//  fTempHistBkgnd(0x0),
  fHistBkgSubtracted(0x0),
  fHistSEOS(0x0),
  fHistSELSbkg(0x0),
  fHistMEbkg(0x0)
  
{
  //
  // Default constructor
  //
  for(Int_t i=0; i<kNVariables; ++i) {
    fUsedVars[i] = kFALSE;
    fVarIndex[i] = -1;
  }
  //fFitRange[0]=0.0; fFitRange[1]=10.0;
  //fExclusionRange[0]=0.0; fExclusionRange[1]=0.0;
  fMassRange[0]=0.0; fMassRange[1]=20.0;
  fSignalRange[0]=0.0; fSignalRange[1]=20.0;
  fCentralityRange[0]=gkDEFAULT_FLOAT_INIT; fCentralityRange[1]=gkDEFAULT_FLOAT_INIT;
  fVertexZRange[0]=gkDEFAULT_FLOAT_INIT; fVertexZRange[1]=gkDEFAULT_FLOAT_INIT; 
  fEP2Range[0]=gkDEFAULT_FLOAT_INIT; fEP2Range[1]=gkDEFAULT_FLOAT_INIT;
  //fPtRange[0]=gkDEFAULT_FLOAT_INIT; fPtRange[1]=gkDEFAULT_FLOAT_INIT;
  for(Int_t i=0; i<kNFitValues; ++i) fFitValues[i] = 0.0;
}

//_______________________________________________________________________________
AliResonanceFits::~AliResonanceFits()
{
  //
  // De-constructor
  //
  //if(fTempHistSignal) delete fTempHistSignal;
  //if(fTempHistBkgnd) delete fTempHistBkgnd;
  if(fHistBkgSubtracted) delete fHistBkgSubtracted;
  if(fHistSEOS) delete fHistSEOS;
  if(fHistSELSbkg) delete fHistSELSbkg;
  if(fHistMEbkg) delete fHistMEbkg;
  if(fSEOS) delete fSEOS;
  if(fMEOS) delete fMEOS;
  if(fSELSleg1) delete fSELSleg1;
  if(fSELSleg2) delete fSELSleg2;
  if(fMELSleg1) delete fMELSleg1;
  if(fMELSleg2) delete fMELSleg2;
  if(fEventVsCent) delete fEventVsCent;
  //if(fEffVsPtCent) delete fEffVsPtCent;
  //if(fWeightVsCent) delete fWeightVsCent;
  //if(fSignalMCshape) delete fSignalMCshape;
}

//_______________________________________________________________________________
void AliResonanceFits::SetHistograms(THnF* seos, THnF* meos, 
		                     THnF* selsLeg1 /*=0x0*/, THnF* selsLeg2 /*=0x0*/, 
		                     THnF* melsLeg1 /*=0x0*/, THnF* melsLeg2 /*=0x0*/) 
{
  //
  // Set the multi-dim histograms
  //
  fSEOS = seos; fMEOS = meos;
  fSELSleg1 = selsLeg1; fSELSleg2 = selsLeg2;
  fMELSleg1 = melsLeg1; fMELSleg2 = melsLeg2;
  SetDefaultRanges();
}

//_______________________________________________________________________________
void AliResonanceFits::SetVars(Int_t* vars)
{
  //
  // Set the mapping of the various observables in the THn histograms
  //
  for(Int_t i=0; i<kNVariables;++i) {
    if(vars[i]<0) continue;   // this variable is not used
    fUsedVars[i] = kTRUE;
    fVarIndex[i] = vars[i];
  }
  SetDefaultRanges();
}

//____________________________________________________________________________________
void AliResonanceFits::SetDefaultRanges()
{
  //
  // Set default ranges using the input histograms
  //
  if(!fSEOS) return;
  Int_t nUsedVars=0;
  for(Int_t i=0; i<kNVariables;++i) {
    if(fUsedVars[i]) ++nUsedVars;
  }
  if(nUsedVars==0) return;
  
  Double_t minMass = fSEOS->GetAxis(fVarIndex[kMass])->GetXmin();
  Double_t maxMass = fSEOS->GetAxis(fVarIndex[kMass])->GetXmax();
  fFitRange[0] = minMass;
  fFitRange[1] = maxMass;
  fExclusionRange[0] = minMass + 0.25*(maxMass-minMass);
  fExclusionRange[1] = minMass + 0.75*(maxMass-minMass);
  fMassRange[0] = minMass;
  fMassRange[1] = maxMass;
  fSignalRange[0] = minMass + 0.25*(maxMass-minMass);
  fSignalRange[1] = minMass + 0.75*(maxMass-minMass);
  if(fUsedVars[kCentrality]) {
    fCentralityRange[0] = fSEOS->GetAxis(fVarIndex[kCentrality])->GetXmin();
    fCentralityRange[1] = fSEOS->GetAxis(fVarIndex[kCentrality])->GetXmax();
  }
  if(fUsedVars[kVertexZ]) {
    fVertexZRange[0] = fSEOS->GetAxis(fVarIndex[kVertexZ])->GetXmin();
    fVertexZRange[1] = fSEOS->GetAxis(fVarIndex[kVertexZ])->GetXmax();
  }
  if(fUsedVars[kEP2]) {
    fEP2Range[0] = fSEOS->GetAxis(fVarIndex[kEP2])->GetXmin();
    fEP2Range[1] = fSEOS->GetAxis(fVarIndex[kEP2])->GetXmax();
  }
  if(fUsedVars[kPt]) {
    fPtRange[0] = fSEOS->GetAxis(fVarIndex[kPt])->GetXmin();
    fPtRange[1] = fSEOS->GetAxis(fVarIndex[kPt])->GetXmax();
  }
}

//____________________________________________________________________________________
void AliResonanceFits::Reset()
{
  //
  // Reset everything so that another independent fit can be done
  //
  for(Int_t i=0; i<kNVariables; ++i) {
    fUsedVars[i] = kFALSE;
    fVarIndex[i] = -1;
  }
  fFitRange[0]=0.0; fFitRange[1]=10.0;
  fExclusionRange[0]=0.0; fExclusionRange[1]=0.0;
  fSignalRange[0]=0.0;fSignalRange[1]=10.0;
  fCentralityRange[0]=gkDEFAULT_FLOAT_INIT; fCentralityRange[1]=gkDEFAULT_FLOAT_INIT;
  fCentralitySelection = kFALSE;
  fVertexZRange[0]=gkDEFAULT_FLOAT_INIT; fVertexZRange[1]=gkDEFAULT_FLOAT_INIT; 
  fVertexSelection = kFALSE;
  fEP2Range[0]=gkDEFAULT_FLOAT_INIT; fEP2Range[1]=gkDEFAULT_FLOAT_INIT;
  fEPSelection = kFALSE;
  fPtRange[0]=gkDEFAULT_FLOAT_INIT; fPtRange[1]=gkDEFAULT_FLOAT_INIT;
  fPtSelection = kFALSE;
  fIsProcessed = kFALSE;
  fNdf = 0;
  fPlottingOption = 0;
  fUse2DMatching = kFALSE;
  fWeightedAveragePower = 2.0;
  fMinuitFitOption = 1.0;
  fFixScale = kFALSE;
  fMatchingOption = 1;
  fMEMatchOption = 2;
  fLSmethod = 1;
  fBkgMethod = 1;
  fTempHistSignal = 0x0;
  fTempHistBkgnd = 0x0;
  for(Int_t i=0; i<kNFitValues; ++i) fFitValues[i] = 0.0;
  if(fHistBkgSubtracted) {delete fHistBkgSubtracted; fHistBkgSubtracted=0x0;}
  if(fHistSEOS) {delete fHistSEOS; fHistSEOS=0x0;}
  if(fHistSELSbkg) {delete fHistSELSbkg; fHistSELSbkg=0x0;}
  if(fHistMEbkg) {delete fHistMEbkg; fHistMEbkg=0x0;}
  if(fSEOS) delete fSEOS;
  if(fMEOS) delete fMEOS;
  if(fSELSleg1) delete fSELSleg1;
  if(fSELSleg2) delete fSELSleg2;
  if(fMELSleg1) delete fMELSleg1;
  if(fMELSleg2) delete fMELSleg2;
  if(fEventVsCent) delete fEventVsCent;
  //if(fEffVsPtCent) delete fEffVsPtCent;
  //if(fWeightVsCent) delete fWeightVsCent;
  //if(fSignalMCshape) delete fSignalMCshape;
  fSEOS = 0x0; fMEOS = 0x0;
  fSELSleg1 = 0x0; fSELSleg2 = 0x0;
  fMELSleg1 = 0x0; fMELSleg2 = 0x0;
  fEffVsPtCent = 0x0;
  fWeightVsCent = 0x0;
  fEventVsCent = 0x0;
  fSignalMCshape = 0x0;
}

//____________________________________________________________________________________
void AliResonanceFits::Process(Int_t bkgMethod, Int_t matchOption, Float_t waPower, Float_t minuitOption, 
	                       Float_t fixBkgScale, Int_t meMatchBkg, Int_t lsMethod) {
  //
  // Set the fitting options and call the main process function
  //
  fBkgMethod = bkgMethod;
  fMatchingOption = matchOption;
  fWeightedAveragePower = waPower;
  fMinuitFitOption = minuitOption;
  if(fixBkgScale>0.0) {fFixScale=kTRUE; fFitValues[kBkgScale]=fixBkgScale;}
  fMEMatchOption = meMatchBkg;
  fLSmethod = lsMethod;
  Process();
}

//____________________________________________________________________________________
void AliResonanceFits::Process() {
  //
  // Main steering function
  //
  if(!fSEOS) {
    cout << "Error: AliResonanceFits needs the SE-OS histogram but its pointer is null !" << endl;
    return;
  }
  if(!fMEOS) {
    cout << "Error: AliResonanceFits needs the ME-OS histogram but its pointer is null !" << endl;
    return;
  }
  
  if((fBkgMethod==1 && fMEMatchOption==2) || fBkgMethod==2) {  
    if(!fSELSleg1) {
      cout << "Error: AliResonanceFits needs the SE-LS leg1 histogram but its pointer is null !" << endl;
      return;
    }
    if(!fSELSleg2) {
      cout << "Error: AliResonanceFits needs the SE-LS leg2 histogram but its pointer is null !" << endl;
      return;
    }
    if(!fMELSleg1) {
      cout << "Error: AliResonanceFits needs the ME-LS leg1 histogram but its pointer is null !" << endl;
      return;
    }
    if(!fMELSleg2) {
      cout << "Error: AliResonanceFits needs the ME-LS leg2 histogram but its pointer is null !" << endl;
      return;
    }
  }
  
  // Create the signal and bkg histograms
  if(fHistSEOS) delete fHistSEOS;
  if(fHistSELSbkg) delete fHistSELSbkg;
  if(fHistMEbkg) delete fHistMEbkg;
  
  fSEOS->GetAxis(fVarIndex[kMass])->SetRangeUser(fMassRange[0]+0.00001, fMassRange[1]-0.00001);
  if(fUse2DMatching && fPtSelection)
     fSEOS->GetAxis(fVarIndex[kPt])->SetRangeUser(fPtRange[0]+0.001, fPtRange[1]-0.001);
  if(fUse2DMatching) {
    TH2D* proj = fSEOS->Projection(fVarIndex[kPt], fVarIndex[kMass]);
    proj->SetName("tempProj");
    fHistSEOS = (TH2D*)proj->Clone("fHistSEOS");
    delete proj;
    fHistSEOS->Reset();
    fHistSELSbkg = (TH2D*)fHistSEOS->Clone("fHistSELSbkg");
    fHistMEbkg = (TH2D*)fHistSEOS->Clone("fHistMEbkg");
  }
  else {
     TH1D* proj = fSEOS->Projection(fVarIndex[kMass]);
     fHistSEOS = (TH1D*)proj->Clone("fHistSEOS");
     delete proj;
     fHistSEOS->Reset();
     fHistSELSbkg = (TH1D*)fHistSEOS->Clone("fHistSELSbkg");
     fHistMEbkg = (TH1D*)fHistSEOS->Clone("fHistMEbkg");
  }
  
  /*
  TAxis* mAxis = (TAxis*)fSEOS->GetAxis(fVarIndex[kMass])->Clone(Form("mAxis%.6f",gRandom->Rndm()));
  Double_t minMass = mAxis->GetXmin(); Double_t maxMass = mAxis->GetXmax();
  if(fMassRange[0]>minMass) minMass = fMassRange[0]+0.00001;
  if(fMassRange[1]<maxMass) maxMass = fMassRange[1]-0.00001;
  mAxis->SetRangeUser(minMass, maxMass);
  
  if(fUse2DMatching) {
     cout << "AliResonanceFits::Process() fUse2DMatching" << endl;
     TAxis* ptAxis = (TAxis*)fSEOS->GetAxis(fVarIndex[kPt])->Clone(Form("ptAxis%.6f",gRandom->Rndm()));
     Double_t minPt = fPtRange[0]; Double_t maxPt = fPtRange[1];
     if(fPtSelection) {
        cout << "AliResonanceFits::Process() fPtSelection" << endl;
        //ptIdx = fVarIndex[kPt];
        minPt = ptAxis->GetXmin(); maxPt = ptAxis->GetXmax();
        if(fPtRange[0]>minPt) minPt = fPtRange[0]+0.001;
        if(fPtRange[1]<maxPt) maxPt = fPtRange[1]-0.001;
        ptAxis->SetRangeUser(minPt, maxPt);
     }
     if(mAxis->IsVariableBinSize() && ptAxis->IsVariableBinSize()) {
       fHistSEOS = new TH2D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                                                                                       ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
       fHistSELSbkg = new TH2D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                            ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
       fHistMEbkg = new TH2D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                            ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
    }
    if(!mAxis->IsVariableBinSize() && ptAxis->IsVariableBinSize()) {
       fHistSEOS = new TH2D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                            ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
       fHistSELSbkg = new TH2D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                               ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
       fHistMEbkg = new TH2D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                             ptAxis->GetNbins(), ptAxis->GetXbins()->GetArray());
    }
    if(mAxis->IsVariableBinSize() && !ptAxis->IsVariableBinSize()) {
       fHistSEOS = new TH2D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                            ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
       fHistSELSbkg = new TH2D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                               ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
       fHistMEbkg = new TH2D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray(),
                             ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
    }
    if(!mAxis->IsVariableBinSize() && !ptAxis->IsVariableBinSize()) {
       fHistSEOS = new TH2D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                            ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
       fHistSELSbkg = new TH2D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                               ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
       fHistMEbkg = new TH2D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax(),
                             ptAxis->GetNbins(), ptAxis->GetXmin(), ptAxis->GetXmax());
    }
  }  // end if(fUse2DMatching)
  else {
     cout << "AliResonanceFits::Process() fUse2DMatching NO" << endl;
     if(mAxis->IsVariableBinSize()) {
        cout << "AliResonanceFits::Process() mAxis is variable bin size" << endl;
        fHistSEOS = new TH1D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXbins()->GetArray());
        fHistSELSbkg = new TH1D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray());
        fHistMEbkg = new TH1D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXbins()->GetArray());
     }
     else {
        cout << "AliResonanceFits::Process() mAxis is not variable bin size" << endl;
       fHistSEOS = new TH1D("fHistSEOS", "fHistSEOS", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax());
       fHistSELSbkg = new TH1D("fHistSELSbkg", "fHistSELSbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax());
       fHistMEbkg = new TH1D("fHistMEbkg", "fHistMEbkg", mAxis->GetNbins(), mAxis->GetXmin(), mAxis->GetXmax());
     }
  }
  */
  
  //MakeSEOS(); fHistSEOS->Draw();
  //MakeLSbkg(); fHistSELSbkg->Draw();
  //MakeMEbkg(); fHistMEbkg->Draw();
  MakeSEOS();
  if(fBkgMethod==1) {MakeMEbkg(); ExtractSignal(fHistSEOS, fHistMEbkg);}
  if(fBkgMethod==2) {MakeLSbkg(); ExtractSignal(fHistSEOS, fHistSELSbkg);}
  fIsProcessed = kTRUE;
}

//____________________________________________________________________________________
Double_t AliResonanceFits::Chi2(TH1* signal, TH1* bkg, Double_t scale, Double_t scaleError) {
  //
  // Compute the chi2 for the difference between the signal and scaled background
  // Assume signal and background uncertainties are uncorrelated
  //
  // NOTE:
  // In case of 2D matching, the signal and bkg are TH2D histograms with mass on X-axis and pt on Y-axis
  
  Float_t chi2 = 0.0;
  fNdf = 0;
  //cout << "Chi2::scale = " << scale << endl;
  //cout << "Chi2::scaleError = " << scaleError << endl;
  for(Int_t mbin=1; mbin<=signal->GetXaxis()->GetNbins(); ++mbin) {
    Float_t mass = signal->GetXaxis()->GetBinCenter(mbin);
    if(mass<fFitRange[0] || mass>fFitRange[1]) continue;
    if(mass>fExclusionRange[0] && mass<fExclusionRange[1]) continue;
    
    for(Int_t ptbin=1; ptbin<=(fUse2DMatching ? signal->GetYaxis()->GetNbins() : 1); ++ptbin) {
      Float_t pt = (fUse2DMatching ? signal->GetYaxis()->GetBinCenter(ptbin) : 0.0);
      if(fUse2DMatching && (pt<fPtRange[0] || pt>fPtRange[1])) continue;      
      
      Float_t sig = (fUse2DMatching ? signal->GetBinContent(mbin,ptbin) : signal->GetBinContent(mbin));
      if(sig<=0.0) continue;
      Float_t bkgVal = (fUse2DMatching ? bkg->GetBinContent(bkg->FindBin(mass,pt)) : bkg->GetBinContent(bkg->FindBin(mass)));
      if(bkgVal<=0.0) continue;
      Float_t bkgValErr = (fUse2DMatching ? bkg->GetBinError(bkg->FindBin(mass,pt)) : bkg->GetBinError(bkg->FindBin(mass)));
      Float_t sigErr = (fUse2DMatching ? signal->GetBinError(mbin,ptbin) : signal->GetBinError(mbin));
      Float_t err = 0.0;
      if(scale>0.0)
        err = bkgVal*scale*TMath::Sqrt(scaleError*scaleError/scale/scale + bkgValErr*bkgValErr/bkgVal/bkgVal);
      err = sigErr*sigErr+err*err;
      
      chi2 += (err>0.0 ? (sig-scale*bkgVal)*(sig-scale*bkgVal)/err : 0.0);
      
      //cout << "mass / pt :: " << mass << "/" << pt << endl;
      //cout << "sig=" << sig << "; bkg=" << bkgVal << "; sigErr=" << sigErr << "; bkgValErr=" << bkgValErr << "; diff2 = " << (sig-bkgVal)*(sig-bkgVal) << ";  err2 = " << err << "; chi2=" << chi2 << endl;
      
      ++fNdf;
    }
  }  // end loop over mass bins
  //cout << "Chi2::ndf  = " << fNdf << endl;
  //cout << "Chi2::chi2 = " << chi2/Float_t(fNdf) << endl;
  return (fNdf>0 ? chi2/Float_t(fNdf) : 1000.);
}

//____________________________________________________________________________________
void AliResonanceFits::Fcn(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // chi2 function used as interface with minuit
  // par[0] - scale factor for the background histogram
  //
  Double_t chi2 = Chi2(fTempHistSignal, fTempHistBkgnd, par[0]);
  f = chi2;  
}

//____________________________________________________________________________________
void AliResonanceFits::ComputeWeightedScale(TH1* signal, TH1* bkg) {
  //
  // Calculate the ME scale by weighted average
  //
  // NOTE:
  // In case of 2D matching, the signal and bkg are TH2D histograms with mass on X-axis and pt on Y-axis
  
  Double_t sweights = 0.0; Double_t avWeights = 0.0; Double_t nMassBins=0;
  Double_t serror = 0.0;
  for(Int_t im=1; im<=signal->GetXaxis()->GetNbins(); ++im) {
    Double_t m = signal->GetXaxis()->GetBinCenter(im);
    if(m<fFitRange[0]) continue; 
    if(m>fFitRange[1]) continue;
    if(m>fExclusionRange[0] && m<fExclusionRange[1]) continue;
    //cout << "m :: " << m << endl;
    
    for(Int_t ipt=1; ipt<=(fUse2DMatching ? signal->GetYaxis()->GetNbins() : 1); ++ipt) {
      Float_t pt = (fUse2DMatching ? signal->GetYaxis()->GetBinCenter(ipt) : 0.0);
      if(fUse2DMatching && (pt<fPtRange[0] || pt>fPtRange[1])) continue;      
      
      //cout << "pt :: " << pt << endl;
      
      Double_t s = (fUse2DMatching ? signal->GetBinContent(im,ipt) : signal->GetBinContent(im)); 
      Double_t sErr = (fUse2DMatching ? signal->GetBinError(im,ipt) : signal->GetBinError(im));
      Double_t b = (fUse2DMatching ? bkg->GetBinContent(im,ipt) : bkg->GetBinContent(im)); 
      Double_t bErr = (fUse2DMatching ? bkg->GetBinError(im,ipt) : bkg->GetBinError(im));
      if(s<=0.0 || b<=0.0) continue;
      Double_t sOverB = s / b;
      Double_t sOverBerr = sOverB * TMath::Sqrt(TMath::Power(sErr/s,2.0) + TMath::Power(bErr/b,2.0));
    
      // weighting using S/B error
      sweights  += 1.0/TMath::Power(sOverBerr, fWeightedAveragePower);
      avWeights += sOverB/TMath::Power(sOverBerr, fWeightedAveragePower);
      serror += TMath::Power(sOverBerr, 2.0-2.0*fWeightedAveragePower);
    
      nMassBins += 1.0;
      /*cout << "mass/s/serr/b/berr/sOverB/sOverBerr/sweights/avWeights/nMassBins: " 
       *        << m << " / " << sigEntries << " / " << sigEntriesErr << " / " << bkgEntries << " / " << bkgEntriesErr
       *        << " / " << sOverB << " / " << sOverBerr << " / " << sweights << " / " << avWeights << " / " << nMassBins << endl;*/
    }
  }
  if(sweights>0.0) avWeights /= sweights;
  fFitValues[kBkgScale] = avWeights;  
  if(sweights>0.0) fFitValues[kBkgScaleErr] = TMath::Sqrt(serror)/sweights;
  else fFitValues[kBkgScaleErr] = TMath::Sqrt(serror);
}

//____________________________________________________________________________________
void AliResonanceFits::FitScale(TH1* signal, TH1* bkg) {
  //
  // Compute the bkg scaling by fitting
  //
  /*if((signal->GetDimension()==1 && signal->Integral()<5) || 
     (signal->GetDimension()==2 && signal->Integral()<25)) {
    cout << "GOTO ComputeWeightedScale()" << endl;
    ComputeWeightedScale(signal, bkg);
    return;
  }*/
  fTempHistSignal = signal;
  fTempHistBkgnd = bkg;
    
  TMinuit* minuit=new TMinuit(1);
  minuit->SetFCN(Fcn);
  
  Double_t arglist[2];
  Int_t ierflg=0;
  arglist[0] = fMinuitFitOption;   // 1-chisquare minimization; 0.5-log likelihood
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  
  // Set starting values and step sizes for parameters
  Double_t vstart[1];
  Double_t step[1];
  vstart[0] = 1;  
  step[0] = 0.1; 
  
  minuit->mnparm(0, "scale", vstart[0], step[0], 0.0, 1000000.0, ierflg);
  if(fFixScale) minuit->FixParameter(fFitValues[kBkgScale]);
  
  arglist[0] = 500;
  arglist[1] = 1.;
  
  minuit->SetMaxIterations(10000);
  minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
  minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  Double_t scale=0.0, error=0.0;
  minuit->GetParameter(0, scale, error);
  //cout << "ME scale factor = " << scale << " +/- " << error << endl;
  delete minuit;
  fFitValues[kBkgScale] = scale;
  fFitValues[kBkgScaleErr] = error;
}

//____________________________________________________________________________________
void AliResonanceFits::ComputeEntryScale(TH1* signal, TH1* bkg) {
  //
  // Compute the bkg scaling based on the number of entries
  //
  // NOTE:
  // In case of 2D matching, the signal and bkg are TH2D histograms with mass on X-axis and pt on Y-axis
  
  Double_t entriesSE = 0.0; Double_t entriesME = 0.0;
  Double_t seError = 0.0; Double_t meError = 0.0;
  //cout << "AliResonanceFits::ComputeEntryScale signal50=" << signal->GetBinContent(50) << endl;
  //cout << "AliResonanceFits::ComputeEntryScale bkg50=" << bkg->GetBinContent(50) << endl;
  for(Int_t im=1;im<signal->GetXaxis()->GetNbins();++im) {
    Double_t mass = signal->GetXaxis()->GetBinCenter(im);
    cout << "mass " << mass << endl;
    if(mass<fFitRange[0]) continue;
    if(mass>fFitRange[1]) continue;
    if(mass>fExclusionRange[0] && mass<fExclusionRange[1]) continue;
    
    for(Int_t ipt=1;ipt<=(fUse2DMatching ? signal->GetYaxis()->GetNbins() : 1);++ipt) {
      Float_t pt = (fUse2DMatching ? signal->GetYaxis()->GetBinCenter(ipt) : 0.0);
      if(fUse2DMatching && (pt<fPtRange[0] || pt>fPtRange[1])) continue;      
      
      entriesSE += (fUse2DMatching ? signal->GetBinContent(im,ipt) : signal->GetBinContent(im));
      seError += (fUse2DMatching ? TMath::Power(signal->GetBinError(im,ipt),2.0) : TMath::Power(signal->GetBinError(im),2.0));
      entriesME += (fUse2DMatching ? bkg->GetBinContent(im,ipt) : bkg->GetBinContent(im));
      meError += (fUse2DMatching ? TMath::Power(signal->GetBinError(im,ipt),2.0) : TMath::Power(signal->GetBinError(im),2.0));
      cout << "ipt / se / me : " << ipt << "/" << entriesSE << "/" << entriesME << endl;
    }
  }
  //cout << "AliResonanceFits::ComputeEntryScale entriesSE=" << entriesSE << endl;
  //cout << "AliResonanceFits::ComputeEntryScale entriesME=" << entriesME << endl;
  seError = TMath::Sqrt(seError); seError /= (entriesSE>1.0e-6 ? entriesSE : 1.);
  meError = TMath::Sqrt(meError); meError /= (entriesME>1.0e-6 ? entriesME : 1.);
  fFitValues[kBkgScale] = (entriesSE>1.0e-6 && entriesME>1.0e-6 ? entriesSE/entriesME : 0.0);
  fFitValues[kBkgScaleErr] = fFitValues[kBkgScale]*TMath::Sqrt(seError*seError+meError*meError);
}

//____________________________________________________________________________________
void AliResonanceFits::ExtractSignal(TH1* signal, TH1* bkg, Bool_t fixScale /*=kFALSE*/) {
  //
  // Scale the bkg by a constant and then subtract it to obtain the signal
  //
  fFixScale = fixScale;
  
  if(fHistBkgSubtracted) {
    delete fHistBkgSubtracted;
    fHistBkgSubtracted = 0x0;
  }
  if(fUse2DMatching) 
    fHistBkgSubtracted = new TH2D(Form("fHistBkgSubtracted_%s", signal->GetName()),"",
                                  signal->GetXaxis()->GetNbins(), 
				  signal->GetXaxis()->GetXmin(), signal->GetXaxis()->GetXmax(),
				  signal->GetYaxis()->GetNbins(), signal->GetYaxis()->GetXbins()->GetArray());
  else  
    fHistBkgSubtracted = new TH1D(Form("fHistBkgSubtracted_%s", signal->GetName()),"",
                                  signal->GetXaxis()->GetNbins(), 
	  	 		  signal->GetXaxis()->GetXmin(), signal->GetXaxis()->GetXmax());
  fHistBkgSubtracted->Sumw2();
  
  if(!fFixScale) {
    if(fMatchingOption==1) ComputeWeightedScale(signal, bkg);
    if(fMatchingOption==2) FitScale(signal, bkg);
    if(fMatchingOption==3) ComputeEntryScale(signal, bkg);
  }
    
  Float_t oldFitRanges[2] = {fFitRange[0], fFitRange[1]};
  fFitRange[0] = fMassRange[0]; fFitRange[1] = fMassRange[1];
  fFitValues[kChisq] = Chi2(signal, bkg, fFitValues[kBkgScale], fFitValues[kBkgScaleErr]);
  cout << "chisq/ndf = " << fFitValues[kChisq]*fNdf << "/" << fNdf << endl;
  fFitRange[0] = oldFitRanges[0]; fFitRange[1] = oldFitRanges[1];
  
  fFitValues[kSig] = 0.0; fFitValues[kSigErr] = 0.0;
  fFitValues[kBkg] = 0.0; fFitValues[kBkgErr] = 0.0;
  fFitValues[kSplusB] = 0.0; fFitValues[kSplusBerr] = 0.0;
  for(Int_t im=1; im<=signal->GetXaxis()->GetNbins(); ++im) {
    Double_t m = signal->GetXaxis()->GetBinCenter(im);
    for(Int_t ipt=1; ipt<=(fUse2DMatching ? signal->GetYaxis()->GetNbins() : 1); ++ipt) {
      Double_t pt   = (fUse2DMatching ? signal->GetYaxis()->GetBinCenter(ipt) : 0.0);
      Double_t s    = (fUse2DMatching ? signal->GetBinContent(im,ipt) : signal->GetBinContent(im));
      if(s<=0.0) continue;
      Double_t sErr = (fUse2DMatching ? signal->GetBinError(im,ipt) : signal->GetBinError(im));
      Double_t b    = fFitValues[kBkgScale]*(fUse2DMatching ? bkg->GetBinContent(bkg->FindBin(m,pt)) : bkg->GetBinContent(bkg->FindBin(m)));
      if(b<=0.0) continue;
      Double_t bErr = fFitValues[kBkgScale]*(fUse2DMatching ? bkg->GetBinError(bkg->FindBin(m,pt)) : bkg->GetBinError(bkg->FindBin(m)));
               bErr = b*fFitValues[kBkgScale]*TMath::Sqrt(fFitValues[kBkgScaleErr]*fFitValues[kBkgScaleErr]/fFitValues[kBkgScale]/fFitValues[kBkgScale] +
                                                          bErr*bErr/b/b);				    
      if(fUse2DMatching)
	fHistBkgSubtracted->SetBinContent(im, ipt, s-b);
      else
	fHistBkgSubtracted->SetBinContent(im, s-b);
      if(fUse2DMatching)
	fHistBkgSubtracted->SetBinError(im, ipt, TMath::Sqrt(sErr*sErr+bErr*bErr));
      else
	fHistBkgSubtracted->SetBinError(im, TMath::Sqrt(sErr*sErr+bErr*bErr));
      
      if(m<fSignalRange[0] || m>fSignalRange[1]) continue;
      fFitValues[kSig] += (s-b);
      fFitValues[kSigErr] += (sErr*sErr+bErr*bErr);
      fFitValues[kBkg] += b;
      fFitValues[kBkgErr] += (bErr*bErr);
      fFitValues[kSplusB] += s;
      fFitValues[kSplusBerr] += (sErr*sErr);
      cout << "bin (mass,pt) " << m << ", " << pt << ";"<< endl;
      cout << "S+B : " << s << " +/- " << sErr << endl;
      cout << "B   : " << b << " +/- " << bErr << endl;
      cout << "S-B : " << (s-b) << " +/- " << TMath::Sqrt(sErr*sErr+bErr*bErr) << endl;
    }
  }
  fFitValues[kSigErr] = TMath::Sqrt(fFitValues[kSigErr]);
  fFitValues[kBkgErr] = TMath::Sqrt(fFitValues[kBkgErr]);  
  fFitValues[kSplusBerr] = TMath::Sqrt(fFitValues[kSplusB]);
  
  Double_t sErr = fFitValues[kSigErr];
  Double_t bErr = fFitValues[kBkgErr];
  Double_t s = fFitValues[kSig]; Double_t b = fFitValues[kBkg];
  fFitValues[kSoverB] = (b>0.0 ? s/b : 0.0);
  fFitValues[kSoverBerr] = (b>0.0 ? TMath::Sqrt((sErr*sErr+fFitValues[kSoverB]*fFitValues[kSoverB]*bErr*bErr)/b/b) : 0.0);
  fFitValues[kSignif] = ((s+b>0.001) ? (s/(TMath::Sqrt(s+b))) : 0.0);
  fFitValues[kSignifErr] = ((s+b>0.001 && s>0.001) ? (fFitValues[kSignif]*TMath::Sqrt(bErr*bErr+TMath::Power(sErr*(s+2.0*b)/s,2))/2.0/(s+b)) : 0.0);
  
  bkg->Scale(fFitValues[kBkgScale]);
}


//_____________________________________________________________________________________________
void AliResonanceFits::MakeSEOS() {
  //
  // Build the invariant mass spectrum for signal (SEOS)
  //
   Int_t zIdx = fVarIndex[kVertexZ];
  Double_t minZ = fVertexZRange[0]; Double_t maxZ = fVertexZRange[1];
  if(fUsedVars[kVertexZ] && fVertexSelection) {
    //zIdx = fVarIndex[kVertexZ];
    minZ = fSEOS->GetAxis(zIdx)->GetXmin(); maxZ = fSEOS->GetAxis(zIdx)->GetXmax();
    if(fVertexZRange[0]>minZ) minZ = fVertexZRange[0]+0.001;
    if(fVertexZRange[1]<maxZ) maxZ = fVertexZRange[1]-0.001;
    fSEOS->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
  }
  Int_t epIdx = fVarIndex[kEP2];
  Double_t minEP = fEP2Range[0]; Double_t maxEP = fEP2Range[1];
  if(fUsedVars[kEP2] && fEPSelection) {
    //epIdx = fVarIndex[kEP2];
    minEP = fSEOS->GetAxis(epIdx)->GetXmin(); maxEP = fSEOS->GetAxis(epIdx)->GetXmax();
    if(fEP2Range[0]>minEP) minEP = fEP2Range[0]+0.001;
    if(fEP2Range[1]<maxEP) maxEP = fEP2Range[1]-0.001;
    fSEOS->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
  }
  Int_t ptIdx = fVarIndex[kPt];
  Double_t minPt = fPtRange[0]; Double_t maxPt = fPtRange[1];
  if(fUsedVars[kPt] && fPtSelection) {
    //ptIdx = fVarIndex[kPt];
    minPt = fSEOS->GetAxis(ptIdx)->GetXmin(); maxPt = fSEOS->GetAxis(ptIdx)->GetXmax();
    if(fPtRange[0]>minPt) minPt = fPtRange[0]+0.001;
    if(fPtRange[1]<maxPt) maxPt = fPtRange[1]-0.001;
    fSEOS->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
  }
  Int_t centIdx = fVarIndex[kCentrality];
  Double_t minCent = fCentralityRange[0]; Double_t maxCent = fCentralityRange[1];
  if(fUsedVars[kCentrality] && fCentralitySelection) {
    //centIdx = fVarIndex[kCentrality];
    minCent = fSEOS->GetAxis(centIdx)->GetXmin(); maxCent = fSEOS->GetAxis(centIdx)->GetXmax();
    if(fCentralityRange[0]>minCent) minCent = fCentralityRange[0]+0.001;
    if(fCentralityRange[1]<maxCent) maxCent = fCentralityRange[1]-0.001;
    fSEOS->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
  }
  Int_t mIdx = fVarIndex[kMass];
  Double_t minMass = fSEOS->GetAxis(mIdx)->GetXmin(); Double_t maxMass = fSEOS->GetAxis(mIdx)->GetXmax();
  if(fMassRange[0]>minMass) minMass = fMassRange[0]+0.00001;
  if(fMassRange[1]<maxMass) maxMass = fMassRange[1]-0.00001;
  fSEOS->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
        
  if(fUsedVars[kCentrality] && (fEffVsPtCent || fWeightVsCent)) {
    for(Int_t icent=1;icent<=fSEOS->GetAxis(fVarIndex[kCentrality])->GetNbins(); ++icent) {
      if(fSEOS->GetAxis(fVarIndex[kCentrality])->GetBinCenter(icent)<minCent) continue;
      if(fSEOS->GetAxis(fVarIndex[kCentrality])->GetBinCenter(icent)>maxCent) continue;
      	
      Double_t evWeight = 1.0;
      if(fWeightVsCent) evWeight = fWeightVsCent->GetBinContent(icent);    //fWeightVsCent has the same binning as the THnF
      
      fSEOS->GetAxis(fVarIndex[kCentrality])->SetRange(icent, icent);
      TH1* tempProj = 0x0;
      if(fUse2DMatching) tempProj = (TH2D*)fSEOS->Projection(fVarIndex[kPt],fVarIndex[kMass]);
      else tempProj = (TH1D*)fSEOS->Projection(fVarIndex[kMass]);
            
      // Apply corrections or weights if applicable		  
      for(Int_t ipt=1; ipt<=(fUse2DMatching ? tempProj->GetYaxis()->GetNbins() : 1); ++ipt) {
	Float_t pt = (fUse2DMatching ? tempProj->GetYaxis()->GetBinCenter(ipt) : 0.0);
	Int_t ptBin = ((fUse2DMatching && fEffVsPtCent) ? fEffVsPtCent->GetYaxis()->FindBin(pt) : 1);
	Int_t centBin = (fEffVsPtCent ? fEffVsPtCent->GetXaxis()->FindBin(fSEOS->GetAxis(fVarIndex[kCentrality])->GetBinCenter(icent)) : 0);
	Double_t eff = 1.0;
        if(fEffVsPtCent) {
	  if(centBin==0) centBin=1;
          if(centBin==fEffVsPtCent->GetXaxis()->GetNbins()+1) centBin = fEffVsPtCent->GetXaxis()->GetNbins();
	  eff = fEffVsPtCent->GetBinContent(centBin, ptBin);
	}
        for(Int_t im=1;im<=tempProj->GetXaxis()->GetNbins();++im) {
          Double_t corrYield = (fUse2DMatching ? tempProj->GetBinContent(im,ipt) : tempProj->GetBinContent(im))/eff*evWeight;
          Double_t corrYieldErr = (fUse2DMatching ? tempProj->GetBinError(im,ipt) : tempProj->GetBinError(im))/eff*evWeight;
          if(fUse2DMatching) {
	    tempProj->SetBinContent(im,ipt, corrYield);
	    tempProj->SetBinError(im,ipt, corrYieldErr);
	  }
	  else {
	    tempProj->SetBinContent(im, corrYield);
	    tempProj->SetBinError(im, corrYieldErr);
	  }
        }
      }  // end loop over pt bins
            
      //cout << "GetJpsiSignal::centrality / eff / evWeight = " << fSEOS->GetAxis(kCentVZERO)->GetBinCenter(icent)
      //     << " / " << eff << " / " << evWeight << endl;
            /*
      if(!fHistSEOS) {
	if(fUse2DMatching)
	  fHistSEOS = (TH2D*)tempProj->Clone(Form("fHistSEOS_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
	            minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
	else
	  fHistSEOS = (TH1D*)tempProj->Clone(Form("fHistSEOS_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
		    minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
      }
      else*/ 
      fHistSEOS->Add(tempProj);
      delete tempProj;
    }  // end loop over centrality bins
  } // end if(fUsedVars[kCentrality] && (fEffVsPtCent || fWeightVsCent))
  else {
    if(fUse2DMatching) fHistSEOS = fSEOS->Projection(fVarIndex[kPt],fVarIndex[kMass]);
    else fHistSEOS = fSEOS->Projection(fVarIndex[kMass]);
    fHistSEOS->SetName(Form("fHistSEOS_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
			     minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
  }
  cout << "MakeSEOS  fHistSEOS(100)" << fHistSEOS->GetBinContent(100) << endl;
  //fHistSEOS->Draw();
}

//_____________________________________________________________________________________________
TH1* AliResonanceFits::GetLScorrected(TH1* selsLeg1, TH1* selsLeg2, TH1* meos, TH1* melsLeg1, TH1* melsLeg2) {
  //
  // Compute the R-factor corrected SE-LS 
  //
  TH1* mels = 0x0; TH1* sels = 0x0;
  if(fUse2DMatching) {
    mels = (TH2D*)meos->Clone(Form("mels%.6f", gRandom->Rndm())); mels->Reset();
    sels = (TH2D*)selsLeg1->Clone(Form("sels%.6f", gRandom->Rndm())); sels->Reset();
  }
  else {
    mels = (TH1D*)meos->Clone(Form("mels%.6f", gRandom->Rndm())); mels->Reset();
    sels = (TH1D*)selsLeg1->Clone(Form("sels%.6f", gRandom->Rndm())); sels->Reset();
  }
  for(Int_t ipt=1; ipt<=(fUse2DMatching ? meos->GetYaxis()->GetNbins() : 1); ++ipt) {
    for(Int_t im=1; im<=meos->GetXaxis()->GetNbins(); ++im) {
      Double_t me1 = (fUse2DMatching ? melsLeg1->GetBinContent(im,ipt) : melsLeg1->GetBinContent(im));
      Double_t me2 = (fUse2DMatching ? melsLeg2->GetBinContent(im,ipt) : melsLeg2->GetBinContent(im));
      Double_t me1Err = (fUse2DMatching ? melsLeg1->GetBinError(im,ipt) : melsLeg1->GetBinError(im));
      Double_t me2Err = (fUse2DMatching ? melsLeg2->GetBinError(im,ipt) : melsLeg2->GetBinError(im));
      // fLSmethod == 1  -> arithmetic mean
      Double_t melsVal = 0.5*(me1+me2);                                    // arithmetic mean
      Double_t melsValErr = 0.5*TMath::Sqrt(me1Err*me1Err+me2Err*me2Err);
      if(fLSmethod==2) {   // geometric mean
	melsVal = 1.0*TMath::Sqrt(me1*me2);
	if(me1>0. && me2>0.)
	  melsValErr = 1.0*TMath::Sqrt((me1/me2)*me2Err*me2Err + (me2/me1)*me1Err*me1Err);
	else melsValErr = 0.0;
      }
      if(fUse2DMatching) {
	mels->SetBinContent(im,ipt,melsVal); 
	mels->SetBinError(im,ipt,melsValErr);
      }
      else {
	mels->SetBinContent(im,melsVal); 
	mels->SetBinError(im,melsValErr);
      }
          
      Double_t se1 = (fUse2DMatching ? selsLeg1->GetBinContent(im,ipt) : selsLeg1->GetBinContent(im));
      Double_t se2 = (fUse2DMatching ? selsLeg2->GetBinContent(im,ipt) : selsLeg2->GetBinContent(im));
      Double_t se1Err = (fUse2DMatching ? selsLeg1->GetBinError(im,ipt) : selsLeg1->GetBinError(im));
      Double_t se2Err = (fUse2DMatching ? selsLeg2->GetBinError(im,ipt) : selsLeg2->GetBinError(im));
      // fLSmethod == 1 -> arithmetic mean
      Double_t selsVal = 1.0*(se1+se2);
      Double_t selsValErr = 1.0*TMath::Sqrt(se1Err*se1Err+se2Err*se2Err);
      if(fLSmethod==2) {   // geometric mean
	selsVal = 2.0*TMath::Sqrt(se1*se2);
	if(se1>0. && se2>0.)
	  selsValErr = 2.0*TMath::Sqrt((se1/se2)*se2Err*se2Err + (se2/se1)*se1Err*se1Err);
	else selsValErr = 0.0;
      }
      if(fUse2DMatching) {
	sels->SetBinContent(im,ipt,selsVal);
	sels->SetBinError(im,ipt,selsValErr);
      }
      else {
	sels->SetBinContent(im,selsVal);
	sels->SetBinError(im,selsValErr);
      }
    }  // end loop over mass bins
  }  // end loop over pt bins
  //meLS->Draw(); return;
  
  TH1* rFactor = 0x0;
  if(fUse2DMatching) 
    rFactor = (TH2D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm()));
  else
    rFactor = (TH1D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm()));
  rFactor->Divide(mels);
      
  TH1* seLScorr = 0x0;
  if(fUse2DMatching) 
    seLScorr = (TH2D*)sels->Clone(Form("seLScorr%.6f", gRandom->Rndm()));
  else
    seLScorr = (TH1D*)sels->Clone(Form("seLScorr%.6f", gRandom->Rndm()));
  seLScorr->Multiply(rFactor);
  delete rFactor;
  return seLScorr;
}

//_____________________________________________________________________________________________
void AliResonanceFits::MakeLSbkg() {
  //
  // Build the invariant mass spectrum for the LS bkg 
  //
  // Apply user requested ranges
   Int_t zIdx = fVarIndex[kVertexZ];
  Double_t minZ = fVertexZRange[0]; Double_t maxZ = fVertexZRange[1];
  if(fUsedVars[kVertexZ] && fVertexSelection) {
    //zIdx = fVarIndex[kVertexZ];
    minZ = fSEOS->GetAxis(zIdx)->GetXmin(); maxZ = fSEOS->GetAxis(zIdx)->GetXmax();
    if(fVertexZRange[0]>minZ) minZ = fVertexZRange[0]+0.001;
    if(fVertexZRange[1]<maxZ) maxZ = fVertexZRange[1]-0.001;
    fSELSleg1->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    fSELSleg2->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    fMEOS->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    fMELSleg1->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    fMELSleg2->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
  }
  Int_t epIdx = fVarIndex[kEP2];
  Double_t minEP = fEP2Range[0]; Double_t maxEP = fEP2Range[1];
  if(fUsedVars[kEP2] && fEPSelection) {
    //epIdx = fVarIndex[kEP2];
    minEP = fSEOS->GetAxis(epIdx)->GetXmin(); maxEP = fSEOS->GetAxis(epIdx)->GetXmax();
    if(fEP2Range[0]>minEP) minEP = fEP2Range[0]+0.001;
    if(fEP2Range[1]<maxEP) maxEP = fEP2Range[1]-0.001;
    fSELSleg1->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    fSELSleg2->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    fMEOS->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    fMELSleg1->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    fMELSleg2->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
  }
  Int_t ptIdx = fVarIndex[kPt];
  Double_t minPt = fPtRange[0]; Double_t maxPt = fPtRange[1];
  if(fUsedVars[kPt] && fPtSelection) {
    //ptIdx = fVarIndex[kPt];
    minPt = fSEOS->GetAxis(ptIdx)->GetXmin(); maxPt = fSEOS->GetAxis(ptIdx)->GetXmax();
    if(fPtRange[0]>minPt) minPt = fPtRange[0]+0.001;
    if(fPtRange[1]<maxPt) maxPt = fPtRange[1]-0.001;
    fSELSleg1->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    fSELSleg2->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    fMEOS->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    fMELSleg1->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    fMELSleg2->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
  }
  Int_t centIdx = fVarIndex[kCentrality];
  Double_t minCent = fCentralityRange[0]; Double_t maxCent = fCentralityRange[1];
  if(fUsedVars[kCentrality] && fCentralitySelection) {
    //centIdx = fVarIndex[kCentrality];
    minCent = fSEOS->GetAxis(centIdx)->GetXmin(); maxCent = fSEOS->GetAxis(centIdx)->GetXmax();
    if(fCentralityRange[0]>minCent) minCent = fCentralityRange[0]+0.001;
    if(fCentralityRange[1]<maxCent) maxCent = fCentralityRange[1]-0.001;
    fSELSleg1->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    fSELSleg2->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    fMEOS->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    fMELSleg1->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    fMELSleg2->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
  }
  Int_t mIdx = fVarIndex[kMass];
  Double_t minMass = fSEOS->GetAxis(mIdx)->GetXmin(); Double_t maxMass = fSEOS->GetAxis(mIdx)->GetXmax();
  if(fMassRange[0]>minMass) minMass = fMassRange[0]+0.00001;
  if(fMassRange[1]<maxMass) maxMass = fMassRange[1]-0.00001;
  fSELSleg1->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  fSELSleg2->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  fMEOS->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  fMELSleg1->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  fMELSleg2->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  
  TH1* selsLeg1_proj=0; TH1* selsLeg2_proj=0; TH1* melsLeg1_proj=0; TH1* melsLeg2_proj=0; TH1* meos_proj=0;
  
  cout << "AliResonanceFits::MakeLSbkg() fWeightVsCent (" << fWeightVsCent << ")  fEffVsPtCent(" << fEffVsPtCent << ")" << endl;
  
  // TODO: Implement properly the fUsedVars[] flags
  //       Right now, all dimensions are assumed to be present.
  for(Int_t ic=1; ic<=fMEOS->GetAxis(centIdx)->GetNbins(); ++ic) {         // ME bkg can have fewer bins in centrality compared to SE
    Double_t centrality = fMEOS->GetAxis(centIdx)->GetBinCenter(ic);
    Double_t centLowEdge = fMEOS->GetAxis(centIdx)->GetBinLowEdge(ic);
    Double_t centUpEdge = fMEOS->GetAxis(centIdx)->GetBinUpEdge(ic);
    if(centrality<minCent) continue;
    if(centrality>maxCent) continue;
    //cout << "centrality " << centrality << endl;
    
    Double_t evWeight = 1.0;
    if(fWeightVsCent) {
      // TODO: Make sure this is synchronized with the binning for the SE
      evWeight = fWeightVsCent->GetBinContent(ic);    //fWeightVsCent has the same binning as the THnF for the bkg
      //cout << "GetJpsiMEbkg::centrality / eff / evWeight = " << centrality << " / " << eff << " / " << evWeight << endl;
    }
        
    fMEOS->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
    fMELSleg1->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
    fMELSleg2->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
    fSELSleg1->GetAxis(centIdx)->SetRangeUser(centLowEdge+0.001, centUpEdge-0.001);
    fSELSleg2->GetAxis(centIdx)->SetRangeUser(centLowEdge+0.001, centUpEdge-0.001);
    
    for(Int_t ipt=1; ipt<=(fUse2DMatching ? fMEOS->GetAxis(ptIdx)->GetNbins() : 1); ++ipt) {
      Float_t pt = (fUse2DMatching ? fMEOS->GetAxis(ptIdx)->GetBinCenter(ipt) : 0.0);
      Int_t ptBin = ((fUse2DMatching && fEffVsPtCent) ? fEffVsPtCent->GetYaxis()->FindBin(pt) : 1);
      Int_t centBin = (fEffVsPtCent ? fEffVsPtCent->GetXaxis()->FindBin(fMEOS->GetAxis(fVarIndex[kCentrality])->GetBinCenter(ic)) : 0);
      
      if(pt<fPtRange[0] || pt>fPtRange[1]) continue;
      
      Double_t eff = 1.0;
      if(fEffVsPtCent) {
	if(centBin==0) centBin=1;
        if(centBin==fEffVsPtCent->GetXaxis()->GetNbins()+1) centBin = fEffVsPtCent->GetXaxis()->GetNbins();
	eff = fEffVsPtCent->GetBinContent(centBin, ptBin);
      }
      
      for(Int_t iz=1; iz<=fMEOS->GetAxis(zIdx)->GetNbins(); ++iz) {
        Double_t vtxz = fMEOS->GetAxis(zIdx)->GetBinCenter(iz);
        cout << "vtxz " << vtxz << endl;
        fMEOS->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        fMELSleg1->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        fMELSleg2->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        fSELSleg1->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        fSELSleg2->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
      
        for(Int_t iep=1; iep<=fMEOS->GetAxis(epIdx)->GetNbins(); ++iep) {
          Double_t ep = fMEOS->GetAxis(epIdx)->GetBinCenter(iep);
          cout << "EP " << ep << endl;
          fMEOS->GetAxis(epIdx)->SetRangeUser(ep, ep);
          fMELSleg1->GetAxis(epIdx)->SetRangeUser(ep, ep);
          fMELSleg2->GetAxis(epIdx)->SetRangeUser(ep, ep);
	  fSELSleg1->GetAxis(epIdx)->SetRangeUser(ep, ep);
	  fSELSleg2->GetAxis(epIdx)->SetRangeUser(ep, ep);
	
	  if(fUse2DMatching) {
	    selsLeg1_proj = fSELSleg1->Projection(ptIdx,mIdx); selsLeg1_proj->SetName(Form("selsLeg1_proj%.6f", gRandom->Rndm()));
	    selsLeg2_proj = fSELSleg2->Projection(ptIdx,mIdx); selsLeg2_proj->SetName(Form("selsLeg2_proj%.6f", gRandom->Rndm()));
	    melsLeg1_proj = fMELSleg1->Projection(ptIdx,mIdx); melsLeg1_proj->SetName(Form("melsLeg1_proj%.6f", gRandom->Rndm()));
	    melsLeg2_proj = fMELSleg2->Projection(ptIdx,mIdx); melsLeg2_proj->SetName(Form("melsLeg2_proj%.6f", gRandom->Rndm()));
	    meos_proj = fMEOS->Projection(ptIdx,mIdx); meos_proj->SetName(Form("meos_proj%.6f", gRandom->Rndm()));
	  }
	  else {
	    selsLeg1_proj = fSELSleg1->Projection(mIdx); selsLeg1_proj->SetName(Form("selsLeg1_proj%.6f", gRandom->Rndm()));
	    selsLeg2_proj = fSELSleg2->Projection(mIdx); selsLeg2_proj->SetName(Form("selsLeg2_proj%.6f", gRandom->Rndm()));
	    melsLeg1_proj = fMELSleg1->Projection(mIdx); melsLeg1_proj->SetName(Form("melsLeg1_proj%.6f", gRandom->Rndm()));
	    melsLeg2_proj = fMELSleg2->Projection(mIdx); melsLeg2_proj->SetName(Form("melsLeg2_proj%.6f", gRandom->Rndm()));
	    meos_proj = fMEOS->Projection(mIdx); meos_proj->SetName(Form("meos_proj%.6f", gRandom->Rndm()));
	  }
	  cout << "sels1/sels2/mels1/mels2/meos : " << selsLeg1_proj->GetEntries() << "/"
                  << selsLeg2_proj->GetEntries() << "/" << melsLeg1_proj->GetEntries() << "/" << melsLeg2_proj->GetEntries() << "/"
                  << meos_proj->GetEntries() << endl;
	  	
	  TH1* seLScorr = GetLScorrected(selsLeg1_proj, selsLeg2_proj, meos_proj, melsLeg1_proj, melsLeg2_proj);
	  cout << "seLScorr/n/low/high/evWeight/eff : " << seLScorr->GetBinContent(2) << "/" << seLScorr->GetXaxis()->GetNbins() 
                  << "/" << seLScorr->GetXaxis()->GetXmin() << "/" << seLScorr->GetXaxis()->GetXmax() << "/" << evWeight << "/" << eff << endl;
          
          
	  /*if(!fHistSELSbkg) {
	    if(fUse2DMatching) {
              fHistSELSbkg = fMEOS->Projection(ptIdx,mIdx);
              fHistSELSbkg->SetName(Form("fHistSELSbkg_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
                                                                minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
               //fHistSELSbkg = (TH2D*)selsLeg1_proj->Clone(Form("fHistSELSbkg_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
               //                                           minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
            }
	    else {
               fHistSELSbkg = fMEOS->Projection(mIdx);
               fHistSELSbkg->SetName(Form("fHistSELSbkg_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
                                          minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
               //fHistSELSbkg = (TH1D*)selsLeg1_proj->Clone(Form("fHistSELSbkg_pt%.2f_%.2f_cent%.1f_%.1f_msig%.2f_%.2f_mbkg%.2f_%.2f", 
                //                                          minPt, maxPt, minCent, maxCent, fSignalRange[0], fSignalRange[1], fFitRange[0], fFitRange[1]));
            }
	    //fHistSELSbkg->Scale(evWeight/eff);
            //fHistSELSbkg->Reset();
            //fHistSELSbkg->Add(seLScorr, evWeight/eff);
	  } */
	  fHistSELSbkg->Add(seLScorr, evWeight/eff);
          //fHistSELSbkg->Draw(); return;
	
	  delete selsLeg1_proj; delete selsLeg2_proj; delete melsLeg1_proj; 
	  delete melsLeg2_proj; delete meos_proj; delete seLScorr;
        } // end loop over EP bins
      }  // end loop over z-vertex bins
    }  // end loop over pt bins
  }  // end loop over centrality bins
  //fHistSELSbkg->Draw();
}

//_____________________________________________________________________________________________
void AliResonanceFits::MakeMEbkg() {
  //
  // Build the invariant mass spectrum for the ME bkg
  // The ME bkg is scaled either to the region in SEOS around the peak, excluding the peak region (fMEMatchOption=1)
  //  or it is scaled to the SELS bkg corrected with the R factor
  //
  // Apply user requested ranges
   Int_t zIdx = fVarIndex[kVertexZ];
  Double_t minZ = fVertexZRange[0]; Double_t maxZ = fVertexZRange[1];
  if(fUsedVars[kVertexZ] && fVertexSelection) {
    //zIdx = fVarIndex[kVertexZ];
    minZ = fSEOS->GetAxis(zIdx)->GetXmin(); maxZ = fSEOS->GetAxis(zIdx)->GetXmax();
    if(fVertexZRange[0]>minZ) minZ = fVertexZRange[0]+0.001;
    if(fVertexZRange[1]<maxZ) maxZ = fVertexZRange[1]-0.001;
    fMEOS->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    if(fMEMatchOption==1) fSEOS->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    else {
      fSELSleg1->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
      fSELSleg2->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
      fMELSleg1->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
      fMELSleg2->GetAxis(zIdx)->SetRangeUser(minZ, maxZ);
    }
  }
  Int_t epIdx = fVarIndex[kEP2];
  Double_t minEP = fEP2Range[0]; Double_t maxEP = fEP2Range[1];
  if(fUsedVars[kEP2] && fEPSelection) {
    //epIdx = fVarIndex[kEP2];
    minEP = fSEOS->GetAxis(epIdx)->GetXmin(); maxEP = fSEOS->GetAxis(epIdx)->GetXmax();
    if(fEP2Range[0]>minEP) minEP = fEP2Range[0]+0.001;
    if(fEP2Range[1]<maxEP) maxEP = fEP2Range[1]-0.001;
    fMEOS->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    if(fMEMatchOption==1) fSEOS->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    else {
      fSELSleg1->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
      fSELSleg2->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
      fMELSleg1->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
      fMELSleg2->GetAxis(epIdx)->SetRangeUser(minEP, maxEP);
    }
  }
  Int_t ptIdx = fVarIndex[kPt];
  Double_t minPt = fPtRange[0]; Double_t maxPt = fPtRange[1];
  if(fUsedVars[kPt] && fPtSelection) {
    ptIdx = fVarIndex[kPt];
    minPt = fSEOS->GetAxis(ptIdx)->GetXmin(); maxPt = fSEOS->GetAxis(ptIdx)->GetXmax();
    if(fPtRange[0]>minPt) minPt = fPtRange[0]+0.001;
    if(fPtRange[1]<maxPt) maxPt = fPtRange[1]-0.001;
    fMEOS->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    if(fMEMatchOption==1) fSEOS->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    else {
      fSELSleg1->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
      fSELSleg2->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
      fMELSleg1->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
      fMELSleg2->GetAxis(ptIdx)->SetRangeUser(minPt, maxPt);
    }
  }
  Int_t centIdx = fVarIndex[kCentrality];
  Double_t minCent = fCentralityRange[0]; Double_t maxCent = fCentralityRange[1];
  if(fUsedVars[kCentrality] && fCentralitySelection) {
    centIdx = fVarIndex[kCentrality];
    minCent = fSEOS->GetAxis(centIdx)->GetXmin(); maxCent = fSEOS->GetAxis(centIdx)->GetXmax();
    if(fCentralityRange[0]>minCent) minCent = fCentralityRange[0]+0.001;
    if(fCentralityRange[1]<maxCent) maxCent = fCentralityRange[1]-0.001;
    fMEOS->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    if(fMEMatchOption==1) fSEOS->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    else {
      fSELSleg1->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
      fSELSleg2->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
      fMELSleg1->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
      fMELSleg2->GetAxis(centIdx)->SetRangeUser(minCent, maxCent);
    }
  }
  Int_t mIdx = fVarIndex[kMass];
  Double_t minMass = fSEOS->GetAxis(mIdx)->GetXmin(); Double_t maxMass = fSEOS->GetAxis(mIdx)->GetXmax();
  if(fMassRange[0]>minMass) minMass = fMassRange[0]+0.00001;
  if(fMassRange[1]<maxMass) maxMass = fMassRange[1]-0.00001;
  fMEOS->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  if(fMEMatchOption==1) fSEOS->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  else {
    fSELSleg1->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
    fSELSleg2->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
    fMELSleg1->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
    fMELSleg2->GetAxis(mIdx)->SetRangeUser(minMass, maxMass);
  }
  
  TH1* selsLeg1_proj=0; TH1* selsLeg2_proj=0; TH1* melsLeg1_proj=0; TH1* melsLeg2_proj=0; 
  TH1* seos_proj=0; TH1* meos_proj=0;
  
  for(Int_t ic=1; ic<=fMEOS->GetAxis(centIdx)->GetNbins(); ++ic) {         // ME bkg has fewer bins in centrality
    Double_t centrality = fMEOS->GetAxis(centIdx)->GetBinCenter(ic);
    Double_t centLowEdge = fMEOS->GetAxis(centIdx)->GetBinLowEdge(ic);
    Double_t centUpEdge = fMEOS->GetAxis(centIdx)->GetBinUpEdge(ic);
    if(centrality<minCent) continue;
    if(centrality>maxCent) continue;
    cout << "centrality " << centrality << endl;
						
    Double_t evWeight = 1.0;
    if(fWeightVsCent) {
      evWeight = fWeightVsCent->GetBinContent(ic);    //hWeights has the same binning as the THnF for the bkg
      //cout << "GetJpsiMEbkg::centrality / eff / evWeight = " << centrality << " / " << eff << " / " << evWeight << endl;
    }
        
    fMEOS->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
    if(fMEMatchOption==1) fSEOS->GetAxis(centIdx)->SetRangeUser(centLowEdge+0.001, centUpEdge-0.001);
    else {
      fMELSleg1->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
      fMELSleg2->GetAxis(centIdx)->SetRangeUser(centrality, centrality);
      fSELSleg1->GetAxis(centIdx)->SetRangeUser(centLowEdge+0.001, centUpEdge-0.001);
      fSELSleg2->GetAxis(centIdx)->SetRangeUser(centLowEdge+0.001, centUpEdge-0.001);
    }
    
    for(Int_t ipt=1; ipt<=(fUse2DMatching ? fMEOS->GetAxis(ptIdx)->GetNbins() : 1); ++ipt) {
      Float_t pt = (fUse2DMatching ? fMEOS->GetAxis(ptIdx)->GetBinCenter(ipt) : 0.0);
      Int_t ptBin = ((fUse2DMatching && fEffVsPtCent) ? fEffVsPtCent->GetYaxis()->FindBin(pt) : 1);
      Int_t centBin = (fEffVsPtCent ? fEffVsPtCent->GetXaxis()->FindBin(fMEOS->GetAxis(fVarIndex[kCentrality])->GetBinCenter(ic)) : 0);
      
      if(pt<fPtRange[0] || pt>fPtRange[1]) continue;
      
      cout << "pt/ipt/npt :: " << pt << "/" << ipt << "/" << fMEOS->GetAxis(ptIdx)->GetNbins() << endl;
      
      // If the efficiency vs (pt,centrality) is provided, the use it to weight the ME 
      Double_t eff = 1.0;
      if(fEffVsPtCent) {
	if(centBin==0) centBin=1;
        if(centBin==fEffVsPtCent->GetXaxis()->GetNbins()+1) centBin = fEffVsPtCent->GetXaxis()->GetNbins();
	eff = fEffVsPtCent->GetBinContent(centBin, ptBin);
      }
    
      for(Int_t iz=1; iz<=fMEOS->GetAxis(zIdx)->GetNbins(); ++iz) {
        Double_t vtxz = fMEOS->GetAxis(zIdx)->GetBinCenter(iz);
        cout << "vtxz " << vtxz << endl;
        fMEOS->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        if(fMEMatchOption==1) fSEOS->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        else {
          fMELSleg1->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
          fMELSleg2->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
          fSELSleg1->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
          fSELSleg2->GetAxis(zIdx)->SetRangeUser(vtxz, vtxz);
        }
      
        for(Int_t iep=1; iep<=fMEOS->GetAxis(epIdx)->GetNbins(); ++iep) {
          Double_t ep = fMEOS->GetAxis(epIdx)->GetBinCenter(iep);
          cout << "EP " << ep << endl;
          fMEOS->GetAxis(epIdx)->SetRangeUser(ep, ep);
          if(fMEMatchOption==1) fSEOS->GetAxis(epIdx)->SetRangeUser(ep, ep);
	  else {
	    fMELSleg1->GetAxis(epIdx)->SetRangeUser(ep, ep);
            fMELSleg2->GetAxis(epIdx)->SetRangeUser(ep, ep);
	    fSELSleg1->GetAxis(epIdx)->SetRangeUser(ep, ep);
	    fSELSleg2->GetAxis(epIdx)->SetRangeUser(ep, ep);
	  }
	  
	  cout << "1" << endl;
	  if(fUse2DMatching) meos_proj = fMEOS->Projection(ptIdx,mIdx);
	  else meos_proj = fMEOS->Projection(mIdx);
	  cout << "2" << endl;
	  meos_proj->SetName(Form("meos_proj%.6f", gRandom->Rndm()));
	  	
	  if(fMEMatchOption==1) {
	    if(fUse2DMatching) seos_proj = fSEOS->Projection(ptIdx,mIdx);
	    else seos_proj = fSEOS->Projection(mIdx);
            seos_proj->SetName(Form("seos_proj%.6f", gRandom->Rndm()));
	    	    
	    fFitValues[kBkgScale]=1.0;
	    if(fMatchingOption==1) ComputeWeightedScale(seos_proj, meos_proj);
	    if(fMatchingOption==2) FitScale(seos_proj, meos_proj);
	    if(fMatchingOption==3) ComputeEntryScale(seos_proj, meos_proj);
	    cout << "fFitValues[kBkgScale] = " << fFitValues[kBkgScale] << endl;
	    meos_proj->Scale(fFitValues[kBkgScale]);
	    delete seos_proj;
	  }
	  if(fMEMatchOption==2) {
	    if(fUse2DMatching) {
	      selsLeg1_proj = fSELSleg1->Projection(ptIdx,mIdx); selsLeg1_proj->SetName(Form("selsLeg1_proj%.6f", gRandom->Rndm()));
	      selsLeg2_proj = fSELSleg2->Projection(ptIdx,mIdx); selsLeg2_proj->SetName(Form("selsLeg2_proj%.6f", gRandom->Rndm()));
	      melsLeg1_proj = fMELSleg1->Projection(ptIdx,mIdx); melsLeg1_proj->SetName(Form("melsLeg1_proj%.6f", gRandom->Rndm()));
	      melsLeg2_proj = fMELSleg2->Projection(ptIdx,mIdx); melsLeg2_proj->SetName(Form("melsLeg2_proj%.6f", gRandom->Rndm()));
	    }
	    else {
	      selsLeg1_proj = fSELSleg1->Projection(mIdx); selsLeg1_proj->SetName(Form("selsLeg1_proj%.6f", gRandom->Rndm()));
	      selsLeg2_proj = fSELSleg2->Projection(mIdx); selsLeg2_proj->SetName(Form("selsLeg2_proj%.6f", gRandom->Rndm()));
	      melsLeg1_proj = fMELSleg1->Projection(mIdx); melsLeg1_proj->SetName(Form("melsLeg1_proj%.6f", gRandom->Rndm()));
	      melsLeg2_proj = fMELSleg2->Projection(mIdx); melsLeg2_proj->SetName(Form("melsLeg2_proj%.6f", gRandom->Rndm()));
	    }	    
	  
	    TH1* seLScorr = GetLScorrected(selsLeg1_proj, selsLeg2_proj, meos_proj, melsLeg1_proj, melsLeg2_proj);
	    fFitValues[kBkgScale] = 1.0;
	    
	    // NOTE: The matching with the SELS corrected bkg should be done over the full mass range
	    Double_t exclRange_backup[2]={fExclusionRange[0],fExclusionRange[1]};
	    fExclusionRange[0]=0.0; fExclusionRange[1]=0.0;
	    if(fMatchingOption==1) ComputeWeightedScale(seLScorr, meos_proj);
	    if(fMatchingOption==2) FitScale(seLScorr, meos_proj);
	    if(fMatchingOption==3) ComputeEntryScale(seLScorr, meos_proj); 
	    fExclusionRange[0] = exclRange_backup[0];
	    fExclusionRange[1] = exclRange_backup[1];
	    //cout << "fFitValues[kBkgScale] = " << fFitValues[kBkgScale] << endl;
	    meos_proj->Scale(fFitValues[kBkgScale]);
	    delete selsLeg1_proj; delete selsLeg2_proj; 
	    delete melsLeg1_proj; delete melsLeg2_proj; 
	    delete seLScorr;
	  }
	
	  /*if(!fHistMEbkg) {
	    if(fUse2DMatching) 
	      fHistMEbkg = (TH2D*)meos_proj->Clone(Form("fHistMEbkg%.2f_%.2f_cent%.1f_%.1f", 
		  				        minPt, maxPt, minCent, maxCent));
	    else
	      fHistMEbkg = (TH1D*)meos_proj->Clone(Form("fHistMEbkg%.2f_%.2f_cent%.1f_%.1f", 
		  				        minPt, maxPt, minCent, maxCent));
	    fHistMEbkg->Scale(evWeight/eff);
	  }
	  else*/
          fHistMEbkg->Add(meos_proj, evWeight/eff);
	
	  delete meos_proj;
        } // end loop over EP bins
      }  // end loop over z-vertex bins
    }  // end loop over pt
  }  // end loop over centrality bins
}

//_____________________________________________________________________________________________
void AliResonanceFits::DrawSignalExtraction(Bool_t save /*=kFALSE*/, const Char_t* name /*=""*/, const Char_t* outputDir /*=""*/,
                                            Bool_t makeNicePlot /*=kFALSE*/, 
                                            TVirtualPad* externalPad /*=0x0*/, Bool_t noYlabels /*=kFALSE*/, Bool_t noLegends /*=kFALSE*/) {
  //
  // Draw the results of signal extraction
  //
  cout << "1" << endl;
  if(!fIsProcessed) Process();

  if(fUse2DMatching) {
    if(fPlottingOption==2) {
      TCanvas *c1=new TCanvas();
      fHistSEOS->DrawClone("colz");
      TCanvas *c2=new TCanvas();
      if(fBkgMethod==1)
        fHistMEbkg->DrawClone("colz");
      else
        fHistSELSbkg->DrawClone("colz");
      TCanvas *c3=new TCanvas();
      fHistBkgSubtracted->DrawClone("colz");
      //TCanvas *c4=new TCanvas();
    }
    if(fPlottingOption==0) {
      TH1D* projSEOS = ((TH2D*)fHistSEOS)->ProjectionX(Form("%s_massProj",fHistSEOS->GetName()),0);
      TH1* bkg = fHistMEbkg; 
      if(fBkgMethod!=1) bkg = fHistSELSbkg;
      TH1D* projBkg = ((TH2D*)bkg)->ProjectionX(Form("%s_massProj",bkg->GetName()),0);
      TH1D* projSignal = ((TH2D*)fHistBkgSubtracted)->ProjectionX(Form("%s_massProj",fHistBkgSubtracted->GetName()),0);
      DrawMassProjection(projSEOS, projBkg, projSignal, save, name, outputDir, makeNicePlot, externalPad, noYlabels, noLegends);
    }
    if(fPlottingOption==1) {
      TH1D* projSEOS = ((TH2D*)fHistSEOS)->ProjectionY(Form("%s_ptProj",fHistSEOS->GetName()),
						       fHistSEOS->GetXaxis()->FindBin(fSignalRange[0]+0.001),
						       fHistSEOS->GetXaxis()->FindBin(fSignalRange[1]-0.001));
      TH1* bkg = fHistMEbkg; 
      if(fBkgMethod!=1) bkg = fHistSELSbkg;
      TH1D* projBkg = ((TH2D*)bkg)->ProjectionY(Form("%s_ptProj",bkg->GetName()),
						fHistSEOS->GetXaxis()->FindBin(fSignalRange[0]+0.001),
						fHistSEOS->GetXaxis()->FindBin(fSignalRange[1]-0.001));
      TH1D* projSignal = ((TH2D*)fHistBkgSubtracted)->ProjectionY(Form("%s_ptProj",fHistBkgSubtracted->GetName()),
								  fHistSEOS->GetXaxis()->FindBin(fSignalRange[0]+0.001),
								  fHistSEOS->GetXaxis()->FindBin(fSignalRange[1]-0.001));
      TCanvas *c1=new TCanvas();
      projSEOS->DrawClone(); projBkg->DrawClone("same");
      TCanvas *c2=new TCanvas();
      projSignal->DrawClone();
      //DrawMassProjection(projSEOS, projBkg, projSignal, save, name, outputDir, makeNicePlot, externalPad, noYlabels, noLegends);
    }
  }   // end if(fUse2DMatching)
  else {
    TH1* bkg=0;
    if(fBkgMethod==1) bkg = fHistMEbkg;
    else bkg = fHistSELSbkg;
    DrawMassProjection(fHistSEOS, bkg, fHistBkgSubtracted, save, name, outputDir, makeNicePlot, externalPad, noYlabels, noLegends);
  }
}


//_____________________________________________________________________________________________
void AliResonanceFits::DrawMassProjection(TH1* seos, TH1* bkg, TH1* signal,
                                          Bool_t save /*=kFALSE*/, const Char_t* name /*=""*/, const Char_t* outputDir /*=""*/,
                                          Bool_t makeNicePlot /*=kFALSE*/, 
                                          TVirtualPad* externalPad /*=0x0*/, Bool_t noYlabels /*=kFALSE*/, Bool_t noLegends /*=kFALSE*/) {
  //
  // Draw the mass projection for the signal extraction procedure
  //
  //bkg->Draw(); return;
   
  fFitValues[kChisqMC]=0.0;
  Double_t ndfMC=0.0;
  if(fSignalMCshape) {
    fSignalMCshape->SetLineWidth(2.0);
    fSignalMCshape->SetLineColor(1);
    
    Float_t exclRange_backup[2] = {fExclusionRange[0], fExclusionRange[1]};
    fExclusionRange[0]=0.0; fExclusionRange[1]=0.0;
    fFitValues[kBkgScale]=1.0;
    if(fMatchingOption==1) ComputeWeightedScale(signal, fSignalMCshape);
    if(fMatchingOption==2) FitScale(signal, fSignalMCshape);
    if(fMatchingOption==3) {
      Int_t bin1 = signal->GetXaxis()->FindBin(fFitRange[0]+1.0e-6);
      Int_t bin2 = signal->GetXaxis()->FindBin(fFitRange[1]-1.0e-6);
      fFitValues[kBkgScale] = (fSignalMCshape->Integral(bin1,bin2)>1.0e-6 ? signal->Integral(bin1,bin2) / fSignalMCshape->Integral(bin1,bin2) : 0.0);
    }
    fSignalMCshape->Scale(fFitValues[kBkgScale]);
    fExclusionRange[0]=exclRange_backup[0]; fExclusionRange[1]=exclRange_backup[1];	
    
    for(Int_t ibin=1; ibin<=signal->GetXaxis()->GetNbins(); ++ibin) {
      Double_t mass = signal->GetXaxis()->GetBinCenter(ibin);
      Double_t s    = signal->GetBinContent(ibin);
      Double_t sErr = signal->GetBinError(ibin);
      Double_t mc   = fSignalMCshape->GetBinContent(fSignalMCshape->FindBin(mass));
      if(mass<fMassRange[0] || mass>fMassRange[1]) continue;
      if(sErr>0.001) {
        fFitValues[kChisqMC] += ((s-mc)*(s-mc)/sErr/sErr);
        ndfMC+=1.0;
      }
    }
    if(ndfMC>0) fFitValues[kChisqMC] /= Double_t(ndfMC);
  }  // end if (fSignalMCshape)
  cout << "2" << endl;
 
  TLatex* lat=new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.05);
  lat->SetTextColor(2);
  TCanvas* cv=0x0;
  TVirtualPad* c1=0x0;
  if(externalPad) c1=externalPad;
  else {
    TCanvas* oldCanvas = (TCanvas*)gROOT->FindObject("AliResonanceFitsCanvas");
    if(oldCanvas) delete oldCanvas; 
    cv=new TCanvas("AliResonanceFitsCanvas", "AliResonanceFits canvas", 720, 900);
    c1=cv->cd();
  }
  c1->SetTopMargin(0.01);
  c1->SetRightMargin(0.005);
  c1->SetBottomMargin(0.01);
  c1->SetLeftMargin(0.005);
  c1->Divide(1,2,0.0,0.0);
  TVirtualPad* pad = c1->cd(1);
  pad->SetTopMargin(0.01);
  pad->SetRightMargin(0.02);
  pad->SetBottomMargin(0.0);
  pad->SetLeftMargin(0.15); if(noYlabels) pad->SetLeftMargin(0.068);
  bkg->SetLineColor(2);
  seos->SetStats(kFALSE);
  if(fCentralityRange[0]>0.0)
    seos->SetTitle(Form("Centrality: %d - %d", TMath::Nint(fCentralityRange[0]), TMath::Nint(fCentralityRange[1])));
  seos->GetXaxis()->SetRangeUser(fMassRange[0], fMassRange[1]);
  seos->GetYaxis()->SetRangeUser(0.01, seos->GetMaximum()*1.20);
  seos->SetLineWidth(3.0);
  seos->GetYaxis()->SetLabelSize(0.05);
  seos->GetYaxis()->SetTitle(Form("entries per %.0f MeV/#it{c}^{2}", seos->GetXaxis()->GetBinWidth(1)*1000.0));
  if(noYlabels) seos->GetYaxis()->SetTitle("");
  seos->GetYaxis()->SetTitleSize(0.08);
  seos->GetYaxis()->SetTitleOffset(0.84);
  seos->GetYaxis()->CenterTitle();
  seos->SetTitle("");
  seos->GetXaxis()->SetTitleFont(42);
  seos->GetYaxis()->SetTitleFont(42);
  seos->GetXaxis()->SetLabelFont(42);
  seos->GetYaxis()->SetLabelFont(42);
  seos->DrawClone("E");
  bkg->SetStats(kFALSE);
  bkg->SetTitle("");
  bkg->SetLineStyle(2.0);
  bkg->SetLineWidth(2.0);
  bkg->DrawClone("sameE");
  
  cout << "3" << endl;
  
  TLegend* legend1 = new TLegend((noYlabels ? 0.08 : 0.16), 0.05, (noYlabels ? 0.32 : 0.40), 0.23);
  legend1->SetTextSize(0.06);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  legend1->AddEntry(seos, "SE opposite-sign", "l");
  legend1->AddEntry(bkg, (fBkgMethod==1 ? "ME bkg." : "SE-LS bkg."), "l");
  if(!noLegends)
    legend1->Draw();
  TLine line;
  line.SetLineColor(2);
  line.SetLineStyle(2);
  line.SetLineWidth(2.0);
  line.DrawLine(fSignalRange[0], 0.0, fSignalRange[0], seos->GetMaximum()*1.2);
  line.DrawLine(fSignalRange[1], 0.0, fSignalRange[1], seos->GetMaximum()*1.2);
  lat->SetTextFont(42);
  lat->SetTextColor(1);
  lat->SetTextSize(0.060);
  lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.92, Form("Centrality: %.0f - %.0f %%", fCentralityRange[0], fCentralityRange[1]));
  lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.82, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fMassRange[0], fMassRange[1], fFitValues[kChisq]));
  if(!makeNicePlot) {
    if(fEventVsCent)
      lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.72, Form("# events = %.2e", fEventVsCent->Integral(fEventVsCent->GetXaxis()->FindBin(fCentralityRange[0]+0.001), fEventVsCent->GetXaxis()->FindBin(fCentralityRange[1]-0.001))));
    lat->DrawLatex(0.56, 0.54, Form("Bkg. scaling range: %.2f-%.2f GeV/c^{2}", fFitRange[0], fFitRange[1]));
    lat->DrawLatex(0.56, 0.48, Form("Sig. range: %.2f-%.2f GeV/c^{2}", fSignalRange[0], fSignalRange[1]));
  }
    
    cout << "4" << endl;
    
  pad = c1->cd(2);
  pad->SetTopMargin(0.0);
  pad->SetRightMargin(0.02);
  pad->SetBottomMargin(0.17);
  pad->SetLeftMargin(0.15); if(noYlabels) pad->SetLeftMargin(0.068);
  signal->SetStats(kFALSE);
  signal->SetTitle("");
  signal->SetMarkerStyle(20); signal->SetMarkerColor(2);
  signal->SetLineWidth(1.0); signal->SetLineColor(2);
  signal->GetXaxis()->SetRangeUser(fMassRange[0], fMassRange[1]);
  signal->GetYaxis()->SetRangeUser(signal->GetMinimum()*1.4, signal->GetMaximum()*1.6);
  signal->SetLineWidth(2.0);
  signal->GetYaxis()->SetLabelSize(0.05);
  signal->GetYaxis()->SetTitleSize(0.08);
  //signal->GetYaxis()->SetTitle(Form("entries per %.0f MeV/#it{c}^{2}",signal->GetXaxis()->GetBinWidth(1)*1000.0));
  if(noYlabels) signal->GetYaxis()->SetTitle("");
  signal->GetYaxis()->SetTitleOffset(1.2);
  signal->GetXaxis()->SetLabelSize(0.07);
  signal->GetXaxis()->SetTitleSize(0.08);
  signal->GetXaxis()->SetTitleOffset(0.84);
  signal->GetXaxis()->SetTitle("#it{m} (GeV/#it{c}^{2})");
  signal->GetXaxis()->SetTitleFont(42);
  signal->GetYaxis()->SetTitleFont(42);
  signal->GetXaxis()->SetLabelFont(42);
  signal->GetYaxis()->SetLabelFont(42);
  signal->GetXaxis()->SetNdivisions(506);
  signal->GetYaxis()->SetNdivisions(506);
  
  cout << "5" << endl;
  
  signal->DrawClone();
  if(fSignalMCshape) fSignalMCshape->DrawClone("sameHISTL");
  line.DrawLine(fMassRange[0],0.0,fMassRange[1],0.0);
  line.DrawLine(fSignalRange[0], signal->GetMinimum()*1.05, fSignalRange[0], signal->GetMaximum()*1.5);
  line.DrawLine(fSignalRange[1], signal->GetMinimum()*1.05, fSignalRange[1], signal->GetMaximum()*1.5);
  lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.92, Form("Signal: %.0f #pm %.0f", fFitValues[kSig], fFitValues[kSigErr]));
  lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.83, Form("S/B: %.2f #pm %.2f %%", 100.0*fFitValues[kSoverB], 100.0*fFitValues[kSoverBerr]));
  lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.74, Form("S/#sqrt{S+B}: %.2f #pm %.2f", fFitValues[kSignif], fFitValues[kSignifErr]));
  if(fSignalMCshape)
    lat->DrawLatex((noYlabels ? 0.52 : 0.56), 0.65, Form("#chi^{2}/NDF (%.1f, %.1f) = %.2f", fMassRange[0], fMassRange[1], fFitValues[kChisqMC]));
  TLegend* legend2 = new TLegend((noYlabels ? 0.08 : 0.16), 0.78, (noYlabels ? 0.37 : 0.43), 0.98);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  legend2->AddEntry(signal, "Data", "lp");
  if(fSignalMCshape) legend2->AddEntry(fSignalMCshape, "MC", "l");
  if(!noLegends)
    legend2->Draw();
  
  cout << "6" << endl;
  
  if(save) {
    c1->SaveAs(Form("%s/%s.png", outputDir, name));
    delete c1;
  }
}
