/************************************************************************* 
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

// Author: Redmer Alexander Bertens, Utrecht University, Utrecht, Netherlands
//         (rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl)
//
// Tools class for Jet Flow Analysis, replaces 'extractJetFlow.C' macro
// Class is designed to manipulate the raw output of the jet flow
// tasks in PWG and PWGJE 
// to use this class, see $ALICE_ROOT/PWGCF/FLOW/macros/jetFlowTools.C

// root includes
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TList.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMath.h"
// aliroot include
#include "AliUnfolding.h"
#include "AliAnaChargedJetResponseMaker.h"
// local include
#include "AliJetFlowTools.h"

using namespace std;

//_____________________________________________________________________________
AliJetFlowTools::AliJetFlowTools() :
    fActiveString       (""),
    fActiveDir          (0x0),
    fInputList          (0x0),
    fOutputFileName     ("UnfoldedSpectra.root"),
    fOutputFile         (0x0),
    fCentralityBin      (0),
    fdPhidPt            (kTRUE),
    fDetectorResponse   (0x0),
    fBeta               (.001),
    fBetaPerDOF         (.0005),
    fAvoidRoundingError (kFALSE),
    fBinsTrue           (0x0),
    fBinsRec            (0x0),
    fDeltaPtDeltaPhi    (0x0),
    fJetPtDeltaPhi      (0x0),
    fSpectrumIn         (0x0),
    fSpectrumOut        (0x0),
    fDptInDist          (0x0),
    fDptOutDist         (0x0),
    fDptIn              (0x0),
    fDptOut             (0x0),
    fFullResponseIn     (0x0),
    fFullResponseOut    (0x0),
    fUnfoldedIn         (0x0),
    fUnfoldedOut        (0x0) { /* constructor */ }
//_____________________________________________________________________________
void AliJetFlowTools::Make() {
    // core function of the class
    // 1) manipulates input  histograms to prepare them for the unfolding procedure
    // 2) calls the unfolding routine
    // 3) writes output to file
    
    // 1) manipulation of input histograms
    // check if the input variables are present
    if(!PrepareForUnfolding()) {
        printf(" AliJetFlowTools::Make() Fatal error \n - couldn't prepare for unfolding ! \n");
        return;
    }
    // the unfolding algorithm needs a starting point.
    // as a starting point, take the input spectrum, rebinned to the
    // dimensions of the desired output spectrum
    TH1D* unfoldingTemplateIn  = GetUnfoldingTemplate(fSpectrumIn, fBinsTrue, TString("in"));   
    TH1D* unfoldingTemplateOut = GetUnfoldingTemplate(fSpectrumOut, fBinsTrue, TString("out"));

     // scale to the number of events
//    NormalizeTH1D(unfoldingTemplateIn, 1000);
//    NormalizeTH1D(unfoldingTemplateOut, 1000);
    
    // resize the jet spectrum
    TH1D* resizedJetPtIn  = ResizeXaxisTH1D(fSpectrumIn, 30, 110, TString("in"));
    TH1D* resizedJetPtOut = ResizeXaxisTH1D(fSpectrumOut, 30, 110, TString("out"));
//    NormalizeTH1D(resizedJetPt, 1000);

    // normalize the detector response matrix (dpt is already normalized by design)
    fDetectorResponse = NormalizeTH2D(fDetectorResponse);

    // get the full response matrix
    fFullResponseIn  = MatrixMultiplicationTH2D(fDptIn, fDetectorResponse);
    fFullResponseOut = MatrixMultiplicationTH2D(fDptOut, fDetectorResponse);

    // normalize the product
    NormalizeTH2D(fFullResponseIn);
    NormalizeTH2D(fFullResponseOut);

    // rebin the full response
    // as input the desired binning scheme is necessary, 
    // supplied as a TArrayD which hold the bin borders
    TH2D* fullResponseRebinnedIn  = RebinTH2DX(fFullResponseIn, fBinsTrue, TString("in"));
    TH2D* fullResponseRebinnedOut = RebinTH2DX(fFullResponseOut, fBinsTrue, TString("out"));

    // perform and check the normalization
    NormalizeTH2D(fFullResponseIn);
    NormalizeTH2D(fFullResponseOut);
    NormalizeTH2D(fullResponseRebinnedIn);
    NormalizeTH2D(fullResponseRebinnedOut);

    // resize the histo
    TH2D* resizedResonseIn  = ResizeYaxisTH2D(fullResponseRebinnedIn, fBinsTrue, fBinsRec, TString("in"));
    TH2D* resizedResonseOut = ResizeYaxisTH2D(fullResponseRebinnedOut, fBinsTrue, fBinsRec, TString("out"));
    TH1D* kinematicEfficiencyIn  = resizedResonseIn->ProjectionX();
    TH1D* kinematicEfficiencyOut = resizedResonseOut->ProjectionX();

    kinematicEfficiencyIn->SetNameTitle("kin_eff_IN","kin_eff_IN");
    kinematicEfficiencyOut->SetNameTitle("kin_eff_OUT", "kin_eff_OUT");

    // 2, 3) call the actual unfolding. results and transient objects are stored in a dedicated TDirectoryFile
    fActiveDir->cd();                   // select active dir
    TDirectoryFile* dirIn = new TDirectoryFile(Form("InPlane___%s", fActiveString.Data()), Form("InPlane___%s", fActiveString.Data()));
    dirIn->cd();                        // select inplane subdir
    Bool_t convergedIn = UnfoldSpectrum(       // do the inplane unfolding
            resizedJetPtIn,
            resizedResonseIn,
            kinematicEfficiencyIn,
            unfoldingTemplateIn,
            fUnfoldedIn, 
            TString("in"));
    fSpectrumIn->SetNameTitle("JetSpectrum", "Jet spectrum, in plane");
    fSpectrumIn->Write();
    fDptInDist->SetNameTitle("DeltaPt", "#delta p_{T} distribution, in plane");
    fDptInDist->Write();
    fDptIn->SetNameTitle("DeltaPtMatrix","#delta p_{T} matrix, in plane");
    fDptIn->Write();
    fFullResponseIn->SetNameTitle("ResponseMatrix", "Response matrix, in plane");
    fFullResponseIn->Write();
    fullResponseRebinnedIn->SetNameTitle("ResponseMatrixRebinned", "Rebinned response matrix, in plane");
    fullResponseRebinnedIn->Write();
    resizedResonseIn->Write();
    kinematicEfficiencyIn->SetNameTitle("KinematicEfficiency","Kinematic efficiency, in plane");
    kinematicEfficiencyIn->Write();
    fDetectorResponse->SetNameTitle("DetectorResponse", "Detector response matrix");
    fDetectorResponse->Write();

    fActiveDir->cd();
    TDirectoryFile* dirOut = new TDirectoryFile(Form("OutOfPlane___%s", fActiveString.Data()), Form("OutOfPlane___%s", fActiveString.Data()));
    dirOut->cd();
    Bool_t convergedOut = UnfoldSpectrum(
            resizedJetPtOut,
            resizedResonseOut,
            kinematicEfficiencyOut,
            unfoldingTemplateOut,
            fUnfoldedOut,
            TString("out"));
    fSpectrumOut->SetNameTitle("JetSpectrum", "Jet spectrum, Out plane");
    fSpectrumOut->Write();
    fDptOutDist->SetNameTitle("DeltaPt", "#delta p_{T} distribution, Out plane");
    fDptOutDist->Write();
    fDptOut->SetNameTitle("DeltaPtMatrix","#delta p_{T} matrix, Out plane");
    fDptOut->Write();
    fFullResponseOut->SetNameTitle("ResponseMatrix", "Response matrix, Out plane");
    fFullResponseOut->Write();
    fullResponseRebinnedOut->SetNameTitle("ResponseMatrixRebinned", "Rebinned response matrix, Out plane");
    fullResponseRebinnedOut->Write();
    resizedResonseOut->Write();
    kinematicEfficiencyOut->SetNameTitle("KinematicEfficiency","Kinematic efficiency, Out plane");
    kinematicEfficiencyOut->Write();
    fDetectorResponse->SetNameTitle("DetectorResponse", "Detector response matrix");
    fDetectorResponse->Write();
    
    fActiveDir->cd();
    TGraphErrors* ratio = GetRatio((TH1D*)fUnfoldedIn->Clone("hUnfolded_in"), (TH1D*)fUnfoldedOut->Clone("hUnfolded_out"));
    ratio->SetNameTitle("RatioInOutPlane", "Ratio in plane, out of plane jet spectrum");
    ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratio->GetYaxis()->SetTitle("yield IN / yield OUT");
    ratio->Write();
    fDeltaPtDeltaPhi->Write();
    fJetPtDeltaPhi->Write();
    TH1F* beta = new TH1F("UnfoldingConfiguration","UnfoldingConfiguration", 6, -.5, 5.5);
    beta->SetBinContent(1, fBeta);
    beta->GetXaxis()->SetBinLabel(1, "fBeta");
    beta->SetBinContent(2, fBetaPerDOF);
    beta->GetXaxis()->SetBinLabel(2, "fBetaPerDOF");
    beta->SetBinContent(3, fCentralityBin);
    beta->GetXaxis()->SetBinLabel(3, "fCentralityBin");
    beta->SetBinContent(4, (int)convergedIn);
    beta->GetXaxis()->SetBinLabel(4, "convergedIn");
    beta->SetBinContent(5, (int)convergedOut);
    beta->GetXaxis()->SetBinLabel(5, "convergedOut");
    beta->SetBinContent(6, (int)fAvoidRoundingError);
    beta->GetXaxis()->SetBinLabel(6, "fAvoidRoundingError");
    beta->Write();
}
//_____________________________________________________________________________
Bool_t AliJetFlowTools::UnfoldSpectrum(
        TH1D* resizedJetPt, 
        TH2D* resizedResonse,
        TH1D* kinematicEfficiency,
        TH1D* unfoldingTemplate,
        TH1D *&unfolded,
        TString suffix)
{
    /* the following is inspired by Marta Verweij's unfolding twiki */

    Bool_t converged = kFALSE;
    // step 1) clone all input histograms
    TH1D *fh1RawJetsCurrentOrig = (TH1D*)resizedJetPt->Clone(Form("fh1RawJetsCurrentOrig_%s", suffix.Data()));
    // full response matrix and kinematic efficiency
    TH2D* h2ResponseMatrixOrig = (TH2D*)resizedResonse->Clone(Form("h2ResponseMatrixOrig_%s", suffix.Data()));
    TH1D* hEfficiency = (TH1D*)kinematicEfficiency->Clone(Form("hEfficiency_%s", suffix.Data()));
    // the initial guess for the unfolded pt spectrum. 
    // equal to the folded spectrum, but different binning and ranges
    TH1D *hUnfoldedOrig = (TH1D*)unfoldingTemplate->Clone(Form("hUnfoldedOrig_%s", suffix.Data()));
    TH1D *hInitial = (TH1D*)hUnfoldedOrig->Clone(Form("hInitial_%s", suffix.Data()));
    // 'bookkeeping' histograms, this is not input
    TH2D *h2ResponseMatrix = (TH2D*)h2ResponseMatrixOrig->Clone(Form("h2ResponseMatrix_%s", suffix.Data()));
    TH1D *hUnfolded = (TH1D*)hUnfoldedOrig->Clone(Form("hUnfolded_%s", suffix.Data()));
    hUnfolded->Reset();
    TH1D *fh1RawJetsCurrent = (TH1D*)fh1RawJetsCurrentOrig->Clone(Form("fh1RawJetsCurrent_%s", suffix.Data()));
    // step 2) start the unfolding routine
    AliUnfolding::SetUnfoldingMethod(AliUnfolding::kChi2Minimization);
    AliUnfolding::SetNbins(fh1RawJetsCurrent->GetNbinsX(), hUnfolded->GetNbinsX());
    AliUnfolding::SetChi2Regularization(AliUnfolding::kLogLog,fBeta);
    AliUnfolding::SetMinuitStepSize(1.);
    AliUnfolding::SetMinuitPrecision(1e-6);
    AliUnfolding::SetMinuitMaxIterations(100000);
    AliUnfolding::SetMinuitStrategy(2.);
    AliUnfolding::SetDebug(1);
    //Unfolding loop starts here
    Int_t iret = -1; Int_t niter = 0;
    while(iret < 0 && niter < 100) {
        if (niter > 0) hInitial = (TH1D*)hUnfolded->Clone(Form("hInitial_%s", suffix.Data()));
        iret = AliUnfolding::Unfold(h2ResponseMatrix, hEfficiency,fh1RawJetsCurrent ,hInitial ,hUnfolded);
        niter++;
    }
    Int_t errMat(gMinuit->fISW[1]);
    if(iret==0 && errMat==3) {
        TVirtualFitter *fitter(TVirtualFitter::GetFitter());
        double arglist[1] = {0.};
        fitter->ExecuteCommand("SHOW COV", arglist, 1);
        TMatrixD covmatrix(hUnfolded->GetNbinsX(),hUnfolded->GetNbinsX(),fitter->GetCovarianceMatrix());
        TMatrixD *pearson((TMatrixD*)CalculatePearsonCoefficients(&covmatrix));
        TH2D *hPearson(new TH2D(*pearson));
        hPearson->SetYTitle("bin number");
        hPearson->SetXTitle("bin number");
        hPearson->SetNameTitle("PearsonCoefficients", Form("Pearson coefficients, %s plane", suffix.Data()));
        hPearson->SetMinimum(-1.);
        hPearson->SetMaximum(1.);
        hPearson->GetZaxis()->SetLabelSize(0.06);
        hPearson->Write();
        converged = kTRUE;
    }  
    AliAnaChargedJetResponseMaker *corr(new AliAnaChargedJetResponseMaker());
    TH1D *hFolded(corr->MultiplyResponseGenerated(hUnfolded,h2ResponseMatrix,hEfficiency));
    hFolded->SetNameTitle("RefoldedSpectrum", Form("Refolded jet spectrum, %s plane", suffix.Data()));
    hFolded->Write();
    hUnfolded->SetNameTitle("UnfoldedSpectrum", Form("Unfolded jet spectrum, %s plane", suffix.Data()));
    hUnfolded->Write();
    unfolded = hUnfolded;
    TGraphErrors* ratio(GetRatio(hFolded,fh1RawJetsCurrent, kTRUE, fBinsTrue->At(0), fBinsTrue->At(fBinsTrue->GetSize()-1)));
    ratio->SetNameTitle("RatioRefoldedMeasured", Form("Ratio refolded and measured spectrum %s plane", suffix.Data()));
    ratio->Write();
    
    return converged;
}
//_____________________________________________________________________________
Bool_t AliJetFlowTools::PrepareForUnfolding()
{
    // prepare for unfolding
    if(!fInputList) {
        printf(" AliJetFlowTools::PrepareForUnfolding() fInputList not found \n - Set a list using AliJetFlowTools::SetInputList() \n");
        return kFALSE;
    }
    if(!fDetectorResponse) {
        printf(" AliJetFlowTools::PrepareForUnfolding() fDetectorResponse not found \n - Set detector response using AliJetFlowTools::SetDetectorResponse() \n ");
        return kFALSE;
    }
    // check if the pt bin for true and rec have been set
    if(!fBinsTrue) {
        printf(" \n AliJetFlowTools::PrepareForUnfolding() \n - warning: fBinsTrue nog set, using defaults ! \n");
        const Double_t ptBins[] = {20, 22.5, 25, 27.5, 30,32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50,52.5, 55, 57.5, 60};
        fBinsTrue = new TArrayD(sizeof(ptBins)/sizeof(ptBins[0]), ptBins);
    }
    if(!fBinsRec) {
        printf(" \n AliJetFlowTools::PrepareForUnfolding() \n - warning: fBinsRec not set, using defaults ! \n");
        Double_t binsY[81];
        for(Int_t i(0); i < 81; i++) binsY[i] = (double)(30+i);
        fBinsRec = new TArrayD(sizeof(binsY)/sizeof(binsY[0]), binsY);
    }
    // extract the spectra
    TString spectrumName(Form("fHistJetPsi2Pt_%i", fCentralityBin));
    fJetPtDeltaPhi = ((TH2D*)fInputList->FindObject(spectrumName.Data()));
    if(!fJetPtDeltaPhi) {
        printf(" Couldn't find spectrum %s ! \n", spectrumName.Data());
        return kFALSE;
    }
    // in plane spectrum
    fSpectrumIn = fJetPtDeltaPhi->ProjectionY(Form("_py_ina_%s", spectrumName.Data()), 1, 10);
    fSpectrumIn->Add(fJetPtDeltaPhi->ProjectionY(Form("_py_inb_%s", spectrumName.Data()), 31, 40));
    // out of plane spectrum
    fSpectrumOut = fJetPtDeltaPhi->ProjectionY(Form("_py_out_%s", spectrumName.Data()), 11, 30);

    // extract the delta pt matrices
    TString deltaptName(Form("fHistDeltaPtDeltaPhi2_%i", fCentralityBin));
    fDeltaPtDeltaPhi = ((TH2D*)fInputList->FindObject(deltaptName.Data()));
    if(!fDeltaPtDeltaPhi) {
        printf(" Couldn't find delta pt matrix %s ! \n", deltaptName.Data());
    }
    if(fdPhidPt) {
        // in plane delta pt distribution
        fDptInDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_ina_%s", deltaptName.Data()), 1, 10);
        fDptInDist->Add(fDeltaPtDeltaPhi->ProjectionY(Form("_py_inb_%s", deltaptName.Data()), 31, 40));
        // out of plane delta pt distribution
        fDptOutDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_out_%s", deltaptName.Data()), 11, 30);
    } else {
        fDptInDist = fDeltaPtDeltaPhi->ProjectionY("_pyFULLa");
        fDptOutDist = fDeltaPtDeltaPhi->ProjectionY("_pyFULLb");
    }
    
    // create a rec - true smeared response matrix
    TMatrixD* rfIn = new TMatrixD(-50, 249, -50, 249);
    for(Int_t j(-50); j < 250; j++) {   // loop on pt true slices j
        Bool_t skip = kFALSE;
        for(Int_t k(-50); k < 250; k++) {       // loop on pt gen slices k
            (*rfIn)(k, j) = (skip) ? 0. : fDptInDist->GetBinContent(fDptInDist->GetXaxis()->FindBin(k-j));
            if(fAvoidRoundingError && k > j && TMath::AreEqualAbs(fDptInDist->GetBinContent(fDptInDist->GetXaxis()->FindBin(k-j)), 0, 1e-8)) skip = kTRUE;
        }
    }
    TMatrixD* rfOut = new TMatrixD(-50, 249, -50, 249);
    for(Int_t j(-50); j < 250; j++) {   // loop on pt true slices j
        Bool_t skip = kFALSE;
        for(Int_t k(-50); k < 250; k++) {       // loop on pt gen slices k
            (*rfOut)(k, j) = (skip) ? 0. : fDptOutDist->GetBinContent(fDptOutDist->GetXaxis()->FindBin(k-j));
            if(fAvoidRoundingError && k > j && TMath::AreEqualAbs(fDptOutDist->GetBinContent(fDptOutDist->GetXaxis()->FindBin(k-j)), 0, 1e-8)) skip = kTRUE;
        }
    }
    fDptIn = new TH2D(*rfIn);
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityBin), Form("dpt_response_INPLANE_%i", fCentralityBin));
    fDptIn->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");
    fDptOut = new TH2D(*rfOut);
    fDptOut->SetNameTitle(Form("dpt_response_OUTOFPLANE_%i", fCentralityBin), Form("dpt_response_OUTOFPLANE_%i", fCentralityBin));
    fDptOut->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptOut->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");

    return kTRUE;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix) {
    // resize the x-axis of a th1d
    if(!histo) {
        printf(" > ResizeXaxisTH!D:: fatal error, NULL pointer passed < \n");
        return NULL;
    } 
    // see how many bins we need to copy
    TH1D* resized = new TH1D(Form("%s_resized_%s", histo->GetName(), suffix.Data()), Form("%s_resized_%s", histo->GetName(), suffix.Data()), up-low, (double)low, (double)up);
    // low is the bin number of the first new bin
    Int_t l = histo->GetXaxis()->FindBin(low);
    // set the values
    Double_t _x(0), _xx(0);
    for(Int_t i(0); i < up-low; i++) {
        _x = histo->GetBinContent(l+i);
        _xx=histo->GetBinError(l+i);
        resized->SetBinContent(i+1, _x);
        resized->SetBinError(i+1, _xx);
    }
    return resized;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::ResizeYaxisTH2D(TH2D* histo, TArrayD* x, TArrayD* y, TString suffix) {
    // resize the y-axis of a th2d
    if(!histo) {
        printf(" > ResizeYaxisTH2D:: fatal error, NULL pointer passed < \n");
        return NULL;
    } 
    // see how many bins we need to copy
    TH2D* resized = new TH2D(Form("%s_resized_%s", histo->GetName(), suffix.Data()), Form("%s_resized_%s", histo->GetName(), suffix.Data()), x->GetSize()-1, x->GetArray(), y->GetSize()-1, y->GetArray());
    // assume only the y-axis has changed
    // low is the bin number of the first new bin
    Int_t low = histo->GetYaxis()->FindBin(y->At(0));
    // set the values
    Double_t _x(0), _xx(0);
    for(Int_t i(0); i < x->GetSize(); i++) {
        for(Int_t j(0); j < y->GetSize(); j++) {
            _x = histo->GetBinContent(i, low+j);
            _xx=histo->GetBinError(i, low+1+j);
            resized->SetBinContent(i, j, _x);
            resized->SetBinError(i, j, _xx);
        }
    }
    return resized;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::NormalizeTH2D(TH2D* histo) {
    // general method to normalize all vertical slices of a th2 to unity
    // i.e. get a probability matrix
    if(!histo) {
        printf(" > NormalizeTH2D:: fatal error, NULL pointer passed < \n");
        return NULL;
    }
    Int_t binsX = histo->GetXaxis()->GetNbins();
    Int_t binsY = histo->GetYaxis()->GetNbins();
    
    // normalize all slices in x
    for(Int_t i(0); i < binsX; i++) {   // for each vertical slice
        Double_t weight = 0;
        for(Int_t j(0); j < binsY; j++) {       // loop over all the horizontal components
            weight+=histo->GetBinContent(i+1, j+1);
        }       // now we know the total weight
        for(Int_t j(0); j < binsY; j++) {
            if (weight <= 0 ) continue;
            histo->SetBinContent(1+i, j+1, histo->GetBinContent(1+i, j+1)/weight);
            histo->SetBinError(  1+i, j+1, histo->GetBinError(  1+i, j+1)/weight);
        }
    }
    return histo;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::GetUnfoldingTemplate(TH1D* histo, TArrayD* bins, TString suffix) {
    // return a TH1D with the supplied histogram rebinned to the supplied bins
    // this histogram will be used as a startng point for the chi^2 minimization
    //
    if(!histo || !bins) {
        printf(" > RebinTH2D:: fatal error, NULL pointer passed < \n");
        return NULL;
    }
    // create the output histo
    TString name = histo->GetName();
    name+="_template";
    name+=suffix;
    TH1D* rebinned = new TH1D(name.Data(), name.Data(), bins->GetSize()-1, bins->GetArray());
    for(Int_t i(0); i < histo->GetXaxis()->GetNbins(); i++) {
        for(Int_t j(0); j < histo->GetBinContent(i+1); j++) {
            if(i-100 > 1) rebinned->Fill(i-100);
        }
    }
    return rebinned;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::RebinTH2DX(TH2D* histo, TArrayD* bins, TString suffix) {
    // rebin a TH2D to a binning scheme supplied in a TArrayD
    // as a rebinning weight a tsallis fit to pythia is used
    // the returned histogram is new

    // see if input makes sense
    if(!histo || !bins) {
        printf(" > RebinTH2D:: fatal error, NULL pointer passed < \n");
        return NULL;
    }
    // set up some common variables
    TString name = histo->GetName();
    name+="_rebinned_";
    name+=suffix;
    Double_t arrY[301];
    TF1* weightFunction = new TF1("weightFunction", "x*TMath::Power(1.+(1./(8.*0.9))*x, -8.)", 0 ,150);

    for(Int_t i(0); i < 301; i++) arrY[i] = i-50;
    // create the output histo
    TH2D* rebinned = new TH2D(name.Data(), name.Data(), bins->GetSize()-1, bins->GetArray(), 300, arrY);
    rebinned->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    rebinned->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    // loop over the original bins and do the rebinning

    for(Int_t y(0); y < 301; y++) {   // loop over pt gen bins, we don't merge these
        for(Int_t newBin(0); newBin < bins->GetSize()-1; newBin++) {  // loop over desired pt true bins
            Int_t newBinMin = histo->GetXaxis()->FindBin(bins->At(newBin));       // bin with lower border in pt
            Int_t newBinMax = histo->GetXaxis()->FindBin(bins->At(newBin+1));     // bin with upper border in pt
            Double_t temp(0), weightsum(0);
            for(Int_t ll(newBinMin); ll < newBinMax; ll++) {       // fin bins in pt rec
                // loop over the fine, original bins in x
                Double_t cen = histo->GetXaxis()->GetBinCenter(ll);   // center of current bin in x
                weightsum+=weightFunction->Eval(cen);
                temp+=weightFunction->Eval(cen)*histo->GetBinContent(ll+1, y+1);
            }
            rebinned->SetBinContent(newBin+1, y+1, temp);
        }
    }
    return rebinned;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::RebinTH2DY(TH2D* histo, TArrayD* bins) {
    // rebin a TH2D to a binning scheme supplied in a TArrayD
    // as a rebinning weight a tsallis fit to pythia is used
    // the returned histogram is new

    // see if input makes sense
    if(!histo || !bins) {
        printf(" > RebinTH2D:: fatal error, NULL pointer passed < \n");
        return NULL;
    }
    // set up some common variables
    TString name = histo->GetName();
    name+="_rebinned_";
    Double_t arrY[301];
    TF1* weightFunction = new TF1("weightFunction", "x*TMath::Power(1.+(1./(8.*0.9))*x, -8.)", 0 ,150);

    for(Int_t i(0); i < 301; i++) arrY[i] = i-50;
    // create the output histo
    TH2D* rebinned = new TH2D(name.Data(), name.Data(), bins->GetSize()-1, bins->GetArray(), 300, arrY);
    rebinned->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    rebinned->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    // loop over the original bins and do the rebinning

    for(Int_t y(0); y < 301; y++) {   // loop over pt gen bins, we don't merge these
        for(Int_t newBin(0); newBin < bins->GetSize()-1; newBin++) {  // loop over desired pt true bins
            Int_t newBinMin = histo->GetXaxis()->FindBin(bins->At(newBin));       // bin with lower border in pt
            Int_t newBinMax = histo->GetXaxis()->FindBin(bins->At(newBin+1));     // bin with upper border in pt
            Double_t temp(0), weightsum(0);
            for(Int_t ll(newBinMin); ll < newBinMax; ll++) {       // fin bins in pt rec
                // loop over the fine, original bins in x
                Double_t cen = histo->GetXaxis()->GetBinCenter(ll);   // center of current bin in x
                weightsum+=weightFunction->Eval(cen);
                temp+=weightFunction->Eval(cen)*histo->GetBinContent(ll+1, y+1);
            }
            rebinned->SetBinContent(newBin+1, y+1, temp);
        }
    }
    return rebinned;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::MatrixMultiplicationTH2D(TH2D* A, TH2D* B, TString name) {
    // general matrix multiplication
    if(!A || !B) {
       printf(" > MatrixMultiplicationTH2D:: fatal error, NULL pointer passed < \n");
       return NULL;
    } 
    // step A -> see if the matrices are both square and of equal dimensions
    Int_t aX(A->GetXaxis()->GetNbins());
    Int_t aY(A->GetYaxis()->GetNbins());
    Int_t bX(B->GetXaxis()->GetNbins());
    Int_t bY(B->GetYaxis()->GetNbins());
    if(!(aX == aY && aX == bX && aX == bY)) {
        printf(" > MatrixMultiplicationTH2D:: Error, matrices with incompatible dimensions passed ! < \n");
        printf("   - matrix A (%i, %i) \n", aX, aY);
        printf("   - matrix B (%i, %i) \n", bX, bY);
        return NULL;
    }

    // step B -> setup a new matrix template to store the results
    TH2D* C = (TH2D*)A->Clone();
    C->SetNameTitle(name.Data(), name.Data());
    C->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    C->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");

    // step C -> nested loop multiplication
    printf ( " >> parsing matrix data, may take some time (depending on matrix dimensions) << \n");
    for(Int_t i(0); i < aX; i++) {      // x coordinate of target matrix
        for(Int_t j(0); j < aY; j++) {  // y coordinate of target matrix

            // matrix multiplications are ugly by definition. 
            // the thing to keep in mind: the coordinate of the target matrix is
            // (x, y) is the dot product of row x of matrix A with column y of matrix B
            // so in this nexted loop, we need to fill point i, j of the target matrix
            // with the dot product of 
            // matrix A) row i \cdotp matrix B) column j
            // looping through a row means that y remains constant
            // looping through a column means that x remains constant

            Double_t value(0);      // placeholder for target value
            for(Int_t x(0); x < aX; x++) value += A->GetBinContent(x, i) * B->GetBinContent(j, x);
            C->SetBinContent(i, j);
        }
    }
    return C;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::NormalizeTH1D(TH1D* histo, Double_t scale) {
    // normalize a th1d to a certain scale
     Double_t integral = histo->Integral()*scale;
     if (integral > 0 && scale == 1.) histo->Scale(1./integral, "width");
     else if (scale != 1.) histo->Scale(1./scale, "width");
     else printf(" > Histogram integral < 0, cannot normalize \n");
     return histo;
}
//_____________________________________________________________________________
TMatrixD* AliJetFlowTools::CalculatePearsonCoefficients(TMatrixD* covmat) 
{
    // Calculate pearson coefficients from covariance matrix
    TMatrixD *pearsonCoefs((TMatrixD*)covmat->Clone("pearsonCoefs"));
    Int_t nrows(covmat->GetNrows()), ncols (covmat->GetNcols());
    Double_t pearson(0.);
    if(nrows==0 && ncols==0) return 0x0;
    for(int row = 0; row<nrows; row++) {
        for(int col = 0; col<ncols; col++) {
        if((*covmat)(row,row)!=0. && (*covmat)(col,col)!=0.) pearson = (*covmat)(row,col)/TMath::Sqrt((*covmat)(row,row)*(*covmat)(col,col));
        (*pearsonCoefs)(row,col) = pearson;
        }
    }
    return pearsonCoefs;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::GetRatio(TH1 *h1, TH1* h2, Bool_t appendFit, Double_t low, Double_t up) 
{
    TGraphErrors *gr(new TGraphErrors());
    Int_t j(0);
    Float_t binCent(0.), ratio(0.), error2(0.), binWidth(0.);
    for(Int_t i(1); i < h1->GetNbinsX(); i++) {
        binCent = h1->GetXaxis()->GetBinCenter(i);
        j = h2->FindBin(binCent);
        binWidth = h1->GetXaxis()->GetBinWidth(i);
        if(h2->GetBinContent(j) > 0.) {
            ratio = h1->GetBinContent(i)/h2->GetBinContent(j);
            Double_t A = 1./h2->GetBinContent(j)*h1->GetBinError(i);
            error2 = A*A;
            gr->SetPoint(gr->GetN(),binCent,ratio);
            gr->SetPointError(gr->GetN()-1,0.5*binWidth,TMath::Sqrt(error2));
        }
    }
    if(appendFit) {
        TF1* fit = new TF1("name", "pol1", low, up);
        gr->Fit(fit);
    }
    return gr;
}
//_____________________________________________________________________________
