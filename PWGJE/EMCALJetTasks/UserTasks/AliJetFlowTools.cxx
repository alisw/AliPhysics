/************************************************************************ 
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
//
// The task uses input from two analysis tasks:
// - $ALICE_PHYSICS/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskJetV2.cxx
//   used to retrieve jet spectra and delta pt distributions
// - $ALICE_PHYSICS/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskJetMatching.cxx
//   used to construct the detector response function
// and unfolds jet spectra with respect to the event plane. The user can choose
// different alrogithms for unfolding which are available in (ali)root. RooUnfold 
// libraries must be present on the system 
// ( see http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html ).
// 
// The weak spot of this class is the function PrepareForUnfolding, which will read
// output from two output files and expects histograms with certain names and binning. 
// Unfolding methods itself are general and should be able to handle any input, therefore one
// can forgo the PrepareForUnfolding method, and supply necessary input information via the 
// SetRawInput() method
//
// to see an example of how to use this class, see $ALICE_PHYSICS/PWGCF/FLOW/macros/jetFlowTools.C

// root includes
#include "TH1.h"
#include "TF2.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TArrayD.h"
#include "TList.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
// aliroot includes
#include "AliUnfolding.h"
#include "AliAnaChargedJetResponseMaker.h"
// class includes
#include "AliJetFlowTools.h"
// roo unfold includes (make sure you have these available on your system)
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBayes.h"
#include "TSVDUnfold.h"

using namespace std;
//_____________________________________________________________________________
AliJetFlowTools::AliJetFlowTools() :
    fListPrefix         (-1),
    fResponseMaker      (new AliAnaChargedJetResponseMaker()),
    fRMS                (kTRUE),
    fSymmRMS            (0),
    fConstantUE         (kTRUE),
    fRho0               (kFALSE),
    fBootstrap          (kFALSE),
    fPower              (new TF1("fPower","[0]*TMath::Power(x,-([1]))",0.,300.)),
    fSaveFull           (kTRUE),
    fActiveString       (""),
    fActiveDir          (0x0),
    fInputList          (0x0),
    fRefreshInput       (kTRUE),
    fOutputFileName     ("UnfoldedSpectra.root"),
    fOutputFile         (0x0),
    fCentralityArray    (0x0),
    fMergeBinsArray     (0x0),
    fCentralityWeights  (0x0),
    fMergeWithList      (0x0),
    fMergeWithCen       (-1),
    fMergeWithWeight    (1.),
    fDetectorResponse   (0x0),
    fJetFindingEff      (0x0),
    fBetaIn             (.1),
    fBetaOut            (.1),
    fBayesianIterIn     (4),
    fBayesianIterOut    (4),
    fBayesianSmoothIn   (0.),
    fBayesianSmoothOut  (0.),
    fAvoidRoundingError (kFALSE),
    fUnfoldingAlgorithm (kChi2),
    fPrior              (kPriorMeasured),
    fPriorUser          (0x0),
    fBinsTrue           (0x0),
    fBinsRec            (0x0),
    fBinsTruePrior      (0x0),
    fBinsRecPrior       (0x0),
    fSVDRegIn           (5),
    fSVDRegOut          (5),
    fSVDToy             (kTRUE),
    fJetRadius          (0.3),
    fEventCount         (-1),
    fNormalizeSpectra   (kFALSE),
    fSmoothenPrior      (kFALSE),
    fFitMin             (60.),
    fFitMax             (300.),
    fFitStart           (75.),
    fSmoothenCounts     (kTRUE),
    fTestMode           (kFALSE),
    fRawInputProvided   (kFALSE),
    fEventPlaneRes      (.63),
    fUseDetectorResponse(kTRUE),
    fUseDptResponse     (kTRUE),
    fTrainPower         (kTRUE),
    fDphiUnfolding      (kTRUE),
    fDphiDptUnfolding   (kFALSE),
    fExLJDpt            (kTRUE),
    fTitleFontSize      (-999.),
    fRMSSpectrumIn      (0x0),
    fRMSSpectrumOut     (0x0),
    fRMSRatio           (0x0),
    fRMSV2              (0x0),
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
    fPivot              (55.),
    fSubdueError        (kTRUE),
    fUnfoldedSpectrumIn (0x0),
    fUnfoldedSpectrumOut(0x0),
    fHarmonic(2) { // class constructor
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
        fResponseMaker->SetRMMergeWeightFunction(new TF1("weightFunction", "x*TMath::Power(1.+(1./(8.*0.9))*x, -8.)", 0, 200));
        for(Int_t i(0); i < fPower->GetNpar(); i++) fPower->SetParameter(i, 0.);
}
//_____________________________________________________________________________
void AliJetFlowTools::Make(TH1* customIn, TH1* customOut) {
    // core function of the class
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   if(fDphiDptUnfolding) {
        // to extract the yield as function of Dphi, Dpt - experimental
        MakeAU();
        return;
    }
    // 1) rebin the raw output of the jet task to the desired binnings
    // 2) calls the unfolding routine
    // 3) writes output to file
    // can be repeated multiple times with different configurations

    // 1) manipulation of input histograms
    // check if the input variables are present
    if(fRefreshInput) {
        if(!PrepareForUnfolding(customIn, customOut)) {
            printf(" AliJetFlowTools::Make() Fatal error \n - couldn't prepare for unfolding ! \n");
            return;
        }
    }
    // 1a) resize the jet spectrum according to the binning scheme in fBinsTrue
    //     parts of the spectrum can end up in over or underflow bins

    // if bootstrap mode is kTRUE, resample the underlying distributions
    // FIXME think about resampling the rebinned results or raw results, could lead to difference
    // in smoothness of tail of spectrum (which is probably not used in any case, but still ... )
/*
    if(fBootstrap) {
        // resample but leave original spectra intact for the next unfolding round
        fSpectrumIn = reinterpret_cast<TH1D*>(Bootstrap(fSpectrumIn, kFALSE));
        fSpectrumOut = reinterpret_cast<TH1D*>(Bootstrap(fSpectrumOut, kFALSE));
    }
*/
    TH1D* measuredJetSpectrumIn  = RebinTH1D(fSpectrumIn, fBinsRec, TString("resized_in_"), kFALSE);
    TH1D* measuredJetSpectrumOut = RebinTH1D(fSpectrumOut, fBinsRec,  TString("resized_out_"), kFALSE);

    if(fBootstrap) {
        measuredJetSpectrumIn = reinterpret_cast<TH1D*>(Bootstrap(measuredJetSpectrumIn, kFALSE));
        measuredJetSpectrumOut = reinterpret_cast<TH1D*>(Bootstrap(measuredJetSpectrumOut, kFALSE));
    }
    // for now do it BEFORE as after gives an issue in Rebin function (counts are wrong)


    
    // 1b) resize the jet spectrum to 'true' bins. can serve as a prior and as a template for unfolding
    // the template will be used as a prior for the chi2 unfolding
    TH1D* measuredJetSpectrumTrueBinsIn  = RebinTH1D(fSpectrumIn, fBinsTrue, TString("in"), kFALSE);
    TH1D* measuredJetSpectrumTrueBinsOut = RebinTH1D(fSpectrumOut, fBinsTrue, TString("out"), kFALSE);
    // get the full response matrix from the dpt and the detector response
    fDetectorResponse = NormalizeTH2D(fDetectorResponse);
    // get the full response matrix. if test mode is chosen, the full response is replace by a unity matrix
    // so that unfolding should return the initial spectrum
    if(!fTestMode) {
        if(fUseDptResponse && fUseDetectorResponse) {
            fFullResponseIn = MatrixMultiplication(fDptIn, fDetectorResponse);
            fFullResponseOut = MatrixMultiplication(fDptOut, fDetectorResponse);
        } else if (fUseDptResponse && !fUseDetectorResponse) {
            fFullResponseIn = fDptIn;
            fFullResponseOut = fDptOut;
        } else if (!fUseDptResponse && fUseDetectorResponse) {
            fFullResponseIn = fDetectorResponse;
            fFullResponseOut = fDetectorResponse;
        } else if (!fUseDptResponse && !fUseDetectorResponse && !fUnfoldingAlgorithm == AliJetFlowTools::kNone) {
            printf(" > No response, exiting ! < \n" );
            return;
        }
    } else {
        fFullResponseIn = GetUnityResponse(fBinsTrue, fBinsRec, TString("in"));
        fFullResponseOut = GetUnityResponse(fBinsTrue, fBinsRec, TString("out"));
    }
    // normalize each slide of the response to one
    NormalizeTH2D(fFullResponseIn);
    NormalizeTH2D(fFullResponseOut);
    // resize to desired binning scheme
    TH2D* resizedResponseIn  = RebinTH2D(fFullResponseIn, fBinsTrue, fBinsRec, TString("in"));
    TH2D* resizedResponseOut = RebinTH2D(fFullResponseOut, fBinsTrue, fBinsRec, TString("out"));
    // get the kinematic efficiency
    TH1D* kinematicEfficiencyIn  = resizedResponseIn->ProjectionX();
    kinematicEfficiencyIn->SetNameTitle("kin_eff_IN","kin_eff_IN");
    TH1D* kinematicEfficiencyOut = resizedResponseOut->ProjectionX();
    kinematicEfficiencyOut->SetNameTitle("kin_eff_OUT", "kin_eff_OUT");
    // suppress the errors 
    for(Int_t i(0); i < kinematicEfficiencyOut->GetXaxis()->GetNbins(); i++) {
        kinematicEfficiencyIn->SetBinError(1+i, 0.);
        kinematicEfficiencyOut->SetBinError(1+i, 0.);
    }
    TH1D* jetFindingEfficiency(0x0);
    if(fJetFindingEff) {
        jetFindingEfficiency = ProtectHeap(fJetFindingEff);
        jetFindingEfficiency->SetNameTitle(Form("%s_coarse", jetFindingEfficiency->GetName()), Form("%s_coarse", jetFindingEfficiency->GetName()));
        jetFindingEfficiency = RebinTH1D(jetFindingEfficiency, fBinsTrue);
    }
    // 2, 3) call the actual unfolding. results and transient objects are stored in a dedicated TDirectoryFile
    TH1D* unfoldedJetSpectrumIn(0x0);
    TH1D* unfoldedJetSpectrumOut(0x0); 
    fActiveDir->cd();                   // select active dir
    TDirectoryFile* dirIn = new TDirectoryFile(Form("InPlane___%s", fActiveString.Data()), Form("InPlane___%s", fActiveString.Data()));
    dirIn->cd();                        // select inplane subdir
    // do the inplane unfolding
    unfoldedJetSpectrumIn = UnfoldWrapper(
        measuredJetSpectrumIn,
        resizedResponseIn,
        kinematicEfficiencyIn,
        measuredJetSpectrumTrueBinsIn,
        TString("in"),
        jetFindingEfficiency);
    resizedResponseIn->SetNameTitle("ResponseMatrixIn", "response matrix in plane");
    resizedResponseIn->SetXTitle("p_{T, jet}^{true} (GeV/#it{c})");
    resizedResponseIn->SetYTitle("p_{T, jet}^{rec} (GeV/#it{c})");
    resizedResponseIn = ProtectHeap(resizedResponseIn);
    resizedResponseIn->Write();
    kinematicEfficiencyIn->SetNameTitle("KinematicEfficiencyIn","Kinematic efficiency, in plane");
    kinematicEfficiencyIn = ProtectHeap(kinematicEfficiencyIn);
    kinematicEfficiencyIn->Write();
    fDetectorResponse->SetNameTitle("DetectorResponse", "Detector response matrix");
    fDetectorResponse = ProtectHeap(fDetectorResponse, kFALSE);
    fDetectorResponse->Write();
    // optional histograms
    if(fSaveFull) {
        fSpectrumIn->SetNameTitle("[ORIG]JetSpectrum", "[INPUT] Jet spectrum, in plane");
        fSpectrumIn->Write();
        fDptInDist->SetNameTitle("[ORIG]DeltaPt", "#delta p_{T} distribution, in plane");
        fDptInDist->Write();
        fDptIn->SetNameTitle("[ORIG]DeltaPtMatrix","#delta p_{T} matrix, in plane");
        fDptIn->Write();
        fFullResponseIn->SetNameTitle("ResponseMatrix", "Response matrix, in plane");
        fFullResponseIn->Write();
    }
    fActiveDir->cd();
    if(fDphiUnfolding) {
        TDirectoryFile* dirOut = new TDirectoryFile(Form("OutOfPlane___%s", fActiveString.Data()), Form("OutOfPlane___%s", fActiveString.Data()));
        dirOut->cd();
        // do the out of plane unfolding
        unfoldedJetSpectrumOut = UnfoldWrapper(
                    measuredJetSpectrumOut,
                    resizedResponseOut,
                    kinematicEfficiencyOut,
                    measuredJetSpectrumTrueBinsOut,
                    TString("out"),
                    jetFindingEfficiency);
        resizedResponseOut->SetNameTitle("ResponseMatrixOut", "response matrix in plane");
        resizedResponseOut->SetXTitle("p_{T, jet}^{true} (GeV/#it{c})");
        resizedResponseOut->SetYTitle("p_{T, jet}^{rec} (GeV/#it{c})");
        resizedResponseOut = ProtectHeap(resizedResponseOut);
        resizedResponseOut->Write();
        kinematicEfficiencyOut->SetNameTitle("KinematicEfficiencyOut","Kinematic efficiency, Out plane");
        kinematicEfficiencyOut = ProtectHeap(kinematicEfficiencyOut);
        kinematicEfficiencyOut->Write();
        fDetectorResponse->SetNameTitle("DetectorResponse", "Detector response matrix");
        fDetectorResponse = ProtectHeap(fDetectorResponse, kFALSE);
        fDetectorResponse->Write();
        if(jetFindingEfficiency) jetFindingEfficiency->Write();
        // optional histograms
        if(fSaveFull) {
            fSpectrumOut->SetNameTitle("[ORIG]JetSpectrum", "[INPUT]Jet spectrum, Out plane");
            fSpectrumOut->Write();
            fDptOutDist->SetNameTitle("[ORIG]DeltaPt", "#delta p_{T} distribution, Out plane");
            fDptOutDist->Write();
            fDptOut->SetNameTitle("[ORIG]DeltaPtMatrix","#delta p_{T} matrix, Out plane");
            fDptOut->Write();
            fFullResponseOut->SetNameTitle("[ORIG]ResponseMatrix", "Response matrix, Out plane");
            fFullResponseOut->Write();
        }

        // write general output histograms to file
        fActiveDir->cd();
        if(unfoldedJetSpectrumIn && unfoldedJetSpectrumOut && unfoldedJetSpectrumIn && unfoldedJetSpectrumOut) {
            TGraphErrors* ratio(GetRatio((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_in"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_out")));
            if(ratio) {
                ratio->SetNameTitle("RatioInOutPlane", "Ratio in plane, out of plane jet spectrum");
                ratio->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                ratio->GetYaxis()->SetTitle("yield IN / yield OUT");
                ratio = ProtectHeap(ratio);
                ratio->Write();
                // write histo values to RMS files if both routines converged
                // input values are weighted by their uncertainty
                for(Int_t i(0); i < ratio->GetXaxis()->GetNbins(); i++) {
                    if(unfoldedJetSpectrumIn->GetBinError(i+1) > 0) fRMSSpectrumIn->Fill(fRMSSpectrumIn->GetBinCenter(i+1), unfoldedJetSpectrumIn->GetBinContent(i+1), 1./TMath::Power(unfoldedJetSpectrumIn->GetBinError(i+1), 2.));
                    if(unfoldedJetSpectrumOut->GetBinError(i+1) > 0) fRMSSpectrumOut->Fill(fRMSSpectrumOut->GetBinCenter(i+1), unfoldedJetSpectrumOut->GetBinContent(i+1), 1./TMath::Power(unfoldedJetSpectrumOut->GetBinError(i+1), 2.));
                    if(unfoldedJetSpectrumOut->GetBinContent(i+1) > 0) fRMSRatio->Fill(fRMSSpectrumIn->GetBinCenter(i+1), unfoldedJetSpectrumIn->GetBinContent(i+1) / unfoldedJetSpectrumOut->GetBinContent(i+1));
               }
            }
            TGraphErrors* v2(GetV2((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_inv2"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_outv2"), fEventPlaneRes));
            if(v2) {
                v2->SetNameTitle("v2", "v_{2} from different in, out of plane yield");
                v2->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                v2->GetYaxis()->SetTitle("v_{2}");
                v2 = ProtectHeap(v2);
                v2->Write();
            }
        } else if (unfoldedJetSpectrumOut && unfoldedJetSpectrumIn) {
            TGraphErrors* ratio(GetRatio((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_in"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_out"), TString(""), kTRUE, fBinsRec->At(fBinsRec->GetSize()-1)));
            if(ratio) {
                ratio->SetNameTitle("[NC]RatioInOutPlane", "[NC]Ratio in plane, out of plane jet spectrum");
                ratio->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                ratio->GetYaxis()->SetTitle("yield IN / yield OUT");
                ratio = ProtectHeap(ratio);
                ratio->Write();
            }
            TGraphErrors* v2(GetV2((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_inv2"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_outv2"), fEventPlaneRes));
             if(v2) {
                v2->SetNameTitle("v2", "v_{2} from different in, out of plane yield");
                v2->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                v2->GetYaxis()->SetTitle("v_{2}");
                v2 = ProtectHeap(v2);
                v2->Write();
            }
        }
    }   // end of if(fDphiUnfolding)
    fDeltaPtDeltaPhi->Write();
    unfoldedJetSpectrumIn->Sumw2();
    fUnfoldedSpectrumIn = ProtectHeap(unfoldedJetSpectrumIn, kFALSE);
    unfoldedJetSpectrumIn->Write();
    unfoldedJetSpectrumOut->Sumw2();
    fUnfoldedSpectrumOut = ProtectHeap(unfoldedJetSpectrumOut, kFALSE);
    unfoldedJetSpectrumOut->Write();
    fJetPtDeltaPhi->Write();
    // save the current state of the unfolding object
    SaveConfiguration(unfoldedJetSpectrumIn ? kTRUE : kFALSE, unfoldedJetSpectrumOut ? kTRUE : kFALSE);
    TH1D* unfoldedJetSpectrumInForSub((TH1D*)unfoldedJetSpectrumIn->Clone("forSubIn"));
    TH1D* unfoldedJetSpectrumOutForSub((TH1D*)unfoldedJetSpectrumOut->Clone("forSubOut"));
    unfoldedJetSpectrumInForSub->Add(unfoldedJetSpectrumOutForSub, -1.);
    unfoldedJetSpectrumInForSub= ProtectHeap(unfoldedJetSpectrumInForSub);
    unfoldedJetSpectrumInForSub->Write();

}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldWrapper(
        const TH1D* measuredJetSpectrum,        // truncated raw jets (same binning as pt rec of response) 
        const TH2D* resizedResponse,            // response matrix
        const TH1D* kinematicEfficiency,        // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,        // unfolding template: same binning is pt gen of response
        const TString suffix,                   // suffix (in or out of plane)
        const TH1D* jetFindingEfficiency)       // jet finding efficiency
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // wrapper function to call specific unfolding routine
    TH1D* (AliJetFlowTools::*myFunction)(const TH1D*, const TH2D*, const TH1D*, const TH1D*, const TString, const TH1D*);
    // initialize functon pointer
    if(fUnfoldingAlgorithm == kChi2)                    myFunction = &AliJetFlowTools::UnfoldSpectrumChi2;
    else if(fUnfoldingAlgorithm == kBayesian)           myFunction = &AliJetFlowTools::UnfoldSpectrumBayesian;
    else if(fUnfoldingAlgorithm == kBayesianAli)        myFunction = &AliJetFlowTools::UnfoldSpectrumBayesianAli;
    else if(fUnfoldingAlgorithm == kSVD)                myFunction = &AliJetFlowTools::UnfoldSpectrumSVD;
    else if(fUnfoldingAlgorithm == kFold)               myFunction = &AliJetFlowTools::FoldSpectrum;
    else if(fUnfoldingAlgorithm == kNone) {
        TH1D* clone((TH1D*)measuredJetSpectrum->Clone("clone"));
        clone->SetNameTitle(Form("MeasuredJetSpectrum%s", suffix.Data()), Form("measuredJetSpectrum %s", suffix.Data()));
        return clone;//RebinTH1D(clone, fBinsTrue, clone->GetName(), kFALSE);
    }
    else return 0x0; 
    // do the actual unfolding with the selected function
    return (this->*myFunction)( measuredJetSpectrum, resizedResponse, kinematicEfficiency, measuredJetSpectrumTrueBins, suffix, jetFindingEfficiency);
} 
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldSpectrumChi2(
        const TH1D* measuredJetSpectrum,      // truncated raw jets (same binning as pt rec of response) 
        const TH2D* resizedResponse,          // response matrix
        const TH1D* kinematicEfficiency,      // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,        // unfolding template: same binning is pt gen of response
        const TString suffix,                 // suffix (in or out of plane)
        const TH1D* jetFindingEfficiency)     // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // unfold the spectrum using chi2 minimization

    // step 0) setup the static members of AliUnfolding
    ResetAliUnfolding();                // reset from previous iteration
                                        // also deletes and re-creates the global TVirtualFitter
    AliUnfolding::SetUnfoldingMethod(AliUnfolding::kChi2Minimization);
    if(!strcmp("in", suffix.Data())) AliUnfolding::SetChi2Regularization(AliUnfolding::kLogLog, fBetaIn);
    else if(!strcmp("out", suffix.Data())) AliUnfolding::SetChi2Regularization(AliUnfolding::kLogLog, fBetaOut);
    if(!strcmp("prior_in", suffix.Data())) AliUnfolding::SetChi2Regularization(AliUnfolding::kLogLog, fBetaIn);
    else if(!strcmp("prior_out", suffix.Data())) AliUnfolding::SetChi2Regularization(AliUnfolding::kLogLog, fBetaOut);
    AliUnfolding::SetNbins(fBinsRec->GetSize()-1, fBinsTrue->GetSize()-1);

    // step 1) clone all input histograms. the histograms are cloned to make sure that the original histograms
    // stay intact. a local copy of a histogram (which only exists in the scope of this function) is 
    // denoted by the suffix 'Local'
    
    // measuredJetSpectrumLocal holds the spectrum that needs to be unfolded
    TH1D *measuredJetSpectrumLocal = (TH1D*)measuredJetSpectrum->Clone(Form("measuredJetSpectrumLocal_%s", suffix.Data()));
    // unfolded local will be filled with the result of the unfolding
    TH1D *unfoldedLocal(new TH1D(Form("unfoldedLocal_%s", suffix.Data()), Form("unfoldedLocal_%s", suffix.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));

    // full response matrix and kinematic efficiency
    TH2D* resizedResponseLocal = (TH2D*)resizedResponse->Clone(Form("resizedResponseLocal_%s", suffix.Data()));
    TH1D* kinematicEfficiencyLocal = (TH1D*)kinematicEfficiency->Clone(Form("kinematicEfficiencyLocal_%s", suffix.Data()));

    // the initial guess for the unfolded pt spectrum, equal to the folded spectrum, but in 'true' bins
    TH1D *priorLocal = (TH1D*)measuredJetSpectrumTrueBins->Clone(Form("priorLocal_%s", suffix.Data()));
    // optionally, the prior can be smoothened by extrapolating the spectrum using a power law fit
    if(fSmoothenPrior) priorLocal = SmoothenPrior(priorLocal, fPower, fFitMin, fFitMax, fFitStart, kTRUE, fSmoothenCounts);

    // step 2) start the unfolding
    Int_t status(-1), i(0);
    while(status < 0 && i < 100) {
        // i > 0 means that the first iteration didn't converge. in that case, the result of the first
        // iteration (stored in unfoldedLocal) is cloned and used as a starting point for the 
        if (i > 0) priorLocal = (TH1D*)unfoldedLocal->Clone(Form("priorLocal_%s_%i", suffix.Data(), i));
        status = AliUnfolding::Unfold(
                resizedResponseLocal,           // response matrix
                kinematicEfficiencyLocal,       // efficiency applied on the unfolded spectrum (can be NULL)
                measuredJetSpectrumLocal,              // measured spectrum
                priorLocal,                     // initial conditions (set NULL to use measured spectrum)
                unfoldedLocal);                 // results
        // status holds the minuit fit status (where 0 means convergence)
        i++;
    }
    // get the status of TMinuit::mnhess(), fISW[1] == 3 means the hessian matrix was calculated succesfully
    TH2D* hPearson(0x0);
    if(status == 0 && gMinuit->fISW[1] == 3) {
        // if the unfolding converged and the hessian matrix is reliable, plot the pearson coefficients
        TVirtualFitter *fitter(TVirtualFitter::GetFitter());
        if(gMinuit) gMinuit->Command("SET COV");
        TMatrixD covarianceMatrix(fBinsTrue->GetSize()-1, fBinsTrue->GetSize()-1, fitter->GetCovarianceMatrix());
        TMatrixD *pearson((TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix));
        pearson->Print();
        hPearson = new TH2D(*pearson);
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients, %s plane", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
        if(fMergeBinsArray) unfoldedLocal = MergeSpectrumBins(fMergeBinsArray, unfoldedLocal, hPearson);
    } else status = -1; 

    // step 3) refold the unfolded spectrum and save the ratio measured / refolded
    TH1D *foldedLocal(fResponseMaker->MultiplyResponseGenerated(unfoldedLocal, resizedResponseLocal,kinematicEfficiencyLocal));
    foldedLocal->SetNameTitle(Form("RefoldedSpectrum_%s", suffix.Data()), Form("Refolded jet spectrum, %s plane", suffix.Data()));
    unfoldedLocal->SetNameTitle(Form("UnfoldedSpectrum_%s", suffix.Data()), Form("Unfolded jet spectrum, %s plane", suffix.Data()));
    TGraphErrors* ratio(GetRatio(foldedLocal, measuredJetSpectrumLocal, TString(""), kTRUE, fBinsTrue->At(fBinsTrue->GetSize()-1)));
    if(ratio) {
        ratio->SetNameTitle("RatioRefoldedMeasured", Form("Ratio measured, re-folded %s ", suffix.Data()));
        ratio->GetYaxis()->SetTitle("ratio measured / re-folded");
        ratio = ProtectHeap(ratio);
        ratio->Write();
    }

    // step 4) write histograms to file. to ensure that these have unique identifiers on the heap, 
    // objects are cloned using 'ProtectHeap()'
    measuredJetSpectrumLocal->SetNameTitle(Form("InputSpectrum_%s", suffix.Data()), Form("InputSpectrum_%s", suffix.Data()));
    measuredJetSpectrumLocal = ProtectHeap(measuredJetSpectrumLocal);
    measuredJetSpectrumLocal->Write(); 

    resizedResponseLocal = ProtectHeap(resizedResponseLocal);
    resizedResponseLocal->Write();

    unfoldedLocal = ProtectHeap(unfoldedLocal);
    if(jetFindingEfficiency) unfoldedLocal->Divide(jetFindingEfficiency);
    unfoldedLocal->Write();

    foldedLocal = ProtectHeap(foldedLocal);
    foldedLocal->Write();
    
    priorLocal = ProtectHeap(priorLocal);
    priorLocal->Write();

    // step 5) save the fit status (penalty value, degrees of freedom, chi^2 value)
    TH1F* fitStatus(new TH1F(Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), 4, -0.5, 3.5));
    fitStatus->SetBinContent(1, AliUnfolding::fChi2FromFit);
    fitStatus->GetXaxis()->SetBinLabel(1, "fChi2FromFit");
    fitStatus->SetBinContent(2, AliUnfolding::fPenaltyVal);
    fitStatus->GetXaxis()->SetBinLabel(2, "fPenaltyVal");
    fitStatus->SetBinContent(3, fBinsRec->GetSize()-fBinsTrue->GetSize());
    fitStatus->GetXaxis()->SetBinLabel(3, "DOF");
    fitStatus->SetBinContent(4, (!strcmp(suffix.Data(), "in")) ? fBetaIn : fBetaOut);
    fitStatus->GetXaxis()->SetBinLabel(4, (!strcmp(suffix.Data(), "in")) ? "fBetaIn" : "fBetaOut");
    fitStatus->Write();

    return unfoldedLocal;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldSpectrumSVD(
        const TH1D* measuredJetSpectrum,              // jet pt in pt rec bins 
        const TH2D* resizedResponse,                   // full response matrix, normalized in slides of pt true
        const TH1D* kinematicEfficiency,              // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        const TString suffix,                         // suffix (in, out)
        const TH1D* jetFindingEfficiency)             // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    TH1D* priorLocal( GetPrior(
        measuredJetSpectrum,              // jet pt in pt rec bins 
        resizedResponse,                  // full response matrix, normalized in slides of pt true
        kinematicEfficiency,              // kinematic efficiency
        measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        suffix,                           // suffix (in, out)
        jetFindingEfficiency));           // jet finding efficiency (optional)
    if(!priorLocal) {
        printf(" > couldn't find prior ! < \n");
        return 0x0;
    } else printf(" 1) retrieved prior \n");

    // go back to the 'root' directory of this instance of the SVD unfolding routine
    (!strcmp(suffix.Data(), "in")) ? fActiveDir->cd(Form("InPlane___%s", fActiveString.Data())) : fActiveDir->cd(Form("OutOfPlane___%s", fActiveString.Data()));

    // 2) setup all the necessary input for the unfolding routine. all input histograms are copied locally
    // measured jets in pt rec binning
    TH1D *measuredJetSpectrumLocal((TH1D*)measuredJetSpectrum->Clone(Form("jets_%s", suffix.Data())));
    // local copie of the response matrix
    TH2D *resizedResponseLocal((TH2D*)resizedResponse->Clone(Form("resizedResponseLocal_%s", suffix.Data())));
    // local copy of response matrix, all true slides normalized to 1 
    // this response matrix will eventually be used in the re-folding routine
    TH2D *resizedResponseLocalNorm((TH2D*)resizedResponse->Clone(Form("resizedResponseLocalNorm_%s", suffix.Data())));
    resizedResponseLocalNorm = NormalizeTH2D(resizedResponseLocalNorm);
    // kinematic efficiency
    TH1D *kinematicEfficiencyLocal((TH1D*)kinematicEfficiency->Clone(Form("kinematicEfficiency_%s", suffix.Data())));
    // place holder histos
    TH1D *unfoldedLocalSVD(0x0);
    TH1D *foldedLocalSVD(0x0);
    cout << " 2) setup necessary input " << endl;
    // 3) configure routine
    RooUnfold::ErrorTreatment errorTreatment = (fSVDToy) ? RooUnfold::kCovToy : RooUnfold::kCovariance;
    cout << " step 3) configured routine " << endl;

    // 4) get transpose matrices
    // a) get the transpose of the full response matrix
    TH2* responseMatrixLocalTransposePrior(fResponseMaker->GetTransposeResponsMatrix(resizedResponseLocal));
    responseMatrixLocalTransposePrior->SetNameTitle(Form("prior_%s_%s", responseMatrixLocalTransposePrior->GetName(), suffix.Data()),Form("prior_%s_%s", responseMatrixLocalTransposePrior->GetName(), suffix.Data()));
    // normalize it with the prior. this will ensure that high statistics bins will constrain the
    // end result most strenuously than bins with limited number of counts
    responseMatrixLocalTransposePrior = fResponseMaker->NormalizeResponsMatrixYaxisWithPrior(responseMatrixLocalTransposePrior, priorLocal);
    cout << " 4) retrieved first transpose matrix " << endl;
 
    // 5) get response for SVD unfolding
    RooUnfoldResponse responseSVD(0, 0, responseMatrixLocalTransposePrior, Form("respCombinedSVD_%s", suffix.Data()), Form("respCombinedSVD_%s", suffix.Data()));
    cout << " 5) retrieved roo unfold response object " << endl;

    // 6) actualy unfolding loop
    RooUnfoldSvd unfoldSVD(&responseSVD, measuredJetSpectrumLocal, (!strcmp(suffix.Data(), "in")) ? fSVDRegIn : fSVDRegOut);
    unfoldedLocalSVD = (TH1D*)unfoldSVD.Hreco(errorTreatment);
    // correct the spectrum for the kinematic efficiency
    unfoldedLocalSVD->Divide(kinematicEfficiencyLocal);

    // get the pearson coefficients from the covariance matrix
    TMatrixD covarianceMatrix = unfoldSVD.Ereco(errorTreatment);
    TMatrixD *pearson = (TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix);
    TH2D* hPearson(0x0);
    if(pearson) {
        hPearson = new TH2D(*pearson);
        pearson->Print();
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients_%s", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
        if(fMergeBinsArray) unfoldedLocalSVD = MergeSpectrumBins(fMergeBinsArray, unfoldedLocalSVD, hPearson);
    } else return 0x0;       // return if unfolding didn't converge

    // plot singular values and d_i vector
    TSVDUnfold* svdUnfold(unfoldSVD.Impl());
    TH1* hSVal(svdUnfold->GetSV());
    TH1D* hdi(svdUnfold->GetD());
    hSVal->SetNameTitle("SingularValuesOfAC", "Singular values of AC^{-1}");
    hSVal->SetXTitle("singular values");
    hSVal->Write();
    hdi->SetNameTitle("dVector", "d vector after orthogonal transformation");
    hdi->SetXTitle("|d_{i}^{kreg}|");
    hdi->Write();
    cout << " plotted singular values and d_i vector " << endl;

    // 7) refold the unfolded spectrum
    foldedLocalSVD = fResponseMaker->MultiplyResponseGenerated(unfoldedLocalSVD, resizedResponseLocalNorm, kinematicEfficiencyLocal);
    TGraphErrors* ratio(GetRatio(measuredJetSpectrumLocal, foldedLocalSVD, "ratio  measured / re-folded", kTRUE));
    ratio->SetNameTitle(Form("RatioRefoldedMeasured_%s", fActiveString.Data()), Form("Ratio measured / re-folded %s", fActiveString.Data()));
    ratio->GetXaxis()->SetTitle("p_{t}^{rec, rec} [GeV/ c]");
    ratio->GetYaxis()->SetTitle("ratio measured / re-folded");
    ratio->Write();
    cout << " 7) refolded the unfolded spectrum " << endl;

    // write the measured, unfolded and re-folded spectra to the output directory
    measuredJetSpectrumLocal->SetNameTitle(Form("InputSpectrum_%s", suffix.Data()), Form("input spectrum (measured) %s", suffix.Data()));
    measuredJetSpectrumLocal = ProtectHeap(measuredJetSpectrumLocal);
    measuredJetSpectrumLocal->SetXTitle("p_{t}^{rec} [GeV/c]");
    measuredJetSpectrumLocal->Write(); // input spectrum
    unfoldedLocalSVD->SetNameTitle(Form("UnfoldedSpectrum_%s",suffix.Data()), Form("unfolded spectrum %s", suffix.Data()));
    unfoldedLocalSVD = ProtectHeap(unfoldedLocalSVD);
    if(jetFindingEfficiency) unfoldedLocalSVD->Divide(jetFindingEfficiency);
    unfoldedLocalSVD->Write();  // unfolded spectrum
    foldedLocalSVD->SetNameTitle(Form("RefoldedSpectrum_%s", suffix.Data()), Form("refoldedSpectrum_%s", suffix.Data()));
    foldedLocalSVD = ProtectHeap(foldedLocalSVD);
    foldedLocalSVD->Write();    // re-folded spectrum

    // save more general bookkeeeping histograms to the output directory
    responseMatrixLocalTransposePrior->SetNameTitle("TransposeResponseMatrix", "Transpose of response matrix, normalize with prior");
    responseMatrixLocalTransposePrior->SetXTitle("p_{T, jet}^{true} [GeV/c]");
    responseMatrixLocalTransposePrior->SetYTitle("p_{T, jet}^{rec} [GeV/c]");
    responseMatrixLocalTransposePrior->Write();
    priorLocal->SetNameTitle("PriorOriginal", "Prior, original");
    priorLocal->SetXTitle("p_{t} [GeV/c]");
    priorLocal = ProtectHeap(priorLocal);
    priorLocal->Write();
    resizedResponseLocalNorm = ProtectHeap(resizedResponseLocalNorm);
    resizedResponseLocalNorm->Write();

    // save some info 
    TH1F* fitStatus(new TH1F(Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), 1, -0.5, 0.5));
    fitStatus->SetBinContent(1, (!strcmp(suffix.Data(), "in")) ? fSVDRegIn : fSVDRegOut);
    fitStatus->GetXaxis()->SetBinLabel(1, (!strcmp(suffix.Data(), "in")) ? "fSVDRegIn" : "fSVDRegOut");
    fitStatus->Write();

    return unfoldedLocalSVD;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldSpectrumBayesianAli(
        const TH1D* measuredJetSpectrum,              // jet pt in pt rec bins 
        const TH2D* resizedResponse,                  // full response matrix, normalized in slides of pt true
        const TH1D* kinematicEfficiency,              // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        const TString suffix,                         // suffix (in, out)
        const TH1D* jetFindingEfficiency)             // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // unfold the spectrum using the bayesian unfolding impelmented in AliUnfolding
    // FIXME careful, not tested yet ! (06122013) FIXME

    // step 0) setup the static members of AliUnfolding
    ResetAliUnfolding();                // reset from previous iteration
                                        // also deletes and re-creates the global TVirtualFitter
    AliUnfolding::SetUnfoldingMethod(AliUnfolding::kBayesian);
    if(!strcmp("in", suffix.Data())) AliUnfolding::SetBayesianParameters(fBayesianSmoothIn, fBayesianIterIn);
    else if(!strcmp("out", suffix.Data())) AliUnfolding::SetBayesianParameters(fBayesianSmoothOut, fBayesianIterOut);
    else if(!strcmp("prior_in", suffix.Data())) AliUnfolding::SetBayesianParameters(fBayesianSmoothIn, fBayesianIterIn);
    else if(!strcmp("prior_out", suffix.Data())) AliUnfolding::SetBayesianParameters(fBayesianSmoothOut, fBayesianIterOut);
    AliUnfolding::SetNbins(fBinsRec->GetSize()-1, fBinsTrue->GetSize()-1);

    // 1) get a prior for unfolding and clone all the input histograms
    TH1D* priorLocal( GetPrior(
        measuredJetSpectrum,              // jet pt in pt rec bins 
        resizedResponse,                  // full response matrix, normalized in slides of pt true
        kinematicEfficiency,              // kinematic efficiency
        measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        suffix,                           // suffix (in, out)
        jetFindingEfficiency));           // jet finding efficiency (optional)
    if(!priorLocal) {
        printf(" > couldn't find prior ! < \n");
        return 0x0;
    } else printf(" 1) retrieved prior \n");
    // switch back to root dir of this unfolding procedure
    (!strcmp(suffix.Data(), "in")) ? fActiveDir->cd(Form("InPlane___%s", fActiveString.Data())) : fActiveDir->cd(Form("OutOfPlane___%s", fActiveString.Data()));

    // measuredJetSpectrumLocal holds the spectrum that needs to be unfolded
    TH1D *measuredJetSpectrumLocal = (TH1D*)measuredJetSpectrum->Clone(Form("measuredJetSpectrumLocal_%s", suffix.Data()));
    // unfolded local will be filled with the result of the unfolding
    TH1D *unfoldedLocal(new TH1D(Form("unfoldedLocal_%s", suffix.Data()), Form("unfoldedLocal_%s", suffix.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));

    // full response matrix and kinematic efficiency
    TH2D* resizedResponseLocal = (TH2D*)resizedResponse->Clone(Form("resizedResponseLocal_%s", suffix.Data()));
    TH1D* kinematicEfficiencyLocal = (TH1D*)kinematicEfficiency->Clone(Form("kinematicEfficiencyLocal_%s", suffix.Data()));

    // step 2) start the unfolding
    Int_t status(-1), i(0);
    while(status < 0 && i < 100) {
        // i > 0 means that the first iteration didn't converge. in that case, the result of the first
        // iteration (stored in unfoldedLocal) is cloned and used as a starting point for the 
        if (i > 0) priorLocal = (TH1D*)unfoldedLocal->Clone(Form("priorLocal_%s_%i", suffix.Data(), i));
        status = AliUnfolding::Unfold(
                resizedResponseLocal,           // response matrix
                kinematicEfficiencyLocal,       // efficiency applied on the unfolded spectrum (can be NULL)
                measuredJetSpectrumLocal,              // measured spectrum
                priorLocal,                     // initial conditions (set NULL to use measured spectrum)
                unfoldedLocal);                 // results
        // status holds the minuit fit status (where 0 means convergence)
        i++;
    }
    // get the status of TMinuit::mnhess(), fISW[1] == 3 means the hessian matrix was calculated succesfully
    TH2D* hPearson(0x0);
    if(status == 0 && gMinuit->fISW[1] == 3) {
        // if the unfolding converged and the hessian matrix is reliable, plot the pearson coefficients
        TVirtualFitter *fitter(TVirtualFitter::GetFitter());
        if(gMinuit) gMinuit->Command("SET COV");
        TMatrixD covarianceMatrix(fBinsTrue->GetSize()-1, fBinsTrue->GetSize()-1, fitter->GetCovarianceMatrix());
        TMatrixD *pearson((TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix));
        pearson->Print();
        hPearson= new TH2D(*pearson);
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients, %s plane", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
    } else status = -1; 
    if(fMergeBinsArray) unfoldedLocal = MergeSpectrumBins(fMergeBinsArray, unfoldedLocal, hPearson);

    // step 3) refold the unfolded spectrum and save the ratio measured / refolded
    TH1D *foldedLocal(fResponseMaker->MultiplyResponseGenerated(unfoldedLocal, resizedResponseLocal,kinematicEfficiencyLocal));
    foldedLocal->SetNameTitle(Form("RefoldedSpectrum_%s", suffix.Data()), Form("Refolded jet spectrum, %s plane", suffix.Data()));
    unfoldedLocal->SetNameTitle(Form("UnfoldedSpectrum_%s", suffix.Data()), Form("Unfolded jet spectrum, %s plane", suffix.Data()));
    TGraphErrors* ratio(GetRatio(foldedLocal, measuredJetSpectrumLocal, TString(""), kTRUE, fBinsTrue->At(fBinsTrue->GetSize()-1)));
    if(ratio) {
        ratio->SetNameTitle("RatioRefoldedMeasured", Form("Ratio measured, re-folded %s ", suffix.Data()));
        ratio->GetYaxis()->SetTitle("ratio measured / re-folded");
        ratio = ProtectHeap(ratio);
        ratio->Write();
    }

    // step 4) write histograms to file. to ensure that these have unique identifiers on the heap, 
    // objects are cloned using 'ProtectHeap()'
    measuredJetSpectrumLocal->SetNameTitle(Form("InputSpectrum_%s", suffix.Data()), Form("InputSpectrum_%s", suffix.Data()));
    measuredJetSpectrumLocal = ProtectHeap(measuredJetSpectrumLocal);
    measuredJetSpectrumLocal->Write(); 

    resizedResponseLocal = ProtectHeap(resizedResponseLocal);
    resizedResponseLocal->Write();

    unfoldedLocal = ProtectHeap(unfoldedLocal);
    if(jetFindingEfficiency) unfoldedLocal->Divide(jetFindingEfficiency);
    unfoldedLocal->Write();

    foldedLocal = ProtectHeap(foldedLocal);
    foldedLocal->Write();
    
    priorLocal = ProtectHeap(priorLocal);
    priorLocal->Write();

    // step 5) save the fit status (penalty value, degrees of freedom, chi^2 value)
    TH1F* fitStatus(new TH1F(Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), 4, -0.5, 3.5));
    fitStatus->SetBinContent(1, AliUnfolding::fChi2FromFit);
    fitStatus->GetXaxis()->SetBinLabel(1, "fChi2FromFit");
    fitStatus->SetBinContent(2, AliUnfolding::fPenaltyVal);
    fitStatus->GetXaxis()->SetBinLabel(2, "fPenaltyVal");
    fitStatus->SetBinContent(3, fBinsRec->GetSize()-fBinsTrue->GetSize());
    fitStatus->GetXaxis()->SetBinLabel(3, "DOF");
    fitStatus->SetBinContent(4, (!strcmp(suffix.Data(), "in")) ? fBetaIn : fBetaOut);
    fitStatus->GetXaxis()->SetBinLabel(4, (!strcmp(suffix.Data(), "in")) ? "fBetaIn" : "fBetaOut");
    fitStatus->Write();

    return unfoldedLocal;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldSpectrumBayesian(
        const TH1D* measuredJetSpectrum,              // jet pt in pt rec bins 
        const TH2D* resizedResponse,                  // full response matrix, normalized in slides of pt true
        const TH1D* kinematicEfficiency,              // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        const TString suffix,                         // suffix (in, out)
        const TH1D* jetFindingEfficiency)             // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // use bayesian unfolding from the RooUnfold package to unfold jet spectra
    
    // 1) get a prior for unfolding.
    TH1D* priorLocal( GetPrior(
        measuredJetSpectrum,              // jet pt in pt rec bins 
        resizedResponse,                  // full response matrix, normalized in slides of pt true
        kinematicEfficiency,              // kinematic efficiency
        measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        suffix,                           // suffix (in, out)
        jetFindingEfficiency));           // jet finding efficiency (optional)
    if(!priorLocal) {
        printf(" > couldn't find prior ! < \n");
        return 0x0;
    } else printf(" 1) retrieved prior \n");
    (!strcmp(suffix.Data(), "in")) ? fActiveDir->cd(Form("InPlane___%s", fActiveString.Data())) : fActiveDir->cd(Form("OutOfPlane___%s", fActiveString.Data()));

    // 2) setup all the necessary input for the unfolding routine. all input histograms are copied locally
    // measured jets in pt rec binning
    TH1D *measuredJetSpectrumLocal((TH1D*)measuredJetSpectrum->Clone(Form("jets_%s", suffix.Data())));
    // local copie of the response matrix
    TH2D *resizedResponseLocal((TH2D*)resizedResponse->Clone(Form("resizedResponseLocal_%s", suffix.Data())));
    // local copy of response matrix, all true slides normalized to 1 
    // this response matrix will eventually be used in the re-folding routine
    TH2D *resizedResponseLocalNorm((TH2D*)resizedResponse->Clone(Form("resizedResponseLocalNorm_%s", suffix.Data())));
    resizedResponseLocalNorm = NormalizeTH2D(resizedResponseLocalNorm);
    // kinematic efficiency
    TH1D *kinematicEfficiencyLocal((TH1D*)kinematicEfficiency->Clone(Form("kinematicEfficiency_%s", suffix.Data())));
    // place holder histos
    TH1D *unfoldedLocalBayes(0x0);
    TH1D *foldedLocalBayes(0x0);
    cout << " 2) setup necessary input " << endl;
    // 4) get transpose matrices
    // a) get the transpose of the full response matrix
    TH2* responseMatrixLocalTransposePrior(fResponseMaker->GetTransposeResponsMatrix(resizedResponseLocal));
    responseMatrixLocalTransposePrior->SetNameTitle(Form("prior_%s_%s", responseMatrixLocalTransposePrior->GetName(), suffix.Data()),Form("prior_%s_%s", responseMatrixLocalTransposePrior->GetName(), suffix.Data()));
    // normalize it with the prior. this will ensure that high statistics bins will constrain the
    // end result most strenuously than bins with limited number of counts
    responseMatrixLocalTransposePrior = fResponseMaker->NormalizeResponsMatrixYaxisWithPrior(responseMatrixLocalTransposePrior, priorLocal);
    // 3) get response for Bayesian unfolding
    RooUnfoldResponse responseBayes(0, 0, responseMatrixLocalTransposePrior, Form("respCombinedBayes_%s", suffix.Data()), Form("respCombinedBayes_%s", suffix.Data()));

    // 4) actualy unfolding loop
    RooUnfoldBayes unfoldBayes(&responseBayes, measuredJetSpectrumLocal, (!strcmp("in", suffix.Data())) ? fBayesianIterIn : fBayesianIterOut);
    RooUnfold::ErrorTreatment errorTreatment = (fSVDToy) ? RooUnfold::kCovToy : RooUnfold::kCovariance;
    unfoldedLocalBayes = (TH1D*)unfoldBayes.Hreco(errorTreatment);
    // correct the spectrum for the kinematic efficiency
    unfoldedLocalBayes->Divide(kinematicEfficiencyLocal);
    // get the pearson coefficients from the covariance matrix
    TMatrixD covarianceMatrix = unfoldBayes.Ereco(errorTreatment);
    TMatrixD *pearson = (TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix);
    TH2D* hPearson(0x0);
    if(pearson) {
        hPearson = new TH2D(*pearson);
        pearson->Print();
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients_%s", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
        if(fMergeBinsArray) unfoldedLocalBayes = MergeSpectrumBins(fMergeBinsArray, unfoldedLocalBayes, hPearson);
    } else return 0x0;       // return if unfolding didn't converge

    // 5) refold the unfolded spectrum
    foldedLocalBayes = fResponseMaker->MultiplyResponseGenerated(unfoldedLocalBayes, resizedResponseLocalNorm, kinematicEfficiencyLocal);
    TGraphErrors* ratio(GetRatio(measuredJetSpectrumLocal, foldedLocalBayes, "ratio  measured / re-folded", kTRUE));
    ratio->SetNameTitle(Form("RatioRefoldedMeasured_%s", fActiveString.Data()), Form("Ratio measured / re-folded %s", fActiveString.Data()));
    ratio->GetXaxis()->SetTitle("p_{t}^{rec, rec} [GeV/ c]");
    ratio->GetYaxis()->SetTitle("ratio measured / re-folded");
    ratio->Write();
    cout << " 7) refolded the unfolded spectrum " << endl;

    // write the measured, unfolded and re-folded spectra to the output directory
    measuredJetSpectrumLocal->SetNameTitle(Form("InputSpectrum_%s", suffix.Data()), Form("input spectrum (measured) %s", suffix.Data()));
    measuredJetSpectrumLocal = ProtectHeap(measuredJetSpectrumLocal);
    measuredJetSpectrumLocal->SetXTitle("p_{t}^{rec} [GeV/c]");
    measuredJetSpectrumLocal->Write(); // input spectrum
    unfoldedLocalBayes->SetNameTitle(Form("UnfoldedSpectrum_%s",suffix.Data()), Form("unfolded spectrum %s", suffix.Data()));
    unfoldedLocalBayes = ProtectHeap(unfoldedLocalBayes);
    if(jetFindingEfficiency) unfoldedLocalBayes->Divide(jetFindingEfficiency);
    unfoldedLocalBayes->Write();  // unfolded spectrum
    foldedLocalBayes->SetNameTitle(Form("RefoldedSpectrum_%s", suffix.Data()), Form("refoldedSpectrum_%s", suffix.Data()));
    foldedLocalBayes = ProtectHeap(foldedLocalBayes);
    foldedLocalBayes->Write();    // re-folded spectrum

    // save more general bookkeeeping histograms to the output directory
    responseMatrixLocalTransposePrior->SetNameTitle("TransposeResponseMatrix", "Transpose of response matrix, normalize with prior");
    responseMatrixLocalTransposePrior->SetXTitle("p_{T, jet}^{true} [GeV/c]");
    responseMatrixLocalTransposePrior->SetYTitle("p_{T, jet}^{rec} [GeV/c]");
    responseMatrixLocalTransposePrior->Write();
    priorLocal->SetNameTitle("PriorOriginal", "Prior, original");
    priorLocal->SetXTitle("p_{t} [GeV/c]");
    priorLocal = ProtectHeap(priorLocal);
    priorLocal->Write();
    resizedResponseLocalNorm = ProtectHeap(resizedResponseLocalNorm);
    resizedResponseLocalNorm->Write();

    // save some info 
    TH1F* fitStatus(new TH1F(Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), Form("fitStatus_%s_%s", fActiveString.Data(), suffix.Data()), 1, -0.5, 0.5));
    fitStatus->SetBinContent(1, (!strcmp(suffix.Data(), "in")) ? fBayesianIterIn : fBayesianIterOut);
    fitStatus->GetXaxis()->SetBinLabel(1, (!strcmp(suffix.Data(), "in")) ? "fBayesianIterIn" : "fBayesianIterOut");
    fitStatus->Write();

    return  unfoldedLocalBayes; 
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::FoldSpectrum(
        const TH1D* measuredJetSpectrum,      // truncated raw jets (same binning as pt rec of response) 
        const TH2D* resizedResponse,          // response matrix
        const TH1D* kinematicEfficiency,      // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,        // not used
        const TString suffix,                 // suffix (in or out of plane)
        const TH1D* jetFindingEfficiency)     // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // simple function to fold the given spectrum with the in-plane and out-of-plane response

    if(!measuredJetSpectrumTrueBins) SquelchWarning();
    // 0) for consistency with the other methods, keep the same nomenclature which admittedly is a bit confusing 
    // what is 'unfolded' here, is just a clone of the input spectrum, binned to the 'unfolded' binning
    TH1D* unfoldedLocal((TH1D*)measuredJetSpectrum->Clone(Form("unfoldedLocal_%s", suffix.Data())));

    // 1) full response matrix and kinematic efficiency
    TH2D* resizedResponseLocal = (TH2D*)resizedResponse->Clone(Form("resizedResponseLocal_%s", suffix.Data()));
    TH1D* kinematicEfficiencyLocal = (TH1D*)kinematicEfficiency->Clone(Form("kinematicEfficiencyLocal_%s", suffix.Data()));

    // step 2) fold the 'unfolded' spectrum and save the ratio measured / refolded
    TH1D *foldedLocal(fResponseMaker->MultiplyResponseGenerated(unfoldedLocal, resizedResponseLocal,kinematicEfficiencyLocal));
    foldedLocal->SetNameTitle(Form("RefoldedSpectrum_%s", suffix.Data()), Form("Refolded jet spectrum, %s plane", suffix.Data()));
    unfoldedLocal->SetNameTitle(Form("UnfoldedSpectrum_%s", suffix.Data()), Form("Unfolded jet spectrum, %s plane", suffix.Data()));
 
    // step 3) write histograms to file. to ensure that these have unique identifiers on the heap, 
    // objects are cloned using 'ProtectHeap()'
    TH1D* measuredJetSpectrumLocal((TH1D*)(measuredJetSpectrum->Clone("tempObject")));
    measuredJetSpectrumLocal->SetNameTitle(Form("InputSpectrum_%s", suffix.Data()), Form("InputSpectrum_%s", suffix.Data()));
    measuredJetSpectrumLocal = ProtectHeap(measuredJetSpectrumLocal);
    measuredJetSpectrumLocal->Write(); 

    resizedResponseLocal = ProtectHeap(resizedResponseLocal);
    resizedResponseLocal->Write();

    unfoldedLocal = ProtectHeap(unfoldedLocal);
    if(jetFindingEfficiency) unfoldedLocal->Divide(jetFindingEfficiency);
    unfoldedLocal->Write();

    foldedLocal = ProtectHeap(foldedLocal);
    foldedLocal->Write();

    // return the folded result
    return foldedLocal;
}
//_____________________________________________________________________________
Bool_t AliJetFlowTools::PrepareForUnfolding(TH1* customIn, TH1* customOut)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // prepare for unfolding. check if all the necessary input data is available, if not, retrieve it
    if(fRawInputProvided) return kTRUE;
    if(!fInputList) {
        printf(" AliJetFlowTools::PrepareForUnfolding() fInputList not found \n - Set a list using AliJetFlowTools::SetInputList() \n");
        return kFALSE;
    }
    if(!fDetectorResponse) printf(" WARNING, no detector response supplied ! May be ok (depending on what you want to do) \n ");
    // check if the pt bin for true and rec have been set
    if(!fBinsTrue || !fBinsRec) {
        printf(" AliJetFlowTools::PrepareForUnfolding() no true or rec bins set, aborting ! \n");
        return kFALSE;
    }
    if(!fRMSSpectrumIn && fDphiUnfolding) { // initialie the profiles which will hold the RMS values. if binning changes in between unfolding
                          // procedures, these profiles will be nonsensical, user is responsible
        fRMSSpectrumIn = new TProfile("fRMSSpectrumIn", "fRMSSpectrumIn", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
        fRMSSpectrumOut = new TProfile("fRMSSpectrumOut", "fRMSSpectrumOut", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
        fRMSRatio = new TProfile("fRMSRatio", "fRMSRatio", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
    }
    if(!fTrainPower) {
        // clear minuit state to avoid constraining the fit with the results of the previous iteration
        for(Int_t i(0); i < fPower->GetNpar(); i++) fPower->SetParameter(i, 0.);
    }
    // extract the spectra 
    TString spectrumName(Form("fHistJetPsi3Pt_%i", fCentralityArray->At(0)));
    if(fRho0) spectrumName = Form("fHistJetPsi3PtRho0_%i", fCentralityArray->At(0));
    if(!fInputList->FindObject(spectrumName.Data())) {
        printf(" Couldn't find spectrum %s ! \n", spectrumName.Data());
        return kFALSE;
    }

    // get the first scaled spectrum
    fJetPtDeltaPhi = (TH2D*)fInputList->FindObject(spectrumName.Data());
    // clone the spectrum on the heap. this is necessary since scale or add change the
    // contents of the original histogram
    fJetPtDeltaPhi = ProtectHeap(fJetPtDeltaPhi, kFALSE);
    fJetPtDeltaPhi->Scale(fCentralityWeights->At(0));
    printf("Extracted %s wight weight %.2f \n", spectrumName.Data(), fCentralityWeights->At(0));
    // merge subsequent bins (if any)
    for(Int_t i(1); i < fCentralityArray->GetSize(); i++) {
        spectrumName = Form("fHistJetPsi3Pt_%i", fCentralityArray->At(i));
        printf( " Merging with %s with weight %.4f \n", spectrumName.Data(), fCentralityWeights->At(i));
        fJetPtDeltaPhi->Add(((TH2D*)fInputList->FindObject(spectrumName.Data())), fCentralityWeights->At(i));
    }
    // if a second list is passed, merge this with the existing one
    if(fMergeWithList) {
        spectrumName = Form("fHistJetPsi3Pt_%i", fMergeWithCen);
        printf( " Adding additional output %s \n", spectrumName.Data());
        fJetPtDeltaPhi->Add(((TH2D*)fMergeWithList->FindObject(spectrumName.Data())), fMergeWithWeight);
    }

    // in plane spectrum
    if(!fDphiUnfolding) {
        fSpectrumIn = fJetPtDeltaPhi->ProjectionY(Form("_py_in_%s", spectrumName.Data()), 1, 40, "e");
        fSpectrumOut = fJetPtDeltaPhi->ProjectionY(Form("_py_out_%s", spectrumName.Data()), 1, 40, "e");
    } else {
        fSpectrumIn = fJetPtDeltaPhi->ProjectionY(Form("_py_ina_%s", spectrumName.Data()), 1, 10, "e");
        fSpectrumIn->Add(fJetPtDeltaPhi->ProjectionY(Form("_py_inb_%s", spectrumName.Data()), 31, 40, "e"));
        fSpectrumIn = ProtectHeap(fSpectrumIn);
        // out of plane spectrum
        fSpectrumOut = fJetPtDeltaPhi->ProjectionY(Form("_py_out_%s", spectrumName.Data()), 11, 30, "e");
        fSpectrumOut = ProtectHeap(fSpectrumOut);
    }
    // if a custom input is passed, overwrite existing histograms
    if(customIn)        fSpectrumIn = dynamic_cast<TH1D*>(customIn);
    if(customOut)       fSpectrumOut = dynamic_cast<TH1D*>(customOut);

    // normalize spectra to event count if requested
    if(fNormalizeSpectra) {
        TH1* rho((TH1*)fInputList->FindObject(Form("fHistRho_%i", fCentralityArray->At(0))));
        rho->Scale(fCentralityWeights->At(0));
        for(Int_t i(1); i < fCentralityArray->GetSize(); i++) {
           rho->Add((TH1*)fInputList->FindObject(Form("fHistRho_%i", fCentralityArray->At(i))), fCentralityWeights->At(i));
        }
        if(!rho) return 0x0;
        Bool_t normalizeToFullSpectrum = (fEventCount < 0) ? kTRUE : kFALSE;
        if (normalizeToFullSpectrum) fEventCount = rho->GetEntries();
        if(fEventCount > 0) {
            fSpectrumIn->Sumw2();       // necessary for correct error propagation of scale
            fSpectrumOut->Sumw2();
            Double_t pt(0);            
            Double_t error(0); // lots of issues with the errors here ...
            for(Int_t i(0); i < fSpectrumIn->GetXaxis()->GetNbins(); i++) {
                pt = fSpectrumIn->GetBinContent(1+i)/fEventCount;       // normalized count
                error = 1./((double)(fEventCount*fEventCount))*fSpectrumIn->GetBinError(1+i)*fSpectrumIn->GetBinError(1+i);
                fSpectrumIn->SetBinContent(1+i, pt);
                if(pt <= 0 ) fSpectrumIn->SetBinError(1+i, 0.);
                if(error > 0) fSpectrumIn->SetBinError(1+i, error);
                else fSpectrumIn->SetBinError(1+i, TMath::Sqrt(pt));
            }
            for(Int_t i(0); i < fSpectrumOut->GetXaxis()->GetNbins(); i++) {
                pt = fSpectrumOut->GetBinContent(1+i)/fEventCount;       // normalized count
                error = 1./((double)(fEventCount*fEventCount))*fSpectrumOut->GetBinError(1+i)*fSpectrumOut->GetBinError(1+i);
                fSpectrumOut->SetBinContent(1+i, pt);
                if( pt <= 0) fSpectrumOut->SetBinError(1+i, 0.);
                if(error > 0) fSpectrumOut->SetBinError(1+i, error);
                else fSpectrumOut->SetBinError(1+i, TMath::Sqrt(pt));
            }
        }
        if(normalizeToFullSpectrum) fEventCount = -1;
    }
    // extract the delta pt matrices
    TString deltaptName("");
    if(!fRho0) {
        deltaptName += (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJ_%i", fCentralityArray->At(0)) : Form("fHistDeltaPtDeltaPhi3_%i", fCentralityArray->At(0));
    } else {
        deltaptName += (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJRho0_%i", fCentralityArray->At(0)) : Form("fHistDeltaPtDeltaPhi3Rho0_%i", fCentralityArray->At(0));
    }
    fDeltaPtDeltaPhi = ((TH2D*)fInputList->FindObject(deltaptName.Data()));
    if(!fDeltaPtDeltaPhi) {
        printf(" Couldn't find delta pt matrix %s ! \n", deltaptName.Data());
        printf(" > may be ok, depending no what you want to do < \n");
        fRefreshInput = kTRUE;
        return kTRUE;
    }

    // clone the distribution on the heap and if requested merge with other centralities
    fDeltaPtDeltaPhi = ProtectHeap(fDeltaPtDeltaPhi, kFALSE);
    fDeltaPtDeltaPhi->Scale(fCentralityWeights->At(0));
    printf("Extracted %s with weight %.2f \n", deltaptName.Data(), fCentralityWeights->At(0));
    for(Int_t i(1); i < fCentralityArray->GetSize(); i++) {
        deltaptName = (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJ_%i", fCentralityArray->At(i)) : Form("fHistDeltaPtDeltaPhi3_%i", fCentralityArray->At(i));
        printf(" Merging with %s with weight %.4f \n", deltaptName.Data(), fCentralityWeights->At(i));
        fDeltaPtDeltaPhi->Add((TH2D*)fInputList->FindObject(deltaptName.Data()), fCentralityWeights->At(i));
    }
    // if a second list is passed, merge this with the existing one
    if(fMergeWithList) {
        deltaptName = (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJ_%i", fMergeWithCen) : Form("fHistDeltaPtDeltaPhi3_%i", fMergeWithCen);
        printf(" Adding additional data %s \n", deltaptName.Data());
        fDeltaPtDeltaPhi->Add((TH2D*)fMergeWithList->FindObject(deltaptName.Data()), fMergeWithWeight);
    }

    // in plane delta pt distribution
    if(!fDphiUnfolding) {
        fDptInDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_in_%s", deltaptName.Data()), 1, 40, "e");
        fDptOutDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_out_%s", deltaptName.Data()), 1, 40, "e");
    } else {
        fDptInDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_ina_%s", deltaptName.Data()), 1, 10, "e");
        fDptInDist->Add(fDeltaPtDeltaPhi->ProjectionY(Form("_py_inb_%s", deltaptName.Data()), 31, 40, "e"));
        // out of plane delta pt distribution
        fDptOutDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_out_%s", deltaptName.Data()), 11, 30, "e");
        fDptInDist = ProtectHeap(fDptInDist);
        fDptOutDist = ProtectHeap(fDptOutDist);
        // TODO get dpt response matrix from ConstructDPtResponseFromTH1D
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
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)), Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)));
    fDptIn->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T, jet}^{rec} [GeV/c]");
    fDptIn = ProtectHeap(fDptIn);
    fDptOut = new TH2D(*rfOut);
    fDptOut->SetNameTitle(Form("dpt_response_OUTOFPLANE_%i", fCentralityArray->At(0)), Form("dpt_response_OUTOFPLANE_%i", fCentralityArray->At(0)));
    fDptOut->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV/c]");
    fDptOut->GetYaxis()->SetTitle("p_{T, jet}^{rec} [GeV/c]");
    fDptOut = ProtectHeap(fDptOut);
    
    fRefreshInput = kTRUE;     // force cloning of the input
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliJetFlowTools::PrepareForUnfolding(Int_t low, Int_t up) {
    // prepare for unfoldingUA - more robust method to extract input spectra from file
    // will replace PrepareForUnfolding eventually (09012014)
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    if(!fInputList) {
        printf(" AliJetFlowTools::PrepareForUnfolding() fInputList not found \n - Set a list using AliJetFlowTools::SetInputList() \n");
        return kFALSE;
    }
    if(!fDetectorResponse) printf(" WARNING, no detector response supplied ! May be ok (depending on what you want to do) \n ");
    // check if the pt bin for true and rec have been set
    if(!fBinsTrue || !fBinsRec) {
        printf(" AliJetFlowTools::PrepareForUnfolding() no true or rec bins set, aborting ! \n");
        return kFALSE;
    }
    if(!fTrainPower) {
        // clear minuit state to avoid constraining the fit with the results of the previous iteration
        for(Int_t i(0); i < fPower->GetNpar(); i++) fPower->SetParameter(i, 0.);
    }
    // extract the spectra 
    TString spectrumName(Form("fHistJetPsi3Pt_%i", fCentralityArray->At(0)));
    if(fRho0) spectrumName = Form("fHistJetPsi3PtRho0_%i", fCentralityArray->At(0));
    fJetPtDeltaPhi = ((TH2D*)fInputList->FindObject(spectrumName.Data()));
    if(!fJetPtDeltaPhi) {
        printf(" Couldn't find spectrum %s ! \n", spectrumName.Data());
        return kFALSE;
    }
    if(fCentralityArray) {
        for(Int_t i(1); i < fCentralityArray->GetSize(); i++) {
            spectrumName = Form("fHistJetPsi3Pt_%i", fCentralityArray->At(i));
            fJetPtDeltaPhi->Add(((TH2D*)fInputList->FindObject(spectrumName.Data())));
        }
    }
    fJetPtDeltaPhi = ProtectHeap(fJetPtDeltaPhi, kFALSE);
    // in plane spectrum
    fSpectrumIn = fJetPtDeltaPhi->ProjectionY(Form("_py_in_%s", spectrumName.Data()), low, up, "e");
    // extract the delta pt matrices
    TString deltaptName("");
    if(!fRho0) {
        deltaptName += (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJ_%i", fCentralityArray->At(0)) : Form("fHistDeltaPtDeltaPhi3_%i", fCentralityArray->At(0));
    } else {
        deltaptName += (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJRho0_%i", fCentralityArray->At(0)) : Form("fHistDeltaPtDeltaPhi3Rho0_%i", fCentralityArray->At(0));
    }
    fDeltaPtDeltaPhi = ((TH2D*)fInputList->FindObject(deltaptName.Data()));
    if(!fDeltaPtDeltaPhi) {
        printf(" Couldn't find delta pt matrix %s ! \n", deltaptName.Data());
        printf(" > may be ok, depending no what you want to do < \n");
        fRefreshInput = kTRUE;
        return kTRUE;
    }
    if(fCentralityArray) {
        for(Int_t i(1); i < fCentralityArray->GetSize(); i++) {
            deltaptName += (fExLJDpt) ? Form("fHistDeltaPtDeltaPhi3ExLJ_%i", fCentralityArray->At(i)) : Form("fHistDeltaPtDeltaPhi3_%i", fCentralityArray->At(i));
            fDeltaPtDeltaPhi->Add(((TH2D*)fInputList->FindObject(deltaptName.Data())));
        }
    }

    fDeltaPtDeltaPhi = ProtectHeap(fDeltaPtDeltaPhi, kFALSE);
    // in plane delta pt distribution
    fDptInDist = fDeltaPtDeltaPhi->ProjectionY(Form("_py_in_%s", deltaptName.Data()), low, up, "e");
    // create a rec - true smeared response matrix
    TMatrixD* rfIn = new TMatrixD(-50, 249, -50, 249);
    for(Int_t j(-50); j < 250; j++) {   // loop on pt true slices j
        Bool_t skip = kFALSE;
        for(Int_t k(-50); k < 250; k++) {       // loop on pt gen slices k
            (*rfIn)(k, j) = (skip) ? 0. : fDptInDist->GetBinContent(fDptInDist->GetXaxis()->FindBin(k-j));
            if(fAvoidRoundingError && k > j && TMath::AreEqualAbs(fDptInDist->GetBinContent(fDptInDist->GetXaxis()->FindBin(k-j)), 0, 1e-8)) skip = kTRUE;
        }
    }
    fDptIn = new TH2D(*rfIn);
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)), Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)));
    fDptIn->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T, jet}^{rec} [GeV/c]");
    fDptIn = ProtectHeap(fDptIn);
    
    return kTRUE;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::GetPrior(
        const TH1D* measuredJetSpectrum,              // jet pt in pt rec bins 
        const TH2D* resizedResponse,                  // full response matrix, normalized in slides of pt true
        const TH1D* kinematicEfficiency,              // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        const TString suffix,                         // suffix (in, out)
        const TH1D* jetFindingEfficiency)             // jet finding efficiency (optional)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // 1) get a prior for unfolding. 
    // this can be either an unfolded spectrum from e.g. chi2 unfolding or the measured spectrum
    TH1D* unfolded(0x0);
    TDirectoryFile* dirOut = new TDirectoryFile(Form("Prior_%s___%s", suffix.Data(), fActiveString.Data()), Form("Prior_%s___%s", suffix.Data(), fActiveString.Data()));
    dirOut->cd();
    switch (fPrior) {    // select the prior for unfolding
        case kPriorTF1 : {
            if(!fPriorUser) {
                printf("> GetPrior:: FATAL ERROR! TF1 prior requested but prior has not been set ! < \n");
                return 0x0;
            }
            return fPriorUser;
        } break;
        case kPriorPythia : {
            if(!fPriorUser) {
                printf("> GetPrior:: FATAL ERROR! pythia prior requested but prior has not been set ! < \n");
                return 0x0;
            }
            // rebin the given prior to the true spectrum (creates a new histo)
            return RebinTH1D(fPriorUser, fBinsTrue, Form("kPriorPythia_%s", suffix.Data()), kFALSE);
        } break;
        case kPriorChi2 : {
            TArrayI* placeHolder(0x0);
            if(fMergeBinsArray) {
                placeHolder = fMergeBinsArray;
                fMergeBinsArray = 0x0;
            }
            if(fBinsTruePrior && fBinsRecPrior) {       // if set, use different binning for the prior
                TArrayD* tempArrayTrue(fBinsTrue);      // temporarily cache the original binning
                fBinsTrue = fBinsTruePrior;             // switch binning schemes (will be used in UnfoldSpectrumChi2())
                TArrayD* tempArrayRec(fBinsRec);
                fBinsRec = fBinsRecPrior;
                // for the prior, do not re-bin the output
                TH1D* measuredJetSpectrumChi2 = RebinTH1D((!strcmp("in", suffix.Data())) ? fSpectrumIn : fSpectrumOut, fBinsRec, TString("resized_chi2"), kFALSE);
                TH1D* measuredJetSpectrumTrueBinsChi2 = RebinTH1D((!strcmp("in", suffix.Data())) ? fSpectrumIn : fSpectrumOut, fBinsTruePrior, TString("out"), kFALSE);
                TH2D* resizedResponseChi2(RebinTH2D((!strcmp("in", suffix.Data())) ? fFullResponseIn : fFullResponseOut,fBinsTruePrior, fBinsRec, TString("chi2")));
                TH1D* kinematicEfficiencyChi2(resizedResponseChi2->ProjectionX());
                kinematicEfficiencyChi2->SetNameTitle("kin_eff_chi2","kin_eff_chi2");
                for(Int_t i(0); i < kinematicEfficiencyChi2->GetXaxis()->GetNbins(); i++) kinematicEfficiencyChi2->SetBinError(1+i, 0.);
                unfolded= UnfoldSpectrumChi2(
                            measuredJetSpectrumChi2,
                            resizedResponseChi2,
                            kinematicEfficiencyChi2,
                            measuredJetSpectrumTrueBinsChi2,    // prior for chi2 unfolding (measured)
                            TString(Form("prior_%s", suffix.Data())));
               if(!unfolded) {
                    printf(" > GetPrior:: panic, couldn't get prior from Chi2 unfolding! \n");
                    printf("               probably Chi2 unfolding did not converge < \n");
                    if(placeHolder) fMergeBinsArray = placeHolder;
                    return 0x0;
                }
                fBinsTrue = tempArrayTrue;  // reset bins borders
                fBinsRec = tempArrayRec;
                // if the chi2 prior has a different binning, rebin to the true binning for the  unfolding
                unfolded = RebinTH1D(unfolded, fBinsTrue, TString(Form("unfoldedChi2Prior_%s", suffix.Data())));     // rebin unfolded
            } else {
                unfolded = UnfoldSpectrumChi2(
                            measuredJetSpectrum,
                            resizedResponse,
                            kinematicEfficiency,
                            measuredJetSpectrumTrueBins,        // prior for chi2 unfolding (measured)
                            TString(Form("prior_%s", suffix.Data())));
                if(!unfolded) {
                    printf(" > GetPrior:: panic, couldn't get prior from Chi2 unfolding! \n");
                    printf("               probably Chi2 unfolding did not converge < \n");
                    if(placeHolder) fMergeBinsArray = placeHolder;
                    return 0x0;
                }
                if(placeHolder) fMergeBinsArray = placeHolder;
            }
            break;

        }
        case kPriorMeasured : { 
            unfolded = (TH1D*)measuredJetSpectrumTrueBins->Clone(Form("kPriorMeasured_%s", suffix.Data()));       // copy template to unfolded to use as prior
        }
        default : break;
    }
    // it can be important that the prior is smooth, this can be achieved by 
    // extrapolating the spectrum with a fitted power law when bins are sparsely filed 
    if(jetFindingEfficiency) unfolded->Divide(jetFindingEfficiency);
    TH1D *priorLocal((TH1D*)unfolded->Clone(Form("priorUnfolded_%s", suffix.Data())));
    if(fSmoothenPrior) priorLocal = SmoothenPrior(priorLocal, fPower, fFitMin, fFitMax, fFitStart, kTRUE, fSmoothenCounts);
    return priorLocal;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix) {
    // resize the x-axis of a th1d
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
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
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
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
TH2D* AliJetFlowTools::NormalizeTH2D(TH2D* histo, Bool_t noError) {
    // general method to normalize all vertical slices of a th2 to unity
    // i.e. get a probability matrix
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    if(!histo) {
        printf(" > NormalizeTH2D:: NULL pointer passed, returning NULL < \n");
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
            if(noError) histo->SetBinError(  1+i, j+1, 0.);
            else histo->SetBinError(  1+i, j+1, histo->GetBinError(  1+i, j+1)/weight);
        }
    }
    return histo;
}
//_____________________________________________________________________________
TH1* AliJetFlowTools::Bootstrap(TH1* hist, Bool_t kill) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // resample a TH1
    // the returned histogram is new, the original is deleted from the heap if kill is true
    if(!hist) {
        printf(" > Bootstrap:: fatal error,NULL pointer passed! \n");
        return 0x0;
    }
    else printf(" > Bootstrap:: resampling, may take some time \n");
    // clone input histo
    TH1* bootstrapped((TH1*)(hist->Clone(Form("%s_bootstrapped", hist->GetName()))));
    bootstrapped->Reset();

    /* OLD method - slightly underestimates fluctuations
    // reset the content
    bootstrapped->Reset();
    // resample the input histogram 
    for(Int_t i(0); i < hist->GetEntries(); i++) bootstrapped->Fill(hist->GetRandom()); */
    
    // new method
    Double_t mean(0), sigma(0);
    Int_t sampledMean(0), entries(0);

    for(Int_t i(0); i < hist->GetXaxis()->GetNbins(); i++) {
        // for each bin, get the value
        mean = hist->GetBinContent(i+1);
        sigma = hist->GetBinError(i+1);
        // draw a new mean 
        sampledMean = TMath::Nint(gRandom->Gaus(mean, sigma));
        printf(" sampled %i from original number %.2f \n", sampledMean, mean);
        // set the new bin content
        bootstrapped->SetBinContent(i+1, sampledMean);
        if(sampledMean > 0) bootstrapped->SetBinError(i+1, TMath::Sqrt(sampledMean));
        entries += sampledMean;
    }
    printf(" Done bootstrapping, setting number of entries to %i \n", entries);
    bootstrapped->SetEntries((double)entries);


    // if requested kill input histo
    if(kill) delete hist;
    // return resampled histogram
    return bootstrapped;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::RebinTH1D(TH1D* histo, TArrayD* bins, TString suffix, Bool_t kill) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // return a TH1D with the supplied histogram rebinned to the supplied bins
    // the returned histogram is new, the original is deleted from the heap if kill is true
    if(!histo || !bins) {
        printf(" > RebinTH1D:: fatal error, NULL pointer passed < \n");
        return NULL;
    }
    // create the output histo
    TString name = histo->GetName();
    name+="_template";
    name+=suffix;
    TH1D* rebinned = 0x0;
    // check if the histogram is normalized. if so, fall back to ROOT style rebinning
    // if not, rebin manually
    if(histo->GetBinContent(1) > 1 ) {
        rebinned = new TH1D(name.Data(), name.Data(), bins->GetSize()-1, bins->GetArray());
        for(Int_t i(0); i < histo->GetXaxis()->GetNbins(); i++) {
            // loop over the bins of the old histo and fill the new one with its data
            rebinned->Fill(histo->GetBinCenter(i+1), histo->GetBinContent(i+1));
        }
    } else rebinned = reinterpret_cast<TH1D*>(histo->Rebin(bins->GetSize()-1, name.Data(), bins->GetArray()));
    if(kill) delete histo;
    return rebinned;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::RebinTH2D(TH2D* rebinMe, TArrayD* binsTrue, TArrayD* binsRec, TString suffix) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // weighted rebinning of a th2d, implementation for function call to AliAnaChargedJetResponseMaker
    // not static as it is just a wrapper for the response maker object
    if(!fResponseMaker || !binsTrue || !binsRec) {
        printf(" > RebinTH2D:: function called with NULL arguments < \n");
        return 0x0;
    }
    TString name(Form("%s_%s", rebinMe->GetName(), suffix.Data()));
    return (TH2D*)fResponseMaker->MakeResponseMatrixRebin(rebinMe, (TH2*)(new TH2D(name.Data(), name.Data(), binsTrue->GetSize()-1, binsTrue->GetArray(), binsRec->GetSize()-1, binsRec->GetArray())), kTRUE);
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::MatrixMultiplication(TH2D* a, TH2D* b, TString name)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // multiply two matrices
    if (a->GetNbinsX() != b->GetNbinsY()) return 0x0;
    TH2D* c = (TH2D*)a->Clone("c");
    for (Int_t y1 = 1; y1 <= a->GetNbinsY(); y1++) {
        for (Int_t x2 = 1; x2 <= b->GetNbinsX(); x2++) {
            Double_t val = 0;
            for (Int_t x1 = 1; x1 <= a->GetNbinsX(); x1++) {
                Int_t y2 = x1;
	        val += a->GetBinContent(x1, y1) * b->GetBinContent(x2, y2);
            }
            c->SetBinContent(x2, y1, val);
            c->SetBinError(x2, y1, 0.);
        }
    }
    if(strcmp(name.Data(), "")) c->SetNameTitle(name.Data(), name.Data());
    return c;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::NormalizeTH1D(TH1D* histo, Double_t scale) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // normalize a th1d to a certain scale
    histo->Sumw2();
    Double_t integral = histo->Integral()*scale;
    if (integral > 0 && scale == 1.) histo->Scale(1./integral, "width");
    else if (scale != 1.) histo->Scale(1./scale, "width");
    else printf(" > Histogram integral < 0, cannot normalize \n");
    return histo;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::MergeSpectrumBins(TArrayI* bins, TH1D* spectrum, TH2D* corr)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // merge a spectrum histogram taking into account the correlation terms
    if(!(bins&&spectrum)) {
        printf(" > NULL pointer passed as argument in MergeSpectrumBins ! < \n");
        return 0x0;
    }
    Double_t sum(0), error(0), pearson(0);
    // take the sum of the bin content
    sum += spectrum->GetBinContent(bins->At(0));
    sum += spectrum->GetBinContent(bins->At(1));
    // quadratically sum the errors
    error += TMath::Power(spectrum->GetBinError(bins->At(0)), 2);
    error += TMath::Power(spectrum->GetBinError(bins->At(1)), 2);
    // add the covariance term
    pearson = corr->GetBinContent(bins->At(0), bins->At(1));
    if(!corr) {
        printf(" > PANIC ! something is wrong with the covariance matrix, assuming full correlation ! < \n ");
        pearson = 1;
    }
    error += 2.*spectrum->GetBinError(bins->At(0))*spectrum->GetBinError(bins->At(1))*pearson;
    spectrum->SetBinContent(1, sum);
    if(error <= 0) return spectrum;
    spectrum->SetBinError(1, TMath::Sqrt(error));
    printf(" > sum is %.2f \t error is %.8f < \n", sum, error);
    printf(" > REPLACING BIN CONTENT OF FIRST BIN (%i) WITH MERGED BINS (%i) and (%i) < \n", 1, bins->At(0), bins->At(1));
    return spectrum;
}
//_____________________________________________________________________________
TMatrixD* AliJetFlowTools::CalculatePearsonCoefficients(TMatrixD* covarianceMatrix) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // Calculate pearson coefficients from covariance matrix
    TMatrixD *pearsonCoefficients((TMatrixD*)covarianceMatrix->Clone("pearsonCoefficients"));
    Int_t nrows(covarianceMatrix->GetNrows()), ncols(covarianceMatrix->GetNcols());
    Double_t pearson(0.);
    if(nrows==0 && ncols==0) return 0x0;
    for(Int_t row = 0; row < nrows; row++) {
        for(Int_t col = 0; col<ncols; col++) {
        if((*covarianceMatrix)(row,row)!=0. && (*covarianceMatrix)(col,col)!=0.) pearson = (*covarianceMatrix)(row,col)/TMath::Sqrt((*covarianceMatrix)(row,row)*(*covarianceMatrix)(col,col));
        (*pearsonCoefficients)(row,col) = pearson;
        }
    }
    return pearsonCoefficients;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::SmoothenPrior(TH1D* spectrum, TF1* function, Double_t min, Double_t max, Double_t start, Bool_t kill, Bool_t counts) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // smoothen the spectrum using a user defined function
    // returns a clone of the original spectrum if fitting failed
    // if kill is kTRUE the input spectrum will be deleted from the heap
    // if 'count' is selected, bins are filled with integers (necessary if the 
    // histogram is interpreted in a routine which accepts only counts)
    if(!spectrum || !function) return 0x0;
    if(start > max) printf(" > cannot extrapolate fit beyond fit range ! < " );
    TH1D* temp = (TH1D*)spectrum->Clone(Form("%s_smoothened", spectrum->GetName()));
    temp->Sumw2();      // if already called on the original, this will give off a warning but do nothing
    TFitResultPtr r = temp->Fit(function, "", "", min, max);
    if((int)r == 0) {   // MINUIT status
        for(Int_t i(0); i < temp->GetNbinsX() + 1; i++) {
            if(temp->GetBinCenter(i) > start) {     // from this pt value use extrapolation
                (counts) ? temp->SetBinContent(i, (int)(function->Integral(temp->GetXaxis()->GetBinLowEdge(i),temp->GetXaxis()->GetBinUpEdge(i))/temp->GetXaxis()->GetBinWidth(i))) : temp->SetBinContent(i, function->Integral(temp->GetXaxis()->GetBinLowEdge(i),temp->GetXaxis()->GetBinUpEdge(i))/temp->GetXaxis()->GetBinWidth(i));
                if(temp->GetBinContent(i) > 0) temp->SetBinError(i, TMath::Sqrt(temp->GetBinContent(i)));
            }
        }
    }
    if(kill) delete spectrum;
    return temp;
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(Bool_t legacy)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // set global style for your current aliroot session
    if(!gStyle) return;
    // legacy style is pleasing to the eye, default is the formal ALICE style
    if(legacy) {
        gStyle->SetCanvasColor(-1); 
        gStyle->SetPadColor(-1); 
        gStyle->SetFrameFillColor(-1); 
        gStyle->SetHistFillColor(-1); 
        gStyle->SetTitleFillColor(-1); 
        gStyle->SetFillColor(-1); 
        gStyle->SetFillStyle(4000); 
        gStyle->SetStatStyle(0); 
        gStyle->SetTitleStyle(0); 
        gStyle->SetCanvasBorderSize(0); 
        gStyle->SetFrameBorderSize(0); 
        gStyle->SetLegendBorderSize(0); 
        gStyle->SetStatBorderSize(0); 
        gStyle->SetTitleBorderSize(0);
    } else {
        gStyle->Reset("Plain");
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        gStyle->SetCanvasColor(10);
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetFrameLineWidth(1);
        gStyle->SetFrameFillColor(kWhite);
        gStyle->SetPadColor(10);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadBottomMargin(0.17);
        gStyle->SetPadLeftMargin(0.15);
        gStyle->SetHistLineWidth(1);
        gStyle->SetHistLineColor(kRed);
        gStyle->SetFuncWidth(2);
        gStyle->SetFuncColor(kGreen);
        gStyle->SetLineWidth(2);
        gStyle->SetLabelSize(0.045,"xyz");
        gStyle->SetLabelOffset(0.01,"y");
        gStyle->SetLabelOffset(0.01,"x");
        gStyle->SetLabelColor(kBlack,"xyz");
        gStyle->SetTitleSize(0.05,"xyz");
        gStyle->SetTitleOffset(1.25,"y");
        gStyle->SetTitleOffset(1.2,"x");
        gStyle->SetTitleFillColor(kWhite);
        gStyle->SetTextSizePixels(26);
        gStyle->SetTextFont(42);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(kWhite);
        gStyle->SetLegendFont(42);
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TCanvas* c, TString style)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    if(c) return;
    // set a default style for a canvas
    if(!strcmp(style.Data(), "PEARSON")) {
        printf(" > style PEARSON canvas < \n");
        gStyle->SetOptStat(0);
        c->SetGridx();
        c->SetGridy();
        c->SetTicks();
        return;
    } else if(!strcmp(style.Data(), "SPECTRUM")) {
        printf(" > style SPECTRUM canvas < \n");
        gStyle->SetOptStat(0);
        c->SetLogy();
        c->SetGridx();
        c->SetGridy();
        c->SetTicks();
        return;
    } else printf(" > Style called with unknown option %s \n    returning < \n", style.Data());
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TVirtualPad* c, TString style, Bool_t legacy)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // set a default style for a canva
    if(!c) return;
    if(legacy) {
        c->SetLeftMargin(.25);
        c->SetBottomMargin(.25);
    }
    else Style();
    if(!strcmp(style.Data(), "PEARSON")) {
        printf(" > style PEARSON pad < \n");
        gStyle->SetOptStat(0);
        c->SetGridx();
        c->SetGridy();
        c->SetTicks();
        return;
    } else if(!strcmp(style.Data(), "SPECTRUM")) {
        printf(" > style SPECTRUM pad < \n");
        gStyle->SetOptStat(0);
        c->SetLogy();
        c->SetGridx();
        c->SetGridy();
        c->SetTicks();
        return;
    } else if (!strcmp(style.Data(), "GRID")) {
        printf(" > style GRID pad < \n");
        gStyle->SetOptStat(0);
        c->SetGridx();
        c->SetGridy();
        c->SetTicks();
    } else printf(" > Style called with unknown option %s \n    returning < \n", style.Data());
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TLegend* l)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // set a default style for a legend
    if(!l) return;
    l->SetFillColor(0);
    l->SetBorderSize(0);
    if(gStyle) l->SetTextSize(gStyle->GetTextSize()*.08);
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TH1* h, EColor col, histoType type, Bool_t legacy)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // style a histo
    if(!h) return;
    h->GetYaxis()->SetNdivisions(505);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2);
    h->SetMarkerSize(1.5);
    if(legacy) {
        h->SetTitle("");
        h->GetYaxis()->SetLabelSize(0.05);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetTitleOffset(1.5);
        h->GetXaxis()->SetTitleOffset(1.5);
        h->GetYaxis()->SetTitleSize(.05);
        h->GetXaxis()->SetTitleSize(.05);
    } else Style();
    switch (type) {
        case kInPlaneSpectrum : {
            h->SetTitle("IN PLANE");
            h->GetXaxis()->SetTitle("#it{p}_{T}^{ch, jet} (GeV/#it{c})");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
        } break;
        case kOutPlaneSpectrum : {
            h->SetTitle("OUT OF PLANE");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
            h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
       } break;
       case kUnfoldedSpectrum : {
            h->SetTitle("UNFOLDED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
            h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
       } break;
       case kFoldedSpectrum : {
            h->SetTitle("FOLDED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
            h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
       } break;
       case kMeasuredSpectrum : {
            h->SetTitle("MEASURED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
            h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{cht} (GeV/#it{c})");
       } break;
       case kBar : {
            h->SetFillColor(col);
            h->SetBarWidth(.6);
            h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
            h->SetBarOffset(0.2);
       }
       case kRatio : {
            h->SetMarkerStyle(8);
            h->SetMarkerSize(1.5);
       } break;
       case kDeltaPhi : {
            h->GetYaxis()->SetTitle("[counts]");
            h->GetXaxis()->SetTitle("#Delta #phi");
       }
       default : break;
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TGraph* h, EColor col, histoType type, Bool_t legacy)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // style a tgraph
    if(!h) return;
    h->GetYaxis()->SetNdivisions(505);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2);
    h->SetMarkerSize(1.5);
    h->SetTitle("");
    h->SetFillColor(kCyan);
    if(legacy) {
        h->GetYaxis()->SetLabelSize(0.05);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->GetXaxis()->SetTitleOffset(1.6);
        h->GetYaxis()->SetTitleSize(.05);
        h->GetXaxis()->SetTitleSize(.05);
    } else Style();
    h->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
    switch (type) {
        case kInPlaneSpectrum : {
            h->SetTitle("IN PLANE");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
        } break;
        case kOutPlaneSpectrum : {
            h->SetTitle("OUT OF PLANE");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
       } break;
       case kUnfoldedSpectrum : {
            h->SetTitle("UNFOLDED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
       } break;
       case kFoldedSpectrum : {
            h->SetTitle("FOLDED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
       } break;
       case kMeasuredSpectrum : {
            h->SetTitle("MEASURED");
            h->GetYaxis()->SetTitle("#frac{d#it{N}}{d#it{p}_{T}}");
       } break;
       case kRatio : {
            h->GetYaxis()->SetTitle("(d#it{N}^{ch, jet}_{in plane}/(d#it{p}_{T}d#eta))/(d#it{N}^{ch,jet}_{out of plane}/(d#it{p}_{T}d#eta))");
       } break;
       case kV2 : {
            h->GetYaxis()->SetTitle("#it{v}_{2}^{ch, jet} {EP, |#Delta#eta|>0.9 } ");
            h->GetYaxis()->SetRangeUser(-.5, 1.);
            h->SetMarkerStyle(8);
            h->SetMarkerSize(1.5);
       } break;
       default : break;
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::GetNominalValues(
        TH1D*& ratio,           // pointer reference, output of this function
        TGraphErrors*& v2,      // pointer reference, as output of this function
        TArrayI* in,
        TArrayI* out,
        TString inFile,
        TString outFile) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // pass clones of the nominal points and nominal v2 values 
    if(fOutputFile && !fOutputFile->IsZombie()) fOutputFile->Close();   // if for some weird reason the unfolding output is still mutable
    TFile* readMe(new TFile(inFile.Data(), "READ"));                    // open unfolding output read-only
    if(readMe->IsZombie()) {
        printf(" > Fatal error, couldn't read %s for post processing ! < \n", inFile.Data());
        return;
    }
    printf("\n\n\n\t\t GetNominalValues \n > Recovered the following file structure : \n <");
    readMe->ls();
    TFile* output(new TFile(outFile.Data(), "RECREATE"));   // create new output file
    if(output) SquelchWarning();
    // get some placeholders, have to be initialized but will be deleted
    ratio = new TH1D("nominal", "nominal", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
    TH1D* nominalIn(new TH1D("nominal in", "nominal in", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* nominalOut(new TH1D("nominal out", "nominal out", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    Int_t iIn[] = {in->At(0), in->At(0)};          // first value is the nominal point
    Int_t iOut[] = {out->At(0), out->At(0)};

    // call the functions and set the relevant pointer references
    TH1D* dud(0x0);
    DoIntermediateSystematics(
            new TArrayI(2, iIn), 
            new TArrayI(2, iOut),
            dud, dud, dud, dud, dud, dud, 
            ratio,              // pointer reference, output of this function 
            nominalIn,
            nominalOut,
            1, 
            fBinsTrue->At(0), 
            fBinsTrue->At(fBinsTrue->GetSize()-1),
            readMe,
            "nominal_values");
    v2 = GetV2(nominalIn, nominalOut, fEventPlaneRes, "nominal v_{2}");

    // close the open files, reclaim ownership of histograms which are necessary outside of the file scope
    ratio->SetDirectory(0);     // disassociate from current gDirectory
    readMe->Close();
}
//_____________________________________________________________________________
void AliJetFlowTools::GetCorrelatedUncertainty(
        TGraphAsymmErrors*& corrRatio,  // correlated uncertainty pointer reference
        TGraphAsymmErrors*& corrV2,     // correlated uncertainty pointer reference
        TArrayI* variationsIn,          // variantions in plane
        TArrayI* variationsOut,         // variantions out of plane
        Bool_t sym,                     // treat as symmmetric
        TArrayI* variations2ndIn,       // second source of variations
        TArrayI* variations2ndOut,      // second source of variations
        Bool_t sym2nd,                  // treat as symmetric
        TString type,                   // type of uncertaitny
        TString type2,
        Int_t columns,                  // divide the output canvasses in this many columns
        Float_t rangeLow,               // lower pt range
        Float_t rangeUp,                // upper pt range
        Float_t corr,                   // correlation strength
        TString in,                     // input file name (created by this unfolding class)
        TString out                     // output file name (which will hold results of the systematic test)
        ) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // do full systematics
    if(fOutputFile && !fOutputFile->IsZombie()) fOutputFile->Close();   // if for some weird reason the unfolding output is still mutable
    TFile* readMe(new TFile(in.Data(), "READ"));                        // open unfolding output read-only
    if(readMe->IsZombie()) {
        printf(" > Fatal error, couldn't read %s for post processing ! < \n", in.Data());
        return;
    }
    printf("\n\n\n\t\t GetCorrelatedUncertainty \n > Recovered the following file structure : \n <");
    readMe->ls();
    TFile* output(new TFile(out.Data(), "RECREATE"));   // create new output file

    if(sym2nd) SquelchWarning();

    // create some null placeholder pointers
    TH1D* relativeErrorVariationInLow(0x0);
    TH1D* relativeErrorVariationInUp(0x0);
    TH1D* relativeErrorVariationOutLow(0x0);
    TH1D* relativeErrorVariationOutUp(0x0);
    TH1D* relativeError2ndVariationInLow(0x0);
    TH1D* relativeError2ndVariationInUp(0x0);
    TH1D* relativeError2ndVariationOutLow(0x0);
    TH1D* relativeError2ndVariationOutUp(0x0);
    TH1D* relativeStatisticalErrorIn(0x0);
    TH1D* relativeStatisticalErrorOut(0x0);
    // histo for the nominal ratio in / out
    TH1D* nominal(new TH1D("ratio in plane, out of plane", "ratio in plane, out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* nominalIn(new TH1D("in plane jet yield", "in plane jet yield", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* nominalOut(new TH1D("out of plane jet yield", "out of plane jet yield", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));

    // call the functions
    if(variationsIn && variationsOut) {
        SystematicsWrapper(
                variationsIn, 
                variationsOut, 
                relativeErrorVariationInUp,        // pointer reference
                relativeErrorVariationInLow,       // pointer reference
                relativeErrorVariationOutUp,       // pointer reference
                relativeErrorVariationOutLow,      // pointer reference
                relativeStatisticalErrorIn,        // pointer reference
                relativeStatisticalErrorOut,       // pointer reference
                nominal,
                nominalIn,
                nominalOut,
                columns, 
                rangeLow, 
                rangeUp,
                readMe,
                type,
                kFALSE,
                kTRUE);
        if(relativeErrorVariationInUp) {
            // canvas with the error from variation strength
            TCanvas* relativeErrorVariation(new TCanvas(Form("relativeError_%s", type.Data()), Form("relativeError_%s", type.Data())));
            relativeErrorVariation->Divide(2);
            relativeErrorVariation->cd(1);
            Style(gPad, "GRID");
            relativeErrorVariationInUp->DrawCopy("b");
            relativeErrorVariationInLow->DrawCopy("same b");
            Style(AddLegend(gPad));
            relativeErrorVariation->cd(2);
            Style(gPad, "GRID");
            relativeErrorVariationOutUp->DrawCopy("b");
            relativeErrorVariationOutLow->DrawCopy("same b");
            SavePadToPDF(relativeErrorVariation);
            Style(AddLegend(gPad));
            relativeErrorVariation->Write();

            // now smoothen the diced response error (as it is expected to be flat)
            // this is done by fitting a constant to the diced resonse error histo
            //
            TF1* lin = new TF1("lin", "[0]", rangeLow, rangeUp);
            relativeErrorVariationInUp->Fit(lin, "L", "", rangeLow, rangeUp);
            if(gMinuit->fISW[1] != 3) printf(" fit is NOT ok ! " );
            for(Int_t i(0); i < relativeErrorVariationInUp->GetNbinsX(); i++) {
                relativeErrorVariationInUp->SetBinContent(i+1, lin->GetParameter(0));
            }
            relativeErrorVariationInLow->Fit(lin, "L", "", rangeLow, rangeUp);
            printf(" > Fit over diced resonse, new value for all bins is %.4f < \n ", lin->GetParameter(0));
            for(Int_t i(0); i < relativeErrorVariationInUp->GetNbinsX(); i++) {
                relativeErrorVariationInLow->SetBinContent(i+1, 0);//lin->GetParameter(0));
            }
            relativeErrorVariationOutUp->Fit(lin, "L", "", rangeLow, rangeUp);
            printf(" > Fit over diced resonse, new value for all bins is %.4f < \n ", lin->GetParameter(0));
            for(Int_t i(0); i < relativeErrorVariationInUp->GetNbinsX(); i++) {
                relativeErrorVariationOutUp->SetBinContent(i+1, lin->GetParameter(0));
            }
            relativeErrorVariationOutLow->Fit(lin, "L", "", rangeLow, rangeUp);
            printf(" > Fit over diced resonse, new value for all bins is %.4f < \n ", lin->GetParameter(0));
            for(Int_t i(0); i < relativeErrorVariationInUp->GetNbinsX(); i++) {
                relativeErrorVariationOutLow->SetBinContent(i+1, 0);//lin->GetParameter(0));
            }
        }
    }
    // call the functions for a second set of systematic sources
    if(variations2ndIn && variations2ndOut) {
        SystematicsWrapper(
                variations2ndIn, 
                variations2ndOut, 
                relativeError2ndVariationInUp,        // pointer reference
                relativeError2ndVariationInLow,       // pointer reference
                relativeError2ndVariationOutUp,       // pointer reference
                relativeError2ndVariationOutLow,      // pointer reference
                relativeStatisticalErrorIn,           // pointer reference
                relativeStatisticalErrorOut,          // pointer reference
                nominal,
                nominalIn,
                nominalOut,
                columns, 
                rangeLow, 
                rangeUp,
                readMe,
                type2,
                kFALSE,
                kTRUE);
        if(relativeError2ndVariationInUp) {
            // canvas with the error from variation strength
            TCanvas* relativeError2ndVariation(new TCanvas(Form("relativeError2nd_%s", type2.Data()), Form("relativeError2nd_%s", type2.Data())));
            relativeError2ndVariation->Divide(2);
            relativeError2ndVariation->cd(1);
            Style(gPad, "GRID");
            relativeError2ndVariationInUp->DrawCopy("b");
            relativeError2ndVariationInLow->DrawCopy("same b");
            Style(AddLegend(gPad));
            relativeError2ndVariation->cd(2);
            Style(gPad, "GRID");
            relativeError2ndVariationOutUp->DrawCopy("b");
            relativeError2ndVariationOutLow->DrawCopy("same b");
            SavePadToPDF(relativeError2ndVariation);
            Style(AddLegend(gPad));
            relativeError2ndVariation->Write();
        }

    }

    // and the placeholder for the final systematic
    TH1D* relativeErrorInUp(new TH1D("max correlated uncertainty in plane", "max correlated uncertainty in plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorOutUp(new TH1D("max correlated uncertainty out of plane", "max correlated uncertainty out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorInLow(new TH1D("min correlated uncertainty in plane", "min correlated uncertainty in plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorOutLow(new TH1D("min correlated uncertainty out of plane", "min correlated uncertainty out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    relativeErrorInUp->SetYTitle("relative uncertainty");
    relativeErrorOutUp->SetYTitle("relative uncertainty");
    relativeErrorInLow->SetYTitle("relative uncertainty");
    relativeErrorOutLow->SetYTitle("relative uncertainty");

    // merge the correlated errors (FIXME) trivial for one set 
    Double_t aInUp(0.), bInUp(0.), cInUp(0.), dInUp(0.);
    Double_t aOutUp(0.), bOutUp(0.), cOutUp(0.), dOutUp(0.);
    Double_t aInLow(0.), bInLow(0.), cInLow(0.), dInLow(0.);
    Double_t aOutLow(0.), bOutLow(0.), cOutLow(0.), dOutLow(0.);

    for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
        // for the upper bound
        if(relativeErrorVariationInUp) aInUp = relativeErrorVariationInUp->GetBinContent(b+1);
        if(relativeErrorVariationOutUp) aOutUp = relativeErrorVariationOutUp->GetBinContent(b+1);
        if(relativeError2ndVariationInUp) bInUp = relativeError2ndVariationInUp->GetBinContent(b+1);
        if(relativeError2ndVariationOutUp) bInLow = relativeError2ndVariationOutUp->GetBinContent(b+1);
        dInUp  = aInUp*aInUp + bInUp*bInUp + cInUp*cInUp;
        // for a symmetric first variation
        if(sym) dInUp += aInLow*aInLow;
        if(dInUp > 0) relativeErrorInUp->SetBinContent(b+1, TMath::Sqrt(dInUp));
        dOutUp = aOutUp*aOutUp + bOutUp*bOutUp + cOutUp*cOutUp;
        if(sym) dOutUp += aOutLow*aOutLow;
        if(dOutUp > 0) relativeErrorOutUp->SetBinContent(b+1, TMath::Sqrt(dOutUp));
        // for the lower bound
        if(relativeErrorVariationInLow) aInLow = relativeErrorVariationInLow->GetBinContent(b+1);
        if(relativeErrorVariationOutLow) aOutLow = relativeErrorVariationOutLow->GetBinContent(b+1);
        if(relativeError2ndVariationInLow) bInLow = relativeError2ndVariationInLow->GetBinContent(b+1);
        if(relativeError2ndVariationOutLow) bOutLow = relativeError2ndVariationOutLow->GetBinContent(b+1);
        dInLow  = aInLow*aInLow + bInLow*bInLow + cInLow*cInLow;
        if(sym) dInLow += aInUp*aInUp;
        if(dInLow > 0) relativeErrorInLow->SetBinContent(b+1, -1*TMath::Sqrt(dInLow));
        dOutLow = aOutLow*aOutLow + bOutLow*bOutLow + cOutLow*cOutLow;
        if(sym) dOutLow += aOutUp*aOutUp;
        if(dOutLow > 0) relativeErrorOutLow->SetBinContent(b+1, -1.*TMath::Sqrt(dOutLow));
    }
    // project the estimated errors on the nominal ratio
    if(nominal) {
        Double_t* ax = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ay = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* axh = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* axl = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ayh = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ayl = new Double_t[fBinsTrue->GetSize()-1];
        Double_t lowErr(0.), upErr(0.);
        for(Int_t i(0); i < fBinsTrue->GetSize()-1; i++) {
            // add the in and out of plane errors in quadrature
            // propagation is tricky because we should consider two cases:
            // [1] the error is uncorrelated, and we propagate upper 1 with lower 2 and vice versa
            // [2] the error is correlated - in this case upper 1 should be propagated with upper 2
            // as the fluctuations are bound to always to of same sign
            if(corr <= 0) {
                lowErr = relativeErrorInLow->GetBinContent(i+1)*relativeErrorInLow->GetBinContent(i+1)+relativeErrorOutUp->GetBinContent(1+i)*relativeErrorOutUp->GetBinContent(i+1);
                upErr = relativeErrorInUp->GetBinContent(i+1)*relativeErrorInUp->GetBinContent(i+1)+relativeErrorOutLow->GetBinContent(i+1)*relativeErrorOutLow->GetBinContent(i+1);
            } else {
                lowErr = relativeErrorInLow->GetBinContent(i+1)*relativeErrorInLow->GetBinContent(i+1)+relativeErrorOutLow->GetBinContent(1+i)*relativeErrorOutLow->GetBinContent(i+1) -2.*corr*relativeErrorInLow->GetBinContent(i+1)*relativeErrorOutLow->GetBinContent(i+1);
                upErr = relativeErrorInUp->GetBinContent(i+1)*relativeErrorInUp->GetBinContent(i+1)+relativeErrorOutUp->GetBinContent(i+1)*relativeErrorOutUp->GetBinContent(i+1) - 2.*corr*relativeErrorInUp->GetBinContent(i+1)*relativeErrorOutUp->GetBinContent(i+1);
            }
            ayl[i] = TMath::Sqrt(TMath::Abs(lowErr))*nominal->GetBinContent(i+1);
            ayh[i] = TMath::Sqrt(TMath::Abs(upErr))*nominal->GetBinContent(i+1);
            // get the bin width (which is the 'error' on x
            Double_t binWidth(nominal->GetBinWidth(i+1));
            axl[i] = binWidth/2.;
            axh[i] = binWidth/2.;
            // now get the coordinate for the point
            ax[i] = nominal->GetBinCenter(i+1);
            ay[i] = nominal->GetBinContent(i+1);
        }
        // save the nominal ratio
        TGraphAsymmErrors* nominalError(new TGraphAsymmErrors(fBinsTrue->GetSize()-1, ax, ay, axl, axh, ayl, ayh));
        corrRatio = (TGraphAsymmErrors*)nominalError->Clone();
        nominalError->SetName("correlated uncertainty");
        TCanvas* nominalCanvas(new TCanvas("nominalCanvas", "nominalCanvas"));
        nominalCanvas->Divide(2);
        nominalCanvas->cd(1);
        Style(nominal, kBlack);
        Style(nominalError, kCyan, kRatio);
        nominalError->SetLineColor(kCyan);
        nominalError->SetMarkerColor(kCyan);
        nominalError->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
        nominalError->GetYaxis()->SetRangeUser(.7, 2.2);
        nominalError->DrawClone("a2");
        nominal->DrawCopy("same E1");
        Style(AddLegend(gPad));
        Style(gPad, "GRID");
        Style(nominalCanvas);
        // save nominal v2 and systematic errors
        TGraphAsymmErrors* nominalV2Error(GetV2WithSystematicErrors(
                    nominalIn,
                    nominalOut,
                    fEventPlaneRes,
                    "v_{2} with shape uncertainty",
                    relativeErrorInUp,
                    relativeErrorInLow,
                    relativeErrorOutUp,
                    relativeErrorOutLow,
                    corr));
        // pass the nominal values to the pointer references
        corrV2 = (TGraphAsymmErrors*)nominalV2Error->Clone();
        TGraphErrors* nominalV2(GetV2(nominalIn, nominalOut, fEventPlaneRes, "v_{2}"));
        nominalCanvas->cd(2);
        Style(nominalV2, kBlack);
        Style(nominalV2Error, kCyan, kV2);
        nominalV2Error->SetLineColor(kCyan);
        nominalV2Error->SetMarkerColor(kCyan);
        nominalV2Error->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
        nominalV2Error->DrawClone("a2");
        nominalV2->DrawClone("same E1");
        Style(AddLegend(gPad));
        Style(gPad, "GRID");
        Style(nominalCanvas);
        SavePadToPDF(nominalCanvas);
        nominalCanvas->Write();
    }

    TCanvas* relativeError(new TCanvas("relativeCorrelatedError"," relativeCorrelatedError"));
    relativeError->Divide(2);
    relativeError->cd(1);
    Style(gPad, "GRID");
    relativeErrorInUp->GetYaxis()->SetRangeUser(-1.5, 3.);
    Style(relativeErrorInUp, kBlue, kBar);
    Style(relativeErrorInLow, kGreen, kBar);
    relativeErrorInUp->DrawCopy("b");
    relativeErrorInLow->DrawCopy("same b");
    Style(relativeStatisticalErrorIn, kRed);
    relativeStatisticalErrorIn->DrawCopy("same");
    Style(AddLegend(gPad));
    relativeError->cd(2);
    Style(gPad, "GRID");
    relativeErrorOutUp->GetYaxis()->SetRangeUser(-1.5, 3.);
    Style(relativeErrorOutUp, kBlue, kBar);
    Style(relativeErrorOutLow, kGreen, kBar);
    relativeErrorOutUp->DrawCopy("b");
    relativeErrorOutLow->DrawCopy("same b");
    Style(relativeStatisticalErrorOut, kRed);
    relativeStatisticalErrorOut->DrawCopy("same");
    Style(AddLegend(gPad));

    // write the buffered file to disk and close the file
    SavePadToPDF(relativeError);
    relativeError->Write();
    output->Write();
}
//_____________________________________________________________________________
void AliJetFlowTools::GetShapeUncertainty(
        TGraphAsymmErrors*& shapeRatio, // pointer reference to final shape uncertainty of ratio
        TGraphAsymmErrors*& shapeV2,    // pointer reference to final shape uncertainty of v2
        TArrayI* regularizationIn,      // regularization strength systematics, in plane
        TArrayI* regularizationOut,     // regularization strenght systematics, out of plane
        TArrayI* recBinIn,              // rec bin ranges, in plane
        TArrayI* recBinOut,             // rec bin ranges, out of plane
        TArrayI* methodIn,              // method in
        TArrayI* methodOut,             // method out
        Int_t columns,                  // divide the output canvasses in this many columns
        Float_t rangeLow,               // lower pt range
        Float_t rangeUp,                // upper pt range
        Float_t corr,                   // correlation strength
        TString in,                     // input file name (created by this unfolding class)
        TString out,                    // output file name (which will hold results of the systematic test)
        Bool_t regularizationOnV2       // get uncertainty on yields separately or v2 directly
        ) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // do full systematics
    if(fOutputFile && !fOutputFile->IsZombie()) fOutputFile->Close();   // if for some weird reason the unfolding output is still mutable
    TFile* readMe(new TFile(in.Data(), "READ"));                        // open unfolding output read-only
    if(readMe->IsZombie()) {
        printf(" > Fatal error, couldn't read %s for post processing ! < \n", in.Data());
        return;
    }
    printf("\n\n\n\t\t DOSYSTEMATICS \n > Recovered the following file structure : \n <");
    readMe->ls();
    TFile* output(new TFile(out.Data(), "RECREATE"));   // create new output file
    // create some null placeholder pointers
    TH1D* relativeErrorRegularizationInLow(0x0);
    TH1D* relativeErrorRegularizationInUp(0x0);
    TH1D* relativeErrorRecBinInLow(0x0);
    TH1D* relativeErrorRecBinInUp(0x0);
    TH1D* relativeErrorMethodInLow(0x0);
    TH1D* relativeErrorMethodInUp(0x0);
    TH1D* relativeErrorRegularizationOutLow(0x0);
    TH1D* relativeErrorRegularizationOutUp(0x0);
    TH1D* relativeErrorRecBinOutLow(0x0);
    TH1D* relativeErrorRecBinOutUp(0x0);
    TH1D* relativeStatisticalErrorIn(0x0);
    TH1D* relativeStatisticalErrorOut(0x0);
    TH1D* relativeErrorMethodOutLow(0x0);
    TH1D* relativeErrorMethodOutUp(0x0);
    // histo for the nominal ratio in / out
    TH1D* nominal(new TH1D("ratio in plane, out of plane", "ratio in plane, out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* nominalIn(new TH1D("in plane jet yield", "in plane jet yield", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* nominalOut(new TH1D("out of plane jet yield", "out of plane jet yield", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));

    // call the functions
    if(regularizationIn && regularizationOut) {
        SystematicsWrapper(
                regularizationIn, 
                regularizationOut, 
                relativeErrorRegularizationInUp,        // pointer reference
                relativeErrorRegularizationInLow,       // pointer reference
                relativeErrorRegularizationOutUp,       // pointer reference
                relativeErrorRegularizationOutLow,      // pointer reference
                relativeStatisticalErrorIn,             // pointer reference
                relativeStatisticalErrorOut,            // pointer reference
                nominal,
                nominalIn,
                nominalOut,
                columns, 
                rangeLow, 
                rangeUp,
                readMe,
                "regularization",
                fRMS,
                !regularizationOnV2);   // set to false means NO undertainty from ratio, but v2
       if(relativeErrorRegularizationInUp && !regularizationOnV2 ) {
            // canvas with the error from regularization strength
            TCanvas* relativeErrorRegularization(new TCanvas("relativeErrorRegularization", "relativeErrorRegularization"));
            relativeErrorRegularization->Divide(2);
            relativeErrorRegularization->cd(1);
            Style(gPad, "GRID");
            relativeErrorRegularizationInUp->DrawCopy("b");
            relativeErrorRegularizationInLow->DrawCopy("same b");
            Style(AddLegend(gPad));
            relativeErrorRegularization->cd(2);
            Style(gPad, "GRID");
            relativeErrorRegularizationOutUp->DrawCopy("b");
            relativeErrorRegularizationOutLow->DrawCopy("same b");
            SavePadToPDF(relativeErrorRegularization);
            Style(AddLegend(gPad));
            relativeErrorRegularization->Write();
       } else if(relativeErrorRegularizationInUp && regularizationOnV2 ) {
            // canvas with the error from regularization strength
            TCanvas* relativeErrorRegularization(new TCanvas("relativeErrorRegularization", "relativeErrorRegularization"));
            Style(gPad, "GRID");
            relativeErrorRegularization->cd(1);
            TH1F* relativeErrorRegularizationInUpErrors = (TH1F*)relativeErrorRegularizationInUp->Clone();
            for(Int_t i(1); i < relativeErrorRegularizationInUp->GetNbinsX() + 1; i++) {
                relativeErrorRegularizationInUpErrors->SetBinContent(i, relativeErrorRegularizationInUp->GetBinError(i));
                relativeErrorRegularizationInUpErrors->SetBinError(i, 0);
            }
            relativeErrorRegularizationInUpErrors->GetYaxis()->SetTitle("absolute error on v_{2} from unfolding");
            relativeErrorRegularizationInUpErrors->GetXaxis()->SetTitle("#it{p}_{T, jet}^{ch} (GeV/#it{c})");
            relativeErrorRegularizationInUpErrors->DrawCopy("b");
            Style(gPad, "GRID");
            relativeErrorRegularization->Write();
        }
    }
    if(recBinIn && recBinOut) {
        SystematicsWrapper(
                recBinIn, 
                recBinOut, 
                relativeErrorRecBinInUp,       // pointer reference
                relativeErrorRecBinInLow,        // pointer reference
                relativeErrorRecBinOutUp,      // pointer reference
                relativeErrorRecBinOutLow,       // pointer reference
                relativeStatisticalErrorIn,
                relativeStatisticalErrorOut,
                nominal,
                nominalIn,
                nominalOut,
                columns, 
                rangeLow, 
                rangeUp,
                readMe,
                "recBin",
                kFALSE,//fRMS,  // was RMS< i think not ok ...
                kTRUE);         // undertainty from ratio
       if(relativeErrorRecBinOutUp) {
            // canvas with the error from regularization strength
            TCanvas* relativeErrorRecBin(new TCanvas("relativeErrorRecBin"," relativeErrorRecBin"));
            relativeErrorRecBin->Divide(2);
            relativeErrorRecBin->cd(1);
            Style(gPad, "GRID");
            relativeErrorRecBinInUp->DrawCopy("b");
            relativeErrorRecBinInLow->DrawCopy("same b");
            Style(AddLegend(gPad));
            relativeErrorRecBin->cd(2);
            Style(gPad, "GRID");
            relativeErrorRecBinOutUp->DrawCopy("b");
            relativeErrorRecBinOutLow->DrawCopy("same b");
            SavePadToPDF(relativeErrorRecBin);
            Style(AddLegend(gPad));
            relativeErrorRecBin->Write();
        }
    }
    if(methodIn && methodOut) {
        SystematicsWrapper(
                methodIn, 
                methodOut, 
                relativeErrorMethodInUp,       // pointer reference
                relativeErrorMethodInLow,      // pointer reference
                relativeErrorMethodOutUp,      // pointer reference
                relativeErrorMethodOutLow,     // pointer reference
                relativeStatisticalErrorIn,
                relativeStatisticalErrorOut,
                nominal,
                nominalIn,
                nominalOut,
                columns, 
                rangeLow, 
                rangeUp,
                readMe,
                "method",
                kFALSE,
                kTRUE);         // undertaitny from ratio
      if(relativeErrorMethodInUp) {
            TCanvas* relativeErrorMethod(new TCanvas("relativeErrorMethod", "relativeErrorMethod"));
            relativeErrorMethod->Divide(2);
            relativeErrorMethod->cd(1);
            Style(gPad, "GRID");
            relativeErrorMethodInUp->DrawCopy("b");
            relativeErrorMethodInLow->DrawCopy("same b");
            Style(AddLegend(gPad));
            relativeErrorMethod->cd(2);
            Style(gPad, "GRID");
            relativeErrorMethodOutUp->DrawCopy("b");
            relativeErrorMethodOutLow->DrawCopy("same b");
            SavePadToPDF(relativeErrorMethod);
            Style(AddLegend(gPad));
            relativeErrorMethod->Write();
        }
    }

    // and the placeholder for the final systematic
    TH1D* relativeErrorInUp(new TH1D("max shape uncertainty in plane", "max shape uncertainty in plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorOutUp(new TH1D("max shape uncertainty out of plane", "max shape uncertainty out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorInLow(new TH1D("min shape uncertainty in plane", "min shape uncertainty in plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    TH1D* relativeErrorOutLow(new TH1D("min shape uncertainty out of plane", "min shape uncertainty out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    relativeErrorInUp->SetYTitle("relative uncertainty");
    relativeErrorOutUp->SetYTitle("relative uncertainty");
    relativeErrorInLow->SetYTitle("relative uncertainty");
    relativeErrorOutLow->SetYTitle("relative uncertainty");

    // sum of squares for the total systematic error
    Double_t aInUp(0.), cInUp(0.), dInUp(0.), eInUp(0.);
    Double_t aOutUp(0.), cOutUp(0.), dOutUp(0.), eOutUp(0.);
    Double_t aInLow(0.), cInLow(0.), dInLow(0.), eInLow(0.);
    Double_t aOutLow(0.), cOutLow(0.), dOutLow(0.), eOutLow(0.);

    GetErrorFromFit(relativeErrorRecBinInUp, relativeErrorRecBinInLow, rangeLow, rangeUp, fPivot, fSubdueError, "measured, in-plane");
    GetErrorFromFit(relativeErrorRecBinOutUp, relativeErrorRecBinOutLow, rangeLow, rangeUp, fPivot, fSubdueError, "measured, out-of-plane");
    GetErrorFromFit(relativeErrorMethodInUp, relativeErrorMethodInLow, rangeLow, rangeUp, fPivot, fSubdueError, "background, in-plane");
    GetErrorFromFit(relativeErrorMethodOutUp, relativeErrorMethodOutLow, rangeLow, rangeUp, fPivot, fSubdueError, "background, out-of-plane");

    for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
        // for the upper bound only add regularization stuff if it is NOT done on v2 directly
        if(relativeErrorRegularizationInUp && !regularizationOnV2) aInUp = relativeErrorRegularizationInUp->GetBinContent(b+1);
        if(relativeErrorRegularizationOutUp && !regularizationOnV2) aOutUp = relativeErrorRegularizationOutUp->GetBinContent(b+1); 
        if(relativeErrorRecBinInUp) cInUp = relativeErrorRecBinInUp->GetBinContent(b+1);
        if(relativeErrorRecBinOutUp) cOutUp = relativeErrorRecBinOutUp->GetBinContent(b+1);
        if(relativeErrorMethodInUp) dInUp = relativeErrorMethodInUp->GetBinContent(b+1);
        if(relativeErrorMethodOutUp) dOutUp = relativeErrorMethodOutUp->GetBinContent(b+1);
        // for the lower bound
        if(relativeErrorRegularizationInLow) aInLow = relativeErrorRegularizationInLow->GetBinContent(b+1);
        if(relativeErrorRegularizationOutLow) aOutLow = relativeErrorRegularizationOutLow->GetBinContent(b+1);
        if(relativeErrorRecBinInLow) cInLow = relativeErrorRecBinInLow->GetBinContent(b+1);
        if(relativeErrorRecBinOutLow) cOutLow = relativeErrorRecBinOutLow->GetBinContent(b+1);
        if(relativeErrorMethodInLow) dInLow = relativeErrorMethodInLow->GetBinContent(b+1);
        if(relativeErrorMethodOutLow) dOutLow = relativeErrorMethodOutLow->GetBinContent(b+1);

        if(fSymmRMS) {  // since there is no a-priori assumption, we take these errors as symmetric
            RemoveSign(aInLow);
            RemoveSign(cInLow);
            RemoveSign(dInLow);
            RemoveSign(aOutLow);
            RemoveSign(cOutLow);
            RemoveSign(dOutLow);
            if(aInLow > aInUp) { aInUp = aInLow;
            } else {aInLow = aInUp;}
            if(aOutLow > aOutUp) { aOutUp = aOutLow;
            } else { aOutLow = aOutUp;}
            if(cInLow > cInUp) { cInUp = cInLow;
            } else {cInLow = cInUp; };
            if(cOutLow > cOutUp ) { cOutUp = cOutLow;   
            } else { cOutLow = cOutUp; }
            // for confusio the logic changes here ...
            if(dInLow < dInUp) {dInLow = dInUp;
            } else { dInUp = dInLow;}
            if(dOutLow < dOutUp) { dOutLow = dOutUp;
            } else { dOutUp = dOutLow; }
        }
        // get the quadratic sums
        eInUp  = aInUp*aInUp + cInUp*cInUp + dInUp*dInUp;
        if(eInUp > 0) relativeErrorInUp->SetBinContent(b+1, 1.*TMath::Sqrt(eInUp));
        eOutUp = aOutUp*aOutUp + cOutUp*cOutUp + dOutUp*dOutUp;
        if(eOutUp > 0) relativeErrorOutUp->SetBinContent(b+1, 1.*TMath::Sqrt(eOutUp));

        eInLow = aInLow*aInLow + cInLow*cInLow + dInLow*dInLow;
        if(eInLow > 0) relativeErrorInLow->SetBinContent(b+1, -1.*TMath::Sqrt(eInLow));
        eOutLow = aOutLow*aOutLow + cOutLow*cOutLow + dOutLow*dOutLow;
        if(eOutLow > 0) relativeErrorOutLow->SetBinContent(b+1, -1.*TMath::Sqrt(eOutLow));

    printf("%i) \teInUp %.4f \t eInLow %.4f \t eOutUp %4.f \t eOutLow %4.f \n", b, eInUp, eInLow, eOutUp, eOutLow);
    }

    // project the estimated errors on the nominal ratio
   if(nominal) {
        Double_t* ax = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ay = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* axh = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* axl = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ayh = new Double_t[fBinsTrue->GetSize()-1];
        Double_t* ayl = new Double_t[fBinsTrue->GetSize()-1];
        Double_t lowErr(0.), upErr(0.);
        for(Int_t i(0); i < fBinsTrue->GetSize()-1; i++) {
            // add the in and out of plane errors in quadrature
            // take special care here: to propagate the assymetric error, we need to correlate the
            // InLow with OutUp (minimum value of ratio) and InUp with OutLow (maximum value of ratio)
            lowErr = relativeErrorInLow->GetBinContent(i+1)*relativeErrorInLow->GetBinContent(i+1)+relativeErrorOutUp->GetBinContent(1+i)*relativeErrorOutUp->GetBinContent(i+1);
            upErr = relativeErrorInUp->GetBinContent(i+1)*relativeErrorInUp->GetBinContent(i+1)+relativeErrorOutLow->GetBinContent(i+1)*relativeErrorOutLow->GetBinContent(i+1);
            // set the errors 
            ayl[i] = TMath::Sqrt(lowErr)*nominal->GetBinContent(i+1);
            ayh[i] = TMath::Sqrt(upErr)*nominal->GetBinContent(i+1);
            // get the bin width (which is the 'error' on x
            Double_t binWidth(nominal->GetBinWidth(i+1));
            axl[i] = binWidth/2.;
            axh[i] = binWidth/2.;
            // now get the coordinate for the point
            ax[i] = nominal->GetBinCenter(i+1);
            ay[i] = nominal->GetBinContent(i+1);
        }

        // save the nominal ratio
        TGraphAsymmErrors* nominalError(new TGraphAsymmErrors(fBinsTrue->GetSize()-1, ax, ay, axl, axh, ayl, ayh));
        shapeRatio = (TGraphAsymmErrors*)nominalError->Clone();
        nominalError->SetName("shape uncertainty");
        TCanvas* nominalCanvas(new TCanvas("nominalCanvas", "nominalCanvas"));
        nominalCanvas->Divide(2);
        nominalCanvas->cd(1);
        Style(nominal, kBlack);
        Style(nominalError, kCyan, kRatio);
        nominalError->SetLineColor(kCyan);
        nominalError->SetMarkerColor(kCyan);
        nominalError->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
        nominalError->GetYaxis()->SetRangeUser(.7, 2.2);
        nominalError->DrawClone("a2");
        nominal->DrawCopy("same E1");
        Style(AddLegend(gPad));
        Style(gPad, "GRID");
        Style(nominalCanvas);
        // save nominal v2 and systematic errors
        TGraphAsymmErrors* nominalV2Error(GetV2WithSystematicErrors(
                    nominalIn,
                    nominalOut,
                    fEventPlaneRes,
                    "v_{2} with shape uncertainty",
                    relativeErrorInUp,
                    relativeErrorInLow,
                    relativeErrorOutUp,
                    relativeErrorOutLow,
                    corr));
        shapeV2 = (TGraphAsymmErrors*)nominalV2Error->Clone();
        if(regularizationOnV2) shapeV2 = AddHistoErrorsToGraphErrors(shapeV2, relativeErrorRegularizationInUp);
        TGraphErrors* nominalV2(GetV2(nominalIn, nominalOut, fEventPlaneRes, "v_{2}"));
        // add in quadratufe (not very nice, rethink this because it may add additional weight to
        // the rms unfolded piece) the additional error
        nominalCanvas->cd(2);
        Style(nominalV2, kBlack);
        Style(nominalV2Error, kCyan, kV2);
        nominalV2Error->SetLineColor(kCyan);
        nominalV2Error->SetMarkerColor(kCyan);
        nominalV2Error->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
        nominalV2Error->DrawClone("a2");
        nominalV2->DrawClone("same E1");
        Style(AddLegend(gPad));
        Style(gPad, "GRID");
        Style(nominalCanvas);
        SavePadToPDF(nominalCanvas);
        nominalCanvas->Write();
    }
    TCanvas* relativeError(new TCanvas("relativeError"," relativeError"));
    relativeError->Divide(2);
    relativeError->cd(1);
    Style(gPad, "GRID");
    relativeErrorInUp->GetYaxis()->SetRangeUser(-1.5, 3.);
    Style(relativeErrorInUp, kBlue, kBar);
    Style(relativeErrorInLow, kGreen, kBar);
    relativeErrorInUp->DrawCopy("b");
    relativeErrorInLow->DrawCopy("same b");
    Style(relativeStatisticalErrorIn, kRed);
    relativeStatisticalErrorIn->DrawCopy("same");
    Style(AddLegend(gPad));
    relativeError->cd(2);
    Style(gPad, "GRID");
    Style(relativeErrorOutUp, kBlue, kBar);
    Style(relativeErrorOutLow, kGreen, kBar);
    relativeErrorOutUp->DrawCopy("b");
    if(relativeErrorOutLow) relativeErrorOutLow->DrawCopy("same b");
    if(relativeStatisticalErrorOut) Style(relativeStatisticalErrorOut, kRed);
    if(relativeStatisticalErrorOut) relativeStatisticalErrorOut->DrawCopy("same");
    Style(AddLegend(gPad));

    // write the buffered file to disk
    SavePadToPDF(relativeError);
    relativeError->Write();
    output->Write();
}
//_____________________________________________________________________________
void AliJetFlowTools::SystematicsWrapper(
            TArrayI* variationsIn,                  // variantions in plane
            TArrayI* variationsOut,                 // variantions out of plane
            TH1D*& relativeErrorInUp,               // pointer reference to minimum relative error histogram in plane
            TH1D*& relativeErrorInLow,              // pointer reference to maximum relative error histogram in plane
            TH1D*& relativeErrorOutUp,              // pointer reference to minimum relative error histogram out of plane
            TH1D*& relativeErrorOutLow,             // pointer reference to maximum relative error histogram out of plane
            TH1D*& relativeStatisticalErrorIn,      // relative systematic error on ratio
            TH1D*& relativeStatisticalErrorOut,     // relative systematic error on ratio
            TH1D*& nominal,                         // clone of the nominal ratio in / out of plane
            TH1D*& nominalIn,                       // clone of the nominal in plane yield
            TH1D*& nominalOut,                      // clone of the nominal out of plane yield
            Int_t columns,                          // divide the output canvasses in this many columns
            Float_t rangeLow,                       // lower pt range
            Float_t rangeUp,                        // upper pt range
            TFile* readMe,                          // input file name (created by this unfolding class)
            TString source,                         // source of the variation
            Bool_t RMS,                             // return RMS of distribution of variations as error
            Bool_t onRatio                          // use ratio or v2 directly for assessing systematic uncertainty
            ) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // wrapper function to call specific systematics function using function pointers
    void (AliJetFlowTools::*myFunction)(TArrayI*, TArrayI*, TH1D*&, TH1D*&, TH1D*&, TH1D*&, TH1D*&, TH1D*&, 
            TH1D*&, TH1D*&, TH1D*&, Int_t, Float_t, Float_t, TFile*, TString, Bool_t) const;

    printf(" >>> systematic wrapper called for %s <<< \n", source.Data());

    // assign function pointer
    myFunction = (onRatio) ? &AliJetFlowTools::DoIntermediateSystematics : &AliJetFlowTools::DoIntermediateSystematicsOnV2;

    // do the actual unfolding with the selected function
    return (this->*myFunction)( 
           variationsIn,             
           variationsOut,            
           relativeErrorInUp,          
           relativeErrorInLow,         
           relativeErrorOutUp,         
           relativeErrorOutLow,        
           relativeStatisticalErrorIn, 
           relativeStatisticalErrorOut,
           nominal,                    
           nominalIn,                  
           nominalOut,                 
           columns,                     
           rangeLow,                  
           rangeUp,                   
           readMe,                     
           source,                    
           RMS);
} 
//_____________________________________________________________________________
void AliJetFlowTools::DoIntermediateSystematics(
            TArrayI* variationsIn,                  // variantions in plane
            TArrayI* variationsOut,                 // variantions out of plane
            TH1D*& relativeErrorInUp,               // pointer reference to minimum relative error histogram in plane
            TH1D*& relativeErrorInLow,              // pointer reference to maximum relative error histogram in plane
            TH1D*& relativeErrorOutUp,              // pointer reference to minimum relative error histogram out of plane
            TH1D*& relativeErrorOutLow,             // pointer reference to maximum relative error histogram out of plane
            TH1D*& relativeStatisticalErrorIn,      // relative systematic error on ratio
            TH1D*& relativeStatisticalErrorOut,     // relative systematic error on ratio
            TH1D*& nominal,                         // clone of the nominal ratio in / out of plane
            TH1D*& nominalIn,                       // clone of the nominal in plane yield
            TH1D*& nominalOut,                      // clone of the nominal out of plane yield
            Int_t columns,                          // divide the output canvasses in this many columns
            Float_t rangeLow,                       // lower pt range
            Float_t rangeUp,                        // upper pt range
            TFile* readMe,                          // input file name (created by this unfolding class)
            TString source,                         // source of the variation
            Bool_t RMS                              // return RMS of distribution of variations as error
            ) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   // intermediate systematic check function. first index of supplied array is nominal value
   TList* listOfKeys((TList*)readMe->GetListOfKeys());
   if(!listOfKeys) {
       printf(" > Fatal error, couldn't retrieve list of keys. Input file might have been corrupted ! < \n");
       return;
   }
   // check input params
   if(variationsIn->GetSize() != variationsOut->GetSize()) {
       printf(" > DoSystematics: fatal error, input arrays have different sizes ! < \n ");
       return;
   }
   TDirectoryFile* defRootDirIn(dynamic_cast<TDirectoryFile*>(readMe->Get(listOfKeys->At(variationsIn->At(0))->GetName())));
   TDirectoryFile* defRootDirOut(dynamic_cast<TDirectoryFile*>(readMe->Get(listOfKeys->At(variationsOut->At(0))->GetName())));
   if(!(defRootDirIn && defRootDirOut)) {
       printf(" > DoSystematics: fatal error, couldn't retrieve nominal values ! < \n ");
       return;
   }
   TString defIn(defRootDirIn->GetName());
   TString defOut(defRootDirOut->GetName());

   // define lines to make the output prettier
   TLine* lineLow(new TLine(rangeLow, 0., rangeLow, 2.));
   TLine* lineUp(new TLine(rangeUp, 0., rangeUp, 2.));
   lineLow->SetLineColor(11);
   lineUp->SetLineColor(11);
   lineLow->SetLineWidth(3);
   lineUp->SetLineWidth(3);

   // define an output histogram with the maximum relative error from this systematic constribution
   // if the option RMS is set to false, sigma is not really a standard deviation but holds the maximum (or minimum) relative value that the data has
   // reached in this function call.
   // if the option RMS is set to true, sigma holds the RMS value (equal to sigma as the mean is zero for relative errors) of the distribution of variations
   // which should correspond to a 68% confidence level
   relativeErrorInUp = new TH1D(Form("relative error (up) from %s", source.Data()), Form("relative error (up) from %s", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
   relativeErrorInLow = new TH1D(Form("relative error (low) from  %s", source.Data()), Form("relative error (low) from %s", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
   relativeErrorOutUp = new TH1D(Form("relative error (up) from  %s", source.Data()), Form("relative error (up)  from %s", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
   relativeErrorOutLow = new TH1D(Form("relative error (low) from %s", source.Data()), Form("relative error (low) from %s", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
   for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
       if(!RMS) {
           relativeErrorInUp->SetBinContent(b+1, 1.);
           relativeErrorInUp->SetBinError(b+1, 0.);
           relativeErrorOutUp->SetBinContent(b+1, 1.);
           relativeErrorOutUp->SetBinError(b+1, .0);
           relativeErrorInLow->SetBinContent(b+1, 1.);
           relativeErrorInLow->SetBinError(b+1, 0.);
           relativeErrorOutLow->SetBinContent(b+1, 1.);
           relativeErrorOutLow->SetBinError(b+1, .0);
       } else if(RMS) {
           relativeErrorInUp->SetBinContent(b+1, 0.);
           relativeErrorInUp->SetBinError(b+1, 0.);
           relativeErrorOutUp->SetBinContent(b+1, 0.);
           relativeErrorOutUp->SetBinError(b+1, 0.);
           relativeErrorInLow->SetBinContent(b+1, 0.);
           relativeErrorInLow->SetBinError(b+1, 0.);
           relativeErrorOutLow->SetBinContent(b+1, 0.);
           relativeErrorOutLow->SetBinError(b+1, 0.);
       } 
   }
   Int_t relativeErrorInUpN[100] = {0};
   Int_t relativeErrorOutUpN[100] = {0};
   Int_t relativeErrorInLowN[100] = {0};
   Int_t relativeErrorOutLowN[100] = {0};
   
   // define an output histogram with the systematic error from this systematic constribution
   if(!relativeStatisticalErrorIn || !relativeStatisticalErrorOut) {
       relativeStatisticalErrorIn = new TH1D("relative statistical error, in plane", "relative statistical error, in plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
       relativeStatisticalErrorOut = new TH1D("relative statistical error, out of plane", "relative statistital error, out of plane", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
   }

   // prepare necessary canvasses
   TCanvas* canvasIn(new TCanvas(Form("SYST_%s_PearsonIn", source.Data()), Form("SYST_%s_PearsonIn", source.Data())));
   TCanvas* canvasOut(0x0);
   if(fDphiUnfolding) canvasOut = new TCanvas(Form("SYST_%s_PearsonOut", source.Data()), Form("SYST_%s_PearsonOut", source.Data()));
   TCanvas* canvasRatioMeasuredRefoldedIn(new TCanvas(Form("SYST_%s_RefoldedIn", source.Data()), Form("SYST_%s_RefoldedIn", source.Data())));
   TCanvas* canvasRatioMeasuredRefoldedOut(0x0);
   if(fDphiUnfolding) canvasRatioMeasuredRefoldedOut = new TCanvas(Form("SYST_%s_RefoldedOut", source.Data()), Form("SYST_%s_RefoldedOut", source.Data()));
   TCanvas* canvasSpectraIn(new TCanvas(Form("SYST_%s_SpectraIn", source.Data()), Form("SYST_%s_SpectraIn", source.Data()))); 
   TCanvas* canvasSpectraOut(0x0);
   if(fDphiUnfolding) canvasSpectraOut = new TCanvas(Form("SYST_%s_SpectraOut", source.Data()), Form("SYST_%s_SpectraOut", source.Data()));
   TCanvas* canvasRatio(0x0);
   if(fDphiUnfolding) canvasRatio = new TCanvas(Form("SYST_%s_Ratio", source.Data()), Form("SYST_%s_Ratio", source.Data()));
   TCanvas* canvasV2(0x0);
   if(fDphiUnfolding) canvasV2 = new TCanvas(Form("SYST_%s_V2", source.Data()), Form("SYST_%s_V2", source.Data()));
   TCanvas* canvasMISC(new TCanvas(Form("SYST_%s_MISC", source.Data()), Form("SYST_%s_MISC", source.Data())));
   TCanvas* canvasMasterIn(new TCanvas(Form("SYST_%s_defaultIn", source.Data()), Form("SYST_%s_defaultIn", source.Data())));
   TCanvas* canvasMasterOut(0x0);
   if(fDphiUnfolding) canvasMasterOut = new TCanvas(Form("SYST_%s_defaultOut", source.Data()), Form("SYST_%s_defaultOut", source.Data()));
   (fDphiUnfolding) ? canvasMISC->Divide(4, 2) : canvasMISC->Divide(4, 1);

   TCanvas* canvasProfiles(new TCanvas(Form("SYST_%s_canvasProfiles", source.Data()), Form("SYST_%s_canvasProfiles", source.Data())));
   canvasProfiles->Divide(2);
   TProfile* ratioProfile(new TProfile(Form("SYST_%s_ratioProfile", source.Data()), Form("SYST_%s_ratioProfile", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   TProfile* v2Profile(new TProfile(Form("SYST_%s_v2Profile", source.Data()), Form("SYST_%s_v2Profile", source.Data()),fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   // get an estimate of the number of outputs and find the default set

   Int_t rows = 1;
   columns = variationsIn->GetSize()-1;
   (TMath::Floor(variationsIn->GetSize()/(float)columns)+((variationsIn->GetSize()%columns)>0));
   canvasIn->Divide(columns, rows);
   if(canvasOut) canvasOut->Divide(columns, rows);
   canvasRatioMeasuredRefoldedIn->Divide(columns, rows);
   if(canvasRatioMeasuredRefoldedOut) canvasRatioMeasuredRefoldedOut->Divide(columns, rows);
   canvasSpectraIn->Divide(columns, rows);
   if(canvasSpectraOut) canvasSpectraOut->Divide(columns, rows);
   if(canvasRatio) canvasRatio->Divide(columns, rows);
   if(canvasV2) canvasV2->Divide(columns, rows);
   canvasMasterIn->Divide(columns, rows);
   if(canvasMasterOut) canvasMasterOut->Divide(columns, rows);

   // prepare a separate set of canvases to hold the nominal points
   TCanvas* canvasNominalIn(new TCanvas(Form("Nominal_%s_PearsonIn", source.Data()), Form("Nominal_%s_PearsonIn", source.Data())));
   TCanvas* canvasNominalOut(0x0);
   if(fDphiUnfolding) canvasNominalOut = new TCanvas(Form("Nominal_%s_PearsonOut", source.Data()), Form("Nominal_%s_PearsonOut", source.Data()));
   TCanvas* canvasNominalRatioMeasuredRefoldedIn(new TCanvas(Form("Nominal_%s_RefoldedIn", source.Data()), Form("Nominal_%s_RefoldedIn", source.Data())));
   TCanvas* canvasNominalRatioMeasuredRefoldedOut(0x0);
   if(fDphiUnfolding) canvasNominalRatioMeasuredRefoldedOut = new TCanvas(Form("Nominal_%s_RefoldedOut", source.Data()), Form("Nominal_%s_RefoldedOut", source.Data()));
   TCanvas* canvasNominalSpectraIn(new TCanvas(Form("Nominal_%s_SpectraIn", source.Data()), Form("Nominal_%s_SpectraIn", source.Data()))); 
   TCanvas* canvasNominalSpectraOut(0x0);
   if(fDphiUnfolding) canvasNominalSpectraOut = new TCanvas(Form("Nominal_%s_SpectraOut", source.Data()),  Form("Nominal_%s_SpectraOut", source.Data()));
   TCanvas* canvasNominalRatio(0x0);
   if(fDphiUnfolding) canvasNominalRatio = new TCanvas(Form("Nominal_%s_Ratio", source.Data()), Form("Nominal_%s_Ratio", source.Data()));
   TCanvas* canvasNominalV2(0x0);
   if(fDphiUnfolding) canvasNominalV2 = new TCanvas(Form("Nominal_%s_V2", source.Data()), Form("Nominal_%s_V2", source.Data()));
   TCanvas* canvasNominalMISC(new TCanvas(Form("Nominal_%s_MISC", source.Data()), Form("Nominal_%s_MISC", source.Data())));
   TCanvas* canvasNominalMasterIn(new TCanvas(Form("Nominal_%s_defaultIn", source.Data()), Form("Nominal_%s_defaultIn", source.Data())));
   TCanvas* canvasNominalMasterOut(0x0);
   if(fDphiUnfolding) canvasNominalMasterOut = new TCanvas(Form("Nominal_%s_defaultOut", source.Data()), Form("Nominal_%s_defaultOut", source.Data()));
   (fDphiUnfolding) ? canvasNominalMISC->Divide(4, 2) : canvasNominalMISC->Divide(4, 1);

   canvasNominalSpectraIn->Divide(2);
   if(canvasNominalSpectraOut) canvasNominalSpectraOut->Divide(2);

   canvasNominalMasterIn->Divide(2);
   if(canvasNominalMasterOut) canvasNominalMasterOut->Divide(2);

   // extract the default output 
   TH1D* defaultUnfoldedJetSpectrumIn(0x0);
   TH1D* defaultUnfoldedJetSpectrumOut(0x0);
   TDirectoryFile* defDirIn = (TDirectoryFile*)defRootDirIn->Get(Form("InPlane___%s", defIn.Data()));
   TDirectoryFile* defDirOut = (TDirectoryFile*)defRootDirOut->Get(Form("OutOfPlane___%s", defOut.Data()));
   if(defDirIn) defaultUnfoldedJetSpectrumIn = (TH1D*)defDirIn->Get(Form("UnfoldedSpectrum_in_%s", defIn.Data()));
   if(defDirOut) defaultUnfoldedJetSpectrumOut = (TH1D*)defDirOut->Get(Form("UnfoldedSpectrum_out_%s", defOut.Data()));
   printf(" > succesfully extracted default results < \n");

   // loop through the directories, only plot the graphs if the deconvolution converged
   TDirectoryFile* tempDirIn(0x0); 
   TDirectoryFile* tempDirOut(0x0);
   TDirectoryFile* tempIn(0x0);
   TDirectoryFile* tempOut(0x0);
   TH1D* unfoldedSpectrumInForRatio(0x0);
   TH1D* unfoldedSpectrumOutForRatio(0x0);
   for(Int_t i(0), j(-1); i < variationsIn->GetSize(); i++) {
       tempDirIn = (dynamic_cast<TDirectoryFile*>(readMe->Get(listOfKeys->At(variationsIn->At(i))->GetName())));
       tempDirOut = (dynamic_cast<TDirectoryFile*>(readMe->Get(listOfKeys->At(variationsOut->At(i))->GetName())));
       if(!(tempDirIn && tempDirOut)) {
           printf(" > DoSystematics: couldn't get a set of variations < \n");
           continue;
       }
       TString dirNameIn(tempDirIn->GetName());
       TString dirNameOut(tempDirOut->GetName());
       // try to read the in- and out of plane subdirs
       tempIn = (TDirectoryFile*)tempDirIn->Get(Form("InPlane___%s", dirNameIn.Data()));
       tempOut = (TDirectoryFile*)tempDirOut->Get(Form("OutOfPlane___%s", dirNameOut.Data()));
       j++;
       if(tempIn) { 
           // to see if the unfolding converged try to extract the pearson coefficients
           TH2D* pIn((TH2D*)tempIn->Get(Form("PearsonCoefficients_in_%s", dirNameIn.Data())));
           if(pIn) {
               printf(" - %s in plane converged \n", dirNameIn.Data());
               canvasIn->cd(j);
               if(i==0) canvasNominalIn->cd(j);
               Style(gPad, "PEARSON");
               pIn->DrawCopy("colz");
               TGraphErrors* rIn((TGraphErrors*)tempIn->Get(Form("RatioRefoldedMeasured_%s", dirNameIn.Data())));
               if(rIn) {
                   Style(rIn);
                   printf(" > found RatioRefoldedMeasured < \n");
                   canvasRatioMeasuredRefoldedIn->cd(j);
                   if(i==0) canvasNominalRatioMeasuredRefoldedIn->cd(j);
                   Style(gPad, "GRID");
                   rIn->SetFillColor(kRed);
                   rIn->Draw("ap");
               }
               TH1D* dvector((TH1D*)tempIn->Get("dVector"));
               TH1D* avalue((TH1D*)tempIn->Get("SingularValuesOfAC"));
               TH2D* rm((TH2D*)tempIn->Get(Form("ResponseMatrixIn_%s", dirNameIn.Data())));
               TH1D* eff((TH1D*)tempIn->Get(Form("KinematicEfficiencyIn_%s", dirNameIn.Data())));
               if(dvector && avalue && rm && eff) {
                   Style(dvector);
                   Style(avalue);
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(1);
                   if(i==0) canvasNominalMISC->cd(1);
                   Style(gPad, "SPECTRUM");
                   dvector->DrawCopy();
                   canvasMISC->cd(2);
                   if(i==0) canvasNominalMISC->cd(2);
                   Style(gPad, "SPECTRUM");
                   avalue->DrawCopy();
                   canvasMISC->cd(3);
                   if(i==0) canvasNominalMISC->cd(3);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(4);
                   if(i==0) canvasNominalMISC->cd(4);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(3);
                   if(i==0) canvasNominalMISC->cd(3);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(4);
                   if(i==0) canvasNominalMISC->cd(4);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               }
           }
           TH1D* inputSpectrum((TH1D*)tempIn->Get(Form("InputSpectrum_in_%s", dirNameIn.Data())));
           TH1D* unfoldedSpectrum((TH1D*)tempIn->Get(Form("UnfoldedSpectrum_in_%s", dirNameIn.Data())));
           unfoldedSpectrumInForRatio = ProtectHeap(unfoldedSpectrum, kFALSE, TString("ForRatio"));
           TH1D* refoldedSpectrum((TH1D*)tempIn->Get(Form("RefoldedSpectrum_in_%s", dirNameIn.Data())));
           if(inputSpectrum && unfoldedSpectrum && refoldedSpectrum) {
               if(defaultUnfoldedJetSpectrumIn) {
                   Style(defaultUnfoldedJetSpectrumIn, kBlue, kUnfoldedSpectrum);
                   TH1D* temp((TH1D*)defaultUnfoldedJetSpectrumIn->Clone(Form("defaultUnfoldedJetSpectrumIn_%s", dirNameIn.Data())));
                   temp->Divide(unfoldedSpectrum);
                   // get the absolute relative error
                   for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
                       if(!RMS) {       // save the maximum deviation that a variation can cause
                           // the variation is HIGHER than the nominal point, so the bar goes UP
                           if( temp->GetBinContent(b+1) < 1 && temp->GetBinContent(b+1) < relativeErrorInUp->GetBinContent(b+1)) {
                               relativeErrorInUp->SetBinContent(b+1, temp->GetBinContent(b+1));
                               relativeErrorInUp->SetBinError(b+1, temp->GetBinError(b+1));
                           }
                           // the variation is LOWER than the nominal point, so the bar goes DOWN
                           else if(temp->GetBinContent(b+1) > 1 && temp->GetBinContent(b+1) > relativeErrorInLow->GetBinContent(b+1)) {
                               relativeErrorInLow->SetBinContent(b+1, temp->GetBinContent(b+1));
                               relativeErrorInLow->SetBinError(b+1, temp->GetBinError(b+1));
                           }
                       } else if (RMS && !fSymmRMS) { // save info necessary for evaluating the RMS of a distribution of variations
                           printf(" oops shouldnt be here \n " );
                           if(temp->GetBinContent(b+1) < 1) {
                               relativeErrorInUp->SetBinContent(b+1, relativeErrorInUp->GetBinContent(b+1)+TMath::Power(1.-temp->GetBinContent(b+1), 2));
                               relativeErrorInUpN[b]++;
                           }
                           // the variation is LOWER than the nominal point, so the bar goes DOWN
                           else if(temp->GetBinContent(b+1) > 1) {
                               relativeErrorInLow->SetBinContent(b+1, relativeErrorInLow->GetBinContent(b+1)+TMath::Power(1.-temp->GetBinContent(b+1), 2));
                               relativeErrorInLowN[b]++;
                           }
                       } else if (fSymmRMS) {
                           // save symmetric sum of square to get a symmetric rms
                           relativeErrorInUp->SetBinContent(b+1, relativeErrorInUp->GetBinContent(b+1)+TMath::Power(temp->GetBinContent(b+1)-1., 2));
                           relativeErrorInUpN[b]++;
                       }
                       if(temp->GetBinError(b+1) > 0) relativeStatisticalErrorIn->SetBinContent(b+1, temp->GetBinError(b+1)/temp->GetBinContent(b+1));
                   }
                   temp->SetTitle(Form("[%s] / [%s]", defIn.Data(), dirNameIn.Data()));
                   temp->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                   temp->GetYaxis()->SetTitle("ratio");
                   canvasMasterIn->cd(j);
                   temp->GetYaxis()->SetRangeUser(0., 2);
                   Style(gPad, "GRID");
                   temp->DrawCopy();
                   canvasNominalMasterIn->cd(1);
                   Style(gPad, "GRID");
                   if(i > 0 ) {
                       TH1D* tempSyst((TH1D*)temp->Clone(Form("%s_syst", temp->GetName())));
                       tempSyst->SetTitle(Form("[%s] / [%s]", defIn.Data(), dirNameIn.Data()));
                       Style(tempSyst, (EColor)(i+2));
                       if(i==1) tempSyst->DrawCopy();
                       else tempSyst->DrawCopy("same");
                   }
               }
               TH1F* fitStatus((TH1F*)tempIn->Get(Form("fitStatus_%s_in", dirNameIn.Data())));
               canvasSpectraIn->cd(j);
               if(i==0) canvasNominalSpectraIn->cd(1);
               Style(gPad);
               Style(unfoldedSpectrum, kRed, kUnfoldedSpectrum);
               unfoldedSpectrum->DrawCopy();
               Style(inputSpectrum, kGreen, kMeasuredSpectrum);
               inputSpectrum->DrawCopy("same");
               Style(refoldedSpectrum, kBlue, kFoldedSpectrum);
               refoldedSpectrum->DrawCopy("same");
               TLegend* l(AddLegend(gPad));
               Style(l);
               if(fitStatus && fitStatus->GetNbinsX() == 4) { // only available in chi2 fit
                   Float_t chi(fitStatus->GetBinContent(1));
                   Float_t pen(fitStatus->GetBinContent(2));
                   Int_t dof((int)fitStatus->GetBinContent(3));
                   Float_t beta(fitStatus->GetBinContent(4));
                   l->AddEntry((TObject*)0, Form("#chi %.2f \tP %.2f \tDOF %i, #beta %.2f", chi, pen, dof, beta), "");
               } else if (fitStatus) { // only available in SVD fit
                   Int_t reg((int)fitStatus->GetBinContent(1));
                   l->AddEntry((TObject*)0, Form("REG %i", reg), "");
               }
               canvasNominalSpectraIn->cd(2);
               TH1D* tempSyst((TH1D*)unfoldedSpectrum->Clone(Form("%s_syst", unfoldedSpectrum->GetName())));
               tempSyst->SetTitle(Form("[%s]", dirNameIn.Data()));
               Style(tempSyst, (EColor)(i+2));
               Style(gPad, "SPECTRUM");
               if(i==0) tempSyst->DrawCopy();
               else tempSyst->DrawCopy("same");
           }
       }
       if(tempOut) {
           TH2D* pOut((TH2D*)tempOut->Get(Form("PearsonCoefficients_out_%s", dirNameOut.Data())));
           if(pOut) {
               printf(" - %s out of plane converged \n", dirNameOut.Data());
               canvasOut->cd(j);
               if(i==0) canvasNominalOut->cd(j);
               Style(gPad, "PEARSON");
               pOut->DrawCopy("colz");
               TGraphErrors* rOut((TGraphErrors*)tempOut->Get(Form("RatioRefoldedMeasured_%s", dirNameOut.Data())));
               if(rOut) {
                   Style(rOut);
                   printf(" > found RatioRefoldedMeasured < \n");
                   canvasRatioMeasuredRefoldedOut->cd(j);
                   if(i==0) canvasNominalRatioMeasuredRefoldedOut->cd(j);
                   Style(gPad, "GRID");
                   rOut->SetFillColor(kRed);
                   rOut->Draw("ap");
               }
               TH1D* dvector((TH1D*)tempOut->Get("dVector"));
               TH1D* avalue((TH1D*)tempOut->Get("SingularValuesOfAC"));
               TH2D* rm((TH2D*)tempOut->Get(Form("ResponseMatrixOut_%s", dirNameOut.Data())));
               TH1D* eff((TH1D*)tempOut->Get(Form("KinematicEfficiencyOut_%s", dirNameOut.Data())));
               if(dvector && avalue && rm && eff) {
                   Style(dvector);
                   Style(avalue);
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(5);
                   if(i==0) canvasNominalMISC->cd(5);
                   Style(gPad, "SPECTRUM");
                   dvector->DrawCopy();
                   canvasMISC->cd(6);
                   if(i==0) canvasNominalMISC->cd(6);
                   Style(gPad, "SPECTRUM");
                   avalue->DrawCopy();
                   canvasMISC->cd(7);
                   if(i==0) canvasNominalMISC->cd(7);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(8);
                   if(i==0) canvasNominalMISC->cd(8);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(7);
                   if(i==0) canvasNominalMISC->cd(7);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(8);
                   if(i==0) canvasNominalMISC->cd(8);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               }
           }
           TH1D* inputSpectrum((TH1D*)tempOut->Get(Form("InputSpectrum_out_%s", dirNameOut.Data())));
           TH1D* unfoldedSpectrum((TH1D*)tempOut->Get(Form("UnfoldedSpectrum_out_%s", dirNameOut.Data())));
           unfoldedSpectrumOutForRatio = ProtectHeap(unfoldedSpectrum, kFALSE, TString("ForRatio"));
           TH1D* refoldedSpectrum((TH1D*)tempOut->Get(Form("RefoldedSpectrum_out_%s", dirNameOut.Data())));
           if(inputSpectrum && unfoldedSpectrum && refoldedSpectrum) {
               if(defaultUnfoldedJetSpectrumOut) {
                   Style(defaultUnfoldedJetSpectrumOut, kBlue, kUnfoldedSpectrum);
                   TH1D* temp((TH1D*)defaultUnfoldedJetSpectrumOut->Clone(Form("defaultUnfoldedJetSpectrumOut_%s", dirNameOut.Data())));
                   temp->Divide(unfoldedSpectrum);
                   // get the absolute relative error 
                   for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
                       if(!RMS) {
                           // check if the error is larger than the current maximum
                           if(temp->GetBinContent(b+1) < 1 && temp->GetBinContent(b+1) < relativeErrorOutUp->GetBinContent(b+1)) {
                               relativeErrorOutUp->SetBinContent(b+1, temp->GetBinContent(b+1));
                               relativeErrorOutUp->SetBinError(b+1, temp->GetBinError(b+1));
                           }
                           // check if the error is smaller than the current minimum
                           else if(temp->GetBinContent(b+1) > 1 && temp->GetBinContent(b+1) > relativeErrorOutLow->GetBinContent(b+1)) {
                               relativeErrorOutLow->SetBinContent(b+1, temp->GetBinContent(b+1));
                               relativeErrorOutLow->SetBinError(b+1, temp->GetBinError(b+1));
                           }
                       } else if (RMS && !fSymmRMS) {
                           printf(" OOps \n ");
                           if(temp->GetBinContent(b+1) < 1) {
                               relativeErrorOutUp->SetBinContent(b+1, relativeErrorOutUp->GetBinContent(b+1)+TMath::Power(1.-temp->GetBinContent(b+1), 2));
                               relativeErrorOutUpN[b]++;
                           }
                           else if(temp->GetBinContent(b+1) > 1) {
                               relativeErrorOutLow->SetBinContent(b+1, relativeErrorOutLow->GetBinContent(b+1)+TMath::Power(1.-temp->GetBinContent(b+1), 2));
                               relativeErrorOutLowN[b]++;
                           }
                       } else if (fSymmRMS) {
                           // save symmetric rms value
                           relativeErrorOutUp->SetBinContent(b+1, relativeErrorOutUp->GetBinContent(b+1)+TMath::Power(temp->GetBinContent(b+1)-1., 2));
                           relativeErrorOutUpN[b]++;
                       }
                       if(temp->GetBinError(b+1) > 0) relativeStatisticalErrorOut->SetBinContent(b+1, temp->GetBinError(b+1)/temp->GetBinContent(b+1));
                   }
                   temp->SetTitle(Form("[%s] / [%s]", defOut.Data(), dirNameOut.Data()));
                   temp->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                   temp->GetYaxis()->SetTitle("ratio");
                   canvasMasterOut->cd(j);
                   temp->GetYaxis()->SetRangeUser(0., 2);
                   Style(gPad, "GRID");
                   temp->DrawCopy();
                   canvasNominalMasterOut->cd(1);
                   Style(gPad, "GRID");
                   if(i > 0 ) {
                       TH1D* tempSyst((TH1D*)temp->Clone(Form("%s_syst", temp->GetName())));
                       tempSyst->SetTitle(Form("[%s] / [%s]", defOut.Data(), dirNameOut.Data()));
                       Style(tempSyst, (EColor)(i+2));
                       if(i==1) tempSyst->DrawCopy();
                       else tempSyst->DrawCopy("same");
                   }
               }
               TH1F* fitStatus((TH1F*)tempOut->Get(Form("fitStatus_%s_out", dirNameOut.Data())));
               canvasSpectraOut->cd(j);
               if(i==0) canvasNominalSpectraOut->cd(1);
               Style(gPad);
               Style(unfoldedSpectrum, kRed, kUnfoldedSpectrum);
               unfoldedSpectrum->DrawCopy();
               Style(inputSpectrum, kGreen, kMeasuredSpectrum);
               inputSpectrum->DrawCopy("same");
               Style(refoldedSpectrum, kBlue, kFoldedSpectrum);
               refoldedSpectrum->DrawCopy("same");
               TLegend* l(AddLegend(gPad));
               Style(l);
               if(fitStatus && fitStatus->GetNbinsX() == 4) { // only available in chi2 fit
                   Float_t chi(fitStatus->GetBinContent(1));
                   Float_t pen(fitStatus->GetBinContent(2));
                   Int_t dof((int)fitStatus->GetBinContent(3));
                   Float_t beta(fitStatus->GetBinContent(4));
                   l->AddEntry((TObject*)0, Form("#chi %.2f \tP %.2f \tDOF %i, #beta %.2f", chi, pen, dof, beta), "");
               } else if (fitStatus) { // only available in SVD fit
                   Int_t reg((int)fitStatus->GetBinContent(1));
                   l->AddEntry((TObject*)0, Form("REG %i", reg), "");
               }
               canvasNominalSpectraOut->cd(2);
               TH1D* tempSyst((TH1D*)unfoldedSpectrum->Clone(Form("%s_syst", unfoldedSpectrum->GetName())));
               tempSyst->SetTitle(Form("[%s]", dirNameOut.Data()));
               Style(tempSyst, (EColor)(i+2));
               Style(gPad, "SPECTRUM");
               if(i==0) tempSyst->DrawCopy();
               else tempSyst->DrawCopy("same");
           }
       }
       if(canvasRatio && canvasV2) {
           canvasRatio->cd(j);
           if(i==0) {
               canvasNominalRatio->cd(j);
               if(nominal && nominalIn && nominalOut) {
                   // if a nominal ratio is requested, delete the placeholder and update the nominal point
                   delete nominal;
                   delete nominalIn;
                   delete nominalOut;
                   nominalIn =  (TH1D*)unfoldedSpectrumInForRatio->Clone("in plane jet yield");
                   nominalOut =  (TH1D*)unfoldedSpectrumOutForRatio->Clone("out of plane jet yield");
                   nominal = (TH1D*)unfoldedSpectrumInForRatio->Clone("ratio in plane out of plane");
                   nominal->Divide(nominal, unfoldedSpectrumOutForRatio);       // standard root divide for uncorrelated histos
               }
           }
           TGraphErrors* ratioYield(GetRatio(unfoldedSpectrumInForRatio, unfoldedSpectrumOutForRatio, TString(Form("ratio [in=%s, out=%s]", dirNameIn.Data(), dirNameOut.Data()))));
           Double_t _tempx(0), _tempy(0);
           if(ratioYield) {
               Style(ratioYield);
               for(Int_t b(0); b < fBinsTrue->GetSize(); b++) {
                   ratioYield->GetPoint(b, _tempx, _tempy);
                   ratioProfile->Fill(_tempx, _tempy);
               }
               ratioProfile->GetYaxis()->SetRangeUser(-0., 2.);
               ratioProfile->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioYield->GetYaxis()->SetRangeUser(-0., 2.);
               ratioYield->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioYield->SetFillColor(kRed);
               ratioYield->Draw("ap");
           }
           canvasV2->cd(j);
           if(i==0) canvasNominalV2->cd(j);
           TGraphErrors* ratioV2(GetV2(unfoldedSpectrumInForRatio,unfoldedSpectrumOutForRatio, fEventPlaneRes, TString(Form("v_{2} [in=%s, out=%s]", dirNameIn.Data(), dirNameOut.Data()))));
           if(ratioV2) {
               Style(ratioV2);
               for(Int_t b(0); b < fBinsTrue->GetSize(); b++) {
                   ratioV2->GetPoint(b, _tempx, _tempy);
                   v2Profile->Fill(_tempx, _tempy);

               }
               v2Profile->GetYaxis()->SetRangeUser(-0., 2.);
               v2Profile->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioV2->GetYaxis()->SetRangeUser(-.25, .75);
               ratioV2->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioV2->SetFillColor(kRed);
               ratioV2->Draw("ap");
           }
       }
       delete unfoldedSpectrumInForRatio;
       delete unfoldedSpectrumOutForRatio;
   }
   // save the canvasses
   canvasProfiles->cd(1);
   Style(ratioProfile);
   ratioProfile->DrawCopy();
   canvasProfiles->cd(2);
   Style(v2Profile);
   v2Profile->DrawCopy();
   SavePadToPDF(canvasProfiles);
   canvasProfiles->Write(); 
   SavePadToPDF(canvasIn);
   canvasIn->Write();
   if(canvasOut) {
       SavePadToPDF(canvasOut);
       canvasOut->Write();
   }
   SavePadToPDF(canvasRatioMeasuredRefoldedIn);
   canvasRatioMeasuredRefoldedIn->Write();
   if(canvasRatioMeasuredRefoldedOut) {
       SavePadToPDF(canvasRatioMeasuredRefoldedOut);
       canvasRatioMeasuredRefoldedOut->Write();
   }
   SavePadToPDF(canvasSpectraIn);
   canvasSpectraIn->Write();
   if(canvasSpectraOut) {
       SavePadToPDF(canvasSpectraOut);
       canvasSpectraOut->Write();
       SavePadToPDF(canvasRatio);
       canvasRatio->Write();
       SavePadToPDF(canvasV2);
       canvasV2->Write();
   }
   SavePadToPDF(canvasMasterIn);
   canvasMasterIn->Write();
   if(canvasMasterOut) {
       SavePadToPDF(canvasMasterOut);
       canvasMasterOut->Write();
   }
   SavePadToPDF(canvasMISC);
   canvasMISC->Write();
   // save the nomial canvasses
   SavePadToPDF(canvasNominalIn);
   canvasNominalIn->Write();
   if(canvasNominalOut) {
       SavePadToPDF(canvasNominalOut);
       canvasNominalOut->Write();
   }
   SavePadToPDF(canvasNominalRatioMeasuredRefoldedIn);
   canvasNominalRatioMeasuredRefoldedIn->Write();
   if(canvasNominalRatioMeasuredRefoldedOut) {
       SavePadToPDF(canvasNominalRatioMeasuredRefoldedOut);
       canvasNominalRatioMeasuredRefoldedOut->Write();
   }
   canvasNominalSpectraIn->cd(2);
   Style(AddLegend(gPad)); 
   SavePadToPDF(canvasNominalSpectraIn);
   canvasNominalSpectraIn->Write();
   if(canvasNominalSpectraOut) {
       canvasNominalSpectraOut->cd(2);
       Style(AddLegend(gPad));
       SavePadToPDF(canvasNominalSpectraOut);
       canvasNominalSpectraOut->Write();
       SavePadToPDF(canvasNominalRatio);
       canvasNominalRatio->Write();
       SavePadToPDF(canvasNominalV2);
       canvasNominalV2->Write();
   }
   canvasNominalMasterIn->cd(1);
   Style(AddLegend(gPad));
   lineUp->DrawClone("same");
   lineLow->DrawClone("same");
   SavePadToPDF(canvasNominalMasterIn);
   canvasNominalMasterIn->Write();
   if(canvasNominalMasterOut) {
       canvasNominalMasterOut->cd(1);
       Style(AddLegend(gPad));
       lineUp->DrawClone("same");
       lineLow->DrawClone("same");
       SavePadToPDF(canvasNominalMasterOut);
       canvasNominalMasterOut->Write();
   }
   SavePadToPDF(canvasNominalMISC);
   canvasNominalMISC->Write();

   // save the relative errors
   for(Int_t b(0); b < fBinsTrue->GetSize()-1; b++) {
       // to arrive at a min and max from here, combine in up and out low
       if(!RMS) {
           relativeErrorInUp->SetBinContent(b+1, -1.*(relativeErrorInUp->GetBinContent(b+1)-1));
//           relativeErrorInUp->SetBinError(b+1, 0.);
           relativeErrorOutUp->SetBinContent(b+1, -1.*(relativeErrorOutUp->GetBinContent(b+1)-1));
//           relativeErrorOutUp->SetBinError(b+1, .0);
           relativeErrorInLow->SetBinContent(b+1, -1.*(relativeErrorInLow->GetBinContent(b+1)-1));
//           relativeErrorInLow->SetBinError(b+1, 0.);
           relativeErrorOutLow->SetBinContent(b+1, -1.*(relativeErrorOutLow->GetBinContent(b+1)-1));
//           relativeErrorOutLow->SetBinError(b+1, .0);
       } else if (RMS) {
           // these guys are already stored as percentages, so no need to remove the offset of 1
           // RMS is defined as sqrt(sum(squared))/N
           // min is defined as negative, max is defined as positive
           if(!fSymmRMS) {
               if(relativeErrorInUpN[b] < 1) relativeErrorInUpN[b] = 1;
               if(relativeErrorInLowN[b] < 1) relativeErrorInLowN[b] = 1;
               if(relativeErrorOutUpN[b] < 1) relativeErrorOutUpN[b] = 1;
               if(relativeErrorOutLowN[b] < 1) relativeErrorOutLowN[b] = 1;
               relativeErrorInUp->SetBinContent(b+1, TMath::Sqrt(relativeErrorInUp->GetBinContent(b+1)/relativeErrorInUpN[b]));
//               relativeErrorInUp->SetBinError(b+1, 0.);
               relativeErrorOutUp->SetBinContent(b+1, TMath::Sqrt(relativeErrorOutUp->GetBinContent(b+1)/relativeErrorOutUpN[b]));
//               relativeErrorOutUp->SetBinError(b+1, .0);
               relativeErrorInLow->SetBinContent(b+1, -1.*TMath::Sqrt(relativeErrorInLow->GetBinContent(b+1)/relativeErrorInLowN[b]));
//               relativeErrorInLow->SetBinError(b+1, 0.);
               relativeErrorOutLow->SetBinContent(b+1, -1.*TMath::Sqrt(relativeErrorOutLow->GetBinContent(b+1)/relativeErrorOutLowN[b]));
//               relativeErrorOutLow->SetBinError(b+1, .0);
           } else if (fSymmRMS) {
               if(relativeErrorInUpN[b] < 1) relativeErrorInUpN[b] = 1;
               if(relativeErrorOutUpN[b] < 1) relativeErrorOutUpN[b] = 1;
               relativeErrorInUp->SetBinContent(b+1, TMath::Sqrt(relativeErrorInUp->GetBinContent(b+1)/relativeErrorInUpN[b]));
               relativeErrorOutUp->SetBinContent(b+1, TMath::Sqrt(relativeErrorOutUp->GetBinContent(b+1)/relativeErrorOutUpN[b]));
           }
       }
   }
   relativeErrorInUp->SetYTitle("relative uncertainty");
   relativeErrorOutUp->SetYTitle("relative uncertainty");
   relativeErrorInLow->SetYTitle("relative uncertainty");
   relativeErrorOutLow->SetYTitle("relative uncertainty");
   relativeErrorInUp->GetYaxis()->SetRangeUser(-1.5, 3.);
   relativeErrorInLow->GetYaxis()->SetRangeUser(-1.5, 3.);
   relativeErrorOutUp->GetYaxis()->SetRangeUser(-1.5, 3.);
   relativeErrorOutLow->GetYaxis()->SetRangeUser(-1.5, 3.);

   canvasNominalMasterIn->cd(2);
   Style(gPad, "GRID");
   Style(relativeErrorInUp, kBlue, kBar);
   Style(relativeErrorInLow, kGreen, kBar);
   relativeErrorInUp->DrawCopy("b");
   relativeErrorInLow->DrawCopy("same b");
   Style(AddLegend(gPad));
   SavePadToPDF(canvasNominalMasterIn);
   canvasNominalMasterIn->Write();
   canvasNominalMasterOut->cd(2);
   Style(gPad, "GRID");
   Style(relativeErrorOutUp, kBlue, kBar);
   Style(relativeErrorOutLow, kGreen, kBar);
   relativeErrorOutUp->DrawCopy("b");
   relativeErrorOutLow->DrawCopy("same b");
   Style(AddLegend(gPad));
   SavePadToPDF(canvasNominalMasterOut);
   canvasNominalMasterOut->Write();
}
//_____________________________________________________________________________
void AliJetFlowTools::DoIntermediateSystematicsOnV2(
            TArrayI* variationsIn,                  // variantions in plane
            TArrayI* variationsOut,                 // variantions out of plane
            TH1D*& relativeErrorInUp,               // will store systematic error on v2 
            TH1D*& relativeErrorInLow,              // remains NULL
            TH1D*& relativeErrorOutUp,              // remains NULL
            TH1D*& relativeErrorOutLow,             // remains NULL
            TH1D*& relativeStatisticalErrorIn,      // remains NULL
            TH1D*& relativeStatisticalErrorOut,     // remains NULL
            TH1D*& nominal,                         // nominal v2
            TH1D*& nominalIn,                       // remains NULL
            TH1D*& nominalOut,                      // remains NULL
            Int_t columns,                          // trivial
            Float_t rangeLow,                       // trivial
            Float_t rangeUp,                        // trivial
            TFile* readMe,                          // trivial
            TString source,                         // trivial
            Bool_t RMS                              // NOT trivial
            ) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // 0) intermediate systematic check function. first make sure a placeholder is ready
    if(relativeErrorInUp) delete relativeErrorInUp;
    relativeErrorInUp = new TH1D(Form("absolute_systematic_uncertainty_%s", source.Data()), Form("absolute_systematic_uncertainty_%s", source.Data()), fBinsTrue->GetSize()-1, fBinsTrue->GetArray());

    // 1) then loop over the variations
    for(Int_t i(0); i < variationsIn->GetSize(); i++) {
        printf(" variation %i of %i \n", i, variationsIn->GetSize());
        // the frist variation will be stored as nominal, others are looped over
        Int_t iIn[] = {variationsIn->At(i), variationsIn->At(i)};
        Int_t iOut[] = {variationsOut->At(i), variationsOut->At(i)};

        // call the functions and set the relevant pointer references
        TH1D* dud(0x0);
        DoIntermediateSystematics(      // i guess this is the point where i retreive the variation
                new TArrayI(2, iIn),    // except for i is one, that's the nominal measurement
                new TArrayI(2, iOut),
                dud, dud, dud, dud, dud, dud, dud,
                nominalIn,
                nominalOut,
                1, 
                fBinsTrue->At(0), 
                fBinsTrue->At(fBinsTrue->GetSize()-1),
                readMe,
                Form("error_on_v2_variation_%i", i));

        printf(" nominal in %i binsX \n", nominalIn->GetNbinsX());
        printf(" nominal out %i binsX \n", nominalOut->GetNbinsX());

        // 2) get a histo with the v2
        TH1D* v2(GetV2Histo(nominalIn, nominalOut, fEventPlaneRes, "nominal v_{2}"));

        printf(" v2 in %i binsX \n", v2->GetNbinsX());
        // 3) bounce if nominal measurement ...
        if(i==0) {
            nominal = v2; 
        } else {  // 4) calculate the uncertainty. for now only RMS style, easier and the first check
            for(Int_t k(0); k < nominal->GetNbinsX(); k++) {
                // 5) store the sum of squares of difference (so orig + square of new point) differentially
                relativeErrorInUp->SetBinContent(k+1, relativeErrorInUp->GetBinContent(k+1) + TMath::Power(v2->GetBinContent(k+1)-nominal->GetBinContent(k+1), 2));
                printf(" nominal bin %i is %.4f \n", k+1, nominal->GetBinContent(k+1));
                printf(" variati bin %i is %.4f \n", k+1, v2->GetBinContent(k+1));
                printf(" sum of squared difference is %.4f \n", relativeErrorInUp->GetBinContent(k+1));
            }
            delete v2;
        }
    }
    // 6) looped over all the variations, now divide by N and take the square for the RMS 
    //    to make life fun, the function that will add all uncertainties uses bin ERRORS, not content, 
    //    so we'll set the error here and set the relative error as bin content
    for(Int_t k(0); k < nominal->GetNbinsX(); k++) {
        relativeErrorInUp->SetBinError(k+1, TMath::Sqrt(relativeErrorInUp->GetBinContent(k+1)/((double)variationsIn->GetSize())));
        printf(" total absolute error in bin %i is %.4f \n", k+1, relativeErrorInUp->GetBinError(k+1));
        relativeErrorInUp->SetBinContent(k+1, relativeErrorInUp->GetBinError(k+1)/nominal->GetBinContent(k+1));
    }
    if(fConstantUE) {
        Int_t ctr(0);
        Double_t av(0), avct(0);
        for(Int_t k(3); k < 8; k++) {
            av+=relativeErrorInUp->GetBinError(k+1);
            avct+=relativeErrorInUp->GetBinContent(k+1);
            ctr++;
        }
        av/=((double)ctr);
        avct/=((double)ctr);
        for(Int_t k(8); k < nominal->GetNbinsX(); k++) {
            relativeErrorInUp->SetBinError(k+1, av);
            relativeErrorInUp->SetBinContent(k+1, avct);
        }
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::PostProcess(TString def, Int_t columns, Float_t rangeLow, Float_t rangeUp, TString in, TString out) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   // go through the output file and perform post processing routines
   // can either be performed in one go with the unfolding, or at a later stage
   if(fOutputFile && !fOutputFile->IsZombie()) fOutputFile->Close();
   TFile readMe(in.Data(), "READ");     // open file read-only
   if(readMe.IsZombie()) {
       printf(" > Fatal error, couldn't read %s for post processing ! < \n", in.Data());
       return;
   }
   printf("\n\n\n\t\t POSTPROCESSING \n > Recovered the following file structure : \n <");
   readMe.ls();
   TList* listOfKeys((TList*)readMe.GetListOfKeys());
   if(!listOfKeys) {
       printf(" > Fatal error, couldn't retrieve list of keys. Input file might have been corrupted ! < \n");
       return;
   }
   // prepare necessary canvasses
   TCanvas* canvasIn(new TCanvas("PearsonIn", "PearsonIn"));
   TCanvas* canvasOut(0x0);
   if(fDphiUnfolding) canvasOut = new TCanvas("PearsonOut", "PearsonOut");
   TCanvas* canvasRatioMeasuredRefoldedIn(new TCanvas("RefoldedIn", "RefoldedIn"));
   TCanvas* canvasRatioMeasuredRefoldedOut(0x0);
   if(fDphiUnfolding) canvasRatioMeasuredRefoldedOut = new TCanvas("RefoldedOut", "RefoldedOut");
   TCanvas* canvasSpectraIn(new TCanvas("SpectraIn", "SpectraIn")); 
   TCanvas* canvasSpectraOut(0x0);
   if(fDphiUnfolding) canvasSpectraOut = new TCanvas("SpectraOut", "SpectraOut");
   TCanvas* canvasRatio(0x0);
   if(fDphiUnfolding) canvasRatio = new TCanvas("Ratio", "Ratio");
   TCanvas* canvasV2(0x0);
   if(fDphiUnfolding) canvasV2 = new TCanvas("V2", "V2");
   TCanvas* canvasMISC(new TCanvas("MISC", "MISC"));
   TCanvas* canvasMasterIn(new TCanvas("defaultIn", "defaultIn"));
   TCanvas* canvasMasterOut(0x0);
   if(fDphiUnfolding) canvasMasterOut = new TCanvas("defaultOut", "defaultOut");
   (fDphiUnfolding) ? canvasMISC->Divide(4, 2) : canvasMISC->Divide(4, 1);
   TDirectoryFile* defDir(0x0);
   TCanvas* canvasProfiles(new TCanvas("canvasProfiles", "canvasProfiles"));
   canvasProfiles->Divide(2);
   TProfile* ratioProfile(new TProfile("ratioProfile", "ratioProfile", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   TProfile* v2Profile(new TProfile("v2Profile", "v2Profile", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   // get an estimate of the number of outputs and find the default set
   Int_t cacheMe(0);
   for(Int_t i(0); i < listOfKeys->GetSize(); i++) {
       if(dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()))) {
           if(!strcmp(listOfKeys->At(i)->GetName(), def.Data())) defDir = dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()));
           cacheMe++;
       }
   }
   Int_t rows(TMath::Floor(cacheMe/(float)columns)+((cacheMe%columns)>0));
   canvasIn->Divide(columns, rows);
   if(canvasOut) canvasOut->Divide(columns, rows);
   canvasRatioMeasuredRefoldedIn->Divide(columns, rows);
   if(canvasRatioMeasuredRefoldedOut) canvasRatioMeasuredRefoldedOut->Divide(columns, rows);
   canvasSpectraIn->Divide(columns, rows);
   if(canvasSpectraOut) canvasSpectraOut->Divide(columns, rows);
   if(canvasRatio) canvasRatio->Divide(columns, rows);
   if(canvasV2) canvasV2->Divide(columns, rows);

   canvasMasterIn->Divide(columns, rows);
   if(canvasMasterOut) canvasMasterOut->Divide(columns, rows);
   // extract the default output 
   TH1D* defaultUnfoldedJetSpectrumIn(0x0);
   TH1D* defaultUnfoldedJetSpectrumOut(0x0);
   if(defDir) {
       TDirectoryFile* defDirIn = (TDirectoryFile*)defDir->Get(Form("InPlane___%s", def.Data()));
       TDirectoryFile* defDirOut = (TDirectoryFile*)defDir->Get(Form("OutOfPlane___%s", def.Data()));
       if(defDirIn) defaultUnfoldedJetSpectrumIn = (TH1D*)defDirIn->Get(Form("UnfoldedSpectrum_in_%s", def.Data()));
       if(defDirOut) defaultUnfoldedJetSpectrumOut = (TH1D*)defDirOut->Get(Form("UnfoldedSpectrum_out_%s", def.Data()));
       printf(" > succesfully extracted default results < \n");
   }

   // loop through the directories, only plot the graphs if the deconvolution converged
   TDirectoryFile* tempDir(0x0); 
   TDirectoryFile* tempIn(0x0);
   TDirectoryFile*  tempOut(0x0);
   for(Int_t i(0), j(0); i < listOfKeys->GetSize(); i++) {
       // read the directory by its name
       tempDir = dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()));
       if(!tempDir) continue;
       TString dirName(tempDir->GetName());
       // try to read the in- and out of plane subdirs
       tempIn = (TDirectoryFile*)tempDir->Get(Form("InPlane___%s", dirName.Data()));
       tempOut = (TDirectoryFile*)tempDir->Get(Form("OutOfPlane___%s", dirName.Data()));
       j++;
       if(!tempIn) {    // attempt to get MakeAU output
           TString stringArray[] = {"a", "b", "c", "d", "e", "f", "g", "h"};
           TCanvas* tempPearson(new TCanvas(Form("pearson_%s", dirName.Data()), Form("pearson_%s", dirName.Data())));
           tempPearson->Divide(4,2);
           TCanvas* tempRatio(new TCanvas(Form("ratio_%s", dirName.Data()), Form("ratio_%s", dirName.Data())));
           tempRatio->Divide(4,2);
           TCanvas* tempResult(new TCanvas(Form("result_%s", dirName.Data()), Form("result_%s", dirName.Data())));
           tempResult->Divide(4,2);
           for(Int_t q(0); q < 8; q++) {
               tempIn = (TDirectoryFile*)tempDir->Get(Form("%s___%s", stringArray[q].Data(), dirName.Data()));
               if(tempIn) {
                       // to see if the unfolding converged try to extract the pearson coefficients
                       TH2D* pIn((TH2D*)tempIn->Get(Form("PearsonCoefficients_in_%s", dirName.Data())));
                       if(pIn) {
                       printf(" - %s in plane converged \n", dirName.Data());
                           tempPearson->cd(1+q);
                            Style(gPad, "PEARSON");
                            pIn->DrawCopy("colz");
                            TGraphErrors* rIn((TGraphErrors*)tempIn->Get(Form("RatioRefoldedMeasured_%s", dirName.Data())));
                            if(rIn) {
                                Style(rIn);
                                printf(" > found RatioRefoldedMeasured < \n");
                                tempRatio->cd(q+1);
                                rIn->SetFillColor(kRed);
                                rIn->Draw("ap");
                            }
                            TH1D* dvector((TH1D*)tempIn->Get("dVector"));
                            TH1D* avalue((TH1D*)tempIn->Get("SingularValuesOfAC"));
                            TH2D* rm((TH2D*)tempIn->Get(Form("ResponseMatrixIn_%s", dirName.Data())));
                            TH1D* eff((TH1D*)tempIn->Get(Form("KinematicEfficiencyIn_%s", dirName.Data())));
                            if(dvector && avalue && rm && eff) {
                                Style(dvector);
                                Style(avalue);
                                Style(rm);
                                Style(eff);
                                canvasMISC->cd(1);
                                Style(gPad, "SPECTRUM");
                                dvector->DrawCopy();
                                canvasMISC->cd(2);
                                Style(gPad, "SPECTRUM");
                                avalue->DrawCopy();
                                canvasMISC->cd(3);
                                Style(gPad, "PEARSON");
                                rm->DrawCopy("colz");
                                canvasMISC->cd(4);
                                Style(gPad, "GRID");
                                eff->DrawCopy();
                            } else if(rm && eff) {
                                Style(rm);
                                Style(eff);
                                canvasMISC->cd(3);
                                Style(gPad, "PEARSON");
                                rm->DrawCopy("colz");
                                canvasMISC->cd(4);
                                Style(gPad, "GRID");
                                eff->DrawCopy();
                            }
                        }
                       TH1D* inputSpectrum((TH1D*)tempIn->Get(Form("InputSpectrum_in_%s", dirName.Data())));
                       TH1D* unfoldedSpectrum((TH1D*)tempIn->Get(Form("UnfoldedSpectrum_in_%s", dirName.Data())));
                       TH1D* refoldedSpectrum((TH1D*)tempIn->Get(Form("RefoldedSpectrum_in_%s", dirName.Data())));
                       if(inputSpectrum && unfoldedSpectrum && refoldedSpectrum) {
                           if(defaultUnfoldedJetSpectrumIn) {
                               Style(defaultUnfoldedJetSpectrumIn, kBlue, kUnfoldedSpectrum);
                               TH1D* temp((TH1D*)defaultUnfoldedJetSpectrumIn->Clone(Form("defaultUnfoldedJetSpectrumIn_%s", dirName.Data())));
                               temp->Divide(unfoldedSpectrum);
                               temp->SetTitle(Form("ratio default unfolded / %s", dirName.Data()));
                               temp->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                               temp->GetYaxis()->SetTitle(Form("%s / %s", def.Data(), dirName.Data()));
                               canvasMasterIn->cd(j);
                               temp->GetYaxis()->SetRangeUser(0., 2);
                               temp->DrawCopy();
                           }
                           TH1F* fitStatus((TH1F*)tempIn->Get(Form("fitStatus_%s_in", dirName.Data())));
                           tempResult->cd(q+1);
                           Style(gPad);
                           Style(unfoldedSpectrum, kRed, kUnfoldedSpectrum);
                           unfoldedSpectrum->DrawCopy();
                           Style(inputSpectrum, kGreen, kMeasuredSpectrum);
                           inputSpectrum->DrawCopy("same");
                           Style(refoldedSpectrum, kBlue, kFoldedSpectrum);
                           refoldedSpectrum->DrawCopy("same");
                           TLegend* l(AddLegend(gPad));
                           Style(l);
                           if(fitStatus && fitStatus->GetNbinsX() == 4) { // only available in chi2 fit
                               Float_t chi(fitStatus->GetBinContent(1));
                               Float_t pen(fitStatus->GetBinContent(2));
                               Int_t dof((int)fitStatus->GetBinContent(3));
                               Float_t beta(fitStatus->GetBinContent(4));
                               l->AddEntry((TObject*)0, Form("#chi %.2f \tP %.2f \tDOF %i, #beta %.2f", chi, pen, dof, beta), "");
                           } else if (fitStatus) { // only available in SVD fit
                               Int_t reg((int)fitStatus->GetBinContent(1));
                               l->AddEntry((TObject*)0, Form("REG %i", reg), "");
                           }
                       }
               }
           }
       }
       if(tempIn) { 
           // to see if the unfolding converged try to extract the pearson coefficients
           TH2D* pIn((TH2D*)tempIn->Get(Form("PearsonCoefficients_in_%s", dirName.Data())));
           if(pIn) {
               printf(" - %s in plane converged \n", dirName.Data());
               canvasIn->cd(j);
               Style(gPad, "PEARSON");
               pIn->DrawCopy("colz");
               TGraphErrors* rIn((TGraphErrors*)tempIn->Get(Form("RatioRefoldedMeasured_%s", dirName.Data())));
               if(rIn) {
                   Style(rIn);
                   printf(" > found RatioRefoldedMeasured < \n");
                   canvasRatioMeasuredRefoldedIn->cd(j);
                   rIn->SetFillColor(kRed);
                   rIn->Draw("ap");
               }
               TH1D* dvector((TH1D*)tempIn->Get("dVector"));
               TH1D* avalue((TH1D*)tempIn->Get("SingularValuesOfAC"));
               TH2D* rm((TH2D*)tempIn->Get(Form("ResponseMatrixIn_%s", dirName.Data())));
               TH1D* eff((TH1D*)tempIn->Get(Form("KinematicEfficiencyIn_%s", dirName.Data())));
               if(dvector && avalue && rm && eff) {
                   Style(dvector);
                   Style(avalue);
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(1);
                   Style(gPad, "SPECTRUM");
                   dvector->DrawCopy();
                   canvasMISC->cd(2);
                   Style(gPad, "SPECTRUM");
                   avalue->DrawCopy();
                   canvasMISC->cd(3);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(4);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(3);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(4);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               }
           }
           TH1D* inputSpectrum((TH1D*)tempIn->Get(Form("InputSpectrum_in_%s", dirName.Data())));
           TH1D* unfoldedSpectrum((TH1D*)tempIn->Get(Form("UnfoldedSpectrum_in_%s", dirName.Data())));
           TH1D* refoldedSpectrum((TH1D*)tempIn->Get(Form("RefoldedSpectrum_in_%s", dirName.Data())));
           if(inputSpectrum && unfoldedSpectrum && refoldedSpectrum) {
               if(defaultUnfoldedJetSpectrumIn) {
                   Style(defaultUnfoldedJetSpectrumIn, kBlue, kUnfoldedSpectrum);
                   TH1D* temp((TH1D*)defaultUnfoldedJetSpectrumIn->Clone(Form("defaultUnfoldedJetSpectrumIn_%s", dirName.Data())));
                   temp->Divide(unfoldedSpectrum);
                   temp->SetTitle(Form("ratio default unfolded / %s", dirName.Data()));
                   temp->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                   temp->GetYaxis()->SetTitle(Form("%s / %s", def.Data(), dirName.Data()));
                   canvasMasterIn->cd(j);
                   temp->GetYaxis()->SetRangeUser(0., 2);
                   temp->DrawCopy();
               }
               TH1F* fitStatus((TH1F*)tempIn->Get(Form("fitStatus_%s_in", dirName.Data())));
               canvasSpectraIn->cd(j);
               Style(gPad);
               Style(unfoldedSpectrum, kRed, kUnfoldedSpectrum);
               unfoldedSpectrum->DrawCopy();
               Style(inputSpectrum, kGreen, kMeasuredSpectrum);
               inputSpectrum->DrawCopy("same");
               Style(refoldedSpectrum, kBlue, kFoldedSpectrum);
               refoldedSpectrum->DrawCopy("same");
               TLegend* l(AddLegend(gPad));
               Style(l);
               if(fitStatus && fitStatus->GetNbinsX() == 4) { // only available in chi2 fit
                   Float_t chi(fitStatus->GetBinContent(1));
                   Float_t pen(fitStatus->GetBinContent(2));
                   Int_t dof((int)fitStatus->GetBinContent(3));
                   Float_t beta(fitStatus->GetBinContent(4));
                   l->AddEntry((TObject*)0, Form("#chi %.2f \tP %.2f \tDOF %i, #beta %.2f", chi, pen, dof, beta), "");
               } else if (fitStatus) { // only available in SVD fit
                   Int_t reg((int)fitStatus->GetBinContent(1));
                   l->AddEntry((TObject*)0, Form("REG %i", reg), "");
               }
           }
       }
       if(tempOut) {
           TH2D* pOut((TH2D*)tempOut->Get(Form("PearsonCoefficients_out_%s", dirName.Data())));
           if(pOut) {
               printf(" - %s out of plane converged \n", dirName.Data());
               canvasOut->cd(j);
               Style(gPad, "PEARSON");
               pOut->DrawCopy("colz");
               TGraphErrors* rOut((TGraphErrors*)tempOut->Get(Form("RatioRefoldedMeasured_%s", dirName.Data())));
               if(rOut) {
                   Style(rOut);
                   printf(" > found RatioRefoldedMeasured < \n");
                   canvasRatioMeasuredRefoldedOut->cd(j);
                   rOut->SetFillColor(kRed);
                   rOut->Draw("ap");
               }
               TH1D* dvector((TH1D*)tempOut->Get("dVector"));
               TH1D* avalue((TH1D*)tempOut->Get("SingularValuesOfAC"));
               TH2D* rm((TH2D*)tempOut->Get(Form("ResponseMatrixOut_%s", dirName.Data())));
               TH1D* eff((TH1D*)tempOut->Get(Form("KinematicEfficiencyOut_%s", dirName.Data())));
               if(dvector && avalue && rm && eff) {
                   Style(dvector);
                   Style(avalue);
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(5);
                   Style(gPad, "SPECTRUM");
                   dvector->DrawCopy();
                   canvasMISC->cd(6);
                   Style(gPad, "SPECTRUM");
                   avalue->DrawCopy();
                   canvasMISC->cd(7);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(8);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(7);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(8);
                   Style(gPad, "GRID");
                   eff->DrawCopy();
               }
           }
           TH1D* inputSpectrum((TH1D*)tempOut->Get(Form("InputSpectrum_out_%s", dirName.Data())));
           TH1D* unfoldedSpectrum((TH1D*)tempOut->Get(Form("UnfoldedSpectrum_out_%s", dirName.Data())));
           TH1D* refoldedSpectrum((TH1D*)tempOut->Get(Form("RefoldedSpectrum_out_%s", dirName.Data())));
           if(inputSpectrum && unfoldedSpectrum && refoldedSpectrum) {
               if(defaultUnfoldedJetSpectrumOut) {
                   Style(defaultUnfoldedJetSpectrumOut, kBlue, kUnfoldedSpectrum);
                   TH1D* temp((TH1D*)defaultUnfoldedJetSpectrumOut->Clone(Form("defaultUnfoldedJetSpectrumOut_%s", dirName.Data())));
                   temp->Divide(unfoldedSpectrum);
                   temp->SetTitle(Form("ratio default unfolded / %s", dirName.Data()));
                   temp->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
                   temp->GetYaxis()->SetTitle(Form("%s / %s", def.Data(), dirName.Data()));
                   canvasMasterOut->cd(j);
                   temp->GetYaxis()->SetRangeUser(0., 2.);
                   temp->DrawCopy();
               }
               TH1F* fitStatus((TH1F*)tempOut->Get(Form("fitStatus_%s_out", dirName.Data())));
               canvasSpectraOut->cd(j);
               Style(gPad);
               Style(unfoldedSpectrum, kRed, kUnfoldedSpectrum);
               unfoldedSpectrum->DrawCopy();
               Style(inputSpectrum, kGreen, kMeasuredSpectrum);
               inputSpectrum->DrawCopy("same");
               Style(refoldedSpectrum, kBlue, kFoldedSpectrum);
               refoldedSpectrum->DrawCopy("same");
               TLegend* l(AddLegend(gPad));
               Style(l);
               if(fitStatus && fitStatus->GetNbinsX() == 4) { // only available in chi2 fit
                   Float_t chi(fitStatus->GetBinContent(1));
                   Float_t pen(fitStatus->GetBinContent(2));
                   Int_t dof((int)fitStatus->GetBinContent(3));
                   Float_t beta(fitStatus->GetBinContent(4));
                   l->AddEntry((TObject*)0, Form("#chi %.2f \tP %.2f \tDOF %i, #beta %.2f", chi, pen, dof, beta), "");
               } else if (fitStatus) { // only available in SVD fit
                   Int_t reg((int)fitStatus->GetBinContent(1));
                   l->AddEntry((TObject*)0, Form("REG %i", reg), "");
               }
           }
       }
       if(canvasRatio && canvasV2) {
           canvasRatio->cd(j);
           TGraphErrors* ratioYield((TGraphErrors*)tempDir->Get(Form("RatioInOutPlane_%s", dirName.Data())));
           Double_t _tempx(0), _tempy(0);
           if(ratioYield) {
               Style(ratioYield);
               for(Int_t b(0); b < fBinsTrue->GetSize(); b++) {
                   ratioYield->GetPoint(b, _tempx, _tempy);
                   ratioProfile->Fill(_tempx, _tempy);
               }
               ratioProfile->GetYaxis()->SetRangeUser(-0., 2.);
               ratioProfile->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioYield->GetYaxis()->SetRangeUser(-0., 2.);
               ratioYield->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioYield->SetFillColor(kRed);
               ratioYield->Draw("ap");
           }
           canvasV2->cd(j);
           TGraphErrors* ratioV2((TGraphErrors*)tempDir->Get(Form("v2_%s", dirName.Data())));
           if(ratioV2) {
               Style(ratioV2);
               for(Int_t b(0); b < fBinsTrue->GetSize(); b++) {
                   ratioV2->GetPoint(b, _tempx, _tempy);
                   v2Profile->Fill(_tempx, _tempy);

               }
               v2Profile->GetYaxis()->SetRangeUser(-0., 2.);
               v2Profile->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioV2->GetYaxis()->SetRangeUser(-.25, .75);
               ratioV2->GetXaxis()->SetRangeUser(rangeLow, rangeUp);
               ratioV2->SetFillColor(kRed);
               ratioV2->Draw("ap");
           }
       }
   }
   TFile output(out.Data(), "RECREATE");
   canvasProfiles->cd(1);
   Style(ratioProfile);
   ratioProfile->DrawCopy();
   canvasProfiles->cd(2);
   Style(v2Profile);
   v2Profile->DrawCopy();
   SavePadToPDF(canvasProfiles);
   canvasProfiles->Write(); 
   SavePadToPDF(canvasIn);
   canvasIn->Write();
   if(canvasOut) {
       SavePadToPDF(canvasOut);
       canvasOut->Write();
   }
   SavePadToPDF(canvasRatioMeasuredRefoldedIn);
   canvasRatioMeasuredRefoldedIn->Write();
   if(canvasRatioMeasuredRefoldedOut) {
       SavePadToPDF(canvasRatioMeasuredRefoldedOut);
       canvasRatioMeasuredRefoldedOut->Write();
   }
   SavePadToPDF(canvasSpectraIn);
   canvasSpectraIn->Write();
   if(canvasSpectraOut) {
       SavePadToPDF(canvasSpectraOut);
       canvasSpectraOut->Write();
       SavePadToPDF(canvasRatio);
       canvasRatio->Write();
       SavePadToPDF(canvasV2);
       canvasV2->Write();
   }
   SavePadToPDF(canvasMasterIn);
   canvasMasterIn->Write();
   if(canvasMasterOut) {
       SavePadToPDF(canvasMasterOut);
       canvasMasterOut->Write();
   }
   SavePadToPDF(canvasMISC);
   canvasMISC->Write();
   output.Write();
   output.Close();
}
//_____________________________________________________________________________
void AliJetFlowTools::BootstrapSpectra(TString def, TString in, TString out) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   // function to interpret results of bootstrapping routine
   // TString def should hold the true emperical distribution
   if(fOutputFile && !fOutputFile->IsZombie()) fOutputFile->Close();
   TFile readMe(in.Data(), "READ");     // open file read-only
   if(readMe.IsZombie()) {
       printf(" > Fatal error, couldn't read %s for post processing ! < \n", in.Data());
       return;
   }
   printf("\n\n\n\t\t BootstrapSpectra \n > Recovered the following file structure : \n <");
   readMe.ls();
   TList* listOfKeys((TList*)readMe.GetListOfKeys());
   if(!listOfKeys) {
       printf(" > Fatal error, couldn't retrieve list of keys. Input file might have been corrupted ! < \n");
       return;
   }
   // get an estimate of the number of outputs and find the default set
   TDirectoryFile* defDir(0x0);
   for(Int_t i(0); i < listOfKeys->GetSize(); i++) {
       if(dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()))) {
           if(!strcmp(listOfKeys->At(i)->GetName(), def.Data())) defDir = dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()));
       }
   }

   // extract the default output this is the 'emperical' distribution 
   // w.r.t. which the other distributions will be evaluated
   TH1D* defaultUnfoldedJetSpectrumIn(0x0);
   TH1D* defaultUnfoldedJetSpectrumOut(0x0);
   TH1D* defaultInputSpectrumIn(0x0);
   TH1D* defaultInputSpectrumOut(0x0);
   TGraphErrors* v2Emperical(0x0);
   if(defDir) {
       TDirectoryFile* defDirIn = (TDirectoryFile*)defDir->Get(Form("InPlane___%s", def.Data()));
       TDirectoryFile* defDirOut = (TDirectoryFile*)defDir->Get(Form("OutOfPlane___%s", def.Data()));
       if(defDirIn) {
           defaultUnfoldedJetSpectrumIn = (TH1D*)defDirIn->Get(Form("UnfoldedSpectrum_in_%s", def.Data()));
           defaultInputSpectrumIn = (TH1D*)defDirIn->Get(Form("InputSpectrum_in_%s", def.Data()));
       }
       if(defDirOut) {
           defaultUnfoldedJetSpectrumOut = (TH1D*)defDirOut->Get(Form("UnfoldedSpectrum_out_%s", def.Data()));
           defaultInputSpectrumOut = (TH1D*)defDirOut->Get(Form("InputSpectrum_out_%s", def.Data()));
       }
   }
   if((!defaultUnfoldedJetSpectrumIn && defaultUnfoldedJetSpectrumOut)) {
       printf(" BootstrapSpectra: couldn't extract default spectra, aborting! \n");
       return;
   }
   else v2Emperical = GetV2(defaultUnfoldedJetSpectrumIn, defaultUnfoldedJetSpectrumOut, fEventPlaneRes);

   // now that we know for sure that the input is in place, reserve the bookkeeping histograms
   TH1F* delta0[fBinsTrue->GetSize()-1];
   for(Int_t i(0); i < fBinsTrue->GetSize()-1; i++) {
       delta0[i] = new TH1F(Form("delta%i_%.2f-%.2f_gev", i, fBinsTrue->At(i), fBinsTrue->At(i+1)), Form("#Delta_{0, %i} p_{T} %.2f-%.2f GeV/c", i, fBinsTrue->At(i), fBinsTrue->At(i+1)), 30, -1., 1.);//.15, .15);
       delta0[i]->Sumw2();
   }
   // and a canvas for illustration purposes only
   TCanvas* spectraIn(new TCanvas("spectraIn", "spectraIn"));
   TCanvas* spectraOut(new TCanvas("spectraOut", "spectraOut"));
   // common reference (in this case the generated v2)
   TF1* commonReference = new TF1("v2_strong_bkg_comb", "(x<=3)*.07+(x>3&&x<5)*(.07-(x-3)*.035)+(x>30&&x<40)*(x-30)*.005+(x>40)*.05", 0, 200);

   // loop through the directories, only plot the graphs if the deconvolution converged
   TDirectoryFile* tempDir(0x0); 
   TDirectoryFile* tempIn(0x0);
   TDirectoryFile*  tempOut(0x0);
   TH1D* unfoldedSpectrumIn(0x0);
   TH1D* unfoldedSpectrumOut(0x0);
   TH1D* measuredSpectrumIn(0x0);
   TH1D* measuredSpectrumOut(0x0);
   for(Int_t i(0), j(0); i < listOfKeys->GetSize(); i++) {
       // read the directory by its name
       tempDir = dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()));
       if(!tempDir) continue;
       TString dirName(tempDir->GetName());
       // read only bootstrapped outputs
       if(!dirName.Contains("bootstrap")) continue;
       // try to read the in- and out of plane subdirs
       tempIn = (TDirectoryFile*)tempDir->Get(Form("InPlane___%s", dirName.Data()));
       tempOut = (TDirectoryFile*)tempDir->Get(Form("OutOfPlane___%s", dirName.Data()));
       j++;
       // extract the unfolded spectra only if both in- and out-of-plane converted (i.e. pearson coefficients were saved)
       if(tempIn) {
           if(!(TH2D*)tempIn->Get(Form("PearsonCoefficients_in_%s", dirName.Data()))) continue;
           unfoldedSpectrumIn = (TH1D*)tempIn->Get(Form("UnfoldedSpectrum_in_%s", dirName.Data()));
           measuredSpectrumIn = (TH1D*)tempIn->Get(Form("InputSpectrum_in_%s", dirName.Data()));
           spectraIn->cd();
           (j==1) ? measuredSpectrumIn->DrawCopy() : measuredSpectrumIn->DrawCopy("same");
       }
       if(tempOut) {
           if(!(TH2D*)tempOut->Get(Form("PearsonCoefficients_out_%s", dirName.Data()))) continue;
           unfoldedSpectrumOut = (TH1D*)tempOut->Get(Form("UnfoldedSpectrum_out_%s", dirName.Data()));
           measuredSpectrumOut = (TH1D*)tempOut->Get(Form("InputSpectrum_out_%s", dirName.Data()));
           spectraOut->cd();
           (j==1) ? measuredSpectrumOut->DrawCopy() : measuredSpectrumOut->DrawCopy("same");
       }
       // get v2 with statistical uncertainties from the extracted spectra
       TGraphErrors* v2Bootstrapped(GetV2(unfoldedSpectrumIn, unfoldedSpectrumOut, fEventPlaneRes));
       // and loop over all bins to fill the bookkeeping histograms
       Double_t yBoot(0), yEmp(0), xDud(0);
       // optional: common reference (in this case the sampled v2 value)
       for(Int_t k(0); k < fBinsTrue->GetSize()-1; k++) {
           // read values point by point (passed by reference)
           v2Bootstrapped->GetPoint(k, xDud, yBoot);
           v2Emperical->GetPoint(k, xDud, yEmp);
           if(commonReference) {
               if(!commonReference->Eval(xDud)<1e-10) {
                   yEmp/=commonReference->Eval(xDud);
                   yBoot/=commonReference->Eval(xDud);
               } else { // if reference equals zero, take emperical distribution as reference
                   yEmp/=yEmp;  // 1
                   yBoot/=yEmp;
               }
           }
           cout << " yEmp " << yEmp << " yBoot  " << yBoot << endl;
           // save relative difference per pt bin
           if(TMath::Abs(yBoot)>1e-10) delta0[k]->Fill(yEmp - yBoot);
       }
   }
   // extracting final results now, as first estimate just a gaus fit to the distributions 
   // (should be changed perhaps to proper rms eventually)
   // attach relevant data to current buffer in the same loop
   TFile output(out.Data(), "RECREATE");
   defaultInputSpectrumIn->SetLineColor(kRed);
   spectraIn->cd();
   defaultInputSpectrumIn->DrawCopy("same");
   defaultInputSpectrumOut->SetLineColor(kRed);
   spectraOut->cd();
   defaultInputSpectrumOut->DrawCopy("same");
   spectraIn->Write();
   spectraOut->Write();

   TH1F* delta0spread(new TH1F("delta0spread", "#sigma(#Delta_{0})", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   TH1F* unfoldingError(new TH1F("unfoldingError", "error from unfolding algorithm", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
   TF1* gaus(new TF1("gaus", "gaus"/*[0]*exp(-0.5*((x-[1])/[2])**2)"*/, -1., 1.));
   Double_t xDud(0), yDud(0);
   for(Int_t i(0); i < fBinsTrue->GetSize()-1; i++) {
       delta0[i]->Fit(gaus);
       delta0[i]->GetYaxis()->SetTitle("counts");
       delta0[i]->GetXaxis()->SetTitle("(v_{2, jet}^{measured} - v_{2, jet}^{generated}) / input v_{2}");
       delta0spread->SetBinContent(i+1, gaus->GetParameter(1)); // mean of gaus
       delta0spread->SetBinError(i+1, gaus->GetParameter(2));   // sigma of gaus
       cout << " mean " << gaus->GetParameter(1) << endl;
       cout << " sigm " << gaus->GetParameter(2) << endl;
       delta0[i]->Write();
       v2Emperical->GetPoint(i, xDud, yDud);
       unfoldingError->SetBinContent(i+1, 1e-10/*gaus->GetParameter(1)*/);
       if(commonReference && !commonReference->Eval(xDud)<1e-10) unfoldingError->SetBinError(i+1, v2Emperical->GetErrorY(i)/(commonReference->Eval(xDud)));
       else if(yDud>10e-10) unfoldingError->SetBinError(i+1, v2Emperical->GetErrorY(i)/yDud);
       else unfoldingError->SetBinError(i+1, 0.);
   }
   delta0spread->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
   delta0spread->GetYaxis()->SetTitle("(mean v_{2, jet}^{measured} - v_{2, jet}^{generated}) / input v_{2}");
   delta0spread->Write();
   unfoldingError->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
   unfoldingError->GetYaxis()->SetTitle("(mean v_{2, jet}^{measured} - v_{2, jet}^{generated}) / input v_{2}");
   unfoldingError->Write();
   // write the buffer and close the file
   output.Write();
   output.Close();
}
//_____________________________________________________________________________
Bool_t AliJetFlowTools::SetRawInput (
        TH2D* detectorResponse,  // detector response matrix
        TH1D* jetPtIn,           // in plane jet spectrum
        TH1D* jetPtOut,          // out of plane jet spectrum
        TH1D* dptIn,             // in plane delta pt distribution
        TH1D* dptOut,            // out of plane delta pt distribution
        Int_t eventCount) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // set input histograms manually
    fDetectorResponse   = detectorResponse;
    fSpectrumIn         = jetPtIn;
    fSpectrumOut        = jetPtOut;
    fDptInDist          = dptIn;
    fDptOutDist         = dptOut;
    // check if all data is provided
    if(!fDetectorResponse) {
        printf(" fDetectorResponse not found \n ");
        return kFALSE;
    }
    // check if the pt bin for true and rec have been set
    if(!fBinsTrue || !fBinsRec) {
        printf(" No true or rec bins set, please set binning ! \n");
        return kFALSE;
    }
    if(!fRMSSpectrumIn) { // initialie the profiles which will hold the RMS values. if binning changes in between unfolding
                        // procedures, these profiles will be nonsensical, user is responsible
        fRMSSpectrumIn = new TProfile("fRMSSpectrumIn", "fRMSSpectrumIn", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
        fRMSSpectrumOut = new TProfile("fRMSSpectrumOut", "fRMSSpectrumOut", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
        fRMSRatio = new TProfile("fRMSRatio", "fRMSRatio", fBinsTrue->GetSize()-1, fBinsTrue->GetArray());
    }
    // normalize spectra to event count if requested
    if(fNormalizeSpectra) {
        fEventCount = eventCount;
        if(fEventCount > 0) {
            fSpectrumIn->Sumw2();       // necessary for correct error propagation of scale
            fSpectrumOut->Sumw2();
            fSpectrumIn->Scale(1./((double)fEventCount));
            fSpectrumOut->Scale(1./((double)fEventCount));
        }
    }
    if(!fNormalizeSpectra && fEventCount > 0) {
        fSpectrumIn->Sumw2();       // necessary for correct error propagation of scale
        fSpectrumOut->Sumw2();
        fSpectrumIn->Scale(1./((double)fEventCount));
        fSpectrumOut->Scale(1./((double)fEventCount));
    }
    fDptIn = ConstructDPtResponseFromTH1D(fDptInDist, fAvoidRoundingError);
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)), Form("dpt_response_INPLANE_%i", fCentralityArray->At(0)));
    fDptIn->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T, jet}^{rec} [GeV/c]");
    fDptOut = ConstructDPtResponseFromTH1D(fDptOutDist, fAvoidRoundingError);
    fDptOut->SetNameTitle(Form("dpt_response_OUTOFPLANE_%i", fCentralityArray->At(0)), Form("dpt_response_OUTOFPLANE_%i", fCentralityArray->At(0)));
    fDptOut->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV/c]");
    fDptOut->GetYaxis()->SetTitle("p_{T, jet}^{rec} [GeV/c]");
    

    // set the flag to let other functions to read data from input file
    fRawInputProvided = kTRUE;
    return fRawInputProvided;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::GetRatio(TH1 *h1, TH1* h2, TString name, Bool_t appendFit, Int_t xmax) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // return ratio of h1 / h2
    // histograms can have different binning. errors are propagated as uncorrelated
    if(!(h1 && h2) ) {
        printf(" GetRatio called with NULL argument(s) \n ");
        return 0x0;
    }
    Int_t j(0);
    TGraphErrors *gr = new TGraphErrors();
    gr->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
    Double_t binCent(0.), ratio(0.), error2(0.), binWidth(0.);
    TH1* dud((TH1*)h1->Clone("dud"));
    dud->Sumw2();
    h1->Sumw2();
    h2->Sumw2();
    if(!dud->Divide(h1, h2)) {
        printf(" > ROOT failed to divide two histogams, dividing manually \n < ");
        for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
            binCent = h1->GetXaxis()->GetBinCenter(i);
            if(xmax > 0. && binCent > xmax) continue;
            j = h2->FindBin(binCent);
            binWidth = h1->GetXaxis()->GetBinWidth(i);
            if(h2->GetBinContent(j) > 0.) {
                ratio = h1->GetBinContent(i)/h2->GetBinContent(j);
                Double_t A = h1->GetBinError(i)/h1->GetBinContent(i);
                Double_t B = h2->GetBinError(i)/h2->GetBinContent(i);
                error2 = ratio*ratio*A*A+ratio*ratio*B*B;
                if(error2 > 0 ) error2 = TMath::Sqrt(error2);
                gr->SetPoint(i-1, binCent, ratio);
                gr->SetPointError(i-1, 0.5*binWidth, error2);
            }
        }
    } else {
        printf( " > using ROOT to divide two histograms < \n");
        for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
            binCent = dud->GetXaxis()->GetBinCenter(i);
            if(xmax > 0. && binCent > xmax) continue;
            binWidth = dud->GetXaxis()->GetBinWidth(i);
            gr->SetPoint(i-1,binCent,dud->GetBinContent(i));
            gr->SetPointError(i-1, 0.5*binWidth,dud->GetBinError(i));
        }
    }

    if(appendFit) {
        TF1* fit(new TF1("lin", "pol0", 10, 100));
        gr->Fit(fit);
    }
    if(strcmp(name, "")) gr->SetNameTitle(name.Data(), name.Data());
    if(dud) delete dud;
    return gr;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::GetV2(TH1 *h1, TH1* h2, Double_t r, TString name) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // get v2 from difference of in plane, out of plane yield
    // h1 must hold the in-plane yield, h2 holds the out of plane  yield
    // r is the event plane resolution for the chosen centrality
    if(!(h1 && h2) ) {
        printf(" GetV2 called with NULL argument(s) \n ");
        return 0x0;
    }
    TGraphErrors *gr = new TGraphErrors();
    gr->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
    Float_t binCent(0.), ratio(0.), error2(0.), binWidth(0.);
    Double_t pre(TMath::Pi()/(4.*r)), in(0.), out(0.), ein(0.), eout(0.);

    for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
        binCent = h1->GetXaxis()->GetBinCenter(i);
        binWidth = h1->GetXaxis()->GetBinWidth(i);
        if(h2->GetBinContent(i) > 0.) {
            in = h1->GetBinContent(i);
            ein = h1->GetBinError(i);
            out = h2->GetBinContent(i);
            eout = h2->GetBinError(i);
            ratio = pre*((in-out)/(in+out));
            error2 = (4.*out*out/(TMath::Power(in+out, 4)))*ein*ein+(4.*in*in/(TMath::Power(in+out, 4)))*eout*eout;
            error2 = error2*pre*pre;
            if(error2 > 0) error2 = TMath::Sqrt(error2);
            gr->SetPoint(i-1,binCent,ratio);
            gr->SetPointError(i-1,0.5*binWidth,error2);
        }
    }
    if(strcmp(name, "")) gr->SetNameTitle(name.Data(), name.Data());
    return gr;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::GetV2Histo(TH1 *h1, TH1* h2, Double_t r, TString name) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // get v2 from difference of in plane, out of plane yield
    // h1 must hold the in-plane yield, h2 holds the out of plane  yield
    // r is the event plane resolution for the chosen centrality
    if(!(h1 && h2) ) {
        printf(" GetV2 called with NULL argument(s) \n ");
        return 0x0;
    }
    TH1D* gr((TH1D*)h1->Clone(name.Data()));
    gr->GetXaxis()->SetTitle("p_{T, jet} [GeV/c]");
    Float_t ratio(0.), error2(0.);
    Double_t pre(TMath::Pi()/(4.*r)), in(0.), out(0.), ein(0.), eout(0.);

    for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
        if(h2->GetBinContent(i) > 0.) {
            in = h1->GetBinContent(i);
            ein = h1->GetBinError(i);
            out = h2->GetBinContent(i);
            eout = h2->GetBinError(i);
            ratio = pre*((in-out)/(in+out));
            error2 = (4.*out*out/(TMath::Power(in+out, 4)))*ein*ein+(4.*in*in/(TMath::Power(in+out, 4)))*eout*eout;
            error2 = error2*pre*pre;
            if(error2 > 0) error2 = TMath::Sqrt(error2);
            gr->SetBinContent(i,ratio);
            gr->SetBinError(i,error2);
        }
    }
    if(strcmp(name, "")) gr->SetNameTitle(name.Data(), name.Data());
    return gr;
}
//_____________________________________________________________________________
TH1F* AliJetFlowTools::ConvertGraphToHistogram(TGraphErrors* g)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // convert a tgrapherrors to a histogram
    if(!g) {
        printf(" > ConvertGraphToHistogram recevied a NULL pointer > \n");
        return 0x0;
    }
    // first get the frame which we'll use to build the histogram
    TH1F* hist(g->GetHistogram());
    Double_t xref(0), yref(0), yerr(0);
    // then copy the points and errors (remember graph starts at bin 0 and hist at 1)
    for(Int_t i(0); i < hist->GetNbinsX(); i++) {
        g->GetPoint(i, xref, yref);
        yerr = g->GetErrorY(i);
        hist->SetBinContent(i+1, yref);
        hist->SetBinError(i+1, yerr);
    }
    return hist;
}
//_____________________________________________________________________________
TGraphAsymmErrors* AliJetFlowTools::AddHistoErrorsToGraphErrors(TGraphAsymmErrors* g, TH1D* h)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // add the errors of a tgra histogram, assuming oncorrelated uncertainties
    // the points on the graph are returned ! 
    if(!(g&&h)) {
        printf(" > ConvertGraphToHistogram recevied a NULL pointer > \n");
        return 0x0;
    }

    printf("\n\n\n adding histo error to graph error \n\n\n\n");


    Double_t yerrL(0), yerrH(0), herr(0);
    // quadratic sum of the errors, assyming here symmetric ones
    for(Int_t i(0); i < h->GetNbinsX(); i++) {
        yerrH = g->GetErrorY(i);//high(i);
        yerrL = g->GetErrorY(i);//low(i);
        herr = h->GetBinError(i+1);
        cout << " graph error H " << yerrH << "\t histo error " << herr << endl;
        cout << " graph error L " << yerrL << "\t histo error " << herr << endl;
        yerrH = TMath::Sqrt(yerrH*yerrH+herr*herr);
        yerrL = TMath::Sqrt(yerrL*yerrL+herr*herr);
        g->SetPointError(i, g->GetErrorX(i), g->GetErrorX(i), yerrL, yerrH);
    }
    return g;
}
//_____________________________________________________________________________
Double_t AliJetFlowTools::GetRMSOfTH1(TH1* h, Double_t a, Double_t b)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // get rms between a and b for x-axis of histo
    if(!h || b < a ) return 0.;

    // a is a lower edge and b an upper edge
    Int_t binA(h->GetXaxis()->FindBin(a+.1));
    Int_t binB(h->GetXaxis()->FindBin(b+.1));   // on purpose
    Double_t RMS(0), mean(0), c(0);

    for(Int_t i(binA); i < binB; i++) {
        c++;
        mean += h->GetBinContent(i);
    }
    if( c > 0. ) {
        mean /= c;
        for(Int_t i(binA); i < binB; i++) RMS += (mean-h->GetBinContent(i))*(mean-h->GetBinContent(i));
        return TMath::Sqrt(RMS/c);
    }
    return 0.;
}
//_____________________________________________________________________________
TF1* AliJetFlowTools::GetErrorFromFit(TH1* h1, TH1* h2, Double_t a, Double_t b, 
        Float_t pivot, Bool_t subdueError, TString str, Bool_t setContent) 
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // return an error from a fit
    TF1* lin = new TF1("lin", Form("(x<%i)*(pol1)+(x>%i)*[2]", (int)pivot, (int)pivot), a, b);
    // clone the input
    TF1* fit_pol0 = new TF1("fit_pol0", "pol0", pivot, b);
    TF1* fit_pol1 = new TF1("fit_pol1", "pol1", a, pivot);

    TH1* h((TH1*)h1->Clone(str.Data()));
    // reset the guy and fill it with the maxima
    h->Reset();
    Double_t _a(0), _b(0);
    for(Int_t i(1); i < h->GetNbinsX() + 1; i++) {
        _a = TMath::Abs(h1->GetBinContent(i));
        _b = TMath::Abs(h2->GetBinContent(i));
        if(_a > _b) {
            h->SetBinContent(i, _a);
            h->SetBinError(i, h1->GetBinError(i));
        } else {
            h->SetBinContent(i, _b);
            h->SetBinError(i, h2->GetBinError(i));
        }
    }



    // fit to full error. root doesn't like fitting a step-function, so fit these
    // two components separately ...
    h->Fit(fit_pol0, "", "", pivot, b);
    lin->SetParameter(2, (subdueError) ? 0. : fit_pol0->GetParameter(0));

    h->Fit(fit_pol1, "", "", a, pivot);

    //check if it makes sense
    if(fit_pol1->GetParameter(1) > 0) {
        fit_pol1 = new TF1("fit_pol2", "pol0", a, pivot);
        h->Fit(fit_pol1, "", "", a, pivot);
        lin->SetParameter(0, fit_pol1->GetParameter(0));
        lin->SetParameter(1, 0);
    } else {
        lin->SetParameter(0, fit_pol1->GetParameter(0));
        lin->SetParameter(1, fit_pol1->GetParameter(1));
    }


    // .. and then add them together
    h->GetListOfFunctions()->Add(lin);

    if(setContent) {
        // update the histos with the fit result
        for(Int_t i(1); i < h->GetNbinsX() + 1; i++) {
            // calculate the integral in a given bin
            Double_t dud(lin->Integral(h->GetXaxis()->GetBinLowEdge(i),h->GetXaxis()->GetBinUpEdge(i))/h->GetXaxis()->GetBinWidth(i));
            // if it's larger than 0, assign an 'up' error
            h1->SetBinContent(i, 0);
            h2->SetBinContent(i, 0);
            h1->SetBinError(i, 0);
            h2->SetBinError(i, 0);
            // dont' show errors outside of range of interest
            if(a > h->GetXaxis()->GetBinLowEdge(i) || b < h->GetXaxis()->GetBinUpEdge(i)) continue;
            // assign them correctly to up or low
            if(dud > 0) {
                h1->SetBinContent(i, dud);
            } else {
                h2->SetBinContent(i, dud);
            }
        }
        h->Write();
    } 
    return lin;
}
//_____________________________________________________________________________
TGraphAsymmErrors* AliJetFlowTools::GetV2WithSystematicErrors(
        TH1* h1, TH1* h2, Double_t r, TString name, 
        TH1* relativeErrorInUp,
        TH1* relativeErrorInLow,
        TH1* relativeErrorOutUp,
        TH1* relativeErrorOutLow,
        Float_t rho) const
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // get v2 with asymmetric systematic error
    // note that this is ONLY the systematic error, no statistical error!
    // rho is the pearson correlation coefficient
    TGraphErrors* tempV2(GetV2(h1, h2, r, name));
    Double_t* ax = new Double_t[fBinsTrue->GetSize()-1];
    Double_t* ay = new Double_t[fBinsTrue->GetSize()-1];
    Double_t* axh = new Double_t[fBinsTrue->GetSize()-1];
    Double_t* axl = new Double_t[fBinsTrue->GetSize()-1];
    Double_t* ayh = new Double_t[fBinsTrue->GetSize()-1];
    Double_t* ayl = new Double_t[fBinsTrue->GetSize()-1];
    Double_t in(0.), out(0.), einUp(0.), einLow(0.), eoutUp(0.), eoutLow(0.), error2Up(0.), error2Low(0.);
    // loop through the bins and do error propagation
    for(Int_t i(0); i < fBinsTrue->GetSize()-1; i++) {
        // extract the absolute errors
        in = h1->GetBinContent(i+1);
        einUp = TMath::Abs(in*relativeErrorInUp->GetBinContent(i+1));
        einLow = TMath::Abs(in*relativeErrorInLow->GetBinContent(1+i));
        out = h2->GetBinContent(i+1);
        eoutUp = TMath::Abs(out*relativeErrorOutUp->GetBinContent(1+i));
        eoutLow = TMath::Abs(out*relativeErrorOutLow->GetBinContent(1+i));
        // get the error squared
        if(rho <= 0) {
            error2Up = TMath::Power(((r*4.)/(TMath::Pi())),-2.)*((4.*out*out/(TMath::Power(in+out, 4)))*einUp*einUp+(4.*in*in/(TMath::Power(in+out, 4)))*eoutLow*eoutLow);
            error2Low =TMath::Power(((r*4.)/(TMath::Pi())),-2.)*((4.*out*out/(TMath::Power(in+out, 4)))*einLow*einLow+(4.*in*in/(TMath::Power(in+out, 4)))*eoutUp*eoutUp);
        } else {
            error2Up = TMath::Power(((r*4.)/(TMath::Pi())),-2.)*((4.*out*out/(TMath::Power(in+out, 4)))*einUp*einUp+(4.*in*in/(TMath::Power(in+out, 4)))*eoutUp*eoutUp-((8.*out*in)/(TMath::Power(in+out, 4)))*rho*einUp*eoutUp);
            error2Low =TMath::Power(((r*4.)/(TMath::Pi())),-2.)*((4.*out*out/(TMath::Power(in+out, 4)))*einLow*einLow+(4.*in*in/(TMath::Power(in+out, 4)))*eoutLow*eoutLow-((8.*out*in)/(TMath::Power(in+out, 4)))*rho*einLow*eoutLow);
        }
        if(error2Up > 0) error2Up = TMath::Sqrt(error2Up);
        if(error2Low > 0) error2Low = TMath::Sqrt(error2Low);
        // set the errors 
        ayh[i] = error2Up;
        ayl[i] = error2Low;
        // get the bin width (which is the 'error' on x)
        Double_t binWidth(h1->GetBinWidth(i+1));
        axl[i] = binWidth/2.;
        axh[i] = binWidth/2.;
        // now get the coordinate for the poin
        tempV2->GetPoint(i, ax[i], ay[i]);
      /*  if(rho<=0) {
            cout << "[print]  pt " << ax[i] << " errorLow" << gReductionFactor*ayl[i] << " errorUp " << gReductionFactor*ayh[i] << endl;
            cout << "[print]   %" << ax[i] << " low " << (gReductionFactor*ayl[i])/ay[i] << " up " << (gReductionFactor*ayh[i])/ay[i] << endl;
        } else {
            cout << "[print]  pt " << ax[i] << " errorLow" << ayl[i]*gReductionFactorCorr << " errorUp " << ayh[i]*gReductionFactorCorr << endl;
            cout << "[print]   %" << ax[i] << " low " << (gReductionFactorCorr*ayl[i])/ay[i] << " up " << (gReductionFactorCorr*ayh[i])/ay[i] << endl;
        }*/
    }
    // save the nominal ratio
    TGraphAsymmErrors* nominalError(new TGraphAsymmErrors(fBinsTrue->GetSize()-1, ax, ay, axl, axh, ayl, ayh));
    nominalError->SetName("v_{2} with shape uncertainty");
    // do some memory management
    delete tempV2;
    delete[] ax;
    delete[] ay;
    delete[] axh;
    delete[] axl;
    delete[] ayh;
    delete[] ayl;

    return nominalError;
}
//_____________________________________________________________________________
void AliJetFlowTools::WriteObject(TObject* object, TString suffix, Bool_t kill) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // write object, if a unique identifier is given the object is cloned
    // and the clone is saved. setting kill to true will delete the original obect from the heap
    if(!object) {
        printf(" > WriteObject:: called with NULL arguments \n ");
        return;
    } else if(!strcmp("", suffix.Data())) object->Write();
    else {
        TObject* newObject(object->Clone(Form("%s_%s", object->GetName(), suffix.Data())));
        newObject->Write();
    }
    if(kill) delete object;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::ConstructDPtResponseFromTH1D(TH1D* dpt, Bool_t AvoidRoundingError) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // construt a delta pt response matrix from supplied dpt distribution
    // binning is fine, set fBinsTrue and fBinsRec and call 'RebinTH2D' to 
    // do a weighted rebinning to a (coarser) dpt distribution
    // be careful with the binning of the dpt response: it should be equal to that
    // of the response matrix, otherwise dpt and response matrices cannot be multiplied
    //
    // the response matrix will be square and have the same binning
    // (min, max and granularity) of the input histogram
    Int_t bins(dpt->GetXaxis()->GetNbins());        // number of bins, will also be no of rows, columns
    Double_t _bins[bins+1];             // prepare array with bin borders
    for(Int_t i(0); i < bins; i++) _bins[i] = dpt->GetBinLowEdge(i+1);
    _bins[bins] = dpt->GetBinLowEdge(bins)+dpt->GetBinWidth(bins+1);    // get upper edge
    TH2D* res(new TH2D(Form("Response_from_%s", dpt->GetName()), Form("Response_from_%s", dpt->GetName()), bins, _bins, bins, _bins));
    for(Int_t j(0); j < bins+1; j++) {   // loop on pt true slices j
        Bool_t skip = kFALSE;
        for(Int_t k(0); k < bins+1; k++) {       // loop on pt gen slices k
            (skip) ? res->SetBinContent(j, k, 0.) : res->SetBinContent(j, k, dpt->GetBinContent(dpt->GetXaxis()->FindBin(k-j)));
            if(AvoidRoundingError && k > j && TMath::AreEqualAbs(dpt->GetBinContent(dpt->GetBinContent(k-j)), 0, 1e-8)) skip = kTRUE;
        }
    }
    return res;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::GetUnityResponse(TArrayD* binsTrue, TArrayD* binsRec, TString suffix) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   if(!binsTrue || !binsRec) {
        printf(" > GetUnityResponse:: function called with NULL arguments < \n");
        return 0x0;
    }
    TString name(Form("unityResponse_%s", suffix.Data()));
    TH2D* unity(new TH2D(name.Data(), name.Data(), binsTrue->GetSize()-1, binsTrue->GetArray(), binsRec->GetSize()-1, binsRec->GetArray()));
    for(Int_t i(0); i < binsTrue->GetSize(); i++) {
        for(Int_t j(0); j < binsRec->GetSize(); j++) {
            if(i==j) unity->SetBinContent(1+i, 1+j, 1.);
        }
    }
    return unity;
}
//_____________________________________________________________________________
void AliJetFlowTools::SaveConfiguration(Bool_t convergedIn, Bool_t convergedOut) const {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // save configuration parameters to histogram
    TH1F* summary = new TH1F("UnfoldingConfiguration","UnfoldingConfiguration", 20, -.5, 19.5);
    summary->SetBinContent(1, fBetaIn);
    summary->GetXaxis()->SetBinLabel(1, "fBetaIn");
    summary->SetBinContent(2, fBetaOut);
    summary->GetXaxis()->SetBinLabel(2, "fBetaOut");
    summary->SetBinContent(3, fCentralityArray->At(0));
    summary->GetXaxis()->SetBinLabel(3, "fCentralityArray[0]");
    summary->SetBinContent(4, (int)convergedIn);
    summary->GetXaxis()->SetBinLabel(4, "convergedIn");
    summary->SetBinContent(5, (int)convergedOut);
    summary->GetXaxis()->SetBinLabel(5, "convergedOut");
    summary->SetBinContent(6, (int)fAvoidRoundingError);
    summary->GetXaxis()->SetBinLabel(6, "fAvoidRoundingError");
    summary->SetBinContent(7, (int)fUnfoldingAlgorithm);
    summary->GetXaxis()->SetBinLabel(7, "fUnfoldingAlgorithm");
    summary->SetBinContent(8, (int)fPrior);
    summary->GetXaxis()->SetBinLabel(8, "fPrior");
    summary->SetBinContent(9, fSVDRegIn);
    summary->GetXaxis()->SetBinLabel(9, "fSVDRegIn");
    summary->SetBinContent(10, fSVDRegOut);
    summary->GetXaxis()->SetBinLabel(10, "fSVDRegOut");
    summary->SetBinContent(11, (int)fSVDToy);
    summary->GetXaxis()->SetBinLabel(11, "fSVDToy");
    summary->SetBinContent(12, fJetRadius);
    summary->GetXaxis()->SetBinLabel(12, "fJetRadius");
    summary->SetBinContent(13, (int)fNormalizeSpectra);
    summary->GetXaxis()->SetBinLabel(13, "fNormalizeSpectra");
    summary->SetBinContent(14, (int)fSmoothenPrior);
    summary->GetXaxis()->SetBinLabel(14, "fSmoothenPrior");
    summary->SetBinContent(15, (int)fTestMode);
    summary->GetXaxis()->SetBinLabel(15, "fTestMode");
    summary->SetBinContent(16, (int)fUseDetectorResponse);
    summary->GetXaxis()->SetBinLabel(16, "fUseDetectorResponse");
    summary->SetBinContent(17, fBayesianIterIn);
    summary->GetXaxis()->SetBinLabel(17, "fBayesianIterIn");
    summary->SetBinContent(18, fBayesianIterOut);
    summary->GetXaxis()->SetBinLabel(18, "fBayesianIterOut");
    summary->SetBinContent(19, fBayesianSmoothIn);
    summary->GetXaxis()->SetBinLabel(19, "fBayesianSmoothIn");
    summary->SetBinContent(20, fBayesianSmoothOut);
    summary->GetXaxis()->SetBinLabel(20, "fBayesianSmoothOut");
}
//_____________________________________________________________________________
void AliJetFlowTools::ResetAliUnfolding() {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
     // ugly function: reset all unfolding parameters 
     TVirtualFitter* fitter(TVirtualFitter::GetFitter());
     if(fitter) {
         printf(" > Found fitter, will delete it < \n");
         delete fitter;
     }
     if(gMinuit) {
         printf(" > Found gMinuit, will re-create it < \n");
         delete gMinuit;
         gMinuit = new TMinuit;
     }
     AliUnfolding::fgCorrelationMatrix = 0;
     AliUnfolding::fgCorrelationMatrixSquared = 0;
     AliUnfolding::fgCorrelationCovarianceMatrix = 0;
     AliUnfolding::fgCurrentESDVector = 0;
     AliUnfolding::fgEntropyAPriori = 0;
     AliUnfolding::fgEfficiency = 0;
     AliUnfolding::fgUnfoldedAxis = 0;
     AliUnfolding::fgMeasuredAxis = 0;
     AliUnfolding::fgFitFunction = 0;
     AliUnfolding::fgMaxInput  = -1;
     AliUnfolding::fgMaxParams = -1;
     AliUnfolding::fgOverflowBinLimit = -1;
     AliUnfolding::fgRegularizationWeight = 10000;
     AliUnfolding::fgSkipBinsBegin = 0;
     AliUnfolding::fgMinuitStepSize = 0.1;
     AliUnfolding::fgMinuitPrecision = 1e-6;
     AliUnfolding::fgMinuitMaxIterations = 1000000;
     AliUnfolding::fgMinuitStrategy = 1.;
     AliUnfolding::fgMinimumInitialValue = kFALSE;
     AliUnfolding::fgMinimumInitialValueFix = -1;
     AliUnfolding::fgNormalizeInput = kFALSE;
     AliUnfolding::fgNotFoundEvents = 0;
     AliUnfolding::fgSkipBin0InChi2 = kFALSE;
     AliUnfolding::fgBayesianSmoothing  = 1;
     AliUnfolding::fgBayesianIterations = 10;
     AliUnfolding::fgDebug = kFALSE;
     AliUnfolding::fgCallCount = 0;
     AliUnfolding::fgPowern = 5;
     AliUnfolding::fChi2FromFit = 0.;
     AliUnfolding::fPenaltyVal  = 0.;
     AliUnfolding::fAvgResidual = 0.;
     AliUnfolding::fgPrintChi2Details = 0;
     AliUnfolding::fgCanvas = 0;
     AliUnfolding::fghUnfolded = 0;     
     AliUnfolding::fghCorrelation = 0;  
     AliUnfolding::fghEfficiency = 0;
     AliUnfolding::fghMeasured = 0;   
     AliUnfolding::SetMinuitStepSize(1.);
     AliUnfolding::SetMinuitPrecision(1e-6);
     AliUnfolding::SetMinuitMaxIterations(100000);
     AliUnfolding::SetMinuitStrategy(2.);
     AliUnfolding::SetDebug(1);
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::ProtectHeap(TH1D* protect, Bool_t kill, TString suffix) const {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // protect heap by adding unique qualifier to name
    if(!protect) return 0x0;
    TH1D* p = (TH1D*)protect->Clone();
    TString tempString(fActiveString);
    tempString+=suffix;
    p->SetName(Form("%s_%s", protect->GetName(), tempString.Data()));
    p->SetTitle(Form("%s_%s", protect->GetTitle(), tempString.Data()));
    if(kill) delete protect;
    return p;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::ProtectHeap(TH2D* protect, Bool_t kill, TString suffix) const {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // protect heap by adding unique qualifier to name
    if(!protect) return 0x0;
    TH2D* p = (TH2D*)protect->Clone();
    TString tempString(fActiveString);
    tempString+=suffix;
    p->SetName(Form("%s_%s", protect->GetName(), tempString.Data()));
    p->SetTitle(Form("%s_%s", protect->GetTitle(), tempString.Data()));
    if(kill) delete protect;
    return p;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::ProtectHeap(TGraphErrors* protect, Bool_t kill, TString suffix) const {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // protect heap by adding unique qualifier to name
    if(!protect) return 0x0;
    TGraphErrors* p = (TGraphErrors*)protect->Clone();
    TString tempString(fActiveString);
    tempString+=suffix;
    p->SetName(Form("%s_%s", protect->GetName(), tempString.Data()));
    p->SetTitle(Form("%s_%s", protect->GetTitle(), tempString.Data()));
    if(kill) delete protect;
    return p;
}
//_____________________________________________________________________________
void AliJetFlowTools::MakeAU() {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // === azimuthal unfolding ===
    // 
    // unfolds the spectrum in delta phi bins, extracts the yield per bin, and does a fit
    // in transverse momentum and azimuthal correlation space to extract v2 params
    // settings are equal to the ones used for 'Make()'
    //
    // basic steps that are followed:
    // 1) rebin the raw output of the jet task to the desired binnings
    // 2) calls the unfolding routine
    // 3) writes output to file
    // can be repeated multiple times with different configurations

    Int_t low[] = {1, 6, 11, 16, 21, 26, 31, 36};
    Int_t up[] = {5, 10, 15, 20, 25, 30, 35, 40};
    TString stringArray[] = {"a", "b", "c", "d", "e", "f", "g", "h"};
    const Int_t ptBins(fBinsTrue->GetSize()-1);
    const Int_t dPhiBins(8);
    TH1D* dPtdPhi[fBinsTrue->GetSize()];
    for(Int_t i(0); i < ptBins; i++) dPtdPhi[i] = new TH1D(Form("dPtdPhi_%i", i), Form("dPtdPhi_%i", i), dPhiBins, 0, TMath::Pi());

    for(Int_t i(0); i < dPhiBins; i++) {
        // 1) manipulation of input histograms
        // check if the input variables are present
        if(!PrepareForUnfolding(low[i], up[i])) return;
        // 1a) resize the jet spectrum according to the binning scheme in fBinsTrue
        //     parts of the spectrum can end up in over or underflow bins
        TH1D* measuredJetSpectrumIn  = RebinTH1D(fSpectrumIn, fBinsRec, Form("resized_%s", stringArray[i].Data()), kFALSE);

        // 1b) resize the jet spectrum to 'true' bins. can serve as a prior and as a template for unfolding
        // the template will be used as a prior for the chi2 unfolding
        TH1D* measuredJetSpectrumTrueBinsIn  = RebinTH1D(fSpectrumIn, fBinsTrue, stringArray[i], kFALSE);

        // get the full response matrix from the dpt and the detector response
        fDetectorResponse = NormalizeTH2D(fDetectorResponse);
        // get the full response matrix. if test mode is chosen, the full response is replace by a unity matrix
        // so that unfolding should return the initial spectrum
        if(fUseDptResponse && fUseDetectorResponse) fFullResponseIn = MatrixMultiplication(fDptIn, fDetectorResponse);
        else if (fUseDptResponse && !fUseDetectorResponse) fFullResponseIn = fDptIn;
        else if (!fUseDptResponse && fUseDetectorResponse) fFullResponseIn = fDetectorResponse;
        else if (!fUseDptResponse && !fUseDetectorResponse && !fUnfoldingAlgorithm == AliJetFlowTools::kNone) return;
        // normalize each slide of the response to one
        NormalizeTH2D(fFullResponseIn);
        // resize to desired binning scheme
        TH2D* resizedResponseIn  = RebinTH2D(fFullResponseIn, fBinsTrue, fBinsRec, stringArray[i]);
        // get the kinematic efficiency
        TH1D* kinematicEfficiencyIn  = resizedResponseIn->ProjectionX();
        kinematicEfficiencyIn->SetNameTitle(Form("kin_eff_%s", stringArray[i].Data()), Form("kin_eff_%s", stringArray[i].Data()));
        // suppress the errors 
        for(Int_t j(0); j < kinematicEfficiencyIn->GetXaxis()->GetNbins(); j++) kinematicEfficiencyIn->SetBinError(1+j, 0.);
        TH1D* jetFindingEfficiency(0x0);
        if(fJetFindingEff) {
            jetFindingEfficiency = ProtectHeap(fJetFindingEff);
            jetFindingEfficiency->SetNameTitle(Form("%s_coarse", jetFindingEfficiency->GetName()), Form("%s_coarse", jetFindingEfficiency->GetName()));
            jetFindingEfficiency = RebinTH1D(jetFindingEfficiency, fBinsTrue);
        }
        // 2, 3) call the actual unfolding. results and transient objects are stored in a dedicated TDirectoryFile
        TH1D* unfoldedJetSpectrumIn(0x0);
        fActiveDir->cd();                   // select active dir
        TDirectoryFile* dirIn = new TDirectoryFile(Form("%s___%s", stringArray[i].Data(), fActiveString.Data()), Form("%s___%s", stringArray[i].Data(), fActiveString.Data()));
        dirIn->cd();                        // select inplane subdir
        // select the unfolding method
        unfoldedJetSpectrumIn = UnfoldWrapper(
            measuredJetSpectrumIn,
            resizedResponseIn,
            kinematicEfficiencyIn,
            measuredJetSpectrumTrueBinsIn,
            TString("dPtdPhiUnfolding"),
            jetFindingEfficiency);
        // arbitrarily save one of the full outputs (same for all dphi bins, avoid duplicates)
        if(i+1 == ptBins) {
            resizedResponseIn->SetNameTitle(Form("ResponseMatrix_%s", stringArray[i].Data()), Form("response matrix %s", stringArray[i].Data()));
            resizedResponseIn->SetXTitle("p_{T, jet}^{true} [GeV/c]");
            resizedResponseIn->SetYTitle("p_{T, jet}^{rec} [GeV/c]");
            resizedResponseIn = ProtectHeap(resizedResponseIn);
            resizedResponseIn->Write();
            kinematicEfficiencyIn->SetNameTitle(Form("KinematicEfficiency_%s", stringArray[i].Data()), Form("Kinematic efficiency, %s", stringArray[i].Data()));
            kinematicEfficiencyIn = ProtectHeap(kinematicEfficiencyIn);
            kinematicEfficiencyIn->Write();
            fDetectorResponse->SetNameTitle("DetectorResponse", "Detector response matrix");
            fDetectorResponse = ProtectHeap(fDetectorResponse, kFALSE);
            fDetectorResponse->Write();
            // optional histograms
            if(fSaveFull) {
                fSpectrumIn->SetNameTitle("[ORIG]JetSpectrum", Form("[INPUT] Jet spectrum, %s", stringArray[i].Data()));
                fSpectrumIn->Write();
                fDptInDist->SetNameTitle("[ORIG]DeltaPt", Form("#delta p_{T} distribution, %s", stringArray[i].Data()));
                fDptInDist->Write();
                fDptIn->SetNameTitle("[ORIG]DeltaPtMatrix", Form("#delta p_{T} matrix, %s", stringArray[i].Data()));
                fDptIn->Write();
                fFullResponseIn->SetNameTitle("ResponseMatrix", Form("Response matrix, %s", stringArray[i].Data()));
                fFullResponseIn->Write();
            }
        }
        fActiveDir->cd();
        fDeltaPtDeltaPhi->Write();
        fJetPtDeltaPhi->Write();
        
        TH1D* dud(ProtectHeap(unfoldedJetSpectrumIn, kTRUE, stringArray[i]));;
        Double_t integralError(0);
        // at this point in the code, the spectrum has been unfolded in a certain region of dPhi space
        // next step is splitting it in pt space as well to estimate the yield differentially in pt
        for(Int_t j(0); j < ptBins; j++) {
            // get the integrated 
            Double_t integral(dud->IntegralAndError(j+1, j+2, integralError));
            dPtdPhi[j]->SetBinContent(i+1, integral);
            dPtdPhi[j]->SetBinError(i+1, integralError);
        }
        dud->Write();
        // save the current state of the unfolding object
        SaveConfiguration(unfoldedJetSpectrumIn ? kTRUE : kFALSE, kFALSE);
    }
    TF1* fourier = new TF1("fourier", "[0]*(1.+0.5*[1]*(TMath::Cos(2.*x)))", 0, TMath::Pi());
    TH1D* v2(new TH1D("v2FromFit", "v2FromFit", fBinsTrue->GetSize()-1, fBinsTrue->GetArray()));
    for(Int_t i(0); i < ptBins; i++) {
        dPtdPhi[i]->Fit(fourier, "VI");
        Style(gPad, "PEARSON");
        Style(dPtdPhi[i], kBlue, kDeltaPhi); 
        dPtdPhi[i]->DrawCopy();
        AliJetFlowTools::AddText(
                TString(Form("%.2f #LT p_{T} #LT %.2f", fBinsTrue->At(i), fBinsTrue->At(i+1))), 
                kBlack,
                .38,
                .56,
                .62,
                .65
                );
        v2->SetBinContent(1+i, fourier->GetParameter(1));
        v2->SetBinError(1+i, fourier->GetParError(1));
        Style(gPad, "PEARSON");
        Style(v2, kBlack, kV2, kTRUE);
        v2->DrawCopy();
        dPtdPhi[i]->Write();
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::ReplaceBins(TArrayI* array, TGraphErrors* graph) {
   // replace bins
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   Double_t x(0), y(0);
   graph->GetPoint(0, x, y);
   graph->SetPoint(array->At(0)-1, fBinsTrue->At(array->At(0)), y);
   graph->SetPointError(array->At(0)-1, 10, graph->GetErrorY(0));
   graph->SetPoint(array->At(1)-1, -5, -5);
}
//_____________________________________________________________________________
void AliJetFlowTools::ReplaceBins(TArrayI* array, TGraphAsymmErrors* graph) {
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
   // replace bins
   Double_t x(0), y(0);
   graph->GetPoint(0, x, y);
   graph->SetPoint(array->At(0)-1, fBinsTrue->At(array->At(0)), y);
   Double_t yl = graph->GetErrorYlow(0);
   Double_t yh = graph->GetErrorYhigh(0);
   graph->SetPointError(array->At(0)-1, 10, 10, yl, yh);
   graph->SetPoint(array->At(1)-1, -5, -5);
}
//_____________________________________________________________________________
void AliJetFlowTools::GetSignificance(
        TGraphErrors* n,                // points with stat error
        TGraphAsymmErrors* shape,       // points with shape error
        TGraphAsymmErrors* corr,        // points with stat error
        Int_t low,                      // lower bin (tgraph starts at 0)
        Int_t up                        // upper bin
        )
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // main use of this function is filling the static buffers
    Double_t statE(0), shapeE(0), corrE(0), y(0), x(0);

    // print some stuff
    printf(" double v2[] = {\n");
    Int_t iterator(0);
    for(Int_t i(low); i < up+1; i++) { 
        n->GetPoint(i, x, y);
        if(i==up) printf("%.4f}; \n\n", y);
        else printf("%.4f, \n", y);
        gV2->SetAt(y, iterator);
        iterator++;
    }
    iterator = 0;
    printf(" double stat[] = {\n");
    for(Int_t i(low); i < up+1; i++) { 
        y = n->GetErrorYlow(i);
        if(i==up) printf("%.4f}; \n\n", y);
        else printf("%.4f, \n", y);
        gStat->SetAt(y, iterator);
        iterator++;
    }
    iterator = 0;
    printf(" double shape[] = {\n");
    for(Int_t i(low); i < up+1; i++) { 
        y = shape->GetErrorYhigh(i);
        if (y < 0.0001) y = 0.0001;     // avoid crash in fit routine
        y*=gReductionFactor;
        shape->SetPointEYhigh(i, y);
        y = shape->GetErrorYlow(i);
        if ( y < 0.0001) y = 0.0001;
        y*=gReductionFactor;
        shape->SetPointEYlow(i, y);
        if(i==up) printf("%.4f}; \n\n", y);
        else printf("%.4f, \n", y);
        gShape->SetAt(y, iterator);
        iterator++;
    }
    iterator = 0;
    printf(" double corr[] = {\n");
    for(Int_t i(low); i < up+1; i++) { 
        y = gReductionFactorCorr*corr->GetErrorYhigh(i);
        corr->SetPointEYhigh(i, y);
        y = gReductionFactorCorr*corr->GetErrorYlow(i);
        corr->SetPointEYlow(i, y);
        if(i==up) printf("%.4f}; \n\n", y);
        else printf("%.4f, \n", y);
        gCorr->SetAt(y, iterator);
        iterator++;
    }

    // to plot the average error as function of number of events
    Float_t ctr(0);
    for(Int_t i(low); i < up+1; i++) {
        ctr = ctr + 1.;
        // set some flags to 0
        x = 0.;
        y = 0.;
        // get the nominal point
        n->GetPoint(i, x, y);
        statE += n->GetErrorY(i);
        shapeE += shape->GetErrorY(i);
        corrE += corr->GetErrorY(i);
        // combine the errors
    }
    printf(" ======================================\n");
    printf("  > between %i and %i GeV/c \n", low, up);
    cout << " AVERAGE_SHAPE " << shapeE/ctr << endl;
    cout << " AVERAGE_CORR " << corrE/ctr << endl;
    cout << " AVERAGE_STAT " << statE/ctr << endl;
}
//_____________________________________________________________________________
void AliJetFlowTools::MinimizeChi2nd()
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // Choose method upon creation between:
    // kMigrad, kSimplex, kCombined, 
    // kScan, kFumili
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.001);
 
    ROOT::Math::Functor f(&PhenixChi2nd,2); 
    double step[] = {0.0000001, 0.0000001};
    double variable[] = {-1., -1.};
        
    min.SetFunction(f);
    // Set the free variables to be minimized!
    min.SetVariable(0,"epsilon_c",variable[0], step[0]);
    min.SetVariable(1,"epsilon_b",variable[1], step[1]);


    min.Minimize(); 
}
//_____________________________________________________________________________
Double_t AliJetFlowTools::PhenixChi2nd(const Double_t *xx )
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // define arrays with results and errors here, see example at PhenixChi2()

    // return the function value at certain epsilon
    const Double_t epsc = xx[0];
    const Double_t epsb = xx[1];
    Double_t chi2(0);
    Int_t counts(gV2->GetSize() + gOffsetStop);

    if(gPwrtTo > -999) {
        // altered implemtation of eq 3 of arXiv:0801.1665v2
        // see analysis note and QM2014 poster for validation
        for(Int_t i(gOffsetStart); i < counts; i++) {

            // quadratic sum of statistical and uncorrelated systematic error
            Double_t e = gStat->At(i);

            // sum of v2 plus epsilon times correlated error minus hypothesis (gPwrtTo)
            // also the numerator of equation 3 of phenix paper
            Double_t numerator = TMath::Power(gV2->At(i)+epsc*gCorr->At(i)+epsb-gPwrtTo, 2);

            // modified denominator of equation 3 of phenix paper
            Double_t denominator = e*e;

            // add to the sum
            chi2 += numerator/denominator;
        }
        // add the square of epsilon to the total chi2 as penalty

        Double_t sumEpsb(0);
        for(Int_t j(gOffsetStart); j < counts; j++) sumEpsb += (epsb*epsb)/(gShape->At(j)*gShape->At(j));
        chi2 += epsc*epsc + sumEpsb/((float)counts);
    } else {
        // otherwise use the array
        for(Int_t i(gOffsetStart); i < counts; i++) {

            // quadratic sum of statistical and uncorrelated systematic error
            Double_t e = TMath::Sqrt(gStat->At(i)*gStat->At(i) + gPwrtToStatArray->At(i)*gPwrtToStatArray->At(i));

            // sum of v2 plus epsilon times correlated error minus hypothesis (gPwrtTo)
            // also the numerator of equation 3 of phenix paper
            Double_t numerator = TMath::Power(gV2->At(i)+epsc*gCorr->At(i)+epsb-gPwrtToArray->At(i), 2);

            // modified denominator of equation 3 of phenix paper
            Double_t denominator = e*e;

            // add to the sum
            chi2 += numerator/denominator;
        }
        // add the square of epsilon to the total chi2 as penalty

        Double_t sumEpsb(0);
        for(Int_t j(gOffsetStart); j < counts; j++) sumEpsb += (epsb*epsb)/(gShape->At(j)*gShape->At(j));
        chi2 += epsc*epsc + sumEpsb/((float)counts);
    }
    return chi2;
}
//_____________________________________________________________________________
Double_t AliJetFlowTools::ConstructFunctionnd(Double_t *x, Double_t *par)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // internal use only: evaluate the function at a given point
    if(par) SquelchWarning(); 
    return AliJetFlowTools::PhenixChi2nd(x);
}
//_____________________________________________________________________________
TF2* AliJetFlowTools::ReturnFunctionnd(Double_t &p)
{
#ifdef ALIJETFLOWTOOLS_DEBUG_FLAG
    printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
#endif
    // return the fitting function, pass the p-value w.r.t. 0 by reference
    const Int_t DOF(7-2);       // dof is n-2
    TF2 *f1 = new TF2("ndhist", AliJetFlowTools::ConstructFunctionnd, -100, 100, -100, 100, 0);
    printf(" > locating minima < \n");
    Double_t x(0), y(0);
    f1->GetMinimumXY(x, y);
    f1->GetXaxis()->SetTitle("#epsilon{c}");
    f1->GetXaxis()->SetTitle("#epsilon_{b}");
    f1->GetZaxis()->SetTitle("#chi^{2}");

    printf(" ===============================================================================\n");
    printf(" > minimal chi2 f(%.8f, %.8f) = %.8f  (i should be ok ... ) \n", x, y, f1->Eval(x, y));
    cout << "  so the probability of finding data at least as imcompatible with " << gPwrtTo << " as the actually" << endl;
    cout << "  observed data is " << TMath::Prob(f1->Eval(x, y), DOF) << endl; 
    cout << "  minimization parameters: EPSILON_B " << y << endl;
    cout << "                           EPSILON_C " << x << endl;
    printf(" ===============================================================================\n");

    // pass the p-value by reference and return the function
    p = TMath::Prob(f1->Eval(x, y), DOF);
    return f1;
}
//_____________________________________________________________________________
