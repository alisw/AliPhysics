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
//
// The task uses input from two analysis tasks:
// - $ALICE_ROOT/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskRhoVnModulation.cxx
//   used to retrieve jet spectra and delta pt distributions
// - $ALICE_ROOT/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskJetMatching.cxx
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
// to see an example of how to use this class, see $ALICE_ROOT/PWGCF/FLOW/macros/jetFlowTools.C

// root includes
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
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
    fResponseMaker      (new AliAnaChargedJetResponseMaker()),
    fPower              (new TF1("fPower","[0]*TMath::Power(x,-([1]))",0.,300.)),
    fSaveFull           (kFALSE),
    fActiveString       (""),
    fActiveDir          (0x0),
    fInputList          (0x0),
    fRefreshInput       (kTRUE),
    fOutputFileName     ("UnfoldedSpectra.root"),
    fOutputFile         (0x0),
    fCentralityBin      (0),
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
    fInOutUnfolding     (kTRUE),
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
    fFullResponseOut    (0x0) { // class constructor
    // create response maker weight function (tuned to PYTHIA spectrum)
    fResponseMaker->SetRMMergeWeightFunction(new TF1("weightFunction", "x*TMath::Power(1.+(1./(8.*0.9))*x, -8.)", 0 ,200));
}
//_____________________________________________________________________________
void AliJetFlowTools::Make() {
    // core function of the class
    // 1) rebin the raw output of the jet task to the desired binnings
    // 2) calls the unfolding routine
    // 3) writes output to file
    // can be repeated multiple times with different configurations

    // 1) manipulation of input histograms
    // check if the input variables are present
    if(fRefreshInput) {
        if(!PrepareForUnfolding()) {
            printf(" AliJetFlowTools::Make() Fatal error \n - couldn't prepare for unfolding ! \n");
            return;
        }
    }
    // 1a) resize the jet spectrum according to the binning scheme in fBinsTrue
    //     parts of the spectrum can end up in over or underflow bins
    TH1D* measuredJetSpectrumIn  = RebinTH1D(fSpectrumIn, fBinsRec, TString("resized_in_"), kFALSE);
    TH1D* measuredJetSpectrumOut = RebinTH1D(fSpectrumOut, fBinsRec, TString("resized_out_"), kFALSE);

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
    // select the unfolding method
    switch (fUnfoldingAlgorithm) {
        case kChi2 : {
            unfoldedJetSpectrumIn = UnfoldSpectrumChi2(       // do the inplane unfolding
                measuredJetSpectrumIn,
                resizedResponseIn,
                kinematicEfficiencyIn,
                measuredJetSpectrumTrueBinsIn,
                TString("in"),
                jetFindingEfficiency);
            printf(" > Spectrum (in plane) unfolded using kChi2 unfolding < \n");
        } break;
        case kBayesian : {
            unfoldedJetSpectrumIn = UnfoldSpectrumBayesian(       // do the inplane unfolding
                measuredJetSpectrumIn,
                resizedResponseIn,
                kinematicEfficiencyIn,
                measuredJetSpectrumTrueBinsIn,
                TString("in"),
                jetFindingEfficiency);
            printf(" > Spectrum (in plane) unfolded using kBayesian unfolding < \n");
        } break;
        case kBayesianAli : {
            unfoldedJetSpectrumIn = UnfoldSpectrumBayesianAli(       // do the inplane unfolding
                measuredJetSpectrumIn,
                resizedResponseIn,
                kinematicEfficiencyIn,
                measuredJetSpectrumTrueBinsIn,
                TString("in"),
                jetFindingEfficiency);
            printf(" > Spectrum (in plane) unfolded using kBayesianAli unfolding < \n");
        } break;
        case kSVD : {
            unfoldedJetSpectrumIn = UnfoldSpectrumSVD(       // do the inplane unfolding
                measuredJetSpectrumIn,
                resizedResponseIn,
                kinematicEfficiencyIn,
                measuredJetSpectrumTrueBinsIn,
                TString("in"),
                jetFindingEfficiency);
            printf(" > Spectrum (in plane) unfolded using kSVD unfolding < \n");
        } break;
        case kNone : {  // do nothing 
            resizedResponseIn->SetNameTitle("measuredSpectrumIn", "measured spectrum, in plane");
            unfoldedJetSpectrumIn = ProtectHeap(measuredJetSpectrumIn, kTRUE, TString("in"));
        } break;
        
        default : {
            printf(" > Selected unfolding method is not implemented yet ! \n");
            return;
        }
    }
    resizedResponseIn->SetNameTitle("ResponseMatrixIn", "response matrix in plane");
    resizedResponseIn->SetXTitle("p_{T}^{true} [GeV/c]");
    resizedResponseIn->SetYTitle("p_{T}^{rec} [GeV/c]");
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
    if(fInOutUnfolding) {
        TDirectoryFile* dirOut = new TDirectoryFile(Form("OutOfPlane___%s", fActiveString.Data()), Form("OutOfPlane___%s", fActiveString.Data()));
        dirOut->cd();
        switch (fUnfoldingAlgorithm) {
            case kChi2 : {
                unfoldedJetSpectrumOut = UnfoldSpectrumChi2(
                    measuredJetSpectrumOut,
                    resizedResponseOut,
                    kinematicEfficiencyOut,
                    measuredJetSpectrumTrueBinsOut,
                    TString("out"),
                    jetFindingEfficiency);
                printf(" > Spectrum (out of plane) unfolded using kChi2 < \n");
            } break;
            case kBayesian : {
                unfoldedJetSpectrumOut = UnfoldSpectrumBayesian(
                    measuredJetSpectrumOut,
                    resizedResponseOut,
                    kinematicEfficiencyOut,
                    measuredJetSpectrumTrueBinsOut,
                    TString("out"),
                    jetFindingEfficiency);
                printf(" > Spectrum (out of plane) unfolded using kBayesian < \n");
            } break;
            case kBayesianAli : {
                unfoldedJetSpectrumOut = UnfoldSpectrumBayesianAli(
                    measuredJetSpectrumOut,
                    resizedResponseOut,
                    kinematicEfficiencyOut,
                    measuredJetSpectrumTrueBinsOut,
                    TString("out"),
                    jetFindingEfficiency);
                printf(" > Spectrum (out of plane) unfolded using kBayesianAli < \n");
            } break;
            case kSVD : {
                unfoldedJetSpectrumOut = UnfoldSpectrumSVD(
                    measuredJetSpectrumOut,
                    resizedResponseOut,
                    kinematicEfficiencyOut,
                    measuredJetSpectrumTrueBinsOut,
                    TString("out"),
                    jetFindingEfficiency);
                printf(" > Spectrum (out of plane) unfolded using kSVD < \n");
            } break;
            case kNone : {  // do nothing
                resizedResponseOut->SetNameTitle("measuredSpectrumOut", "measured spectrum, out plane");
                unfoldedJetSpectrumOut = ProtectHeap(measuredJetSpectrumOut, kTRUE, TString("out"));
            } break;
            default : {
                printf(" > Selected unfolding method is not implemented yet ! \n");
                return;
            }
        }
        resizedResponseOut->SetNameTitle("ResponseMatrixOut", "response matrix in plane");
        resizedResponseOut->SetXTitle("p_{T}^{true} [GeV/c]");
        resizedResponseOut->SetYTitle("p_{T}^{rec} [GeV/c]");
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
                ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
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
            TGraphErrors* v2(GetV2((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_inv2"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_outv2")));
            if(v2) {
                v2->SetNameTitle("v2", "v_{2} from different in, out of plane yield");
                v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
                v2->GetYaxis()->SetTitle("v_{2}");
                v2 = ProtectHeap(v2);
                v2->Write();
            }
        } else if (unfoldedJetSpectrumOut && unfoldedJetSpectrumIn) {
            TGraphErrors* ratio(GetRatio((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_in"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_out"), TString(""), kTRUE, fBinsRec->At(fBinsRec->GetSize()-1)));
            if(ratio) {
                ratio->SetNameTitle("[NC]RatioInOutPlane", "[NC]Ratio in plane, out of plane jet spectrum");
                ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
                ratio->GetYaxis()->SetTitle("yield IN / yield OUT");
                ratio = ProtectHeap(ratio);
                ratio->Write();
            }
            TGraphErrors* v2(GetV2((TH1D*)unfoldedJetSpectrumIn->Clone("unfoldedLocal_inv2"), (TH1D*)unfoldedJetSpectrumOut->Clone("unfoldedLocal_outv2")));
             if(v2) {
                v2->SetNameTitle("v2", "v_{2} from different in, out of plane yield");
                v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
                v2->GetYaxis()->SetTitle("v_{2}");
                v2 = ProtectHeap(v2);
                v2->Write();
            }
        }
    }   // end of if(fInOutUnfolding)
    fDeltaPtDeltaPhi->Write();
    fJetPtDeltaPhi->Write();
    // save the current state of the unfolding object
    SaveConfiguration(unfoldedJetSpectrumIn ? kTRUE : kFALSE, unfoldedJetSpectrumOut ? kTRUE : kFALSE);
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::UnfoldSpectrumChi2(
        const TH1D* measuredJetSpectrum,      // truncated raw jets (same binning as pt rec of response) 
        const TH2D* resizedResponse,           // response matrix
        const TH1D* kinematicEfficiency,      // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,        // unfolding template: same binning is pt gen of response
        const TString suffix,                 // suffix (in or out of plane)
        const TH1D* jetFindingEfficiency)     // jet finding efficiency (optional)
{
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
    if(status == 0 && gMinuit->fISW[1] == 3) {
        // if the unfolding converged and the hessian matrix is reliable, plot the pearson coefficients
        TVirtualFitter *fitter(TVirtualFitter::GetFitter());
        if(gMinuit) gMinuit->Command("SET COV");
        TMatrixD covarianceMatrix(fBinsTrue->GetSize()-1, fBinsTrue->GetSize()-1, fitter->GetCovarianceMatrix());
        TMatrixD *pearson((TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix));
        pearson->Print();
        TH2D *hPearson(new TH2D(*pearson));
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients, %s plane", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
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
    if(pearson) {
        TH2D* hPearson(new TH2D(*pearson));
        pearson->Print();
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients_%s", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
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
    responseMatrixLocalTransposePrior->SetXTitle("p_{T}^{true} [GeV/c]");
    responseMatrixLocalTransposePrior->SetYTitle("p_{T}^{rec} [GeV/c]");
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
    if(status == 0 && gMinuit->fISW[1] == 3) {
        // if the unfolding converged and the hessian matrix is reliable, plot the pearson coefficients
        TVirtualFitter *fitter(TVirtualFitter::GetFitter());
        if(gMinuit) gMinuit->Command("SET COV");
        TMatrixD covarianceMatrix(fBinsTrue->GetSize()-1, fBinsTrue->GetSize()-1, fitter->GetCovarianceMatrix());
        TMatrixD *pearson((TMatrixD*)CalculatePearsonCoefficients(&covarianceMatrix));
        pearson->Print();
        TH2D *hPearson(new TH2D(*pearson));
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients, %s plane", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
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
TH1D* AliJetFlowTools::UnfoldSpectrumBayesian(
        const TH1D* measuredJetSpectrum,              // jet pt in pt rec bins 
        const TH2D* resizedResponse,                  // full response matrix, normalized in slides of pt true
        const TH1D* kinematicEfficiency,              // kinematic efficiency
        const TH1D* measuredJetSpectrumTrueBins,      // jet pt in pt true bins, also the prior when measured is chosen as prior
        const TString suffix,                         // suffix (in, out)
        const TH1D* jetFindingEfficiency)             // jet finding efficiency (optional)
{
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
    if(pearson) {
        TH2D* hPearson(new TH2D(*pearson));
        pearson->Print();
        hPearson->SetNameTitle(Form("PearsonCoefficients_%s", suffix.Data()), Form("Pearson coefficients_%s", suffix.Data()));
        hPearson = ProtectHeap(hPearson);
        hPearson->Write();
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
    responseMatrixLocalTransposePrior->SetXTitle("p_{T}^{true} [GeV/c]");
    responseMatrixLocalTransposePrior->SetYTitle("p_{T}^{rec} [GeV/c]");
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
Bool_t AliJetFlowTools::PrepareForUnfolding()
{
    // prepare for unfolding
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
    if(!fRMSSpectrumIn && fInOutUnfolding) { // initialie the profiles which will hold the RMS values. if binning changes in between unfolding
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
    TString spectrumName(Form("fHistJetPsi2Pt_%i", fCentralityBin));
    fJetPtDeltaPhi = ((TH2D*)fInputList->FindObject(spectrumName.Data()));
    if(!fJetPtDeltaPhi) {
        printf(" Couldn't find spectrum %s ! \n", spectrumName.Data());
        return kFALSE;
    }
    fJetPtDeltaPhi = ProtectHeap(fJetPtDeltaPhi, kFALSE);
    // in plane spectrum
    if(!fInOutUnfolding) {
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
    // normalize spectra to event count if requested
    if(fNormalizeSpectra) {
        TH1* rho((TH1*)fInputList->FindObject(Form("fHistRho_%i", fCentralityBin)));
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
    TString deltaptName(Form("fHistDeltaPtDeltaPhi2ExLJ_%i", fCentralityBin));
    fDeltaPtDeltaPhi = ((TH2D*)fInputList->FindObject(deltaptName.Data()));
    if(!fDeltaPtDeltaPhi) {
        printf(" Couldn't find delta pt matrix %s ! \n", deltaptName.Data());
        printf(" > may be ok, depending no what you want to do < \n");
        fRefreshInput = kTRUE;
        return kTRUE;
    }
    fDeltaPtDeltaPhi = ProtectHeap(fDeltaPtDeltaPhi, kFALSE);
    // in plane delta pt distribution
    if(!fInOutUnfolding) {
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
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityBin), Form("dpt_response_INPLANE_%i", fCentralityBin));
    fDptIn->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");
    fDptIn = ProtectHeap(fDptIn);
    fDptOut = new TH2D(*rfOut);
    fDptOut->SetNameTitle(Form("dpt_response_OUTOFPLANE_%i", fCentralityBin), Form("dpt_response_OUTOFPLANE_%i", fCentralityBin));
    fDptOut->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptOut->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");
    fDptOut = ProtectHeap(fDptOut);
    
    fRefreshInput = kTRUE;     // force cloning of the input
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
    // 1) get a prior for unfolding. 
    // this can be either an unfolded spectrum from e.g. chi2 unfolding or the measured spectrum
    TH1D* unfolded(0x0);
    TDirectoryFile* dirOut = new TDirectoryFile(Form("Prior_%s___%s", suffix.Data(), fActiveString.Data()), Form("Prior_%s___%s", suffix.Data(), fActiveString.Data()));
    dirOut->cd();
    switch (fPrior) {    // select the prior for unfolding
        case kPriorChi2 : {
            if(fBinsTruePrior && fBinsRecPrior) {       // if set, use different binning for the prior
                TArrayD* tempArrayTrue(fBinsTrue);      // temporarily cache the original binning
                fBinsTrue = fBinsTruePrior;             // switch binning schemes (will be used in UnfoldSpectrumChi2())
                TArrayD* tempArrayRec(fBinsRec);
                fBinsRec = fBinsRecPrior;
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
                    return 0x0;
                }
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
TH2D* AliJetFlowTools::NormalizeTH2D(TH2D* histo, Bool_t noError) {
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
            if(noError) histo->SetBinError(  1+i, j+1, 0.);
            else histo->SetBinError(  1+i, j+1, histo->GetBinError(  1+i, j+1)/weight);
        }
    }
    return histo;
}
//_____________________________________________________________________________
TH1D* AliJetFlowTools::RebinTH1D(TH1D* histo, TArrayD* bins, TString suffix, Bool_t kill) {
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
    TH1D* rebinned = new TH1D(name.Data(), name.Data(), bins->GetSize()-1, bins->GetArray());
    for(Int_t i(0); i < histo->GetXaxis()->GetNbins(); i++) {
        // loop over the bins of the old histo and fill the new one with its data
        rebinned->Fill(histo->GetBinCenter(i+1), histo->GetBinContent(i+1));
    }
    if(kill) delete histo;
    return rebinned;
}
//_____________________________________________________________________________
TH2D* AliJetFlowTools::RebinTH2D(TH2D* rebinMe, TArrayD* binsTrue, TArrayD* binsRec, TString suffix) {
    // weighted rebinning of a th2d, implementation for function call to AliAnaChargedJetResponseMaker
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
    // normalize a th1d to a certain scale
    histo->Sumw2();
    Double_t integral = histo->Integral()*scale;
    if (integral > 0 && scale == 1.) histo->Scale(1./integral, "width");
    else if (scale != 1.) histo->Scale(1./scale, "width");
    else printf(" > Histogram integral < 0, cannot normalize \n");
    return histo;
}
//_____________________________________________________________________________
TMatrixD* AliJetFlowTools::CalculatePearsonCoefficients(TMatrixD* covarianceMatrix) 
{
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
    // smoothen the spectrum using a user defined function
    // returns a clone of the original spectrum if fitting failed
    // if kill is kTRUE the input spectrum will be deleted from the heap
    // if 'count' is selected, bins are filled with integers (necessary if the 
    // histogram is interpreted in a routine which accepts only counts)
    if(!spectrum || !function) return 0x0;
    if(start > max) printf(" > cannot extrapolate fit beyond fit range ! < " );
    TH1D* temp = (TH1D*)spectrum->Clone(Form("%s_smoothened", spectrum->GetName()));
    temp->Sumw2();      // if already called on the original, this will give off a warning but do nothing
    TFitResultPtr r = temp->Fit(function, "WLIS", "", min, max);
    if((int)r == 0) {   // MINUIT status
        for(Int_t i(0); i < temp->GetNbinsX() + 1; i++) {
            if(temp->GetBinCenter(i) > start) {     // from this pt value use extrapolation
                if(counts) temp->SetBinContent(i, (int)function->Integral(temp->GetXaxis()->GetBinLowEdge(i),temp->GetXaxis()->GetBinUpEdge(i))/temp->GetXaxis()->GetBinWidth(i));
                else temp->SetBinContent(i, function->Integral(temp->GetXaxis()->GetBinLowEdge(i),temp->GetXaxis()->GetBinUpEdge(i))/temp->GetXaxis()->GetBinWidth(i));
                if(temp->GetBinContent(i) > 0) temp->SetBinError(i, TMath::Sqrt(temp->GetBinContent(i)));
            }
        }
    }
    if(kill) delete spectrum;
    return temp;
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TCanvas* c, TString style)
{
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
void AliJetFlowTools::Style(TVirtualPad* c, TString style)
{
    // set a default style for a canvas
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
    } else printf(" > Style called with unknown option %s \n    returning < \n", style.Data());
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TLegend* l)
{
    // set a default style for a legend
    l->SetTextSize(.06);
    l->SetFillColor(0);
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TH1* h, EColor col, histoType type)
{
    // style a histo
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2.);
    h->SetMarkerSize(2.);
    h->SetTitleSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetXaxis()->SetTitleOffset(.9);
    switch (type) {
        case kInPlaneSpectrum : {
            h->SetTitle("IN PLANE");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
        } break;
        case kOutPlaneSpectrum : {
            h->SetTitle("OUT OF PLANE");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kUnfoldedSpectrum : {
            h->SetTitle("UNFOLDED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kFoldedSpectrum : {
            h->SetTitle("FOLDED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kMeasuredSpectrum : {
            h->SetTitle("MEASURED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       default : break;
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::Style(TGraph* h, EColor col, histoType type)
{
    // style a histo
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2.);
    h->SetMarkerSize(2.);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetXaxis()->SetTitleOffset(.9);
    switch (type) {
        case kInPlaneSpectrum : {
            h->SetTitle("IN PLANE");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
        } break;
        case kOutPlaneSpectrum : {
            h->SetTitle("OUT OF PLANE");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kUnfoldedSpectrum : {
            h->SetTitle("UNFOLDED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kFoldedSpectrum : {
            h->SetTitle("FOLDED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       case kMeasuredSpectrum : {
            h->SetTitle("MEASURED");
            h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
            h->GetYaxis()->SetTitle("counts");
       } break;
       default : break;
    }
}
//_____________________________________________________________________________
void AliJetFlowTools::PostProcess(TString def, TString in, TString out, Int_t columns) const
{
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
   if(fInOutUnfolding) canvasOut = new TCanvas("PearsonOut", "PearsonOut");
   TCanvas* canvasRatioMeasuredRefoldedIn(new TCanvas("RefoldedIn", "RefoldedIn"));
   TCanvas* canvasRatioMeasuredRefoldedOut(0x0);
   if(fInOutUnfolding) canvasRatioMeasuredRefoldedOut = new TCanvas("RefoldedOut", "RefoldedOut");
   TCanvas* canvasSpectraIn(new TCanvas("SpectraIn", "SpectraIn")); 
   TCanvas* canvasSpectraOut(0x0);
   if(fInOutUnfolding) canvasSpectraOut = new TCanvas("SpectraOut", "SpectraOut");
   TCanvas* canvasRatio(0x0);
   if(fInOutUnfolding) canvasRatio = new TCanvas("Ratio", "Ratio");
   TCanvas* canvasV2(0x0);
   if(fInOutUnfolding) canvasV2 = new TCanvas("V2", "V2");
   TCanvas* canvasMISC(new TCanvas("MISC", "MISC"));
   TCanvas* canvasMasterIn(new TCanvas("defaultIn", "defaultIn"));
   TCanvas* canvasMasterOut(0x0);
   if(fInOutUnfolding) canvasMasterOut = new TCanvas("defaultOut", "defaultOut");
   (fInOutUnfolding) ? canvasMISC->Divide(4, 2) : canvasMISC->Divide(4, 1);
   TDirectoryFile* defDir(0x0);
   
   // get an estimate of the number of outputs and find the default set
   Int_t cacheMe(0);
   for(Int_t i(0); i < listOfKeys->GetSize(); i++) {
       if(dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()))) {
           if(!strcmp(listOfKeys->At(i)->GetName(), def.Data())) defDir = dynamic_cast<TDirectoryFile*>(readMe.Get(listOfKeys->At(i)->GetName()));
           cacheMe++;
       }
   }
   Int_t rows(TMath::Floor(cacheMe/(float)columns)+((cacheMe%4)>0));
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
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(3);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(4);
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
                   temp->GetXaxis()->SetTitle("p_{T} [GeV/c]");
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
                   eff->DrawCopy();
               } else if(rm && eff) {
                   Style(rm);
                   Style(eff);
                   canvasMISC->cd(7);
                   Style(gPad, "PEARSON");
                   rm->DrawCopy("colz");
                   canvasMISC->cd(8);
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
                   temp->GetXaxis()->SetTitle("p_{T} [GeV/c]");
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
           if(ratioYield) {
               Style(ratioYield);
               ratioYield->GetYaxis()->SetRangeUser(-1., 3.);
               ratioYield->SetFillColor(kRed);
               ratioYield->Draw("ap");
           }
           canvasV2->cd(j);
           TGraphErrors* ratioV2((TGraphErrors*)tempDir->Get(Form("v2_%s", dirName.Data())));
           if(ratioV2) {
               Style(ratioV2);
               ratioV2->GetYaxis()->SetRangeUser(-.25, .75);
               ratioV2->SetFillColor(kRed);
               ratioV2->Draw("ap");
           }
       }
   }
   TFile output(out.Data(), "RECREATE");
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
Bool_t AliJetFlowTools::SetRawInput (
        TH2D* detectorResponse,  // detector response matrix
        TH1D* jetPtIn,           // in plane jet spectrum
        TH1D* jetPtOut,          // out of plane jet spectrum
        TH1D* dptIn,             // in plane delta pt distribution
        TH1D* dptOut,            // out of plane delta pt distribution
        Int_t eventCount) {
    // set input histograms manually
    fDetectorResponse   = detectorResponse;
    fSpectrumIn         = jetPtIn;
    fSpectrumOut        = jetPtOut;
    fDptInDist          = dptIn;
    fDptOutDist         = dptOut;
    fRawInputProvided   = kTRUE;
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
    fDptIn->SetNameTitle(Form("dpt_response_INPLANE_%i", fCentralityBin), Form("dpt_response_INPLANE_%i", fCentralityBin));
    fDptIn->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptIn->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");
    fDptOut = ConstructDPtResponseFromTH1D(fDptOutDist, fAvoidRoundingError);
    fDptOut->SetNameTitle(Form("dpt_response_OUTOFPLANE_%i", fCentralityBin), Form("dpt_response_OUTOFPLANE_%i", fCentralityBin));
    fDptOut->GetXaxis()->SetTitle("p_{T}^{gen} [GeV/c]");
    fDptOut->GetYaxis()->SetTitle("p_{T}^{rec} [GeV/c]");
    
    return kTRUE;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::GetRatio(TH1 *h1, TH1* h2, TString name, Bool_t appendFit, Int_t xmax) 
{
    // return ratio of h1 / h2
    // histograms can have different binning. errors are propagated as uncorrelated
    if(!(h1 && h2) ) {
        printf(" GetRatio called with NULL argument(s) \n ");
        return 0x0;
    }
    Int_t j(0);
    TGraphErrors *gr = new TGraphErrors();
    gr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    Double_t binCent(0.), ratio(0.), error2(0.), binWidth(0.);
    for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
        binCent = h1->GetXaxis()->GetBinCenter(i);
        if(xmax > 0. && binCent > xmax) continue;
        j = h2->FindBin(binCent);
        binWidth = h1->GetXaxis()->GetBinWidth(i);
        if(h2->GetBinContent(j) > 0.) {
            ratio = h1->GetBinContent(i)/h2->GetBinContent(j);
            /* original propagation of uncertainty
            Double_t A = 1./h2->GetBinContent(j)*h1->GetBinError(i);
            Double_t B = 0.;
            if(h2->GetBinError(j)>0.) {
                B = -1.*h1->GetBinContent(i)/(h2->GetBinContent(j)*h2->GetBinContent(j))*h2->GetBinError(j);
                error2 = A*A + B*B;
            } else error2 = A*A;        */
            Double_t A = h1->GetBinError(i)/h1->GetBinContent(i);
            Double_t B = h2->GetBinError(i)/h2->GetBinContent(i);
            error2 = ratio*ratio*A*A+ratio*ratio*B*B;
            if(error2 > 0 ) error2 = TMath::Sqrt(error2);
            gr->SetPoint(gr->GetN(),binCent,ratio);
            gr->SetPointError(gr->GetN()-1,0.5*binWidth,error2);
        }
    }
    if(appendFit) {
        TF1* fit(new TF1("lin", "pol0", 10, 100));
        gr->Fit(fit);
    }
    if(strcmp(name, "")) gr->SetNameTitle(name.Data(), name.Data());
    return gr;
}
//_____________________________________________________________________________
TGraphErrors* AliJetFlowTools::GetV2(TH1 *h1, TH1* h2, Double_t r, TString name) 
{
    // get v2 from difference of in plane, out of plane yield
    // h1 must hold the in-plane yield, h2 holds the out of plane  yield
    // different binning is allowed but will mean that the error
    // propagation is unreliable
    // r is the event plane resolution for the chosen centrality
    if(!(h1 && h2) ) {
        printf(" GetV2 called with NULL argument(s) \n ");
        return 0x0;
    }
    Int_t j(0);
    TGraphErrors *gr = new TGraphErrors();
    gr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    Float_t binCent(0.), ratio(0.), error2(0.), binWidth(0.);
    Double_t pre(TMath::Pi()/(4.*r)), in(0.), out(0.), ein(0.), eout(0.);
    for(Int_t i(1); i <= h1->GetNbinsX(); i++) {
        binCent = h1->GetXaxis()->GetBinCenter(i);
        j = h2->FindBin(binCent);
        binWidth = h1->GetXaxis()->GetBinWidth(i);
        if(h2->GetBinContent(j) > 0.) {
            in = h1->GetBinContent(i);
            ein = h1->GetBinError(i);
            out = h2->GetBinContent(j);
            eout = h2->GetBinError(j);
            ratio = pre*((in-out)/(in+out));
            error2 =((r*4.)/(TMath::Pi()))*((4.*out*out/(TMath::Power(in+out, 4)))*ein*ein+(4.*in*in/(TMath::Power(in+out, 4)))*eout*eout);
            if(error2 > 0) error2 = TMath::Sqrt(error2);
            gr->SetPoint(gr->GetN(),binCent,ratio);
            gr->SetPointError(gr->GetN()-1,0.5*binWidth,error2);
        }
    }
    if(strcmp(name, "")) gr->SetNameTitle(name.Data(), name.Data());
    return gr;
}
//_____________________________________________________________________________
void AliJetFlowTools::WriteObject(TObject* object, TString suffix, Bool_t kill) {
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
    // save configuration parameters to histogram
    TH1F* summary = new TH1F("UnfoldingConfiguration","UnfoldingConfiguration", 20, -.5, 19.5);
    summary->SetBinContent(1, fBetaIn);
    summary->GetXaxis()->SetBinLabel(1, "fBetaIn");
    summary->SetBinContent(2, fBetaOut);
    summary->GetXaxis()->SetBinLabel(2, "fBetaOut");
    summary->SetBinContent(3, fCentralityBin);
    summary->GetXaxis()->SetBinLabel(3, "fCentralityBin");
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
TH1D* AliJetFlowTools::ProtectHeap(TH1D* protect, Bool_t kill, TString suffix) {
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
TH2D* AliJetFlowTools::ProtectHeap(TH2D* protect, Bool_t kill, TString suffix) {
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
TGraphErrors* AliJetFlowTools::ProtectHeap(TGraphErrors* protect, Bool_t kill, TString suffix) {
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
