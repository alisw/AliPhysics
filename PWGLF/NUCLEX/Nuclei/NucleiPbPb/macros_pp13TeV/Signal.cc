#include "src/Common.h"
#include "src/FitModules.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <memory>
#include <functional>
using std::function;
#include <utility>
using std::pair;
#include <vector>
using std::vector;

#include <TAxis.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;

using namespace RooFit;

enum RangeVar{
  kWidenRange,
  kShiftRange
};

void BinCountingVariations(RangeVar rv, FitModule& sigFunction, FitModule& bkgFunction, RooRealVar& sig_norm_var, RooRealVar& bkg_norm_var, TH1* dat, bool isTPC, TDirectory* dir, int iBin, int iS, int iC, const char* list_name, TH1* hRawCounts, TH1* hChi2, TH1* hSignificance, TH1* hVarSyst, TH1* hVarSystJoined, TH1D** hVarCounts);

void reset_fit_parameters(FitExpExpTailGaus& sigTOF, FitExpExpTailGaus& bkgTOF, FitExpTailTailGaus& sigTPC, bool isMB);

void Signal(bool useMBsignal=false, bool use_extended=true, bool isMC = false) {

  gStyle->SetOptStat(0);

  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel=kError; // Suppressing warning outputs

  /// Taking all the histograms from the MC file
  const char* input_name = (isMC) ? kMCfilename.data() : kDataFilename.data();
  const char* input_name_MB = (isMC) ? kMCfilenameMB.data() : kDataFilenameMB.data();
  const char* output_name = (isMC) ? kSignalMCOutput.data() : kSignalOutput.data();
  TFile input_file(input_name);
  TFile input_file_MB(input_name_MB);
  TFile output_file(output_name,"recreate");

  string kRefListNames = (isMC) ? kFilterListNamesMCasData.data() : kFilterListNames.data();

  int iList=0;

  for (auto list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kRefListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;

    if(string(list_key->GetName())!=kRefListNames.data() && isMC){
      printf("list_key: %s found: %s\n",list_key->GetName(),kRefListNames.data());
      continue;
    }

    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    Requires(list, list_key->GetName());
    string base_list = list_key->GetName();
    string list_name_MB = base_list;
    if(!isMC) list_name_MB.insert(kRefListNames.size()-1,"MB");
    TTList* listMB = (TTList*)input_file_MB.Get(list_name_MB.data());
    Requires(listMB,list_name_MB.data());
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Taking all the necessary histogram to perform the analysis
    TH3F *fATOFsignal = (TH3F*)list->Get("fATOFsignal");
    TH3F *fMTOFsignal = (TH3F*)list->Get("fMTOFsignal");
    TH3F *fATPCcounts = (TH3F*)list->Get("fATPCcounts");
    TH3F *fMTPCcounts = (TH3F*)list->Get("fMTPCcounts");
    // MB histograms
    TH3F *fATOFsignalMB = (TH3F*)listMB->Get("fATOFsignal");
    TH3F *fMTOFsignalMB = (TH3F*)listMB->Get("fMTOFsignal");
    TH3F *fATPCcountsMB = (TH3F*)listMB->Get("fATPCcounts");
    TH3F *fMTPCcountsMB = (TH3F*)listMB->Get("fMTPCcounts");

    /// Taking information about centrality bins
    const int n_centralities = fATOFsignal->GetNbinsX();
    auto cent_labels = *(fATOFsignal->GetXaxis()->GetXbins());

    /// Now it comes a bit of a complication. To improve fit quality
    /// and to minimize the manual interventions one should do the fits to
    /// the same pT bin for all the species and all the centrality classes
    /// before moving to another pT bin. This requires to build all the
    /// arrays before starting the actual analysis.

    /// Build arrays to analyse both deuteron and anti-deuterons
    /// with the same code
    TH3F* tof_histo[2] = {fMTOFsignal,fATOFsignal};
    TH3F* tpc_histo[2] = {fMTPCcounts,fATPCcounts};

    TH3F* tof_histoMB[2] = {fMTOFsignalMB,fATOFsignalMB};
    TH3F* tpc_histoMB[2] = {fMTPCcountsMB,fATPCcountsMB};

    /// Build arrays to analyse all the centrality classes for
    /// both deuteron and anti-deuterons. Complicate stuff just for fun.
    TH1D* hRawCounts[2][kCentLength];
    TH1D* hSignalTailGaus[2][kCentLength];
    TH1D* hSystFit[2][kCentLength];
    TH1D* hSignificance[2][kCentLength];
    TH1D* hChiSquare[2][kCentLength];
    TH1D* hChiSquareTPC[2][kCentLength];
    TH1D* hRawCountsTPC[2][kCentLength];
    TH1D* hSignificanceTPC[2][kCentLength];

    TCanvas* cComparison[2][kCentLength];
    TH1D* hRatioTPC[2][kCentLength];

    int kNewGreen = kGreen + 3;
    int color_vector[] = {kBlack,kBlue,kNewGreen,kOrange,kRed};
    TH1D* hWidenRangeCounts[2][kCentLength][kNsigmaVar];
    TH1D* hShiftRangeCounts[2][kCentLength][kNshiftVar];
    TH1D* hWidenRangeCountsTPC[2][kCentLength][kNsigmaVar];
    TH1D* hShiftRangeCountsTPC[2][kCentLength][kNshiftVar];
    TH1D* hWidenRangeSyst[2][kCentLength];
    TH1D* hShiftRangeSyst[2][kCentLength];
    TH1D* hWidenRangeSystTPC[2][kCentLength];
    TH1D* hShiftRangeSystTPC[2][kCentLength];
    TH1D* hWidenRangeSystJoined[2][kCentLength];
    TH1D* hShiftRangeSystJoined[2][kCentLength];

    /// Creating the directories to be used to store the results
    for (int iS = 1; iS>=0; iS--) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      TDirectory* sig_dir= dir->mkdir("Fits");
      TDirectory *root_sig_dir = dir->mkdir("root_TOF");
      TDirectory* side_dir = dir->mkdir("Sidebands");
      TDirectory* TPCOnly_dir = dir->mkdir("TPConly");
      TDirectory* significance_dir = dir->mkdir("Significance");
      TDirectory* significance_TPC_dir = dir->mkdir("SignificanceTPC");
      TDirectory* chi_dir = dir->mkdir("ChiSquare");
      TDirectory* chiTPC_dir = dir->mkdir("ChiSquareTPC");
      if(string(list_key->GetName())==kFilterListNames.data()){
        TDirectory* syst_dir = dir->mkdir("Systematic");
        TDirectory* range_width_dir = dir->mkdir("Range_width");
        TDirectory* range_shift_dir = dir->mkdir("Range_shift");
        TDirectory* range_width_TPC_dir = dir->mkdir("Range_width_TPC");
        TDirectory* range_shift_TPC_dir = dir->mkdir("Range_shift_TPC");
        for(int iC = 0; iC < kCentLength; iC++){
          range_width_dir->mkdir(Form("C_%d",iC));
          range_shift_dir->mkdir(Form("C_%d",iC));
          range_width_TPC_dir->mkdir(Form("C_%d",iC));
          range_shift_TPC_dir->mkdir(Form("C_%d",iC));
          syst_dir->mkdir(Form("C_%d",iC));
        }
      }
      for(int iC = 0; iC < kCentLength; iC++){
        sig_dir->mkdir(Form("C_%d",iC));
        side_dir->mkdir(Form("C_%d",iC));
        TPCOnly_dir->mkdir(Form("C_%d",iC));
        root_sig_dir->mkdir(Form("C_%d", iC));
      }
    }

    for (int iS = 1; iS>=0; iS--) {
      for (int iC = 0; iC < kCentLength; ++iC) {
        hRawCounts[iS][iC] = new TH1D(Form("hRawCounts%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RawCounts",kNPtBins,kPtBins);
        plotting::SetHistStyle(hRawCounts[iS][iC],plotting::kListColors[iList]);
        hRawCountsTPC[iS][iC] = new TH1D(Form("hRawCountsTPC%c%i",kLetter[iS],iC),";p_{T} GeV/c; TPC raw counts",kNPtBins,kPtBins);
        plotting::SetHistStyle(hRawCountsTPC[iS][iC],plotting::kListColors[iList]);
        hSignificance[iS][iC] = new TH1D(Form("hSignificance_%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}",kNPtBins,kPtBins);
        hSignificanceTPC[iS][iC] = new TH1D(Form("hSignificanceTPC_%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}",kNPtBins,kPtBins);
        hChiSquare[iS][iC] = new TH1D(Form("hChiSquare%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #chi^{2}/NDF",kNPtBins,kPtBins);
        hChiSquareTPC[iS][iC] = new TH1D(Form("hChiSquareTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #chi^{2}/NDF",kNPtBins,kPtBins);
        hSignalTailGaus[iS][iC] = new TH1D(Form("hSignalTailGaus%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RawCounts",kNPtBins,kPtBins);
        plotting::SetHistStyle(hSignalTailGaus[iS][iC],kBlack);
        // Bin-counting related histograms
        if(string(list_key->GetName())==kFilterListNames.data()){
          hWidenRangeSyst[iS][iC] = new TH1D(Form("hWidenRangeSyst%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          hWidenRangeSystTPC[iS][iC] = new TH1D(Form("hWidenRangeSystTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          hShiftRangeSyst[iS][iC] = new TH1D(Form("hShiftRangeSyst%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          hShiftRangeSystTPC[iS][iC] = new TH1D(Form("hShiftRangeSystTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          hWidenRangeSystJoined[iS][iC] = new TH1D(Form("hWidenRangeSystJoined%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          hShiftRangeSystJoined[iS][iC] = new TH1D(Form("hShiftRangeSystJoined%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",kNPtBins,kPtBins);
          for(int iCut=0; iCut<kNsigmaVar; iCut++){
            hWidenRangeCounts[iS][iC][iCut] = new TH1D(Form("hWidenRangeCounts%c%i_%i",kLetter[iS],iC,iCut),"; p_{T}(GeV/c); raw counts",kNPtBins,kPtBins);
            hWidenRangeCountsTPC[iS][iC][iCut] = new TH1D(Form("hWidenRangeCountsTPC%c%i_%i",kLetter[iS],iC,iCut),"; p_{T}(GeV/c); raw counts",kNPtBins,kPtBins);
          }
          for(int iCut=0; iCut<kNshiftVar; iCut++){
            hShiftRangeCounts[iS][iC][iCut] = new TH1D(Form("hShiftRangeCounts%c%i_%i",kLetter[iS],iC,iCut),"; p_{T}(GeV/c); raw counts",kNPtBins,kPtBins);
            hShiftRangeCountsTPC[iS][iC][iCut] = new TH1D(Form("hShiftRangeCountsTPC%c%i_%i",kLetter[iS],iC,iCut),"; p_{T}(GeV/c); raw counts",kNPtBins,kPtBins);
          }
        }
      }
    } 

    for (int iB = 1; iB <= kNPtBins; ++iB) {
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (bin_center < kPtRange[0] || bin_center > kPtRange[1]) continue;
      for (int iS = 1; iS>=0; iS--) {
        /// Setting up the fitting environment for TOF analysis
        float low_pt_safe = (bin_center<kTOFlowPt) ? 1.04 : 0.;
        RooRealVar m("dm2", "m^{2} - m^{2}_{d}", kLimitsTOF[0], kLimitsTOF[1]+low_pt_safe, "GeV^{2}/#it{c}^{4}");
        m.setBins(1000,"cache");
        m.setRange("Full", kLimitsTOF[0],kLimitsTOF[1]);

        FitExpExpTailGaus fModelTOF(&m,use_extended);
        fModelTOF.mMu->setUnit("GeV^{2}/#it{c}^{4}");
        fModelTOF.mSigma->setUnit("GeV^{2}/#it{c}^{4}");
        fModelTOF.mTau0->setUnit("GeV^{-2}#it{c}^{4}");
        fModelTOF.mTau1->setUnit("GeV^{-2}#it{c}^{4}");

        //Background
        RooRealVar m_bis("dm2_bis","m^{2} - m^{2}_{d}",kLimitsTOFbkg[0],kLimitsTOFbkg[1],"GeV^{2}/#it{c}^{4}");
        m_bis.setBins(1000,"cache");
        m_bis.setRange("Full", kLimitsTOFbkg[0],kLimitsTOFbkg[1]);
        FitExpExpTailGaus fBkgTOF(&m_bis,use_extended);
        fBkgTOF.UseSignal(false);
        fBkgTOF.mTau0->setUnit("GeV^{-2}#it{c}^{4}");
        fBkgTOF.mTau1->setUnit("GeV^{-2}#it{c}^{4}");

        // Setting up the fitting environment for the TPC analysis
        RooRealVar ns("ns", "n#sigma_{d}", kLimitsTPC[0], kLimitsTPC[1], "a. u.");
        ns.setBins(1000,"cache");
        ns.setRange("Full", kLimitsTPC[0], kLimitsTPC[1]);

        // TPC analysis
        FitExpTailTailGaus fModelTPC(&ns);
        fModelTPC.mSigma->setUnit("a. u.");
        fModelTPC.mMu->setUnit("a. u.");
        fModelTPC.mAlpha0->setUnit("a. u.");
        fModelTPC.mAlpha1->setUnit("a. u.");

        for (int iC = kCentLength; iC--;)
        {
          // Reset fit parameters
          if(useMBsignal){
            reset_fit_parameters(fModelTOF,fBkgTOF,fModelTPC,iC==kCentLength-1);
          } else {
            reset_fit_parameters(fModelTOF,fBkgTOF,fModelTPC,true);
          }
          // TOF analysis
          if(bin_center > kTOFminPt){
            if(bin_center > kCentPtLimits[iC]) continue;
            TString iTitle = Form("%s, Class %s, %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", kSymbols[iS], kRomanLabels[iC], kPtBins[iB-1], kPtBins[iB]);
            TString iName = Form("d%i_%i", iC, iB);
            int iBin = tof_histo[iS]->GetYaxis()->FindBin(bin_center);
            TH1D *dat = (iC==9) ? tof_histoMB[iS]->ProjectionZ(Form("data%i_%i", iC, iB), 1, 13, iBin, iBin) : tof_histo[iS]->ProjectionZ(Form("data%i_%i", iC, iB), kCentBinsArray[iC][0], kCentBinsArray[iC][1], iBin, iBin);
            /// Fits
            base_dir->cd(Form("%s/Fits/C_%d",kNames[iS].data(),iC));
            if (bin_center < kNoBkgTOF) {
              fModelTOF.UseBackground(false);
              fBkgTOF.UseBackground(false);
            } else {
              fModelTOF.UseBackground(true);
              fBkgTOF.UseBackground(true);
              if (bin_center < kSingleExpBkg){
                fModelTOF.mKbkg->setVal(0.);
                fModelTOF.mKbkg->setConstant(true);
                fModelTOF.mTau0->setConstant(true);
                fModelTOF.mTau1->setRange(-10,0.2);
                fModelTOF.mTau1->setVal(-0.5);
              }
              else{
                fModelTOF.mKbkg->setVal(0.5);
                fModelTOF.mKbkg->setConstant(false);
                fModelTOF.mTau0->setConstant(false);
              }
              if (bin_center < kSingleExpSideBkg){
                fBkgTOF.mKbkg->setVal(0.);
                fBkgTOF.mKbkg->setConstant(true);
                fBkgTOF.mTau0->setConstant(true);
                fBkgTOF.mTau1->setRange(-10,0.2);
                fBkgTOF.mTau1->setVal(-0.5);
              }
              else{
                fBkgTOF.mKbkg->setVal(0.5);
                fBkgTOF.mKbkg->setConstant(false);
                fBkgTOF.mTau0->setConstant(false);
              }
            }
            float first_mu = 0.;
            float first_sigma = 0.;
            float second_mu = 0.;
            float second_sigma = 0.;
            // Control with normal histograms
            dat->SetMarkerStyle(20);
            dat->SetMarkerSize(0.7);
            dat->SetTitle(Form("%1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", kPtBins[iB - 1], kPtBins[iB]));
            dat->SetOption("PEx0");
            TF1 *func = nullptr;
            RooPlot *TOFmodelPlot;
            // printf("list name: %s iS: %d iC: %d iB: %d\n",list_key->GetName(),iS,iC,iB);
            if(bin_center < kTOFlowPt){
              fModelTOF.mMu->setConstant(false);
              fModelTOF.mSigma->setConstant(false);
              fModelTOF.mAlpha0->setConstant(false);
              // fModelTOF.mAlpha1->setConstant(false);
              fModelTOF.mAlpha0->setRange(0.5,4.);
              fModelTOF.mAlpha0->setVal(1.1);
              // fModelTOF.mAlpha1->setRange(-4.,-0.8);
              // fModelTOF.mAlpha1->setVal(-1.1);
              TOFmodelPlot = fModelTOF.FitData(dat, iName, iTitle, "Full", "Full");
              first_mu = fModelTOF.mMu->getVal();
              int mu_bin = dat->FindBin(first_mu);
              float right_edge_gaus_fit = dat->GetBinLowEdge(mu_bin+1+3);
              float left_edge_gaus_fit = dat->GetBinLowEdge(mu_bin-3);
              func = new TF1(Form("gaus_fit_%d_%d", iC, iB), "gaus",left_edge_gaus_fit,right_edge_gaus_fit);
              func->SetParLimits(1, -0.1, 0.7);
              func->SetParLimits(2, 0.09, 0.4);
              dat->Fit(Form("gaus_fit_%d_%d", iC, iB), "RQ");
              second_mu = func->GetParameter(1);
              second_sigma = func->GetParameter(2);
              fModelTOF.mMu->setVal(second_mu);
              fModelTOF.mMu->setConstant(true);
              fModelTOF.mSigma->setVal(second_sigma);
              fModelTOF.mSigma->setConstant(true);
              m.setRange(Form("Fit_%i_%c",iB,kLetter[iS]), fModelTOF.mMu->getVal()-3.*fModelTOF.mSigma->getVal(), fModelTOF.mMu->getVal()+5.*fModelTOF.mSigma->getVal());
              TOFmodelPlot = fModelTOF.FitData(dat, iName, iTitle, Form("Fit_%i_%c",iB,kLetter[iS]), Form("Fit_%i_%c",iB,kLetter[iS]));
              fModelTOF.mMu->setConstant(false);
              fModelTOF.mSigma->setConstant(false);
            }
            else{
              if(useMBsignal){
                if(iC==kCentLength-1){
                  fModelTOF.mMu->setConstant(false);
                  fModelTOF.mSigma->setConstant(false);
                  fModelTOF.mAlpha0->setConstant(false);
                }
              }
              TOFmodelPlot = fModelTOF.FitData(dat, iName, iTitle, "Full", "Full");
              if(bin_center < kNoBkgTOF){
                // printf("#1) fModelTOF.mMu: %f second_mu: %f fModelTOF.mSigma: %f second_sigma: %f\n",fModelTOF.mMu->getVal(),second_mu,fModelTOF.mSigma->getVal(),second_sigma);
                first_mu = fModelTOF.mMu->getVal();
                int mu_bin = dat->FindBin(first_mu);
                float right_edge_gaus_fit = dat->GetBinLowEdge(mu_bin+1+1);
                float left_edge_gaus_fit = dat->GetBinLowEdge(mu_bin-3);
                func = new TF1(Form("gaus_fit_%d_%d", iC, iB), "gaus",left_edge_gaus_fit,right_edge_gaus_fit);
                func->SetParLimits(1, -0.1, 0.7);
                func->SetParLimits(2, 0.09, 0.4);
                dat->Fit(Form("gaus_fit_%d_%d", iC, iB), "RQ");
                second_mu = func->GetParameter(1);
                second_sigma = func->GetParameter(2);
                fModelTOF.mMu->setVal(second_mu);
                fModelTOF.mMu->setConstant(true);
                fModelTOF.mSigma->setVal(second_sigma);
                fModelTOF.mSigma->setConstant(true);
                m.setRange(Form("Fit_%i_%c",iB,kLetter[iS]),
                fModelTOF.mMu->getVal()-3.*fModelTOF.mSigma->getVal(), fModelTOF.mMu->getVal()+5.*fModelTOF.mSigma->getVal());
                fModelTOF.mMu->setVal(second_mu);
                fModelTOF.mSigma->setVal(second_sigma);
                // printf("#2) fModelTOF.mMu: %f second_mu: %f fModelTOF.mSigma: %f second_sigma: %f\n",fModelTOF.mMu->getVal(),second_mu,fModelTOF.mSigma->getVal(),second_sigma);
                if(useMBsignal){
                  reset_fit_parameters(fModelTOF,fBkgTOF,fModelTPC,iC==kCentLength-1);
                } else {
                  reset_fit_parameters(fModelTOF,fBkgTOF,fModelTPC,true);
                }
                fModelTOF.UseBackground(false);
                fModelTOF.mMu->setVal(second_mu);
                fModelTOF.mSigma->setVal(second_sigma);
                // printf("#3) fModelTOF.mMu: %f second_mu: %f fModelTOF.mSigma: %f second_sigma: %f\n",fModelTOF.mMu->getVal(),second_mu,fModelTOF.mSigma->getVal(),second_sigma);
                TOFmodelPlot = fModelTOF.FitData(dat, iName, iTitle, Form("Fit_%i_%c",iB,kLetter[iS]), Form("Fit_%i_%c",iB,kLetter[iS]));
                fModelTOF.mMu->setConstant(false);
                fModelTOF.mSigma->setConstant(false);
                // printf("#4) fModelTOF.mMu: %f second_mu: %f fModelTOF.mSigma: %f second_sigma: %f\n",fModelTOF.mMu->getVal(),second_mu,fModelTOF.mSigma->getVal(),second_sigma);
              } else if (bin_center<kFixSigma){
                first_mu = fModelTOF.mMu->getVal();
                int mu_bin = dat->FindBin(first_mu);
                float right_edge_gaus_fit = dat->GetBinLowEdge(mu_bin+1+1);
                float left_edge_gaus_fit = dat->GetBinLowEdge(mu_bin-3);
                func = new TF1(Form("gaus_fit_%d_%d", iC, iB), "gaus",left_edge_gaus_fit,right_edge_gaus_fit);
                func->SetParLimits(1, -0.1, 0.7);
                func->SetParLimits(2, 0.09, 0.4);
                dat->Fit(Form("gaus_fit_%d_%d", iC, iB), "RQ");
                second_mu = func->GetParameter(1);
                second_sigma = func->GetParameter(2);
                fModelTOF.mMu->setVal(second_mu);
                fModelTOF.mMu->setConstant(true);
                fModelTOF.mSigma->setVal(second_sigma);
                fModelTOF.mSigma->setConstant(true);
                TOFmodelPlot = fModelTOF.FitData(dat, iName, iTitle, "Full", "Full");
                fModelTOF.mMu->setConstant(false);
                fModelTOF.mSigma->setConstant(true);
              }
              if(useMBsignal){
                if(iC==kCentLength-1){
                  fModelTOF.mMu->setConstant(true);
                  fModelTOF.mSigma->setConstant(true);
                  fModelTOF.mAlpha0->setConstant(true);
                  // fModelTOF.mAlpha1->setConstant(true);
                }
              }
            } 
            //fModelTOF.mSigma->setConstant(false);
            if(bin_center > kTOFminPt){
              base_dir->cd(Form("%s/Fits/C_%d", kNames[iS].data(), iC));
              TOFmodelPlot->Write();
              base_dir->cd(Form("%s/root_TOF/C_%d", kNames[iS].data(), iC));
              dat->Write();
            }
            if(use_extended){
              hSignalTailGaus[iS][iC]->SetBinContent(iB,fModelTOF.mSigCounts->getVal());
              hSignalTailGaus[iS][iC]->SetBinError(iB,fModelTOF.mSigCounts->getError());
            }
            else{
              hSignalTailGaus[iS][iC]->SetBinContent(iB,fModelTOF.mFraction->getVal()*fModelTOF.mNentries);
              hSignalTailGaus[iS][iC]->SetBinError(iB,fModelTOF.mFraction->getError()*fModelTOF.mNentries);
            } 
            if(bin_center > kNoBkgTOF && bin_center > kBinCountingCut) {
              if(use_extended){
                hRawCounts[iS][iC]->SetBinContent(iB,fModelTOF.mSigCounts->getVal());
                hRawCounts[iS][iC]->SetBinError(iB,fModelTOF.mSigCounts->getError());
              }
              else{
                hRawCounts[iS][iC]->SetBinContent(iB,fModelTOF.mFraction->getVal()*fModelTOF.mNentries);
                hRawCounts[iS][iC]->SetBinError(iB,fModelTOF.mFraction->getError()*fModelTOF.mNentries);
              }
              hChiSquare[iS][iC]->SetBinContent(iB, fModelTOF.mChi2);
              hChiSquare[iS][iC]->SetBinError(iB, 0);
            } 
            /// TOF bin counting: range widening
            BinCountingVariations(RangeVar::kWidenRange,fModelTOF,fBkgTOF,m,m_bis,dat,false,base_dir,iB,iS,iC,list_key->GetName(),
              hRawCounts[iS][iC],hChiSquare[iS][iC],hSignificance[iS][iC],hWidenRangeSyst[iS][iC],hWidenRangeSystJoined[iS][iC],hWidenRangeCounts[iS][iC]
            );
            /// TOF bin counting: range shifting
            BinCountingVariations(RangeVar::kShiftRange,fModelTOF,fBkgTOF,m,m_bis,dat,false,base_dir,iB,iS,iC,list_key->GetName(),
              hRawCounts[iS][iC],hChiSquare[iS][iC],hSignificance[iS][iC],hShiftRangeSyst[iS][iC],hShiftRangeSystJoined[iS][iC],hShiftRangeCounts[iS][iC]
            );
          }

          /// TPC analysis
          if(bin_center < kTPCmaxPt){
            TString iTitle = Form("%s, Class %s, %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", kSymbols[iS], kRomanLabels[iC], kPtBins[iB-1], kPtBins[iB]);
            TString iName = Form("TPC_d%i_%i", iC, iB);
            int iBin_tpc = tpc_histo[iS]->GetYaxis()->FindBin(bin_center);
            TH1D *tpc_dat = (iC==9) ? tpc_histoMB[iS]->ProjectionZ(Form("tpc_data%i_%i", iC, iB), 1, 13, iBin_tpc, iBin_tpc) : tpc_histo[iS]->ProjectionZ(Form("tpc_data%i_%i", iC, iB), kCentBinsArray[iC][0], kCentBinsArray[iC][1], iBin_tpc, iBin_tpc);
            tpc_dat->Rebin(4);
            /// Fits
            base_dir->cd(Form("%s/TPConly/C_%i",kNames[iS].data(),iC));
            if(bin_center < kNoBkgTPC) {
              fModelTPC.UseBackground(false);
            } else {
              fModelTPC.UseBackground(true);
            }
            RooPlot *expTailTailGausPlot;
            if(useMBsignal){
              if(iC==kCentLength-1){
                fModelTPC.mMu->setConstant(false);
                fModelTPC.mSigma->setConstant(false);
                fModelTPC.mAlpha0->setConstant(false);
                fModelTPC.mAlpha1->setConstant(false);
              }
            }
            expTailTailGausPlot = fModelTPC.FitData(tpc_dat, iName, iTitle, "Full", "Full");
            if(useMBsignal){
              if(iC==kCentLength-1){
                fModelTPC.mMu->setConstant(true);
                fModelTPC.mSigma->setConstant(true);
                fModelTPC.mAlpha0->setConstant(true);
                fModelTPC.mAlpha1->setConstant(true);
              }
            }
            expTailTailGausPlot->Write();

            hRawCountsTPC[iS][iC]->SetBinContent(iB, fModelTPC.mSigCounts->getVal());
            hRawCountsTPC[iS][iC]->SetBinError(iB, fModelTPC.mSigCounts->getError());
            hChiSquareTPC[iS][iC]->SetBinContent(iB,fModelTPC.mChi2);
            hChiSquareTPC[iS][iC]->SetBinError(iB,0.);
            
            /// TPC bin counting: range widening
            BinCountingVariations(RangeVar::kWidenRange,fModelTPC,fModelTPC,ns,ns,tpc_dat,true,base_dir,iB,iS,iC,list_key->GetName(),
              hRawCountsTPC[iS][iC],hChiSquareTPC[iS][iC],hSignificanceTPC[iS][iC],hWidenRangeSystTPC[iS][iC],hWidenRangeSystJoined[iS][iC],hWidenRangeCountsTPC[iS][iC]
            );
            /// TPC bin counting: range shifting
            BinCountingVariations(RangeVar::kShiftRange,fModelTPC,fModelTPC,ns,ns,tpc_dat,true,base_dir,iB,iS,iC,list_key->GetName(),
              hRawCountsTPC[iS][iC],hChiSquareTPC[iS][iC],hSignificanceTPC[iS][iC],hShiftRangeSystTPC[iS][iC],hShiftRangeSystJoined[iS][iC],hShiftRangeCountsTPC[iS][iC]
            );
          }
        } // End loop iC
      } // End loop iS
    } // End loop iB

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < kCentLength; ++iC) {
        base_dir->cd(Form("%s/Fits",kNames[iS].data()));
        hRawCounts[iS][iC]->Write();
        hSignalTailGaus[iS][iC]->Write();
        if(string(list_key->GetName())==kFilterListNames.data()){
          base_dir->cd(Form("%s/Systematic/C_%i",kNames[iS].data(),iC));
          hShiftRangeSyst[iS][iC]->Write();
          hWidenRangeSyst[iS][iC]->Write();
          hWidenRangeSystTPC[iS][iC]->Write();
          hShiftRangeSystTPC[iS][iC]->Write();
          hWidenRangeSystJoined[iS][iC]->Write();
          hShiftRangeSystJoined[iS][iC]->Write();
          base_dir->cd(Form("%s/Range_width/C_%i",kNames[iS].data(),iC));
          for(int iCut=0; iCut<kNsigmaVar; iCut++ ){
            plotting::SetHistStyle(hWidenRangeCounts[iS][iC][iCut],kRed,21);
            hWidenRangeCounts[iS][iC][iCut]->Write();
          }
          base_dir->cd(Form("%s/Range_shift/C_%i",kNames[iS].data(),iC));
          for(int iCut=0; iCut<kNshiftVar; iCut++ ){
            plotting::SetHistStyle(hShiftRangeCounts[iS][iC][iCut],kRed,21);
            hShiftRangeCounts[iS][iC][iCut]->Write();
          }
          base_dir->cd(Form("%s/Range_width_TPC/C_%i",kNames[iS].data(),iC));
          for(int iCut=0; iCut<kNsigmaVar; iCut++ ){
            plotting::SetHistStyle(hWidenRangeCountsTPC[iS][iC][iCut],kRed,21);
            hWidenRangeCountsTPC[iS][iC][iCut]->Write();
          }
          base_dir->cd(Form("%s/Range_shift_TPC/C_%i",kNames[iS].data(),iC));
          for(int iCut=0; iCut<kNshiftVar; iCut++ ){
            plotting::SetHistStyle(hShiftRangeCountsTPC[iS][iC][iCut],kRed,21);
            hShiftRangeCountsTPC[iS][iC][iCut]->Write();
          }
          //
          cComparison[iS][iC] = new TCanvas(Form("cComparison_%c_%i",kLetter[iS],iC),Form("cComparison_%c_%i",kLetter[iS],iC));
          
          cComparison[iS][iC]->cd();
          hSignalTailGaus[iS][iC]->Draw("PE");
          hWidenRangeCounts[iS][iC][0]->Draw("PE SAME");
          TLegend leg(0.65,0.65, 0.90, 0.82,Form("%s, Class %s",kNames[iS].data(),kRomanLabels[iC]),"brNDC");
          leg.SetBorderSize(0);
          leg.AddEntry(hSignalTailGaus[iS][iC],"Fit","pe");
          leg.AddEntry(hWidenRangeCounts[iS][iC][0],"Bin counting","pe");
          leg.Draw();
          base_dir->cd(Form("%s/Fits",kNames[iS].data()));
          cComparison[iS][iC]->Write();
          if(iS==0 && iC==0){
            cComparison[iS][iC]->Print(Form("%s/plots/comparison.pdf[",kBaseOutputDir.data()),"pdf");
          }
          //} else if (!(iS==0 && iC==0)) {
            cComparison[iS][iC]->Print(Form("%s/plots/comparison.pdf",kBaseOutputDir.data()),"pdf");
          //} else {
          if(iS==1 && iC==kCentLength-1){
            cComparison[iS][iC]->Print(Form("%s/plots/comparison.pdf]",kBaseOutputDir.data()),"pdf");
          }
        }
        base_dir->cd(Form("%s/Significance",kNames[iS].data()));
        hSignificance[iS][iC]->Write();
        base_dir->cd(Form("%s/SignificanceTPC",kNames[iS].data()));
        hSignificanceTPC[iS][iC]->Write();
        base_dir->cd(Form("%s/TPConly",kNames[iS].data()));
        hRawCountsTPC[iS][iC]->Write();
        base_dir->cd(Form("%s/ChiSquare",kNames[iS].data()));
        hChiSquare[iS][iC]->Write();
        base_dir->cd(Form("%s/ChiSquareTPC",kNames[iS].data()));
        hChiSquareTPC[iS][iC]->Write();
      }
    }
    base_dir->Close();
    iList++;
  } 
  output_file.Close();
}

void BinCountingVariations(RangeVar rv, FitModule& sigFunction, FitModule& bkgFunction, RooRealVar& sig_norm_var, RooRealVar& bkg_norm_var, TH1* dat, bool isTPC, TDirectory* dir, int iBin, int iS, int iC, const char* list_name, TH1* hRawCounts, TH1* hChi2, TH1* hSignificance, TH1* hVarSyst, TH1* hVarSystJoined, TH1D** hVarCounts){
  float bin_center = (kPtBins[iBin]+kPtBins[iBin-1])/2;
  if (!isTPC && sigFunction.mSigma->getVal()>0.39995) printf("PAY ATTENTION! mSigma to the limit for list %s at cent %i and bin_center %f (TOF)\n",list_name,iC,bin_center);
  if (isTPC && sigFunction.mSigma->getVal()>1.1995) printf("PAY ATTENTION! mSigma to the limit for list %s at cent %i and bin_center %f (TPC)\n",list_name,iC,bin_center);
  int simga_limit = 0;
  std::vector<float> residual_vector;
  if(rv == RangeVar::kWidenRange) simga_limit = kNsigmaVar;
  else if (rv == RangeVar::kShiftRange) simga_limit = kNshiftVar;
  else{
    std::cout<<"Wrong value of bin counting variations. ABORT." << std::endl;
    exit(1);
  }
  float bin_width = dat->GetBinLowEdge(2) - dat->GetBinLowEdge(1);
  for(int iSigma=0; iSigma < simga_limit; iSigma++){
    float left_sigma = sigFunction.mMu->getVal() - 3*sigFunction.mSigma->getVal();
    float right_sigma = (isTPC) ? (sigFunction.mMu->getVal() + 3*sigFunction.mSigma->getVal()) : (sigFunction.mMu->getVal() + 5*sigFunction.mSigma->getVal());
    int left_edge_bin = dat->FindBin(left_sigma);
    int right_edge_bin = dat->FindBin(right_sigma);
    if(rv == RangeVar::kWidenRange){
      left_edge_bin -= vSigmaWidth[iSigma];
      right_edge_bin += vSigmaWidth[iSigma];
    } else {
      left_edge_bin += vSigmaShift[iSigma];
      right_edge_bin += vSigmaShift[iSigma];
    }
    float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
    float left_edge_float_fit_bkg = left_edge_float;
    if(!isTPC && left_edge_float < -0.33){ 
      left_edge_float_fit_bkg = -0.32;
    }
    float right_edge_float = dat->GetBinLowEdge(right_edge_bin+1);
    float right_edge_float_fit_bkg = right_edge_float;
    if(!isTPC && right_edge_float > 0.73){
      right_edge_float_fit_bkg =  0.72;
    }
    float low_pt_safe = (bin_center<kTOFlowPt) ? 1. : 0.;
    if (isTPC) {
        if (left_edge_float<kLimitsTPC[0]){
          printf("WARNING: left limit for bin-counting below histogram range!\n TPC, list: %s iB: %d iS: %d iC %d \n value: %f limit: %f\n Study: %d iSigma: %d\n",list_name,iBin,iS,iC,left_edge_float,kLimitsTPC[0],rv,iSigma);
        }
        if (right_edge_float>kLimitsTPC[1]){
          printf("WARNING: right limit for bin-counting above histogram range!\n TPC, list: %s iB: %d iS: %d iC %d \n value: %f limit: %f\n Study: %d iSigma: %d\n",list_name,iBin,iS,iC,right_edge_float,kLimitsTPC[1],rv,iSigma);
        }
    } else {
        if (left_edge_float<kLimitsTOFbkg[0]){
          printf("WARNING: left limit for bin-counting below histogram range!\n TOF, list: %s iB: %d iS: %d iC %d \n value: %f limit: %f\n Study: %d iSigma: %d\n",list_name,iBin,iS,iC,left_edge_float,kLimitsTOFbkg[0],rv,iSigma);
        }
        if (right_edge_float>kLimitsTOFbkg[1]+low_pt_safe){
          printf("WARNING: right limit for bin-counting above histogram range!\n TOF, list: %s iB: %d iS: %d iC %d \n value: %f limit: %f\n Study: %d iSigma: %d\n",list_name,iBin,iS,iC,right_edge_float,kLimitsTOFbkg[1],rv,iSigma);
        }
    }
    bkgFunction.mX->setRange("signal",left_edge_float,right_edge_float);
    sigFunction.mX->setRange("signal",left_edge_float,right_edge_float);
    if (iSigma==0) {
      TString iTitle = Form("%s, Class %s, %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", kSymbols[iS], kRomanLabels[iC], kPtBins[iBin-1], kPtBins[iBin]);
      TString iName = Form("d%i_%i", iC, iBin);
      if(!isTPC && rv==RangeVar::kWidenRange){
        bkgFunction.mX->setRange("left",kLimitsTOFbkg[0],left_edge_float_fit_bkg);
        bkgFunction.mX->setRange("right",right_edge_float_fit_bkg,kLimitsTOFbkg[1]);
        dir->cd(Form("%s/Sidebands/C_%d",kNames[iS].data(),iC));
        RooPlot *bkgPlot; 
        if(bin_center > kNoBkgTOF) {
            bkgPlot = bkgFunction.FitData(dat, Form("%s_sideband",iName.Data()), iTitle, "left,right","Full",true);
            bkgPlot->Write();
        }
      }
    }
    float bkg_integral = ((!isTPC && bin_center > kNoBkgTOF) || (isTPC && bin_center > kNoBkgTPC)) ? bkgFunction.mBackground->createIntegral(bkg_norm_var, NormSet(bkg_norm_var), Range("signal"))->getVal() * bkgFunction.mBkgCounts->getVal() : 0;
    float tot_integral = dat->Integral(left_edge_bin,right_edge_bin);
    float sig_integral = tot_integral - bkg_integral;
    float sig_err = TMath::Sqrt(tot_integral+bkg_integral);
    if(string(list_name)==kFilterListNames.data()){
      hVarCounts[iSigma]->SetBinContent(iBin, sig_integral);
      hVarCounts[iSigma]->SetBinError(iBin, sig_err);
    }
    if(rv==RangeVar::kWidenRange && iSigma==0){
      if(!isTPC){
        if(bin_center < kNoBkgTOF || bin_center < kBinCountingCut){
          hChi2->SetBinContent(iBin, bkgFunction.mChi2);
          hChi2->SetBinError(iBin, 0);
          hRawCounts->SetBinContent(iBin, sig_integral);
          hRawCounts->SetBinError(iBin, sig_err);
          hSignificance->SetBinContent(iBin, sig_integral / TMath::Sqrt(tot_integral));
          hSignificance->SetBinError(iBin,0);
        } else {
          sig_integral = sigFunction.mSignal->createIntegral(sig_norm_var, NormSet(sig_norm_var), Range("signal"))->getVal() * sigFunction.mSigCounts->getVal();
          bkg_integral = sigFunction.mBackground->createIntegral(sig_norm_var, NormSet(sig_norm_var), Range("signal"))->getVal() * sigFunction.mBkgCounts->getVal();
          hChi2->SetBinContent(iBin, sigFunction.mChi2);
          hChi2->SetBinError(iBin, 0);
          hSignificance->SetBinContent(iBin, sig_integral / TMath::Sqrt(sig_integral+bkg_integral));
          hSignificance->SetBinError(iBin,0);
        }
      } else {
        float sig_integral_fit = sigFunction.mSignal->createIntegral(sig_norm_var, NormSet(sig_norm_var), Range("signal"))->getVal() * sigFunction.mSigCounts->getVal();
        float bkg_integral_fit = sigFunction.mBackground->createIntegral(sig_norm_var, NormSet(sig_norm_var), Range("signal"))->getVal() * sigFunction.mBkgCounts->getVal();
        hSignificance->SetBinContent(iBin, sig_integral_fit / TMath::Sqrt(sig_integral_fit+bkg_integral_fit));
        hSignificance->SetBinError(iBin,0);
      }
    }
    if(iSigma==0 && rv==RangeVar::kWidenRange) continue;
    if(string(list_name)!=kFilterListNames.data()) break; 
    residual_vector.push_back(sig_integral);
  }
  if(string(list_name)!=kFilterListNames.data()) return;
  float var_range_syst = TMath::RMS(residual_vector.size(),residual_vector.data());
  var_range_syst /= hRawCounts->GetBinContent(iBin);
  hVarSyst->SetBinContent(iBin, var_range_syst);
  if(isTPC){
    if(bin_center<1.) hVarSystJoined->SetBinContent(iBin,var_range_syst);
  } else {
    if(bin_center>1.) hVarSystJoined->SetBinContent(iBin,var_range_syst);
  }
}

void reset_fit_parameters(FitExpExpTailGaus &sigTOF, FitExpExpTailGaus &bkgTOF, FitExpTailTailGaus &sigTPC, bool isMB)
{
  // TOF signal function
  if(isMB){
    sigTOF.mMu->setRange(0.00001,0.7);
    sigTOF.mMu->setVal(0.1);
    sigTOF.mSigma->setRange(0.09,0.30);
    sigTOF.mSigma->setVal(0.2);
    // sigTOF.mAlpha1->setRange(-4.,-0.8);
    // sigTOF.mAlpha1->setVal(-1.3);
    sigTOF.mAlpha0->setRange(0.8,4.);
    sigTOF.mAlpha0->setVal(1.5);
  }
  sigTOF.mSigCounts->setRange(9.,100000.);
  sigTOF.mSigCounts->setVal(50000);
  sigTOF.mBkgCounts->setRange(0.,1.e10);
  sigTOF.mBkgCounts->setVal(20000);
  sigTOF.mTau0->setRange(-10.,-0.8);
  sigTOF.mTau0->setVal(-5.);
  sigTOF.mTau1->setRange(-0.8,0.);
  sigTOF.mTau1->setVal(-0.3);
  sigTOF.mKbkg->setRange(0.,1.);
  sigTOF.mKbkg->setVal(0.5);

  // TOF background function
  bkgTOF.mBkgCounts->setRange(0.,1.e10);
  bkgTOF.mBkgCounts->setVal(20000);
  bkgTOF.mTau0->setRange(-10.,-0.8);
  bkgTOF.mTau0->setVal(-5.);
  bkgTOF.mTau1->setRange(-0.8,0.);
  bkgTOF.mTau1->setVal(-0.3);
  bkgTOF.mKbkg->setRange(0.,1.);
  bkgTOF.mKbkg->setVal(0.5);

  // TPC signal function
  if(isMB){
    sigTPC.mSigma->setRange(0.2,1.2);
    sigTPC.mSigma->setVal(0.75);
    sigTPC.mMu->setRange(-1.,1.);
    sigTPC.mMu->setVal(0.);
    sigTPC.mAlpha0->setRange(1.,3.);
    sigTPC.mAlpha0->setVal(1.1);
    sigTPC.mAlpha1->setRange(-3.,-1.);
    sigTPC.mAlpha1->setVal(-1.1);
  }
  sigTPC.mSigCounts->setRange(9.,1000000.);
  sigTPC.mSigCounts->setVal(50000);
  sigTPC.mBkgCounts->setRange(0.,1.e10);
  sigTPC.mBkgCounts->setVal(20000);
  sigTPC.mTau0->setRange(-10.,-0.01);
  sigTPC.mTau0->setVal(-5.);
}