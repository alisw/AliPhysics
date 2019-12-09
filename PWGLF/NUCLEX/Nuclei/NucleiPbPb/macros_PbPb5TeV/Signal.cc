#include "src/Common.h"
#include "src/FitModules.h"
#include "src/Utils.h"
using namespace utils;

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
#include <TStyle.h>
#include <TError.h>

#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooMsgService.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

using namespace RooFit;

void Signal() {

  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel=kError; // Suppressing warning outputs

  /// Taking all the histograms from the MC file
  TFile input_file(kDataFilename.data());
  TFile output_file(kSignalOutput.data(),"recreate");

  /// Setting up the fitting environment for TOF analysis
  RooRealVar m("dm2","m^{2} - m^{2}_{d}",kFitminPt,kFitmaxPt,"GeV^{2}/#it{c}^{4}");
  m.setBins(1000,"cache");
  m.setRange("Full", kFitminPt,kFitmaxPt);

  FitExpExpTailGaus fExpExpTailGaus(&m);
  fExpExpTailGaus.mMu->setRange(0.00001,0.5);
  fExpExpTailGaus.mMu->setVal(0.1);
  fExpExpTailGaus.mMu->setUnit("GeV^{2}/#it{c}^{4}");
  fExpExpTailGaus.mSigma->setRange(0.05,0.15);
  fExpExpTailGaus.mSigma->setVal(0.1);
  fExpExpTailGaus.mSigma->setUnit("GeV^{2}/#it{c}^{4}");
  fExpExpTailGaus.mAlpha0->setRange(1.1,3.);
  fExpExpTailGaus.mAlpha0->setVal(1.2);
  fExpExpTailGaus.mAlpha0->setUnit("GeV^{2}/#it{c}^{4}");
  fExpExpTailGaus.mSigCounts->setRange(0.,kFitmaxNBkg);
  fExpExpTailGaus.mTau0->setUnit("GeV^{-2}#it{c}^{4}");
  fExpExpTailGaus.mTau1->setUnit("GeV^{-2}#it{c}^{4}");

  //Background
  RooRealVar m_bis("dm2_bis","m^{2} - m^{2}_{d}",kFitminPt,kFitmaxPt,"GeV^{2}/#it{c}^{4}");
  m_bis.setBins(1000,"cache");
  m_bis.setRange("Full", kFitminPt,kFitmaxPt);
  FitExpExpTailGaus fBkg(&m_bis);
  fBkg.UseSignal(false);
  fBkg.mTau0->setUnit("GeV^{-2}#it{c}^{4}");
  fBkg.mTau1->setUnit("GeV^{-2}#it{c}^{4}");

  // Low pTbins
  FitExpPolDSCrystalBall fExpPolDSCrystalBall(&m);
  fExpPolDSCrystalBall.mMu->setRange(0.00001,0.5);
  fExpPolDSCrystalBall.mMu->setVal(0.1);
  fExpPolDSCrystalBall.mMu->setUnit("GeV^{2}/#it{c}^{4}");
  fExpPolDSCrystalBall.mSigma->setRange(0.05,0.15);
  fExpPolDSCrystalBall.mSigma->setVal(0.1);
  fExpPolDSCrystalBall.mSigma->setUnit("GeV^{2}/#it{c}^{4}");
  fExpPolDSCrystalBall.mSigCounts->setRange(0.,kFitmaxNBkg);
  fExpPolDSCrystalBall.mTau0->setUnit("GeV^{-2}#it{c}^{4}");
  fExpPolDSCrystalBall.mTau1->setUnit("GeV^{-2}#it{c}^{4}");
  

  // Setting up the fitting environment for the TPC analysis
  RooRealVar ns("ns","n#sigma_{d}",-3.,3,"a. u.");
  ns.setBins(1000,"cache");
  ns.setRange("Full", -3., 3.);
  ns.setRange("Special", -3., 3.);

  // TPC analysis
  FitGausGaus fGausGaus(&ns);
  fGausGaus.mSigma->setRange(0.2,0.9);
  fGausGaus.mSigma->setVal(0.75);
  fGausGaus.mSigma->setUnit("a. u.");
  fGausGaus.mMu->setRange(-0.5,0.5);
  fGausGaus.mMu->setUnit("a. u.");
  fGausGaus.mMuBkg->setRange(-10.,-2.);
  fGausGaus.mMuBkg->setVal(-7);
  fGausGaus.mMuBkg->setUnit("a. u.");
  fGausGaus.mSigmaBkg->setRange(0.2,4.);
  fGausGaus.mSigmaBkg->setUnit("a. u.");

  for (auto list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Taking all the necessary histogram to perform the analysis
    TH3F *fATOFsignal = (TH3F*)list->Get("fATOFsignal");
    TH3F *fMTOFsignal = (TH3F*)list->Get("fMTOFsignal");
    TH3F *fATPCcounts = (TH3F*)list->Get("fATPCcounts");
    TH3F *fMTPCcounts = (TH3F*)list->Get("fMTPCcounts");

    /// Taking information about centrality bins
    const int n_centralities = fATOFsignal->GetNbinsX();
    auto cent_labels = *(fATOFsignal->GetXaxis()->GetXbins());
    /// Taking information about \f$p_{\mathrm{T}}\f$ bins
    const int n_pt_bins = fATOFsignal->GetNbinsY();
    auto pt_axis = fATOFsignal->GetYaxis();
    auto pt_labels = *(pt_axis->GetXbins());

    /// Now it comes a bit of a complication. To improve fit quality
    /// and to minimize the manual interventions one should do the fits to
    /// the same pT bin for all the species and all the centrality classes
    /// before moving to another pT bin. This requires to build all the
    /// arrays before starting the actual analysis.

    /// Build arrays to analyse both deuteron and anti-deuterons
    /// with the same code
    TH3F* tof_histo[2] = {fMTOFsignal,fATOFsignal};
    TH3F* tpc_histo[2] = {fMTPCcounts,fATPCcounts};

    /// Build arrays to analyse all the centrality classes for
    /// both deuteron and anti-deuterons. Complicate stuff just for fun.
    TH1D* hRawCounts[2][kCentLength];
    TH1D* hSignalTailTailGaus[2][kCentLength];
    TH1D* hSystFit[2][kCentLength];
    TH1D* hSignificance[2][kCentLength];
    TH1D* hChiSquare[2][kCentLength];
    TH1D* hChiSquareTPC[2][kCentLength];
    TH1D* hTPConly[2][kCentLength];

    vector<float> n_sigma_vec = {3.0,3.1,3.2,3.3,3.4,3.5};
    vector<float> v_shift = {-0.1,0.05,0.,0.05,0.1};
    int n_shifts = v_shift.size();
    int kNewGreen = kGreen + 3;
    int color_vector[] = {kBlack,kBlue,kNewGreen,kOrange,kRed};
    TH1D* hWidenRangeSyst[2][kCentLength];
    TH1D* hShiftRangeSyst[2][kCentLength];
    TH1D* hWidenRangeSystTPC[2][kCentLength];
    TH1D* hShiftRangeSystTPC[2][kCentLength];

    /// Creating the directories to be used to store the results
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      //dir->mkdir("Fits");
      TDirectory* sig_dir = dir->mkdir("TailTail");
      sig_dir->cd();
      for (int iC = 0; iC < kCentLength; ++iC) sig_dir->mkdir(Form("C_%d",iC));
      dir->cd();
      TDirectory* side_dir = dir->mkdir("Sidebands");
      side_dir->cd();
      for (int iC = 0; iC < kCentLength; ++iC) side_dir->mkdir(Form("C_%d",iC));
      dir->cd();
      dir->mkdir("Significance");
      dir->mkdir("Systematic");
      dir->mkdir("TPConly");
      dir->mkdir("ChiSquare");
    }

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < kCentLength; ++iC) {
        hTPConly[iS][iC] = new TH1D(Form("hTPConly%c%i",kLetter[iS],iC),";p_{T} GeV/c; TPC raw counts",n_pt_bins,pt_labels.GetArray());
        hSignificance[iS][iC] = new TH1D(Form("hSignificance%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}",n_pt_bins,pt_labels.GetArray());
        hChiSquare[iS][iC] = new TH1D(Form("hChiSquare%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #chi^{2}/NDF",n_pt_bins,pt_labels.GetArray());
        hChiSquareTPC[iS][iC] = new TH1D(Form("hChiSquareTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); #chi^{2}/NDF",n_pt_bins,pt_labels.GetArray());
        hRawCounts[iS][iC] = new TH1D(Form("hRawCounts%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RawCounts",n_pt_bins,pt_labels.GetArray());
        hSignalTailTailGaus[iS][iC] = new TH1D(Form("hSignalTailTailGaus%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RawCounts",n_pt_bins,pt_labels.GetArray());
        hWidenRangeSyst[iS][iC] = new TH1D(Form("hWidenRangeSyst%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",n_pt_bins,pt_labels.GetArray());
        hShiftRangeSyst[iS][iC] = new TH1D(Form("hShiftRangeSyst%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",n_pt_bins,pt_labels.GetArray());
        hWidenRangeSystTPC[iS][iC] = new TH1D(Form("hWidenRangeSystTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",n_pt_bins,pt_labels.GetArray());
        hShiftRangeSystTPC[iS][iC] = new TH1D(Form("hShiftRangeSystTPC%c%i",kLetter[iS],iC),"; p_{T}(GeV/c); RMS",n_pt_bins,pt_labels.GetArray());
      }
    }

    float width_range_syst = 0.;
    float width_range_syst_tpc = 0.;
    float pos_range_syst = 0.;
    float pos_range_syst_tpc = 0.;

    for (int iB = 0; iB < n_pt_bins; ++iB) {
      if (pt_axis->GetBinCenter(iB+1) < kPtRange[0] || pt_axis->GetBinCenter(iB+1) > kPtRange[1]) continue;
      float sigma_deut[kCentLength];
      float sigma_deut_tpc[kCentLength];
      for (int iS = 0; iS < 2; ++iS) {
        for (int iC = 0; iC < kCentLength; ++iC) {
          // TOF analysis
          if (pt_axis->GetBinCenter(iB+1) > kCentPtLimits[iC]) continue;
          // std::cout <<"test" <<std::endl;
          TString iTitle = Form(" %1.0f - %1.0f %% %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", cent_labels[kCentBinsArray[iC][0]], cent_labels[kCentBinsArray[iC][1]], pt_labels[iB], pt_labels[iB + 1]);
          TString iName = Form("d%i_%i",iC,iB);
          TH1D *dat = tof_histo[iS]->ProjectionZ(Form("data%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iB + 1,iB + 1);
          // std::cout <<"test2" <<std::endl;
          // if (pt_axis->GetBinCenter(iB+1) > kPtRebin[iC]){
          //   dat->Rebin();
          //   fExpExpTailGaus.mTau0->setRange(-10.,-3.);
          //   fExpExpTailGaus.mTau1->setRange(-2.,-0.01);
          //   fExpExpTailGaus.mKbkg->setRange(0.2,1.);
          // }
          // else{
          //   fExpExpTailGaus.mTau0->setRange(-10.,-1.5);
          //   fExpExpTailGaus.mTau1->setRange(-0.5,-0.01);
          //   fExpExpTailGaus.mKbkg->setRange(0.,1.);
          // }
          fExpExpTailGaus.mTau0->setVal(-5.);
          fExpExpTailGaus.mTau1->setVal(-0.3);
          fExpExpTailGaus.mTau0->setVal(0.5);
          RooDataHist data("data","data",RooArgList(m),Import(*dat));

          /// TailTail
          base_dir->cd(Form("%s/TailTail/C_%d",kNames[iS].data(),iC));
          if(iB<=8){
            fExpExpTailGaus.UseBackground(true);
            fExpExpTailGaus.mMu->setRange(-0.5,1.5);
            fExpExpTailGaus.mSigma->setRange(0.05,0.3);
            fExpExpTailGaus.mSigma->setVal(0.1);

            fExpPolDSCrystalBall.UseBackground(true);
            fExpPolDSCrystalBall.mMu->setRange(-0.5,1.5);
            fExpPolDSCrystalBall.mSigma->setRange(0.05,0.3);
            fExpPolDSCrystalBall.mSigma->setVal(0.1);
          }
          else{
            fExpExpTailGaus.UseBackground(true);
            fExpExpTailGaus.mSigma->setRange(0.05,0.3);
            fExpExpTailGaus.mSigma->setVal(0.1);
            if(iB>=9 && iB<=10){
              fExpExpTailGaus.mKbkg->setVal(0.);
              fExpExpTailGaus.mKbkg->setConstant(true);
              fExpExpTailGaus.mTau0->setConstant(true);
              fBkg.mKbkg->setVal(0.);
              fBkg.mKbkg->setConstant(true);
              fBkg.mTau0->setConstant(true);
            }
            else{
              fExpExpTailGaus.mKbkg->setVal(0.5);
              fExpExpTailGaus.mKbkg->setConstant(false);
              fExpExpTailGaus.mTau0->setConstant(false);
              fBkg.mKbkg->setVal(0.5);
              fBkg.mKbkg->setConstant(false);
              fBkg.mTau0->setConstant(false);
            }
          }
          // if(iS==1){
          //   fExpExpTailGaus.mSigma->setVal(sigma_deut[iC]);
          //   fExpExpTailGaus.mSigma->setConstant(true);
          // }
          if(iB<=8){
              RooPlot* expPolDSCrystalBall = fExpPolDSCrystalBall.FitData(dat, iName, iTitle, "Full", "Full",false,kFitminPt,kFitmaxPt);
            fExpPolDSCrystalBall.mSigma->setConstant(false);
            if(iS==0) sigma_deut[iC] = fExpPolDSCrystalBall.mSigma->getVal();
            if(pt_axis->GetBinCenter(iB+1) > kTOFminPt) expPolDSCrystalBall->Write();
            hSignalTailTailGaus[iS][iC]->SetBinContent(iB+1,fExpPolDSCrystalBall.mSigCounts->getVal());
            hSignalTailTailGaus[iS][iC]->SetBinError(iB+1,fExpPolDSCrystalBall.mSigCounts->getError());
            hRawCounts[iS][iC]->SetBinContent(iB + 1, fExpPolDSCrystalBall.mSigCounts->getVal());
            hRawCounts[iS][iC]->SetBinError(iB + 1, fExpPolDSCrystalBall.mSigCounts->getError());
          }
          else{
            RooPlot* expExpTailTailGausPlot = fExpExpTailGaus.FitData(dat, iName, iTitle, "Full", "Full",false,kFitminPt,kFitmaxPt);
            fExpExpTailGaus.mSigma->setConstant(false);
            if(iS==0) sigma_deut[iC] = fExpExpTailGaus.mSigma->getVal();
            if(pt_axis->GetBinCenter(iB+1) > kTOFminPt) expExpTailTailGausPlot->Write();
            hSignalTailTailGaus[iS][iC]->SetBinContent(iB+1,fExpExpTailGaus.mSigCounts->getVal());
            hSignalTailTailGaus[iS][iC]->SetBinError(iB+1,fExpExpTailGaus.mSigCounts->getError());
            hRawCounts[iS][iC]->SetBinContent(iB + 1, fExpExpTailGaus.mSigCounts->getVal());
            hRawCounts[iS][iC]->SetBinError(iB + 1, fExpExpTailGaus.mSigCounts->getError());
          }
          /// Bin counting TOF
          float residual_vector[n_sigma_vec.size()];
          for(size_t iSigma=0; iSigma < n_sigma_vec.size(); iSigma++){
            float left_sigma = fExpExpTailGaus.mMu->getVal() - n_sigma_vec[iSigma]*fExpExpTailGaus.mSigma->getVal();
            float right_sigma = fExpExpTailGaus.mMu->getVal() + (float(n_sigma_vec[iSigma])+2.)*fExpExpTailGaus.mSigma->getVal();
            int left_edge_bin = dat->FindBin(left_sigma);
            float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
            int right_edge_bin = dat->FindBin(right_sigma);
            float right_edge_float = dat->GetBinLowEdge(right_edge_bin+1);
            fBkg.mX->setRange("signal",left_edge_float,right_edge_float);
            if (iSigma==0) {
              fBkg.mX->setRange("left",kFitminPt,left_edge_float);
              fBkg.mX->setRange("right",right_edge_float,kFitmaxPt);
              RooPlot* bkgPlot = fBkg.FitData(dat, Form("%s_sideband",iName.Data()), iTitle, "left,right","Full");
              base_dir->cd(Form("%s/Sidebands/C_%d",kNames[iS].data(),iC));
              bkgPlot->Write();

            }
            float bkg_integral = (iB>8) ? fBkg.mBackground->createIntegral(m_bis,NormSet(m_bis),Range("signal"))->getVal() * fBkg.mBkgCounts->getVal() : 0;
            if(iB>8){
              hChiSquare[iS][iC]->SetBinContent(iB+1, fBkg.mChi2);
              hChiSquare[iS][iC]->SetBinError(iB+1, 0);
            }
            float tot_integral = dat->Integral(left_edge_bin,right_edge_bin);
            float sig_integral = tot_integral - bkg_integral;
            float sig_err = TMath::Sqrt(tot_integral+bkg_integral);
            if(iSigma==0){
              // hRawCounts[iS][iC]->SetBinContent(iB + 1, sig_integral);
              // hRawCounts[iS][iC]->SetBinError(iB + 1, sig_err);
              hSignificance[iS][iC]->SetBinContent(iB + 1, sig_integral / TMath::Sqrt(tot_integral));
            }
            residual_vector[iSigma] = sig_integral;
          }
          width_range_syst = TMath::RMS(n_sigma_vec.size(),residual_vector);
          width_range_syst /= hRawCounts[iS][iC]->GetBinContent(iB + 1);
          hWidenRangeSyst[iS][iC]->SetBinContent(iB + 1, width_range_syst);
          // Moving the counting range
          float shift_vector[n_shifts];
          for(int iShift=0; iShift<n_shifts; iShift++){
            float left_sigma = fExpExpTailGaus.mMu->getVal()-3.*fExpExpTailGaus.mSigma->getVal()-v_shift[iShift];
            float right_sigma = fExpExpTailGaus.mMu->getVal()+5.*fExpExpTailGaus.mSigma->getVal()-v_shift[iShift];
            int left_edge_bin = dat->FindBin(left_sigma);
            float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
            int right_edge_bin = dat->FindBin(right_sigma);
            float right_edge_float = dat->GetBinLowEdge(right_edge_bin+1);
            fBkg.mX->setRange("signal",left_edge_float,right_edge_float);
            float bkg_integral = (iB>7) ? fBkg.mBackground->createIntegral(m_bis,NormSet(m_bis),Range("signal"))->getVal() * fBkg.mBkgCounts->getVal() : 0;
            float tot_integral = dat->Integral(left_edge_bin,right_edge_bin);
            float sig_integral = tot_integral - bkg_integral;
            float sig_err = TMath::Sqrt(tot_integral+bkg_integral);
            shift_vector[iShift] = sig_integral;
          }
          pos_range_syst = TMath::RMS(n_shifts,shift_vector);
          pos_range_syst /= hRawCounts[iS][iC]->GetBinContent(iB + 1);
          hShiftRangeSyst[iS][iC]->SetBinContent(iB + 1, pos_range_syst);

          /// TPC analysis
          if(pt_axis->GetBinCenter(iB+1) < kTPCmaxPt){
            base_dir->cd(Form("%s/TPConly",kNames[iS].data()));
            TH1D *tpc_dat = tpc_histo[iS]->ProjectionZ(Form("tpc_data%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1], iB + 1, iB + 1);
            RooDataHist tpc_data("tpc_data","tpc_data",RooArgList(ns),Import(*tpc_dat));

            if(iB<=5) fGausGaus.UseBackground(false);
            else fGausGaus.UseBackground(true);
            // if(iS==1){
            //   fGausGaus.mSigma->setVal(sigma_deut_tpc[iC]);
            //   fGausGaus.mSigma->setConstant(true);
            // }
            RooPlot* gausGausPlot = fGausGaus.FitData(tpc_dat,Form("TPC_d_%i_%i",iC,iB),iTitle,"Full","Full");
            fGausGaus.mSigma->setConstant(false);
            if(iS==0) sigma_deut_tpc[iC] = fGausGaus.mSigma->getVal();
            gausGausPlot->Write();
            float bin_width = tpc_dat->GetBinWidth(1);

            /// Bin Counting TPC
            float count_tpc_vector[n_sigma_vec.size()];
            float count_tpc_err_vector[n_sigma_vec.size()];
            float TpcChi=0.;

            for(size_t iSigma=0; iSigma < n_sigma_vec.size(); iSigma++){
              float left_sigma_tpc = fGausGaus.mMu->getVal()-n_sigma_vec[iSigma]*fGausGaus.mSigma->getVal();
              float right_sigma_tpc = fGausGaus.mMu->getVal()+n_sigma_vec[iSigma]*fGausGaus.mSigma->getVal();
              int left_edge_bin_tpc = tpc_dat->FindBin(left_sigma_tpc);
              float left_edge_float_tpc = tpc_dat->GetBinLowEdge(left_edge_bin_tpc);
              int right_edge_bin_tpc = tpc_dat->FindBin(right_sigma_tpc);
              float right_edge_float_tpc = tpc_dat->GetBinLowEdge(right_edge_bin_tpc+1);
              fGausGaus.mX->setRange("signal",left_edge_float_tpc,right_edge_float_tpc);
              float integral_tpc = tpc_dat->Integral(left_edge_bin_tpc,right_edge_bin_tpc);
              float err2 = integral_tpc;
              if(iB>5){
                integral_tpc -= fGausGaus.mBackground->createIntegral(ns,NormSet(ns),Range("signal"))->getVal() * fGausGaus.mBkgCounts->getVal();
                err2 += fGausGaus.mBackground->createIntegral(ns,NormSet(ns),Range("signal"))->getVal() * fGausGaus.mBkgCounts->getVal();
                if(iSigma==0) TpcChi= fGausGaus.mChi2;
              }
              count_tpc_vector[iSigma] = integral_tpc;
              count_tpc_err_vector[iSigma] = TMath::Sqrt(err2);
            }
            if(iB>5 && iB<=8){
              hChiSquareTPC[iS][iC]->SetBinContent(iB+1,TpcChi);
              hChiSquareTPC[iS][iC]->SetBinError(iB+1,0.);
            }
            float sValTPC = count_tpc_vector[0];
            float sErrTPC = count_tpc_err_vector[0];
            hTPConly[iS][iC]->SetBinContent(iB + 1, fGausGaus.mSigCounts->getVal());
            hTPConly[iS][iC]->SetBinError(iB + 1, fGausGaus.mSigCounts->getError());
            // hTPConly[iS][iC]->SetBinContent(iB + 1, sValTPC);
            // hTPConly[iS][iC]->SetBinError(iB + 1, sErrTPC);
            width_range_syst_tpc = TMath::RMS(n_sigma_vec.size(),count_tpc_vector);
            width_range_syst_tpc /= hTPConly[iS][iC]->GetBinContent(iB + 1);
            hWidenRangeSystTPC[iS][iC]->SetBinContent(iB + 1, width_range_syst_tpc);

            float shift_vector_tpc[n_shifts];
            for(int iShift=0; iShift<n_shifts; iShift++){
              float left_sigma = fGausGaus.mMu->getVal()-3.*fGausGaus.mSigma->getVal()-v_shift[iShift];
              float right_sigma = fGausGaus.mMu->getVal()+3.*fGausGaus.mSigma->getVal()-v_shift[iShift];
              int left_edge_bin = dat->FindBin(left_sigma);
              float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
              int right_edge_bin = dat->FindBin(right_sigma);
              float right_edge_float = dat->GetBinLowEdge(right_edge_bin+1);
              fGausGaus.mX->setRange("signal",left_edge_float,right_edge_float);
              float bkg_integral = (iB>5) ? fGausGaus.mBackground->createIntegral(ns,NormSet(ns),Range("signal"))->getVal() * fGausGaus.mBkgCounts->getVal() : 0;
              float tot_integral = tpc_dat->Integral(left_edge_bin,right_edge_bin);
              float sig_integral = tot_integral - bkg_integral;
              float sig_err = TMath::Sqrt(tot_integral+bkg_integral);
              shift_vector[iShift] = sig_integral;
            }
            pos_range_syst_tpc = TMath::RMS(n_shifts,shift_vector_tpc);
            pos_range_syst_tpc /= hTPConly[iS][iC]->GetBinContent(iB + 1);
            hShiftRangeSystTPC[iS][iC]->SetBinContent(iB + 1, pos_range_syst_tpc);
          }
        }
      }
    }

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < kCentLength; ++iC) {
        base_dir->cd(Form("%s/TailTail",kNames[iS].data()));
        hRawCounts[iS][iC]->Write();
        hSignalTailTailGaus[iS][iC]->Write();
        base_dir->cd(Form("%s/Systematic",kNames[iS].data()));
        hShiftRangeSyst[iS][iC]->Write();
        hWidenRangeSyst[iS][iC]->Write();
        hWidenRangeSystTPC[iS][iC]->Write();
        hShiftRangeSystTPC[iS][iC]->Write();
        base_dir->cd(Form("%s/Significance",kNames[iS].data()));
        hSignificance[iS][iC]->Write();
        base_dir->cd(Form("%s/TPConly",kNames[iS].data()));
        hTPConly[iS][iC]->Write();
        base_dir->cd(Form("%s/ChiSquare",kNames[iS].data()));
        hChiSquare[iS][iC]->Write();
        hChiSquareTPC[iS][iC]->Write();
      }
    }
    base_dir->Close();
    break;
  }
  output_file.Close();
}
