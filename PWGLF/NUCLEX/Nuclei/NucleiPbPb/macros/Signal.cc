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
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  /// Taking all the histograms from the MC file
  TFile input_file(kDataFilename.data());
  TFile output_file(kSignalOutput.data(),"recreate");

  /// Setting up the fitting environment for TOF analysis
  RooRealVar m("dm2","m_{TOF}^{2} - m^{2}_{PDG}",-2.,2.5,"GeV^{2}/c^{4}");
  m.setBins(10000,"cache");
  m.setRange("Full", -2., 2.5);

  FitExpExpTailTailGaus fExpExpTailTailGaus(&m);
  fExpExpTailTailGaus.mMu->setRange(-.1,1);
  fExpExpTailTailGaus.mSigma->setRange(0.05,0.3);
  fExpExpTailTailGaus.mSigma->setVal(0.1);
  fExpExpTailTailGaus.mAlpha0->setRange(-3.,-1.);
  fExpExpTailTailGaus.mAlpha1->setRange(0.5,3.);
  fExpExpTailTailGaus.mSigCounts->setRange(0.,4000.);

  //Background
  RooRealVar m_bis("dm2_bis","m_{TOF}^{2} - m^{2}_{PDG}",-1.2,1.5,"GeV^{2}/c^{4}");
  m_bis.setBins(10000,"cache");
  m_bis.setRange("Full", -1.2, 1.5);
  FitExpExpTailGaus fBkg(&m_bis);
  fBkg.UseSignal(false);

  // Setting up the fitting environment for the TPC analysis
  RooRealVar ns("ns","n#sigma_{d}",-3.,3,"a. u.");
  ns.setBins(1000,"cache");
  ns.setRange("Full", -3., 3.);
  ns.setRange("Special", -3., 3.);

  // TPC analysis
  FitGausGaus fGausGaus(&ns);
  fGausGaus.mSigma->setRange(0.4,0.9);
  fGausGaus.mSigma->setVal(0.75);
  fGausGaus.mMu->setRange(-0.5,0.5);
  fGausGaus.mSigma->setRange(0.02,0.9);
  fGausGaus.mMuBkg->setRange(-10.,-2.);
  fGausGaus.mMuBkg->setVal(-7);
  fGausGaus.mSigmaBkg->setRange(0.2,4.);

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
    TH1D* hRawCounts[2][n_centralities];
    TH1D* hSignalTailTailGaus[2][n_centralities];
    TH1D* hSystFit[2][n_centralities];
    TH1D* hSignificance[2][n_centralities];
    TH1D* hChiSquare[2][n_centralities];
    TH1D* hChiSquareTPC[2][n_centralities];
    TH1D* hTPConly[2][n_centralities];

    vector<float> n_sigma_vec={3.0,3.1,3.2,3.3,3.4,3.5};
    int n_vec_sigma = n_sigma_vec.size();
    vector<float> v_shift = {-0.1,0.05,0.,0.05,0.1};
    int n_shifts = v_shift.size();
    int kNewGreen = kGreen + 3;
    int color_vector[] = {kBlack,kBlue,kNewGreen,kOrange,kRed};
    TH1D* hWidenRangeSyst[2][n_centralities];
    TH1D* hShiftRangeSyst[2][n_centralities];
    TH1D* hWidenRangeSystTPC[2][n_centralities];
    TH1D* hShiftRangeSystTPC[2][n_centralities];

    /// Creating the directories to be used to store the results
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      //dir->mkdir("Fits");
      dir->mkdir("TailTail");
      dir->mkdir("Sidebands");
      dir->mkdir("Significance");
      dir->mkdir("Systematic");
      dir->mkdir("TPConly");
      dir->mkdir("ChiSquare");
    }

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < n_centralities; ++iC) {
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
      if ( pt_axis->GetBinCenter(iB+1) < kPtRange[0] || pt_axis->GetBinCenter(iB+1) > kPtRange[1]) continue;
      float sigma_deut[n_centralities];
      float sigma_deut_tpc[n_centralities];
      for (int iS = 0; iS < 2; ++iS) {
        for (int iC = 0; iC < n_centralities; ++iC) {
          // TOF analysis
          TString iTitle = Form("%1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}",pt_labels[iB],pt_labels[iB + 1]);
          TString iName = Form("d%i_%i",iC,iB);
          TH1D *dat = tof_histo[iS]->ProjectionZ(Form("data%i_%i",iC,iB),iC + 1,iC + 1, iB + 1, iB + 1);
          RooDataHist data("data","data",RooArgList(m),Import(*dat));

          /// TailTail
          base_dir->cd(Form("%s/TailTail",kNames[iS].data()));
          if(iB<=7){
            fExpExpTailTailGaus.UseBackground(false);
            fExpExpTailTailGaus.mSigma->setRange(0.05,0.13);
          }
          else  fExpExpTailTailGaus.UseBackground(true);
          if(iS==1){
            fExpExpTailTailGaus.mSigma->setVal(sigma_deut[iC]);
            fExpExpTailTailGaus.mSigma->setConstant(true);
          }
          RooPlot* expExpTailTailGausPlot = fExpExpTailTailGaus.FitData(dat, iName, iTitle, "Full", "Full",true,-1.2,1.5);
          fExpExpTailTailGaus.mSigma->setConstant(false);
          if(iS==0) sigma_deut[iC] = fExpExpTailTailGaus.mSigma->getVal();
          if(pt_axis->GetBinCenter(iB+1) > kTOFminPt) expExpTailTailGausPlot->Write();
          hSignalTailTailGaus[iS][iC]->SetBinContent(iB+1,fExpExpTailTailGaus.mSigCounts->getVal());
          hSignalTailTailGaus[iS][iC]->SetBinError(iB+1,fExpExpTailTailGaus.mSigCounts->getError());

          /// Bin counting TOF
          float residual_vector[n_vec_sigma];
          for(size_t iSigma=0; iSigma < n_sigma_vec.size(); iSigma++){
            float left_sigma = fExpExpTailTailGaus.mMu->getVal()-n_sigma_vec[iSigma]*fExpExpTailTailGaus.mSigma->getVal();
            float right_sigma = fExpExpTailTailGaus.mMu->getVal()+(n_sigma_vec[iSigma])*fExpExpTailTailGaus.mSigma->getVal();
            int left_edge_bin = dat->FindBin(left_sigma);
            float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
            int right_edge_bin = dat->FindBin(right_sigma);
            float right_edge_float = dat->GetBinLowEdge(right_edge_bin+1);
            fBkg.mX->setRange("signal",left_edge_float,right_edge_float);
            if (!iSigma) {
              fBkg.mX->setRange("left",-1.2,left_edge_float);
              fBkg.mX->setRange("right",right_edge_float,1.5);
              RooPlot* bkgPlot = fBkg.FitData(dat, Form("%s_%d_sideband",iName.Data(),(int)iSigma), iTitle, "left,right","Full",true);
              base_dir->cd(Form("%s/Sidebands",kNames[iS].data()));
              bkgPlot->Write();

            }
            float bkg_integral = (iB>7) ? fBkg.mBackground->createIntegral(m_bis,NormSet(m_bis),Range("signal"))->getVal() * fBkg.mBkgCounts->getVal() : 0;
            if(iB>7){
              hChiSquare[iS][iC]->SetBinContent(iB+1, fBkg.mChi2);
              hChiSquare[iS][iC]->SetBinError(iB+1, 0);
            }
            float tot_integral = dat->Integral(left_edge_bin,right_edge_bin);
            float sig_integral = tot_integral - bkg_integral;
            float sig_err = TMath::Sqrt(tot_integral+bkg_integral);
            if(iSigma==0){
              hRawCounts[iS][iC]->SetBinContent(iB + 1, sig_integral);
              hRawCounts[iS][iC]->SetBinError(iB + 1, sig_err);
              hSignificance[iS][iC]->SetBinContent(iB + 1, sig_integral / TMath::Sqrt(tot_integral));
            }
            residual_vector[iSigma] = sig_integral;
          }
          width_range_syst = TMath::RMS(n_vec_sigma,residual_vector);
          width_range_syst /= hRawCounts[iS][iC]->GetBinContent(iB + 1);
          hWidenRangeSyst[iS][iC]->SetBinContent(iB + 1, width_range_syst);
          // Moving the counting range
          float shift_vector[n_shifts];
          for(int iShift=0; iShift<n_shifts; iShift++){
            float left_sigma = fExpExpTailTailGaus.mMu->getVal()-3.*fExpExpTailTailGaus.mSigma->getVal()-v_shift[iShift];
            float right_sigma = fExpExpTailTailGaus.mMu->getVal()+3.*fExpExpTailTailGaus.mSigma->getVal()-v_shift[iShift];
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
            TH1D *tpc_dat = tpc_histo[iS]->ProjectionZ(Form("tpc_data%i_%i",iC,iB),iC + 1,iC + 1, iB + 1, iB + 1);
            RooDataHist tpc_data("tpc_data","tpc_data",RooArgList(ns),Import(*tpc_dat));

            if(iB<=5) fGausGaus.UseBackground(false);
            else fGausGaus.UseBackground(true);
            if(iS==1){
              fGausGaus.mSigma->setVal(sigma_deut_tpc[iC]);
              fGausGaus.mSigma->setConstant(true);
            }
            RooPlot* gausGausPlot = fGausGaus.FitData(tpc_dat,Form("TPC_d_%i_%i",iC,iB),iTitle,"Full","Full");
            fGausGaus.mSigma->setConstant(false);
            if(iS==0) sigma_deut_tpc[iC] = fGausGaus.mSigma->getVal();
            gausGausPlot->Write();
            float bin_width = tpc_dat->GetBinWidth(1);

            /// Bin Counting TPC
            float count_tpc_vector[n_vec_sigma];
            float count_tpc_err_vector[n_vec_sigma];
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
            if(iB>5 && iB<=7){
              hChiSquareTPC[iS][iC]->SetBinContent(iB+1,TpcChi);
              hChiSquareTPC[iS][iC]->SetBinError(iB+1,0.);
            }
            float sValTPC = count_tpc_vector[0];
            float sErrTPC = count_tpc_err_vector[0];
            hTPConly[iS][iC]->SetBinContent(iB + 1, sValTPC);
            hTPConly[iS][iC]->SetBinError(iB + 1, sErrTPC);
            width_range_syst_tpc = TMath::RMS(n_vec_sigma,count_tpc_vector);
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
      for (int iC = 0; iC < n_centralities; ++iC) {
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
  }
  output_file.Close();
}
