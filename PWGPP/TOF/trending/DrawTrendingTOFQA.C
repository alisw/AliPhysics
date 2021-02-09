/*
  fbellini@cern.ch - last update on 27/02/2014
  Macro to draw the TOF QA trending plots by accessing the std tree.
  To be mainly used with the automatic scripts to fill the QA repository.
  Launch with 
  aliroot -x -b -q "DrawTrendingTOFQA.C" 
  The macro produces one png file for each trending variables
  and a .root file with the histograms
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGrid.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#endif

//Utility variables
TString plotDir = "./";
TString plotExt = "png";

//Utility functions
void DrawStdCanvas(TH1* h1, TH1* h2, TH1* h3, TH1* h4, TH1* h5, Double_t yMin, Double_t yMax, TLegend* leg = nullptr, const TString cname = "", Bool_t Gridy = kFALSE);
void DrawStdCanvas(TH1* h1, TH1* h2, TH1* h3, Double_t yMin, Double_t yMax, TLegend* leg = nullptr, TString cname = "");
void DrawStdCanvas(TH1* h1, TH1* h2, Double_t yMin, Double_t yMax, TLegend* leg = nullptr, TString cname = "");
void DrawStdCanvas(TH1* h1, Double_t yMin, Double_t yMax, TString cname = "");
void SetPlotExt();
TLegend* GetLeg(Double_t x1, Double_t y1, Double_t x2, Double_t y2);

//Main function
Int_t DrawTrendingTOFQA(TString mergedTrendFile = "trending.root", // trending tree file name
    Bool_t displayAll = kFALSE)                                    //set to kTRUE to display trending for expert plots
{
  //
  //reads merged trending.root file and draws trending plots from tree
  //
  if (mergedTrendFile.IsNull() || mergedTrendFile.IsWhitespace() || !mergedTrendFile) {
    Printf("TOF QA input file name not provided");
    return 1;
  }
  const TString outfilename = "ProductionQA.hist.root"; //Output file name

  //Tree variables
  Int_t runNumber = 0;
  Double_t avTime = 0., peakTime = 0., spreadTime = 0., peakTimeErr = 0., spreadTimeErr = 0., negTimeRatio = 0.;
  Double_t avRawTime = 0., peakRawTime = 0., spreadRawTime = 0., peakRawTimeErr = 0., spreadRawTimeErr = 0.;
  Double_t avTot = 0., peakTot = 0., spreadTot = 0., peakTotErr = 0., spreadTotErr = 0.;
  Double_t meanResTOF = 0., spreadResTOF = 0., meanResTOFerr = 0., spreadResTOFerr = 0.;
  Double_t orphansRatio = 0., avL = 0., negLratio = 0.;
  Double_t matchEffIntegratedErr = -9999., matchEffIntegrated = -9999., matchEffLinFit1Gev = 0., matchEffLinFit1GevErr = 0.;
  Double_t avPiDiffTime = 0., peakPiDiffTime = 0., spreadPiDiffTime = 0., peakPiDiffTimeErr = 0., spreadPiDiffTimeErr = 0., avT0fillRes = 0.;

  Double_t avT0A = 0., peakT0A = 0., spreadT0A = 0., peakT0AErr = 0., spreadT0AErr = 0.;
  Double_t avT0C = 0., peakT0C = 0., spreadT0C = 0., peakT0CErr = 0., spreadT0CErr = 0.;
  Double_t avT0AC = 0., peakT0AC = 0., spreadT0AC = 0., peakT0ACErr = 0., spreadT0ACErr = 0.;
  Double_t avT0res = 0., peakT0res = 0., spreadT0res = 0., peakT0resErr = 0., spreadT0resErr = 0.;

  Double_t StartTime_pBestT0 = 0.0, StartTime_pBestT0Err = 0.0, StartTime_pFillT0 = 0.0, StartTime_pFillT0Err = 0.0;
  Double_t StartTime_pTOFT0 = 0.0, StartTime_pTOFT0Err = 0.0, StartTime_pT0ACT0 = 0.0, StartTime_pT0ACT0Err = 0.0, StartTime_pT0AT0 = 0.0, StartTime_pT0AT0Err = 0.0, StartTime_pT0CT0 = 0.0, StartTime_pT0CT0Err = 0.0;
  Double_t StartTime_pBestT0_Res = 0.0, StartTime_pFillT0_Res = 0.0, StartTime_pTOFT0_Res = 0.0, StartTime_pT0ACT0_Res = 0.0, StartTime_pT0AT0_Res = 0.0, StartTime_pT0CT0_Res = 0.0;

  Float_t avMulti = 0;
  Float_t fractionEventsWHits = 0.0;
  Double_t goodChannelsRatio = 0.0;
  Double_t goodChannelsRatioInAcc = 0.0;

  TFile fin(mergedTrendFile, "READ");
  if (!fin.IsOpen()) {
    ::Error("DrawTrendingTOFQA", "input tree file %s not found!!", mergedTrendFile.Data());
    return 3;
  }
  TTree* ttree = nullptr;
  fin.GetObject("trending", ttree);
  if (!ttree) {
    fin.ls();
    ::Error("DrawTrendingTOFQA", " Invalid trending tree.");
    return 2;
  }
  //Set branch addresses
#define SetAddress(var)                                                                             \
  if (0 != ttree->SetBranchAddress(#var, &var)) {                                                   \
    ::Warning("DrawTrendingTOFQA", "Address for variable %s is not found in input TTree !!", #var); \
    var = 0;                                                                                        \
  }
  ttree->SetBranchAddress("run", &runNumber); // Run number
  SetAddress(avMulti);                        // Average TOF hit multiplicity
  SetAddress(goodChannelsRatio);              // Fraction of active TOF channels
  SetAddress(goodChannelsRatioInAcc);         // Fraction of active TOF channels in acceptance
  SetAddress(avTime);                         // mean time
  SetAddress(peakTime);                       // main peak time after fit
  SetAddress(spreadTime);                     // spread of main peak of time after fit
  SetAddress(peakTimeErr);                    // main peak time after fit error
  SetAddress(spreadTimeErr);                  // spread of main peak of time after fit error
  SetAddress(negTimeRatio);                   // negative time ratio
  SetAddress(avRawTime);                      // mean raw time
  SetAddress(peakRawTime);                    // mean peak of raw time after fit
  SetAddress(spreadRawTime);                  // spread of main peak of raw time after fit
  SetAddress(peakRawTimeErr);                 // main peak raw  time after fit error
  SetAddress(spreadRawTimeErr);               // spread of  raw main peak of time after fit error
  SetAddress(avTot);                          // main peak tot
  SetAddress(peakTot);                        // main peak of tot after fit
  SetAddress(spreadTot);                      // spread of main peak of tot after fit
  SetAddress(peakTotErr);                     // main peak of tot after fit
  SetAddress(spreadTotErr);                   // spread of main peak of tot after fit
  SetAddress(orphansRatio);                   // orphans ratio
  SetAddress(avL);                            // mean track length
  SetAddress(negLratio);                      // ratio of tracks with track length <350 cm
  SetAddress(matchEffIntegrated);             // matching eff. integrated in pt (1-10GeV/c)
  SetAddress(matchEffIntegratedErr);          // matching eff. error integrated in pt (1-10GeV/c)
  SetAddress(matchEffLinFit1Gev);             // matching eff fit param
  SetAddress(matchEffLinFit1GevErr);          // matching eff fit param error
  SetAddress(avPiDiffTime);                   // mean t-texp_pi
  SetAddress(peakPiDiffTime);                 // main peak t-texp_pi after fit
  SetAddress(spreadPiDiffTime);               // spread of main peak t-texp_pi after fit
  SetAddress(peakPiDiffTimeErr);              // main peak t-texp_pi after fit error
  SetAddress(spreadPiDiffTimeErr);            // spread of main peak of t-texp_pi after fit error
  SetAddress(meanResTOF);                     // mean of t-texp_pi-t0_TOF
  SetAddress(spreadResTOF);                   // spread of t-texp_pi-t0_TOF, ie. resolution
  SetAddress(meanResTOFerr);                  // error on mean of t-texp_pi-t0_TOF
  SetAddress(spreadResTOFerr);                // error on the spread of t-texp_pi-t0_TOF
  SetAddress(avT0A);                          // main peak t0A
  SetAddress(peakT0A);                        // main peak of t0A after fit
  SetAddress(spreadT0A);                      // spread of main peak of t0A after fit
  SetAddress(peakT0AErr);                     // main peak of t0A after fit
  SetAddress(spreadT0AErr);                   // spread of main peak of t0A after fit
  SetAddress(avT0C);                          // main peak t0C
  SetAddress(peakT0C);                        // main peak of t0C after fit
  SetAddress(spreadT0C);                      // spread of main peak of t0C after fit
  SetAddress(peakT0CErr);                     // main peak of t0C after fit
  SetAddress(spreadT0CErr);                   // spread of main peak of t0C after fit
  SetAddress(avT0AC);                         // main peak t0AC
  SetAddress(peakT0AC);                       // main peak of t0AC after fit
  SetAddress(spreadT0AC);                     // spread of main peak of t0AC after fit
  SetAddress(peakT0ACErr);                    // main peak of t0AC after fit
  SetAddress(spreadT0ACErr);                  // spread of main peak of t0AC after fit
  SetAddress(avT0res);                        // main peak t0AC
  SetAddress(peakT0res);                      // main peak of t0AC after fit
  SetAddress(spreadT0res);                    // spread of main peak of t0AC after fit
  SetAddress(peakT0resErr);                   // main peak of t0AC after fit
  SetAddress(spreadT0resErr);                 // spread of main peak of t0AC after fit
  SetAddress(avT0fillRes);                    // t0 fill res
  SetAddress(StartTime_pBestT0);              // T0Best
  SetAddress(StartTime_pFillT0);              // T0Fill
  SetAddress(StartTime_pTOFT0);               // T0TOF
  SetAddress(StartTime_pT0ACT0);              // T0AC
  SetAddress(StartTime_pT0AT0);               // T0A
  SetAddress(StartTime_pT0CT0);               // T0C
  SetAddress(StartTime_pBestT0_Res);          // T0Best
  SetAddress(StartTime_pFillT0_Res);          // T0Fill res
  SetAddress(StartTime_pTOFT0_Res);           // T0TOF res
  SetAddress(StartTime_pT0ACT0_Res);          // T0AC res
  SetAddress(StartTime_pT0AT0_Res);           // T0A res
  SetAddress(StartTime_pT0CT0_Res);           // T0C res
  SetAddress(StartTime_pBestT0Err);           // T0Best
  SetAddress(StartTime_pFillT0Err);           // T0Fill
  SetAddress(StartTime_pTOFT0Err);            // T0TOF
  SetAddress(StartTime_pT0ACT0Err);           // T0AC
  SetAddress(StartTime_pT0AT0Err);            // T0A
  SetAddress(StartTime_pT0CT0Err);            // T0C
#undef SetAddress
  //Fetch period-integrated PID plots
  //Pions
  TH2F* hDiffTimePi = (TH2F*)fin.Get("hExpTimePiVsP_all");
  hDiffTimePi->SetTitle("PIONS t-t_{exp,#pi} (from tracking) vs. P");
  //Kaon
  TH2F* hDiffTimeKa = (TH2F*)fin.Get("hExpTimeKaVsP_all");
  hDiffTimeKa->SetTitle("KAONS t-t_{exp,K} (from tracking) vs. P");
  //Protons
  TH2F* hDiffTimePro = (TH2F*)fin.Get("hExpTimeProVsP_all");
  hDiffTimePro->SetTitle("PROTONS t-t_{exp,p} (from tracking) vs. P");

  //Create trending plots
  const Int_t nRuns = ttree->GetEntries();
  TList lista;
  //Parameters for histogram creation
  Int_t lineWidth = 1;
  Int_t lineStyle = 0;
  Int_t markerStyle = 20;
  Color_t markerColor = kBlue;
  //Macro to create histos
#define HCreate(H, t)                          \
  TH1F* H = new TH1F(#H, t, nRuns, 0., nRuns); \
  H->SetMarkerStyle(markerStyle);              \
  H->SetMarkerColor(markerColor);              \
  H->SetLineColor(markerColor);                \
  H->SetLineWidth(lineWidth);                  \
  H->SetLineStyle(lineStyle);                  \
  lista.Add(H);

  HCreate(hAvMulti, "Average multiplicity of matched tracks <N_{TOF}>;; <N_{TOF}>");
  HCreate(hAvDiffTimeVsRun, "Mean t-t_{exp} (no fit);run;<t^{TOF}-t_{exp,#pi}> (ps)");
  HCreate(hPeakDiffTimeVsRun, "t-t_{exp} (gaussian fit) ;; <t^{TOF}-t_{exp,#pi}> (ps)");
  HCreate(hSpreadDiffTimeVsRun, "#sigma(t-t_{exp}) (gaussian fit);; #sigma(t^{TOF}-t_{exp,#pi}) (ns)");
  HCreate(hMeanTOFResVsRun, "Mean value of t-t_{exp,#pi}-t0_{TOF} (ps);;<t^{TOF}-t_{exp,#pi}-t_{0,TOF}> (ps)");
  HCreate(hSigmaTOFResVsRun, "Spread of t-t_{exp,#pi}-t0_{TOF} (ps);;#sigma(t^{TOF}-t_{exp,#pi}-t_{0,TOF}) (ps)");
  HCreate(hAvTimeVsRun, "<t^{TOF}>;;<t^{TOF}> (ns)");
  HCreate(hPeakTimeVsRun, "Peak value of t^{TOF} (landau fit);;t_{peak}^{TOF} (ns)");
  HCreate(hSpreadTimeVsRun, "Spread of t^{TOF} (landau fit);; #sigma(t^{TOF}) (ns)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerStyle = 21;
  markerColor = kGreen;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hAvRawTimeVsRun, "Peak value of raw t^{TOF};;<t_{raw}^{TOF}> (ns)");
  HCreate(hPeakRawTimeVsRun, "Peak value of raw t^{TOF} (landau fit);;t_{peak,raw}^{TOF} (ns)");
  HCreate(hSpreadRawTimeVsRun, "Spread of raw t^{TOF} (landau fit);;#sigma(t_{raw}^{TOF}) (ns)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerStyle = 22;
  markerColor = kBlack;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hAvTotVsRun, "<ToT> (no fit);run;<ToT> (ns)");
  HCreate(hPeakTotVsRun, "<ToT> (gaussian fit);;ToT_{peak} (ns)");
  HCreate(hSpreadTotVsRun, "#sigma(ToT) (gaussian fit);#sigma(ToT) (ns)");

  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hNegTimeRatioVsRun, "Ratio of tracks with t^{TOF}<12.5 ns; ; ratio of tracks with t^{TOF}<12.5 ns (%)");
  HCreate(hOrphansRatioVsRun, "Ratio of orphans (hits with ToT=0); ; ratio of orphans (%)");
  HCreate(hMeanLVsRun, "Average track length;; <L> (cm)");
  HCreate(hNegLRatioVsRun, "Ratio of tracks with L<350 cm;; ratio of tracks with L<350 cm (%)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kRed;
  lineWidth = 2;
  markerStyle = 0;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hMatchEffVsRun, "#epsilon_{match} (linear fit for p_{T}>1.0 GeV/c);;#epsilon_{match} (p_{T}>1.0 GeV/c)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kCyan + 2;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hMatchEffVsRunNormToGoodCh, "#epsilon_{match} normalized to percentage of TOF good channels;;#epsilon_{match}(p_{T}>1.0 GeV/c,|#eta|<0.8)/f_{all good}");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kBlue;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hMatchEffVsRunNormToGoodChInAcc, "#epsilon_{match} normalized to TOF good channels in |#eta|<0.8;;#epsilon_{match}(p_{T}>1.0 GeV/c,|#eta|<0.8/f_{good}(|#eta|<0.8)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kOrange;
  lineStyle = 7;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hMatchEffIntegratedVsRun, "#it{p}_{T} integrated #epsilon_{match}; ; #epsilon_{match} (1 < p_{T} < 10 GeV/c)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kBlack;
  lineStyle = 1;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hPeakT0AVsRun, "Peak value of T0A (gaussian fit);;t0A (ps)");
  HCreate(hPeakT0CVsRun, "Peak value of T0C (gaussian fit);;t0AC (ps)");
  HCreate(hPeakT0ACVsRun, "Peak value of T0AC (gaussian fit);;t0AC (ps)");
  HCreate(hT0fillResVsRun, "t0_fill spread;;t0_spread (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerStyle = 20;
  markerColor = kOrange;
  lineWidth = 2;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0BestVsRun, "start time by best_t0;;t0 Best (ps)");
  HCreate(hT0BestVsRunRes, "#sigma of best_t0;; #sigma t0 Best (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kBlack;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0FillVsRun, "start time by fill_t0;;t0 Fill (ps)");
  HCreate(hT0FillVsRunRes, "fill_t0;; #sigmat0 Fill (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kBlue;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0TOFVsRun, "start time by TOF_t0;;t0 TOF (ps)");
  HCreate(hT0TOFVsRunRes, "TOF_t0;; #sigma t0 TOF (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kRed;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0T0ACVsRun, "start time by T0AC;;t0 T0AC (ps)");
  HCreate(hT0T0ACVsRunRes, "T0AC_t0;; #sigma t0 T0AC (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kGreen + 2;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0T0AVsRun, "start time by T0A;;t0 T0A (ps)");
  HCreate(hT0T0AVsRunRes, "T0A_t0;; #sigma t0 T0A (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kMagenta;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hT0T0CVsRun, "start time by T0C;;t0 T0C (ps)");
  HCreate(hT0T0CVsRunRes, "T0C_t0;; #sigma t0 T0C (ps)");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerStyle = 1;
  markerColor = kCyan - 1;
  lineWidth = 2;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hGoodChannelsRatio, "Fraction of TOF good channels;;fraction of good channels");

  //~~~~~~~~~~~~~~~~~~~~~//
  markerColor = kBlue + 2;
  //~~~~~~~~~~~~~~~~~~~~~//

  HCreate(hGoodChannelsRatioInAcc, "Fraction of TOF good channels in |#eta|<0.8;;fraction of good channels in |#eta|<0.8");

#undef HCreate

  for (Int_t irun = 0; irun < nRuns; irun++) {
    ttree->GetEntry(irun);
    TString runlabel = Form("%i", runNumber);

    // Macro to fill the histograms
#define HFill(H, y, ye)          \
  H->SetBinContent(irun + 1, y); \
  H->SetBinError(irun + 1, ye);  \
  H->GetXaxis()->SetBinLabel(irun + 1, runlabel);

    HFill(hAvMulti, avMulti, 0);
    HFill(hAvDiffTimeVsRun, avPiDiffTime, 0);
    HFill(hPeakDiffTimeVsRun, peakPiDiffTime, peakPiDiffTimeErr);
    HFill(hSpreadDiffTimeVsRun, spreadPiDiffTime, spreadPiDiffTimeErr);
    HFill(hMeanTOFResVsRun, meanResTOF, meanResTOFerr);
    HFill(hSigmaTOFResVsRun, spreadResTOF, spreadResTOFerr);
    HFill(hAvTimeVsRun, avTime, 0);
    HFill(hPeakTimeVsRun, peakTime, peakTimeErr);
    HFill(hSpreadTimeVsRun, spreadTime, spreadTimeErr);
    HFill(hAvRawTimeVsRun, avRawTime, 0);
    HFill(hPeakRawTimeVsRun, peakRawTime, peakRawTimeErr);
    HFill(hSpreadRawTimeVsRun, spreadRawTime, spreadRawTimeErr);
    HFill(hAvTotVsRun, avTot, 0);
    HFill(hPeakTotVsRun, peakTot, peakTotErr);
    HFill(hSpreadTotVsRun, spreadTot, spreadTotErr);
    HFill(hNegTimeRatioVsRun, negTimeRatio, 0);
    HFill(hOrphansRatioVsRun, orphansRatio, 0);
    HFill(hMeanLVsRun, avL, 0);
    HFill(hNegLRatioVsRun, negLratio, 0);
    HFill(hMatchEffVsRun, matchEffLinFit1Gev, matchEffLinFit1GevErr);
    HFill(hMatchEffIntegratedVsRun, matchEffIntegrated, matchEffIntegratedErr);
    HFill(hMatchEffVsRunNormToGoodCh, goodChannelsRatio > 0 ? matchEffLinFit1Gev / goodChannelsRatio : 0.0, matchEffLinFit1GevErr);
    HFill(hMatchEffVsRunNormToGoodChInAcc, goodChannelsRatioInAcc > 0 ? matchEffLinFit1Gev / goodChannelsRatioInAcc : 0.0, matchEffLinFit1GevErr);
    HFill(hGoodChannelsRatio, goodChannelsRatio, 0.0);
    HFill(hGoodChannelsRatioInAcc, goodChannelsRatioInAcc, 0.0);
    HFill(hPeakT0AVsRun, peakT0A, spreadT0A);
    HFill(hPeakT0CVsRun, peakT0C, spreadT0C);
    HFill(hPeakT0ACVsRun, peakT0AC, spreadT0AC);
    HFill(hT0fillResVsRun, avT0fillRes, 0);
    HFill(hT0BestVsRun, StartTime_pBestT0, StartTime_pBestT0Err);
    HFill(hT0FillVsRun, StartTime_pFillT0, StartTime_pFillT0Err);
    HFill(hT0TOFVsRun, StartTime_pTOFT0, StartTime_pTOFT0Err);
    HFill(hT0T0ACVsRun, StartTime_pT0ACT0, StartTime_pT0ACT0Err);
    HFill(hT0T0AVsRun, StartTime_pT0AT0, StartTime_pT0AT0Err);
    HFill(hT0T0CVsRun, StartTime_pT0CT0, StartTime_pT0CT0Err);
    HFill(hT0BestVsRunRes, StartTime_pBestT0_Res, 1.e-5);
    HFill(hT0FillVsRunRes, StartTime_pFillT0_Res, 1.e-5);
    HFill(hT0TOFVsRunRes, StartTime_pTOFT0_Res, 1.e-5);
    HFill(hT0T0ACVsRunRes, StartTime_pT0ACT0_Res, 1.e-5);
    HFill(hT0T0AVsRunRes, StartTime_pT0AT0_Res, 1.e-5);
    HFill(hT0T0CVsRunRes, StartTime_pT0CT0_Res, 1.e-5);
#undef HFill
  }

  TFile* fout = new TFile(outfilename, "RECREATE");
  fout->cd();
  lista.Write();
  fout->Close();

  gStyle->SetOptStat(10);
  //Plot t-texp trend
  DrawStdCanvas(hPeakDiffTimeVsRun, -50, 50);
  DrawStdCanvas(hSpreadDiffTimeVsRun, 0., 400.);

  //Plot average of t-texp-t0tof and resolution trend
  DrawStdCanvas(hMeanTOFResVsRun, -50., 50.);
  DrawStdCanvas(hSigmaTOFResVsRun, 0., 200.);

  //Plot matching efficiency trend
  TLegend* leg = GetLeg(0.5095602, 0.1206897, 0.8891013, 0.3314176);
  //
  leg->AddEntry(hMatchEffVsRun, "#epsilon_{match} (linear fit for p_{T}>1.0 GeV/c)", "lpf");
  leg->AddEntry(hMatchEffIntegratedVsRun, "#epsilon_{match} (integrated for p_{T}>1.0 GeV/c)", "lpf");
  //
  DrawStdCanvas(hMatchEffVsRun, hMatchEffIntegratedVsRun, 0., 1., leg);
  DrawStdCanvas(hMatchEffVsRunNormToGoodChInAcc, 0., 1.);
  DrawStdCanvas(hMatchEffVsRunNormToGoodCh, 0., 1.);

  leg = GetLeg(0.5095602, 0.1206897, 0.8891013, 0.3314176);
  //
  leg->AddEntry(hMatchEffVsRun, "#epsilon_{match} (linear fit for p_{T}>1.0 GeV/c)", "lpf");
  leg->AddEntry(hMatchEffVsRunNormToGoodCh, "#epsilon_{match} norm. to fraction of TOF good channels", "lpf");
  leg->AddEntry(hMatchEffVsRunNormToGoodChInAcc, "#epsilon_{match} norm. to fraction of TOF good channels in |#eta|<0.8", "lpf");
  //
  DrawStdCanvas(hMatchEffVsRun, hMatchEffVsRunNormToGoodCh, hMatchEffVsRunNormToGoodChInAcc, 0.4, 0.8, leg, "cMatchEffSummary");

  //Plot start time trend
  leg = GetLeg(0.6, 0.7, 0.9, 0.9);
  leg->SetNColumns(2);
  //
  leg->AddEntry(hT0TOFVsRun, "", "lp");
  leg->AddEntry(hT0T0ACVsRun, "", "lp");
  leg->AddEntry(hT0T0AVsRun, "", "lp");
  leg->AddEntry(hT0T0CVsRun, "", "lp");
  leg->AddEntry(hT0BestVsRun, "", "lp");
  //
  DrawStdCanvas(hT0TOFVsRun, hT0T0ACVsRun, hT0T0AVsRun, hT0T0CVsRun, hT0BestVsRun, -100., 100., leg, "cStartTimeSummary;Start Time by different methods;Start Time (ps)", kTRUE);

  leg = GetLeg(0.6, 0.7, 0.9, 0.9);
  leg->SetNColumns(2);
  //
  leg->AddEntry(hT0TOFVsRunRes, "TOF_T0 res.", "lp");
  leg->AddEntry(hT0T0ACVsRunRes, "T0AC_T0 res.", "lp");
  leg->AddEntry(hT0T0AVsRunRes, "T0A_T0 res.", "lp");
  leg->AddEntry(hT0T0CVsRunRes, "T0C_T0 res.", "lp");
  leg->AddEntry(hT0BestVsRunRes, "Best_T0 res.", "lp");
  //
  DrawStdCanvas(hT0TOFVsRunRes, hT0T0ACVsRunRes, hT0T0AVsRunRes, hT0T0CVsRunRes, hT0BestVsRunRes, 0., 200., leg, "cStartTimeResolutionSummary;Start Time Resolution by different methods;#sigma Start Time (ps)", kTRUE);

  //Plot active channels trend
  DrawStdCanvas(hGoodChannelsRatio, 0.75, 1., "cGoodCh");
  DrawStdCanvas(hGoodChannelsRatioInAcc, 0.75, 1., "cGoodChInAcc");

  //Plot PID performance
  TCanvas* cPidPerformance = new TCanvas("cPidPerformance", "summary of PID performance", 1200, 500);
  cPidPerformance->Divide(3, 1);
  cPidPerformance->cd(1);
  gPad->SetLogz();
  hDiffTimePi->Draw("colz");
  cPidPerformance->cd(2);
  gPad->SetLogz();
  hDiffTimeKa->Draw("colz");
  cPidPerformance->cd(3);
  gPad->SetLogz();
  hDiffTimePro->Draw("colz");
  cPidPerformance->Print(Form("%s/cPIDExpTimes.%s", plotDir.Data(), plotExt.Data()));

  if (displayAll) {
    DrawStdCanvas(hPeakT0AVsRun, -300, 300);
    DrawStdCanvas(hPeakT0CVsRun, -300, 300);
    DrawStdCanvas(hPeakT0ACVsRun, -300, 300);
    DrawStdCanvas(hT0fillResVsRun, -10000, 10000);

    //Plot TOF signal trend
    DrawStdCanvas(hAvDiffTimeVsRun, 0, 2000);
    DrawStdCanvas(hAvTimeVsRun, 0, 25);
    DrawStdCanvas(hPeakTimeVsRun, 0, 30);
    DrawStdCanvas(hSpreadTimeVsRun, -2, 2);
    DrawStdCanvas(hAvRawTimeVsRun, 0, 300);
    DrawStdCanvas(hPeakRawTimeVsRun, 0, 300);
    DrawStdCanvas(hSpreadRawTimeVsRun, 0, 30);
    DrawStdCanvas(hAvTotVsRun, 0, 20);
    DrawStdCanvas(hPeakTotVsRun, 0, 20);
    DrawStdCanvas(hSpreadTotVsRun, 0, 2);
    DrawStdCanvas(hNegTimeRatioVsRun, 0, 0.05);
    DrawStdCanvas(hOrphansRatioVsRun, 0, 0.05);
    DrawStdCanvas(hMeanLVsRun, 0, 500);
    DrawStdCanvas(hNegLRatioVsRun, 0, 1.1);
  }
  return 0;
}
//
void DrawStdCanvas(TH1* h1, TH1* h2, TH1* h3, TH1* h4, TH1* h5, Double_t yMin, Double_t yMax, TLegend* leg, const TString cname, Bool_t Gridy)
{
  const Bool_t verbose = kFALSE;
  TString hname = h1->GetName();
  hname.Replace(0, 1, "c");
  TString ctitle = hname;
  if (!cname.IsNull() && !cname.IsWhitespace() && cname.Sizeof() > 2) {
    hname = cname;
    if (cname.Contains(";")) {
      TObjArray* arr = cname.Tokenize(";");
      hname = arr->At(0)->GetName();
      ctitle = arr->At(1)->GetName();
      if (arr->GetEntries() > 2) {
        const TString ytitle = arr->At(2)->GetName();
        if (!ytitle.IsNull())
          ctitle += ";" + ytitle;
      }
    }
  }
  TCanvas* cDisplay = new TCanvas(hname, ctitle, 50, 50, 1050, 550);
  if (verbose)
    ::Info("DrawStdCanvas", "\nDrawing in standard canvas '%s' (title: '%s')\n with asked custom name '%s': \n\th1:'%s' \n\th2:'%s' \n\th3:'%s' \n\th4:'%s' \n\th5:'%s'", hname.Data(), ctitle.Data(), cname.Data(), h1->GetName(), h2 ? h2->GetName() : "", h3 ? h3->GetName() : "", h4 ? h4->GetName() : "", h5 ? h5->GetName() : "");
  //
  h1->Draw();
  if (yMin < yMax)
    h1->GetYaxis()->SetRangeUser(yMin, yMax);
  if (h2)
    h2->Draw("same");
  if (h3)
    h3->Draw("same");
  if (h4)
    h4->Draw("same");
  if (h5)
    h5->Draw("same");
  if (leg)
    leg->Draw("same");
  if (Gridy)
    cDisplay->SetGridy();
  //
  cDisplay->Modified();
  cDisplay->Update();
  cDisplay->Print(Form("%s/%s.%s", plotDir.Data(), cDisplay->GetName(), plotExt.Data()));
}
//
void DrawStdCanvas(TH1* h1, TH1* h2, TH1* h3, Double_t yMin, Double_t yMax, TLegend* leg, TString cname)
{
  DrawStdCanvas(/*TH1*/ h1, /*TH1*/ h2, /*TH1*/ h3, /*TH1*/ nullptr, /*TH1*/ nullptr, /*Double_t*/ yMin, /*Double_t*/ yMax, /*TLegend*/ leg, /*TString*/ cname, /*Bool_t*/ kFALSE);
}
//
void DrawStdCanvas(TH1* h1, TH1* h2, Double_t yMin, Double_t yMax, TLegend* leg, TString cname)
{
  DrawStdCanvas(/*TH1*/ h1, /*TH1*/ h2, /*TH1*/ nullptr, /*Double_t*/ yMin, /*Double_t*/ yMax, /*TLegend*/ leg, /*TString*/ cname);
}
//
void DrawStdCanvas(TH1* h1, Double_t yMin, Double_t yMax, TString cname)
{
  DrawStdCanvas(/*TH1*/ h1, /*TH1*/ nullptr, /*Double_t*/ yMin, /*Double_t*/ yMax, /*TLegend*/ nullptr, /*TString*/ cname);
}
//
void SetPlotExt()
{
  const TString desiredext = gSystem->Getenv("TOFQAPLOTEXTENSION");
  if (desiredext.EqualTo("pdf") || desiredext.EqualTo("root"))
    plotExt = desiredext;
  else if (!desiredext.IsNull())
    ::Error("SetPlotExt", "Unrecognized extension: '%s'", desiredext.Data());
}
//
TLegend* GetLeg(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  return leg;
}
