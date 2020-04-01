/*
 * AliFemtoDreamTrackHist.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamTrackHist.h"
#include "TMath.h"
ClassImp(AliFemtoDreamTrackHist)
AliFemtoDreamTrackHist::AliFemtoDreamTrackHist()
    : fpTmin(0),
      fpTmax(0),
      fMinimalBooking(false),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fHistList(0),
      fConfig(0),
      fCutCounter(0),
      fDCAXYPtBins(0),
      fTOFMass(0),
      fNSigCom(0) {
  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = nullptr;
    fpTDist[i] = nullptr;
    fpDist[i] = nullptr;
    fpTPCDist[i] = nullptr;
    fetaDist[i] = nullptr;
    fphiDist[i] = nullptr;
    fTPCCls[i] = nullptr;
    fTPCClsS[i] = nullptr;
    fTrackChi2[i] = nullptr;
    fShrdClsITS[i] = nullptr;
    fDCAxy[i] = nullptr;
    fDCAz[i] = nullptr;
    fDCAxyProp[i] = nullptr;
    fDCAzProp[i] = nullptr;
    fTPCCrossedRows[i] = nullptr;
    fTPCRatio[i] = nullptr;
    fITSdedx[i] = nullptr;
    fTPCdedx[i] = nullptr;
    fTOFbeta[i] = nullptr;
    fNSigITS[i] = nullptr;
    fNSigITSMod[i] = nullptr;
    fNSigTPC[i] = nullptr;
    fNSigTPCMod[i] = nullptr;
    fNSigTOF[i] = nullptr;
    fITSStatus[i] = nullptr;
    fTPCStatus[i] = nullptr;
    fTOFStatus[i] = nullptr;
    fNSigComITSTPC[i] = nullptr;
    fNSigComTPCTOF[i] = nullptr;
    fTPCClsCPiluUp[i] = nullptr;
    fITShrdClsPileUp[i] = nullptr;
  }
}
AliFemtoDreamTrackHist::AliFemtoDreamTrackHist(bool DCADist, bool CombSig, bool TOFM, float pTmin, float pTmax)
    : fpTmin(pTmin),
      fpTmax(pTmax),
      fMinimalBooking(false),
      fMultRangeLow(27),
      fMultRangeHigh(55) {
  TString sName[2] = { "before", "after" };
  float ptmin = 0;
  float ptmax = 6.0;
  int ptBins = 120;
  int twoDBins = 400;
  fHistList = new TList();
  fHistList->SetName("TrackCuts");
  fHistList->SetOwner();

  fConfig = new TProfile("TrackCutConfig", "Track Cut Config", 30, 0, 30);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1, "pT_{min}");
  fConfig->GetXaxis()->SetBinLabel(2, "pT_{max}");
  fConfig->GetXaxis()->SetBinLabel(3, "pT_{excl,min}");
  fConfig->GetXaxis()->SetBinLabel(4, "pT_{excl,max}");
  fConfig->GetXaxis()->SetBinLabel(5, "#eta_{min}");
  fConfig->GetXaxis()->SetBinLabel(6, "#eta_{max}");
  fConfig->GetXaxis()->SetBinLabel(7, "Charge");
  fConfig->GetXaxis()->SetBinLabel(8, "ClsTPC_{min}");
  fConfig->GetXaxis()->SetBinLabel(9, "Filter Bit");
  fConfig->GetXaxis()->SetBinLabel(10, "Is MC");
  fConfig->GetXaxis()->SetBinLabel(11, "DCA_{xy}");
  fConfig->GetXaxis()->SetBinLabel(12, "DCA_{z}");
  fConfig->GetXaxis()->SetBinLabel(13, "Max Shared Cls");
  fConfig->GetXaxis()->SetBinLabel(14, "Shared Cls");
  fConfig->GetXaxis()->SetBinLabel(15, "TPC Crossed Rows");
  fConfig->GetXaxis()->SetBinLabel(16, "Ratio TPC Crossed Rows");
  fConfig->GetXaxis()->SetBinLabel(17, "PID Mom Thresh");
  fConfig->GetXaxis()->SetBinLabel(18, "PID nsig");
  fConfig->GetXaxis()->SetBinLabel(19, "PID ITS nsig");
  fConfig->GetXaxis()->SetBinLabel(20, "Reject Pions");
  fConfig->GetXaxis()->SetBinLabel(21, "Smallest Sig");
  fConfig->GetXaxis()->SetBinLabel(22, "ITS Hit");
  fConfig->GetXaxis()->SetBinLabel(23, "SPD Hit TOF Timing");
  fConfig->GetXaxis()->SetBinLabel(24, "TOF Timing");
  fConfig->GetXaxis()->SetBinLabel(25, "Pile Up Rej");
  fConfig->GetXaxis()->SetBinLabel(26, "TPC Refit");
  fConfig->GetXaxis()->SetBinLabel(27, "#chi^{2} min");
  fConfig->GetXaxis()->SetBinLabel(28, "#chi^{2} max");
  fConfig->GetXaxis()->SetBinLabel(29, "ESDFiltering");



  fHistList->Add(fConfig);

  fCutCounter = new TH1F("CutCounter", "Cut Counter", 35, 0, 35);
  fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
  fCutCounter->GetXaxis()->SetBinLabel(2, "Filter Bit");
  fCutCounter->GetXaxis()->SetBinLabel(3, "p_{T} Cut");
  fCutCounter->GetXaxis()->SetBinLabel(4, "p_{T, excl} Cut");
  fCutCounter->GetXaxis()->SetBinLabel(5, "#eta Cut");
  fCutCounter->GetXaxis()->SetBinLabel(6, "Charge");
  fCutCounter->GetXaxis()->SetBinLabel(7, "PileUpITS");
  fCutCounter->GetXaxis()->SetBinLabel(8, "PileUpSPDTOF");
  fCutCounter->GetXaxis()->SetBinLabel(9, "PileUpTOF");
  fCutCounter->GetXaxis()->SetBinLabel(10, "PileUp");
  fCutCounter->GetXaxis()->SetBinLabel(11, "nClsTPC");
  fCutCounter->GetXaxis()->SetBinLabel(12, "Max Shared Cls TPC");
  fCutCounter->GetXaxis()->SetBinLabel(13, "Shared Cls");
  fCutCounter->GetXaxis()->SetBinLabel(14, "TPC Refit");
  fCutCounter->GetXaxis()->SetBinLabel(15, "TPC Crossed Rows");
  fCutCounter->GetXaxis()->SetBinLabel(16, "TPC Row Ratio");
  fCutCounter->GetXaxis()->SetBinLabel(17, "#chi^{2} OK");
  fCutCounter->GetXaxis()->SetBinLabel(18, "ITS OK");
  fCutCounter->GetXaxis()->SetBinLabel(19, "TPC OK");
  fCutCounter->GetXaxis()->SetBinLabel(20, "TOF OK");
  fCutCounter->GetXaxis()->SetBinLabel(21, "TPC TOF OK");
  fCutCounter->GetXaxis()->SetBinLabel(22, "ITS PID");
  fCutCounter->GetXaxis()->SetBinLabel(23, "TPC PID");
  fCutCounter->GetXaxis()->SetBinLabel(24, "TPC TOF PID");
  fCutCounter->GetXaxis()->SetBinLabel(25, "Reject Pions");
  fCutCounter->GetXaxis()->SetBinLabel(26, "Smallest Sig");
  fCutCounter->GetXaxis()->SetBinLabel(27, "Passes PID");
  fCutCounter->GetXaxis()->SetBinLabel(28, "Passes ITS_d PID");
  fCutCounter->GetXaxis()->SetBinLabel(29, "DCA_{Z}");
  fCutCounter->GetXaxis()->SetBinLabel(30, "DCA_{XY}");
  fCutCounter->GetXaxis()->SetBinLabel(31, "ITS_d status Ok");
  fCutCounter->GetXaxis()->SetBinLabel(32, "ITS_d PID");//_d stands for deuteron analysis
  fCutCounter->GetYaxis()->SetTitle("Entries");

  fHistList->Add(fCutCounter);

  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = new TList();
    fTrackCutQA[i]->SetName(sName[i].Data());
    fTrackCutQA[i]->SetOwner();
    fHistList->Add(fTrackCutQA[i]);

    TString ptName = Form("pTDist_%s", sName[i].Data());
    fpTDist[i] = new TH1F(ptName.Data(), ptName.Data(), 2. * ptBins, ptmin, ptmax);
    fpTDist[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fpTDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fpTDist[i]);

    TString pITSName = Form("pDist_%s", sName[i].Data());
    fpDist[i] = new TH1F(pITSName.Data(), pITSName.Data(), ptBins, ptmin,
                            ptmax);
    fpDist[i]->GetXaxis()->SetTitle("#it{p}_{} (GeV/#it{c})");
    fpDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fpDist[i]);

    TString pTPCName = Form("pTPCDist_%s", sName[i].Data());
    fpTPCDist[i] = new TH1F(pTPCName.Data(), pTPCName.Data(), ptBins, ptmin,
                            ptmax);
    fpTPCDist[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fpTPCDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fpTPCDist[i]);

    TString etaName = Form("EtaDist_%s", sName[i].Data());
    fetaDist[i] = new TH1F(etaName.Data(), etaName.Data(), 200, -1.5, 1.5);
    fetaDist[i]->GetXaxis()->SetTitle("#eta");
    fetaDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fetaDist[i]);

    TString phiName = Form("phiDist_%s", sName[i].Data());
    fphiDist[i] = new TH1F(phiName.Data(), phiName.Data(), 200, 0.,
                           2 * TMath::Pi());
    fphiDist[i]->GetXaxis()->SetTitle("#phi");
    fphiDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fphiDist[i]);

    TString TPCName = Form("TPCCls_%s", sName[i].Data());
    fTPCCls[i] = new TH1F(TPCName.Data(), TPCName.Data(), 160, 0, 160);
    fTPCCls[i]->GetXaxis()->SetTitle("#it{n}_{cls} (TPC)");
    fTPCCls[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fTPCCls[i]);

    TString DCAXYName = Form("DCAXY_%s", sName[i].Data());
    fDCAxy[i] = new TH2F(DCAXYName.Data(), DCAXYName.Data(), ptBins, ptmin,
                         ptmax, 2.5 * twoDBins, -5., 5.);
    fDCAxy[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fDCAxy[i]->GetYaxis()->SetTitle("DCA_{xy}");

    fTrackCutQA[i]->Add(fDCAxy[i]);

    TString DCAZName = Form("DCAZ_%s", sName[i].Data());
    fDCAz[i] = new TH2F(DCAZName.Data(), DCAZName.Data(), ptBins, ptmin, ptmax,
                        2.5 * twoDBins, -5., 5.);
    fDCAz[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fDCAz[i]->GetYaxis()->SetTitle("DCA_{z}");
    fTrackCutQA[i]->Add(fDCAz[i]);

    TString DCAXYPropName = Form("DCAXYProp_%s", sName[i].Data());
    fDCAxyProp[i] = new TH2F(DCAXYPropName.Data(), DCAXYPropName.Data(), ptBins,
                             ptmin, ptmax, 2.5 * twoDBins, -5., 5.);
    fDCAxyProp[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fDCAxyProp[i]->GetYaxis()->SetTitle("DCA_{xy,Prop}");

    fTrackCutQA[i]->Add(fDCAxyProp[i]);

    TString DCAZPropName = Form("DCAZProp_%s", sName[i].Data());
    fDCAzProp[i] = new TH2F(DCAZPropName.Data(), DCAZPropName.Data(), ptBins,
                            ptmin, ptmax, 2.5 * twoDBins, -5., 5.);
    fDCAzProp[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fDCAzProp[i]->GetYaxis()->SetTitle("DCA_{z,Prop}");
    fTrackCutQA[i]->Add(fDCAzProp[i]);

    TString TPCCRName = Form("CrossedRows_%s", sName[i].Data());
    fTPCCrossedRows[i] = new TH1F(TPCCRName.Data(), TPCCRName.Data(), 160, 0,
                                  160.);
    fTPCCrossedRows[i]->GetXaxis()->SetTitle("#it{n}_{crossed} (TPC)");
    fTPCCrossedRows[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fTPCCrossedRows[i]);

    TString TPCratioName = Form("TPCRatio_%s", sName[i].Data());
    fTPCRatio[i] = new TH1F(TPCratioName.Data(), TPCratioName.Data(), 150, 0.,
                            1.5);
    fTPCRatio[i]->GetXaxis()->SetTitle("#it{n}_{cluster}/#it{n}_{findable} (TPC)");
    fTPCRatio[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fTPCRatio[i]);

    TString TPCClsSName = Form("TPCClsS_%s", sName[i].Data());
    fTPCClsS[i] = new TH1F(TPCClsSName.Data(), TPCClsSName.Data(), 160, 0, 160);
    fTPCClsS[i]->GetXaxis()->SetTitle("#it{n}_{shared} (TPC)");
    fTPCClsS[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fTPCClsS[i]);

    TString ChiSquareName = Form("TrackChi2_%s", sName[i].Data());
    fTrackChi2[i] = new TH2F(ChiSquareName.Data(), ChiSquareName.Data(), ptBins,
                             ptmin, ptmax, 100, 0, 20);
    fTrackChi2[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fTrackChi2[i]->GetYaxis()->SetTitle("#chi^{2}/NDF");
    fTrackCutQA[i]->Add(fTrackChi2[i]);

    TString ShrdClsITSName = Form("SharedClsITS_%s", sName[i].Data());
    fShrdClsITS[i] = new TH2F(ShrdClsITSName.Data(), ShrdClsITSName.Data(), 6,
                              1, 7, 2, 0, 2);
    fShrdClsITS[i]->GetXaxis()->SetTitle("ITS Layer");
    fShrdClsITS[i]->GetYaxis()->SetBinLabel(1, "Has shared Cls");
    fShrdClsITS[i]->GetYaxis()->SetBinLabel(2, "Has no Shard Cls");
    fTrackCutQA[i]->Add(fShrdClsITS[i]);

    TString ITSdedxName = Form("ITSdedx_%s", sName[i].Data());
    fITSdedx[i] = new TH2F(ITSdedxName.Data(), ITSdedxName.Data(), ptBins,
                           ptmin, ptmax, 2 * twoDBins, 0., 400);
    fITSdedx[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fITSdedx[i]->GetYaxis()->SetTitle("d#it{E}/d#it{x} (arb. Units)");

    fTrackCutQA[i]->Add(fITSdedx[i]);

    TString TPCdedxName = Form("TPCdedx_%s", sName[i].Data());
    fTPCdedx[i] = new TH2F(TPCdedxName.Data(), TPCdedxName.Data(), ptBins,
                           ptmin, ptmax, 2 * twoDBins, 0., 400);
    fTPCdedx[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fTPCdedx[i]->GetYaxis()->SetTitle("d#it{E}/d#it{x} (arb. Units)");

    fTrackCutQA[i]->Add(fTPCdedx[i]);

    TString TOFbetaName = Form("TOFbeta_%s", sName[i].Data());
    fTOFbeta[i] = new TH2F(TOFbetaName.Data(), TOFbetaName.Data(), ptBins,
                           ptmin, ptmax, 12 * twoDBins, -0.1, 1.1);
    fTOFbeta[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fTOFbeta[i]->GetYaxis()->SetTitle("#beta_{TOF}");

    fTrackCutQA[i]->Add(fTOFbeta[i]);

    TString NSigITSName = Form("NSigITS_%s", sName[i].Data());
    fNSigITS[i] = new TH2F(NSigITSName.Data(), NSigITSName.Data(), ptBins,
                           ptmin, ptmax, 3. * twoDBins, -60., 60.);
    fNSigITS[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fNSigITS[i]->GetYaxis()->SetTitle("n#sigma_{ITS}");

    fTrackCutQA[i]->Add(fNSigITS[i]);

    TString NSigITSModName = Form("NSigITSMod_%s", sName[i].Data());
    fNSigITSMod[i] = new TH2F(NSigITSModName.Data(), NSigITSModName.Data(),
                              ptBins, ptmin, 1., 3. * twoDBins, -1., 5.);
    fNSigITSMod[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fNSigITSMod[i]->GetYaxis()->SetTitle("|n#sigma_{ITS}|");

    fTrackCutQA[i]->Add(fNSigITSMod[i]);


    TString NSigTPCName = Form("NSigTPC_%s", sName[i].Data());
    fNSigTPC[i] = new TH2F(NSigTPCName.Data(), NSigTPCName.Data(), ptBins,
                           ptmin, ptmax, 3. * twoDBins, -60., 60.);
    fNSigTPC[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fNSigTPC[i]->GetYaxis()->SetTitle("n#sigma_{TPC}");

    fTrackCutQA[i]->Add(fNSigTPC[i]);

    TString NSigTPCModName = Form("NSigTPCMod_%s", sName[i].Data());
    fNSigTPCMod[i] = new TH2F(NSigTPCModName.Data(), NSigTPCModName.Data(),
                              ptBins, ptmin, 1., 3. * twoDBins, -1., 5.);
    fNSigTPCMod[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fNSigTPCMod[i]->GetYaxis()->SetTitle("|n#sigma_{TPC}|");

    fTrackCutQA[i]->Add(fNSigTPCMod[i]);

    TString NSigTOFName = Form("NSigTOF_%s", sName[i].Data());
    fNSigTOF[i] = new TH2F(NSigTOFName.Data(), NSigTOFName.Data(), ptBins,
                           ptmin, ptmax, 3. * twoDBins, -60., 60.);
    fNSigTOF[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fNSigTOF[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");

    fTrackCutQA[i]->Add(fNSigTOF[i]);

    TString ITSstatusname = Form("ITSStatus_%s", sName[i].Data());
    fITSStatus[i] = new TH1F(ITSstatusname.Data(), ITSstatusname.Data(), 4, 0,
                             4);
    fITSStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
    fITSStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
    fITSStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
    fITSStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
    fITSStatus[i]->GetYaxis()->SetTitle("Entries");

    fTrackCutQA[i]->Add(fITSStatus[i]);

    TString TPCstatusname = Form("TPCStatus_%s", sName[i].Data());
    fTPCStatus[i] = new TH1F(TPCstatusname.Data(), TPCstatusname.Data(), 4, 0,
                             4);
    fTPCStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
    fTPCStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
    fTPCStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
    fTPCStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
    fTPCStatus[i]->GetYaxis()->SetTitle("Entries");

    fTrackCutQA[i]->Add(fTPCStatus[i]);

    TString TOFstatusname = Form("TOFStatus_%s", sName[i].Data());
    fTOFStatus[i] = new TH1F(TOFstatusname.Data(), TOFstatusname.Data(), 4, 0,
                             4);
    fTOFStatus[i]->GetXaxis()->SetBinLabel(1, "kDetNoSignal");
    fTOFStatus[i]->GetXaxis()->SetBinLabel(2, "kDetPidOk");
    fTOFStatus[i]->GetXaxis()->SetBinLabel(3, "kDetMismatch");
    fTOFStatus[i]->GetXaxis()->SetBinLabel(4, "kDetNoParams");
    fTOFStatus[i]->GetYaxis()->SetTitle("Entries");

    fTrackCutQA[i]->Add(fTOFStatus[i]);

    TString NSigComITSTPCName = Form("NSigComITSTPC_%s", sName[i].Data());
    fNSigComITSTPC[i] = new TH2F(NSigComITSTPCName.Data(),
                                 NSigComITSTPCName.Data(), ptBins, ptmin, 7,
                                 3. * twoDBins, -1., 7.);
    fNSigComITSTPC[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fNSigComITSTPC[i]->GetYaxis()->SetTitle(
        "n#sigma_{comb}=#sqrt{n#sigma_{ITS}^{2}+n#sigma_{TPC}^{2}}");

    fTrackCutQA[i]->Add(fNSigComITSTPC[i]);

    TString NSigComTPCTOFName = Form("NSigComTPCTOF_%s", sName[i].Data());
    fNSigComTPCTOF[i] = new TH2F(NSigComTPCTOFName.Data(),
                                 NSigComTPCTOFName.Data(), ptBins, ptmin, 7,
                                 3. * twoDBins, -1., 7.);
    fNSigComTPCTOF[i]->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
    fNSigComTPCTOF[i]->GetYaxis()->SetTitle(
        "n#sigma_{comb}=#sqrt{n#sigma_{TPC}^{2}+n#sigma_{TOF}^{2}}");

    fTrackCutQA[i]->Add(fNSigComTPCTOF[i]);

    TString TPCClsCPileUpName = Form("TPCClsCPileUp_%s", sName[i].Data());
    fTPCClsCPiluUp[i] = new TH2F(TPCClsCPileUpName.Data(),
                                 TPCClsCPileUpName.Data(), 15, 0, 15, 160, 0,
                                 160);
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(1, "Hit ITS Layer 1");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(2, "Hit ITS Layer 2");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(3, "Hit ITS Layer 3");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(4, "Hit ITS Layer 4");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(5, "Hit ITS Layer 5");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(6, "Hit ITS Layer 6");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(7, "TOF Bunch.Xing");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(8, "no Hit ITS Layer 1");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(9, "no Hit ITS Layer 2");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(10, "no Hit ITS Layer 3");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(11, "no Hit ITS Layer 4");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(12, "no Hit ITS Layer 5");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(13, "no Hit ITS Layer 6");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(14, "no TOF Bunch.Xing");
    fTPCClsCPiluUp[i]->GetXaxis()->SetBinLabel(15, "all ITS-TOF Req false");
    fTrackCutQA[i]->Add(fTPCClsCPiluUp[i]);

    TString ITSShrdClsPileUpName = Form("ITSShrdClsPileUp_%s", sName[i].Data());
    fITShrdClsPileUp[i] = new TH2F(ITSShrdClsPileUpName.Data(),
                                   ITSShrdClsPileUpName.Data(), 15, 0, 15, 12,
                                   0, 12);
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(1, "Hit ITS Layer 1");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(2, "Hit ITS Layer 2");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(3, "Hit ITS Layer 3");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(4, "Hit ITS Layer 4");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(5, "Hit ITS Layer 5");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(6, "Hit ITS Layer 6");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(7, "TOF Bunch.Xing");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(8, "no Hit ITS Layer 1");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(9, "no Hit ITS Layer 2");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(10, "no Hit ITS Layer 3");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(11, "no Hit ITS Layer 4");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(12, "no Hit ITS Layer 5");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(13, "no Hit ITS Layer 6");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(14, "no Hit TOF Bunch.Xing");
    fITShrdClsPileUp[i]->GetXaxis()->SetBinLabel(15,
                                                 "all Hit ITS-TOF Req false");

    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(1, "Has shared Cls Layer 1");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(2, "Has shared Cls Layer 2");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(3, "Has shared Cls Layer 3");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(4, "Has shared Cls Layer 4");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(5, "Has shared Cls Layer 5");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(6, "Has shared Cls Layer 6");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(7, "Has no Shard Cls 1");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(8, "Has no Shard Cls 2");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(9, "Has no Shard Cls 3");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(10, "Has no Shard Cls 4");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(11, "Has no Shard Cls 5");
    fITShrdClsPileUp[i]->GetYaxis()->SetBinLabel(12, "Has no Shard Cls 6");
    fTrackCutQA[i]->Add(fITShrdClsPileUp[i]);
  }
  if (CombSig) {
    fNSigCom = new TH3F("NSigComb", "NSigComb", 14, 0.5, 4.0, 300, -30, 30, 300,
                        -30, 30);
    fHistList->Add(fNSigCom);
  } else {
    fNSigCom = 0;
  }
  if (DCADist) {
    TString dcaPtBinName = Form("DCAXYPtBinningTot");
    fDCAXYPtBins = new TH2F(dcaPtBinName.Data(), dcaPtBinName.Data(), 20, pTmin,
                            pTmax, 500, -5, 5);
    fDCAXYPtBins->GetXaxis()->SetTitle("P#_{T}");
    fDCAXYPtBins->GetYaxis()->SetTitle("dca_{XY}");
    fHistList->Add(fDCAXYPtBins);

    dcaPtBinName = Form("DCAXYPtMult_0_%i", fMultRangeLow);
    fDCAXYPtBinsMult[0] = new TH2F(
        dcaPtBinName.Data(),
        Form("0 < mult < %i;P#_{T};dca_{XY}", fMultRangeLow), 20, 0.5, 4.05,
        500, -5, 5);

    dcaPtBinName = Form("DCAXYPtMult_%i_%i", fMultRangeLow, fMultRangeHigh);
    fDCAXYPtBinsMult[1] = new TH2F(
        dcaPtBinName.Data(),
        Form("%i < mult < %i;P#_{T};dca_{XY}", fMultRangeLow, fMultRangeHigh),
        20, 0.5, 4.05, 500, -5, 5);

    dcaPtBinName = Form("DCAXYPtMult_%i_inf", fMultRangeHigh);
    fDCAXYPtBinsMult[2] = new TH2F(
        dcaPtBinName.Data(), Form("mult > %i;P#_{T};dca_{XY}", fMultRangeHigh),
        20, 0.5, 4.05, 500, -5, 5);

    fHistList->Add(fDCAXYPtBinsMult[0]);
    fHistList->Add(fDCAXYPtBinsMult[1]);
    fHistList->Add(fDCAXYPtBinsMult[2]);
  } else {
    fDCAXYPtBins = 0;
  }

  if (TOFM)
  {
  TString TOFMassName = Form("TOFMass");
  fTOFMass = new TH2F(TOFMassName.Data(), TOFMassName.Data(), ptBins,
                         ptmin, ptmax, 1400, 0., 1.1);
  fTOFMass->GetXaxis()->SetTitle("p_{primary}");
  fTOFMass->GetYaxis()->SetTitle("m_{TOF}");
  fHistList->Add(fTOFMass);
  } else {
      fTOFMass = nullptr;
  }

}

AliFemtoDreamTrackHist::AliFemtoDreamTrackHist(TString MinimalBooking)
    : fMinimalBooking(true),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fConfig(0),
      fCutCounter(0),
      fDCAXYPtBins(0),
      fTOFMass(0),
      fNSigCom(0) {
  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = nullptr;
    fpDist[i] = nullptr;
    fpTPCDist[i] = nullptr;
    fetaDist[i] = nullptr;
    fphiDist[i] = nullptr;
    fTPCCls[i] = nullptr;
    fTPCClsS[i] = nullptr;
    fShrdClsITS[i] = nullptr;
    fDCAxy[i] = nullptr;
    fDCAz[i] = nullptr;
    fDCAxyProp[i] = nullptr;
    fDCAzProp[i] = nullptr;
    fTPCCrossedRows[i] = nullptr;
    fTPCRatio[i] = nullptr;
    fITSdedx[i] = nullptr;
    fTPCdedx[i] = nullptr;
    fTOFbeta[i] = nullptr;
    fNSigITS[i] = nullptr;
    fNSigITSMod[i] = nullptr;
    fNSigTPC[i] = nullptr;
    fNSigTPCMod[i] = nullptr;
    fNSigTOF[i] = nullptr;
    fITSStatus[i] = nullptr;
    fTPCStatus[i] = nullptr;
    fTOFStatus[i] = nullptr;
    fNSigComITSTPC[i] = nullptr;
    fNSigComTPCTOF[i] = nullptr;
    fTPCClsCPiluUp[i] = nullptr;
    fITShrdClsPileUp[i] = nullptr;
  }
  fHistList = new TList();
  fHistList->SetOwner();
  fHistList->SetName(MinimalBooking.Data());

  fpTDist[0] = 0;
  TString ptName = Form("pTDist_%s", "after");
  fpTDist[1] = new TH1F(ptName.Data(), ptName.Data(), 200, 0, 5);
  fpTDist[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fpTDist[1]->GetYaxis()->SetTitle("Entries");
  fHistList->Add(fpTDist[1]);
}

AliFemtoDreamTrackHist::~AliFemtoDreamTrackHist() {
  delete fHistList;
}

void AliFemtoDreamTrackHist::FillNSigComb(float pT, float nSigTPC,
                                          float nSigTOF) {
  if (!fMinimalBooking)
    fNSigCom->Fill(pT, nSigTPC, nSigTOF);
}

void AliFemtoDreamTrackHist::FillDCAXYPtBins(float pT, float dcaxy,
                                             int multiplicity) {
  if (!fMinimalBooking) {
    fDCAXYPtBins->Fill(pT, dcaxy);
    if (multiplicity < fMultRangeLow) {
      fDCAXYPtBinsMult[0]->Fill(pT, dcaxy);
    } else if (multiplicity >= fMultRangeLow && multiplicity < fMultRangeHigh) {
      fDCAXYPtBinsMult[1]->Fill(pT, dcaxy);
    } else {
      fDCAXYPtBinsMult[2]->Fill(pT, dcaxy);
    }
  }
}
