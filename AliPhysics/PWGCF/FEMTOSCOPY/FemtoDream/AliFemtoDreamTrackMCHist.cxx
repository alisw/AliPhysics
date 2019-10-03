/*
 * AliFemtoDreamTrackMCHist.cxx
 *
 *  Created on: Nov 14,2017
 *      Author: gu74req
 */

#include "AliFemtoDreamTrackMCHist.h"
#include "AliLog.h"
AliFemtoDreamTrackMCHist::AliFemtoDreamTrackMCHist()
    : fpTmin(0.5),
      fpTmax(4.05),
      fpTbins(20),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fDoSplitting(false),
      fDoDCAPlots(false),
      fDoMultiplicityBinning(false),
      fMCList(nullptr),
      fDCAPlots(nullptr),
      fMCCorrPt(nullptr),
      fMCIdentPt(nullptr),
      fMCGenPt(nullptr),
      fMCContPt(nullptr),
      fMCUnknownPt(nullptr),
      fMCPrimaryPt(nullptr),
      fMCMaterialPt(nullptr),
      fMCFeeddownWeakPt(nullptr),
      fHistMCMother(nullptr),
      fHistMCMotherPDG(nullptr),
      fMCPrimDCAXYPtBins(nullptr),
      fMCMaterialDCAXYPtBins(nullptr),
      fMCSecondaryDCAXYPtBins(nullptr),
      fMCSecLambdaDCAXYPtBins(nullptr),
      fMCSecSigmaDCAXYPtBins(nullptr),
      fPtResolution(nullptr),
      fThetaResolution(nullptr),
      fPhiResolution(nullptr) {
  for (int i = 0; i < 4; ++i) {
    if (i < 3) {
      fMCPrimDCAXYPtBinsMult[i] = 0;
      fMCMaterialDCAXYPtBinsMult[i] = 0;
      fMCSecondaryDCAXYPtBinsMult[i] = 0;
      fMCSecLambdaDCAXYPtBinsMult[i] = 0;
      fMCSecSigmaDCAXYPtBinsMult[i] = 0;
    }
    fMCQAPlots[i] = 0;
    fMCpTPCDist[i] = 0;
    fMCetaDist[i] = 0;
    fMCphiDist[i] = 0;
    fMCTPCCls[i] = 0;
    fMCDCAxy[i] = 0;
    fMCDCAz[i] = 0;
    fMCTPCCrossedRows[i] = 0;
    fMCTPCRatio[i] = 0;
    fMCTPCdedx[i] = 0;
    fMCTOFbeta[i] = 0;
    fMCNSigTPC[i] = 0;
    fMCNSigTOF[i] = 0;
  }
}

AliFemtoDreamTrackMCHist::~AliFemtoDreamTrackMCHist() {
  if (fMCList) {
    delete fMCList;
  }
}

AliFemtoDreamTrackMCHist::AliFemtoDreamTrackMCHist(bool contribSplitting,
                                                   bool DCADist,
                                                   bool DoMultBinning,
                                                   bool checkMother)
    : fpTmin(0.5),
      fpTmax(4.05),
      fpTbins(20),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fDoSplitting(contribSplitting),
      fDoDCAPlots(DCADist),
      fDoMultiplicityBinning(DoMultBinning),
      fMCList(nullptr),
      fDCAPlots(nullptr),
      fMCCorrPt(nullptr),
      fMCIdentPt(nullptr),
      fMCGenPt(nullptr),
      fMCContPt(nullptr),
      fMCUnknownPt(nullptr),
      fMCPrimaryPt(nullptr),
      fMCMaterialPt(nullptr),
      fMCFeeddownWeakPt(nullptr),
      fHistMCMother(nullptr),
      fHistMCMotherPDG(nullptr),
      fMCPrimDCAXYPtBins(nullptr),
      fMCMaterialDCAXYPtBins(nullptr),
      fMCSecondaryDCAXYPtBins(nullptr),
      fMCSecLambdaDCAXYPtBins(nullptr),
      fMCSecSigmaDCAXYPtBins(nullptr),
      fPtResolution(nullptr),
      fThetaResolution(nullptr),
      fPhiResolution(nullptr) {
  float ptmin = 0;
  float ptmax = 5;
  float ptBins = 200;
  float twoDBins = 400;
  fMCList = new TList();
  fMCList->SetName("MonteCarlo");
  fMCList->SetOwner();

  fMCCorrPt = new TH1F("CorrParPt", "Correct Particles Pt", ptBins, ptmin,
                       ptmax);
  fMCCorrPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCCorrPt);

  fMCIdentPt = new TH1F("IdentPartPt", "Ident Particles Pt", ptBins, ptmin,
                        ptmax);
  fMCIdentPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCIdentPt);

  fMCGenPt = new TH1F("GenPartPt", "Gen Particles Pt", ptBins, ptmin, ptmax);
  fMCGenPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCGenPt);

  fPtResolution = new TH2F("DeltaPtRecoTruevsPtReco", "DeltaPtRecoTruevsPtReco",
                           100, 0, 5, 750, -3, 1);
  fPtResolution->GetXaxis()->SetTitle("P_{T,True}");
  fPtResolution->GetYaxis()->SetTitle("(P_{T,True}-P_{T,Reco})/P_{T,True}");
  fMCList->Add(fPtResolution);

  fThetaResolution = new TH2F("DeltaThetaRecoTruevsPtReco",
                              "DeltaThetaRecoTruevsPtReco", 100, 0, 5, 400,
                              -0.25, 0.25);
  fThetaResolution->GetXaxis()->SetTitle("P_{T,True}");
  fThetaResolution->GetYaxis()->SetTitle("#Theta_{T,True}-#Theta_{T,Reco}");
  fMCList->Add(fThetaResolution);

  fPhiResolution = new TH2F("DeltaPhiRecoTruevsPtReco",
                            "DeltaPhiRecoTruevsPtReco", 100, 0, 5, 200, -0.2,
                            0.2);
  fPhiResolution->GetXaxis()->SetTitle("P_{T,True}");
  fPhiResolution->GetYaxis()->SetTitle("(#Phi_{T,True}-#Phi_{T,Reco})");
  fMCList->Add(fPhiResolution);

  if (checkMother) {
    fHistMCMother = new TH2F("fHistMCMother",
                             "; #it{p}_{T} (GeV/#it{c}^{2}); PDG code mother",
                             100, 0., 10., 4000, 0, 4000);
    fMCList->Add(fHistMCMother);

    fHistMCMotherPDG = new TH1I("fHistMCMotherPDG",
                                "; Entries; PDG code mother", 10000000, 0,
                                10000000);
    fMCList->Add(fHistMCMotherPDG);
  }

  if (contribSplitting) {
    fMCContPt = new TH1F("ContPt", "ContPt", ptBins, ptmin, ptmax);
    fMCContPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCContPt);

    fMCUnknownPt = new TH1F("UnknPt", "UnknPt", ptBins, ptmin, ptmax);
    fMCUnknownPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCUnknownPt);

    fMCPrimaryPt = new TH1F("PrimaryPt", "PrimaryPt", ptBins, ptmin, ptmax);
    fMCPrimaryPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCPrimaryPt);

    fMCMaterialPt = new TH1F("MatPt", "MatPT", ptBins, ptmin, ptmax);
    fMCMaterialPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCMaterialPt);

    fMCFeeddownWeakPt = new TH2F("FeeddownPt", "Feeddown Pt", ptBins, ptmin,
                                 ptmax, 213, 3121, 3334);
    fMCFeeddownWeakPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCFeeddownWeakPt);

    TString MCModes[4] = { "Primary", "Secondary", "Material", "Contamination" };
    for (int i = 0; i < 4; ++i) {
      fMCQAPlots[i] = new TList();
      fMCQAPlots[i]->SetName(MCModes[i].Data());
      fMCQAPlots[i]->SetOwner();
      fMCList->Add(fMCQAPlots[i]);

      TString MCpTPCName = Form("MCpTPCDist%s", MCModes[i].Data());
      TString MCetaName = Form("MCEtaDist%s", MCModes[i].Data());
      TString MCphiName = Form("MCphiDist%s", MCModes[i].Data());
      TString MCTPCName = Form("MCTPCCls%s", MCModes[i].Data());
      TString MCDCAXYName = Form("MCMCDCAXY%s", MCModes[i].Data());
      TString MCDCAZName = Form("MCDCAZ%s", MCModes[i].Data());
      TString MCTPCCRName = Form("MCCrossedRows%s", MCModes[i].Data());
      TString MCTPCratioName = Form("MCTPCRatio%s", MCModes[i].Data());
      TString MCTPCdedxName = Form("MCTPCdedx%s", MCModes[i].Data());
      TString MCTOFbetaName = Form("MCTOFbeta%s", MCModes[i].Data());
      TString MCNSigTPCName = Form("MCNSigTPC%s", MCModes[i].Data());
      TString MCNSigTOFName = Form("MCNSigTOF%s", MCModes[i].Data());

      fMCpTPCDist[i] = new TH1F(MCpTPCName.Data(), MCpTPCName.Data(), ptBins,
                                ptmin, ptmax);
      fMCQAPlots[i]->Add(fMCpTPCDist[i]);

      fMCetaDist[i] = new TH1F(MCetaName.Data(), MCetaName.Data(), 200, -10.,
                               10.);
      fMCQAPlots[i]->Add(fMCetaDist[i]);

      fMCphiDist[i] = new TH1F(MCphiName.Data(), MCphiName.Data(), 200, 0.,
                               2 * TMath::Pi());
      fMCQAPlots[i]->Add(fMCphiDist[i]);

      fMCTPCCls[i] = new TH2F(MCTPCName.Data(), MCTPCName.Data(), ptBins, ptmin,
                              ptmax, 100, 0, 200.);
      fMCQAPlots[i]->Add(fMCTPCCls[i]);

      fMCDCAxy[i] = new TH2F(MCDCAXYName.Data(), MCDCAXYName.Data(), ptBins,
                             ptmin, ptmax, 2.5 * twoDBins, -5., 5.);
      fMCQAPlots[i]->Add(fMCDCAxy[i]);

      fMCDCAz[i] = new TH2F(MCDCAZName.Data(), MCDCAZName.Data(), ptBins, ptmin,
                            ptmax, 2.5 * twoDBins, -5., 5.);
      fMCQAPlots[i]->Add(fMCDCAz[i]);

      fMCTPCCrossedRows[i] = new TH2F(MCTPCCRName.Data(), MCTPCCRName.Data(),
                                      ptBins, ptmin, ptmax, 100, 0, 200.);
      fMCQAPlots[i]->Add(fMCTPCCrossedRows[i]);

      fMCTPCRatio[i] = new TH2F(MCTPCratioName.Data(), MCTPCratioName.Data(),
                                ptBins, ptmin, ptmax, 100, 0., 2.);
      fMCQAPlots[i]->Add(fMCTPCRatio[i]);
      fMCTPCdedx[i] = new TH2F(MCTPCdedxName.Data(), MCTPCdedxName.Data(),
                               ptBins, ptmin, ptmax, 2 * twoDBins, 0., 400);
      fMCQAPlots[i]->Add(fMCTPCdedx[i]);
      fMCTOFbeta[i] = new TH2F(MCTOFbetaName.Data(), MCTOFbetaName.Data(),
                               ptBins, ptmin, ptmax, 3.5 * twoDBins, 0.4, 1.1);
      fMCQAPlots[i]->Add(fMCTOFbeta[i]);
      fMCNSigTPC[i] = new TH2F(MCNSigTPCName.Data(), MCNSigTPCName.Data(),
                               ptBins, ptmin, ptmax, twoDBins, -20., 20.);
      fMCQAPlots[i]->Add(fMCNSigTPC[i]);
      fMCNSigTOF[i] = new TH2F(MCNSigTOFName.Data(), MCNSigTOFName.Data(),
                               ptBins, ptmin, ptmax, twoDBins, -20., 20.);
      fMCQAPlots[i]->Add(fMCNSigTOF[i]);
    }
  } else {
    fMCContPt = 0;
    fMCUnknownPt = 0;
    fMCPrimaryPt = 0;
    fMCMaterialPt = 0;
    fMCFeeddownWeakPt = 0;
    for (int i = 0; i < 4; ++i) {
      fMCQAPlots[i] = 0;
      fMCpTPCDist[i] = 0;
      fMCetaDist[i] = 0;
      fMCphiDist[i] = 0;
      fMCTPCCls[i] = 0;
      fMCDCAxy[i] = 0;
      fMCDCAz[i] = 0;
      fMCTPCCrossedRows[i] = 0;
      fMCTPCRatio[i] = 0;
      fMCTPCdedx[i] = 0;
      fMCTOFbeta[i] = 0;
      fMCNSigTPC[i] = 0;
      fMCNSigTOF[i] = 0;
    }
  }

  if (DCADist) {
    fDCAPlots = new TList();
    fDCAPlots->SetName("DCAPtBinning");
    fDCAPlots->SetOwner();
    fMCList->Add(fDCAPlots);
    TString MCPridcaPtBinName = Form("DCAPtBinningPri");
    TString MCMatdcaPtBinName = Form("DCAPtBinningMat");
    TString MCSecdcaPtBinName = Form("DCAPtBinningSec");
    TString MCSecLamdcaPtBinName = Form("DCAPtBinningSecLam");
    TString MCSecSigdcaPtBinName = Form("DCAPtBinningSecSig");

    fMCPrimDCAXYPtBins = new TH2F(MCPridcaPtBinName.Data(),
                                  MCPridcaPtBinName.Data(), fpTbins, fpTmin,
                                  fpTmax, 500, -5, 5);
    fMCPrimDCAXYPtBins->GetXaxis()->SetTitle("dca_{XY}");
    fMCPrimDCAXYPtBins->GetYaxis()->SetTitle("dca_{Z}");
    fDCAPlots->Add(fMCPrimDCAXYPtBins);

    fMCMaterialDCAXYPtBins = new TH2F(MCMatdcaPtBinName.Data(),
                                      MCMatdcaPtBinName.Data(), fpTbins, fpTmin,
                                      fpTmax, 500, -5, 5);
    fMCMaterialDCAXYPtBins->GetXaxis()->SetTitle("dca_{XY}");
    fMCMaterialDCAXYPtBins->GetYaxis()->SetTitle("dca_{Z}");
    fDCAPlots->Add(fMCMaterialDCAXYPtBins);

    fMCSecondaryDCAXYPtBins = new TH2F(MCSecdcaPtBinName.Data(),
                                       MCSecdcaPtBinName.Data(), fpTbins,
                                       fpTmin, fpTmax, 500, -5, 5);
    fMCSecondaryDCAXYPtBins->GetXaxis()->SetTitle("dca_{XY}");
    fMCSecondaryDCAXYPtBins->GetYaxis()->SetTitle("dca_{Z}");
    fDCAPlots->Add(fMCSecondaryDCAXYPtBins);

    fMCSecLambdaDCAXYPtBins = new TH2F(MCSecLamdcaPtBinName.Data(),
                                       MCSecLamdcaPtBinName.Data(), fpTbins,
                                       fpTmin, fpTmax, 500, -5, 5);
    fMCSecLambdaDCAXYPtBins->GetXaxis()->SetTitle("dca_{XY}");
    fMCSecLambdaDCAXYPtBins->GetYaxis()->SetTitle("dca_{Z}");
    fDCAPlots->Add(fMCSecLambdaDCAXYPtBins);

    fMCSecSigmaDCAXYPtBins = new TH2F(MCSecSigdcaPtBinName.Data(),
                                      MCSecSigdcaPtBinName.Data(), fpTbins,
                                      fpTmin, fpTmax, 500, -5, 5);
    fMCSecSigmaDCAXYPtBins->GetXaxis()->SetTitle("dca_{XY}");
    fMCSecSigmaDCAXYPtBins->GetYaxis()->SetTitle("dca_{Z}");
    fDCAPlots->Add(fMCSecSigmaDCAXYPtBins);

    if (fDoMultiplicityBinning) {
      TString name1 = "DCAPtBinningPriMult_0_";
      name1 += fMultRangeLow;
      TString name2 = "DCAPtBinningPriMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      TString name3 = "DCAPtBinningPriMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString primPtBinName[3] = { name1, name2, name3 };

      name1 = "DCAPtBinningMatMult_0_";
      name1 += fMultRangeLow;
      name2 = "DCAPtBinningMatMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "DCAPtBinningMatMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString matPtBinName[3] = { name1, name2, name3 };

      name1 = "DCAPtBinningSecMult_0_";
      name1 += fMultRangeLow;
      name2 = "DCAPtBinningSecMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "DCAPtBinningSecMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString secPtBinName[3] = { name1, name2, name3 };

      name1 = "DCAPtBinningSecLambdaMult_0_";
      name1 += fMultRangeLow;
      name2 = "DCAPtBinningSecLambdaMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "DCAPtBinningSecLambdaMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString secLambdaPtBinName[3] = { name1, name2, name3 };

      name1 = "DCAPtBinningSecSigmaMult_0_";
      name1 += fMultRangeLow;
      name2 = "DCAPtBinningSecSigmaMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "DCAPtBinningSecSigmaMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString secSigmaPtBinName[3] = { name1, name2, name3 };

      name1 = "0 < mult < ";
      name1 += fMultRangeLow;
      name1 += ";P#_{T};dca_{XY}";
      name2 = "";
      name2 += fMultRangeLow;
      name2 += " < mult < ";
      name2 += fMultRangeHigh;
      name2 += ";P#_{T};dca_{XY}";
      name3 = "mult > ";
      name3 += fMultRangeHigh;
      name3 += ";P#_{T};dca_{XY}";
      TString axisRange[3] = { name1, name2, name3 };

      for (int i = 0; i < 3; ++i) {
        fMCPrimDCAXYPtBinsMult[i] = new TH2F(primPtBinName[i].Data(),
                                             axisRange[i].Data(), fpTbins,
                                             fpTmin, fpTmax, 500, -5, 5);
        fDCAPlots->Add(fMCPrimDCAXYPtBinsMult[i]);

        fMCMaterialDCAXYPtBinsMult[i] = new TH2F(matPtBinName[i].Data(),
                                                 axisRange[i].Data(), fpTbins,
                                                 fpTmin, fpTmax, 500, -5, 5);
        fDCAPlots->Add(fMCMaterialDCAXYPtBinsMult[i]);

        fMCSecondaryDCAXYPtBinsMult[i] = new TH2F(secPtBinName[i].Data(),
                                                  axisRange[i].Data(), fpTbins,
                                                  fpTmin, fpTmax, 500, -5, 5);
        fDCAPlots->Add(fMCSecondaryDCAXYPtBinsMult[i]);

        fMCSecLambdaDCAXYPtBinsMult[i] = new TH2F(secLambdaPtBinName[i].Data(),
                                                  axisRange[i].Data(), fpTbins,
                                                  fpTmin, fpTmax, 500, -5, 5);
        fDCAPlots->Add(fMCSecLambdaDCAXYPtBinsMult[i]);

        fMCSecSigmaDCAXYPtBinsMult[i] = new TH2F(secSigmaPtBinName[i].Data(),
                                                 axisRange[i].Data(), fpTbins,
                                                 fpTmin, fpTmax, 500, -5, 5);
        fDCAPlots->Add(fMCSecSigmaDCAXYPtBinsMult[i]);
      }
    }
  } else {
    fDCAPlots = 0;
    fMCPrimDCAXYPtBins = 0;
    fMCMaterialDCAXYPtBins = 0;
    fMCSecondaryDCAXYPtBins = 0;
    fMCSecLambdaDCAXYPtBins = 0;
    fMCSecSigmaDCAXYPtBins = 0;
  }
}
void AliFemtoDreamTrackMCHist::FillMCDCAXYPtBins(
    AliFemtoDreamBasePart::PartOrigin org, int PDGCodeMoth, float pT,
    float dcaxy, int multiplicity) {
  if (!fDoDCAPlots) {
    AliFatal("FullBooking not set for SPCutHistograms! Cannot use this method");
  }
  if (org == AliFemtoDreamBasePart::kPhysPrimary) {
    fMCPrimDCAXYPtBins->Fill(pT, dcaxy);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, dcaxy, fMCPrimDCAXYPtBinsMult[0],
                             fMCPrimDCAXYPtBinsMult[1],
                             fMCPrimDCAXYPtBinsMult[2]);
    }
  } else if (org == AliFemtoDreamBasePart::kWeak) {
    fMCSecondaryDCAXYPtBins->Fill(pT, dcaxy);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, dcaxy,
                             fMCSecondaryDCAXYPtBinsMult[0],
                             fMCSecondaryDCAXYPtBinsMult[1],
                             fMCSecondaryDCAXYPtBinsMult[2]);
    }
    if (TMath::Abs(PDGCodeMoth) == 3222) {
      fMCSecSigmaDCAXYPtBins->Fill(pT, dcaxy);
      if (fDoMultiplicityBinning) {
        FillMultiplicityHistos(multiplicity, pT, dcaxy,
                               fMCSecSigmaDCAXYPtBinsMult[0],
                               fMCSecSigmaDCAXYPtBinsMult[1],
                               fMCSecSigmaDCAXYPtBinsMult[2]);
      }
    } else if (TMath::Abs(PDGCodeMoth) == 3122) {
      fMCSecLambdaDCAXYPtBins->Fill(pT, dcaxy);
      if (fDoMultiplicityBinning) {
        FillMultiplicityHistos(multiplicity, pT, dcaxy,
                               fMCSecLambdaDCAXYPtBinsMult[0],
                               fMCSecLambdaDCAXYPtBinsMult[1],
                               fMCSecLambdaDCAXYPtBinsMult[2]);
      }
    } else {
      TString ErrHistSP = Form("Feeddown for %d not implemented", PDGCodeMoth);
      AliWarning(ErrHistSP.Data());
    }
  } else if (org == AliFemtoDreamBasePart::kMaterial) {
    fMCMaterialDCAXYPtBins->Fill(pT, dcaxy);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, dcaxy,
                             fMCMaterialDCAXYPtBinsMult[0],
                             fMCMaterialDCAXYPtBinsMult[1],
                             fMCMaterialDCAXYPtBinsMult[2]);
    }
  } else {
    AliFatal("Particle Origin not implemented");
  }
  return;
}

void AliFemtoDreamTrackMCHist::FillMultiplicityHistos(int multiplicity,
                                                      float pT, float dcaxy,
                                                      TH2F *histo1,
                                                      TH2F *histo2,
                                                      TH2F *histo3) {
  if (multiplicity < fMultRangeLow) {
    histo1->Fill(pT, dcaxy);
  } else if (multiplicity >= fMultRangeLow && multiplicity < fMultRangeHigh) {
    histo2->Fill(pT, dcaxy);
  } else {
    histo3->Fill(pT, dcaxy);
  }
}

void AliFemtoDreamTrackMCHist::FillMCPtResolution(float pTTrue, float pTReco) {
  fPtResolution->Fill(pTTrue, (pTTrue - pTReco) / pTTrue);
}

void AliFemtoDreamTrackMCHist::FillMCThetaResolution(float ThetaTrue,
                                                     float ThetaReco,
                                                     float pTTrue) {
  fThetaResolution->Fill(pTTrue, ThetaTrue - ThetaReco);
}

void AliFemtoDreamTrackMCHist::FillMCPhiResolution(float PhiTrue, float PhiReco,
                                                   float pTTrue) {
  fPhiResolution->Fill(pTTrue, PhiTrue - PhiReco);
}
