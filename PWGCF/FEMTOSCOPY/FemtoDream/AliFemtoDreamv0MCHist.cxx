/*
 * AliFemtoDreamv0MCHist.cxx
 *
 *  Created on: Dec 12,2017
 *      Author: gu74req
 */

#include "AliFemtoDreamv0MCHist.h"
#include "AliFemtoDreamBasePart.h"
ClassImp(AliFemtoDreamv0MCHist)
AliFemtoDreamv0MCHist::AliFemtoDreamv0MCHist()
    : fMCList(0),
      fCPAPlots(nullptr),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fDoMultiplicityBinning(false),
      fMCCorrPt(nullptr),
      fMCIdentPt(nullptr),
      fMCGenPt(nullptr),
      fMCContPt(nullptr),
      fMCPrimaryPt(nullptr),
      fMCMaterialPt(nullptr),
      fMCFeeddownWeakPt(nullptr),
      fHistMCMother(nullptr),
      fHistMCMotherPDG(nullptr),
      fMCPrimCPAPtBins(nullptr),
      fMCMaterialCPAPtBins(nullptr),
      fMCSecondaryCPAPtBins(nullptr),
      fMCContCPAPtBins(nullptr),
      fPtResolution(nullptr),
      fThetaResolution(nullptr),
      fPhiResolution(nullptr) {
  for (int i = 0; i < 5; ++i) {
    if (i < 3) {
      fMCPrimCPAPtBinsMult[i] = 0;
      fMCMaterialCPAPtBinsMult[i] = 0;
      fMCSecondaryCPAPtBinsMult[i] = 0;
      fMCContCPAPtBinsMult[i] = 0;
    }
    fMCQAPlots[i] = 0;
    fMCpTDist[i] = 0;
    fMCetaDist[i] = 0;
    fMCphiDist[i] = 0;
    fMCDecayVtxv0X[i] = 0;
    fMCDecayVtxv0Y[i] = 0;
    fMCDecayVtxv0Z[i] = 0;
    fMCTransverseRadius[i] = 0;
    fMCDCAPosDaugToPV[i] = 0;
    fMCDCANegDaugToPV[i] = 0;
    fMCDCADaugToVtx[i] = 0;
    fMCCosPointing[i] = 0;
    fMCInvMass[i] = 0;
    fMCBachDCAToPV[i] = 0;
    fMCv0DecayLength[i] = 0;
    fMCDecayLength[i] = 0;
    fMCXiRapidity[i] = 0;
    fMCOmegaRapidity[i] = 0;
    fMCPodolanski[i] = 0;
    fMCXiInvMass[i] = 0;
    fMCOmegaInvMass[i] = 0;
  }
}

AliFemtoDreamv0MCHist::AliFemtoDreamv0MCHist(int MassNBins, float MassMin,
                                             float MassMax,
                                             bool contribSplitting,
                                             bool CPADist, bool DoMultBinning,
                                             bool checkMother)
    : fMCList(nullptr),
      fCPAPlots(nullptr),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fDoMultiplicityBinning(DoMultBinning),
      fMCCorrPt(nullptr),
      fMCIdentPt(nullptr),
      fMCGenPt(nullptr),
      fMCContPt(nullptr),
      fMCPrimaryPt(nullptr),
      fMCMaterialPt(nullptr),
      fMCFeeddownWeakPt(nullptr),
      fHistMCMother(nullptr),
      fHistMCMotherPDG(nullptr),
      fMCPrimCPAPtBins(nullptr),
      fMCMaterialCPAPtBins(nullptr),
      fMCSecondaryCPAPtBins(nullptr),
      fMCContCPAPtBins(nullptr),
      fPtResolution(nullptr),
      fThetaResolution(nullptr),
      fPhiResolution(nullptr) {
  float ptmin = -0.2;
  float ptmax = 6.3;
  int ptBins = 13;

  fMCList = new TList();
  fMCList->SetName("v0MonteCarlo");
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
                           100, 0, 5, 500, -1, 1);
  fPtResolution->GetXaxis()->SetTitle("P_{T,True}");
  fPtResolution->GetYaxis()->SetTitle("(P_{T,True}-P_{T,Reco})/P_{T,True}");
  fMCList->Add(fPtResolution);

  fThetaResolution = new TH2F("DeltaThetaRecoTruevsPtReco",
                              "DeltaThetaRecoTruevsPtReco", 100, 0, 5, 500,
                              -0.3, 0.3);
  fThetaResolution->GetXaxis()->SetTitle("P_{T,True}");
  fThetaResolution->GetYaxis()->SetTitle("#Theta_{T,True}-#Theta_{T,Reco}");
  fMCList->Add(fThetaResolution);

  fPhiResolution = new TH2F("DeltaPhiRecoTruevsPtReco",
                            "DeltaPhiRecoTruevsPtReco", 100, 0, 5, 500, -0.3,
                            0.3);
  fPhiResolution->GetXaxis()->SetTitle("P_{T,True}");
  fPhiResolution->GetYaxis()->SetTitle("#Phi_{T,True}-#Phi_{T,Reco}");
  fMCList->Add(fPhiResolution);

  if (checkMother) {
    fHistMCMother = new TH2F("fHistMCMother",
                             "; #it{p}_{T} (GeV/#it{c}^{2}); PDG code mother",
                             100, 0., 10., 4000, 0, 4000);
    fMCList->Add(fHistMCMother);

    fHistMCMotherPDG = new TH1I("fHistMCMotherPDG", ";Entries; PDG code mother",
                                10000000, 0, 10000000);
    fMCList->Add(fHistMCMotherPDG);
  }

  if (contribSplitting) {

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

    fMCContPt = new TH1F("ContPt", "ContPt", ptBins, ptmin, ptmax);
    fMCContPt->GetXaxis()->SetTitle("p_{T}");
    fMCList->Add(fMCContPt);

    TString MCModes[5] = { "Primary", "Secondary", "Material", "Contamination",
        "Fake" };
    for (int i = 0; i < 5; ++i) {
      fMCQAPlots[i] = new TList();
      fMCQAPlots[i]->SetOwner();
      fMCQAPlots[i]->SetName(MCModes[i].Data());
      fMCList->Add(fMCQAPlots[i]);

      TString MCPtDist = Form("MCPt%s", MCModes[i].Data());
      fMCpTDist[i] = new TH1F(MCPtDist.Data(), MCPtDist.Data(), ptBins, ptmin,
                              ptmax);
      fMCQAPlots[i]->Add(fMCpTDist[i]);

      TString MCEtaDist = Form("MCEta%s", MCModes[i].Data());
      fMCetaDist[i] = new TH1F(MCEtaDist.Data(), MCEtaDist.Data(), 200, -10.,
                               10.);
      fMCQAPlots[i]->Add(fMCetaDist[i]);

      TString MCPhiDist = Form("MCPhi%s", MCModes[i].Data());
      fMCphiDist[i] = new TH1F(MCPhiDist.Data(), MCPhiDist.Data(), 100, 0.,
                               2 * TMath::Pi());
      fMCQAPlots[i]->Add(fMCphiDist[i]);

      TString MCDecayVtxv0XDist = Form("MCDecayVtxXPV%s", MCModes[i].Data());
      fMCDecayVtxv0X[i] = new TH2F(MCDecayVtxv0XDist.Data(),
                                   MCDecayVtxv0XDist.Data(), ptBins, ptmin,
                                   ptmax, 400, 0, 200);
      fMCQAPlots[i]->Add(fMCDecayVtxv0X[i]);

      TString MCDecayVtxv0YDist = Form("MCDecayVtxYPV%s", MCModes[i].Data());
      fMCDecayVtxv0Y[i] = new TH2F(MCDecayVtxv0YDist.Data(),
                                   MCDecayVtxv0YDist.Data(), ptBins, ptmin,
                                   ptmax, 400, 0, 200);
      fMCQAPlots[i]->Add(fMCDecayVtxv0Y[i]);

      TString MCDecayVtxv0ZDist = Form("MCDecayVtxZPV%s", MCModes[i].Data());
      fMCDecayVtxv0Z[i] = new TH2F(MCDecayVtxv0ZDist.Data(),
                                   MCDecayVtxv0ZDist.Data(), ptBins, ptmin,
                                   ptmax, 400, 0, 200);
      fMCQAPlots[i]->Add(fMCDecayVtxv0Z[i]);

      TString MCTransverseRadius = Form("MCTransverseRadius%s",
                                        MCModes[i].Data());
      fMCTransverseRadius[i] = new TH2F(MCTransverseRadius.Data(),
                                        MCTransverseRadius.Data(), ptBins,
                                        ptmin, ptmax, 750, 0, 150);
      fMCQAPlots[i]->Add(fMCTransverseRadius[i]);

      TString MCDCADaugPVP = Form("MCDCADauPToPV%s", MCModes[i].Data());
      fMCDCAPosDaugToPV[i] = new TH2F(MCDCADaugPVP.Data(), MCDCADaugPVP.Data(),
                                      ptBins, ptmin, ptmax, 500, 0, 100);
      fMCQAPlots[i]->Add(fMCDCAPosDaugToPV[i]);

      TString MCDCADaugPVN = Form("MCDCADauNToPV%s", MCModes[i].Data());
      fMCDCANegDaugToPV[i] = new TH2F(MCDCADaugPVN.Data(), MCDCADaugPVN.Data(),
                                      ptBins, ptmin, ptmax, 500, 0, 100);
      fMCQAPlots[i]->Add(fMCDCANegDaugToPV[i]);

      TString MCDCADaugToVTX = Form("MCDCADauToVtx%s", MCModes[i].Data());
      fMCDCADaugToVtx[i] = new TH2F(MCDCADaugToVTX.Data(),
                                    MCDCADaugToVTX.Data(), ptBins, ptmin, ptmax,
                                    100, 0, 10);
      fMCQAPlots[i]->Add(fMCDCADaugToVtx[i]);

      TString MCCosPointing = Form("MCPointingAngle%s", MCModes[i].Data());
      fMCCosPointing[i] = new TH2F(MCCosPointing.Data(), MCCosPointing.Data(),
                                   ptBins, ptmin, ptmax, 400, 0.85, 1.001);
      ;
      fMCQAPlots[i]->Add(fMCCosPointing[i]);

      TString MCInvMass = Form("MCInvMass%s", MCModes[i].Data());
      fMCInvMass[i] = new TH1F(MCInvMass.Data(), MCInvMass.Data(), MassNBins,
                               MassMin, MassMax);
      fMCQAPlots[i]->Add(fMCInvMass[i]);

      TString MCBachDCAPV = Form("MCBachDCAPV%s", MCModes[i].Data());
      fMCBachDCAToPV[i] = new TH2F(MCBachDCAPV.Data(), MCBachDCAPV.Data(),
                                   ptBins, ptmin, ptmax, 50, 0, 10);
      fMCQAPlots[i]->Add(fMCBachDCAToPV[i]);

      TString MCv0DecayLength = Form("MCv0DecayLength%s", MCModes[i].Data());
      fMCv0DecayLength[i] = new TH2F(MCv0DecayLength.Data(),
                                     MCv0DecayLength.Data(), ptBins, ptmin,
                                     ptmax, 200, 0, 10);
      fMCQAPlots[i]->Add(fMCv0DecayLength[i]);

      TString MCv0CPA = Form("MCv0CPA%s", MCModes[i].Data());
      fMCv0CPA[i] = new TH2F(MCv0CPA.Data(), MCv0CPA.Data(), ptBins, ptmin,
                             ptmax, 400, 0.85, 1.001);
      fMCQAPlots[i]->Add(fMCv0CPA[i]);

      TString MCXiDecayLength = Form("MCDecayLength%s", MCModes[i].Data());
      fMCDecayLength[i] = new TH2F(MCXiDecayLength.Data(),
                                   MCXiDecayLength.Data(), ptBins, ptmin, ptmax,
                                   200, 0, 10);
      fMCQAPlots[i]->Add(fMCDecayLength[i]);

      TString MCXiRapidity = Form("MCXiRapidity%s", MCModes[i].Data());
      fMCXiRapidity[i] = new TH2F(MCXiRapidity.Data(), MCXiRapidity.Data(), 50,
                                  -2, 2, ptBins, ptmin, ptmax);
      fMCQAPlots[i]->Add(fMCXiRapidity[i]);

      TString MCOmegaRapidity = Form("MCOmegaRapidity%s", MCModes[i].Data());
      fMCOmegaRapidity[i] = new TH2F(MCOmegaRapidity.Data(),
                                     MCOmegaRapidity.Data(), 50, -2, 2, ptBins,
                                     ptmin, ptmax);
      fMCQAPlots[i]->Add(fMCOmegaRapidity[i]);

      TString MCPodolanski = Form("MCPodolanski%s", MCModes[i].Data());
      fMCPodolanski[i] = new TH2F(MCPodolanski.Data(), MCPodolanski.Data(), 50,
                                  -1, 1, 50, 0, 1);
      fMCQAPlots[i]->Add(fMCPodolanski[i]);

      TString MCXiInvMass = Form("MCXiInvMass%s", MCModes[i].Data());
      fMCXiInvMass[i] = new TH2F(MCXiInvMass.Data(), MCXiInvMass.Data(), ptBins,
                                 ptmin, ptmax, 500, 1.31486 * 0.9,
                                 1.31486 * 1.3);
      fMCQAPlots[i]->Add(fMCXiInvMass[i]);

      TString MCOmegaInvMass = Form("MCOmegaInvMass%s", MCModes[i].Data());
      fMCOmegaInvMass[i] = new TH2F(MCOmegaInvMass.Data(),
                                    MCOmegaInvMass.Data(), ptBins, ptmin, ptmax,
                                    500, 1.672 * 0.9, 1.672 * 1.3);
      fMCQAPlots[i]->Add(fMCOmegaInvMass[i]);
    }
  } else {
    fMCPrimaryPt = 0;
    fMCMaterialPt = 0;
    fMCFeeddownWeakPt = 0;
    fMCContPt = 0;
    for (int i = 0; i < 4; ++i) {
      fMCQAPlots[i] = 0;
      fMCpTDist[i] = 0;
      fMCetaDist[i] = 0;
      fMCphiDist[i] = 0;
      fMCDecayVtxv0X[i] = 0;
      fMCDecayVtxv0Y[i] = 0;
      fMCDecayVtxv0Z[i] = 0;
      fMCTransverseRadius[i] = 0;
      fMCDCAPosDaugToPV[i] = 0;
      fMCDCANegDaugToPV[i] = 0;
      fMCDCADaugToVtx[i] = 0;
      fMCCosPointing[i] = 0;
      fMCInvMass[i] = 0;
      fMCBachDCAToPV[i] = 0;
      fMCv0DecayLength[i] = 0;
      fMCDecayLength[i] = 0;
      fMCXiRapidity[i] = 0;
      fMCOmegaRapidity[i] = 0;
      fMCPodolanski[i] = 0;
      fMCXiInvMass[i] = 0;
      fMCOmegaInvMass[i] = 0;
    }
  }

  if (CPADist) {
    fCPAPlots = new TList();
    fCPAPlots->SetName("CPAPtBinning");
    fCPAPlots->SetOwner();
    fMCList->Add(fCPAPlots);
    TString MCPricpaPtBinName = Form("CPAPtBinningPri");
    TString MCMatcpaPtBinName = Form("CPAPtBinningMat");
    TString MCSeccpaPtBinName = Form("CPAPtBinningSec");
    TString MCConcpaPtBinName = Form("CPAPtBinningCont");

    fMCPrimCPAPtBins = new TH2F(MCPricpaPtBinName.Data(),
                                MCPricpaPtBinName.Data(), 8, 0.3, 4.3, 1000,
                                0.90, 1.);
    fMCPrimCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCPrimCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCPrimCPAPtBins);

    fMCMaterialCPAPtBins = new TH2F(MCMatcpaPtBinName.Data(),
                                    MCMatcpaPtBinName.Data(), 8, 0.3, 4.3, 1000,
                                    0.90, 1.);
    fMCMaterialCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCMaterialCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCMaterialCPAPtBins);

    fMCSecondaryCPAPtBins = new TH2F(MCSeccpaPtBinName.Data(),
                                     MCSeccpaPtBinName.Data(), 8, 0.3, 4.3,
                                     1000, 0.90, 1.);
    fMCSecondaryCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCSecondaryCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCSecondaryCPAPtBins);

    fMCContCPAPtBins = new TH2F(MCConcpaPtBinName.Data(),
                                MCConcpaPtBinName.Data(), 8, 0.3, 4.3, 1000,
                                0.90, 1.);
    fMCContCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCContCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCContCPAPtBins);

    if (fDoMultiplicityBinning) {
      TString name1 = "CPAPtBinningPriMult_0_";
      name1 += fMultRangeLow;
      TString name2 = "CPAPtBinningPriMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      TString name3 = "CPAPtBinningPriMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString primPtBinName[3] = { name1, name2, name3 };

      name1 = "CPAPtBinningMatMult_0_";
      name1 += fMultRangeLow;
      name2 = "CPAPtBinningMatMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "CPAPtBinningMatMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString matPtBinName[3] = { name1, name2, name3 };

      name1 = "CPAPtBinningSecMult_0_";
      name1 += fMultRangeLow;
      name2 = "CPAPtBinningSecMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "CPAPtBinningSecMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString secPtBinName[3] = { name1, name2, name3 };

      name1 = "CPAPtBinningContMult_0_";
      name1 += fMultRangeLow;
      name2 = "CPAPtBinningContMult_";
      name2 += fMultRangeLow;
      name2 += "_";
      name2 += fMultRangeHigh;
      name3 = "CPAPtBinningContMult_";
      name3 += fMultRangeHigh;
      name3 += "_inf";
      TString contPtBinName[3] = { name1, name2, name3 };

      name1 = "0 < mult < ";
      name1 += fMultRangeLow;
      name1 += ";P#_{T};CPA";
      name2 = "";
      name2 += fMultRangeLow;
      name2 += " < mult < ";
      name2 += fMultRangeHigh;
      name2 += ";P#_{T};CPA";
      name3 = "mult > ";
      name3 += fMultRangeHigh;
      name3 += ";P#_{T};CPA";
      TString axisRange[3] = { name1, name2, name3 };

      for (int i = 0; i < 3; ++i) {
        fMCPrimCPAPtBinsMult[i] = new TH2F(primPtBinName[i].Data(),
                                           axisRange[i].Data(), 8, 0.3, 4.3,
                                           1000, 0.90, 1.);
        fCPAPlots->Add(fMCPrimCPAPtBinsMult[i]);

        fMCMaterialCPAPtBinsMult[i] = new TH2F(matPtBinName[i].Data(),
                                               axisRange[i].Data(), 8, 0.3, 4.3,
                                               1000, 0.90, 1.);
        fCPAPlots->Add(fMCMaterialCPAPtBinsMult[i]);

        fMCSecondaryCPAPtBinsMult[i] = new TH2F(secPtBinName[i].Data(),
                                                axisRange[i].Data(), 8, 0.3,
                                                4.3, 1000, 0.90, 1.);
        fCPAPlots->Add(fMCSecondaryCPAPtBinsMult[i]);

        fMCContCPAPtBinsMult[i] = new TH2F(contPtBinName[i].Data(),
                                           axisRange[i].Data(), 8, 0.3, 4.3,
                                           1000, 0.90, 1.);
        fCPAPlots->Add(fMCContCPAPtBinsMult[i]);
      }
    }
  } else {
    fCPAPlots = 0;
    fMCPrimCPAPtBins = 0;
    fMCMaterialCPAPtBins = 0;
    fMCSecondaryCPAPtBins = 0;
    fMCContCPAPtBins = 0;
  }
}

AliFemtoDreamv0MCHist::~AliFemtoDreamv0MCHist() {
  if (fMCList) {
    delete fMCList;
  }
}

void AliFemtoDreamv0MCHist::FillMCCPAPtBins(
    AliFemtoDreamBasePart::PartOrigin org, float pT, float cpa,
    int multiplicity) {
  if (org == AliFemtoDreamBasePart::kPhysPrimary) {
    fMCPrimCPAPtBins->Fill(pT, cpa);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, cpa, fMCPrimCPAPtBinsMult[0],
                             fMCPrimCPAPtBinsMult[1], fMCPrimCPAPtBinsMult[2]);
    }
  } else if (org == AliFemtoDreamBasePart::kWeak) {
    fMCSecondaryCPAPtBins->Fill(pT, cpa);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, cpa,
                             fMCSecondaryCPAPtBinsMult[0],
                             fMCSecondaryCPAPtBinsMult[1],
                             fMCSecondaryCPAPtBinsMult[2]);
    }
  } else if (org == AliFemtoDreamBasePart::kMaterial) {
    fMCMaterialCPAPtBins->Fill(pT, cpa);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, cpa, fMCMaterialCPAPtBinsMult[0],
                             fMCMaterialCPAPtBinsMult[1],
                             fMCMaterialCPAPtBinsMult[2]);
    }
  } else if (org == AliFemtoDreamBasePart::kFake) {
    fMCContCPAPtBins->Fill(pT, cpa);
    if (fDoMultiplicityBinning) {
      FillMultiplicityHistos(multiplicity, pT, cpa, fMCContCPAPtBinsMult[0],
                             fMCContCPAPtBinsMult[1], fMCContCPAPtBinsMult[2]);
    }
  } else {
    AliFatal("Particle Origin not implemented");
  }
}

void AliFemtoDreamv0MCHist::FillMultiplicityHistos(int multiplicity, float pT,
                                                   float cpa, TH2F *histo1,
                                                   TH2F *histo2, TH2F *histo3) {
  if (multiplicity < fMultRangeLow) {
    histo1->Fill(pT, cpa);
  } else if (multiplicity >= fMultRangeLow && multiplicity < fMultRangeHigh) {
    histo2->Fill(pT, cpa);
  } else {
    histo3->Fill(pT, cpa);
  }
}

void AliFemtoDreamv0MCHist::FillMCPtResolution(float pTTrue, float pTReco) {
  fPtResolution->Fill(pTTrue, (pTTrue - pTReco) / pTTrue);
}

void AliFemtoDreamv0MCHist::FillMCThetaResolution(float ThetaTrue,
                                                  float ThetaReco,
                                                  float pTTrue) {
  fThetaResolution->Fill(pTTrue, ThetaTrue - ThetaReco);
}

void AliFemtoDreamv0MCHist::FillMCPhiResolution(float PhiTrue, float PhiReco,
                                                float pTTrue) {
  fPhiResolution->Fill(pTTrue, PhiTrue - PhiReco);
}
