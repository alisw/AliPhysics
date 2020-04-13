/*
 * AliFemtoDreamv0Hist.cxx
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamv0Hist.h"
#include "TMath.h"
#include <iostream>
ClassImp(AliFemtoDreamv0Hist)

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist()
    : fMinimalBooking(false),
      fHistList(0),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fConfig(0),
      fCutCounter(0),
      fInvMassBefKaonRej(0),
      fInvMassKaon(0),
      fInvMassBefSelection(0),
      fInvMassPt(0),
      fCPAPtBins(0),
      fCPAPtBinsMult(),
      fInvMassPerRunNumber() {
  for (int i = 0; i < 2; ++i) {
    fv0CutQA[i] = nullptr;
    fOnFly[i] = nullptr;
    fpTDist[i] = nullptr;
    fetaDist[i] = nullptr;
    fDecayVtxv0X[i] = nullptr;
    fDecayVtxv0Y[i] = nullptr;
    fDecayVtxv0Z[i] = nullptr;
    fTransRadius[i] = nullptr;
    fDCAPosDaugToPrimVtx[i] = nullptr;
    fDCANegDaugToPrimVtx[i] = nullptr;
    fDCADaugToVtx[i] = nullptr;
    fCPA[i] = nullptr;
    fInvMass[i] = nullptr;
  }
}

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist(int MassNBins, float MassMin,
                                         float MassMax, bool CPAPlots,
                                         bool perRunnumber, int iRunMin = 0,
                                         int iRunMax = 0)
    : fMinimalBooking(false),
      fMultRangeLow(27),
      fMultRangeHigh(55) {
  TString sName[2] = { "before", "after" };

  fHistList = new TList();
  fHistList->SetName("v0Cuts");
  fHistList->SetOwner();

  fConfig = new TProfile("TrackCutConfig", "Track Cut Config", 20, 0, 20);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1, "OnFly Status");
  fConfig->GetXaxis()->SetBinLabel(2, "Charge");
  fConfig->GetXaxis()->SetBinLabel(3, "p_{T, Min}");
  fConfig->GetXaxis()->SetBinLabel(4, "p_{T, Max}");
  fConfig->GetXaxis()->SetBinLabel(5, "K0 Rejection");
  fConfig->GetXaxis()->SetBinLabel(6, "Max Decay Vertex XYZ");
  fConfig->GetXaxis()->SetBinLabel(7, "Min Transverse Radius");
  fConfig->GetXaxis()->SetBinLabel(8, "Max Transverse Radius");
  fConfig->GetXaxis()->SetBinLabel(9, "Min DCA to PV of Daug");
  fConfig->GetXaxis()->SetBinLabel(10, "Max Distance Daug to Vertex");
  fConfig->GetXaxis()->SetBinLabel(11, "Inv Mass Cut down");
  fConfig->GetXaxis()->SetBinLabel(12, "Inv Mass Cut up");
  fConfig->GetXaxis()->SetBinLabel(13, "Min Cos Pointing Angle");
  fConfig->GetXaxis()->SetBinLabel(14, "Armenteros q_{T} low");
  fConfig->GetXaxis()->SetBinLabel(15, "Armenteros q_{T} up");
  fConfig->GetXaxis()->SetBinLabel(16, "Armenteros #alpha low");
  fConfig->GetXaxis()->SetBinLabel(17, "Armenteros #alpha up");
  fConfig->GetXaxis()->SetBinLabel(18, "Pileup requirement");

  fHistList->Add(fConfig);

  fCutCounter = new TH1F("CutCounter", "Cut Counter", 20, 0, 20);
  fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
  fCutCounter->GetXaxis()->SetBinLabel(2, "Has Daug");
  fCutCounter->GetXaxis()->SetBinLabel(3, "Daugters pass Cuts");
  fCutCounter->GetXaxis()->SetBinLabel(4, "On Fly Status");
  fCutCounter->GetXaxis()->SetBinLabel(5, "Charge");
  fCutCounter->GetXaxis()->SetBinLabel(6, "p_{T}");
  fCutCounter->GetXaxis()->SetBinLabel(7, "v0 DCA PV");
  fCutCounter->GetXaxis()->SetBinLabel(8, "Transverse Radius");
  fCutCounter->GetXaxis()->SetBinLabel(9, "Daug PV");
  fCutCounter->GetXaxis()->SetBinLabel(10, "Daug Vtx");
  fCutCounter->GetXaxis()->SetBinLabel(11, "Pileup");
  fCutCounter->GetXaxis()->SetBinLabel(12, "Armenteros-Podolandski");
  fCutCounter->GetXaxis()->SetBinLabel(13, "K0 Rejection");
  fCutCounter->GetXaxis()->SetBinLabel(14, "Inv Mass Cut");
  fCutCounter->GetXaxis()->SetBinLabel(15, "Cos Pointing Angle");
  fCutCounter->GetXaxis()->SetBinLabel(16, "D1&D2 right");
  fCutCounter->GetXaxis()->SetBinLabel(17, "D1&D2 pass cuts");
  fCutCounter->GetXaxis()->SetBinLabel(18, "D1&D2 wrong");
  fCutCounter->GetXaxis()->SetBinLabel(19, "D1&D2 pass cuts");
  fCutCounter->GetYaxis()->SetTitle("Entries");

  fHistList->Add(fCutCounter);

  fInvMassBefKaonRej = new TH1F("InvMassBefK0Rej", "InvMassBefK0Rej", MassNBins,
                                MassMin, MassMax);
  fInvMassBefKaonRej->GetXaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
  fInvMassBefKaonRej->GetYaxis()->SetTitle("Entries");
  fHistList->Add(fInvMassBefKaonRej);

  fInvMassKaon = new TH1F("InvMassKaon", "InvMassKaon", 400, 0.4, 0.6);
  fInvMassKaon->GetXaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
  fInvMassKaon->GetYaxis()->SetTitle("Entries");
  fHistList->Add(fInvMassKaon);

  fInvMassBefSelection = new TH1F("InvMasswithCuts", "InvMasswithCuts",
                                  MassNBins, MassMin, MassMax);
  fInvMassBefSelection->GetXaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
  fInvMassBefSelection->GetYaxis()->SetTitle("Entries");
  fHistList->Add(fInvMassBefSelection);

  fInvMassPt = new TH2F("InvMassPt", "Invariant Mass in Pt Bins", 8, 0.3, 4.3,
                        MassNBins, MassMin, MassMax);
  fInvMassPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fInvMassPt->GetYaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
  fHistList->Add(fInvMassPt);

  for (int i = 0; i < 2; ++i) {
    fv0CutQA[i] = new TList();
    fv0CutQA[i]->SetName(sName[i].Data());
    fv0CutQA[i]->SetOwner();
    fHistList->Add(fv0CutQA[i]);

    TString OnFlyName = Form("OnFly_%s", sName[i].Data());
    fOnFly[i] = new TH1F(OnFlyName.Data(), OnFlyName.Data(), 2, 0, 2.);
    fOnFly[i]->GetXaxis()->SetBinLabel(1, "Online");
    fOnFly[i]->GetXaxis()->SetBinLabel(2, "Offline");
    fOnFly[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fOnFly[i]);

    TString ptname = Form("pTDist_%s", sName[i].Data());
    fpTDist[i] = new TH1F(ptname.Data(), ptname.Data(), 200, 0, 10.);
    fpTDist[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fpTDist[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fpTDist[i]);

    TString etaname = Form("EtaDist_%s", sName[i].Data());
    fetaDist[i] = new TH1F(etaname.Data(), etaname.Data(), 200, -2., 2.);
    fetaDist[i]->GetXaxis()->SetTitle("#eta");
    fetaDist[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fetaDist[i]);

    TString phiname = Form("PhiDist_%s", sName[i].Data());
    fPhiDist[i] = new TH1F(phiname.Data(), phiname.Data(), 100, 0.,
                           2 * TMath::Pi());
    fPhiDist[i]->GetXaxis()->SetTitle("#phi");
    fPhiDist[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fPhiDist[i]);

    TString decayVtxXname = Form("DecayVtxXPV_%s", sName[i].Data());
    fDecayVtxv0X[i] = new TH1F(decayVtxXname.Data(), decayVtxXname.Data(), 400,
                               0., 200);
    fDecayVtxv0X[i]->GetXaxis()->SetTitle("Decay Vtx To PV X (cm)");
    fDecayVtxv0X[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDecayVtxv0X[i]);

    TString decayVtxYname = Form("DecayVtxYPV_%s", sName[i].Data());
    fDecayVtxv0Y[i] = new TH1F(decayVtxYname.Data(), decayVtxYname.Data(), 400,
                               0., 200);
    fDecayVtxv0Y[i]->GetXaxis()->SetTitle("Decay Vtx To PV Y (cm)");
    fDecayVtxv0Y[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDecayVtxv0Y[i]);

    TString decayVtxZname = Form("DecayVtxZPV_%s", sName[i].Data());
    fDecayVtxv0Z[i] = new TH1F(decayVtxZname.Data(), decayVtxZname.Data(), 400,
                               0., 200);
    fDecayVtxv0Z[i]->GetXaxis()->SetTitle("Decay Vtx To PV Z (cm)");
    fDecayVtxv0Z[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDecayVtxv0Z[i]);

    TString transverseRadname = Form("TransverseRadius_%s", sName[i].Data());
    fTransRadius[i] = new TH1F(transverseRadname.Data(),
                               transverseRadname.Data(), 750, 0, 150);
    fTransRadius[i]->GetXaxis()->SetTitle("#it{r}_{xy} (cm)");
    fTransRadius[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fTransRadius[i]);

    TString DCADauPVPname = Form("DCADauPToPV_%s", sName[i].Data());
    fDCAPosDaugToPrimVtx[i] = new TH1F(DCADauPVPname.Data(),
                                       DCADauPVPname.Data(), 500, 0, 100);
    fDCAPosDaugToPrimVtx[i]->GetXaxis()->SetTitle("DCADaugther P to PV");
    fDCAPosDaugToPrimVtx[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDCAPosDaugToPrimVtx[i]);

    TString DCADauPVNname = Form("DCADauNToPV_%s", sName[i].Data());
    fDCANegDaugToPrimVtx[i] = new TH1F(DCADauPVNname.Data(),
                                       DCADauPVNname.Data(), 500, 0, 100);
    fDCANegDaugToPrimVtx[i]->GetXaxis()->SetTitle("DCADaugther N to PV");
    fDCANegDaugToPrimVtx[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDCANegDaugToPrimVtx[i]);

    TString DCADaugVtxname = Form("DCADauToVtx_%s", sName[i].Data());
    fDCADaugToVtx[i] = new TH1F(DCADaugVtxname.Data(), DCADaugVtxname.Data(),
                                100, 0, 5);
    fDCADaugToVtx[i]->GetXaxis()->SetTitle("DCA Daug to Vtx");
    fDCADaugToVtx[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fDCADaugToVtx[i]);

    TString cosPointName = Form("PointingAngle_%s", sName[i].Data());
    fCPA[i] = new TH1F(cosPointName.Data(), cosPointName.Data(), 500, 0.9,
                       1.001);
    fCPA[i]->GetXaxis()->SetTitle("cos(#alpha)");
    fCPA[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fCPA[i]);

    TString armenterosName = Form("ArmenterosPodolandski_%s", sName[i].Data());
    fArmenterosPodolandski[i] = new TH2F(
        armenterosName.Data(), armenterosName.Data(), 200, -1, 1, 100, 0, 0.25);
    fArmenterosPodolandski[i]->GetXaxis()->SetTitle("#alpha");
    fArmenterosPodolandski[i]->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");
    fv0CutQA[i]->Add(fArmenterosPodolandski[i]);

    TString invMassName = Form("InvariantMass_%s", sName[i].Data());
    fInvMass[i] = new TH1F(invMassName.Data(), invMassName.Data(), MassNBins,
                           MassMin, MassMax);
    fInvMass[i]->GetXaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
    fInvMass[i]->GetYaxis()->SetTitle("Entries");
    fv0CutQA[i]->Add(fInvMass[i]);
  }

  if (CPAPlots) {
    fCPAPtBins = new TH2F("CPAPtBinsTot", "CPAPtBinsTot", 8, 0.3, 4.3, 1000,
                          0.90, 1.);
    fCPAPtBins->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fCPAPtBins->GetYaxis()->SetTitle("cos(#alpha)");
    fHistList->Add(fCPAPtBins);

    TString cpaPtBinName = "CPAPtBinsMult_0_";
    cpaPtBinName += fMultRangeLow;
    TString cpaAxisName = "0 < mult < ";
    cpaAxisName += fMultRangeLow;
    cpaAxisName += ";#it{p}_{T} (GeV/#it{c});cos(#alpha)";
    fCPAPtBinsMult[0] = new TH2F(cpaPtBinName.Data(), cpaAxisName.Data(), 8,
                                 0.3, 4.3, 1000, 0.90, 1.);

    cpaPtBinName = "CPAPtBinsMult_";
    cpaPtBinName += fMultRangeLow;
    cpaPtBinName += "_";
    cpaPtBinName += fMultRangeHigh;
    cpaAxisName = fMultRangeLow;
    cpaAxisName += "0 < mult < ";
    cpaAxisName += fMultRangeHigh;
    cpaAxisName += ";#it{p}_{T} (GeV/#it{c});cos(#alpha)";
    fCPAPtBinsMult[1] = new TH2F(cpaPtBinName.Data(), cpaAxisName.Data(), 8,
                                 0.3, 4.3, 1000, 0.90, 1.);

    cpaPtBinName = "CPAPtBinsMult_";
    cpaPtBinName += fMultRangeHigh;
    cpaPtBinName += "_inf";
    cpaAxisName = "mult > ";
    cpaAxisName += fMultRangeHigh;
    cpaAxisName += ";#it{p}_{T} (GeV/#it{c});cos(#alpha)";
    fCPAPtBinsMult[2] = new TH2F(cpaPtBinName.Data(), cpaAxisName.Data(), 8,
                                 0.3, 4.3, 1000, 0.90, 1.);

    fHistList->Add(fCPAPtBinsMult[0]);
    fHistList->Add(fCPAPtBinsMult[1]);
    fHistList->Add(fCPAPtBinsMult[2]);

  } else {
    fCPAPtBins = nullptr;
  }
  if (perRunnumber) {
    int nBins = iRunMax - iRunMin;
    if (nBins > 2000) {
        std::cout << "Grouping Run Numbers in the Run Number vs. Invariant Mass Plots\n";
      nBins/=10;
    }
    TString InvMassRunNumbName = "InvMassPerRunnumber";
    fInvMassPerRunNumber = new TH2F(InvMassRunNumbName.Data(),
                                    InvMassRunNumbName.Data(), nBins, iRunMin,
                                    iRunMax, 200, 1.08,1.16);
    fHistList->Add(fInvMassPerRunNumber);
  } else {
    fInvMassPerRunNumber = nullptr;
  }
}

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist(TString MinimalBooking, int MassNBins,
                                         float MassMin, float MassMax)
    : fMinimalBooking(true),
      fMultRangeLow(27),
      fMultRangeHigh(55),
      fConfig(0),
      fCutCounter(0),
      fInvMassBefKaonRej(0),
      fInvMassKaon(0),
      fInvMassBefSelection(0),
      fInvMassPt(0),
      fCPAPtBins(0),
      fCPAPtBinsMult(),
      fInvMassPerRunNumber() {
  for (int i = 0; i < 2; ++i) {
    fv0CutQA[i] = nullptr;
    fOnFly[i] = nullptr;
    fpTDist[i] = nullptr;
    fetaDist[i] = nullptr;
    fDecayVtxv0X[i] = nullptr;
    fDecayVtxv0Y[i] = nullptr;
    fDecayVtxv0Z[i] = nullptr;
    fTransRadius[i] = nullptr;
    fDCAPosDaugToPrimVtx[i] = nullptr;
    fDCANegDaugToPrimVtx[i] = nullptr;
    fDCADaugToVtx[i] = nullptr;
    fCPA[i] = nullptr;
    fInvMass[i] = nullptr;
  }
  fHistList = new TList();
  fHistList->SetOwner();
  fHistList->SetName(MinimalBooking.Data());

  fInvMassPt = new TH2F("InvMassPt", "Invariant Mass in Pt Bins", 8, 0.3, 4.3,
                        MassNBins, MassMin, MassMax);
  fInvMassPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c}^{2})");
  fInvMassPt->GetYaxis()->SetTitle("#it{M}_{Pair} (GeV/#it{c}^{2})");
  fHistList->Add(fInvMassPt);
}

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist(const AliFemtoDreamv0Hist& hists)
    : fMinimalBooking(hists.fMinimalBooking),
      fHistList(hists.fHistList),
      fMultRangeLow(hists.fMultRangeLow),
      fMultRangeHigh(hists.fMultRangeHigh),
      fConfig(hists.fConfig),
      fCutCounter(hists.fCutCounter),
      fInvMassBefKaonRej(hists.fInvMassBefKaonRej),
      fInvMassKaon(hists.fInvMassKaon),
      fInvMassBefSelection(hists.fInvMassBefSelection),
      fInvMassPt(hists.fInvMassPt),
      fCPAPtBins(hists.fCPAPtBins),
      fInvMassPerRunNumber(hists.fInvMassPerRunNumber) {
  for (int i = 0; i < 2; ++i) {
    fv0CutQA[i] = hists.fv0CutQA[i];
    fOnFly[i] = hists.fOnFly[i];
    fpTDist[i] = hists.fpTDist[i];
    fetaDist[i] = hists.fetaDist[i];
    fDecayVtxv0X[i] = hists.fDecayVtxv0X[i];
    fDecayVtxv0Y[i] = hists.fDecayVtxv0Y[i];
    fDecayVtxv0Z[i] = hists.fDecayVtxv0Z[i];
    fTransRadius[i] = hists.fTransRadius[i];
    fDCAPosDaugToPrimVtx[i] = hists.fDCAPosDaugToPrimVtx[i];
    fDCANegDaugToPrimVtx[i] = hists.fDCANegDaugToPrimVtx[i];
    fDCADaugToVtx[i] = hists.fDCADaugToVtx[i];
    fCPA[i] = hists.fCPA[i];
    fInvMass[i] = hists.fInvMass[i];
  }
  for (int i = 0; i < 3; ++i) {
    fCPAPtBinsMult[i] = hists.fCPAPtBinsMult[i];
  }
}

void AliFemtoDreamv0Hist::FillCPAPtBins(float pT, float cpa, int multiplicity) {
  if (!fMinimalBooking) {
    fCPAPtBins->Fill(pT, cpa);
    if (multiplicity < fMultRangeLow) {
      fCPAPtBinsMult[0]->Fill(pT, cpa);
    } else if (multiplicity >= fMultRangeLow && multiplicity < fMultRangeHigh) {
      fCPAPtBinsMult[1]->Fill(pT, cpa);
    } else {
      fCPAPtBinsMult[2]->Fill(pT, cpa);
    }
  }
}

AliFemtoDreamv0Hist& AliFemtoDreamv0Hist::operator=(
    const AliFemtoDreamv0Hist& hists) {
  if (this != &hists) {
    this->fMinimalBooking = hists.fMinimalBooking;
    this->fHistList = hists.fHistList;
    this->fMultRangeLow = hists.fMultRangeLow;
    this->fMultRangeHigh = hists.fMultRangeHigh;
    this->fConfig = hists.fConfig;
    this->fCutCounter = hists.fCutCounter;
    this->fInvMassBefKaonRej = hists.fInvMassBefKaonRej;
    this->fInvMassKaon = hists.fInvMassKaon;
    this->fInvMassBefSelection = hists.fInvMassBefSelection;
    this->fInvMassPt = hists.fInvMassPt;
    this->fCPAPtBins = hists.fCPAPtBins;
    this->fInvMassPerRunNumber = hists.fInvMassPerRunNumber;
    for (int i = 0; i < 2; ++i) {
      this->fv0CutQA[i] = hists.fv0CutQA[i];
      this->fOnFly[i] = hists.fOnFly[i];
      this->fpTDist[i] = hists.fpTDist[i];
      this->fetaDist[i] = hists.fetaDist[i];
      this->fDecayVtxv0X[i] = hists.fDecayVtxv0X[i];
      this->fDecayVtxv0Y[i] = hists.fDecayVtxv0Y[i];
      this->fDecayVtxv0Z[i] = hists.fDecayVtxv0Z[i];
      this->fTransRadius[i] = hists.fTransRadius[i];
      this->fDCAPosDaugToPrimVtx[i] = hists.fDCAPosDaugToPrimVtx[i];
      this->fDCANegDaugToPrimVtx[i] = hists.fDCANegDaugToPrimVtx[i];
      this->fDCADaugToVtx[i] = hists.fDCADaugToVtx[i];
      this->fCPA[i] = hists.fCPA[i];
      this->fInvMass[i] = hists.fInvMass[i];
    }
    for (int i = 0; i < 3; ++i) {
      fCPAPtBinsMult[i] = hists.fCPAPtBinsMult[i];
    }
  }
  return *this;
}

AliFemtoDreamv0Hist::~AliFemtoDreamv0Hist() {
  delete fHistList;
}

