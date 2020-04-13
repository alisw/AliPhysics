/*
 * AliFemtoDreamEventHist.cxx
 *
 *  Created on: Nov 22,2017
 *      Author: gu74req
 */

#include "AliFemtoDreamEventHist.h"
ClassImp(AliFemtoDreamEventHist)

AliFemtoDreamEventHist::AliFemtoDreamEventHist()
    : fEventCutList(0),
      fEvtCounter(0),
      fCutConfig(0),
      fCentVsMultPlots(false),
      fCentVsV0A(0),
      fCentVsV0M(0),
      fCentVsV0C(0),
      fCentVsRefMult(0) {
  for (int i = 0; i < 2; ++i) {
    fEvtNCont[i] = nullptr;
    fEvtVtxX[i] = nullptr;
    fEvtVtxY[i] = nullptr;
    fEvtVtxZ[i] = nullptr;
    fSPDTrklClsLy0[i] = nullptr;
    fSPDTrklClsLy1[i] = nullptr;
    fSPDTrklClsLySum[i] = nullptr;
    fSPDTrklClsLySum[i] = nullptr;
    fMultDistSPD[i] = nullptr;
    fMultDistV0A[i] = nullptr;
    fMultDistV0C[i] = nullptr;
    fMultDistV0M[i] = nullptr;
    fMultDistRef08[i] = nullptr;
    fMultPercentV0[i] = nullptr;
    fEvtSpher[i] = nullptr;
    fEvtSphero[i] = nullptr;
    fPileUpVZEROTime[i] = nullptr;
  }
}
AliFemtoDreamEventHist::AliFemtoDreamEventHist(bool centVsMultPlot) {
  fEventCutList = new TList();
  fEventCutList->SetName("Event Cuts");
  fEventCutList->SetOwner();

  fEvtCounter = new TH1F("EventCounter", "Event Counter", 13, 0, 13);
  fEvtCounter->GetXaxis()->SetBinLabel(1, "Events");
  fEvtCounter->GetXaxis()->SetBinLabel(2, "AliEventCuts");
  fEvtCounter->GetXaxis()->SetBinLabel(3, "Phys. Sel.");
  fEvtCounter->GetXaxis()->SetBinLabel(4, "nContrib");
  fEvtCounter->GetXaxis()->SetBinLabel(5, "zVtx");
  fEvtCounter->GetXaxis()->SetBinLabel(6, "PileUp");
  fEvtCounter->GetXaxis()->SetBinLabel(7, "SPD Mult Cleanup");
  fEvtCounter->GetXaxis()->SetBinLabel(8, "V0A Mult Cleanup");
  fEvtCounter->GetXaxis()->SetBinLabel(9, "V0C Mult Cleanup");
  fEvtCounter->GetXaxis()->SetBinLabel(10, "RefMult08 Cleanup");
  fEvtCounter->GetXaxis()->SetBinLabel(11, "Sphericity");
  fEvtCounter->GetXaxis()->SetBinLabel(12, "Mult. percentile");
  fEvtCounter->GetXaxis()->SetBinLabel(13, "Spherocity");
  fEvtCounter->GetYaxis()->SetTitle("Entries");

  fEventCutList->Add(fEvtCounter);

  fCutConfig = new TProfile("CutConfig", "Cut Config", 21, 0, 21);
  fCutConfig->GetXaxis()->SetBinLabel(1, "Min Contrib");
  fCutConfig->GetXaxis()->SetBinLabel(2, "CutZvtx");
  fCutConfig->GetXaxis()->SetBinLabel(3, "Min Zvtx");
  fCutConfig->GetXaxis()->SetBinLabel(4, "Max ZVtx");
  fCutConfig->GetXaxis()->SetBinLabel(5, "PileUp Rejection");
  fCutConfig->GetXaxis()->SetBinLabel(6, "MV PileUp Rejection");
  fCutConfig->GetXaxis()->SetBinLabel(7, "SPD Mult");
  fCutConfig->GetXaxis()->SetBinLabel(8, "V0A Mult");
  fCutConfig->GetXaxis()->SetBinLabel(9, "V0C Mult");
  fCutConfig->GetXaxis()->SetBinLabel(10, "RefMult08");
  fCutConfig->GetXaxis()->SetBinLabel(11, "AliEvtCuts");
  fCutConfig->GetXaxis()->SetBinLabel(12, "Low Spher");
  fCutConfig->GetXaxis()->SetBinLabel(13, "Up Spher");
  fCutConfig->GetXaxis()->SetBinLabel(14, "Mult. percentile");
  fCutConfig->GetXaxis()->SetBinLabel(15, "Low Sphero");
  fCutConfig->GetXaxis()->SetBinLabel(16, "Up Sphero");
  fEventCutList->Add(fCutConfig);

  fCentVsMultPlots = centVsMultPlot;
  if (fCentVsMultPlots) {
    TString vsV0AName = Form("CentvsV0A");
    fCentVsV0A = new TH2F(vsV0AName.Data(), vsV0AName.Data(), 200, -0.5, 99.5,
                          300, 0.5, 600.5);
    fEventCutList->Add(fCentVsV0A);

    TString vsV0MName = Form("CentvsV0M");
    fCentVsV0M = new TH2F(vsV0MName.Data(), vsV0MName.Data(), 200, -0.5, 99.5,
                          300, 0.5, 600.5);
    fEventCutList->Add(fCentVsV0M);

    TString vsV0CName = Form("CentvsV0C");
    fCentVsV0C = new TH2F(vsV0CName.Data(), vsV0CName.Data(), 200, -0.5, 99.5,
                          300, 0.5, 600.5);
    fEventCutList->Add(fCentVsV0C);

    TString vsV0RefName = Form("CentvsRefMult");
    fCentVsRefMult = new TH2F(vsV0RefName.Data(), vsV0RefName.Data(), 200, -0.5, 99.5, 100, 0.5, 200.5);
    fEventCutList->Add(fCentVsRefMult);
  } else {
    fCentVsV0A = 0;
    fCentVsV0M = 0;
    fCentVsV0C = 0;
    fCentVsRefMult = 0;
  }

  TString sName[2] = { "before", "after" };

  for (int i = 0; i < 2; ++i) {

    fEvtCutQA[i] = new TList();
    fEvtCutQA[i]->SetName(sName[i].Data());
    fEvtCutQA[i]->SetOwner();

    fEventCutList->Add(fEvtCutQA[i]);

    TString nEvtNContName = Form("nContributors_%s", sName[i].Data());
    fEvtNCont[i] = new TH1F(nEvtNContName.Data(), nEvtNContName.Data(), 350.,
                            -0.5, 349.5);
    fEvtNCont[i]->GetXaxis()->SetTitle("Number of Contributors");
    fEvtNCont[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtNCont[i]);

    TString EvtVtxXName = Form("VtxX_%s", sName[i].Data());
    fEvtVtxX[i] = new TH1F(EvtVtxXName.Data(), EvtVtxXName.Data(), 100, -0.15,
                           0.35);
    fEvtVtxX[i]->GetXaxis()->SetTitle("vtx_{x} (cm)");
    fEvtVtxX[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtVtxX[i]);

    TString EvtVtxYName = Form("VtxY_%s", sName[i].Data());
    fEvtVtxY[i] = new TH1F(EvtVtxYName.Data(), EvtVtxYName.Data(), 100, 0.1,
                           0.6);
    fEvtVtxY[i]->GetXaxis()->SetTitle("vtx_{y} (cm)");
    fEvtVtxY[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtVtxY[i]);

    TString EvtVtxZName = Form("VtxZ_%s", sName[i].Data());
    fEvtVtxZ[i] = new TH1F(EvtVtxZName.Data(), EvtVtxZName.Data(), 300, -12.5,
                           12.5);
    fEvtVtxZ[i]->GetXaxis()->SetTitle("vtx_{z} (cm)");
    fEvtVtxZ[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtVtxZ[i]);

    TString MultNameSPD = Form("MultiplicitySPD_%s", sName[i].Data());
    fMultDistSPD[i] = new TH1F(MultNameSPD.Data(), MultNameSPD.Data(), 600, 0.,
                               600.);
    fMultDistSPD[i]->GetXaxis()->SetTitle("Multiplicity (SPD)");
    fMultDistSPD[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultDistSPD[i]);

    TString MultNameV0A = Form("MultiplicityV0A_%s", sName[i].Data());
    fMultDistV0A[i] = new TH1F(MultNameV0A.Data(), MultNameV0A.Data(), 600, 0.,
                               600.);
    fMultDistV0A[i]->GetXaxis()->SetTitle("Multiplicity (V0A)");
    fMultDistV0A[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultDistV0A[i]);

    TString MultNameV0C = Form("MultiplicityV0C_%s", sName[i].Data());
    fMultDistV0C[i] = new TH1F(MultNameV0C.Data(), MultNameV0C.Data(), 600, 0.,
                               600.);
    fMultDistV0C[i]->GetXaxis()->SetTitle("Multiplicity (V0C)");
    fMultDistV0C[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultDistV0C[i]);

    TString MultNameV0M = Form("MultiplicityV0M_%s", sName[i].Data());
    fMultDistV0M[i] = new TH1F(MultNameV0M.Data(), MultNameV0M.Data(), 600, 0.,
                               600.);
    fMultDistV0M[i]->GetXaxis()->SetTitle("Multiplicity (V0M)");
    fMultDistV0M[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultDistV0M[i]);

    TString MultNameRefMult08 = Form("MultiplicityRef08_%s", sName[i].Data());
    fMultDistRef08[i] = new TH1F(MultNameRefMult08.Data(),
                                 MultNameRefMult08.Data(), 600, 0., 600.);
    fMultDistRef08[i]->GetXaxis()->SetTitle("Multiplicity (RefMult08)");
    fMultDistRef08[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultDistRef08[i]);

    TString MultNameMultPercentV0 = Form("MultiplicityPercentileV0_%s", sName[i].Data());
    fMultPercentV0[i] = new TH1F(MultNameMultPercentV0.Data(),
                                 MultNameMultPercentV0.Data(), 10000, 0., 100.);
    fMultPercentV0[i]->GetXaxis()->SetTitle("Multiplicity precentile V0");
    fMultPercentV0[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fMultPercentV0[i]);

    TString SPDtrklClsLy0Name = Form("SPDTrackletsVsClusterL0_%s",
                                     sName[i].Data());
    fSPDTrklClsLy0[i] = new TH2F(SPDtrklClsLy0Name.Data(),
                                 SPDtrklClsLy0Name.Data(), 250, 0, 250, 1000, 0,
                                 1000);
    fSPDTrklClsLy0[i]->GetXaxis()->SetTitle("SPD Tracklets");
    fSPDTrklClsLy0[i]->GetYaxis()->SetTitle("SPD Cluster L0");
    fEvtCutQA[i]->Add(fSPDTrklClsLy0[i]);

    TString SPDtrklClsLy1Name = Form("SPDTrackletsVsClusterL1_%s",
                                     sName[i].Data());
    fSPDTrklClsLy1[i] = new TH2F(SPDtrklClsLy1Name.Data(),
                                 SPDtrklClsLy1Name.Data(), 250, 0, 250, 1000, 0,
                                 1000);
    fSPDTrklClsLy1[i]->GetXaxis()->SetTitle("SPD Tracklets");
    fSPDTrklClsLy1[i]->GetYaxis()->SetTitle("SPD Cluster L1");
    fEvtCutQA[i]->Add(fSPDTrklClsLy1[i]);

    TString SPDtrklClsSumName = Form("SPDTrackletsVsClusterL01Sum_%s",
                                     sName[i].Data());
    fSPDTrklClsLySum[i] = new TH2F(SPDtrklClsSumName.Data(),
                                   SPDtrklClsSumName.Data(), 250, 0, 250, 1000,
                                   0, 1000);
    fSPDTrklClsLySum[i]->GetXaxis()->SetTitle("SPD Tracklets");
    fSPDTrklClsLySum[i]->GetYaxis()->SetTitle("SPD Cluster Sum (L0+L1)");
    fEvtCutQA[i]->Add(fSPDTrklClsLySum[i]);

    TString SPDvsTrkZVtxName = Form("SPDvsTrackZVtxPos_%s", sName[i].Data());
    fSPDTrackZVtx[i] = new TH2F(SPDvsTrkZVtxName.Data(),
                                SPDvsTrkZVtxName.Data(), 300, -15, 15, 300, -15,
                                15);
    fSPDTrackZVtx[i]->GetXaxis()->SetTitle("vtx_{z, SPD} (cm)");
    fSPDTrackZVtx[i]->GetYaxis()->SetTitle("vtx_{z, Tracks} (cm)");
    fEvtCutQA[i]->Add(fSPDTrackZVtx[i]);

    TString SPDTrkZVtxDisplName = Form("SPDTrackZVtxDisplacement%s",
                                       sName[i].Data());
    fSPDTrkZVtxDispl[i] = new TH1F(SPDTrkZVtxDisplName.Data(),
                                   SPDTrkZVtxDisplName.Data(), 300, 0, 1.5);
    fSPDTrkZVtxDispl[i]->GetXaxis()->SetTitle("|vtx_{z, SPD} - vtx_{z, Tracks}| (cm)");
    fSPDTrkZVtxDispl[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fSPDTrkZVtxDispl[i]);

    TString BFieldName = Form("MagneticFieldkGauss_%s", sName[i].Data());
    fBField[i] = new TH1F(BFieldName.Data(), BFieldName.Data(), 20, -10, 10);
    fBField[i]->GetXaxis()->SetTitle("B (G)");
    fBField[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fBField[i]);

    TString EvtSpherName = Form("Sphericity_%s", sName[i].Data());
    fEvtSpher[i] = new TH1F(EvtSpherName.Data(), EvtSpherName.Data(), 50, 0.,
                            1.);
    fEvtSpher[i]->GetXaxis()->SetTitle("Sphericity #it{S}_{T}");
    fEvtSpher[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtSpher[i]);

    TString EvtSpheroName = Form("Spherocity_%s", sName[i].Data());
    fEvtSphero[i] = new TH1F(EvtSpheroName.Data(), EvtSpheroName.Data(), 50, 0.,
                            1.);
    fEvtSphero[i]->GetXaxis()->SetTitle("Spherocity #it{S}_{0}");
    fEvtSphero[i]->GetYaxis()->SetTitle("Entries");
    fEvtCutQA[i]->Add(fEvtSphero[i]);

    TString PileupName = Form("VZEROtiming_%s", sName[i].Data());
    fPileUpVZEROTime[i] = new TH2F(PileupName.Data(), PileupName.Data(), 500,
                                   -20, 30, 500, -20, 30);
    fPileUpVZEROTime[i]->GetXaxis()->SetTitle("#it{t}_{V0A} - #it{t}_{V0C}");
    fPileUpVZEROTime[i]->GetYaxis()->SetTitle("#it{t}_{V0A} + #it{t}_{V0C}");
    fEvtCutQA[i]->Add(fPileUpVZEROTime[i]);
  }
}

AliFemtoDreamEventHist::AliFemtoDreamEventHist(
    const AliFemtoDreamEventHist& hists)
    : fEventCutList(hists.fEventCutList),
      fEvtCounter(hists.fEvtCounter),
      fCutConfig(hists.fCutConfig),
      fCentVsMultPlots(hists.fCentVsMultPlots),
      fCentVsV0A(hists.fCentVsV0A),
      fCentVsV0M(hists.fCentVsV0M),
      fCentVsV0C(hists.fCentVsV0C),
      fCentVsRefMult(hists.fCentVsRefMult) {
  for (int i = 0; i < 2; ++i) {
    fEvtNCont[i] = hists.fEvtNCont[i];
    fEvtVtxX[i] = hists.fEvtVtxX[i];
    fEvtVtxY[i] = hists.fEvtVtxY[i];
    fEvtVtxZ[i] = hists.fEvtVtxZ[i];
    fSPDTrklClsLy0[i] = hists.fSPDTrklClsLy0[i];
    fSPDTrklClsLy1[i] = hists.fSPDTrklClsLy1[i];
    fSPDTrklClsLySum[i] = hists.fSPDTrklClsLySum[i];
    fMultDistSPD[i] = hists.fMultDistSPD[i];
    fMultDistV0A[i] = hists.fMultDistV0A[i];
    fMultDistV0C[i] = hists.fMultDistV0C[i];
    fMultDistV0M[i] = hists.fMultDistV0M[i];
    fMultDistRef08[i] = hists.fMultDistRef08[i];
    fMultPercentV0[i] = hists.fMultPercentV0[i];
    fBField[i] = hists.fBField[i];
    fEvtSpher[i] = hists.fEvtSpher[i];
    fEvtSphero[i] = hists.fEvtSphero[i];
    fPileUpVZEROTime[i] = hists.fPileUpVZEROTime[i];
  }
}
AliFemtoDreamEventHist& AliFemtoDreamEventHist::operator=(
    const AliFemtoDreamEventHist& hists) {
  if (this != &hists) {
    this->fEventCutList = hists.fEventCutList;
    this->fEvtCounter = hists.fEvtCounter;
    this->fCutConfig = hists.fCutConfig;
    this->fCentVsMultPlots = hists.fCentVsMultPlots;
    this->fCentVsV0A = hists.fCentVsV0A;
    this->fCentVsV0M = hists.fCentVsV0M;
    this->fCentVsV0C = hists.fCentVsV0C;
    this->fCentVsRefMult = hists.fCentVsRefMult;
    for (int i = 0; i < 2; ++i) {
      this->fEvtNCont[i] = hists.fEvtNCont[i];
      this->fEvtVtxX[i] = hists.fEvtVtxX[i];
      this->fEvtVtxY[i] = hists.fEvtVtxY[i];
      this->fEvtVtxZ[i] = hists.fEvtVtxZ[i];
      this->fSPDTrklClsLy0[i] = hists.fSPDTrklClsLy0[i];
      this->fSPDTrklClsLy1[i] = hists.fSPDTrklClsLy1[i];
      this->fSPDTrklClsLySum[i] = hists.fSPDTrklClsLySum[i];
      this->fMultDistSPD[i] = hists.fMultDistSPD[i];
      this->fMultDistV0A[i] = hists.fMultDistV0A[i];
      this->fMultDistV0C[i] = hists.fMultDistV0C[i];
      this->fMultDistV0M[i] = hists.fMultDistV0M[i];
      this->fMultDistRef08[i] = hists.fMultDistRef08[i];
      this->fMultPercentV0[i] = hists.fMultPercentV0[i];
      this->fBField[i] = hists.fBField[i];
      this->fEvtSpher[i] = hists.fEvtSpher[i];
      this->fEvtSphero[i] = hists.fEvtSphero[i];
      this->fPileUpVZEROTime[i] = hists.fPileUpVZEROTime[i];
    }
  }
  return *this;
}

AliFemtoDreamEventHist::~AliFemtoDreamEventHist() {

}
