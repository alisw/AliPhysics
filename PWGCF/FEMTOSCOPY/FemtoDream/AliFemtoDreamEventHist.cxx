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
    fSPDTrklCls[i] = nullptr;
    fMultDistSPD[i] = nullptr;
    fMultDistV0A[i] = nullptr;
    fMultDistV0C[i] = nullptr;
    fMultDistRef08[i] = nullptr;
  }
}
AliFemtoDreamEventHist::AliFemtoDreamEventHist(bool centVsMultPlot) {
  fEventCutList = new TList();
  fEventCutList->SetName("Event Cuts");
  fEventCutList->SetOwner();

  fEvtCounter = new TH1F("EventCounter", "Event Counter", 10, 0, 10);
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
  fEventCutList->Add(fEvtCounter);

  fCutConfig = new TProfile("CutConfig", "Cut Config", 20, 0, 20);
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
  fEventCutList->Add(fCutConfig);

  fCentVsMultPlots = centVsMultPlot;
  if (fCentVsMultPlots) {
    TString vsV0AName = Form("CentvsV0A");
    fCentVsV0A = new TH2F(vsV0AName.Data(), vsV0AName.Data(), 100, 0.5, 100,
                          300, 0.5, 600.5);
    fCentVsV0A->Sumw2();
    fEventCutList->Add(fCentVsV0A);

    TString vsV0MName = Form("CentvsV0M");
    fCentVsV0M = new TH2F(vsV0MName.Data(), vsV0MName.Data(), 100, 0.5, 100,
                          300, 0.5, 600.5);
    fCentVsV0M->Sumw2();
    fEventCutList->Add(fCentVsV0M);

    TString vsV0CName = Form("CentvsV0C");
    fCentVsV0C = new TH2F(vsV0CName.Data(), vsV0CName.Data(), 100, 0.5, 100,
                          300, 0.5, 600.5);
    fCentVsV0C->Sumw2();
    fEventCutList->Add(fCentVsV0C);

    TString vsV0RefName = Form("CentvsRefMult");
    fCentVsRefMult = new TH2F(vsV0RefName.Data(), vsV0RefName.Data(), 100, 0.5,
                              100, 100, 0.5, 200.5);
    fCentVsRefMult->Sumw2();
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
    fEvtNCont[i]->Sumw2();
    fEvtNCont[i]->GetXaxis()->SetTitle("Number of Contributors");
    fEvtCutQA[i]->Add(fEvtNCont[i]);

    TString EvtVtxXName = Form("VtxX_%s", sName[i].Data());
    fEvtVtxX[i] = new TH1F(EvtVtxXName.Data(), EvtVtxXName.Data(), 50, -2.5,
                           2.5);
    fEvtVtxX[i]->Sumw2();
    fEvtVtxX[i]->GetXaxis()->SetTitle("DCA_{x}");
    fEvtCutQA[i]->Add(fEvtVtxX[i]);

    TString EvtVtxYName = Form("VtxY_%s", sName[i].Data());
    fEvtVtxY[i] = new TH1F(EvtVtxYName.Data(), EvtVtxYName.Data(), 50, -2.5,
                           2.5);
    fEvtVtxY[i]->Sumw2();
    fEvtVtxY[i]->GetXaxis()->SetTitle("DCA_{y}");
    fEvtCutQA[i]->Add(fEvtVtxY[i]);

    TString EvtVtxZName = Form("VtxZ_%s", sName[i].Data());
    fEvtVtxZ[i] = new TH1F(EvtVtxZName.Data(), EvtVtxZName.Data(), 300, -15.,
                           15.);
    fEvtVtxZ[i]->Sumw2();
    fEvtVtxZ[i]->GetXaxis()->SetTitle("DCA_{z}");
    fEvtCutQA[i]->Add(fEvtVtxZ[i]);

    TString MultNameSPD = Form("MultiplicitySPD_%s", sName[i].Data());
    fMultDistSPD[i] = new TH1F(MultNameSPD.Data(), MultNameSPD.Data(), 600, 0.,
                               600.);
    fMultDistSPD[i]->Sumw2();
    fMultDistSPD[i]->GetXaxis()->SetTitle("Multiplicity (SPD)");
    fEvtCutQA[i]->Add(fMultDistSPD[i]);

    TString MultNameV0A = Form("MultiplicityV0A_%s", sName[i].Data());
    fMultDistV0A[i] = new TH1F(MultNameV0A.Data(), MultNameV0A.Data(), 600, 0.,
                               600.);
    fMultDistV0A[i]->Sumw2();
    fMultDistV0A[i]->GetXaxis()->SetTitle("Multiplicity (V0A)");
    fEvtCutQA[i]->Add(fMultDistV0A[i]);

    TString MultNameV0C = Form("MultiplicityV0C_%s", sName[i].Data());
    fMultDistV0C[i] = new TH1F(MultNameV0C.Data(), MultNameV0C.Data(), 600, 0.,
                               600.);
    fMultDistV0C[i]->Sumw2();
    fMultDistV0C[i]->GetXaxis()->SetTitle("Multiplicity (V0C)");
    fEvtCutQA[i]->Add(fMultDistV0C[i]);

    TString MultNameRefMult08 = Form("MultiplicityRef08_%s", sName[i].Data());
    fMultDistRef08[i] = new TH1F(MultNameRefMult08.Data(),
                                 MultNameRefMult08.Data(), 600, 0., 600.);
    fMultDistRef08[i]->Sumw2();
    fMultDistRef08[i]->GetXaxis()->SetTitle("Multiplicity (RefMult08)");
    fEvtCutQA[i]->Add(fMultDistRef08[i]);

    TString SPDtrklClsName = Form("SPDTrackletsVsCluster_%s", sName[i].Data());
    fSPDTrklCls[i] = new TH2F(SPDtrklClsName.Data(), SPDtrklClsName.Data(), 250,
                              0, 250, 1000, 0, 1000);
    fSPDTrklCls[i]->Sumw2();
    fSPDTrklCls[i]->GetXaxis()->SetTitle("SPD Tracklets");
    fSPDTrklCls[i]->GetYaxis()->SetTitle("SPD Cluster");
    fEvtCutQA[i]->Add(fSPDTrklCls[i]);

    TString SPDvsTrkZVtxName = Form("SPDvsTrackZVtxPos_%s", sName[i].Data());
    fSPDTrackZVtx[i] = new TH2F(SPDvsTrkZVtxName.Data(),
                                SPDvsTrkZVtxName.Data(), 300, -15, 15, 300, -15,
                                15);
    fSPDTrackZVtx[i]->Sumw2();
    fSPDTrackZVtx[i]->GetXaxis()->SetTitle("zVtx Position SPD");
    fSPDTrackZVtx[i]->GetYaxis()->SetTitle("zVtx Position Tracks");
    fEvtCutQA[i]->Add(fSPDTrackZVtx[i]);

    TString SPDTrkZVtxDisplName = Form("SPDTrackZVtxDisplacement%s",
                                       sName[i].Data());
    fSPDTrkZVtxDispl[i] = new TH1F(SPDTrkZVtxDisplName.Data(),
                                   SPDTrkZVtxDisplName.Data(), 300, 0, 1.5);
    fSPDTrkZVtxDispl[i]->Sumw2();
    fSPDTrkZVtxDispl[i]->GetXaxis()->SetTitle("zVtx Position |SPD - Tracks|");
    fEvtCutQA[i]->Add(fSPDTrkZVtxDispl[i]);
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
    fSPDTrklCls[i] = hists.fSPDTrklCls[i];
    fMultDistSPD[i] = hists.fMultDistSPD[i];
    fMultDistV0A[i] = hists.fMultDistV0A[i];
    fMultDistV0C[i] = hists.fMultDistV0C[i];
    fMultDistRef08[i] = hists.fMultDistRef08[i];
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
      this->fSPDTrklCls[i] = hists.fSPDTrklCls[i];
      this->fMultDistSPD[i] = hists.fMultDistSPD[i];
      this->fMultDistV0A[i] = hists.fMultDistV0A[i];
      this->fMultDistV0C[i] = hists.fMultDistV0C[i];
      this->fMultDistRef08[i] = hists.fMultDistRef08[i];
    }
  }
  return *this;
}

AliFemtoDreamEventHist::~AliFemtoDreamEventHist() {

}

