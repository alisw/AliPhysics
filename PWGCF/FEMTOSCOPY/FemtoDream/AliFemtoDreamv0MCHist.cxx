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
:fMCList(0)
,fCPAPlots(0)
,fMCCorrPt(0)
,fMCIdentPt(0)
,fMCGenPt(0)
,fMCContPt(0)
,fMCUnknownPt(0)
,fMCPrimaryPt(0)
,fMCMaterialPt(0)
,fMCFeeddownWeakPt(0)
,fMCPrimCPAPtBins(0)
,fMCMaterialCPAPtBins(0)
,fMCSecondaryCPAPtBins(0)
,fMCContCPAPtBins(0)
{
  for (int i=0;i<4;++i) {
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
  }
}

AliFemtoDreamv0MCHist::AliFemtoDreamv0MCHist(
    int MassNBins,double MassMin,double MassMax,bool contribSplitting,
    bool CPADist)
{
  double ptmin=0.;
  double ptmax=5.;
  int ptBins=50;

  fMCList=new TList();
  fMCList->SetName("v0MonteCarlo");
  fMCList->SetOwner();

  fMCCorrPt = new TH1F("CorrParPt","Correct Particles Pt",ptBins,ptmin,ptmax);
  fMCCorrPt->Sumw2();
  fMCCorrPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCCorrPt);

  fMCIdentPt = new TH1F("IdentPartPt","Ident Particles Pt",ptBins,ptmin,ptmax);
  fMCIdentPt->Sumw2();
  fMCIdentPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCIdentPt);

  fMCGenPt = new TH1F("GenPartPt","Gen Particles Pt",ptBins,ptmin,ptmax);
  fMCGenPt->Sumw2();
  fMCGenPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCGenPt);

  fMCContPt = new TH1F("ContPt","ContPt",ptBins,ptmin,ptmax);
  fMCContPt->Sumw2();
  fMCContPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCContPt);

  fMCUnknownPt = new TH1F("UnknPt","UnknPt",ptBins,ptmin,ptmax);
  fMCUnknownPt->Sumw2();
  fMCUnknownPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCUnknownPt);

  fMCPrimaryPt = new TH1F("PrimaryPt","PrimaryPt",ptBins,ptmin,ptmax);
  fMCPrimaryPt->Sumw2();
  fMCPrimaryPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCPrimaryPt);

  fMCMaterialPt = new TH1F("MatPt","MatPT",ptBins,ptmin,ptmax);
  fMCMaterialPt->Sumw2();
  fMCMaterialPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCMaterialPt);

  fMCFeeddownWeakPt = new TH2F(
      "FeeddownPt","Feeddown Pt",ptBins,ptmin,ptmax,213,3121,3334);
  fMCFeeddownWeakPt->Sumw2();
  fMCFeeddownWeakPt->GetXaxis()->SetTitle("p_{T}");
  fMCList->Add(fMCFeeddownWeakPt);

  if (contribSplitting) {
    TString MCModes[4] = {"Primary","Secondary","Material","Contamination"};
    for (int i=0;i<4;++i) {
      fMCQAPlots[i] = new TList();
      fMCQAPlots[i]->SetOwner();
      fMCQAPlots[i]->SetName(MCModes[i].Data());
      fMCList->Add(fMCQAPlots[i]);

      TString MCPtDist = Form("MCPt%s",MCModes[i].Data());
      fMCpTDist[i] = new TH1F(MCPtDist.Data(),MCPtDist.Data(),
                              ptBins,ptmin,ptmax);
      fMCpTDist[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCpTDist[i]);

      TString MCEtaDist = Form("MCEta%s",MCModes[i].Data());
      fMCetaDist[i] = new TH1F(MCEtaDist.Data(),MCEtaDist.Data(),
                               200,-10.,10.);
      fMCetaDist[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCetaDist[i]);

      TString MCPhiDist = Form("MCPhi%s",MCModes[i].Data());
      fMCphiDist[i] = new TH1F(MCPhiDist.Data(),MCPhiDist.Data(),
                               100,0.,2*TMath::Pi());
      fMCphiDist[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCphiDist[i]);

      TString MCDecayVtxv0XDist = Form("MCDecayVtxXPV%s",MCModes[i].Data());
      fMCDecayVtxv0X[i] = new TH2F(
          MCDecayVtxv0XDist.Data(),MCDecayVtxv0XDist.Data(),
          ptBins,ptmin,ptmax,400,0,200);
      fMCDecayVtxv0X[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDecayVtxv0X[i]);

      TString MCDecayVtxv0YDist = Form("MCDecayVtxYPV%s",MCModes[i].Data());
      fMCDecayVtxv0Y[i] = new TH2F(
          MCDecayVtxv0YDist.Data(),MCDecayVtxv0YDist.Data(),
          ptBins,ptmin,ptmax,400,0,200);
      fMCDecayVtxv0Y[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDecayVtxv0Y[i]);

      TString MCDecayVtxv0ZDist = Form("MCDecayVtxZPV%s",MCModes[i].Data());
      fMCDecayVtxv0Z[i] = new TH2F(
          MCDecayVtxv0ZDist.Data(),MCDecayVtxv0ZDist.Data(),
          ptBins,ptmin,ptmax,400,0,200);
      fMCDecayVtxv0Z[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDecayVtxv0Z[i]);

      TString MCTransverseRadius =
          Form("MCTransverseRadius%s",MCModes[i].Data());
      fMCTransverseRadius[i] = new TH2F(
          MCTransverseRadius.Data(),MCTransverseRadius.Data(),
          ptBins,ptmin,ptmax,750,0,150);
      fMCTransverseRadius[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCTransverseRadius[i]);

      TString MCDCADaugPVP = Form("MCDCADauPToPV%s",MCModes[i].Data());
      fMCDCAPosDaugToPV[i] = new TH2F(MCDCADaugPVP.Data(),MCDCADaugPVP.Data(),
                                      ptBins,ptmin,ptmax,500,0,100);
      fMCDCAPosDaugToPV[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDCAPosDaugToPV[i]);

      TString MCDCADaugPVN = Form("MCDCADauNToPV%s",MCModes[i].Data());
      fMCDCANegDaugToPV[i] = new TH2F(MCDCADaugPVN.Data(),MCDCADaugPVN.Data(),
                                      ptBins,ptmin,ptmax,500,0,100);
      fMCDCANegDaugToPV[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDCANegDaugToPV[i]);

      TString MCDCADaugToVTX = Form("MCDCADauToVtx%s",MCModes[i].Data());
      fMCDCADaugToVtx[i] = new TH2F(
          MCDCADaugToVTX.Data(),MCDCADaugToVTX.Data(),
          ptBins,ptmin,ptmax,100,0,10);
      fMCDCADaugToVtx[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCDCADaugToVtx[i]);

      TString MCCosPointing = Form("MCPointingAngle%s",MCModes[i].Data());
      fMCCosPointing[i] = new TH2F(
          MCCosPointing.Data(),MCCosPointing.Data(),
          ptBins,ptmin,ptmax,400,0.85,1.001);;
      fMCCosPointing[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCCosPointing[i]);

      TString MCInvMass = Form("MCInvMass%s",MCModes[i].Data());
      fMCInvMass[i] = new TH1F(MCInvMass.Data(),MCInvMass.Data(),
                               MassNBins,MassMin,MassMax);
      fMCInvMass[i]->Sumw2();
      fMCQAPlots[i]->Add(fMCInvMass[i]);
    }
  } else {
    for (int i=0;i<4;++i) {
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
    }
  }

  if (CPADist) {
    fCPAPlots=new TList();
    fCPAPlots->SetName("CPAPtBinning");
    fCPAPlots->SetOwner();
    fMCList->Add(fCPAPlots);
    TString MCPricpaPtBinName = Form("CPAPtBinningPri");
    TString MCMatcpaPtBinName = Form("CPAPtBinningMat");
    TString MCSeccpaPtBinName = Form("CPAPtBinningSec");
    TString MCConcpaPtBinName = Form("CPAPtBinningCont");

    fMCPrimCPAPtBins = new TH2F(
        MCPricpaPtBinName.Data(),MCPricpaPtBinName.Data(),
        8,0.3,4.3,1000,0.90,1.);
    fMCPrimCPAPtBins->Sumw2();
    fMCPrimCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCPrimCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCPrimCPAPtBins);

    fMCMaterialCPAPtBins = new TH2F(
        MCMatcpaPtBinName.Data(),MCMatcpaPtBinName.Data(),
        8,0.3,4.3,1000,0.90,1.);
    fMCMaterialCPAPtBins->Sumw2();
    fMCMaterialCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCMaterialCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCMaterialCPAPtBins);

    fMCSecondaryCPAPtBins = new TH2F(
        MCSeccpaPtBinName.Data(),MCSeccpaPtBinName.Data(),
        8,0.3,4.3,1000,0.90,1.);
    fMCSecondaryCPAPtBins->Sumw2();
    fMCSecondaryCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCSecondaryCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCSecondaryCPAPtBins);

    fMCContCPAPtBins = new TH2F(
        MCConcpaPtBinName.Data(),MCConcpaPtBinName.Data(),
        8,0.3,4.3,1000,0.90,1.);
    fMCContCPAPtBins->Sumw2();
    fMCContCPAPtBins->GetXaxis()->SetTitle("P_{T}");
    fMCContCPAPtBins->GetYaxis()->SetTitle("CPA");
    fCPAPlots->Add(fMCContCPAPtBins);
  } else {
    fCPAPlots=0;
    fMCPrimCPAPtBins=0;
    fMCMaterialCPAPtBins=0;
    fMCSecondaryCPAPtBins=0;
    fMCContCPAPtBins=0;
  }
}

AliFemtoDreamv0MCHist::~AliFemtoDreamv0MCHist() {
  if (fMCList) {
    delete fMCList;
  }
}

void AliFemtoDreamv0MCHist::FillMCCPAPtBins(
    AliFemtoDreamBasePart::PartOrigin org,double pT,double cpa) {
  if (org==AliFemtoDreamBasePart::kPhysPrimary) {
    fMCPrimCPAPtBins->Fill(pT,cpa);
  } else if(org == AliFemtoDreamBasePart::kWeak) {
    fMCSecondaryCPAPtBins->Fill(pT,cpa);
  } else if(org == AliFemtoDreamBasePart::kMaterial) {
    fMCMaterialCPAPtBins->Fill(pT,cpa);
  } else if(org == AliFemtoDreamBasePart::kContamination) {
    fMCContCPAPtBins->Fill(pT,cpa);
  } else {
    AliFatal("Particle Origin not implemented");
  }
}
