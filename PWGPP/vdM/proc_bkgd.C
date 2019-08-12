// -*- C++ -*-

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TString.h>

#include "AliTriggerBCMask.h"

#include "AliVdMPileup.h"
#include "proc_bkgd.h"

template<typename T>
T* SetAttr(T* h, Int_t color) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(kFullDotMedium);
  h->SetLineWidth(1);
  return h;
}
TH1* FitPol0(TH1* h) {
  h->Fit("pol0");
  dynamic_cast<TF1*>(h->FindObject("pol0"))->SetLineColor(h->GetLineColor());
  return h;
}

void proc_bkgd(const AliVdMMetaData& vdmMetaData,
               AliVdMScanData& allData,
               const char* className)
{
  TFile::Open(AliVdMMetaData::GetFileName(TString::Format("%d/%d.root",
                                                          vdmMetaData.GetFillNumber(),
                                                          vdmMetaData.GetFillNumber())));
  TTree *tOffline = nullptr;
  gFile->GetObject("Vertex_Performance/cOutputVtxESD", tOffline);

  const std::vector<std::string> neededBranches({"timestamp", "bx", "isT0", "timeAD", "timeV0"});
  tOffline->SetBranchStatus("*", 0);
  for (const auto& branchName : neededBranches)
    tOffline->SetBranchStatus(branchName.c_str(), 1);

  struct TreeData {
    UInt_t timeStamp;
    UInt_t bx;
    Bool_t isT0;
    Float_t timeV0[2];
    Float_t timeAD[2];
    Short_t decV0online[2];
    Short_t decADonline[2];
  } treeData;
  tOffline->SetBranchAddress("timeStamp", &treeData.timeStamp);
  tOffline->SetBranchAddress("bx",        &treeData.bx);
  tOffline->SetBranchAddress("isT0",      &treeData.isT0);
  tOffline->SetBranchAddress("timeV0",    &treeData.timeStamp);
  tOffline->SetBranchAddress("timeAD",    &treeData.timeAD);

  typedef std::map<Short_t, TGraphErrors> map_t; // BCID -> TGraphErrors
  map_t gAnotC, gCnotA;  // one-arm/two-arm ratios

  // (1) fill one-arm/two-arm ratio graphs for all BCIDs
  for (Int_t iScan=0; iScan<4; ++iScan) {
    AliVdMTree& vt  = allData.GetMap(iScan)[className];
    vtAND.Zip3([&gAnotC,&gCnotA](const AliVdMTree::DefaultBranchData& d,
                                 AliVdMTree::branchMap_t& mapAC,
                                 AliVdMTree::branchMap_t& mapAnotC,
                                 AliVdMTree::branchMap_t& mapCnotA)
               {
                 AliVdMTree::ValErr v1 = mapAnotC["rate"];
                 v1 /= mapAC["rate"];
                 // v1 /= mapAC["relBkgd"];
                 v1 *= mapAnotC["relBkgd"];
                 if (!v1.isInf() && v1.val()) {
                   const Int_t m1 = gAnotC[d.BCID()].GetN();
                   gAnotC[d.BCID()].SetPoint     (m1, mapAC["mu"].val(), v1.val());
                   gAnotC[d.BCID()].SetPointError(m1, mapAC["mu"].err(), v1.err());
                 }
                 AliVdMTree::ValErr v2 = mapCnotA["rate"];
                 v2 /= mapAC["rate"];
                 // v2 /= mapAC["relBkgd"];
                 v2 *= mapCnotA["relBkgd"];
                 if (!v2.isInf() && v2.val()) {
                   const Int_t m2 = gCnotA[d.BCID()].GetN();
                   gCnotA[d.BCID()].SetPoint     (m2, mapAC["mu"].val(), v2.val());
                   gCnotA[d.BCID()].SetPointError(m2, mapAC["mu"].err(), v2.err());
                 }
               },
               vtAnotC, vtCnotA);
  }

  // (2) fit model
  AliVdMPileup pileupModel;

  TString pn = TString::Format("pileup_fill%d_%s.pdf",
                               vdmMetaData.GetFillNumber(),
                               classAC);
  TCanvas *c1 = new TCanvas;
  c1->SaveAs(pn+"[");

  const AliTriggerBCMask& bcMask = vdmMetaData.GetTriggerBCMask();
  const Int_t             nBCs   = bcMask.GetNUnmaskedBCs();

  const TString histTitleBase = TString::Format("fill %d;BCID;", vdmMetaData.GetFillNumber());
  TH1 *hPar[5] = {
    SetAttr(new TH1D("hrA",      histTitleBase+"r_{A}",           nBCs,0,nBCs), kRed),
    SetAttr(new TH1D("hrC",      histTitleBase+"r_{C}",           nBCs,0,nBCs), kBlue),
    SetAttr(new TH1D("hbkgdA",   histTitleBase+"bkgd_{A}",        nBCs,0,nBCs), kRed),
    SetAttr(new TH1D("hbkgdC",   histTitleBase+"bkgd_{C}",        nBCs,0,nBCs), kBlue),
    SetAttr(new TH1D("hChi2NDF", histTitleBase+"#chi^{2}/n.d.f.", nBCs,0,nBCs), kBlue)
  };

  for (Int_t bc=0, counter=0; bc<3564; ++bc) {
    if (bcMask.GetMask(bc))
      continue;

    const TString binLabel = TString::Format("%d", bc);
    for (Int_t i=0; i<5; ++i)
      hPar[i]->GetXaxis()->SetBinLabel(1+counter, binLabel);

    c1->Clear();
    c1->SetLogx();
    c1->SetLogy();
    TH1 *hf = c1->DrawFrame(1e-6, 0.01, 0.5, 20);
    hf->SetTitle(TString::Format("fill %d BCID=%d %s;#mu_{0} [rate(%s,BC=%d)/11245Hz];one-arm/two-arm",
                                 vdmMetaData.GetFillNumber(), bc, classAC, classAC, bc));
    hf->GetXaxis()->SetTitleOffset(1.3);
    pileupModel.DoFit(&gAnotC[bc], &gCnotA[bc], &par0[0]);

    SetAttr(&gAnotC[bc], kBlue);
    SetAttr(&gCnotA[bc], kRed);
    gAnotC[bc].Draw("PE");
    gCnotA[bc].Draw("PE");

    TF1 *fAnotC = SetAttr(new TF1("fAnotC", &pileupModel, &AliVdMPileup::fcnAnotC, 1e-6, 0.5, 5), kBlue);
    fAnotC->SetParameters(pileupModel.GetPar());
    fAnotC->SetNpx(1000);
    fAnotC->Draw("same");

    TF1 *fCnotA = SetAttr(new TF1("fCnotA", &pileupModel, &AliVdMPileup::fcnCnotA, 1e-6, 0.5, 5), kRed);
    fCnotA->SetParameters(pileupModel.GetPar());
    fCnotA->SetNpx(1000);
    fCnotA->Draw("same");

    TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg->AddEntry(&gCnotA[bc], "AnotC/AandC", "PEL");
    leg->AddEntry(&gAnotC[bc], "CnotA/AandC", "PEL");
    leg->Draw();

    TPaveText *pt = new TPaveText(0.6, 0.4, 0.9, 0.7, "NDC NB");
    pt->SetFillStyle(0);
    pt->AddText(TString::Format("#chi^{2}/n.d.f = %.0f/%.0f = %.2f", pileupModel.GetChi2(), pileupModel.GetNDF(), pileupModel.GetChi2()/pileupModel.GetNDF()));
    {
      double curval,err, lowlim, uplim;
      int iuint;
      TString name;
      for (Int_t ivar=0; ivar<4; ++ivar) {
        gMinuit->mnpout(ivar, name, curval, err, lowlim, uplim,iuint);
        hPar[ivar]->SetBinContent(1+counter, curval);
        hPar[ivar]->SetBinError(1+counter, err);
        if (ivar==0) {
          hf->SetMinimum(0.5*curval);
        }
        if (ivar==1) {
          hf->SetMinimum(TMath::Min(hf->GetMinimum(), 0.5*curval));
        }
        if (ivar < 2)
          pt->AddText(TString::Format("%s = %.4f#pm%.4f", name.Data(), curval, err));
        else
          pt->AddText(TString::Format("%s = %.1e#pm%.1e", name.Data(), curval, err));
        pt->GetLine(1+ivar)->SetTextColor(ivar%2 ? kRed : kBlue);
      }
      hPar[4]->SetBinContent(1+counter, pileupModel.GetChi2()/pileupModel.GetNDF());
    }
    pt->Draw();
    c1->SaveAs(pn);
    Printf("%f / %f", pileupModel.GetChi2(), pileupModel.GetNDF());
    ++counter;
  }
  gStyle->SetOptStat("n");
  gStyle->SetOptFit(111);
  TCanvas *c2 = new TCanvas;
  for (Int_t i=0; i<4; ++i) {
    FitPol0(hPar[i])->Draw();
    c2->SaveAs(pn);
  }
  hPar[4]->SetMinimum(0);
  hPar[4]->Draw();
  c2->SaveAs(pn+")");
}

