// -*- C++ -*-

#include <vector>
#include <string>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <TCanvas.h>

#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"

#include "AliVdMPileup.h"

#include "proc_pileup.h"

// -----------------------------------------------------------------------------
//  analysis for fill 4937

AliVdMScanData allData;

void proc()
{
  const std::vector<std::string> triggerNames{
    "c2TVX",
    "c2VBAandVBC",
    "c2UBAandUBC",

    "c2T0AandNotT0C",
    "c2T0CandNotT0A",

    "c2VBAandNotVBC",
    "c2VBCandNotVBA",

    "c2UBAandNotUBC",
    "c2UBCandNotUBA"
};

  const AliVdMMetaData vdmMetaData(AliVdMMetaData::GetFileName("4937/4937.xml"));

  TFile::Open(AliVdMMetaData::GetFileName("4937/vdm_time_4937_6m11_12p17_1_v3.root"), "READ");
  TTree * VdM = nullptr;
  gFile->GetObject("VdM", VdM);
  VdM->AddFriend("DDL2", AliVdMMetaData::GetFileName("4937/vdm_DDL2_4937.root"));

  allData.FillDefaultBranches(vdmMetaData, VdM, triggerNames);

  Printf("#scans: %ld", std::distance(vdmMetaData.GetScansBegin(), vdmMetaData.GetScansEnd()));

  // (1) compute uncorrected rates
  for (Int_t iScan=0; iScan<6; ++iScan) {
    for (const std::string& triggerName : triggerNames) {
      AliVdMTree& vt = allData.GetMap(iScan)[triggerName];
      vt.AddBranch("rate", [](const AliVdMTree::DefaultBranchData& d,
                              AliVdMTree::branchMap_t& map) {
                     map["rate"].val() = d.Counter() / d.DeltaT();
                     map["rate"].err() = TMath::Sqrt(d.Counter()) / d.DeltaT();
                   });
    }
  }
  // (2) compute mu for (T0,V0,AD)_and triggers
  for (Int_t iScan=0; iScan<6; ++iScan) {
    for (const std::string& triggerName : triggerNames) {
      if (TString(triggerName.c_str()).Contains("Not")) // skip one-arm triggers
        continue;
      AliVdMTree& vt = allData.GetMap(iScan)[triggerName];
      vt.AddBranch("mu", [](const AliVdMTree::DefaultBranchData& d,
                            AliVdMTree::branchMap_t& map) {
                     map["mu"].val() = map["rate"].val() / 11245.0;
                     map["mu"].err() = map["rate"].err() / 11245.0;
                   });
    }
  }

  // (3) compute bkgd due to from tails from previous interactions
  for (Int_t iScan=0; iScan<6; ++iScan) {
    for (const std::string& triggerName : triggerNames) {
      AliVdMTree& vt = allData.GetMap(iScan)[triggerName];
      vt.AddBranch("relBkgd", [](const AliVdMTree::DefaultBranchData& d,
                                 AliVdMTree::branchMap_t& map)
                   {
                     const Double_t *y = d.Counters().GetMatrixArray();
                     const Int_t     n = d.Counters().GetNoElements();
                     TVectorD x(n), ey(n);
                     for (Int_t l=0; l<n; ++l) {
                       x[l]  = l;
                       ey[l] = TMath::Sqrt(d.Counter(l));
                     }
                     TGraphErrors tmp(n, x.GetMatrixArray(), y, nullptr, ey.GetMatrixArray());
                     TF1 f1("f1", "[0]*exp([1]*x)");
                     f1.SetParameters(TMath::Mean(n-5, y+5), 1e-3);
                     f1.SetParLimits(0,0,1000);
                     f1.SetParLimits(1,0,0.2);
                     tmp.Fit(&f1, "Q", "", 2.5,n+0.5);
                     if (f1.GetParameter(0)>10 && kFALSE) {
                       tmp.Draw("APL");
                       char _c;
                       gets(&_c);
                       gPad->Update();
                     }
                     map["relBkgd"].val() = (d.Counter(0) ? (d.Counter(0)-f1.Eval(0))/d.Counter(0) : 1.0);
                     map["relBkgd"].err() = (d.Counter(0) && f1.Eval(0)
                                             ? f1.GetParError(0)/d.Counter(0)
                                             : 0.0);
                   });
    }
  }

  // (3) determine pile-up (work in progress)
  proc_pileup(vdmMetaData, allData,
              "c2UBAandUBC", "c2UBAandNotUBC", "c2UBCandNotUBA",
              {0.24,0.13, 4.4e-5,2.3e-5});
  proc_pileup(vdmMetaData, allData,
              "c2VBAandVBC", "c2VBAandNotVBC", "c2VBCandNotVBA",
              {0.061,0.076, 14e-5,3.5e-5});
  proc_pileup(vdmMetaData, allData,
              "c2TVX", "c2T0AandNotT0C", "c2T0CandNotT0A",
              {0.39,0.44, 13e-5,3.3e-5});
}

