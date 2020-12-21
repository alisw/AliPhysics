// -*- C++ -*-

#include <TPad.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "proc.h"

void proc(const AliVdMMetaData& vdmMetaData,
          AliVdMScanData& vdmScanData,
          const std::vector<std::string>& triggerNames,
          Bool_t computeBkgd)
{

  const Int_t nScans = vdmScanData.GetNScans();
  Printf("nScans = %d", nScans);

  // (1) compute uncorrected rates
  for (Int_t iScan=0; iScan<nScans; ++iScan) {
    Printf("%d", iScan);
    for (const std::string& triggerName : triggerNames) {
      AliVdMTree& vt = vdmScanData.GetMap(iScan)[triggerName];
      vt.AddBranch("rate", [](const AliVdMTree::DefaultBranchData& d,
                              AliVdMTree::branchMap_t& map) {
                     map["rate"].val() = d.Counter() / d.DeltaT();
                     map["rate"].err() = TMath::Sqrt(d.Counter()) / d.DeltaT();
                     Printf("rate: %f +- %f (%.0f %f)", map["rate"].val(), map["rate"].err(), d.Counter(), d.DeltaT());
                   });
    }
  }

  // (2) compute mu for (T0,V0,AD)_and triggers
  for (Int_t iScan=0; iScan<nScans; ++iScan) {
    for (const std::string& triggerName : triggerNames) {
      if (TString(triggerName.c_str()).Contains("Not") ||
          TString(triggerName.c_str()).Contains("NOT")) // skip one-arm triggers
        continue;
      AliVdMTree& vt = vdmScanData.GetMap(iScan)[triggerName];
      vt.AddBranch("mu", [](const AliVdMTree::DefaultBranchData& d,
                            AliVdMTree::branchMap_t& map) {
                     map["mu"].val() = map["rate"].val() / 11245.0;
                     map["mu"].err() = map["rate"].err() / 11245.0;
                     Printf("mu=%f", map["mu"].val());
                   });
    }
  }

  // (3) compute bkgd due to from tails from previous interactions
  if (computeBkgd) {
    for (Int_t iScan=0; iScan<nScans; ++iScan) {
      for (const std::string& triggerName : triggerNames) {
        AliVdMTree& vt = vdmScanData.GetMap(iScan)[triggerName];
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
  }
}

