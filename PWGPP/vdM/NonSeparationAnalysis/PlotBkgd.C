// -*- C++ -*-
// $Id$

#include <sstream>

#include <TH1.h>
#include <TCut.h>
#include <TH2.h>
#include <TTree.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TText.h>

#include "AliVdMData.h"


void PlotBkgd(TCanvas *c1, Int_t fillNumber, TTree *TE, const AliXMLEngine::Node& n, const TCut& vtxCuts) {
  TTree tSep;
  std::istringstream iss(n.GetData());
  tSep.ReadStream(iss);
  Int_t m = tSep.Draw("timeStart:timeEnd:sep", "", "GOFF");
  const Double_t *timeStart = tSep.GetV1();
  const Double_t *timeEnd   = tSep.GetV2();

  Printf("%s %d %.0f %.0f",
         AliVdMData::GetScanName(n).Data(),
         AliVdMData::GetScanType(n),
         timeStart[0], timeEnd[m-1]);

  TCut timeCut(TString::Format("timestamp>%.0f && timestamp<%.0f", timeStart[0], timeEnd[m-1]));

  c1->cd();
  c1->Clear();
  TPad *p1 = new TPad("p1", "", 0,0, 1,0.95);
  p1->Draw();
  p1->Divide(2,2);
  c1->cd();

  TText tt;
  tt.SetTextSize(0.05);
  tt.SetTextAlign(21);
  tt.DrawTextNDC(0.5, 0.96, TString::Format("Fill %d ", fillNumber)+AliVdMData::GetScanName(n));

  TString binsX = TString::Format("%.0f,%.0f,%.0f", timeEnd[m-1]-timeStart[0], timeStart[0], timeEnd[m-1]);

  // TE->SetBranchStatus("*", 0);
  // TE->SetBranchStatus("timestamp", 1);
  // TE->SetBranchStatus("timeAD.*", 1);
  // TE->SetBranchStatus("timeV0.*", 1);

  p1->cd(1);
  TE->Draw("timeAD.A:timestamp>>hAD_A("+binsX+",128,50.75,63.25)", timeCut*vtxCuts, "GOFF");
  TE->GetHistogram()->SetTitle("ADA;;time ADA (ns)");
  TE->GetHistogram()->SetStats(0);
  TE->GetHistogram()->GetXaxis()->SetTimeDisplay(1);
  TE->GetHistogram()->GetYaxis()->SetRangeUser(0,-1);
  TE->GetHistogram()->Draw("COLZ");
  p1->cd(2);
  TE->Draw("timeAD.C:timestamp>>hAD_C("+binsX+",128,58.75,71.25)", timeCut*vtxCuts, "GOFF");
  TE->GetHistogram()->SetTitle("ADC;;time ADC (ns)");
  TE->GetHistogram()->SetStats(0);
  TE->GetHistogram()->GetXaxis()->SetTimeDisplay(1);
  TE->GetHistogram()->GetYaxis()->SetRangeUser(0,-1);
  TE->GetHistogram()->Draw("COLZ");
  p1->cd(3);
  TE->Draw("timeV0.A:timestamp>>hV0_A("+binsX+",128,4.75,17.25)", timeCut*vtxCuts, "GOFF");
  TE->GetHistogram()->SetTitle("V0A;;time V0A (ns)");
  TE->GetHistogram()->SetStats(0);
  TE->GetHistogram()->GetXaxis()->SetTimeDisplay(1);
  TE->GetHistogram()->GetYaxis()->SetRangeUser(0,-1);
  TE->GetHistogram()->Draw("COLZ");
  p1->cd(4);
  TE->Draw("timeV0.C:timestamp>>hV0_C("+binsX+",128,-3.25,9.25)", timeCut*vtxCuts, "GOFF");
  TE->GetHistogram()->SetTitle("V0C;;time V0C (ns)");
  TE->GetHistogram()->SetStats(0);
  TE->GetHistogram()->GetXaxis()->SetTimeDisplay(1);
  TE->GetHistogram()->GetYaxis()->SetRangeUser(0,-1);
  TE->GetHistogram()->Draw("COLZ");
}
