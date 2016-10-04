/**
 * @file   MakeEposWeight.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:11:23 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
/** 
 * @{ 
 * @name Calculate epos-lhc weights 
 *
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * Get a profile from file 
 * 
 * @param filename File to get it from 
 * 
 * @return Pointer or null
 */
TProfile* GetOne(const char* filename)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) return 0;

  Bool_t mc  = !TString(filename).Contains("real");  
  TList* top = static_cast<TList*>(file->Get(Form("MidRapidity%sResults",
						   mc ? "MC" : "")));
  if (!top) {
    Error("GetOne", "Couldn't get result container in %s", filename);
    file->Close();
    return 0;
  }

  TProfile* prf = static_cast<TProfile*>(top->FindObject("centTracklets"));
  if (!prf) 
    Error("GetOne", "Couldn't get profile from %s", filename);

  file->Close();
  return prf;
}

/** 
 * Make EPOS weights 
 * 
 * @param file1 
 * @param file2 
 * @relates AliTrackletBaseWeights
 */
void
MakeEposWeight(const char* file1, const char* file2)
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include "
			  "-I${ALICE_PHYSICS}/include "
			  "-I${ANA_SRC}/dndeta/traclets3");
  gROOT->SetMacroPath(Form("%s:$ANA_SRC/dndeta/tracklets3",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  AliTrackletPtPidStrWeights* weights = new AliTrackletPtPidStrWeights("weights");

  // --- pT weight ---------------------------------------------------
  // Int_t    nBins  = 10;
  // Double_t bins[] = { 0,    5,   10, 20, 30, 40, 50, 60, 70, 80, 90 };
  // Double_t facs[] = { 1.12, 1.01,1.,.99,.98,.98,.98,.97,.95,.88 };  
  // TH2D* cPt = new TH2D("cPt", "cPt", nBins, bins, 1, 0, 100);
  // for (Int_t i = 1; i <= nBins; i++) 
  //   cPt->SetBinContent(i, 1, facs[i-1]);
  TProfile* p1 = GetOne(file1);
  TProfile* p2 = GetOne(file2);
  if (!p1 || !p2) return;

  TH1::SetDefaultSumw2();
  TProfile* r  = static_cast<TProfile*>(p1->Clone("ratio"));
  r->SetErrorOption(" ");
  r->Divide(p2);
  TF1* f = new TF1("rf", "pol6", 0, 90);
  r->Fit(f, "Q0R+", "");

  TH2D* cPT = new TH2D("cPt", "cPt", 100, 0, 100, 1, 0, 100);
  for (Int_t i = 1; i <= cPt->GetNbinsX(); i++) 
    cPt->SetBinContent(i, 1, f->Eval(i-.5));
  
  TCanvas* cp = new TCanvas("p","p");
  cp->Divide(2,1);
  TVirtualPad* pp = cp->cd(1);
  p1->SetTitle(file1);
  p2->SetTitle(file2);
  p1->SetMarkerColor(kGreen+2); p2->SetMarkerColor(kRed+2);  
  p1->SetLineColor  (kGreen+2); p2->SetLineColor  (kRed+2);
  p1->Draw();
  p2->Draw("same");
  pp->BuildLegend();
  pp = cp->cd(2);
  // r->Draw("i");
  cPt->ProjectionX()->Draw("");
  f->Draw("same");
  cp->Modified();
  cp->Update();
  cp->cd();
  
  weights->SetPtWeight(cPt);

  // --- Abundance weights -------------------------------------------
  // We do not add any!

  // --- Strangeness weights -----------------------------------------
  // We do not add any!

  new TBrowser;
  
  // --- Write to file -----------------------------------------------
  TFile* outW = TFile::Open("epos.root", "RECREATE");  
  weights->Write();
  outW->Write();


  TCanvas* c = new TCanvas("w","w");
  weights->Draw();
  weights->Print();
}
/* @} */
//
// EOF
//
