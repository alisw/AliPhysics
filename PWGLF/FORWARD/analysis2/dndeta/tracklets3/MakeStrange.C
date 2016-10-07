/**
 * @file   MakeStrange.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:14:22 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
/** 
 * @{
 * @name Make strangeness weights 
 *
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * Make strangeness weights
 * 
 * @relates AliTrackletBaseWeights
 */
void
MakeStrange()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  AliTrackletPtPidStrWeights* weights = new AliTrackletPtPidStrWeights("weights");

  Double_t c2Bins[] = { 0, 2.5, 5, 7.5, 10, 20, 30,
			40, 50, 60, 70, 80, 90, 100 };

  TFile* pidFile = TFile::Open("../tracklets/REWEIGHTpid.root","READ");
  Double_t c2Bins[] = { 0, 2.5, 5, 7.5, 10, 20, 30,
			40, 50, 60, 70, 80, 90, 100 };
  Short_t pids[] = { 321 };
  for (Int_t i = 0; i < 1; i++) {
    TString histName;
    Short_t pdg = pids[i];
    histName.Form("pidWeight_%s", (pdg == 211 ? "pi" :
				   pdg == 321 ? "ka" :
				   "pr"));
    TH1* hist = static_cast<TH1*>(pidFile->Get(histName));
    if (!hist) {
      Warning("MakeWeight", "pid histogram %s not found", histName.Data());
      continue;
    }
    TH1D* out = new TH1D(Form("w%d", pdg), Form("PID %d weight", pdg),
			 13, c2Bins);
    for (Int_t j = 1; j <= out->GetNbinsX(); j++) {
      out->SetBinContent(j, hist->GetBinContent(j));
      out->SetBinError  (j, hist->GetBinError  (j));
    }
    weights->AddAbundanceWeight(pdg, out);
  }
  
  TFile* strFile = TFile::Open("../tracklets/REWEIGHTstr.root","READ");
  Short_t strs[] = { 321, 310, 3122, 3112, 3222, 3312 /* 3212, 3322 */};
  for (Int_t i = 0; i < 6; i++) {
    TString histName;
    Short_t pdg = strs[i];
    histName.Form("strWeight_%s", (pdg == 211  ? "pi" :
				   pdg == 321  ? "ka" :
				   pdg == 2212 ? "pr" :
				   pdg == 310  ? "k0" :
				   pdg == 3122 ? "la" :
				   pdg == 3112 ? "si" :
				   pdg == 3222 ? "si" : 
				   pdg == 3212 ? "si" :
				   pdg == 3312 ? "xi" : 
				   pdg == 3322 ? "xi" : ""));
    TH1* hist = static_cast<TH1*>(strFile->Get(histName));
    if (!hist) {
      Warning("MakeWeight", "strangeness histogram %s not found",
	      histName.Data());
      continue;
    }
    TH1D* out = new TH1D(Form("w%d", pdg), Form("Strangeness %d weight", pdg),
			 13, c2Bins);
    for (Int_t j = 1; j <= out->GetNbinsX(); j++) {
      out->SetBinContent(j, hist->GetBinContent(j));
      out->SetBinError  (j, hist->GetBinError  (j));
    }
    weights->AddStrangenessWeight(pdg, out);
  }
  

  
  
  new TBrowser;

  TFile* outW = TFile::Open("str.root", "RECREATE");
  weights->Write();
  outW->Write();
  weights->Draw();
  weights->Print();
}
/* @} */
//
// EOF
//
