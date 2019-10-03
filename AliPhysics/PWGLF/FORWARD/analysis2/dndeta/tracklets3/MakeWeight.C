/**
 * @file   MakeWeight.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:15:00 2016
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * @{ 
 * @name Make Roberto's weights 
 *
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * Extract weights from Roberto's files 
 * 
 * @relates AliTrackletBaseWeights
 */
void
MakeWeight()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  AliTrackletPtPidStrWeights* weights = new AliTrackletPtPidStrWeights("weights");

  TFile* ptFile = TFile::Open("../tracklets/REWEIGHTpt.root","READ");

  Double_t c1Bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80 };
  TH2D*    cPt      = new TH2D("cPt", "cPt", 9, c1Bins, 
			       1000, 0, 5);
  for (Int_t i = 1; i <= cPt->GetNbinsX(); i++) {
    TString histName;
    histName.Form("ptWeight_c%d_%d",
		  Int_t(cPt->GetXaxis()->GetBinLowEdge(i)),
		  Int_t(cPt->GetXaxis()->GetBinUpEdge(i)));
    TH1* hist = static_cast<TH1*>(ptFile->Get(histName));
    if (!hist) {
      Warning("MakeWeight", "pt histogram %s not found", histName);
      continue;
    }
    for (Int_t j = 1; j <= cPt->GetNbinsY(); j++) {
      cPt->SetBinContent(i, j, hist->GetBinContent(j));
      cPt->SetBinError  (i, j, hist->GetBinError  (j));
    }
  }
  weights->SetPtWeight(cPt);
      
  TFile* pidFile = TFile::Open("../tracklets/REWEIGHTpid.root","READ");
  Double_t c2Bins[] = { 0, 2.5, 5, 7.5, 10, 20, 30,
			40, 50, 60, 70, 80, 90, 100 };
  Short_t pids[] = { 211, 321, 2212 };
  for (Int_t i = 0; i < 3; i++) {
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
  Short_t strs[] = { 211, 321, 2212, 310, 3122, 3112, 3222, 3312 /*3212,3322*/};
  for (Int_t i = 0; i < 8; i++) {
    TString histName;
    Short_t pdg = strs[i];
    histName.Form("strWeight_%s", (pdg == 211  ? "pi" :
				   pdg == 321  ? "ka" :
				   pdg == 2212 ? "pr" :
				   pdg == 310  ? "k0" :
				   pdg == 3122 ? "la" :
				   pdg == 3222 ? "si" :
				   pdg == 3112 ? "si" :
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

  TFile* outW = TFile::Open("weights.root", "RECREATE");
  weights->Write();
  outW->Write();
  weights->Draw();
  weights->Print();
}
/* @} */
//
// EOF
// 
