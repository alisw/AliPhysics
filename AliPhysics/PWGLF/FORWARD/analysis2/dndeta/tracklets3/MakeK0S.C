/**
 * @file   MakeK0S.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:12:09 2016
 * 
 * @brief  Make K0short weights
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */

/** 
 * @{ 
 * @name Name @f$ K_S^0@f$ weights
 *
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * Make @f$ K^0_S@f$ weights
 * 
 * @relates AliTrackletBaseWeights
 */
void
MakeK0S()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  AliTrackletPtPidStrWeights* weights = new AliTrackletPtPidStrWeights("weights");

  Double_t c2Bins[] = { 0, 2.5, 5, 7.5, 10, 20, 30,
			40, 50, 60, 70, 80, 90, 100 };
  
  TFile* strFile = TFile::Open("../tracklets/REWEIGHTstr.root","READ");
  Short_t strs[] = { 310 };
  for (Int_t i = 0; i < 1; i++) {
    TString histName;
    Short_t pdg = strs[i];
    histName.Form("strWeight_%s", (pdg == 211  ? "pi" :
				   pdg == 321  ? "ka" :
				   pdg == 2212 ? "pr" :
				   pdg == 310  ? "k0" :
				   pdg == 3122 ? "la" :
				   pdg == 3212 ? "si" :
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

  TFile* outW = TFile::Open("k0s.root", "RECREATE");
  weights->Write();
  outW->Write();
  weights->Draw();
  weights->Print();
}
/* @} */
//
// EOF
//
