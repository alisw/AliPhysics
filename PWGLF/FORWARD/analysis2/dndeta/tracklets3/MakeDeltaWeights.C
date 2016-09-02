/**
 * @file   MakeDeltaWeights.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 31 10:56:38 2016
 * 
 * @brief  Make weight based on delta, eta, ipz
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#include <TSystem.h>
#include "AliTrackletAODUtils.C"
#include "AliTrackletWeights.C"

/**
 * Structure to make weights for delta 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct MakeDeltaWeights : public AliTrackletAODUtils
{
  MakeDeltaWeights() {}

  AliTrackletDeltaWeights*Run(const char* realFileName,
			      const char* simFileName)
  {
    TFile* realFile = 0;
    TFile* simFile  = 0;
    if (!(realFile = OpenFile(realFileName))) return 0;
    if (!(simFile  = OpenFile(simFileName)))  return 0;

    TString     outName(gSystem->BaseName(simFileName));    
    const char* base    = "MidRapidity";
    Container*  realTop = GetC(realFile, Form("%sResults",   base));
    Container*  simTop  = GetC(simFile,  Form("%sMCResults", base));
    if (!realTop) {
      realTop = GetC(realFile, Form("%sMCResults",   base));
      if (realTop)
	Warning("Run","\n"
		"*********************************************\n"
		"* Testing MC vs. MC:                        *\n"
		"*  'Data' file:      %23s *\n"
		"*  Simulation file:  %23s *\n"
		"*********************************************\n",
		realFileName, simFileName);
      // fRealIsSim = true;
      outName.Append("_");
      outName.Append(gSystem->BaseName(realFileName));
      outName.ReplaceAll(".root","");
    }
    outName.ReplaceAll(".root","");
    outName.Prepend("delta_");
    outName.Append(".root");
    
    TH1* realCent = GetH1(realTop, "cent");
    TH1* simCent  = GetH1(simTop,  "cent");
    TH1* realIPz  = GetH1(realTop, "ipz");
    TH1* simIPz   = GetH1(simTop,  "ipz");

    // Check consistency of found histograms 
    if (!CheckConsistency(realCent, simCent)) {
      Warning("Post", "Centrality bins are incompatible, giving up");
      return 0;
    }
    if (!CheckConsistency(realIPz, simIPz)) {
      Warning("Post", "IPz bins are incompatible, giving up");
      return 0;
    }

    TFile* out = TFile::Open(outName, "RECREATE");
    out->cd();

    AliTrackletDeltaWeights* weights =
      new AliTrackletDeltaWeights("weights", "w_{#Delta}");
    weights->SetCentAxis(*(realCent->GetXaxis()));
    // Loop over defined centrality bins 
    for (Int_t i = 1; i <= realCent->GetNbinsX(); i++) {
      Double_t c1 = realCent->GetXaxis()->GetBinLowEdge(i);
      Double_t c2 = realCent->GetXaxis()->GetBinUpEdge (i);
      
      ProcessBin(i, c1, c2, realTop, simTop, weights);
    }

    weights->Write();
    out->Write();
    Printf("Output writtten to %s", outName.Data());

    return weights;
  }
  void Normalize(TH3* h)
  {
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      for (Int_t j = 1; j <= h->GetNbinsZ(); j++) {
	// Project out delta distribution for this bin
	TH1*     tmp  = h->ProjectionY(Form("%s_%02d_%02d",
					    h->GetName(), i, j),
				       i, i, j, j);
	Double_t intg = tmp->Integral();
	tmp->SetDirectory(0);
	// Scale to average 
	if (intg > 1e-3) tmp->Scale(tmp->GetNbinsX()*1./intg);
	else             tmp->Reset();
	for (Int_t k = 1; k <= h->GetNbinsY(); k++) {
	  h->SetBinContent(i, k, j, tmp->GetBinContent(k));
	  h->SetBinError  (i, k, j, tmp->GetBinError  (k));
	}
	delete tmp;
      }
    }
  }
  /** 
   * Process a centrality bin 
   * 
   * @param bin      Bin number 
   * @param c1       Lowest centrality 
   * @param c2       Highest centrality 
   * @param realTop  Top container of the real data 
   * @param simTop   Top container of the simulated data
   * @param weights  The weight structure to fill 
   * 
   * @return true on success. 
   */
  Bool_t ProcessBin(Int_t                    bin,
		    Double_t                 c1,
		    Double_t                 c2,
		    Container*               realTop,
		    Container*               simTop,
		    AliTrackletDeltaWeights* weights)
  {
    TString    name      = CentName(c1, c2);
    Color_t    color     = CentColor(bin);
    Container* realBin   = GetC(realTop,   name);
    Container* simBin    = GetC(simTop,    name);
    if (!realBin || !simBin) return false;
    Container* realMeas  = GetC(realBin,   "measured");
    Container* simMeas   = GetC(simBin,    "measured");
    if (!realMeas || !simMeas) return false;
    TH3*       realDelta = GetH3(realMeas, "etaDeltaIPz");
    TH3*       simDelta  = GetH3(simMeas,  "etaDeltaIPz");
    if (!realDelta || !simDelta) return false;
    // Normalize(realDelta);
    // Normalize(simDelta);
    
    TH3*       ratio     = static_cast<TH3*>(realDelta->Clone(name));
    ratio->Divide(simDelta);
    ratio->SetDirectory(0);    
    ratio->SetMarkerColor(color);
    ratio->SetLineColor(color);
    ratio->SetFillColor(color);
    // Normalize(ratio);
    
    weights->SetHisto(bin, ratio);
    return true;
  }
};

//
// EOF
//

