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
 * @relates AliTrackletDeltaWeights
 * @ingroup pwglf_forward_tracklets
 */
struct MakeDeltaWeights : public AliTrackletAODUtils
{
  MakeDeltaWeights() {}
  /** 
   * Create @f$\Delta@f$ dependent weights.  
   * 
   * If @a dimen is set to 3, then the weights depend on the tracklet 
   * @f$(\Delta,\eta)@f$ and event @f$\mathrm{IP}_z@f$. 
   *
   * If @a dimen is set to 2, then the weights depend on the tracklet 
   * @f$(\Delta,\eta)@f$ - this is the default. 
   *
   * If @a dimen is set to 1, then the weights depend on the tracklet 
   * @f$\Delta@f$ only
   * 
   * @param realFileName File with real data (to correct to)
   * @param simFileName  File with simulated data (to correct)
   * @param dimen        Number of dimensions 
   * @param scaleToTail  If we should scale the simulated distribution to 
   *                     the tail of the real distribution 
   * 
   * @return The weight structure on success, null otherwise 
   */
  AliTrackletDeltaWeights* Run(const char* realFileName,
			       const char* simFileName,
			       Int_t       dimen=2,
			       Bool_t      scaleToTail=true)
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
      
      ProcessBin(i, dimen, scaleToTail, c1, c2, realTop, simTop, weights);
    }

    weights->Write();
    out->Write();
    Printf("Output writtten to %s", outName.Data());

    return weights;
  }
  /** 
   * @deprecated 
   * 
   * @param h Histogram to normalize 
   */
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
   * Calculate 1-D weights i.e., @f$\Delta@f$ dependent
   * 
   * @param realMeas     Real container 
   * @param simMeas      Simulation container 
   * @param scaleToTail  If true, scale simulation to real tail first 
   * 
   * @return Weights histogram
   */
  TH1* Weights1(Container* realMeas,
		Container* simMeas,
		Bool_t     scaleToTail)
  {
    TH1*      realDelta = GetH1(realMeas, "delta");
    TH1*      simDelta  = GetH1(simMeas,  "delta");
    if (scaleToTail) {
      Double_t realTail      = GetD(realMeas,  "deltaTail");
      Double_t realTailE     = GetD(realMeas,  "deltaTailError");
      Double_t simTail       = GetD(simMeas,   "deltaTail");
      Double_t simTailE      = GetD(simMeas,   "deltaTailError");
      Double_t scaleE, scale = RatioE(realTail, realTailE,
				      simTail,  simTailE,
				      scaleE);
      Scale(simDelta, scale, scaleE);
    }

    TH1*       ratio     = static_cast<TH1*>(realDelta->Clone("tmp"));
    ratio->Divide(simDelta);
    ratio->Smooth(3);
    
    return ratio;
  }
  /** 
   * Calculate 2-D weights i.e., @f$\Delta,\eta@f$ dependent
   * 
   * @param realMeas     Real container 
   * @param simMeas      Simulation container 
   * @param scaleToTail  If true, scale simulation to real tail first 
   * 
   * @return Weights histogram
   */
  TH2* Weights2(Container* realMeas,
		Container* simMeas,
		Bool_t     scaleToTail)
  {
    TH2*  realDelta = GetH2(realMeas, "etaDelta");
    TH2*  simDelta  = GetH2(simMeas,  "etaDelta");
    if (scaleToTail) {
      TH1*  realTail  = GetH1(realMeas, "etaDeltaTail");
      TH1*  simTail   = GetH1(simMeas,  "etaDeltaTail");
      TH1*  scale     = static_cast<TH1*>(realTail->Clone("scale"));
      scale->Divide(simTail);
      scale->SetDirectory(0);
      Scale(simDelta, scale);
      delete scale;
    }

    TH2*       ratio     = static_cast<TH2*>(realDelta->Clone("tmp"));
    ratio->Divide(simDelta);
    ratio->Smooth(1);
    
    return ratio;
  }
  /** 
   * Calculate 3-D weights i.e., @f$\Delta,\eta,\mathrm{IP}_{z}@f$ dependent
   * 
   * @param realMeas     Real container 
   * @param simMeas      Simulation container 
   * @param scaleToTail  If true, scale simulation to real tail first 
   * 
   * @return Weights histogram
   */
  TH3* Weights3(Container* realMeas,
		Container* simMeas,
		Bool_t     scaleToTail)
  {
    TH3*  realDelta = GetH3(realMeas, "etaDeltaIPz");
    TH3*  simDelta  = GetH3(simMeas,  "etaDeltaIPz");
    if (scaleToTail) {
      TH2*  realTail  = GetH2(realMeas, "etaIPzDeltaTail");
      TH2*  simTail   = GetH2(simMeas,  "etaIPzDeltaTail");
      TH2*  scale     = static_cast<TH2*>(realTail->Clone("scale"));
      scale->Divide(simTail);
      scale->SetDirectory(0);
      ScaleDelta(simDelta, scale);
      delete scale;
    }

    TH3*       ratio     = static_cast<TH3*>(realDelta->Clone("tmp"));
    ratio->Divide(simDelta);
    
    return ratio;
  }
  /** 
   * Process a centrality bin 
   * 
   * @param bin         Bin number 
   * @param dimen       Dimensions
   * @param scaleToTail If true, scale simulation to real tail first 
   * @param c1          Lowest centrality 
   * @param c2          Highest centrality 
   * @param realTop     Top container of the real data 
   * @param simTop      Top container of the simulated data
   * @param weights     The weight structure to fill 
   * 
   * @return true on success. 
   */
  Bool_t ProcessBin(Int_t                    bin,
		    Int_t                    dimen,
		    Bool_t                   scaleToTail,
		    Double_t                 c1,
		    Double_t                 c2,
		    Container*               realTop,
		    Container*               simTop,
		    AliTrackletDeltaWeights* weights)
  {
    Printf(" %5.2f-%5.2f%%", c1, c2);
    TString    name      = CentName(c1, c2);
    Color_t    color     = CentColor(bin);
    Container* realBin   = GetC(realTop,   name);
    Container* simBin    = GetC(simTop,    name);
    if (!realBin || !simBin) return false;
    Container* realMeas  = GetC(realBin,   "measured");
    Container* simMeas   = GetC(simBin,    "measured");
    if (!realMeas || !simMeas) return false;
    TH1*       ratio     = 0;
    switch (dimen) {
    case 2:
      ratio = Weights2(realMeas, simMeas, scaleToTail);
      break;
    case 3:
      ratio = Weights3(realMeas, simMeas, scaleToTail);
      break;
    default:
      ratio = Weights1(realMeas, simMeas, scaleToTail);
      break;
    }  
    if (!ratio) return false;

    ratio->SetName(name);
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

