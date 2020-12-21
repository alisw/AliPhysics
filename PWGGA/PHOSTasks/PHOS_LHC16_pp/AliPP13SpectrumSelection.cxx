// --- Custom header files ---
#include "AliPP13SpectrumSelection.h"

// --- ROOT system ---
#include <TH2F.h>
#include <TH3F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13SpectrumSelection);


//________________________________________________________________
void AliPP13SpectrumSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	//

	for (Int_t i = 0; i < 2; ++i)
	{
		const char * s = (i == 0) ? "" : "Mix";
		TH1 * hist = new TH2F(
			Form("h%sMassPt", s),
			"(M,p_{T})_{#gamma#gamma}, ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})",
			fLimits.nM, fLimits.mMin, fLimits.mMax,
			fLimits.nPt, fLimits.ptMin, fLimits.ptMax
		);
		fInvariantMass[i] = new AliPP13DetectorHistogram(hist, fListOfHistos);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}

	// These histograms are needed only to check the performance
	// Don't do any analysis with these histograms.
	//

	fClusters = new TH1F("hClusterPt_SM0", "Cluster p_{T} spectrum with default cuts, all modules; p_{T} (GeV/#it{c})", fLimits.nPt, fLimits.ptMin, fLimits.ptMax);	
	fListOfHistos->Add(fClusters);
}

//________________________________________________________________
void AliPP13SpectrumSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1 = ClusterMomentum(c1, eflags);
	TLorentzVector p2 = ClusterMomentum(c2, eflags);
	TLorentzVector psum = p1 + p2;

	// Pair cuts can be applied here
	if (psum.M2() < 0)  return;

	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	Double_t ma12 = psum.M();
	Double_t pt12 = psum.Pt();

	Float_t eff = fWeights->TofEfficiency(p1.E()) * fWeights->TofEfficiency(p2.E());
	fInvariantMass[Int_t(eflags.isMixing)]->FillAll(sm1, sm2, ma12, pt12, 1. / eff);
}

//________________________________________________________________
void AliPP13SpectrumSelection::FillClusterHistograms(const AliVCluster * clus, const EventFlags & eflags)
{
	TLorentzVector p = ClusterMomentum(clus, eflags);
	
	if(fClusters)
		fClusters->Fill(p.Pt());
}
