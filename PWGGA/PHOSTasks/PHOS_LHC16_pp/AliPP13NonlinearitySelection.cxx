// --- Custom header files ---
#include "AliPP13NonlinearitySelection.h"
#include "AliPP13DetectorHistogram.h"

// --- ROOT system ---
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13NonlinearitySelection);

//________________________________________________________________
void AliPP13NonlinearitySelection::SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & eflags)
{
	// NB: Nonlinearity is a function of photon energy
	//     therefore the histograms should be filled for each photon.

	// Int_t counter = 0;
	// Consider N^2 - N combinations, excluding only same-same clusters.
	for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); ++i)
	{
		AliVCluster * first = dynamic_cast<AliVCluster *> (photonCandidates.At(i));

		for (Int_t j = 0; j < photonCandidates.GetEntriesFast(); ++j)
		{
			if (j == i) // Skip the same clusters
				continue;

			AliVCluster * second = dynamic_cast<AliVCluster *> (photonCandidates.At(j));

			if (!fCuts.AcceptPair(first, second, eflags))
				continue;

			ConsiderPair(first, second, eflags);
		} // second cluster loop
	} // cluster loop
	// Int_t Nn = photonCandidates.GetEntriesFast();
}

//________________________________________________________________
void AliPP13NonlinearitySelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1, p2, psum;
	c1->GetMomentum(p1, eflags.vtxBest);
	c2->GetMomentum(p2, eflags.vtxBest);
	psum = p1 + p2;
	Float_t pt12 = psum.Pt();
	Float_t m12 = psum.M();


	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok


	// NB: This is the data cut
	if (dynamic_cast<AliPP13SelectionWeightsTOF *>(fWeights))
	{
		Float_t eff = fWeights->TofEfficiency(p1.E()) * fWeights->TofEfficiency(p2.E());
		fMassPt[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, p1.E(), 1. / eff);	
		return;
	}

	// NB: Weight by meson spectrum, but fill only for the first photon
	Float_t weight = fWeights->Weights(pt12, eflags);
	fMassPt[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, p1.E(), weight);	
}

//________________________________________________________________
void AliPP13NonlinearitySelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = fLimits.nM;
	Double_t mMin  = fLimits.mMin;
	Double_t mMax  = fLimits.mMax;
	Int_t nPt      = fLimits.nPt;
	Double_t ptMin = fLimits.ptMin;
	Double_t ptMax = fLimits.ptMax;


	for (Int_t i = 0; i < 2; ++i)
	{
		const char * sf = (i == 0) ? "" : "Mix";
		TH2F * hist = new TH2F(Form("h%sMassPt_", sf), "(M_{#gamma#gamma}, pT_{#gamma}) ; M_{#gamma#gamma} (GeV/#it{c}^{2}); E_{#gamma}, GeV", nM, mMin, mMax, nPt, ptMin, ptMax);
		fMassPt[i] = new AliPP13DetectorHistogram(hist, fListOfHistos, AliPP13DetectorHistogram::kModules);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}

}
