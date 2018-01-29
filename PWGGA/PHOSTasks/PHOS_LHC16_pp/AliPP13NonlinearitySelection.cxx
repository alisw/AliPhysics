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
void AliPP13NonlinearitySelection::FillPi0Mass(TObjArray * clusArray, TList * pool, const EventFlags & eflags)
{
	(void) pool;
	// Ensure that we are not doing mixing
	EventFlags flags = eflags;
	flags.isMixing = kFALSE;

	// Select photons
	TObjArray photonCandidates;
	SelectPhotonCandidates(clusArray, &photonCandidates, flags);

	// Consider N^2 - N combinations, excluding only same-same clusters.
	for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); ++i)
	{
		AliVCluster * first = dynamic_cast<AliVCluster *> (photonCandidates.At(i));

		for (Int_t j = 0; j < photonCandidates.GetEntriesFast(); ++j)
		{
			if (j == i) // Skip the same clusters
				continue;

			AliVCluster * second = dynamic_cast<AliVCluster *> (photonCandidates.At(j));
			ConsiderPair(first, second, flags);
		} // second cluster loop
	} // cluster loop

	MixPhotons(photonCandidates, pool, flags);
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


	if (dynamic_cast<AliPP13SelectionWeights *>(fWeights))
	{
		Float_t eff = fWeights->Weight(p1.E()) * fWeights->Weight(p2.E());
		fMassPt[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, p1.Pt(), 1. / eff);	
		return;
	}

	Float_t weight = fWeights->Weight(pt12);
	// NB: Weight by meson spectrum, but fill only for the first photon
	// fMassPt[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, pt12, weight);
	
	fMassPt[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, p1.Pt(), weight);	
}


//________________________________________________________________
void AliPP13NonlinearitySelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 750;
	Double_t mMin  = 0.0;
	Double_t mMax  = 1.5;
	Int_t nPt      = 400;
	Double_t ptMin = 0;
	Double_t ptMax = 20;


	for (Int_t i = 0; i < 2; ++i)
	{
		const char * sf = (i == 0) ? "" : "Mix";
		TH2F * hist = new TH2F(Form("h%sMassPt_", sf), "(M_{#gamma#gamma}, pT_{#gamma}) ; M_{#gamma#gamma}, GeV; pT_{#gamma}, GeV", nM, mMin, mMax, nPt, ptMin, ptMax);
		fMassPt[i] = new AliPP13DetectorHistogram(hist, fListOfHistos, AliPP13DetectorHistogram::kModules);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}

}
