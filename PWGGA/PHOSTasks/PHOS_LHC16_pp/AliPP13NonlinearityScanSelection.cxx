// --- Custom header files ---
#include "AliPP13NonlinearityScanSelection.h"

// --- AliRoot header files ---
#include <AliPHOSAodCluster.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13NonlinearityScanSelection);


//________________________________________________________________
TLorentzVector AliPP13NonlinearityScanSelection::ClusterMomentumBinned(const AliVCluster * c1, const EventFlags & eflags, Int_t ia, Int_t ib) const
{
	TLorentzVector p;
	c1->GetMomentum(p, eflags.vtxBest);

	Float_t energy = c1->E();
	p *= fWeightsScan[ia][ib].Nonlinearity(energy);
	return p;
}

//________________________________________________________________
void AliPP13NonlinearityScanSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = fLimits.nM;
	Double_t mMin  = fLimits.mMin;
	Double_t mMax  = fLimits.mMax;
	Int_t nPt      = fLimits.nPt;
	Double_t ptMin = fLimits.ptMin;
	Double_t ptMax = fLimits.ptMax;


	for (Int_t ia = 0; ia < kNbinsA; ++ia)
	{
		for (Int_t ib = 0; ib < kNbinsB; ++ib)
		{
			Float_t a = fWeightsScan[ia][ib].fE;
			Float_t b = fWeightsScan[ia][ib].fD;

			fInvariantMass[ia][ib] = new TH2F(Form("hMassPt_%d_%d", ia, ib), Form("%f %f; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", a, b), nM, mMin, mMax, nPt, ptMin, ptMax);
			fMixInvariantMass[ia][ib] = new TH2F(Form("hMixMassPt_%d_%d", ia, ib), Form("%f %f; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", a, b), 10, mMin, mMax, 10, ptMin, ptMax);

			fListOfHistos->Add(fInvariantMass[ia][ib]);
			fListOfHistos->Add(fMixInvariantMass[ia][ib]);
		}
	}
	fPtPrimaryPi0 = new TH1F(
	    "hPt_primary_#pi^{0}_",
	    "Generated p_{T} spectrum of primary #pi^{0}s; p_{T} (GeV/#it{c})",
	    nPt, ptMin, ptMax);


	// NB: Reduce the selection size
	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}


// NB: We need scan to test all possible nonlinearities
//________________________________________________________________
void AliPP13NonlinearityScanSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	for (Int_t ia = 0; ia < kNbinsA; ++ia)
	{
		for (Int_t ib = 0; ib < kNbinsB; ++ib)
		{
			TLorentzVector p1 = ClusterMomentumBinned(c1, eflags, ia, ib);
			TLorentzVector p2 = ClusterMomentumBinned(c2, eflags, ia, ib);
			TLorentzVector psum = p1 + p2;

			if (psum.M2() < 0)
				return;

			Double_t m12 = psum.M();
			Double_t pt12 = psum.Pt();
			TH2 * hist = dynamic_cast<TH2 *> ((!eflags.isMixing) ? fInvariantMass[ia][ib] : fMixInvariantMass[ia][ib]);

			Float_t weight = fWeights->Weights(pt12, eflags);
			hist->Fill(m12, pt12, weight);
		}
	}
}

//________________________________________________________________
void AliPP13NonlinearityScanSelection::ConsiderGeneratedParticles(const EventFlags & eflags)
{
	if (!eflags.fMcParticles)
		return;

	for (Int_t i = 0; i < eflags.fMcParticles->GetEntriesFast(); i++)
	{
		AliAODMCParticle * particle = ( AliAODMCParticle *) eflags.fMcParticles->At(i);
		Int_t code = TMath::Abs(particle->GetPdgCode());

		// NB: replace this condition by find, if the number of particles will grow
		//
		if (code != kPi0) // Only neutral pions
			continue;


		Double_t pt = particle->Pt();
		Double_t w = fWeights->Weights(pt, eflags);

		// Use this to remove forward photons that can modify our true efficiency
		if (TMath::Abs(particle->Y()) > 0.5) // NB: Use rapidity instead of pseudo rapidity!
			continue;

		Bool_t primary = IsPrimary(particle);
		if (!primary)
			continue;
		fPtPrimaryPi0->Fill(pt, w);
	}
}
