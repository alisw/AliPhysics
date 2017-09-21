// --- Custom header files ---
#include "NonlinearityScanSelection.h"

// --- AliRoot header files ---
#include <AliPHOSAodCluster.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(NonlinearityScanSelection);


//________________________________________________________________
TLorentzVector NonlinearityScanSelection::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags, Int_t ia, Int_t ib) const
{
	Float_t energy = c1->E();

	TLorentzVector p;
	c1->GetMomentum(p, eflags.vtxBest);
	p *= Nonlinearity(energy, ia, ib);
	return p;
}

//________________________________________________________________
Float_t NonlinearityScanSelection::Nonlinearity(Float_t x, Int_t ia, Int_t ib) const
{
	Float_t non_a = GetA(ia);
	Float_t non_sigma = GetSigma(ib);

	return fGlobalEnergyScale * (1. + non_a * TMath::Exp(-x / 2. * x / non_sigma / non_sigma));
}

//________________________________________________________________
void NonlinearityScanSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 750;
	Double_t mMin  = 0.0;
	Double_t mMax  = 1.5;
	Int_t nPt      = 400;
	Double_t ptMin = 0;
	Double_t ptMax = 20;

	for (Int_t ia = 0; ia < kNbinsA; ++ia)
	{
		for (Int_t ib = 0; ib < kNbinsSigma; ++ib)
		{
			Float_t a = GetA(ia);
			Float_t b = GetSigma(ib);

			fInvariantMass[ia][ib]    = new TH2F(Form("hMassPt_%d_%d", ia, ib), Form("%f %f; M_{#gamma#gamma}, GeV; p_{T}, GeV/c", a, b), nM, mMin, mMax, nPt, ptMin, ptMax);
			fMixInvariantMass[ia][ib] = new TH2F(Form("hMixMassPt_%d_%d", ia, ib), Form("%f %f; M_{#gamma#gamma}, GeV; p_{T}, GeV/c", a, b), nM, mMin, mMax, nPt, ptMin, ptMax);

			fListOfHistos->Add(fInvariantMass[ia][ib]);
			fListOfHistos->Add(fMixInvariantMass[ia][ib]);
		}
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
}


void NonlinearityScanSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	for (Int_t ia = 0; ia < kNbinsA; ++ia)
	{
		for (Int_t ib = 0; ib < kNbinsSigma; ++ib)
		{
			TLorentzVector p1 = ClusterMomentum(c1, eflags, ia, ib);
			TLorentzVector p2 = ClusterMomentum(c2, eflags, ia, ib);
			TLorentzVector psum = p1 + p2;

			if (psum.M2() < 0)
				return;

			Double_t ma12 = psum.M();
			Double_t pt12 = psum.Pt();
			TH1 * hist = (!eflags.isMixing) ? fInvariantMass[ia][ib] : fMixInvariantMass[ia][ib];

			hist->Fill(ma12, pt12);
		}
	}
}

