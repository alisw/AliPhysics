// --- Custom header files ---
#include "AliPP13SpectrumSelectionSimple.h"

// --- ROOT system ---
#include <TH2F.h>
#include <TH3F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13SpectrumSelectionSimple);


//________________________________________________________________
void AliPP13SpectrumSelectionSimple::InitSelectionHistograms()
{
	fMassPt = new TH2F(
		"hMassPtSM0",
		"(M,p_{T})_{#gamma#gamma}, ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})",
		fLimits.nM, fLimits.mMin, fLimits.mMax,
		fLimits.nPt, fLimits.ptMin, fLimits.ptMax
	);
	fMassPt->Sumw2();


	fMixMassPt = new TH2F(
		"hMixMassPtSM0",
		"(M,p_{T})_{#gamma#gamma}, ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})",
		fLimits.nM, fLimits.mMin, fLimits.mMax,
		fLimits.nPt, fLimits.ptMin, fLimits.ptMax
	);
	fMixMassPt->Sumw2();

	fListOfHistos->Add(fMassPt);
	fListOfHistos->Add(fMixMassPt);
}

//________________________________________________________________
void AliPP13SpectrumSelectionSimple::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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
	if(!eflags.isMixing)
	{
		fMassPt->Fill(ma12, pt12, 1. / eff);
	}
	else
	{
		fMixMassPt->Fill(ma12, pt12, 1. / eff);
	}

}
