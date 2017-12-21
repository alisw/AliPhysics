// --- Custom header files ---
#include "AliPP13TagAndProbeSelection.h"
#include "AliPP13DetectorHistogram.h"

// --- ROOT system ---
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13TagAndProbeSelection);


//________________________________________________________________
void AliPP13TagAndProbeSelection::FillPi0Mass(TObjArray * clusArray, TList * pool, const EventFlags & eflags)
{
	(void) pool;
	// Ensure that we are not doing mixing
	EventFlags flags = eflags;
	flags.isMixing = kFALSE;

	// Select photons
	TObjArray photonCandidates;
	SelectPhotonCandidates(clusArray, &photonCandidates, flags);

	// Consider N^2 - N combinations, excluding only same-same clusters.
	for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); i++)
	{
		AliVCluster * tag = dynamic_cast<AliVCluster *> (photonCandidates.At(i));

		if (TMath::Abs(tag->GetTOF()) > fTimingCut)
			continue;

		for (Int_t j = 0; j < photonCandidates.GetEntriesFast(); j++)
		{
			if (i == j) // Skip the same clusters
				continue;

			AliVCluster * proble = dynamic_cast<AliVCluster *> (photonCandidates.At(j));
			ConsiderPair(tag, proble, flags);
		} // second cluster loop
	} // cluster loop

	MixPhotons(photonCandidates, pool, flags);
}

//________________________________________________________________
void AliPP13TagAndProbeSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1, p2, psum;
	c1->GetMomentum(p1, eflags.vtxBest);
	c2->GetMomentum(p2, eflags.vtxBest);
	psum = p1 + p2;
	Float_t energy = p2.E();
	Float_t m12 = psum.M();


	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	fMassEnergyAll[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, energy);

	if (TMath::Abs(c2->GetTOF()) > fTimingCut)
		return;

	fMassEnergyTOF[int(eflags.isMixing)]->FillAll(sm1, sm2, m12, energy);
}


//________________________________________________________________
void AliPP13TagAndProbeSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 250;
	Double_t mMin  = 0.0;
	Double_t mMax  = 0.3;

	Int_t nE      = 2000;
	Double_t eMin = 0;
	Double_t eMax = 20;


	for (Int_t i = 0; i < 2; ++i)
	{
		const char * sf = (i == 0) ? "" : "Mix";
		TH2F * hist1 = new TH2F(Form("h%sMassEnergyAll_", sf), "(M_{#gamma#gamma}, E_{probe}) ; M_{#gamma#gamma}, GeV; E_{proble}, GeV", nM, mMin, mMax, nE, eMin, eMax);
		TH2F * hist2 = new TH2F(Form("h%sMassEnergyTOF_", sf), "(M_{#gamma#gamma}, E_{probe}) ; M_{#gamma#gamma}, GeV; E_{proble}, GeV", nM, mMin, mMax, nE, eMin, eMax);

		fMassEnergyAll[i] = new AliPP13DetectorHistogram(hist1, fListOfHistos, AliPP13DetectorHistogram::kModules);
		fMassEnergyTOF[i] = new AliPP13DetectorHistogram(hist2, fListOfHistos, AliPP13DetectorHistogram::kModules);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}

}