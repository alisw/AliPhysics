// --- Custom header files ---
#include "AliPP13TriggerEfficiency.h"
#include "AliPP13DetectorHistogram.h"
#include <AliPP13AnalysisCluster.h>

// --- ROOT system ---
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13TriggerEfficiency);

//________________________________________________________________
void AliPP13TriggerEfficiency::SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & eflags)
{
	// NB: Trigger efficiency is a function of a photon registration efficiency
	//     therefore the histograms should be filled for each photon.

	// Consider N^2 - N combinations, excluding only same-same clusters.
	for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); i++)
	{
		AliPP13AnalysisCluster * tag = dynamic_cast<AliPP13AnalysisCluster *> (photonCandidates.At(i));

		// TODO: Implement the trigger cluster selection
		if (!tag->IsTrigger())
			continue;

		for (Int_t j = 0; j < photonCandidates.GetEntriesFast(); j++)
		{
			if (i == j) // Skip the same clusters
				continue;

			AliPP13AnalysisCluster * probe = dynamic_cast<AliPP13AnalysisCluster *> (photonCandidates.At(j));

			if (!fCuts.AcceptPair(tag, probe, eflags))
				continue;

			ConsiderPair(tag, probe, eflags);
		} // second cluster loop
	} // cluster loop}
}


//________________________________________________________________
void AliPP13TriggerEfficiency::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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


	const AliPP13AnalysisCluster * tag = dynamic_cast<const AliPP13AnalysisCluster * >(c1);
	const AliPP13AnalysisCluster * probe = dynamic_cast<const AliPP13AnalysisCluster * >(c2);

	// Filter out clusters from the same 4x4 patch
	//

	Bool_t close = TMath::Abs(x1 - x2) < 4 && TMath::Abs(z1 - z2) < 4;
	if ((probe->TRU() == tag->TRU()) &&  (sm1 == sm2) && close)
		return;

	// Int_t tru = probe->TRU();
	Int_t mix = int(eflags.isMixing);

	fTotalMassEnergyAll[mix]->FillAll(sm1, sm1, m12, energy);
	// fMassEnergyAll[tru][mix]->FillAll(sm1, sm1, m12, energy);

	if (!probe->IsTrigger())
		return;

	fTotalMassEnergyTrigger[mix]->FillAll(sm1, sm1, m12, energy);
	// fMassEnergyTrigger[tru][mix]->FillAll(sm1, sm1, m12, energy);
}


//________________________________________________________________
void AliPP13TriggerEfficiency::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 250;
	Double_t mMin  = 0.0;
	Double_t mMax  = 0.3;

	Int_t nE      = 2000;
	Double_t eMin = 0;
	Double_t eMax = 20;

	fNevents = new TH1F("hNevents", "Number of events", 2, 0.5, 2.5);
	fListOfHistos->Add(fNevents);

	for (Int_t i = 0; i < 2; ++i)
	{
		const char * sf = (i == 0) ? "" : "Mix";
		const char * title = "(M_{#gamma#gamma}, E_{probe}); M_{#gamma#gamma} (GeV/#it{c}^{2}); E_{probe}, GeV";

		TH1 * hist1 = new TH2F(Form("h%sMassEnergyAll_", sf), title, nM, mMin, mMax, nE, eMin, eMax);
		TH1 * hist2 = new TH2F(Form("h%sMassEnergyTrg_", sf), title, nM, mMin, mMax, nE, eMin, eMax);

		fTotalMassEnergyAll[i] = new AliPP13DetectorHistogram(hist1, fListOfHistos);
		fTotalMassEnergyTrigger[i] = new AliPP13DetectorHistogram(hist2, fListOfHistos);
	}


	// for (Int_t tru = 0; tru < kTRUs; ++tru)
	// {
	// 	for (Int_t i = 0; i < 2; ++i)
	// 	{
	// 		const char * sf = (i == 0) ? "" : "Mix";
	// 		const char * title = Form("(M_{#gamma#gamma}, E_{probe}) TRU #%d ; M_{#gamma#gamma} (GeV/#it{c}^{2}); E_{probe}, GeV", tru);
	// 		TH2F * hist1 = new TH2F(Form("h%sMassEnergyAll_TRU_%d_", sf, tru), title, nM, mMin, mMax, nE, eMin, eMax);
	// 		TH2F * hist2 = new TH2F(Form("h%sMassEnergyTrg_TRU_%d_", sf, tru), title, nM, mMin, mMax, nE, eMin, eMax);

	// 		fMassEnergyAll[tru][i] = new AliPP13DetectorHistogram(hist1, fListOfHistos, AliPP13DetectorHistogram::kSingleHist);
	// 		fMassEnergyTrigger[tru][i] = new AliPP13DetectorHistogram(hist2, fListOfHistos, AliPP13DetectorHistogram::kSingleHist);
	// 	}
	// }

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}

}

void AliPP13TriggerEfficiency::SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags)
{
	AliPP13PhysicsSelection::SelectPhotonCandidates(clusArray, candidates, eflags);
	fNevents->Fill(1);

	if (eflags.fTriggerEvent)
		fNevents->Fill(2);
}
