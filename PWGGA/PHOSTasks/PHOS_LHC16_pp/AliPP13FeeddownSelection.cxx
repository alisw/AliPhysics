
// #include "iterator"

// --- Custom header files ---
#include "AliPP13FeeddownSelection.h"

// --- ROOT system ---
#include <TParticle.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliVCluster.h>
#include <AliAnalysisManager.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13FeeddownSelection);


//________________________________________________________________
void AliPP13FeeddownSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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

	fInvMass[eflags.isMixing]->Fill(ma12, pt12);

	if (eflags.isMixing)
		return;

	AliAODMCParticle * gamma1_mother = GetMother(c1, eflags.fMcParticles);
	AliAODMCParticle * gamma2_mother = GetMother(c2, eflags.fMcParticles);

	if (!gamma1_mother || !gamma2_mother)
		return;

	if (gamma1_mother != gamma2_mother)
		return;

	if (gamma1_mother->GetPdgCode() != kPi0)
		return;

	// Looking at the source of pi0
	//

	AliAODMCParticle * hadron = GetMother(gamma1_mother, eflags.fMcParticles);
	if (!hadron)
		return;

	if (hadron->GetPdgCode() != kK0s)
		return;

	Double_t weight = fWeights->Weights(hadron->Pt(), eflags);
	dynamic_cast<TH2F *>(fFeedownK0s[0])->Fill(ma12, pt12, weight);
}


//________________________________________________________________
void AliPP13FeeddownSelection::InitSelectionHistograms()
{
	Int_t nM       = fLimits.nM;
	Double_t mMin  = fLimits.mMin;
	Double_t mMax  = fLimits.mMax;
	Int_t nPt      = fLimits.nPt;
	Double_t ptMin = fLimits.ptMin;
	Double_t ptMax = fLimits.ptMax;

	for (Int_t i = 0; i < 2; ++i)
	{
		fInvMass[i] = new TH2F(Form("h%sMassPt", i == 0 ? "" : "Mix") , "(M,p_{T})_{#gamma#gamma}, N_{cell}>2; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax);
		fListOfHistos->Add(fInvMass[i]);
	}
	fFeedownK0s[0] = new TH2F("hMassPt_#pi^{0}_feeddown_K^{s}_{0}", "(M,p_{T})_{#gamma#gamma} originating form K^{s}_{0}; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax);
	fListOfHistos->Add(fFeedownK0s[0]);

	fFeedownK0s[1] = new TH1F("hMassPt_#pi^{0}_feeddown_K^{s}_{0}_generated", "(M,p_{T})_{#gamma#gamma} originating form K^{s}_{0}; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nPt, ptMin, ptMax);
	fListOfHistos->Add(fFeedownK0s[1]);

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}

//________________________________________________________________
void AliPP13FeeddownSelection::ConsiderGeneratedParticles(const EventFlags & eflags)
{
	if (!eflags.fMcParticles)
		return;

	for (Int_t i = 0; i < eflags.fMcParticles->GetEntriesFast(); i++)
	{
		AliAODMCParticle * particle = (AliAODMCParticle *) eflags.fMcParticles->At(i);
		Int_t code = TMath::Abs(particle->GetPdgCode());

		// NB: replace this condition by find, if the number of particles will grow
		//
		if (code != kPi0)
			continue;

		Double_t pt = particle->Pt();

		// Use this to remove forward photons that can modify our true efficiency
		if (TMath::Abs(particle->Y()) > 0.5) // NB: Use rapidity instead of pseudo rapidity!
			continue;

		AliAODMCParticle * hadron = GetMother(particle, eflags.fMcParticles);

		if (!hadron)
			return;

		if (hadron->GetPdgCode() == kK0s)
			fFeedownK0s[1]->Fill(pt, fWeights->Weights(pt, eflags));
	}
}

//________________________________________________________________
AliAODMCParticle * AliPP13FeeddownSelection::GetMother(const AliAODMCParticle * particle, TClonesArray * particles) const
{
	Int_t plabel = particle->GetMother();

	if (plabel <= -1)
		return 0;

	AliAODMCParticle * parent = dynamic_cast<AliAODMCParticle * >(particles->At(plabel));
	return parent;
}
