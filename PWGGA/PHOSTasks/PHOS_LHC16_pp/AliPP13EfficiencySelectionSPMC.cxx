
// #include "iterator"

// --- Custom header files ---
#include "AliPP13EfficiencySelectionSPMC.h"

// --- ROOT system ---
#include <TParticle.h>
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13EfficiencySelectionSPMC);


//________________________________________________________________
void AliPP13EfficiencySelectionSPMC::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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

    // AliAODMCParticle * origin = (AliAODMCParticle*)eflags.fMcParticles->At(0);//0 is always generated particle by AliGenBox.
	// Double_t w = fWeights->Weight(origin->Pt());
	Double_t w = fWeights->Weights(pt12, eflags);
	TH2 * hist = dynamic_cast<TH2 *> (fInvMass[eflags.isMixing]);
	hist->Fill(ma12, pt12, w);

	if (eflags.isMixing)
		return;
}

void AliPP13EfficiencySelectionSPMC::ConsiderGeneratedParticles(const EventFlags & eflags)
{
	if (!eflags.fMcParticles)
		return;

	for (Int_t i = 0; i < eflags.fMcParticles->GetEntriesFast(); i++)
	{
		AliAODMCParticle * particle = ( AliAODMCParticle *) eflags.fMcParticles->At(i);
		Int_t code = TMath::Abs(particle->GetPdgCode());

		// NB: replace this condition by find, if the number of particles will grow
		//
		if (code != kGamma && code != kPi0 && code != kEta)
			continue;


		Double_t pt = particle->Pt();

	    // AliAODMCParticle * origin = (AliAODMCParticle*)eflags.fMcParticles->At(0);//0 is always generated particle by AliGenBox.
		// Double_t w = fWeights->Weight(origin->Pt());
		Double_t w = fWeights->Weights(pt, eflags);


		// Use this to remove forward photons that can modify our true efficiency
		if (TMath::Abs(particle->Y()) > 0.5) // NB: Use rapidity instead of pseudo rapidity!
			continue;

		Double_t r = TMath::Sqrt(particle->Xv() * particle->Xv() + particle->Yv() * particle->Yv());

		fSpectrums[code]->fPt->Fill(pt, w);
		fSpectrums[code]->fPtRadius->Fill(pt, r, w);

		Bool_t primary = IsPrimary(particle);


		// Tese conditions are just for QA purpose
		if (primary && particle->E() > 0.3)
		{
			fSpectrums[code]->fPtLong->Fill(pt, w);
			fSpectrums[code]->fPtAllRange->Fill(pt, w);
			fSpectrums[code]->fEtaPhi->Fill(particle->Phi(), particle->Y());
		}

		fSpectrums[code]->fPtPrimaries[Int_t(primary)]->Fill(pt, w);
		fSpectrums[code]->fPtPrimariesStandard[Int_t(primary)]->Fill(pt, w);
		ConsiderGeneratedParticle(i, pt, primary, eflags);
	}
}