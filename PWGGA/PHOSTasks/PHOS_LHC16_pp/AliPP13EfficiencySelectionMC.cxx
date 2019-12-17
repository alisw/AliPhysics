// --- Custom header files ---
#include "AliPP13EfficiencySelectionMC.h"

// --- ROOT system ---
#include <TParticle.h>
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13EfficiencySelectionMC);


//________________________________________________________________
void AliPP13EfficiencySelectionMC::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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
	Double_t w = fWeights->Weights(pt12, eflags);

	TH2 * hist = dynamic_cast<TH2 *> (fInvMass[eflags.isMixing]);
	hist->Fill(ma12, pt12, w);
}



//________________________________________________________________
void AliPP13EfficiencySelectionMC::InitSelectionHistograms()
{
	Int_t nM       = fLimits.nM;
	Double_t mMin  = fLimits.mMin;
	Double_t mMax  = fLimits.mMax;
	Int_t nPt      = fLimits.nPt;
	Double_t ptMin = fLimits.ptMin;
	Double_t ptMax = fLimits.ptMax;

	const char * rtitle = "(M,p_{T})_{#gamma#gamma}, N_{cell}>2; M_{#gamma#gamma}";
	for (Int_t i = 0; i < 2; ++i)
	{
		const char * rname = Form("h%sMassPt", i == 0 ? "" : "Mix");
		fInvMass[i] = new TH2F(rname, rtitle, nM, mMin, mMax, nPt, ptMin, ptMax);
		fListOfHistos->Add(fInvMass[i]);
	}

	Float_t ptbins[] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0};
	Int_t ptsize = sizeof(ptbins) / sizeof(Float_t);

	for (EnumNames::iterator i = fPartNames.begin(); i != fPartNames.end(); ++i)
	{
		const char * n = (const char *) i->second.Data();
		fSpectrums[i->first] = new ParticleSpectrum(n, fListOfHistos, ptsize - 1, ptbins);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}

//________________________________________________________________
void AliPP13EfficiencySelectionMC::ConsiderGeneratedParticles(const EventFlags & eflags)
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
	}
}
