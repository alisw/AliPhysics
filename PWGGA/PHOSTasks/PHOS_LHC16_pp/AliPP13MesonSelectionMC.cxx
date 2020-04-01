
// #include "iterator"

// --- Custom header files ---
#include "AliPP13MesonSelectionMC.h"

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


ClassImp(AliPP13MesonSelectionMC);


//________________________________________________________________
void AliPP13MesonSelectionMC::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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

	Int_t label1 = c1->GetLabelAt(0) ;
	Int_t label2 = c2->GetLabelAt(0) ;

	AliAODMCParticle * mother1 = GetParent(label1, eflags.fMcParticles);
	AliAODMCParticle * mother2 = GetParent(label2, eflags.fMcParticles);

	if (!mother1 || !mother2)
		return;

	if (mother1 != mother2)
		return;

	if (mother1->GetPdgCode() != kPi0)
		return;

	// Check if the selected \pi^{0} is primary
	//
	Bool_t primary = IsPrimary(mother1);

	if (primary)
		fPrimaryPi0[kReconstructed]->FillS(ma12, pt12);
	else
		fSecondaryPi0[kReconstructed]->FillS(ma12, pt12);

	// Looking at the source of pi0
	Int_t source_label = mother1->GetMother();

	// It's not decay pi0
	if (source_label == -1)
		return;

	AliAODMCParticle * hadron = dynamic_cast<AliAODMCParticle *> (eflags.fMcParticles->At(source_label));

	if (!hadron)
		return;

	Int_t hcode = hadron->GetPdgCode();
	fFeedDownPi0[kReconstructed]->FillAll(hcode, ma12, pt12);
}


//________________________________________________________________
void AliPP13MesonSelectionMC::InitSelectionHistograms()
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


	Float_t ptbins[] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0};
	Int_t ptsize = sizeof(ptbins) / sizeof(Float_t);

	// Sources of neutral pions, as a histogram
	for (Int_t i = 0; i < 2; ++i)
	{
		Int_t sstart = -10000;
		Int_t sstop = 10000 + 1;
		Int_t sbins = sstop - sstart;

		const char * s = (i == 0) ? "secondary" : "primary";
		fPi0Sources[i] = new TH1F(Form("hMC_%s_sources_%s", fPartNames[kPi0].Data(), s), Form("Sources of %s %ss ; PDG code", s, fPartNames[kPi0].Data()), sbins, sstart, sstop);
		fListOfHistos->Add(fPi0Sources[i]);
	}

	// Fill Generated histograms
	const char * np = fPartNames[kPi0];
	TH1 * hist1 = new TH1F(Form("hPt_%s_primary_", np), "Distribution of primary #pi^{0}s from primary ; p_{T} (GeV/#it{c})", ptsize - 1, ptbins);
	TH1 * hist2 = new TH1F(Form("hPt_%s_secondary_", np), "Distribution of secondary #pi^{0}s from secondary ; p_{T} (GeV/#it{c})", ptsize - 1, ptbins);
	TH1 * hist3 = new TH1F(Form("hPt_%s_feeddown_", np), "Distribution of primary #pi^{0}s from secondary ; p_{T} (GeV/#it{c})", ptsize - 1, ptbins);
	TH1 * hist4 = new TH2F(Form("hMassPt_%s_primary_", np), "(M,p_{T})_{#gamma#gamma} from primary ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax);
	TH1 * hist5 = new TH2F(Form("hMassPt_%s_secondary_", np), "(M,p_{T})_{#gamma#gamma} from secondary ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax);
	TH1 * hist6 = new TH2F(Form("hMassPt_%s_feeddown_", np), "(M,p_{T})_{#gamma#gamma} from secondary ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax);

	fPrimaryPi0[kGenerated]       = new AliPP13ParticlesHistogram(hist1, fListOfHistos, fPi0SourcesNames);
	fSecondaryPi0[kGenerated]     = new AliPP13ParticlesHistogram(hist2, fListOfHistos, fPi0SourcesNames);
	fFeedDownPi0[kGenerated]      = new AliPP13ParticlesHistogram(hist3, fListOfHistos, fPi0SourcesNames);
	fPrimaryPi0[kReconstructed]   = new AliPP13ParticlesHistogram(hist4, fListOfHistos, fPi0SourcesNames);
	fSecondaryPi0[kReconstructed] = new AliPP13ParticlesHistogram(hist5, fListOfHistos, fPi0SourcesNames);
	fFeedDownPi0[kReconstructed]  = new AliPP13ParticlesHistogram(hist6, fListOfHistos, fPi0SourcesNames);

	
	for (EnumNames::iterator i = fPartNames.begin(); i != fPartNames.end(); ++i)
	{
		const char * n = (const char *) i->second.Data();
		fSpectrums[i->first] = new ParticleSpectrum(n, fListOfHistos, ptsize - 1, ptbins, i->first != kPi0);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}


void AliPP13MesonSelectionMC::ConsiderGeneratedParticles(const EventFlags & flags)
{
	if (!flags.fMcParticles)
		return;

	for (Int_t i = 0; i < flags.fMcParticles->GetEntriesFast(); i++)
	{
		AliAODMCParticle * particle = ( AliAODMCParticle *) flags.fMcParticles->At(i);
		Int_t code = TMath::Abs(particle->GetPdgCode());

		// NB: replace this condition by find, if the number of particles will grow
		//
		if (code != kGamma && code != kPi0 && code != kEta)
			continue;


		Double_t pt = particle->Pt();


		// Use this to remove forward photons that can modify our true efficiency
		if (TMath::Abs(particle->Y()) > 0.5) // NB: Use rapidity instead of pseudo rapidity!
			continue;

		Double_t r = TMath::Sqrt(particle->Xv() * particle->Xv() + particle->Yv() * particle->Yv());

		fSpectrums[code]->fPt->Fill(pt);
		fSpectrums[code]->fPtRadius->Fill(pt, r);

		Bool_t primary = IsPrimary(particle);


		if (primary && particle->E() > 0.3)
		{
			fSpectrums[code]->fPtLong->Fill(pt);
			fSpectrums[code]->fPtAllRange->Fill(pt);
			fSpectrums[code]->fEtaPhi->Fill(particle->Phi(), particle->Y());
		}

		if (code != kPi0)
		{
			fSpectrums[code]->fPtPrimaries[Int_t(primary)]->Fill(pt);
			continue;
		}

		// TODO: Scale input distribution
		ConsiderGeneratedPi0(i, pt, primary, flags);
	}
}

void AliPP13MesonSelectionMC::ConsiderGeneratedPi0(Int_t i, Double_t pt, Bool_t primary, const EventFlags & flags)
{
	// Reject MIPS and count again
	if (pt < 0.3)
		return;

	if (primary)
		fPrimaryPi0[kGenerated]->FillS(pt);
	else
		fSecondaryPi0[kGenerated]->FillS(pt);

	AliAODMCParticle * parent = GetParent(i, flags.fMcParticles);

	if (!parent)
		return;

	Int_t pcode = parent->GetPdgCode();

	fPi0Sources[Int_t(primary)]->Fill(pcode);

	if (primary)
	{
		fPrimaryPi0[kGenerated]->Fill(pcode, pt);
		return;
	}

	// Only for decay pi0s
	//
	if (!IsPrimary(parent))
	{
		fSecondaryPi0[kGenerated]->Fill(pcode, pt);
		return;
	}

	fFeedDownPi0[kGenerated]->FillAll(pcode, pt);
}

//________________________________________________________________
AliAODMCParticle * AliPP13MesonSelectionMC::GetParent(Int_t label, Int_t & plabel, TClonesArray * particles) const
{
	if (label <= -1)
		return 0;

	// Int_t primLabel = cluster->GetLabelAt(0) ;
	// Particle # reached PHOS front surface
	AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle * >(particles->At(label));

	if (!particle)
		return 0;

	plabel = particle->GetMother();

	if (plabel <= -1)
		return 0;

	AliAODMCParticle * parent = dynamic_cast<AliAODMCParticle * >(particles->At(plabel));
	return parent;
}
