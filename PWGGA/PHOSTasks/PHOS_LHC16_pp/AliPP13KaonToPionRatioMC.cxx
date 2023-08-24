
// #include "iterator"

// --- Custom header files ---
#include "AliPP13KaonToPionRatioMC.h"

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


ClassImp(AliPP13KaonToPionRatioMC);


//________________________________________________________________
void AliPP13KaonToPionRatioMC::InitSelectionHistograms()
{
	for (EnumNames::iterator i = fPartNames.begin(); i != fPartNames.end(); ++i)
	{
		fAll[i->first] = new TH1F(Form("hPt_%s_", i->second.Data()), Form("Generated Spectrum of %s; p_{T} (GeV/#it{c})", i->second.Data()), 2000, 0, 20);
		fListOfHistos->Add(fAll[i->first]);

		fPrimary[i->first] = new TH1F(Form("hPt_%s_primary", i->second.Data()), Form("Generated Spectrum of primary %s; p_{T} (GeV/#it{c})", i->second.Data()), 2000, 0, 20);
		fListOfHistos->Add(fPrimary[i->first]);
	}

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}


void AliPP13KaonToPionRatioMC::ConsiderGeneratedParticles(const EventFlags & flags)
{
	if (!flags.fMcParticles)
		return;

	for (Int_t i = 0; i < flags.fMcParticles->GetEntriesFast(); i++)
	{
		AliAODMCParticle * particle = ( AliAODMCParticle *) flags.fMcParticles->At(i);
		Int_t code = particle->GetPdgCode();

		// NB: replace this condition by find, if the number of particles will grow
		//
		if (fAll.find(code) == fAll.end())
			continue;

		Double_t pt = particle->Pt();
		fAll[code]->Fill(pt);

		if (IsPrimary(particle))
			fPrimary[code]->Fill(pt);
	}
}
