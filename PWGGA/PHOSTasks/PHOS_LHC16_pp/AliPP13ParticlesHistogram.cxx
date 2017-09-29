// --- Custom header files ---
#include "AliPP13ParticlesHistogram.h"

// --- ROOT system ---
#include <TString.h>


#include <iostream>
using namespace std;


//________________________________________________________________
AliPP13ParticlesHistogram::AliPP13ParticlesHistogram(TH1 * hist, TList * owner, EnumNames & sources):
	fSources(sources)
{
	TString name = hist->GetName();
	TString title = hist->GetTitle();

	// Add all histogram
	fSources[kAll] = "";
	for (EnumNames::iterator s = fSources.begin(); s != fSources.end(); ++s)
	{
		const char * ns = (const char *) s->second.Data();
		fHistograms[s->first] = ns ? dynamic_cast<TH1 *>(hist->Clone(name + ns)) : hist;
		fHistograms[s->first]->SetTitle(title + ns);
		owner->Add(fHistograms[s->first]);
	}
}



//________________________________________________________________
void AliPP13ParticlesHistogram::FillS(Float_t x, Float_t y)
{
	fHistograms[kAll]->Fill(x, y);
}


//________________________________________________________________
void AliPP13ParticlesHistogram::Fill(Int_t pdg, Float_t x, Float_t y)
{
	EnumNames::iterator s = fSources.find(pdg);
	if (s == fSources.end())
		return;

	fHistograms[pdg]->Fill(x, y);
}


//________________________________________________________________
void AliPP13ParticlesHistogram::FillAll(Int_t pdg, Float_t x, Float_t y)
{
	FillS(x, y);
	Fill(pdg, x, y);
}
