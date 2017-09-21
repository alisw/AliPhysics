// --- Custom header files ---
#include "ParticlesHistogram.h"
// #include "AliAnalysisTaskPP.h"

// --- ROOT system ---
#include <TString.h>


#include <iostream>
using namespace std;


//________________________________________________________________
ParticlesHistogram::ParticlesHistogram(TH1 * hist, TList * owner, EnumNames & sources):
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
void ParticlesHistogram::FillS(Float_t x, Float_t y)
{
	fHistograms[kAll]->Fill(x, y);
}


//________________________________________________________________
void ParticlesHistogram::Fill(Int_t pdg, Float_t x, Float_t y)
{
	EnumNames::iterator s = fSources.find(pdg);
	if (s == fSources.end())
		return;

	fHistograms[pdg]->Fill(x, y);
}


//________________________________________________________________
void ParticlesHistogram::FillAll(Int_t pdg, Float_t x, Float_t y)
{
	FillS(x, y);
	Fill(pdg, x, y);
}
