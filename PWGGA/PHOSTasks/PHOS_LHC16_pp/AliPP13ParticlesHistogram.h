#ifndef ALIPP13PARTICLESHISTOGRAM_H 
#define ALIPP13PARTICLESHISTOGRAM_H 

#include <map>

// --- ROOT system ---
#include <TH1.h>
#include <TList.h>


class AliPP13ParticlesHistogram
{
public:
	typedef std::map<Int_t, TString> EnumNames;
	typedef std::map<Int_t, TH1 *> HistogramList;
	AliPP13ParticlesHistogram(TH1 * hist, TList * owner, EnumNames & sources);

	void FillS(Float_t x, Float_t y = 1.0);
	void Fill(Int_t pdg, Float_t x, Float_t y = 1.0);
	void FillAll(Int_t pdg, Float_t x, Float_t y = 1.0);

private:
	AliPP13ParticlesHistogram(const AliPP13ParticlesHistogram &);

	EnumNames fSources;
	HistogramList fHistograms;
	const Int_t kAll = -9999;
};
#endif
