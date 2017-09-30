// --- Custom header files ---
#include "AliPP13DetectorHistogram.h"

// --- ROOT system ---
#include <TString.h>
#include <TH2.h>


#include <iostream>
using namespace std;


//
//________________________________________________________________
AliPP13DetectorHistogram::AliPP13DetectorHistogram(TH1 * hist, TList * owner, AliPP13DetectorHistogram::Mode mode):
	fHistogram123(0),
	fHistograms(),
	fInterModuleHistograms(),
	fMode(mode)

{
	if (!owner->IsOwner())
		cout << "Warning: You are adding histograms to the list that doesn't have ownership" << endl;

	// TODO: Use factory method here once it will be clear what are the requiremets for this class
	//

	TString name = hist->GetName();
    TString title = hist->GetTitle();

	fHistogram123 = dynamic_cast<TH1 *>( hist->Clone(name + "SM123") );
	fHistogram123->SetTitle(Title(title.Data(), -1));
	owner->Add(fHistogram123);

	for (Int_t sm = 0; sm < (kLastModule + 1); ++sm)
	{
		if (sm > 1 && (fMode == kSingleHist || fMode == kInterModules))
			continue;

		fHistograms[sm] = (sm == 0) ? hist : dynamic_cast<TH1 *>( hist->Clone(name + Form("SM%d", sm)) );
		if (sm == 0 && fMode == kModules)
			fHistograms[sm]->SetName(name + Form("SM%d", sm));

		fHistograms[sm]->SetTitle(Title(title.Data(), sm));
		owner->Add(fHistograms[sm]);
	}

	if (fMode != kInterModules)
		return;

	for (Int_t sm = kFirstModule; sm < (kLastModule + 1); ++sm)
	{
		for (Int_t sm2 = sm; sm2 < (kLastModule + 1); ++sm2)
		{
			fInterModuleHistograms[Index(sm, sm2)] = dynamic_cast<TH1 *>(hist->Clone(name + Form("SM%dSM%d", sm, sm2)));
			fInterModuleHistograms[Index(sm, sm2)]->SetTitle(Title(title, sm, sm2));
			owner->Add(fInterModuleHistograms[Index(sm, sm2)]);
		}
	}
}


//________________________________________________________________
void AliPP13DetectorHistogram::FillAll(Int_t sm1, Int_t sm2, Float_t x, Float_t y)
{
	if (sm1 < kFirstModule || sm1 > kLastModule)
	{
		cout << "Illegal module: " << sm1 << endl;
		return;
	}

	if (sm2 < kFirstModule || sm2 > kLastModule)
	{
		cout << "Illegal module: " << sm2 << endl;
		return;
	}

	// Fill histograms for all the modules
	//
	fHistograms[0]->Fill(x, y);

	// Fill histograms without specific Module
	//
	if (sm1 != 4 &&  sm2 != 4)
		fHistogram123->Fill(x, y);

	FillModules(sm1, sm2, x, y);
}

//________________________________________________________________
void AliPP13DetectorHistogram::FillAll(Int_t sm1, Int_t sm2, Float_t x, Float_t y, Float_t z)
{
	if (sm1 < kFirstModule || sm1 > kLastModule)
	{
		cout << "Illegal module: " << sm1 << endl;
		return;
	}

	if (sm2 < kFirstModule || sm2 > kLastModule)
	{
		cout << "Illegal module: " << sm2 << endl;
		return;
	}

	// Fill histograms for all the modules
	//
	dynamic_cast<TH2 *>(fHistograms[0])->Fill(x, y, z);

	// Fill histograms without specific Module
	//
	if (sm1 != 4 &&  sm2 != 4)
		dynamic_cast<TH2 *>(fHistogram123)->Fill(x, y, z);

	FillModules(sm1, sm2, x, y, z);
}

//________________________________________________________________
void AliPP13DetectorHistogram::FillModules(Int_t sm1, Int_t sm2, Float_t x, Float_t y, Float_t z)
{

	if (sm1 == sm2 && fMode == kModules)
		dynamic_cast<TH2 *>(fHistograms[sm1])->Fill(x, y, z);

	if (fMode == kInterModules)
		dynamic_cast<TH2 *>(fInterModuleHistograms[Index(sm1, sm2)])->Fill(x, y, z);
}



//________________________________________________________________
void AliPP13DetectorHistogram::FillModules(Int_t sm1, Int_t sm2, Float_t x, Float_t y)
{

	if (sm1 == sm2 && fMode == kModules)
		fHistograms[sm1]->Fill(x, y);

	if (fMode == kInterModules)
		fInterModuleHistograms[Index(sm1, sm2)]->Fill(x, y);
}

//________________________________________________________________
TString AliPP13DetectorHistogram::Title(TString title, Int_t i) const
{
	if (i == -1)
		return title + "SM123";

	TString s = (i == 0) ? "all modules" : Form("SM%d", i);
	return title + s;
}

//________________________________________________________________
TString AliPP13DetectorHistogram::Title(TString title, Int_t i, Int_t j) const
{
	TString s = (i == j) ? Form("SM%d", i) : Form("SM%dSM%d", i, j);
	return title + s;
}

//________________________________________________________________
Int_t AliPP13DetectorHistogram::Index(Int_t sm1, Int_t sm2) const
{
	if (sm1 > sm2)
		swap(sm1, sm2);

	return sm2 * (sm2 + 1) / 2 + sm1;
}