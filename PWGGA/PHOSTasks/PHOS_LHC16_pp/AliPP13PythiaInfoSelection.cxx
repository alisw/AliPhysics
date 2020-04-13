// --- Custom header files ---
#include "AliPP13PythiaInfoSelection.h"

// --- ROOT system ---
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliAnalysisManager.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13PythiaInfoSelection);

//________________________________________________________________
void AliPP13PythiaInfoSelection::CountMBEvent()
{

	AliPP13PhysicsSelection::CountMBEvent();

	// Fetch the histgram file
	TTree * tree = AliAnalysisManager::GetAnalysisManager()->GetTree();

	if (!tree)
	{
		AliError(Form("%s - UserNotify: No current tree!", GetName()));
		return;
	}

	TFile * curfile = tree->GetCurrentFile();
	if (!curfile)
	{
		AliError(Form("%s - UserNotify: No current file!", GetName()));
		return;
	}

	TString file(curfile->GetName());
	file.ReplaceAll(gSystem->BaseName(file.Data()), "");

	TFile * xsecFile = TFile::Open(Form("%s%s", file.Data(), "pyxsec_hists.root"));
	if (!xsecFile)
	{
		AliError(Form("There is no pyxsec_hists.root in this directory."));
		return;
	}

	// find the tlist we want to be independtent of the name so use the Tkey
	TKey * key = (TKey *)xsecFile->GetListOfKeys()->At(0);
	if (!key)
		return;

	TList * list = dynamic_cast<TList *>(key->ReadObj());
	if (!list)
		return;

	Float_t xsec    = ((TProfile *)list->FindObject("h1Xsec"))  ->GetBinContent(1);
	Float_t trials  = ((TH1F *)    list->FindObject("h1Trials"))->GetBinContent(1);

	xsecFile->Close();

	fXsec->Fill(0.5, xsec);
	fTrials->Fill(0.5, trials);
}

//________________________________________________________________
void AliPP13PythiaInfoSelection::FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags)
{
	(void) clusArray;
	(void) pool;
	(void) eflags;

}

//________________________________________________________________
void AliPP13PythiaInfoSelection::InitSelectionHistograms()
{

	fXsec = new TH1F("hXsec", "xsec from pyxsec.root", 1, 0, 1);
	fXsec->GetXaxis()->SetBinLabel(1, "<#sigma>");
	fListOfHistos->Add(fXsec);

	fTrials = new TH1F("hTrials", "trials root file", 1, 0, 1);
	fTrials->GetXaxis()->SetBinLabel(1, "#sum_{ntrials}");
	fListOfHistos->Add(fTrials);

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}
