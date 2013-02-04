/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// efficiency correction
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TF1.h>

#include "AliLnEfficiency.h"
#include "B2.h"

ClassImp(AliLnEfficiency)

AliLnEfficiency::AliLnEfficiency(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag)
: TObject()
, fParticle(particle)
, fSimuFilename(simuFilename)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
, fG3Fluka(0)
, fAddFakeTracks(1)
{
//
// constructor
//
	TH1::SetDefaultSumw2();
}

AliLnEfficiency::~AliLnEfficiency()
{
//
// destructor
//
}

Int_t AliLnEfficiency::Exec()
{
//
// extract the efficiency for the given simulation
//
	using namespace std;
	
	TFile* fsimu = new TFile(fSimuFilename.Data(),"read");
	if (fsimu->IsZombie()) exit(1);
	
	TFile* foutput = new TFile(fOutputFilename.Data(),"recreate");
	if(fOutputTag != "")
	{
		foutput->mkdir(fOutputTag.Data());
		foutput->cd(fOutputTag.Data());
	}
	
	TH1D* hGenPhSpPt = (TH1D*)FindObj(fsimu, fParticle + "_Gen_PhS_Prim_Pt");
	TH1D* hGenTrigPt = (TH1D*)FindObj(fsimu, fParticle + "_Gen_Trig_Prim_Pt");
	TH1D* hGenVtxPt  = (TH1D*)FindObj(fsimu, fParticle + "_Gen_Vtx_Prim_Pt");
	TH1D* hGenAccPt  = (TH1D*)FindObj(fsimu, fParticle + "_Gen_Acc_Prim_Pt");
	TH1D* hPrimRecPt = (TH1D*)FindObj(fsimu, fParticle + "_Sim_Prim_Pt");
	TH1D* hFakePrimRecPt = (TH1D*)FindObj(fsimu, fParticle + "_Sim_Fake_Prim_Pt");
	
	if(fAddFakeTracks) hPrimRecPt->Add(hFakePrimRecPt);
	
	TH1D* hEffTrigPt   = Divide(hGenTrigPt, hGenPhSpPt, fParticle + "_Eff_Trig_Pt");
	TH1D* hEffVtxPt    = Divide(hGenVtxPt, hGenTrigPt, fParticle + "_Eff_Vtx_Pt");
	TH1D* hEffAccPt    = Divide(hGenAccPt, hGenVtxPt, fParticle + "_Eff_Acc_Pt");
	TH1D* hEffTrkPt    = Divide(hPrimRecPt, hGenAccPt, fParticle + "_Eff_Trk_Pt");
	TH1D* hEffAccTrkPt = Divide(hPrimRecPt, hGenVtxPt, fParticle + "_Eff_AccTrk_Pt");
	
	if(fG3Fluka)
	{
		TF1* fncG3Fluka = 0;
		if(fParticle == "Proton") fncG3Fluka = GetProtonG3FlukaCor("Proton_G3FLUKA_Corr_Pt");
		if(fParticle == "AntiProton") fncG3Fluka = GetAntiProtonG3FlukaCor("AntiProton_G3FLUKA_Corr_Pt");
		
		if(fncG3Fluka != 0)
		{
			hEffTrkPt->Divide(fncG3Fluka);
			hEffAccTrkPt->Divide(fncG3Fluka);
		}
		
		delete fncG3Fluka;
	}
	
	if(!fParticle.Contains("Anti"))
	{
		TH1D* hStats = (TH1D*)FindObj(fsimu, fParticle + "_Stats");
		hStats->Write();
	}
	
	hEffAccTrkPt->SetTitle(fParticle.Data());
	hEffTrkPt->SetTitle(fParticle.Data());
	hEffAccPt->SetTitle(fParticle.Data());
	hEffVtxPt->SetTitle(fParticle.Data());
	hEffTrigPt->SetTitle(fParticle.Data());
	
	hEffAccTrkPt->SetYTitle("Acceptance #times Track Efficiency");
	hEffAccPt->SetYTitle("Acceptance Efficiency");
	hEffTrkPt->SetYTitle("Track Efficiency");
	hEffVtxPt->SetYTitle("Vertex Efficiency");
	hEffTrigPt->SetYTitle("Trigger Efficiency");
	
	hEffAccTrkPt->Write();
	hEffAccPt->Write();
	hEffTrkPt->Write();
	hEffVtxPt->Write();
	hEffTrigPt->Write();
	
	delete hEffAccTrkPt;
	delete hEffAccPt;
	delete hEffTrkPt;
	delete hEffVtxPt;
	delete hEffTrigPt;
	
	delete foutput;
	delete fsimu;
	
	return 0;
}

TF1* AliLnEfficiency::GetProtonG3FlukaCor(const TString& name) const
{
//
// GEANT3/FLUKA correction for protons (Baryon Twiki)
// (TPC tracks)
//
	TF1* fnc = new TF1(name.Data(), "1-[0]*exp([1]*x)+[2]", 0.001, 10);
	fnc->SetParameters(4113, -32.9, -0.009383);
	Double_t err[] = {1640.4, 1.5, 0.000917};
	fnc->SetParErrors(err);
	
	return fnc;
}

TF1* AliLnEfficiency::GetAntiProtonG3FlukaCor(const TString& name) const
{
//
// GEANT3/FLUKA correction for antiprotons (Baryon Twiki)
// (TPC tracks)
//
	TF1* fnc = new TF1(name.Data(), "1-[0]*exp([1]*x)+[2]+[3]*log(x)/pow(x,0.2)", 0.001, 10);
	fnc->SetParameters(177.4, -22.03, -0.06538, 0.04431 );
	Double_t err[] = {73., 1.5, 0.00221, 0.00555};
	fnc->SetParErrors(err);
	
	return fnc;
}

