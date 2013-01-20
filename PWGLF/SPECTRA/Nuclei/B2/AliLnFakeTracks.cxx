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

// fake track fraction
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TF1.h>

#include "AliLnFakeTracks.h"
#include "B2.h"

ClassImp(AliLnFakeTracks)

AliLnFakeTracks::AliLnFakeTracks(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag)
: TObject()
, fParticle(particle)
, fSimuFilename(simuFilename)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
{
//
// constructor
//
	TH1::SetDefaultSumw2();
}

AliLnFakeTracks::~AliLnFakeTracks()
{
//
// destructor
//
}

Int_t AliLnFakeTracks::Exec()
{
//
// extract the fraction of fake tracks for the given fParticle.Data()
//
	using namespace std;
	
	TFile* fsimu = new TFile(fSimuFilename.Data(), "read");
	if (fsimu->IsZombie()) exit(1);
	
	TFile* foutput = new TFile(fOutputFilename.Data(), "recreate");
	if(fOutputTag != "")
	{
		foutput->mkdir(fOutputTag.Data());
		foutput->cd(fOutputTag.Data());
	}
	
	TH1D* hFakeTracks = (TH1D*)FindObj(fsimu, fParticle + "_Sim_Fake_Pt");
	TH1D* hAllTracks  = (TH1D*)FindObj(fsimu, fParticle + "_Sim_Pt");
	
	TH1D* hFracFakePt = Divide(hFakeTracks, hAllTracks, fParticle + "_Frac_Fake_Pt");
	
	hFracFakePt->SetYTitle("Fake Tracks/All Tracks");
	hFracFakePt->Write();
	
	delete hFracFakePt;
	
	delete foutput;
	delete fsimu;
	
	return 0;
}
