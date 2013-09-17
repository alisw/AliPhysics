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

// antiparticle / particle ratio
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TString.h>
#include <TROOT.h>

#include "AliLnRatio.h"
#include "B2.h"

ClassImp(AliLnRatio)

AliLnRatio::AliLnRatio(const TString& species, const TString& ptFilename, const TString& tag, const TString& outputFilename, const TString& otag)
: TObject()
, fSpecies(species)
, fPtFilename(ptFilename)
, fPtTag(tag)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
{
//
// constructor
//
}

AliLnRatio::~AliLnRatio()
{
//
// destructor
//
	TH1::SetDefaultSumw2();
}

Int_t AliLnRatio::Exec()
{
//
// antiparticle/particle ratio
//
	using namespace std;
	
	TFile* finput = new TFile(fPtFilename.Data(), "read");
	if (finput->IsZombie()) exit(1);
	
	const TString kPrefix[] = { "", "Anti" };
	
	TH1D* hPt[2];
	
	for(Int_t i=0; i<2; ++i)
	{
		hPt[i] = FindObj<TH1D>(finput, kPrefix[i] + fSpecies + "_Pt");
	}
	
	TH1D* hRatioPt = Divide(hPt[1],hPt[0],Form("Anti%s%s_Ratio_Pt",fSpecies.Data(),fSpecies.Data()));
	
	TFile* foutput = new TFile(fOutputFilename.Data(), "recreate");
	
	foutput->mkdir(fOutputTag.Data());
	foutput->cd(fOutputTag.Data());
	
	hRatioPt->Write();
	
	delete hRatioPt;
	delete foutput;
	delete finput;
	
	return 0;
}
