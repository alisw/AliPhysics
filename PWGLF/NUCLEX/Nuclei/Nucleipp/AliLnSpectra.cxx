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

// (invariant) differential yield
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "AliLnSpectra.h"
#include "B2.h"

ClassImp(AliLnSpectra)

AliLnSpectra::AliLnSpectra(const TString& particle, const TString& ptFilename, const TString& tag, const TString& outputFilename, const TString& otag)
: TObject()
, fParticle(particle)
, fPtFilename(ptFilename)
, fTag(tag)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
, fYMin(-0.5)
, fYMax(0.5)
, fINEL(1)
{
//
// constructor
//
}

AliLnSpectra::~AliLnSpectra()
{
//
// destructor
//
}

Int_t AliLnSpectra::Exec()
{
//
// (invariant) differential yield projection in pt for |y| < ymax-ymin
//
	using namespace std;
	
	TFile* finput = new TFile(fPtFilename.Data(), "read");
	if(finput->IsZombie()) exit(1);
	
	// pt and number of events
	
	TH1D* hPt = FindObj<TH1D>(finput, fParticle + "_Pt");
	TH1D* hStats = FindObj<TH1D>(finput, fParticle + "_Stats");
	
	Double_t nEvent = hStats->Integral(3,3);
	if(fINEL) nEvent = hStats->Integral(4,4);
	
	// ouputfile
	
	TFile* foutput = new TFile(fOutputFilename.Data(),"recreate");
	
	foutput->mkdir(fOutputTag.Data());
	foutput->cd(fOutputTag.Data());
	
	// differential yield
	
	TGraphErrors* grDYieldPt = this->GetDiffYieldPt(hPt, nEvent, fParticle + "_DiffYield_Pt");
	grDYieldPt->Write();
	
	// invariant differential yield
	
	TGraphErrors* grInvDYieldPt = this->GetInvDiffYieldPt(grDYieldPt, fParticle + "_InvDiffYield_Pt");
	grInvDYieldPt->Write();
	
	// clean
	
	delete grDYieldPt;
	delete grInvDYieldPt;
	
	delete foutput;
	delete finput;
	
	return 0;
}

TGraphErrors* AliLnSpectra::GetDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const
{
//
// projection of the differential yield in pt
// only statistical error of N is taken into account
//
	TGraphErrors* grDYieldPt = new TGraphErrors();
	grDYieldPt->SetName(name.Data());
	
	Double_t dy = fYMax - fYMin;
	
	TAxis* xaxis = hPt->GetXaxis();
	
	for(Int_t i=1, j=0; i < xaxis->GetNbins()+1; ++i)
	{
		Double_t pt  = xaxis->GetBinCenter(i);
		Double_t dpt = xaxis->GetBinWidth(i);
		
		Double_t n = hPt->GetBinContent(i);
		Double_t eN = hPt->GetBinError(i);
		
		if(n==0) continue;
		
		Double_t yield = n/(nEvent*dpt*dy);
		Double_t yError = yield*eN/n;
		
		grDYieldPt->SetPoint(j, pt, yield);
		grDYieldPt->SetPointError(j++, dpt/2., yError);
	}
	
	return grDYieldPt;
}

TGraphErrors* AliLnSpectra::GetInvDiffYieldPt(const TGraphErrors* grDYieldPt, const TString& name) const
{
//
// projection of the invariant differential yield in pt
// only statistical error of N is taken into account
//
	TGraphErrors* grInvDYieldPt = new TGraphErrors();
	grInvDYieldPt->SetName(name.Data());
	
	for(Int_t i=0; i < grDYieldPt->GetN(); ++i)
	{
		Double_t pt, y;
		grDYieldPt->GetPoint(i,pt,y);
		
		Double_t errpt = grDYieldPt->GetErrorX(i);
		Double_t erry = grDYieldPt->GetErrorY(i);
		
		y = y/(2.*TMath::Pi()*pt);
		erry = erry / (2.*TMath::Pi()*pt);
		
		grInvDYieldPt->SetPoint(i, pt, y);
		grInvDYieldPt->SetPointError(i, errpt, erry);
	}
	
	return grInvDYieldPt;
}

TGraphErrors* AliLnSpectra::GetInvDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const
{
//
// projection of the invariant differential yield in pt
// only statistical error of N is taken into account
//
	TGraphErrors* grDYieldPt = this->GetDiffYieldPt(hPt, nEvent, "_tmp_diff_yield_");
	TGraphErrors* grInvDYieldPt = this->GetInvDiffYieldPt(grDYieldPt, name);
	
	delete grDYieldPt;
	return grInvDYieldPt;
}
