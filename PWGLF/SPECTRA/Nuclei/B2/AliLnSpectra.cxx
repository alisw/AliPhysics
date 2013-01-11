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

// invariant differential yields and cross sections
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

AliLnSpectra::AliLnSpectra(const TString& particle, const TString& ptFilename, const TString& tag, const TString& outputFilename, const TString& otag, const Double_t xsec[3])
: TObject()
, fParticle(particle)
, fPtFilename(ptFilename)
, fTag(tag)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
, fYMin(-0.5)
, fYMax(0.5)
, fNormToInel(1)
, fIsOnlyGen(0)
, fSysErr(1)
{
//
// constructor
//
	for(Int_t i=0; i<3; ++i) fInelXsec[i] = xsec[i];
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
// (invariant) differential yield projection in pt for |y| < 0.5
//
	using namespace std;
	
	TFile* finput = new TFile(fPtFilename.Data(), "read");
	if(finput->IsZombie()) exit(1);
	
	// pt and number of events
	
	TH1D* hPt = (TH1D*)FindObj(finput, fParticle + "_Pt");
	TH1D* hStats = (TH1D*)FindObj(finput, fParticle + "_Stats");
	
	Double_t nEvent = hStats->Integral(3,3);
	if(fNormToInel) nEvent = hStats->Integral(4,4);
	
	// ouputfile
	
	TFile* foutput = new TFile(fOutputFilename.Data(),"recreate");
	
	foutput->mkdir(fOutputTag.Data());
	foutput->cd(fOutputTag.Data());
	
	// differential yield
	
	TGraphErrors* grDYieldPt = this->GetDiffYieldPt(hPt, nEvent, fParticle + "_DiffYield_Pt");
	grDYieldPt->Write();
	
	TGraphErrors* grSysErrDYieldPt = this->AddSystError(grDYieldPt, fSysErr, fParticle + "_SysErr_DiffYield_Pt");
	grSysErrDYieldPt->Write();
	
	TF1* fncTsallis0 = this->TsallisDiffYield(GetMass(fParticle), fParticle + "_Fit_DiffYield_Pt");
	fncTsallis0->SetParLimits(0, 0, 1);
	fncTsallis0->SetParLimits(1, 4, 50);
	fncTsallis0->SetParLimits(2, 0.01, 10);
	
	grSysErrDYieldPt->Fit(fncTsallis0,"RNQ");
	fncTsallis0->Write();
	
	// invariant differential yield
	
	TGraphErrors* grInvDYieldPt = this->GetInvDiffYieldPt(grDYieldPt, fParticle + "_InvDiffYield_Pt");
	grInvDYieldPt->Write();
	
	TGraphErrors* grSysErrInvDYieldPt = this->AddSystError(grInvDYieldPt, fSysErr, fParticle + "_SysErr_InvDiffYield_Pt");
	grSysErrInvDYieldPt->Write();
	
	TF1* fncTsallis1 = this->Tsallis(GetMass(fParticle), fParticle + "_Fit_InvDiffYield_Pt");
	fncTsallis1->SetParLimits(0, 0, 1);
	fncTsallis1->SetParLimits(1, 4, 50);
	fncTsallis1->SetParLimits(2, 0.01, 10);
	
	grSysErrInvDYieldPt->Fit(fncTsallis1,"RNQ");
	fncTsallis1->Write();
	
	// invariant differential cross section
	
	TGraphErrors* grInvXsectPt = GetInvDiffXsectionPt(grInvDYieldPt, fInelXsec, fParticle + "_InvDiffXSection_Pt");
	grInvXsectPt->Write();
	
	TGraphErrors* grSysErrInvXsectPt = GetInvDiffXsectionPt(grSysErrInvDYieldPt, fInelXsec, fParticle + "_SysErr_InvDiffXSection_Pt");
	grSysErrInvXsectPt->Write();
	
	TF1* fncTsallis2 = this->Tsallis(GetMass(fParticle), fInelXsec[0], fParticle + "_Fit_InvDiffXSection_Pt");
	fncTsallis2->SetParLimits(0, 0, 1);
	fncTsallis2->SetParLimits(1, 4, 50);
	fncTsallis2->SetParLimits(2, 0.01, 10);
	
	grSysErrInvXsectPt->Fit(fncTsallis2,"RNQ");
	fncTsallis2->Write();
	
	// clean
	
	delete fncTsallis0;
	delete fncTsallis1;
	delete fncTsallis2;
	
	delete grDYieldPt;
	delete grInvDYieldPt;
	delete grSysErrDYieldPt;
	delete grSysErrInvDYieldPt;
	delete grInvXsectPt;
	
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
		
		grDYieldPt->SetPoint(j,pt,yield);
		grDYieldPt->SetPointError(j++,0,yError);
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

TGraphErrors* AliLnSpectra::GetInvDiffXsectionPt(const TGraphErrors* grInvDYieldPt, const Double_t* sigma, const TString& name) const
{
//
// invariant differential cross section, only stat. errors
// (multiply each value for the total cross section)
//
	Double_t xsec = sigma[0];
	Double_t errXSec = sigma[1];
	
	TGraphErrors* grInvDXsecPt = new TGraphErrors();
	grInvDXsecPt->SetName(name.Data());
	
	for(Int_t i=0; i < grInvDYieldPt->GetN(); ++i)
	{
		Double_t pt, y;
		grInvDYieldPt->GetPoint(i,pt,y);
		
		Double_t errpt = grInvDYieldPt->GetErrorX(i);
		Double_t erry = grInvDYieldPt->GetErrorY(i);
		
		Double_t invxsec = xsec*y;
		Double_t err = invxsec*TMath::Sqrt(TMath::Power(errXSec/xsec,2) + TMath::Power(erry/y,2));
		
		grInvDXsecPt->SetPoint(i, pt, invxsec);
		grInvDXsecPt->SetPointError(i, errpt, err);
	}
	
	return grInvDXsecPt;
}

TGraphErrors* AliLnSpectra::AddSystError(const TGraphErrors* gr, Double_t percent, const TString& name) const
{
//
// set the error of h as the given percent of its value
//
	TGraphErrors* grSyst = new TGraphErrors();
	grSyst->SetName(name.Data());
	
	for(Int_t i=0; i < gr->GetN(); ++i)
	{
		Double_t x, y;
		gr->GetPoint(i,x,y);
		Double_t err = percent*y;
		
		grSyst->SetPoint(i, x, y);
		grSyst->SetPointError(i, 0.025, err);
	}
	
	return grSyst;
}

TF1* AliLnSpectra::Tsallis(Double_t m0, const TString& name, Double_t xmin, Double_t xmax) const
{
//
// Tsallis distribution
// Phys. Rev. C 83, 064903 (2011)
// Phys. Rev. C 75, 064901 (2007)
//
	TF1* fnc = new TF1(name.Data(), Form("[0]*([1]-1)*([1]-2)*TMath::Power(1+(sqrt(x*x+%f*%f)-%f)/([1]*[2]),-[1])/(2*TMath::Pi()*[1]*[2]*([1]*[2]+%f*([1]-2)))", m0, m0, m0, m0), xmin, xmax);
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.1, 7, 0.2);
	
	return fnc;
}

TF1* AliLnSpectra::Tsallis(Double_t m0, Double_t xsect, const TString& name, Double_t xmin, Double_t xmax) const
{
//
// Tsallis distribution to fit to invariant cross section
// Phys. Rev. C 83, 064903 (2011)
// Phys. Rev. C 75, 064901 (2007)
//
	TF1* fnc = new TF1(name.Data(), Form("%f*[0]*([1]-1)*([1]-2)*TMath::Power(1+(sqrt(x*x+%f*%f)-%f)/([1]*[2]),-[1])/(2*TMath::Pi()*[1]*[2]*([1]*[2]+%f*([1]-2)))",xsect, m0, m0, m0, m0), xmin, xmax);
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.1, 7, 0.2);
	
	return fnc;
}

TF1* AliLnSpectra::TsallisDiffYield(Double_t m0, const TString& name, Double_t xmin, Double_t xmax) const
{
//
// Tsallis distribution to fit differential yield
// Phys. Rev. C 83, 064903 (2011)
// Phys. Rev. C 75, 064901 (2007)
//
	TF1* fnc = new TF1(name.Data(), Form("x*[0]*([1]-1)*([1]-2)*TMath::Power(1+(sqrt(x*x+%f*%f)-%f)/([1]*[2]),-[1])/([1]*[2]*([1]*[2]+%f*([1]-2)))", m0, m0, m0, m0), xmin, xmax);
	fnc->SetParNames("dN/dy","n","C");
	fnc->SetParameters(0.1, 7, 0.2);
	
	return fnc;
}
