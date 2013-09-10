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

// coalescence parameter
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "AliLnB2.h"
#include "B2.h"

ClassImp(AliLnB2)

AliLnB2::AliLnB2(const TString& protonSpectra, const TString& protonTag, const TString& nucleusSpectra, const TString& nucleusTag, const TString& outputFilename, const TString& otag, Int_t a, Int_t z)
: TObject()
, fProtonSpectra(protonSpectra)
, fProtonTag(protonTag)
, fNucleusSpectra(nucleusSpectra)
, fNucleusTag(nucleusTag)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
, fA(2)
, fZ(1)
, fNucleusName("Deuteron")
, fCd(0.15)
{
//
// constructor
//
	this->SetNucleus(a,z);
}

AliLnB2::~AliLnB2()
{
//
// destructor
//
}

void AliLnB2::SetNucleus(Int_t a, Int_t z)
{
//
// set nucleus mass and name
//
	fA = a;
	fZ = z;
	
	Int_t zz = TMath::Abs(z);
	
	if(a==1 && zz==1)      fNucleusName = "Proton";
	else if(a==2 && zz==1) fNucleusName = "Deuteron";
	else if(a==3 && zz==1) fNucleusName = "Triton";
	else if(a==3 && zz==2) fNucleusName = "He3";
	else if(a==4 && zz==2) fNucleusName = "Alpha";
	else
	{
		this->Warning("SetNucleus", "unknown nucleus A = %d, Z = %d", a, z);
		fNucleusName = "Unknown";
	}
}

Int_t AliLnB2::Run()
{
//
// coalescence parameter
//
	TFile* finput1 = new TFile(fProtonSpectra.Data(), "read");
	if (finput1->IsZombie()) exit(1);
	
	TFile* finputA = new TFile(fNucleusSpectra.Data(), "read");
	if (finputA->IsZombie()) exit(1);
	
	TString prefix = "";
	TString suffix = "";
	if(fZ < 0)
	{
		prefix = "Anti";
		suffix = "bar";
	}
	
	// invariant differential yields
	
	TGraphErrors* grPrtInvDYieldPt = FindObj<TGraphErrors>(finput1, fProtonTag, prefix + "Proton_InvDiffYield_Pt");
	TGraphErrors* grNucInvDYieldPt = FindObj<TGraphErrors>(finputA, fNucleusTag, prefix + fNucleusName + "_InvDiffYield_Pt");
	
	TFile* foutput = new TFile(fOutputFilename.Data(),"recreate");
	
	foutput->mkdir(fOutputTag.Data());
	foutput->cd(fOutputTag.Data());
	
	// coalescence parameter
	
	TGraphErrors* grB2Pt = this->GetBAPt(grPrtInvDYieldPt, grNucInvDYieldPt, Form("B2%s_Pt",suffix.Data()));
	
	
	grB2Pt->Write();
	
	// homogeneity volume
	
	TGraphErrors* grR3Pt = this->Rside2Rlong(grB2Pt, Form("R3%s_Pt", suffix.Data()), fCd);
	
	grR3Pt->Write();
	
	// propagate systematic errors if possible
	
	TString grPrtName = fProtonTag +  "/" + prefix + "Proton_SystErr_InvDiffYield_Pt;1";
	TString grNucName = fNucleusTag + "/" + prefix + fNucleusName + "_SystErr_InvDiffYield_Pt;1";
	
	TGraphErrors* grSystErrPrtInvDYieldPt = dynamic_cast<TGraphErrors*>(finput1->Get(grPrtName.Data()));
	TGraphErrors* grSystErrNucInvDYieldPt = dynamic_cast<TGraphErrors*>(finput1->Get(grNucName.Data()));
	
	if( (grSystErrPrtInvDYieldPt != 0) && (grSystErrNucInvDYieldPt != 0) )
	{
		TGraphErrors* grSystErrB2Pt = this->GetBAPt(grSystErrPrtInvDYieldPt, grSystErrNucInvDYieldPt, Form("SystErr_B2%s_Pt",suffix.Data()));
		
		grSystErrB2Pt->Write();
		
		TGraphErrors* grSystErrR3Pt = this->Rside2Rlong(grSystErrB2Pt, Form("SystErr_R3%s_Pt", suffix.Data()), fCd);
		grSystErrR3Pt->Write();
		
		delete grSystErrB2Pt;
		delete grSystErrR3Pt;
	}
	else
	{
		this->Warning("Run", "systematic errors are not propagated");
	}
	
	delete grB2Pt;
	delete grR3Pt;
	
	delete foutput;
	delete finput1;
	delete finputA;
	
	return 0;
}

TGraphErrors* AliLnB2::GetBAPt(const TGraphErrors* grPrtInvDYieldPt, const TGraphErrors* grNucInvDYieldPt, const TString& name) const
{
//
// coalescence parameter
//
	TGraphErrors* grBAPt = new TGraphErrors();
	grBAPt->SetName(name.Data());
	
	for(Int_t i=0, j=0; i < grNucInvDYieldPt->GetN(); ++i)
	{
		Double_t ptNuc, yNuc;
		
		grNucInvDYieldPt->GetPoint(i, ptNuc, yNuc);
		
		if(ptNuc<0.8) continue; // acceptance
		
		Double_t yPrt = grPrtInvDYieldPt->Eval(ptNuc/fA); // interpolate
		
		if(yPrt == 0 || yNuc == 0 ) continue;
		
		Double_t bA = yNuc/TMath::Power(yPrt,fA);
		
		// error
		Double_t ePrt = this->GetErrorY(grPrtInvDYieldPt, ptNuc/fA);
		Double_t eNuc = grNucInvDYieldPt->GetErrorY(i);
		
		Double_t errPt = grNucInvDYieldPt->GetErrorX(i)/fA;
		Double_t errBA = bA*TMath::Sqrt(TMath::Power(eNuc/yNuc,2) + TMath::Power(fA*ePrt/yPrt,2));
		
		grBAPt->SetPoint(j, ptNuc/fA, bA);
		grBAPt->SetPointError(j++, errPt, errBA);
	}
	
	return grBAPt;
}

Double_t AliLnB2::GetErrorY(const TGraphErrors* gr, Double_t x0) const
{
//
// estimate error of gr(x0) with the closest point to x0
//
	const Double_t kEpsilon  = 1.e-6;
	const Double_t kUnknownError = 1.e+6;
	
	for(Int_t i=0; i<gr->GetN(); ++i)
	{
		Double_t x = gr->GetX()[i];
		
		if( TMath::Abs(x-x0) < kEpsilon) return gr->GetErrorY(i);
		
		if( ((i == 0) && (x > x0)) || ((i == (gr->GetN()-1)) && (x < x0)) )
		{
			this->Warning("GetErrorY", "%f is out of bounds",x);
			return kUnknownError;
		}
		
		if( x > x0 )
		{
			this->Warning("GetErrorY", "Interpolating error at %f",x0);
			return (gr->GetErrorY(i)+gr->GetErrorY(i-1))/2;
		}
	}
	
	return 0;
}

Double_t AliLnB2::Rside2Rlong(Double_t pt, Double_t B2, Double_t Cd) const
{
//
// Rside^2*Rlong from B2 value
// Phys. Rev. C, vol. 59, no. 3, pp. 1585--1602, 1999
// (for pp ignore the exponential term)
//
	Double_t hbarc = 0.1973269631; // GeV*fm
	Double_t m = 0.938272013;  // GeV/c^2
	Double_t pi = TMath::Pi();
	
	if(B2==0) return 0;
	
	Double_t mt = TMath::Sqrt(pt*pt + m*m);
	Double_t r3 = 3.*TMath::Power(pi,3./2.)*TMath::Power(hbarc,3.)*Cd/(2.*mt*B2);
	
	return r3;
}

TGraphErrors* AliLnB2::Rside2Rlong(const TGraphErrors* grB2, const TString& name, Double_t Cd) const
{
//
// Rside^2*Rlong from B2 value
// Phys. Rev. C, vol. 59, no. 3, pp. 1585--1602, 1999
// (for pp ignore the exponential term)
//
	TGraphErrors* grR3 = new TGraphErrors(*grB2);
	grR3->SetName(name.Data());
	
	Double_t pt, b2;
	
	for(Int_t i=0; i < grB2->GetN(); ++i)
	{
		grB2->GetPoint(i, pt, b2);
		
		if(b2==0) continue;
		Double_t r3 = Rside2Rlong(pt, b2, Cd);
		
		grR3->SetPoint(i, pt, r3);
		
		// only statistical error propagation of B2
		Double_t errPt = grB2->GetErrorX(i);
		Double_t errB2 = grB2->GetErrorY(i);
		Double_t errR = r3*errB2/b2;
		
		grR3->SetPointError(i, errPt, errR);
	}
	
	return grR3;
}
