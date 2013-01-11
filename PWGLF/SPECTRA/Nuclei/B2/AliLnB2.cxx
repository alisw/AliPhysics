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
	
	TGraphErrors* grPrtInvDYieldPt = (TGraphErrors*)FindObj(finput1, fProtonTag, prefix + "Proton_InvDiffYield_Pt");
	TGraphErrors* grNucInvDYieldPt = (TGraphErrors*)FindObj(finputA, fNucleusTag, prefix + fNucleusName + "_InvDiffYield_Pt");
	
	TGraphErrors* grSysErrPrtInvDYieldPt = (TGraphErrors*)FindObj(finput1, fProtonTag, prefix + "Proton_SysErr_InvDiffYield_Pt");
	TGraphErrors* grSysErrNucInvDYieldPt = (TGraphErrors*)FindObj(finputA, fNucleusTag, prefix + fNucleusName + "_SysErr_InvDiffYield_Pt");
	
	TFile* foutput = new TFile(fOutputFilename.Data(),"recreate");
	
	foutput->mkdir(fOutputTag.Data());
	foutput->cd(fOutputTag.Data());
	
	// coalescence parameter
	
	TGraphErrors* grB2Pt = this->GetBAPt(grPrtInvDYieldPt, grNucInvDYieldPt, Form("B2%s_Pt",suffix.Data()));
	
	TGraphErrors* grSysErrB2Pt = this->GetBAPt(grSysErrPrtInvDYieldPt, grSysErrNucInvDYieldPt, Form("B2%s_SysErr_Pt",suffix.Data()), 0.025);
	
	grB2Pt->Write();
	grSysErrB2Pt->Write();
	
	// homogeneity volume
	
	TGraphErrors* grR3Pt = this->Rside2Rlong(grB2Pt, Form("R3%s_Pt", suffix.Data()), fCd);
	TGraphErrors* grSysErrR3Pt = this->Rside2Rlong(grSysErrB2Pt, Form("R3%s_SysErr_Pt", suffix.Data()), fCd, 0.025);
	
	grR3Pt->Write();
	grSysErrR3Pt->Write();
	
	// clean
	
	delete grB2Pt;
	delete grSysErrB2Pt;
	delete grR3Pt;
	delete grSysErrR3Pt;
	
	delete foutput;
	delete finput1;
	delete finputA;
	
	return 0;
}

TGraphErrors* AliLnB2::GetBAPt(const TGraphErrors* grPrtInvDYieldPt, const TGraphErrors* grNucInvDYieldPt, const TString& name, Double_t errPt) const
{
//
// coalescence parameter
//
	TGraphErrors* grBAPt = new TGraphErrors();
	grBAPt->SetName(name.Data());
	
	Int_t offset = 0;
	Double_t pt, yPrt, ptNuc, yNuc;
	
	grPrtInvDYieldPt->GetPoint(0,pt,yPrt);
	
	// find offset
	for(Int_t i=0; i < grNucInvDYieldPt->GetN(); ++i)
	{
		grNucInvDYieldPt->GetPoint(i, ptNuc, yNuc);
		if(TMath::Abs(fA*pt-ptNuc)<0.001)
		{
			offset = i;
			break;
		}
	}
	
	for(Int_t i=0, j=0; i < grPrtInvDYieldPt->GetN() && i+offset < grNucInvDYieldPt->GetN(); ++i)
	{
		grPrtInvDYieldPt->GetPoint(i,pt,yPrt);
		grNucInvDYieldPt->GetPoint(i+offset,ptNuc, yNuc);
		
		if(yPrt == 0 || yNuc == 0 ) continue;
		if(TMath::Abs(fA*pt-ptNuc)>0.001) continue; // compare at equal momentum per nucleon
		
		Double_t bA = yNuc/TMath::Power(yPrt,fA);
		
		// error
		Double_t ePrt = grPrtInvDYieldPt->GetErrorY(i);
		Double_t eNuc = grNucInvDYieldPt->GetErrorY(i);
		
		Double_t errBA = bA*TMath::Sqrt(TMath::Power(eNuc/yNuc,2) + TMath::Power(fA*ePrt/yPrt,2));
		
		grBAPt->SetPoint(j, pt, bA);
		grBAPt->SetPointError(j++, errPt, errBA);
	}
	
	return grBAPt;
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

TGraphErrors* AliLnB2::Rside2Rlong(const TGraphErrors* grB2, const TString& name, Double_t Cd, Double_t errPt) const
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
		Double_t errB2 = grR3->GetErrorY(i);
		Double_t errR = r3*errB2/b2;
		
		grR3->SetPointError(i, errPt, errR);
	}
	
	return grR3;
}
