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

// B2 as a function of multiplicity
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFileMerger.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>

#include "AliLnB2.h"
#include "B2.h"
#include "Config.h"

Double_t GetCd(Double_t z)
{
//
// parameterization of <Cd> as a function of multiplicity
// from ALICE Rlong and Rside measurements
//
	return 0.046133 + 0.0484458*z;
}

Int_t B2Mult(  const TString& pSpectra     = "~/alice/output/Proton-lhc10d-nsd-Mult-Spectra.root"
             , const TString& ptag         = "lhc10d-nsd"
             , const TString& dSpectra     = "~/alice/output/Deuteron-lhc10bcde-nsd-Mult-Spectra.root"
             , const TString& dtag         = "lhc10bcde-nsd"
             , const TString& outputMultPt = "~/alice/output/B2-Mult-Pt.root"
             , const TString& outputPtMult = "~/alice/output/B2-Pt-Mult.root"
             , const TString& otag         = "pp-nsd")
{
//
// B2 as a function of multiplicity
//
	using namespace B2mult;
	
	const Int_t kNpart = 2;
	
	const Int_t kNMinPt = 0;
	const Int_t kNMaxPt = 6;
	
	const TString kPrefix[] = { "", "Anti"};
	const TString kSuffix[] = { "", "bar" };
	Int_t kCharge[]         = {1, -1};
	
	// B2 as a function of pt for each multiplicity class
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		TFileMerger m;
		
		for(Int_t j=0; j<kNpart; ++j)
		{
			TString b2file = kPrefix[j] + "B2.root";
			
			AliLnB2 b2(pSpectra, ptag + "-" + kMultTag[i], dSpectra, dtag + "-" + kMultTag[i], b2file, otag + "-" + kMultTag[i], 2, kCharge[j]);
			
			b2.SetCd(GetCd(kKNOmult[i]));
			
			b2.Run();
			
			m.AddFile(b2file.Data(),0);
		}
		
		// merge B2 and B2bar
		
		TString outputfile = otag + "-" + kMultTag[i] + "-B2.root";
		
		m.OutputFile(outputfile.Data());
		m.Merge();
		
		gSystem->Exec("rm -f B2.root AntiB2.root");
	}
	
	// merge multiplicity classes
	
	TFileMerger m;
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		TString b2 = otag + "-" + kMultTag[i] + "-B2.root";
		m.AddFile(b2.Data(),0);
	}
	
	m.OutputFile(outputMultPt.Data());
	m.Merge();
	
	// delete tmp files
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		gSystem->Exec(Form("rm -f %s-%s-B2.root",otag.Data(),kMultTag[i].Data()));
	}
	
	// B2 as a function of multiplicity for each pt
	
	TFile* finput = new TFile(outputMultPt.Data());
	if (finput->IsZombie()) exit(1);
	
	TGraphErrors* grB2pt[kNpart][kNmult];
	TGraphErrors* grR3pt[kNpart][kNmult];
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		for(Int_t j=0; j<kNmult; ++j)
		{
			grB2pt[i][j] = (TGraphErrors*)FindObj(finput, otag + "-" + kMultTag[j], Form("B2%s_Pt", kSuffix[i].Data()));
			grR3pt[i][j] = (TGraphErrors*)FindObj(finput, otag + "-" + kMultTag[j], Form("R3%s_Pt", kSuffix[i].Data()));
		}
	}
	
	TFile* foutput = new TFile(outputPtMult.Data(),"recreate");
	
	Double_t* pt = grB2pt[0][0]->GetX();
	TString ptLabel[kNMaxPt];
	
	if(kNMaxPt > grB2pt[0][0]->GetN())
	{
		std::cerr << "max pt too big" << std::endl;
		exit(1);
	}
	
	for(Int_t i=kNMinPt; i<kNMaxPt; ++i)
	{
		ptLabel[i] = Form("pT%.02fA",pt[i]);
		foutput->mkdir(ptLabel[i].Data());
	}
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		for(Int_t j=kNMinPt; j<kNMaxPt; ++j)
		{
			Double_t B2[kNmult];
			Double_t B2StatErr[kNmult];
			
			Double_t R3[kNmult];
			Double_t R3StatErr[kNmult];
			
			for(Int_t k=0; k<kNmult; ++k)
			{
				Double_t x, y, ey;
				grB2pt[i][k]->GetPoint(j,x,y);
				ey = grB2pt[i][k]->GetErrorY(j);
				
				B2[k] = y;
				B2StatErr[k] = ey;
				
				grR3pt[i][k]->GetPoint(j,x,y);
				ey = grR3pt[i][k]->GetErrorY(j);
				
				R3[k] = y;
				R3StatErr[k] = ey;
			}
			
			TGraphErrors* grB2Mult = new TGraphErrors(kNmult, kKNOmult, B2, kKNOmultErr, B2StatErr);
			grB2Mult->SetName(Form("B2%s_Zmult", kSuffix[i].Data()));
			
			TGraphErrors* grR3Mult = new TGraphErrors(kNmult, kKNOmult, R3, kKNOmultErr, R3StatErr);
			grR3Mult->SetName(Form("R3%s_Zmult", kSuffix[i].Data()));
			
			foutput->cd(ptLabel[j].Data());
			
			grB2Mult->Write();
			grR3Mult->Write();
			
			delete grB2Mult;
			delete grR3Mult;
		}
	}
	
	delete foutput;
	delete finput;
	
	// particle ratios
	
	gROOT->ProcessLine(Form(".x RatioMult.C+g(\"%s\",\"%s\",\"%s\",\"%s\")", pSpectra.Data(), dSpectra.Data(), ptag.Data(), dtag.Data()));
	
	// draw B2 as a function of pt
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"B2%s_Pt\",\"\",0,2, 1.e-3, 7.e-2,\"p_{T}/A (GeV/c)\",\"B_{2} (GeV^{2}/c^{3})\", 0,\"c%d.B2pt\",\"B2%spt\")", outputMultPt.Data(), kSuffix[i].Data(), i, kSuffix[i].Data()));
		
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"R3%s_Pt\",\"\",0,2, 0, 1.7,\"p_{T}/A (GeV/c)\",\"R_{side}^{2} R_{long} (fm^{3})\", 0,\"c%d.R3pt\",\"R3%spt\")", outputMultPt.Data(), kSuffix[i].Data(), i, kSuffix[i].Data()));
	}
	
	// draw B2 as a function of z
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"B2%s_Zmult\",\"\",0,5, 3.e-3, 6.e-2,\"z\",\"B_{2} (GeV^{2}/c^{3})\", 0,\"c%d.B2z\",\"B2%sZ\")", outputPtMult.Data(), kSuffix[i].Data(), i, kSuffix[i].Data()));
		
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"R3%s_Zmult\",\"\",0,5, 0, 4,\"z\",\"R_{side}^{2} R_{long} (fm^{3})\", 0,\"c%d.R3z\",\"R3%sZ\")", outputPtMult.Data(), kSuffix[i].Data(), i, kSuffix[i].Data()));
	}
	
	return 0;
}
