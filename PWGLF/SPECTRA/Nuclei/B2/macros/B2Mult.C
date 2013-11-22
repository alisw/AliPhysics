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

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFileMerger.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include "AliLnBA.h"
#endif

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

extern void RatioMult(  const TString&, const TString&, const TString&, const TString&, const Int_t, const TString*, const Double_t*, const Double_t*, const TString&, const Bool_t, const TString& , const TString& );

Int_t B2Mult(  const TString& pSpectra     = "~/alice/output/Proton-lhc10d-Mult-Spectra.root"
             , const TString& ptag         = "lhc10d"
             , const TString& dSpectra     = "~/alice/output/Deuteron-lhc10d-Mult-Spectra.root"
             , const TString& dtag         = "lhc10d"
             , const TString& outputMultPt = "~/alice/output/B2-lhc10d-MultPt.root"
             , const TString& outputPtMult = "~/alice/output/B2-lhc10d-PtMult.root"
             , const TString& otag         = "lhc10d"
             , const Bool_t qTsallis       = 0
             , const TString& tsallisTag   = "Tsallis"
             , const TString& oratio       = "~/alice/output/Particle-Ratios-lhc10d-Tsallis.root")
{
//
// B2 as a function of multiplicity
//
	using namespace B2mult;
	using namespace std;
	
	const Int_t kNpart = 2;
	
	const Int_t kNMinPt = 0;
	const Int_t kNMaxPt = 6;
	
	const TString kNucleus[kNpart] = { "Deuteron", "AntiDeuteron" };
	
	// B2 as a function of pt for each multiplicity class
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		TFileMerger m;
		
		const Int_t kZ[kNpart] = { 1, -1 };
		const TString kB2File[kNpart] = {"B2.root", "AntiB2.root" };
		
		for(Int_t j=0; j<kNpart; ++j)
		{
			cout << kMultTag[i] << endl;
			
			AliLnBA b2(pSpectra, ptag + kMultTag[i], dSpectra, dtag + kMultTag[i], kB2File[j], otag + kMultTag[i], 2, kZ[j]);
			
			b2.SetCd(GetCd(kKNOmult[i]));
			
			b2.Run();
			
			m.AddFile(kB2File[j].Data(),0);
		}
		
		// merge B2 and B2bar
		
		TString outputfile = otag + kMultTag[i] + "-B2.root";
		
		m.OutputFile(outputfile.Data());
		m.Merge();
		
		gSystem->Exec(Form("rm -f %s %s",kB2File[0].Data(),kB2File[1].Data()));
	}
	
	// merge multiplicity classes
	
	TFileMerger m;
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		TString b2 = otag + kMultTag[i] + "-B2.root";
		m.AddFile(b2.Data(),0);
	}
	
	m.OutputFile(outputMultPt.Data());
	m.Merge();
	
	// delete tmp files
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		gSystem->Exec(Form("rm -f %s%s-B2.root",otag.Data(),kMultTag[i].Data()));
	}
	
	// B2 as a function of multiplicity for each pt
	
	TFile* finput = new TFile(outputMultPt.Data());
	if (finput->IsZombie()) exit(1);
	
	TGraphErrors* grB2pt[kNpart][kNmult];
	TGraphErrors* grR3pt[kNpart][kNmult];
	TGraphErrors* grSysB2pt[kNpart][kNmult];
	TGraphErrors* grSysR3pt[kNpart][kNmult];
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		for(Int_t j=0; j<kNmult; ++j)
		{
			grB2pt[i][j] = FindObj<TGraphErrors>(finput, otag + kMultTag[j], kNucleus[i] + "_B2_Pt");
			grR3pt[i][j] = FindObj<TGraphErrors>(finput, otag + kMultTag[j], kNucleus[i] + "_R3_Pt");
			
			grSysB2pt[i][j] = FindObj<TGraphErrors>(finput, otag + kMultTag[j], kNucleus[i] + "_SystErr_B2_Pt");
			grSysR3pt[i][j] = FindObj<TGraphErrors>(finput, otag + kMultTag[j], kNucleus[i] + "_SystErr_R3_Pt");
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
		ptLabel[i] = Form("pT %.02fA",pt[i]);
		foutput->mkdir(ptLabel[i].Data());
	}
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		for(Int_t j=kNMinPt; j<kNMaxPt; ++j)
		{
			Double_t B2[kNmult];
			Double_t B2StatErr[kNmult];
			Double_t B2SystErr[kNmult];
			
			Double_t R3[kNmult];
			Double_t R3StatErr[kNmult];
			Double_t R3SystErr[kNmult];
			
			for(Int_t k=0; k<kNmult; ++k)
			{
				Double_t x, y, staterr, systerr;
				grB2pt[i][k]->GetPoint(j,x,y);
				staterr = grB2pt[i][k]->GetErrorY(j);
				systerr = grSysB2pt[i][k]->GetErrorY(j);
				
				B2[k] = y;
				B2StatErr[k] = staterr;
				B2SystErr[k] = systerr;
				
				grR3pt[i][k]->GetPoint(j,x,y);
				staterr = grR3pt[i][k]->GetErrorY(j);
				systerr = grSysR3pt[i][k]->GetErrorY(j);
				
				R3[k] = y;
				R3StatErr[k] = staterr;
				R3SystErr[k] = systerr;
			}
			
			TGraphErrors* grB2Mult = new TGraphErrors(kNmult, kKNOmult, B2, kKNOmultErr, B2StatErr);
			grB2Mult->SetName(Form("%s_B2_Zmult", kNucleus[i].Data()));
			
			TGraphErrors* grR3Mult = new TGraphErrors(kNmult, kKNOmult, R3, kKNOmultErr, R3StatErr);
			grR3Mult->SetName(Form("%s_R3_Zmult", kNucleus[i].Data()));
			
			Double_t zMultSystErr[kNmult];
			for(Int_t k=0; k<kNmult; ++k) zMultSystErr[k]=0.07;
			
			TGraphErrors* grSysB2Mult = new TGraphErrors(kNmult, kKNOmult, B2, zMultSystErr, B2SystErr);
			grSysB2Mult->SetName(Form("%s_SystErr_B2_Zmult", kNucleus[i].Data()));
			
			TGraphErrors* grSysR3Mult = new TGraphErrors(kNmult, kKNOmult, R3, zMultSystErr, R3SystErr);
			grSysR3Mult->SetName(Form("%s_SystErr_R3_Zmult", kNucleus[i].Data()));
			
			foutput->cd(ptLabel[j].Data());
			
			grB2Mult->Write();
			grR3Mult->Write();
			grSysB2Mult->Write();
			grSysR3Mult->Write();
			
			delete grB2Mult;
			delete grR3Mult;
			delete grSysB2Mult;
			delete grSysR3Mult;
		}
	}
	
	delete foutput;
	delete finput;
	
	//
	// Particle ratios
	// -----------------------------------
	
	RatioMult(pSpectra, dSpectra, ptag, dtag, kNmult, kMultTag, kKNOmult, kKNOmultErr, kKNOmultName, qTsallis, oratio, tsallisTag);
	
	
	// draw B2 as a function of pt
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_B2_Pt\",\"\",0,2, 1.e-3, 7.e-2,\"p_{T}/A (GeV/c)\",\"B_{2} (GeV^{2}/c^{3})\", 0,\"c%d.B2pt\",\"%s_B2_Pt\")", outputMultPt.Data(), kNucleus[i].Data(), i, kNucleus[i].Data()));
		
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_R3_Pt\",\"\",0,2, 0, 1.7,\"p_{T}/A (GeV/c)\",\"R_{side}^{2} R_{long} (fm^{3})\", 0,\"c%d.R3pt\",\"%s_R3_Pt\")", outputMultPt.Data(), kNucleus[i].Data(), i, kNucleus[i].Data()));
	}
	
	// draw B2 as a function of z
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_B2_Zmult\",\"\",0,5, 3.e-3, 6.e-2,\"z\",\"B_{2} (GeV^{2}/c^{3})\", 0,\"c%d.B2z\",\"%s_B2_Zmult\")", outputPtMult.Data(), kNucleus[i].Data(), i, kNucleus[i].Data()));
		
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_R3_Zmult\",\"\",0,5, 0, 4,\"z\",\"R_{side}^{2} R_{long} (fm^{3})\", 0,\"c%d.R3z\",\"%s_R3_Zmult\")", outputPtMult.Data(), kNucleus[i].Data(), i, kNucleus[i].Data()));
	}
	
	return 0;
}
