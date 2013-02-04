#ifndef CONFIG_H
#define CONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// functions for configs
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <cstdlib>

namespace B2mult
{
//
// multiplicity classes
//
	const Int_t kNmult = 6;
	const TString kMultClass[kNmult]   = { "ntrk0002", "ntrk0204", "ntrk0408", "ntrk0811", "ntrk1120", "ntrk20xx" };
	const Double_t kKNOmult[kNmult]    = { 0.20, 0.60, 1.01, 1.60, 2.60, 4.35 };
	const Double_t kKNOmultErr[kNmult] = { 0, 0, 0, 0, 0, 0 };
};

TString GetCollSystem(const TString& period)
{
//
// translate period name into colliding system name
//
	TString name = period;
	name.ToLower();
	
	if(name == "lhc10c900")    return "pp0.9TeV";
	if(name == "lhc10b_pass2") return "pp7TeV";
	if(name == "lhc10c_pass2") return "pp7TeV";
	if(name == "lhc10b")       return "pp7TeV";
	if(name == "lhc10c")       return "pp7TeV";
	if(name == "lhc10d")       return "pp7TeV";
	if(name == "lhc10e")       return "pp7TeV";
	if(name.Contains("lhc10")) return "pp7TeV";
	if(name == "lhc11a_wsdd")  return "pp2.76TeV";
	if(name == "lhc11a_wosdd") return "pp2.76TeV";
	if(name.Contains("900"))   return "pp0.9TeV";
	if(name.Contains("7000"))  return "pp7TeV";
	if(name.Contains("2760"))  return "pp2.76TeV";
	
	return "unknown";
}

void GetInelXSection(Double_t xsection[3], const TString& period)
{
//
// inelastic cross section in mb and largest stat. and syst error
// measured by ALICE for the given colliding system
// http://arxiv.org/abs/1208.4968
//
	TString collsystem = GetCollSystem(period);
	
	if( collsystem == "pp0.9TeV" )
	{
		xsection[0] = 50.3;
		xsection[1] = 0.4;
		xsection[2] = 1.0;
	}
	else if( collsystem == "pp2.76TeV" )
	{
		xsection[0] = 62.8;
		xsection[1] = 4.0;
		xsection[2] = 4.6;
	}
	else if( collsystem == "pp7TeV" )
	{
		xsection[0] = 73.2;
		xsection[1] = 1.2;
		xsection[2] = 2.6;
	}
	else
	{
		std::cerr << "Warning: unknown colliding system " << collsystem << " for period " << period << std::endl;
	}
}

void GetTriggerEfficiency(Double_t trigeff[3], const TString& trigname, const TString& period)
{
//
// trigger efficiency for the given colliding system and largest stat. and syst. errors
// http://arxiv.org/abs/1208.4968
//
	TString collsystem = GetCollSystem(period);
	TString trigger = trigname;
	trigger.ToLower();
	
	if( trigger == "mbor"  && collsystem == "pp0.9TeV" )
	{
		trigeff[0] = 0.910;
		trigeff[1] = 0; // stat. error
		trigeff[2] = 0.032; // syst. error
	}
	else if( trigger == "mbor"  && collsystem == "pp2.76TeV")
	{
		trigeff[0] = 0.881;
		trigeff[1] = 0;
		trigeff[2] = 0.059;
	}
	else if( trigger == "mbor"  && collsystem == "pp7TeV")
	{
		trigeff[0] = 0.852;
		trigeff[1] = 0;
		trigeff[2] = 0.062;
	}
	else if( trigger == "mband" && collsystem == "pp0.9TeV" )
	{
		trigeff[0] = 0.763;
		trigeff[1] = 0;
		trigeff[2] = 0.022;
	}
	else if( trigger == "mband" && collsystem == "pp2.76TeV")
	{
		trigeff[0] = 0.760;
		trigeff[1] = 0;
		trigeff[2] = 0.052;
	}
	else if( trigger == "mband" && collsystem == "pp7TeV")
	{
		trigeff[0] = 0.742;
		trigeff[1] = 0;
		trigeff[2] = 0.050;
	}
	else
	{
		std::cerr << "Warning: unknown trigger/colliding system " << trigname << "/" << collsystem << " for period " << period << std::endl;
	}
}

TString GetSimuPeriod(const TString& period)
{
//
// simulation code for the given period
//
	TString name = period;
	name.ToLower();
	
	if(period=="lhc10c900")   return "lhc10e13";
	if(period=="lhc10b")      return "lhc10d1";
	if(period=="lhc10c")      return "lhc10d4";
	if(period=="lhc10d")      return "lhc10f6a";
	if(period=="lhc10e")      return "lhc10e21";
	if(period=="lhc11a_wsdd") return "lhc11e3a_wsdd";
	
	return "";
}

TString GetSimuFixPeriod(const TString& period)
{
//
// simulation fix code for the given period
//
	TString name = period;
	name.ToLower();
	
	if(period=="lhc10c900")   return "lhc12a5a";
	if(period=="lhc10b")      return "lhc12a5bb";
	if(period=="lhc10c")      return "lhc12a5bc";
	if(period=="lhc10d")      return "lhc12a5bd";
	if(period=="lhc10e")      return "lhc12a5be";
	if(period=="lhc10cde")    return "lhc12a5bcde";
	if(period=="lhc10bcde")   return "lhc12a5bbcde";
	if(period=="lhc11a_wsdd") return "lhc12a5c_wsdd";
	
	return "";
}

TString MakeInputName(const TString& species, const TString& period, const TString& trksel)
{
//
// make input name for data
//
	return species + "-" + period + "-" + trksel;
}

TString MakeSimuName(const TString& species, const TString& period, const TString& trksel)
{
//
// make input name for simulation
//
	TString simu = (species == "Proton") ? GetSimuPeriod(period) : GetSimuFixPeriod(period);
	
	if(simu=="")
	{
		std::cerr << "No simulation for period: " << period << std::endl;
		exit(1);
	}
	
	return species + "-" + simu + "-" + trksel;
}

TString MakeSimuFixName(const TString& species, const TString& period, const TString& trksel, Bool_t g3Fluka=0)
{
//
// make input name for simulation fix
//
	TString simufix = (species == "Proton" && g3Fluka) ? GetSimuPeriod(period) : GetSimuFixPeriod(period);
	
	if(simufix=="")
	{
		std::cerr << "No simulation fix for period: " << period << std::endl;
		exit(1);
	}
	
	return species + "-" + simufix + "-" + trksel;
}

TString MakeOutputName(const TString& species, const TString& outputTag)
{
//
// make output name
//
	return species + "-" + outputTag;
}

TStyle* GetDrawingStyle()
{
//
// define a default style for drawing
//
	TStyle* st = new TStyle();
	
	//st->SetPalette(51);
	st->SetPalette(1);
	st->SetPadTickX(1);
	st->SetPadTickY(1);
	st->SetPadGridX(1);
	st->SetPadGridY(1);
	
	st->SetCanvasColor(0);
	st->SetFrameBorderMode(0);
	st->SetStatBorderSize(1);
	st->SetStatColor(0);
	st->SetFrameFillColor(0);
	st->SetTitleFillColor(0);
	st->SetLabelFont(62,"XYZ");
	st->SetTitleFont(62,"XYZ");
	
	//st->SetOptTitle(1);
	st->SetOptStat(0);
	
	return st;
}

void DrawOutputCorr(const TString& species, const TString& corr, const TString& tag="")
{
//
// draw corrections for the given species
//
	gROOT->ProcessLine(Form(".x DrawCorr.C+g(\"%s\",\"%s\",\"%s\")", species.Data(), corr.Data(), tag.Data()));
}

void DrawCorrDebug(const TString& sec, const TString& tag, const TString& species, Int_t lowbin=1, Int_t hibin=10, Double_t dcaxyMin=-1.5, Double_t dcaxyMax=1.5)
{
//
// draw correction for secondaries
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawSec.C+g(\"%s\",\"%s\",\"%s\", %d, %d, %f, %f)", sec.Data(), tag.Data(), kParticle[i].Data(), lowbin, hibin, dcaxyMin, dcaxyMax));
	}
}

void DrawPtDebug(const TString& pt, const TString& tag, const TString& species, Bool_t m2pid=0, Int_t hiptbin=17, Int_t lowm2bin=9, Int_t him2bin=17)
{
//
// draw pt debug for the particle species
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawPt.C+g(\"%s\",\"%s\",\"%s\",%d, %d, %d, %d)", pt.Data(), tag.Data(), kParticle[i].Data(), hiptbin, m2pid, lowm2bin, him2bin));
	}
}

void DrawOutputRatio(const TString& ratio, const TString& tag, const TString& species)
{
//
// draw ratio antiparticle/particle
//
	gROOT->ProcessLine(Form(".x DrawRatio.C+g(\"%s\",\"%s\",\"%s\")", ratio.Data(), tag.Data(), species.Data()));
}

void DrawOutputSpectra(const TString& spectra, const TString& tag, const TString& species)
{
//
// draw spectra for the particle species
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawSpectra.C+g(\"%s\",\"%s\",\"%s\")", spectra.Data(), tag.Data(), kParticle[i].Data()));
	}
}

void DrawOutputSpectraMult(const TString& spectra, const TString& species, Double_t ymin, Double_t ymax, Int_t option=0, const TString& refdir="")
{
//
// draw spectra for each multiplicity class
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_InvDiffYield_Pt\",\"%s\",0,4.5, %g, %g, \"p_{T} (GeV/c)\", \"#frac{1}{2#piN_{inel}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c^{3})\", %d, \"c%d\",\"%s\")", spectra.Data(), kParticle[i].Data(), refdir.Data(),ymin, ymax, option, i, kParticle[i].Data()));
	}
}

#endif // CONFIG_H
