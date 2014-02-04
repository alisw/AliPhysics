#ifndef CONFIG_H
#define CONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// common functions for config macros
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <cstdlib>
#endif

namespace CollSystem
{
//
// pp collision energies
//
	const Int_t    kNener               = 3;
	const TString  kEnergyTag[kNener]   = { "-900GeV", "-2.76TeV", "-7TeV" };
	const TString  kEnergyLabel[kNener] = { "900 GeV", "2.76 TeV", "7 TeV" };
	const Double_t kEnergy[kNener]      = { 0.9, 2.76, 7 };
	const Double_t kEnergyError[kNener] = { 0 };
	const TString  kEnergyName          = "#sqrt{s} (TeV)";
};

namespace B2mult
{
//
// multiplicity classes
//

	const Int_t    kNmult              = 5;
	const TString  kMultTag[kNmult]    = { "-ntrk0103", "-ntrk0405", "-ntrk0608", "-ntrk0914", "-ntrk15xx" };
	const Double_t kKNOmult[kNmult]    = { 0.35, 0.74, 1.13, 1.82, 3.14 };
	const Double_t kDtrPtMax[kNmult]   = { 1.8, 2.0, 2.2, 2.3, 2.5 };
	const Double_t kKNOmultErr[kNmult] = { 0 };
	const TString  kKNOmultName        = "z";
};

namespace B2HiLowMult
{
//
// high and low multiplicity classes
//
	const Int_t    kNmult              = 2;
	const TString  kMultTag[kNmult]    = { "-ntrk0108", "-ntrk09xx" };
	const Double_t kKNOmult[kNmult]    = { 0.64, 2.23 };
	const Double_t kDtrPtMax[kNmult]   = { 2.0, 2.5 };
	const Double_t kKNOmultErr[kNmult] = { 0 };
	const TString  kKNOmultName        = "z";
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
// trigger efficiency for the given colliding system and largest systematic error
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
	if(period=="lhc10de")     return "lhc10f6ade";
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
	if(period=="lhc10bc")     return "lhc12a5bbc";
	if(period=="lhc10de")     return "lhc12a5bde";
	if(period=="lhc10cde")    return "lhc12a5bcde";
	if(period=="lhc10bcd")   return "lhc12a5bbcd";
	if(period=="lhc10bcde")   return "lhc12a5bbcde";
	if(period=="lhc11a_wsdd") return "lhc12a5c_wsdd";
	
	return "";
}

TString MakeInputName(const TString& species, const TString& period, const TString& trksel)
{
//
// make input name for data
//
	if(trksel == "")           return species + "-" + period;
	if(trksel.BeginsWith("-")) return species + "-" + period + trksel;
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

void DrawOutputCorr(const TString& species, const TString& corr, const TString& tag="")
{
//
// draw corrections for the given species
//
	gROOT->ProcessLine(Form(".x DrawCorr.C+g(\"%s\",\"%s\",\"%s\")", species.Data(), corr.Data(), tag.Data()));
}

void DrawCorrDebug(const TString& sec, const TString& tag, const TString& species, Double_t ptmin=0.5, Double_t ptmax=3., Double_t dcaxyMin=-1.5, Double_t dcaxyMax=1.5)
{
//
// draw correction for secondaries
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawSec.C+g(\"%s\",\"%s\",\"%s\", %f, %f, %f, %f)", sec.Data(), tag.Data(), kParticle[i].Data(), ptmin, ptmax, dcaxyMin, dcaxyMax));
	}
}

void DrawPtDebug(const TString& pt, const TString& tag, const TString& species, Bool_t m2pid=0, Double_t ptmax=3., Double_t ptpid=1.2)
{
//
// draw pt debug for the particle species
//
	const TString kParticle[] = { species, Form("Anti%s", species.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawPt.C+g(\"%s\",\"%s\",\"%s\",%f, %d, %f)", pt.Data(), tag.Data(), kParticle[i].Data(), ptmax, m2pid, ptpid));
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
		gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"%s_InvDiffYield_Pt\",\"%s\",0,3.5, %g, %g, \"p_{T} (GeV/c)\", \"#frac{1}{2#piN_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c^{3})\", %d, \"%s\",\"%s\")", spectra.Data(), kParticle[i].Data(), refdir.Data(),ymin, ymax, option, kParticle[i].Data(), kParticle[i].Data()));
	}
}

#endif // CONFIG_H
