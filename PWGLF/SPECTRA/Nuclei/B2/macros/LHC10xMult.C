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

// LHC10x multiplicity config for protons and deuterons
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <TFileMerger.h>
#include "AliLnDriver.h"
#endif

#include "Config.h"

Int_t LHC10xMult(  const TString& species   = "Deuteron"
                 , const TString& inputDir  = "~/alice/input"
                 , const TString& outputDir = "~/alice/output"
                 , const TString& period    = "lhc10d"
                 , const TString& otag      = "lhc10d"
                 , const TString& trkselTag = "-tpc3-nsd-moc"
                 , Double_t       ymax      = 0.5
                 , Int_t          option    = 2)
{
//
// lhc10x multiplicity config
// call Config_XXX for each multiplicity class
//
// if option = 0 then use Config_XXX_TPC
// if option = 1 then use Config_XXX_TOF
// if option = 2 then use Config_TPCTOF
//
	using namespace B2mult;
	using namespace std;
	
	const Bool_t  kINEL[kNmult] = { 0 };
	
	Double_t ptmin   = (species=="Proton") ? 0.4 : 0.7;
	Double_t ptjoint = (species=="Proton") ? 1.0 : 1.0;
	Double_t ptmax   = (species=="Proton") ? 2.0 : 2.5;
	Double_t ptpid   = (species=="Proton") ? 0.4 : 1.3;
	
	if( (option<0) || (option>2) || ((species != "Proton") && (species != "Deuteron")))
	{
		cerr << "unknown species/option: " << species << "/" << option << endl;
		cerr << "species: Proton or Deuteron, options: 0 (TPC), 1 (TOF), 2 (TPCTOF)" << endl;
		exit(1);
	}
	
	TFileMerger m1,m2;
	
	TString ratio[kNmult];
	TString spectra[kNmult];
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		// limit the ptmax value for deuterons
		if(species == "Deuteron") ptmax = kDtrPtMax[i];
		
		cout << endl;
		cout << "Multiplicity class : " << kMultTag[i] << endl;
		cout << "Period             : " << period << endl;
		
		TString outputTag = otag + kMultTag[i];
		
		TString arg =          inputDir    + "\","
			      + "\"" + outputDir   + "\","
			      + "\"" + period      + "\","
			      + "\"" + outputTag   + "\","
		              + "\"" + trkselTag   + "\","
			      + "\"" + kMultTag[i] + "\"," // data
			      + "\"" + "";          // same simulations for all mult
			
		switch(option)
		{
			case 0:
				cout << "Config_" << species << "_TPC_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_%s_TPC_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f)", species.Data(), arg.Data(), ymax, kINEL[i], ptmin, ptmax));
				break;
			case 1:
				cout << "Config_" << species << "_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_%s_TOF_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f, %f)", species.Data(), arg.Data(), ymax, kINEL[i], ptmin, ptmax, ptpid));
				break;
			case 2:
				cout << "Config_TPCTOF_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_TPCTOF_LHC10x.C+g(\"%s\", %f, %d, 0, \"%s\", %f, %f, %f, %f)", arg.Data(), ymax, kINEL[i], species.Data(), ptmin, ptjoint, ptmax, ptpid));
				break;
		}
		
		ratio[i]   = outputDir + "/" + species + "-" + outputTag + "-Ratio.root";
		spectra[i] = outputDir + "/" + species + "-" + outputTag + "-Spectra.root";
		
		m1.AddFile(ratio[i].Data(),0);
		m2.AddFile(spectra[i].Data(),0);
	}
	
	// merge
	
	TString allRatios  = outputDir + "/" + species + "-" + otag + "-Mult-Ratio.root";
	TString allSpectra = outputDir + "/" + species + "-" + otag + "-Mult-Spectra.root";
	
	m1.OutputFile(allRatios.Data());
	m2.OutputFile(allSpectra.Data());
	
	m1.Merge();
	m2.Merge();
	
	// delete tmp files
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		gSystem->Exec(Form("rm -f %s", ratio[i].Data()));
		gSystem->Exec(Form("rm -f %s", spectra[i].Data()));
	}
	
	// draw output
	
	gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"Anti%s%s_Ratio_Pt\",\"\",0,4.5, 0., 1.8, \"#it{p}_{T} (GeV/c)\", \"#bar{p}/p\", 0, \"cRatio\",\"Particle ratio\")", allRatios.Data(),species.Data(),species.Data()));
	
	Double_t minYield = (species=="Proton") ? 1.1e-6 : 1.1e-8;
	Double_t maxYield = (species=="Proton") ? 4.e-1  : 7.e-4;
	
	DrawOutputSpectraMult(allSpectra, species, minYield, maxYield);
	
	return 0;
}
