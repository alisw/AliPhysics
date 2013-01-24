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

#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <TFileMerger.h>

#include "AliLnDriver.h"
#include "Config.h"

Int_t LHC10xMult(const TString& species   = "Proton",
                 const TString& inputDir  = "~/alice/input",
                 const TString& outputDir = "~/alice/output",
                 const TString& period    = "lhc10d",
                 Int_t option = 2)
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
	
	const Bool_t  kNormToInel[kNmult] = { 1, 0 }; // only normalize first bin
	
	Double_t ymin = (species=="Proton") ? 1.1e-6 : 1.1e-8;
	Double_t ymax = (species=="Proton") ? 4.e-1 : 4.e-4;
	
	TFileMerger m1,m2;
	
	TString ratio[kNmult];
	TString spectra[kNmult];
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		cout << endl;
		cout << "Multiplicity class : " << kMultClass[i] << endl;
		cout << "Period             : " << period << endl;
		
		TString outputTag  = period + "-" + kMultClass[i];
		
		TString arg =          inputDir            + "\","
			      + "\"" + outputDir           + "\","
			      + "\"" + period              + "\","
			      + "\"" + outputTag           + "\","
			      + "\"" + "-" + kMultClass[i] + "\"," // data
			      + "\"" + "";                         // same simulations for all mult
			
		if(species=="Proton" && option==0)
		{
			cout << "Config_Proton_TPC_LHC10x.C" << endl << endl;
			gROOT->ProcessLine(Form(".x Config_Proton_TPC_LHC10x.C+g(\"%s\",%d,0)", arg.Data(), kNormToInel[i]));
		}
		else if(species=="Proton" && option==1)
		{
			cout << "Config_Proton_TOF_LHC10x.C" << endl << endl;
			gROOT->ProcessLine(Form(".x Config_Proton_TOF_LHC10x.C+g(\"%s\",%d,0)", arg.Data(), kNormToInel[i]));
		}
		else if(species=="Deuteron" && option==0)
		{
			cout << "Config_Deuteron_TPC_LHC10x.C" << endl << endl;
			gROOT->ProcessLine(Form(".x Config_Deuteron_TPC_LHC10x.C+g(\"%s\",%d,0)", arg.Data(), kNormToInel[i]));
		}
		else if(species=="Deuteron" && option==1)
		{
			cout << "Config_Deuteron_TOF_LHC10x.C" << endl << endl;
			gROOT->ProcessLine(Form(".x Config_Deuteron_TOF_LHC10x.C+g(\"%s\",%d,0)", arg.Data(), kNormToInel[i]));
		}
		else if((species=="Proton" || species=="Deuteron") && option==2)
		{
			cout << "Config_TPCTOF_LHC10x.C" << endl << endl;
			gROOT->ProcessLine(Form(".x Config_TPCTOF_LHC10x.C+g(\"%s\",%d,0,\"%s\")", arg.Data(), kNormToInel[i], species.Data()));
		}
		else
		{
			cerr << "unknown species/option: " << species << "/" << option << endl;
			cerr << "usage: Proton/0, Proton/1, Proton/2 or Deuteron" << endl;
			exit(1);
		}
		
		ratio[i]   = outputDir + "/" + species + "-" + outputTag + "-Ratio.root";
		spectra[i] = outputDir + "/" + species + "-" + outputTag + "-Spectra.root";
		
		m1.AddFile(ratio[i].Data(),0);
		m2.AddFile(spectra[i].Data(),0);
	}
	
	// merge
	
	TString allRatios  = outputDir + "/" + species + "-" + period + "-Mult-Ratio.root";
	TString allSpectra = outputDir + "/" + species + "-" + period + "-Mult-Spectra.root";
	
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
	
	gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"Anti%s%s_Ratio_Pt\",\"\",0,4.5, 0., 1.8, \"p_{T} (GeV/c)\", \"#bar{p}/p\", 0, \"cRatio\",\"Particle ratio\")", allRatios.Data(),species.Data(),species.Data()));
	
	DrawOutputSpectraMult(allSpectra, species, ymin, ymax);
	
	return 0;
}
