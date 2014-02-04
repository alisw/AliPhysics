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

// call Config_XXX for each period
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFileMerger.h>
#include "AliLnDriver.h"
#endif

#include "Config.h"

Int_t LHC10bcde(  const TString& species    = "Deuteron"
                , const TString& inputDir   = "~/alice/input"
                , const TString& outputDir  = "~/alice/output"
                , const TString& outputTag  = "lhc10bcde"
                , const TString& trkselTag  = "-tpc3-nsd-moc"
                , const TString& multTag    = ""
                , const TString& multCorTag = ""
                , Double_t ymax             = 0.5
                , Bool_t inel               = 0  // for mult
                , Bool_t drawOutput         = 1  // for batch
                , Double_t ptmin            = 0.8
                , Double_t ptmax            = 3.2
                , Double_t ptjoint          = 1.0
                , Double_t ptpid            = 1.4
                , Int_t option              = 2)
{
//
// call Config_XXX for each period, merge the corrected pt and then get the results
//
	const Int_t kNper = 4;
	const TString kPeriod[kNper]    = { "lhc10b", "lhc10c", "lhc10d", "lhc10e" };
	const TString kOutputTag[kNper] = { "lhc10b", "lhc10c", "lhc10d", "lhc10e" };
	
	using namespace std;
	
	if( (option<0) || (option>2))
	{
		cerr << "unknown option: " << option << endl;
		cerr << "valid options : 0 (TPC), 1 (TOF), 2 (TPCTOF)" << endl;
		exit(1);
	}
	
	TFileMerger m;
	
	for(Int_t i=0; i<kNper; ++i)
	{
		cout << endl << "Period: " << kPeriod[i] << endl;
		
		TString arg =           inputDir       + "\","
			      + "\""  + outputDir      + "\","
			      + "\""  + kPeriod[i]     + "\","
			      + "\""  + kOutputTag[i]  + "\","
		              + "\""  + trkselTag      + "\","
			      + "\""  + multTag        + "\","
			      + "\""  + multCorTag;
			
		switch(option)
		{
			case 0:
				cout << "Config_" << species << "_TPC_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_%s_TPC_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f,1,1,1,1,1)", species.Data(), arg.Data(), ymax, inel, ptmin, ptmax));
				break;
			case 1:
				cout << "Config_" << species << "_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_%s_TOF_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f, %f,1,1,1,1,1)", species.Data(), arg.Data(), ymax, inel, ptmin, ptmax, ptpid));
				break;
			case 2:
				cout << "Config_TPCTOF_LHC10x.C" << endl << endl;
				gROOT->ProcessLine(Form(".x Config_TPCTOF_LHC10x.C+g(\"%s\", %f, %d, 0, \"%s\", %f, %f, %f, %f)", arg.Data(), ymax, inel, species.Data(), ptmin, ptjoint, ptmax, ptpid));
				break;
		}
		
		TString ptfile = outputDir + "/" + species + "-" + kOutputTag[i] + "-Pt.root";
		m.AddFile(ptfile,0);
	}
	
	TString outputPt      = outputDir + "/" + species + "-" + outputTag + multTag + "-Pt.root";
	TString outputRatio   = outputDir + "/" + species + "-" + outputTag + multTag + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + species + "-" + outputTag + multTag + "-Spectra.root";
	
	// pt
	
	m.OutputFile(outputPt.Data());
	m.Merge();
	
	// spectra
	
	AliLnDriver driver;
	
	driver.SetSpecies(species);
	
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	driver.SetOutputTag(outputTag);
	
	driver.SetRapidityInterval(-ymax,ymax);
	driver.SetExtrapolateToINEL(inel);
	
	driver.SetMakeCorrections(0);
	driver.SetMakePt(0);
	driver.SetMakeRatio(1);
	driver.SetMakeSpectra(1);
	
	driver.Run();
	
	// merge all results for comparison
	
	TFileMerger m2, m3;
	
	for(Int_t i=0; i<kNper; ++i)
	{
		m2.AddFile((outputDir + "/" + species + "-" + kPeriod[i] + "-Ratio.root").Data(),0);
		m3.AddFile((outputDir + "/" + species + "-" + kPeriod[i] + "-Spectra.root").Data(),0);
	}
	
	m2.AddFile(outputRatio.Data(),0);
	m3.AddFile(outputSpectra.Data(),0);
	
	TString allRatios = outputDir + "/" + species + "-" + outputTag + "-2" + multTag + "-Ratio.root";
	TString allSpectra = outputDir + "/" + species + "-" + outputTag + "-2" + multTag + "-Spectra.root";
	
	m2.OutputFile(allRatios.Data());
	m3.OutputFile(allSpectra.Data());
	
	m2.Merge();
	m3.Merge();
	
	// compare periods
	
	if(!drawOutput) return 0;
	
	gROOT->ProcessLine(Form(".x DrawDir.C+g(\"%s\",\"Anti%s%s_Ratio_Pt\",\"%s\",0,4.5, 0., 1.8, \"p_{T} (GeV/c)\", \"neg/pos\", 2, \"cRatio\",\"Particle ratio\")", allRatios.Data(), species.Data(), species.Data(),outputTag.Data()));
	
	Double_t minYield  = (species=="Proton") ? 1.1e-6 : 2.e-8;
	Double_t maxYield  = (species=="Proton") ? 4.e-1  : 9.e-4;
	
	DrawOutputSpectraMult(allSpectra, species, minYield, maxYield, 2, outputTag);
	
	return 0;
}
