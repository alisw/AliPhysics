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

// TPC+TOF config
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TString.h>
#include <TFileMerger.h>
#include "AliLnDriver.h"
#endif

#include "Config.h"

Int_t Config_TPCTOF_LHC10x(  const TString& inputDir   = "~/alice/input"
                           , const TString& outputDir  = "~/alice/output"
                           , const TString& period     = "lhc10d"
                           , const TString& outputTag  = "lhc10d"
                           , const TString& trkselTag  = "-bayes-nsd"
                           , const TString& multTag    = ""
                           , const TString& multCorTag = ""
                           , Double_t ymax             = 0.5
                           , Bool_t inel               = 0
                           , Bool_t drawOutput         = 1
                           , const TString& species    = "Proton"
                           , Double_t ptmin            = 0.4
                           , Double_t ptjoint          = 1.0
                           , Double_t ptmax            = 3.2
                           , Double_t ptpid            = 4)
{
//
// combine TPC and TOF for protons and deuterons
//
	using namespace std;
	
	if((species != "Proton") && (species != "Deuteron"))
	{
		cerr << "Particle species " << species << " not implemented, only 'Proton' and 'Deuteron'." << endl;
		exit(1);
	}
	
	const TString kOutputTagTPC = outputTag + "-tpc";
	const TString kOutputTagTOF = outputTag + "-tof";
	
	const TString kArgTPC =        inputDir        + "\","
		              + "\"" + outputDir       + "\","
		              + "\"" + period          + "\","
		              + "\"" + kOutputTagTPC   + "\","
		              + "\"" + trkselTag       + "\","
		              + "\"" + multTag         + "\","
		              + "\"" + multCorTag;
		
	const TString kArgTOF =        inputDir        + "\","
		              + "\"" + outputDir       + "\","
		              + "\"" + period          + "\","
		              + "\"" + kOutputTagTOF   + "\","
		              + "\"" + trkselTag       + "\","
		              + "\"" + multTag         + "\","
		              + "\"" + multCorTag;
	
	cout << "Config_" << species << "_TPC_LHC10x.C" << endl << endl;
	gROOT->ProcessLine(Form(".x Config_%s_TPC_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f, 0,1,1,0,0)"
				, species.Data()
				, kArgTPC.Data()
				, ymax
				, inel
				, ptmin
				, ptjoint));
		
		
	cout << "Config_" << species << "_TOF_LHC10x.C" << endl << endl;
	gROOT->ProcessLine(Form(".x Config_%s_TOF_LHC10x.C+g(\"%s\", %f, %d, 0, %f, %f, %f, 1,1,1,0,0)"
				, species.Data()
				, kArgTOF.Data()
				, ymax
				, inel
				, ptjoint
				, ptmax
				, ptpid));
		
	TString outputPtTPC = outputDir + "/" + MakeOutputName(species, kOutputTagTPC) + "-Pt.root";
	TString outputPtTOF = outputDir + "/" + MakeOutputName(species, kOutputTagTOF) + "-Pt.root";
	
	TString outputPtTPCdbg = outputDir + "/" + MakeOutputName(species, kOutputTagTPC) + "-Pt-debug.root";
	TString outputPtTOFdbg = outputDir + "/" + MakeOutputName(species, kOutputTagTOF) + "-Pt-debug.root";
	
	// combine TPC and TOF pt
	
	TString outputPt      = outputDir + "/" + MakeOutputName(species, outputTag) + "-Pt.root";
	TString outputRatio   = outputDir + "/" + MakeOutputName(species, outputTag) + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + MakeOutputName(species, outputTag) + "-Spectra.root";
	
	TFileMerger m;
	
	m.AddFile(outputPtTPC.Data(),0);
	m.AddFile(outputPtTOF.Data(),0);
	
	m.OutputFile(outputPt.Data());
	
	m.Merge();
	
	// remove tmp files
	gSystem->Exec(Form("rm -f %s %s %s %s", outputPtTPC.Data(), outputPtTOF.Data(), outputPtTPCdbg.Data(), outputPtTOFdbg.Data()));
	
	// make ratio and spectra
	
	AliLnDriver driver;
	
	driver.SetSpecies(species);
	
	driver.SetRapidityInterval(-ymax,ymax);
	
	driver.SetExtrapolateToINEL(inel);
	
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	
	driver.SetOutputTag(outputTag);
	
	driver.SetMakeCorrections(0);
	driver.SetMakePt(0);
	driver.SetMakeRatio(1);
	driver.SetMakeSpectra(1);
	
	driver.Run();
	
	if(!drawOutput) return 0;
	
	DrawOutputRatio(outputRatio, outputTag, species);
	DrawOutputSpectra(outputSpectra, outputTag, species);
	
	return 0;
}
