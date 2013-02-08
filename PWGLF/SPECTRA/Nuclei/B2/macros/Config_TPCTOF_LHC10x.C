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

#include <Riostream.h>
#include <TSystem.h>
#include <TString.h>
#include <TFileMerger.h>

#include "AliLnDriver.h"
#include "Config.h"

Int_t Config_TPCTOF_LHC10x(const TString& inputDir   = "~/alice/input",
                           const TString& outputDir  = "~/alice/output",
                           const TString& period     = "lhc10d",
                           const TString& outputTag  = "lhc10d",
                           const TString& multTag    = "",
                           const TString& multCorTag = "",
                           Bool_t inel               = 1,  // for mult
                           Bool_t drawOutput         = 1,  // for batch
                           const TString& species    = "Proton",
                           Int_t lowPtBin            = 5,
                           Int_t jointPtBin          = 11,
                           Int_t hiPtBin             = 36)
{
//
// combine TPC and TOF for protons and deuterons
//
	const Double_t kProtonSysErr[2]   = {0.08, 0.08} ;
	const Double_t kDeuteronSysErr[2] = {0.10, 0.11} ;
	
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
		              + "\"" + multTag         + "\","
		              + "\"" + multCorTag;
		
	const TString kArgTOF =        inputDir        + "\","
		              + "\"" + outputDir       + "\","
		              + "\"" + period          + "\","
		              + "\"" + kOutputTagTOF   + "\","
		              + "\"" + multTag         + "\","
		              + "\"" + multCorTag;
	
	cout << "Config_" << species << "_TPC_LHC10x.C" << endl << endl;
	gROOT->ProcessLine(Form(".x Config_%s_TPC_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 0,1,1,0,0)"
				, species.Data()
				, kArgTPC.Data()
				, inel
				, lowPtBin
				, jointPtBin));
		
	cout << "Config_" << species << "_TOF_LHC10x.C" << endl << endl;
	gROOT->ProcessLine(Form(".x Config_%s_TOF_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 1,1,1,0,0)"
				, species.Data()
				, kArgTOF.Data()
				, inel
				, jointPtBin
				, hiPtBin));
		
	TString outputPtTPC = outputDir + "/" + MakeOutputName(species, kOutputTagTPC) + "-Pt.root";
	TString outputPtTOF = outputDir + "/" + MakeOutputName(species, kOutputTagTOF) + "-Pt.root";
	
	// combine TPC and TOF pt
	
	TString outputPt      = outputDir + "/" + MakeOutputName(species, outputTag) + "-Pt.root";
	TString outputRatio   = outputDir + "/" + MakeOutputName(species, outputTag) + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + MakeOutputName(species, outputTag) + "-Spectra.root";
	
	TFileMerger m;
	
	m.AddFile(outputPtTPC.Data(),0);
	m.AddFile(outputPtTOF.Data(),0);
	
	m.OutputFile(outputPt.Data());
	
	m.Merge();
	
	// make ratio and spectra
	
	AliLnDriver driver;
	
	driver.SetSpecies(species);
	
	Double_t xsec[3];
	GetInelXSection(xsec, period);
	
	driver.SetInelXSection(xsec);
	driver.SetExtrapolateToINEL(inel);

	if(species == "Proton") driver.SetSysErr(kProtonSysErr[0],kProtonSysErr[1]);
	if(species == "Deuteron") driver.SetSysErr(kDeuteronSysErr[0],kDeuteronSysErr[1]);
	
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	
	driver.SetOutputTag(outputTag);
	
	driver.SetMakeCorrections(0);
	driver.SetMakePt(0);
	driver.SetMakeRatio(1);
	driver.SetMakeSpectra(1);
	
	driver.Run();
	
	if(!drawOutput) return 0;
	
	TStyle* st = GetDrawingStyle();
	st->cd();
	gROOT->ForceStyle();
	
	DrawOutputRatio(outputRatio, outputTag, species);
	DrawOutputSpectra(outputSpectra, outputTag, species);
	
	return 0;
}
