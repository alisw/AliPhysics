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
                           const TString& species    = "Proton")
{
//
// combine TPC and TOF for protons and deuterons
//
	const Int_t kProtonLowPtBin    = 5 ;
	const Int_t kProtonJointBin    = 11;
	const Int_t kProtonHiPtBin     = 36;
	
	const Int_t kDeuteronLowPtBin  = 4;
	const Int_t kDeuteronJointBin  = 6;
	const Int_t kDeuteronHiPtBin   = 13;
	
	const Double_t kProtonSysErr[2]   = {0.08, 0.08} ;
	const Double_t kDeuteronSysErr[2] = {0.10, 0.11} ;
	
	// -------------
	
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
	
	using namespace std;
	
	if(species == "Proton")
	{
		cout << "Config_Proton_TPC_LHC10x.C" << endl << endl;
		gROOT->ProcessLine(Form(".x Config_Proton_TPC_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 0,1,1,0,0)"
					, kArgTPC.Data()
					, inel
					, kProtonLowPtBin
					, kProtonJointBin));
		
		cout << "Config_Proton_TOF_LHC10x.C" << endl << endl;
		gROOT->ProcessLine(Form(".x Config_Proton_TOF_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 1,1,1,0,0)"
					, kArgTOF.Data()
					, inel
					, kProtonJointBin
					, kProtonHiPtBin));
	}
	else if (species == "Deuteron")
	{
		cout << "Config_Deuteron_TPC_LHC10x.C" << endl << endl;
		gROOT->ProcessLine(Form(".x Config_Deuteron_TPC_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 0,1,1,0,0)"
					, kArgTPC.Data()
					, inel
					, kDeuteronLowPtBin
					, kDeuteronJointBin));
		
		cout << "Config_Deuteron_TOF_LHC10x.C" << endl << endl;
		gROOT->ProcessLine(Form(".x Config_Deuteron_TOF_LHC10x.C+g(\"%s\", %d, 0, %d, %d, 1,1,1,0,0)"
					, kArgTOF.Data()
					, inel
					, kDeuteronJointBin
					, kDeuteronHiPtBin));
	}
	else
	{
		cerr << "Particle species " << species << " not implemented." << endl;
		exit(1);
	}
	
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
