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

// LHC10x config for deuterons and antideuterons
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TROOT.h>
#include <TString.h>
#include "AliLnDriver.h"
#include "Config.h"

Int_t Config_Deuteron_TOF_LHC10x(const TString& inputDir   = "~/alice/input",
                                 const TString& outputDir  = "~/alice/output",
                                 const TString& period     = "lhc10d",
                                 const TString& outputTag  = "lhc10d",
                                 const TString& multTag    = "",
                                 const TString& multCorTag = "",
                                 Bool_t inel               = 1,  // for mult
                                 Bool_t drawOutput         = 1,  // for batch
                                 Int_t  lowPtBin           = 5,  // 0.8 Gev/c
                                 Int_t  hiPtBin            = 17, // 3.2 GeV/c
                                 Bool_t makeStats          = 1,
                                 Bool_t makeCor            = 1,
                                 Bool_t makePt             = 1,
                                 Bool_t makeRatio          = 1,
                                 Bool_t makeSpectra        = 1 )
{
//
// lhc10b, lhc10c, lhc10d, lhc10e config for deuterons and antideuterons
//
	const TString  kSpecies         = "Deuteron";
	const TString  kTrkSel          = "its_tpc_tof_dca-tpc3";
	const TString  kTrigName        = "mbor";
	const Bool_t   kMCtoINEL        = 1;
	const Int_t    kM2Bin[2]        = {8,18};
	const Bool_t   kPidM2           = 1;
	const Bool_t   kUnfolding       = 0;
	const Int_t    kIter            = 5;
	const Bool_t   kSecondaries     = 1;
	const Int_t    kSecProd         = 0; // 0 tff, 1 roofit, 2 mc
	const Int_t    kMatDCAxyMod     = 1; // 0 geant, 1 flat
	const Bool_t   kAntiNucTemplate = 0;
	const Int_t    kNbin            = 5;
	const Double_t kDCAxy[2]        = {-0.2,0.2};
	const Double_t kM2Bkg[2]        = {2.2,5.};
	const Double_t kM2tpc[2]        = {2.,6.5};
	const Bool_t   kEfficiency      = 1;
	const Bool_t   kFitFrac         = 0;
	const Double_t kSysErr[2]       = {0.10, 0.11} ;
	
	Double_t xsec[3];
	GetInelXSection(xsec, period);
	
	Double_t trigEff[3];
	GetTriggerEfficiency(trigEff, kTrigName, period);
	
	// input and output filenames
	
	TString inputData     = inputDir + "/" + period + "/" + MakeInputName(kSpecies, period, kTrkSel+multTag) + ".root";
	TString inputSimu     = inputDir + "/" + period + "/" + MakeSimuName(kSpecies, period, kTrkSel+multCorTag) + ".root";
	TString inputSimuFix  = inputDir + "/" + period + "/" + MakeSimuFixName(kSpecies, period, kTrkSel+multCorTag) + ".root";
	TString inputCorr     = inputDir + "/" + period + "/" + MakeInputName(kSpecies, period, kTrkSel+multTag) + "-corr.root";
	
	TString outputPt      = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Pt.root";
	TString outputRatio   = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Spectra.root";
	
	// configure the driver and run
	
	AliLnDriver driver;
	
	driver.SetSpecies(kSpecies);
	
	driver.SetInputFilenames(inputData, inputSimu, inputSimuFix, inputCorr);
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	
	driver.SetOutputTag(outputTag);
	driver.SetTriggerEfficiency(trigEff);
	driver.SetInelXSection(xsec);
	driver.SetExtrapolateToINEL(inel);
	driver.SetMCtoINEL(kMCtoINEL);
	driver.SetPtBinInterval(lowPtBin, hiPtBin);
	driver.SetPidM2(kPidM2);
	driver.SetM2BinInterval(kM2Bin[0], kM2Bin[1]);
	driver.SetM2BkgInterval(kM2Bkg[0], kM2Bkg[1]);
	driver.SetM2TPCInterval(kM2tpc[0], kM2tpc[1]);
	driver.SetUnfolding(kUnfolding, kIter);
	driver.SetSecondaries(kSecondaries);
	driver.SetSecProd(kSecProd);
	driver.SetMatDCAxyModel(kMatDCAxyMod);
	driver.SetAntiNucleusAsTemplate(kAntiNucTemplate);
	driver.SetNBin(kNbin);
	driver.SetDCAxyInterval(kDCAxy[0], kDCAxy[1]);
	driver.SetEfficiency(kEfficiency,0);
	driver.SetFitFractionCorr(kFitFrac);
	driver.SetSysErr(kSysErr[0],kSysErr[1]);
	
	
	driver.SetMakeStats(makeStats);
	driver.SetMakeCorrections(makeCor);
	driver.SetMakePt(makePt);
	driver.SetMakeRatio(makeRatio);
	driver.SetMakeSpectra(makeSpectra);
	
	driver.Run();
	
	// draw output
	
	if(!drawOutput) return 0;
	
	TStyle* st = GetDrawingStyle();
	st->cd();
	gROOT->ForceStyle();
	
	DrawOutputCorr(kSpecies, inputCorr, driver.GetOutputCorrTag());
	
	if(kSecProd != 2) gROOT->ProcessLine(Form(".x DrawSec.C+g(\"%s\",\"\",\"Deuteron\", %d, %d, %f, %f)", driver.GetPtCorrDebugFilename().Data(), lowPtBin, hiPtBin, kDCAxy[0], kDCAxy[1]));
	
	DrawPtDebug(driver.GetPtDebugFilename(), outputTag, kSpecies, kPidM2, hiPtBin, kM2Bin[0], kM2Bin[1]);
	DrawOutputRatio(outputRatio, outputTag, kSpecies);
	DrawOutputSpectra(outputSpectra, outputTag, kSpecies);
	
	return 0;
}
