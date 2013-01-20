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

#include <Riostream.h>
#include <TSystem.h>
#include <TString.h>
#include <TFileMerger.h>

#include "AliLnDriver.h"
#include "Config.h"

Int_t Config_Deuteron_TPCTOF_LHC10x(const TString& inputDir   = "~/alice/input",
                                    const TString& outputDir  = "~/alice/output",
                                    const TString& period     = "lhc10d",
                                    const TString& outputTag  = "lhc10d",
                                    const TString& multTag    = "",
                                    const TString& multCorTag = "",
                                    Bool_t normToInel         = 1,  // for mult
                                    Bool_t drawOutput         = 1)  // for batch
{
//
// lhc10b, lhc10c, lhc10d, lhc10e config for deuterons and antideuterons
// (combine TPC and TOF)
//
	const TString  kSpecies         = "Deuteron";
	
	const TString  kTrkSelTPC       = "its_tpc_dca_spd-tpc3";
	const TString  kTrkSelTOF       = "its_tpc_tof_dca-tpc3";
	
	const TString  kOutputTagTPC    = outputTag + "-tpc";
	const TString  kOutputTagTOF    = outputTag + "-tof";
	
	const Int_t    kLowPtBin        = 4; // bin width = 0.2 GeV/c
	const Int_t    kJointBin        = 6;
	const Int_t    kHiPtBin         = 18;
	
	const TString  kTrigName        = "mbor";
	const Bool_t   kVtxCorr         = 0;
	const Double_t kVtxCorrVal      = GetVertexCorrection(period);
	
	const Bool_t   kUnfolding       = 1;
	const Int_t    kIter            = 5;
	const Bool_t   kFakeTracks      = 0;
	const Bool_t   kSecondaries     = 1;
	const Int_t    kSecProd         = 0; // 0 tff, 1 roofit, 2 mc
	const Int_t    kMatDCAxyMod     = 1; // 0 geant, 1 flat
	const Bool_t   kAntiNucTemplate = 0;
	const Int_t    kNbin            = 5;
	const Double_t kDCAxy[2]        = {-0.2,0.2};
	const Bool_t   kEfficiency      = 1;
	const Bool_t   kFitFrac         = 0;
	
	// for TOF
	const Int_t    kM2Bin[2]        = {8,18};
	const Double_t kM2Bkg[2]        = {2.2,5.};
	const Double_t kM2tpc[2]        = {2.,6.5};
	
	const Double_t kSysErr[2]       = {0.10, 0.11} ;
	
	Double_t xsec[3];
	GetInelXSection(xsec, period);
	
	Double_t trigEff[3];
	GetTriggerEfficiency(trigEff, kTrigName, period);
	
	// common options
	
	AliLnDriver driver;
	
	driver.SetSpecies(kSpecies);
	
	driver.SetTriggerEfficiency(trigEff);
	driver.SetInelXSection(xsec);
	driver.SetNormalizeToINEL(normToInel);
	driver.SetVertexCorrection(kVtxCorr, kVtxCorrVal);
	
	driver.SetPidM2(0);
	driver.SetUnfolding(kUnfolding, kIter);
	driver.SetFakeTracks(kFakeTracks);
	driver.SetSecondaries(kSecondaries);
	driver.SetSecProd(kSecProd);
	driver.SetMatDCAxyModel(kMatDCAxyMod);
	driver.SetAntiNucleusAsTemplate(kAntiNucTemplate);
	driver.SetNBin(kNbin);
	driver.SetDCAxyInterval(kDCAxy[0], kDCAxy[1]);
	driver.SetEfficiency(kEfficiency,0);
	driver.SetFitFractionCorr(kFitFrac);
	
	driver.SetSysErr(kSysErr[0],kSysErr[1]);
	
	// get the pt with TPC up to kJointBin
	
	TString inputDataTPC     = inputDir + "/" + period + "/"
	                         + MakeInputName(kSpecies, period, kTrkSelTPC+multTag) + ".root";
	
	TString inputSimuTPC     = inputDir + "/" + period + "/"
	                         + MakeSimuName(kSpecies, period, kTrkSelTPC+multCorTag) + ".root";
	
	TString inputSimuFixTPC  = inputDir + "/" + period + "/"
	                         + MakeSimuFixName(kSpecies, period, kTrkSelTPC+multCorTag, 0) + ".root";
	
	TString inputCorrTPC     = inputDir + "/" + period + "/"
	                         + MakeInputName(kSpecies, period, kTrkSelTPC+multTag) + "-corr.root";
	
	TString outputPtTPC      = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTPC) + "-Pt.root";
	TString outputRatioTPC   = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTPC) + "-Ratio.root";
	TString outputSpectraTPC = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTPC) + "-Spectra.root";
	
	
	driver.SetInputFilenames(inputDataTPC, inputSimuTPC, inputSimuFixTPC, inputCorrTPC);
	driver.SetOutputFilenames(outputPtTPC, outputRatioTPC, outputSpectraTPC);
	
	driver.SetOutputTag(kOutputTagTPC);
	
	driver.SetMakeStats(0);
	driver.SetMakeCorrections(1);
	driver.SetMakePt(1);
	driver.SetMakeRatio(0);
	driver.SetMakeSpectra(0);
	
	driver.SetPtBinInterval(kLowPtBin, kJointBin);
	
	TString outputCorrDebugTPC = driver.GetPtCorrDebugFilename();
	TString outputPtDebugTPC   = driver.GetPtDebugFilename();
	
	driver.Run();
	
	// get the pt with TOF from kJointBin to kHiPtBin
	
	TString inputDataTOF     = inputDir + "/" + period + "/"
	                         + MakeInputName(kSpecies, period, kTrkSelTOF+multTag) + ".root";
	
	TString inputSimuTOF     = inputDir + "/" + period + "/"
	                         + MakeSimuName(kSpecies, period, kTrkSelTOF+multCorTag) + ".root";
	
	TString inputSimuFixTOF  = inputDir + "/" + period + "/"
	                         + MakeSimuFixName(kSpecies, period, kTrkSelTOF+multCorTag) + ".root";
	
	TString inputCorrTOF     = inputDir + "/" + period + "/"
	                         + MakeInputName(kSpecies, period, kTrkSelTOF+multTag) + "-corr.root";
	
	TString outputPtTOF      = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTOF) + "-Pt.root";
	TString outputRatioTOF   = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTOF) + "-Ratio.root";
	TString outputSpectraTOF = outputDir + "/" + MakeOutputName(kSpecies, kOutputTagTOF) + "-Spectra.root";
	
	
	driver.SetInputFilenames(inputDataTOF, inputSimuTOF, inputSimuFixTOF, inputCorrTOF);
	driver.SetOutputFilenames(outputPtTOF, outputRatioTOF, outputSpectraTOF);
	
	driver.SetOutputTag(kOutputTagTOF);
	
	driver.SetMakeStats(1);
	driver.SetMakeCorrections(1);
	driver.SetMakePt(1);
	driver.SetMakeRatio(0);
	driver.SetMakeSpectra(0);
	
	driver.SetPtBinInterval(kJointBin, kHiPtBin);
	driver.SetUnfolding(0);
	driver.SetEfficiency(kEfficiency,0);
	driver.SetPidM2(1);
	driver.SetM2BinInterval(kM2Bin[0], kM2Bin[1]);
	driver.SetM2BkgInterval(kM2Bkg[0], kM2Bkg[1]);
	driver.SetM2TPCInterval(kM2tpc[0], kM2tpc[1]);
	
	TString outputCorrDebugTOF = driver.GetPtCorrDebugFilename();
	TString outputPtDebugTOF   = driver.GetPtDebugFilename();
	
	driver.Run();
	
	// combine TPC and TOF pt
	
	TString outputPt      = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Pt.root";
	TString outputRatio   = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Spectra.root";
	
	TFileMerger m;
	
	m.AddFile(outputPtTPC.Data(),0);
	m.AddFile(outputPtTOF.Data(),0);
	
	m.OutputFile(outputPt.Data());
	
	m.Merge();
	
	// make ratio and spectra
	
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	
	driver.SetOutputTag(outputTag);
	
	driver.SetMakeCorrections(0);
	driver.SetMakePt(0);
	driver.SetMakeRatio(1);
	driver.SetMakeSpectra(1);
	
	driver.Run();
	
	// delete tmp files
	
	gSystem->Exec(Form("rm -f %s %s %s %s", inputCorrTPC.Data(), outputPtTPC.Data(), outputCorrDebugTPC.Data(), outputPtDebugTPC.Data()));
	gSystem->Exec(Form("rm -f %s %s %s %s", inputCorrTOF.Data(), outputPtTOF.Data(), outputCorrDebugTOF.Data(), outputPtDebugTOF.Data()));
	
	// draw output
	
	if(!drawOutput) return 0;
	
	TStyle* st = GetDrawingStyle();
	st->cd();
	gROOT->ForceStyle();
	
	DrawOutputRatio(outputRatio, outputTag, kSpecies);
	DrawOutputSpectra(outputSpectra, outputTag, kSpecies);
	
	return 0;
}
