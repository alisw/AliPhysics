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

// LHC10x config for protons and antiprotons
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TString.h>
#include "AliLnDriver.h"
#include "Config.h"

Int_t Config_Proton_TOF_LHC10x(const TString& inputDir   = "~/alice/input",
                               const TString& outputDir  = "~/alice/output",
                               const TString& period     = "lhc10d",
                               const TString& outputTag  = "lhc10d",
                               const TString& multTag    = "",
                               const TString& multCorTag = "",
                               Bool_t inel               = 1,  // for mult
                               Bool_t drawOutput         = 1,  // for batch
                               Int_t  lowPtBin           = 11, // 1.0 GeV/c
                               Int_t  hiPtBin            = 36, // 3.5 GeV/c
                               Bool_t makeStats          = 1,
                               Bool_t makeCor            = 1,
                               Bool_t makePt             = 1,
                               Bool_t makeRatio          = 1,
                               Bool_t makeSpectra        = 1 )
{
//
// lhc10b, lhc10c, lhc10d, lhc10e config for protons and antiprotons
// (TOF)
//
	const TString  kSpecies     = "Proton";
	const TString  kTrkSel      = "its_tpc_tof_dca_spd-bayes";
	const TString  kTrigName    = "mbor";
	const Bool_t   kMCtoINEL    = 1;
	const Bool_t   kUnfolding   = 0;
	const Bool_t   kSecondaries = 1;
	const Int_t    kSecProd     = 0; // 0 tff, 1 roofit, 2 mc
	const Int_t    kMatDCAxyMod = 1; // 0 geant, 1 flat
	const Int_t    kNbin        = 15;
	const Double_t kDCAxy[2]    = {-1.,1.};
	const Bool_t   kEfficiency  = 1;
	const Double_t kMatScaling  = 1.9;
	const Double_t kFdwnScaling = 1.9;
	const Bool_t   kFitFrac     = 0;
	const Bool_t   kSameFdwn    = 1;
	const Double_t kSysErr[2]   = {0.08,0.08} ;
	
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
	driver.SetPidM2(0);
	driver.SetUnfolding(kUnfolding);
	driver.SetSecondaries(kSecondaries);
	driver.SetSecProd(kSecProd);
	driver.SetMatDCAxyModel(kMatDCAxyMod);
	driver.SetNBin(kNbin);
	driver.SetDCAxyInterval(kDCAxy[0], kDCAxy[1]);
	driver.SetEfficiency(kEfficiency);
	driver.SetScalingFactors(kMatScaling, kFdwnScaling);
	driver.SetFitFractionCorr(kFitFrac);
	driver.SetSameFeedDownCorr(kSameFdwn);
	
	driver.SetSysErr(kSysErr[0],kSysErr[1]);
	
	//driver.PrintFilenames();
	
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
	
	if(kSecProd != 2) DrawCorrDebug(driver.GetPtCorrDebugFilename(), driver.GetOutputCorrTag(), kSpecies, lowPtBin, hiPtBin, kDCAxy[0], kDCAxy[1]);
	
	DrawPtDebug(driver.GetPtDebugFilename(), outputTag, kSpecies, 0);
	DrawOutputRatio(outputRatio, outputTag, kSpecies);
	DrawOutputSpectra(outputSpectra, outputTag, kSpecies);
	
	return 0;
}
