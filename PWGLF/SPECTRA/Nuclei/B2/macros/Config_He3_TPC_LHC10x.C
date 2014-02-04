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

// LHC10x config for He3 and AntiHe3
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TString.h>
#include "AliLnDriver.h"
#endif

#include "Config.h"

Int_t Config_He3_TPC_LHC10x(  const TString& inputDir   = "~/alice/input"
                            , const TString& outputDir  = "~/alice/output"
                            , const TString& period     = "lhc10bcde"
                            , const TString& outputTag  = "pp-mband-7TeV"
                            , const TString& trkselTag  = "-tpc3-nsd-moc-vbin"
                            , const TString& multTag    = ""
                            , const TString& multCorTag = ""
                            , Double_t ymax             = 0.5
                            , Bool_t   inel             = kFALSE
                            , Bool_t   drawOutput       = kTRUE
                            , Double_t ptmin            = 0.4  // GeV/c
                            , Double_t ptmax            = 10.  // GeV/c
                            , Bool_t   makeStats        = kTRUE
                            , Bool_t   makeCor          = kTRUE
                            , Bool_t   makePt           = kTRUE
                            , Bool_t   makeRatio        = kTRUE
                            , Bool_t   makeSpectra      = kTRUE )
{
//
// lhc10b, lhc10c, lhc10d, lhc10e config for He3 and AntiHe3
// (TPC)
//
	const TString  kSpecies         = "He3";
	const TString  kTrkSel          = "its_tpc_dca";
	const TString  kTrigName        = "mband";
	const Bool_t   kMCtoINEL        = kTRUE;
	const Double_t kPidEff          = 1.;
	const Bool_t   kSecondaries     = kFALSE;
	const Int_t    kSecProc         = 1; // 0 tff, 1 mc
	const Int_t    kMatDCAxyMod     = 1; // 0 geant, 1 flat
	const Bool_t   kAntiNucTemplate = kFALSE;
	const Int_t    kNbin            = 5;
	const Double_t kDCAxy[2]        = {-0.2, 0.2};
	const Bool_t   kEfficiency      = kTRUE;
	const Int_t    kDebugLevel      = 1;
	
	Double_t trigEff[3];
	GetTriggerEfficiency(trigEff, kTrigName, period);
	
	// input and output filenames
	
	TString inputData     = inputDir + "/" + period + "/" + MakeInputName(kSpecies, period, kTrkSel+trkselTag+multTag) + ".root";
	TString inputSimu     = inputDir + "/" + period + "/" + MakeSimuName(kSpecies, period, kTrkSel+trkselTag+multCorTag) + ".root";
	TString inputSimuFix  = inputDir + "/" + period + "/" + MakeSimuFixName(kSpecies, period, kTrkSel+trkselTag+multCorTag) + ".root";
	TString inputCorr     = inputDir + "/" + period + "/" + MakeInputName(kSpecies, period, kTrkSel+trkselTag+multTag) + "-corr.root";
	
	TString outputPt      = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Pt.root";
	TString outputRatio   = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Ratio.root";
	TString outputSpectra = outputDir + "/" + MakeOutputName(kSpecies, outputTag) + "-Spectra.root";
	
	// configure the driver and run
	
	AliLnDriver driver;
	
	driver.SetSpecies(kSpecies);
	
	driver.SetInputFilenames(inputData, inputSimu, inputSimuFix, inputCorr);
	driver.SetOutputFilenames(outputPt, outputRatio, outputSpectra);
	
	driver.SetRapidityInterval(-ymax,ymax);
	
	driver.SetOutputTag(outputTag);
	driver.SetOutputCorTag(outputTag);
	
	driver.SetTriggerEfficiency(trigEff);
	driver.SetExtrapolateToINEL(inel);
	driver.SetMCtoINEL(kMCtoINEL);
	driver.SetPtInterval(ptmin, ptmax);
	driver.SetPid(0);
	driver.SetPidEfficiency(kPidEff);
	driver.SetSecondaries(kSecondaries);
	driver.SetSecProcedure(kSecProc);
	driver.SetMatDCAxyModel(kMatDCAxyMod);
	driver.SetAntiNucleusAsTemplate(kAntiNucTemplate);
	driver.SetNBin(kNbin);
	driver.SetDCAxyInterval(kDCAxy[0], kDCAxy[1]);
	driver.SetEfficiency(kEfficiency,0);
	
	driver.SetMakeStats(makeStats);
	driver.SetMakeCorrections(makeCor);
	driver.SetMakePt(makePt);
	driver.SetMakeRatio(makeRatio);
	driver.SetMakeSpectra(makeSpectra);
	
	driver.SetDebugLevel(kDebugLevel);
	
	driver.Run();
	
	// draw output
	
	if(!drawOutput) return 0;
	
	DrawOutputCorr(kSpecies, inputCorr, driver.GetOutputCorrTag());
	
	DrawPtDebug(driver.GetPtDebugFilename(), outputTag, kSpecies, 0);
	DrawOutputRatio(outputRatio, outputTag, kSpecies);
	DrawOutputSpectra(outputSpectra, outputTag, kSpecies);
	
	return 0;
}
