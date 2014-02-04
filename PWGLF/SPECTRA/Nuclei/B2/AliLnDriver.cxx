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

// driver for computing the pt and spectra
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TObject.h>
#include <TSystem.h>
#include <TFileMerger.h>
#include <TString.h>
#include <TFile.h>
#include <TError.h>

#include "AliLnCorr.h"
#include "AliLnPt.h"
#include "AliLnSecondaries.h"
#include "AliLnEfficiency.h"
#include "AliLnRatio.h"
#include "AliLnSpectra.h"
#include "AliLnDriver.h"
#include "B2.h"


ClassImp(AliLnDriver)

AliLnDriver::AliLnDriver()
: TObject()
, fSpecies("Deuteron")
, fOutputTag("test")
, fOutputCorTag("")
, fIsOnlyGen(0)
, fINEL(1)
, fMakeCorr(1)
, fMakePt(1)
, fMakeRatio(1)
, fMakeSpectra(1)
, fMakeStats(1)
, fPtMin(3)
, fPtMax(15)
, fPid(0.4)
, fPidPt(1.)
, fSecondaries(1)
, fSecProc(AliLnSecondaries::kTFractionFitter)
, fMatDCAxyMod(AliLnSecondaries::kGeantDCAxy)
, fANucTemplate(0)
, fNbin(1)
, fYMin(-0.5)
, fYMax(0.5)
, fMinDCAxy(-1.5)
, fMaxDCAxy(1.5)
, fBkgMin(2.2)
, fBkgMax(5.)
, fIntMin(2.)
, fIntMax(6.5)
, fEfficiency(1)
, fG3Fluka(0)
, fScMat(1)
, fScFd(1)
, fInputData()
, fInputSimu()
, fInputSimuFix()
, fOutputPtCorr()
, fOutputPt()
, fOutputRatio()
, fOutputSpectra()
, fOutputPtCorrDebug()
, fOutputPtDebug()
, fFitFrac(0)
, fFdwnCorr(1)
, fSameFdwn(0)
, fMCtoINEL(0)
, fAddFakeTracks(1)
, fPidProc(AliLnPt::kMassSquared)
, fPidEff(1)
, fDebugLevel(0)
{
//
// constructor
//
	for(Int_t i=0; i<3; ++i)
	{
		fTrigEff[i] = 1;
	}
	
	gErrorIgnoreLevel = kWarning;
	// gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
}

AliLnDriver::~AliLnDriver()
{
//
// destructor
//
}

void AliLnDriver::SetInputFilenames(const TString& data, const TString& simu, const TString& simuFix, const TString& ptcorr)
{
//
// input filenames to get the pt
//
	fInputData = data;
	fInputSimu = simu;
	fInputSimuFix = simuFix;
	fOutputPtCorr = ptcorr;
	
	fOutputPtCorrDebug = ptcorr;
	fOutputPtCorrDebug.Replace(ptcorr.Length()-5,5,"-debug.root");
}

void AliLnDriver::SetOutputFilenames(const TString& pt, const TString& ratio, const TString& spectra)
{
//
// output filenames
//
	fOutputPt = pt;
	fOutputRatio=ratio;
	fOutputSpectra=spectra;
	
	fOutputPtDebug = pt;
	fOutputPtDebug.Replace(pt.Length()-5,5,"-debug.root");
}

void AliLnDriver::SetTriggerEfficiency(Double_t eff[3])
{
//
// set trigger efficiency
//
	for(Int_t i=0; i<3; ++i) fTrigEff[i] = eff[i];
}

Int_t AliLnDriver::Run() const
{
//
// run script
//
	using namespace std;
	
	if(!fIsOnlyGen && fMakeCorr)
	{
		if(fDebugLevel > 0 )
		{
			cout << endl;
			cout << "********* Rebuild corrections ***********" << endl;
			cout << "Species               : " << fSpecies << endl;
			cout << "Simulation file       : " << fInputSimu << endl;
			cout << "Simulation file (fix) : " << fInputSimuFix << endl;
			cout << "Output file           : " << fOutputPtCorr << endl;
			cout << "Output file (debug)   : " << fOutputPtCorrDebug << endl;
			cout << endl;
		}
		
		this->MakePtCorr();
	}
	
	if(fMakePt)
	{
		if(fDebugLevel > 0 )
		{
			cout << endl;
			cout << "********* Make pt ***********" << endl;
			cout << "Species             : " << fSpecies << endl;
			cout << "pt interval (GeV/c) : (" << fPtMin << ", " << fPtMax << ")" << endl;
			cout << "pid pt (GeV/c)      : " << fPidPt <<  endl;
			cout << "Input file          : " << fInputData << endl;
			cout << "Correction file     : " << fOutputPtCorr << endl;
			cout << "Output file         : " << fOutputPt << endl;
			cout << "Output file (debug) : " << fOutputPtDebug << endl;
			cout << endl;
		}
		
		this->MakePt();
	}
	
	if(fMakeRatio)
	{
		if(fDebugLevel > 0 )
		{
			cout << endl;
			cout << "********* Make ratio ***********" << endl;
			cout << "Species             : " << fSpecies << endl;
			cout << "Input file          : " << fOutputPt << endl;
			cout << "Output file         : " << fOutputRatio << endl;
			cout << endl;
		}
		
		this->MakeRatio();
	}
	
	if(fMakeSpectra)
	{
		if(fDebugLevel > 0 )
		{
			cout << endl;
			cout << "********* Make spectra ***********" << endl;
			cout << "Species             : " << fSpecies << endl;
			cout << "INEL extrapolation  : " << fINEL << endl;
			cout << "Rapidity interval   : (" << fYMin << ", " << fYMax << ")" << endl;
			cout << "Input file          : " << fOutputPt << endl;
			cout << "Output file         : " << fOutputSpectra << endl;
			cout << endl;
		}
		
		this->MakeSpectra();
	}
	
	return 0;
}

void AliLnDriver::MakePtCorr() const
{
//
// make pt correction file
//
	const Int_t kNpart = 2;
	const TString kParticle[kNpart] = {fSpecies, Form("Anti%s",fSpecies.Data())};
	
	TFileMerger m1, m2;
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		TString outputfile1 = kParticle[i] + "-corrections.root";
		TString outputfile2 = Form("debug-%s", outputfile1.Data());
		
		AliLnCorr lncorr(kParticle[i], fInputData, fInputSimu, fInputSimuFix, outputfile1, fOutputCorTag);
		
		lncorr.GetLnSecondaries()->SetCorBins(fPtMin, fPtMax);
		lncorr.GetLnSecondaries()->SetProcedure(fSecProc);
		lncorr.GetLnSecondaries()->SetMatDCAxyModel(fMatDCAxyMod);
		lncorr.GetLnSecondaries()->SetAntiNucleusAsTemplate(fANucTemplate);
		lncorr.GetLnSecondaries()->SetDCAxyInterval(fMinDCAxy, fMaxDCAxy);
		lncorr.GetLnSecondaries()->SetNBin(fNbin);
		lncorr.GetLnSecondaries()->SetScalingFactors(fScMat, fScFd);
		lncorr.GetLnSecondaries()->SetAddFakeTracks(fAddFakeTracks);
		
		lncorr.GetLnEfficiency()->SetG3Fluka(fG3Fluka);
		lncorr.GetLnEfficiency()->SetAddFakeTracks(fAddFakeTracks);
		
		lncorr.Exec();
		
		m1.AddFile(outputfile1.Data(),0);
		m2.AddFile(outputfile2.Data(),0);
	}
	
	// merge and remove tmp files
	
	m1.OutputFile(fOutputPtCorr.Data());
	m1.Merge();
	
	gSystem->Exec(Form("rm -f %s-corrections.root Anti%s-corrections.root", fSpecies.Data(), fSpecies.Data()));
	
	m2.OutputFile(fOutputPtCorrDebug.Data());
	m2.Merge();
	
	gSystem->Exec(Form("rm -f debug-%s-corrections.root debug-Anti%s-corrections.root", fSpecies.Data(), fSpecies.Data()));
}

void AliLnDriver::MakePt() const
{
//
// make pt from data and correction file
//
	const Int_t kNpart = 2;
	
	const TString kParticle[kNpart] = { fSpecies, Form("Anti%s",fSpecies.Data())};
	
	TFileMerger m, m2;
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		TString ptfile   = kParticle[i] + "-Pt.root";
		
		AliLnPt lnpt(kParticle[i], fTrigEff[0], fInputData, ptfile, fOutputTag, fOutputPtCorr, fOutputCorTag);
		
		lnpt.SetOnlyGeneration(fIsOnlyGen);
		
		lnpt.SetPtInterval(fPtMin, fPtMax);
		lnpt.SetPid(fPid);
		lnpt.SetPidProcedure(fPidProc);
		lnpt.SetPidPt(fPidPt);
		lnpt.SetBkgInterval(fBkgMin, fBkgMax);
		lnpt.SetPidInterval(fIntMin, fIntMax);
		lnpt.SetSecondaries(fSecondaries);
		lnpt.SetEfficiency(fEfficiency);
		lnpt.SetMakeStats(fMakeStats);
		lnpt.SetMCtoINEL(fMCtoINEL);
		lnpt.SetFitFractionCorr(fFitFrac);
		lnpt.SetFeedDownCorr(fFdwnCorr);
		lnpt.SetSameFeedDownCorr(fSameFdwn);
		lnpt.SetPidEfficiency(fPidEff);
		lnpt.SetDebugLevel(fDebugLevel);
		
		lnpt.Exec();
		
		m.AddFile(ptfile.Data(),0);
		if(!fIsOnlyGen) m2.AddFile(Form("debug-%s",ptfile.Data()),0);
	}
	
	// merge and remove tmp files
	
	m.OutputFile(fOutputPt.Data());
	m.Merge();
	
	gSystem->Exec(Form("rm -f %s-Pt.root Anti%s-Pt.root", fSpecies.Data(), fSpecies.Data()));
	
	if(!fIsOnlyGen)
	{
		m2.OutputFile(fOutputPtDebug.Data());
		m2.Merge();
		
		gSystem->Exec(Form("rm -f debug-%s-Pt.root debug-Anti%s-Pt.root", fSpecies.Data(),fSpecies.Data()));
	}
}

void AliLnDriver::MakeRatio() const
{
//
// make antiparticle/particle ratio
//
	AliLnRatio lnr(fSpecies, fOutputPt, fOutputTag, fOutputRatio, fOutputTag);
	
	lnr.Exec();
}

void AliLnDriver::MakeSpectra() const
{
//
// make differential yields
//
	TFileMerger m;
	
	const TString kParticle[2] = { fSpecies, Form("Anti%s",fSpecies.Data())};
	
	for(Int_t i=0; i<2; ++i)
	{
		TString spectrafile = kParticle[i] + "-Spectra.root";
		
		AliLnSpectra lnspectra(kParticle[i], fOutputPt, fOutputTag, spectrafile, fOutputTag);
		
		lnspectra.SetRapidityInterval(fYMin, fYMax);
		lnspectra.SetExtrapolateToINEL(fINEL);
		
		lnspectra.Exec();
		
		m.AddFile(spectrafile.Data(),0);
	}
	
	// merge and remove tmp files
	
	m.OutputFile(fOutputSpectra.Data());
	m.Merge();
	
	gSystem->Exec(Form("rm -f %s-Spectra.root Anti%s-Spectra.root", fSpecies.Data(), fSpecies.Data()));
}
