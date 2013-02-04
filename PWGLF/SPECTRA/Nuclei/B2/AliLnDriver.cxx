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
#include "AliLnB2.h"
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
, fLowPtBin(3)
, fHiPtBin(15)
, fLowM2Bin(9)
, fHighM2Bin(17)
, fUnfolding(0)
, fNIter(4)
, fSecondaries(1)
, fSecProd(AliLnSecondaries::kTFractionFitter)
, fMatDCAxyMod(AliLnSecondaries::kGeantDCAxy)
, fANucTemplate(0)
, fNbin(1)
, fYMin(-0.5)
, fYMax(0.5)
, fMinDCAxy(-1.5)
, fMaxDCAxy(1.5)
, fMinM2Bkg(2.2)
, fMaxM2Bkg(5.)
, fMinM2tpc(2.)
, fMaxM2tpc(6.5)
, fEfficiency(1)
, fG3Fluka(0)
, fScMat(1)
, fScFd(1)
, fSysPos(1)
, fSysNeg(1)
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
{
//
// constructor
//
	for(Int_t i=0; i<3; ++i)
	{
		fTrigEff[i] = 1;
		fXsec[i] = 1;
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

void AliLnDriver::PrintFilenames() const
{
//
// print filenames to stdout
//
	using namespace std;
	
	cout << endl;
	cout << "input data        : " << fInputData << endl;
	cout << "input simulation  : " << fInputSimu << endl;
	cout << "input simu. fix   : " << fInputSimuFix << endl;
	
	cout << "corrections       : " << fOutputPtCorr << endl;
	cout << "corrections debug : " << fOutputPtCorrDebug << endl;
	
	cout << "output pt         : " << fOutputPt << endl;
	cout << "output pt debug   : " << fOutputPtDebug << endl;
	cout << "output ratio      : " << fOutputRatio << endl;
	cout << "output spectra    : " << fOutputSpectra << endl;
	cout << endl;
}

void AliLnDriver::SetTriggerEfficiency(Double_t eff[3])
{
//
// set trigger efficiency
//
	for(Int_t i=0; i<3; ++i) fTrigEff[i] = eff[i];
}

void AliLnDriver::SetInelXSection(Double_t xsec[3])
{
//
// set total cross section, stat. and syst. errors
//
	for(Int_t i=0; i<3; ++i) fXsec[i] = xsec[i];
}

Int_t AliLnDriver::Run() const
{
//
// run script
//
	if(!fIsOnlyGen && fMakeCorr) this->MakePtCorr();
	
	if(fMakePt) this->MakePt();
	
	if(fMakeRatio) this->MakeRatio();
	
	if(fMakeSpectra) this->MakeSpectra();
	
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
		
		lncorr.GetLnSecondaries()->SetCorBins(fLowPtBin, fHiPtBin);
		lncorr.GetLnSecondaries()->SetProcedure(fSecProd);
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
		lnpt.SetRapidityInterval(fYMin, fYMax);
		
		lnpt.SetPtBinInterval(fLowPtBin, fHiPtBin);
		lnpt.SetM2BinInterval(fLowM2Bin, fHighM2Bin);
		lnpt.SetM2BkgInterval(fMinM2Bkg, fMaxM2Bkg);
		lnpt.SetM2TPCInterval(fMinM2tpc, fMaxM2tpc);
		lnpt.SetPidM2(fPidM2);
		lnpt.SetUnfolding(fUnfolding, fNIter);
		lnpt.SetSecondaries(fSecondaries);
		lnpt.SetEfficiency(fEfficiency);
		lnpt.SetMakeStats(fMakeStats);
		lnpt.SetMCtoINEL(fMCtoINEL);
		lnpt.SetFitFractionCorr(fFitFrac);
		lnpt.SetFeedDownCorr(fFdwnCorr);
		lnpt.SetSameFeedDownCorr(fSameFdwn);
		
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
	
	const TString kParticle[]  = { fSpecies, Form("Anti%s",fSpecies.Data())};
	const Double_t kSysErr[] = { fSysPos, fSysNeg };
	
	for(Int_t i=0; i<2; ++i)
	{
		TString spectrafile = kParticle[i] + "-Spectra.root";
		
		AliLnSpectra lnspectra(kParticle[i], fOutputPt, fOutputTag, spectrafile, fOutputTag, fXsec);
		
		lnspectra.SetRapidityInterval(fYMin, fYMax);
		lnspectra.SetExtrapolateToINEL(fINEL);
		lnspectra.SetOnlyGeneration(fIsOnlyGen);
		lnspectra.SetScalingFactor(kSysErr[i]);
		
		lnspectra.Exec();
		
		m.AddFile(spectrafile.Data(),0);
	}
	
	// merge and remove tmp files
	
	m.OutputFile(fOutputSpectra.Data());
	m.Merge();
	
	gSystem->Exec(Form("rm -f %s-Spectra.root Anti%s-Spectra.root", fSpecies.Data(), fSpecies.Data()));
}
