/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
#include <TVirtualFitter.h>
#include <TStopwatch.h>
#include "AliLog.h"
#include "AliDielectronBtoJPSItoEleCDFfitHandler.h"
#include "AliDielectronBtoJPSItoEleCDFfitFCN.h"

//-------------------------------------------------------------------------
//                      Class AliDielectronBtoJPSItoEleCDFfitHandler
//            Class to perform unbinned log-likelihood fit
//      
//                        Origin: C. Di Giglio
//     Contact: carmelo.digiglio@ba.infn.it , giuseppe.bruno@ba.infn.it
//-------------------------------------------------------------------------

void CDFFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag);

ClassImp(AliDielectronBtoJPSItoEleCDFfitHandler)

	//______________________________________________________________________________
void CDFFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	// This function is called by minuit
	// The corresponding member method is called
	// using SetObjectFit/GetObjectFit methods of TMinuit
	AliDielectronBtoJPSItoEleCDFfitHandler* dummy = (AliDielectronBtoJPSItoEleCDFfitHandler *)TVirtualFitter::GetFitter()->GetObjectFit();
	dummy->CdfFCN(npar, gin, f, par, iflag);
}


//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler::AliDielectronBtoJPSItoEleCDFfitHandler():
	fIsParamFixed(20),
	fPrintStatus(kFALSE),
	fUp(0),
	fX(0x0),
	fM(0x0),
	fLikely(0x0),
	fNcand(0)
{
	//
	// default constructor
	//
}
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler::AliDielectronBtoJPSItoEleCDFfitHandler(Double_t* decaytime, 
		Double_t* invariantmass, Int_t ncand) :
	fIsParamFixed(20),
	fPrintStatus(kFALSE),
	fUp(0),
	fX(decaytime),
	fM(invariantmass),
	fLikely(0x0),
	fNcand(ncand)
{
	//
	// constructor
	//
	AliInfo("\n+++\n+++ Minimization object AliDielectronBtoJPSItoEleCDFfitHandler created\n+++\n");
	fLikely = new AliDielectronBtoJPSItoEleCDFfitFCN();
	AliInfo("\n+++\n+++ CDF fit function object AliDielectronBtoJPSItoEleCDFfitFCN created\n+++\n");
	AliInfo("Parameter 0  ----> fWeightRes");
	AliInfo("Parameter 1  ----> fPos");
	AliInfo("Parameter 2  ----> fNeg");
	AliInfo("Parameter 3  ----> fSym");
	AliInfo("Parameter 4  ----> fOneOvLamPlus");
	AliInfo("Parameter 5  ----> fOneOvLamMinus");
	AliInfo("Parameter 6  ----> fOneOvLamSym");
	AliInfo("Parameter 7  ----> fB");
	AliInfo("Parameter 8  ----> fFsig");
	AliInfo("Parameter 9  ----> fMmean");
	AliInfo("Parameter 10 ----> fNexp");
	AliInfo("Parameter 11 ----> fSigma");
	AliInfo("Parameter 12 ----> fAlpha");
	AliInfo("Parameter 13 ----> fNorm");
	AliInfo("Parameter 14 ----> fBkgNorm");
	AliInfo("Parameter 15 ----> fBkgMean");
	AliInfo("Parameter 16 ----> fBkgSlope");
	AliInfo("Parameter 17 ----> fBkgConst");
	AliInfo("Parameter 18 ----> fNormGaus1");
	AliInfo("Parameter 19 ----> fNormGaus2");

	AliInfo(Form("\n+++\n+++ Number of candidates ---> %d\n+++\n ", ncand));
}
//___________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler& AliDielectronBtoJPSItoEleCDFfitHandler::operator=(const AliDielectronBtoJPSItoEleCDFfitHandler& c)
{
	//
	// Assignment operator
	//
	if (this!=&c) {
		fIsParamFixed = c.fIsParamFixed;
		fPrintStatus  = c.fPrintStatus;
		fUp           = c.fUp;
		fX            = c.fX;
		fM            = c.fM;
		fLikely       = c.fLikely;
		fNcand        = c.fNcand;
	}
	return *this;
}

//_______________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler::AliDielectronBtoJPSItoEleCDFfitHandler(const AliDielectronBtoJPSItoEleCDFfitHandler& c) :
	TNamed(c),
	fIsParamFixed(c.fIsParamFixed),
	fPrintStatus(c.fPrintStatus),
	fUp(c.fUp),
	fX(c.fX),
	fM(c.fM),
	fLikely(c.fLikely),
	fNcand(c.fNcand)
{
	//
	// Copy Constructor
	//
}
//_______________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler::~AliDielectronBtoJPSItoEleCDFfitHandler()
{
	//
	//destructor
	//
	delete fLikely;
}
//_______________________________________________________________________________________
Int_t AliDielectronBtoJPSItoEleCDFfitHandler::DoMinimization()
{
	//
	// performs the minimization
	//
	static TVirtualFitter *fitter = TVirtualFitter::Fitter(this,20);
	fitter->SetFCN(CDFFunction);

	fitter->SetParameter(0,"fWeightRes",fParamStartValues[0], 1.e-08, 0., 1.e+06);
	fitter->SetParameter(1,"fPos",fParamStartValues[1], 1.e-08, 0.,1.e+06);
	fitter->SetParameter(2,"fNeg",fParamStartValues[2], 1.e-08, 0.,1.e+06);
	fitter->SetParameter(3,"fSym",fParamStartValues[3], 1.e-08, 0.,1.e+06);
	fitter->SetParameter(4,"fOneOvLamPlus",fParamStartValues[4], 1.e-10, 0.0000001, 5.e+01);
	fitter->SetParameter(5,"fOneOvLamMinus",fParamStartValues[5], 1.e-10, 0.00000001, 5.e+01);
	fitter->SetParameter(6,"fOneOvLamSym",fParamStartValues[6], 1.e-10, 0.00000001, 5.e+01);
	fitter->SetParameter(7,"fB",fParamStartValues[7], 1.e-10, 0., 1.);
	fitter->SetParameter(8,"fFsig",fParamStartValues[8], 1.e-10, 0., 1.);
	fitter->SetParameter(9,"fMmean",fParamStartValues[9], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(10,"fNexp",fParamStartValues[10], 1.e-08, 0., 1.e+02);
	fitter->SetParameter(11,"fSigma",fParamStartValues[11], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(12,"fAlpha",fParamStartValues[12], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(13,"fNorm",fParamStartValues[13], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(14,"fBkgNorm",fParamStartValues[14], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(15,"fBkgMean",fParamStartValues[15], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(16,"fBkgSlope",fParamStartValues[16], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(17,"fBkgSlope",fParamStartValues[17], 1.e-08, 0., 1.e+04);
	fitter->SetParameter(18,"fNormGaus1",fParamStartValues[18], 1.e-08, 0., 1.e+05);
	fitter->SetParameter(19,"fNormGaus2",fParamStartValues[19], 1.e-08, 0., 1.e+05); 

	for(UInt_t indexparam = 0; indexparam < 20; indexparam++){
		if(IsParamFixed(indexparam))fitter->FixParameter((Int_t)indexparam);
	}

	Double_t arglist[2]={10000,1.0}; 
	Int_t iret=fitter->ExecuteCommand("MIGRAD", arglist ,2);
	fitter->PrintResults(4,0);

	AliInfo("Minimization procedure finished\n");

	return iret;
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::CdfFCN(Int_t & /* npar */, 
		Double_t * /* gin */,Double_t &f,Double_t *par,Int_t /* iflag */)
{
	// 
	// Definition of the FCN to be used by minuit
	//
	fLikely->SetAllParameters(par);
	fLikely->ComputeMassIntegral();
	if(fPrintStatus)fLikely->PrintStatus();

	TStopwatch t;
	t.Start();

	f = fLikely->EvaluateLikelihood(fX,fM,fNcand);

	t.Stop();
	AliDebug(2,Form("Real time spent to calculate function == %f \n", t.RealTime()));
	AliDebug(2,Form("CPU time spent to calculate function == %f \n", t.CpuTime()));
	AliDebug(2,Form("Actual value of the AliDielectronBtoJPSItoEleCDFfitFCN function == %f \n", f));

	return;
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetParamStartValues(Double_t inputparamvalues[20])
{
	for(Int_t index=0; index < 20; index++) fParamStartValues[index] = inputparamvalues[index];
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetResolutionConstants(Double_t* resolutionConst)
{
	//
	// Sets constants for the resolution function
	//
	fLikely->SetResolutionConstants(resolutionConst);

}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetCrystalBallFunction(Bool_t okCB) 
{
	//
	// Sets the CB as the parametrization for the signal invariant mass spectrum 
	// (otherwise Landau is chosen)
	//
	fLikely->SetCrystalBallFunction(okCB);
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetMassWndHigh(Double_t limit) 
{ 
	//
	// Sets upper limit for the invariant mass window (under J/PSI mass peak)
	//
	fLikely->SetMassWndHigh(limit);
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetMassWndLow(Double_t limit)
{
	//
	// Sets lower limit for the invariant mass window (under J/PSI mass peak)
	//
	fLikely->SetMassWndLow(limit);
}

