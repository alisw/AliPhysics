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
#include <TFitter.h>
#include <TMinuit.h>
#include <TStopwatch.h>
#include <TCanvas.h>
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
	fIsParamFixed(45),
	fPrintStatus(kFALSE),
	fUp(0),
	fX(0x0),
	fM(0x0),
	fType(0x0),
        fLikely(0x0),
	fNcand(0),
	fContPlot1(0x0),
	fContPlot2(0x0),
	fContPlot3(0x0),
	fitter(0)
{
	//
	// default constructor
	//
	for (Int_t i=0; i<45; ++i) fParamStartValues[i]=0;
}
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitHandler::AliDielectronBtoJPSItoEleCDFfitHandler(Double_t* decaytime, 
		Double_t* invariantmass, Int_t *type, Int_t ncand) :
	fIsParamFixed(45),
	fPrintStatus(kFALSE),
	fUp(0),
	fX(decaytime),
	fM(invariantmass),
        fType(type),
	fLikely(0x0),
	fNcand(ncand),
	fContPlot1(0x0),
	fContPlot2(0x0),
	fContPlot3(0x0),
	fitter(0)
{
	//
	// constructor
	//
        for (Int_t i=0; i<45; ++i) fParamStartValues[i]=0;
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
	AliInfo("Parameter 18 ----> fNormGaus1 (FF)");
	AliInfo("Parameter 19 ----> fNormGaus2 (FF)");
        AliInfo("Parameter 20 ----> fMean1Res (FF)");
        AliInfo("Parameter 21 ----> fsigma1Res (FF)"); 
        AliInfo("Parameter 22 ----> fMean2Res (FF)");
        AliInfo("Parameter 23 ----> fsigma2Res (FF)");
        AliInfo("Parameter 24 ----> fAlfaRes (FF)");
        AliInfo("Parameter 25 ----> fLambdaRes (FF)");
        AliInfo("Parameter 26 ----> fNormResExp (FF)");
        AliInfo("Parameter 27 ----> fNormGaus1 (FS)");
        AliInfo("Parameter 28 ----> fNormGaus2 (FS)");
        AliInfo("Parameter 29 ----> fMean1Res (FS)");
        AliInfo("Parameter 30 ----> fsigma1Res (FS)");
        AliInfo("Parameter 31 ----> fMean2Res (FS)");
        AliInfo("Parameter 32 ----> fsigma2Res (FS)");
        AliInfo("Parameter 33 ----> fAlfaRes (FS)");
        AliInfo("Parameter 34 ----> fLambdaRes (FS)");
        AliInfo("Parameter 35 ----> fNormResExp (FS)");
        AliInfo("Parameter 36 ----> fNormGaus1 (SS)");
        AliInfo("Parameter 37 ----> fNormGaus2 (SS)");
        AliInfo("Parameter 38 ----> fMean1Res (SS)");
        AliInfo("Parameter 39 ----> fsigma1Res (SS)");
        AliInfo("Parameter 40 ----> fMean2Res (SS)");
        AliInfo("Parameter 41 ----> fsigma2Res (SS)");
        AliInfo("Parameter 42 ----> fAlfaRes (SS)");
        AliInfo("Parameter 43 ----> fLambdaRes (SS)");
        AliInfo("Parameter 44 ----> fNormResExp (SS)");

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
		fType         = c.fType;
                fLikely       = c.fLikely;
		fNcand        = c.fNcand;
		fContPlot1    = c.fContPlot1;
		fContPlot2    = c.fContPlot2;
		fContPlot3    = c.fContPlot3;
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
	fType(c.fType),
        fLikely(c.fLikely),
	fNcand(c.fNcand),
	fContPlot1(c.fContPlot1),
	fContPlot2(c.fContPlot2),
	fContPlot3(c.fContPlot3),
	fitter(c.fitter)
{
	//
	// Copy Constructor
	//
        for (Int_t i=0; i<45; ++i) fParamStartValues[i]=c.fParamStartValues[i];
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
Int_t AliDielectronBtoJPSItoEleCDFfitHandler::DoMinimization(Int_t step)
{
	//
	// performs the minimization
	//
	if(step == 0){ 
		//fitter = TVirtualFitter::Fitter(this,20);
		fitter = (TFitter*)TVirtualFitter::Fitter(this,45);
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
		fitter->SetParameter(17,"fBkgConst",fParamStartValues[17], 1.e-08, 0., 1.e+04);
		fitter->SetParameter(18,"fNormGaus1FF",fParamStartValues[18], 1.e-08, 0., 1.e+05);
		fitter->SetParameter(19,"fNormGaus2FF",fParamStartValues[19], 1.e-08, 0., 1.e+05); 
	        fitter->SetParameter(20,"fMean1ResFF",fParamStartValues[20], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(21,"fSigma1ResFF",fParamStartValues[21], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(22,"fMean2ResFF",fParamStartValues[22], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(23,"fSigma2ResFF",fParamStartValues[23], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(24,"fAlfaResFF",fParamStartValues[24], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(25,"fLambdaResFF",fParamStartValues[25], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(26,"fResNormExpFF",fParamStartValues[26], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(27,"fNormGaus1FS",fParamStartValues[27], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(28,"fNormGaus2FS",fParamStartValues[28], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(29,"fMean1ResFS",fParamStartValues[29], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(30,"fSigma1ResFS",fParamStartValues[30], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(31,"fMean2ResFS",fParamStartValues[31], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(32,"fSigma2ResFS",fParamStartValues[32], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(33,"fAlfaResFS",fParamStartValues[33], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(34,"fLambdaResFS",fParamStartValues[34], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(35,"fResNormExpFS",fParamStartValues[35], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(36,"fNormGaus1SS",fParamStartValues[36], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(37,"fNormGaus2SS",fParamStartValues[37], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(38,"fMean1ResSS",fParamStartValues[38], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(39,"fSigma1ResSS",fParamStartValues[39], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(40,"fMean2ResSS",fParamStartValues[40], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(41,"fSigma2ResSS",fParamStartValues[41], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(42,"fAlfaResSS",fParamStartValues[42], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(43,"fLambdaResSS",fParamStartValues[43], 1.e-08, 0., 1.e+05);
                fitter->SetParameter(44,"fResNormExpSS",fParamStartValues[44], 1.e-08, 0., 1.e+05);
                }

	for(UInt_t indexparam = 0; indexparam < 45; indexparam++){
		if(IsParamFixed(indexparam)) fitter->FixParameter((Int_t)indexparam); 
		else fitter->ReleaseParameter((Int_t)indexparam);
	}
	Double_t arglist[2]={10000,0.1};
	if(step == 2) {Int_t  iret1 = fitter->ExecuteCommand("MINOS", arglist ,1); return iret1;}
	Int_t iret=fitter->ExecuteCommand("MIGRAD", arglist ,2);
	fitter->PrintResults(4,0);

	if(step == 3) {

		TMinuit* minuitPoint = fitter->GetMinuit(); 

		TCanvas *c2 = new TCanvas("c2","contours",800,800);

		//68.27% (1 sigma) confidence level for 2 parameters (fSIG versus fB)   
		minuitPoint->SetErrorDef(1.15); // 2.3/2
		fContPlot1 = (TGraph*)minuitPoint->Contour(100,7,8);
		fContPlot1->GetXaxis()->SetRange(0,1);
                fContPlot1->GetYaxis()->SetRange(0,1);
                fContPlot1->SetLineColor(42);
		fContPlot1->SetLineWidth(3);

		//95% (2 sigma) confidence level for 2 parameters (fSIG versus fB)   
		minuitPoint->SetErrorDef(2.995); // 5.99/2
		fContPlot2 = (TGraph*)minuitPoint->Contour(100,7,8);
                fContPlot2->GetXaxis()->SetRange(0,1);
                fContPlot2->GetYaxis()->SetRange(0,1);
		fContPlot2->SetLineColor(38);
		fContPlot2->SetLineWidth(3);

		//99.73% (3 sigma) confidence level for 2 parameters (fSIG versus fB)   
		minuitPoint->SetErrorDef(5.915); // 11.83/2
		fContPlot3 = (TGraph*)minuitPoint->Contour(100,7,8);
                fContPlot3->GetXaxis()->SetRange(0,1);
                fContPlot3->GetXaxis()->SetTitle("f_{B}");
                fContPlot3->GetYaxis()->SetTitle("f_{Sig}[2.4-4]");
                fContPlot3->GetYaxis()->SetRange(0,1);
  		fContPlot3->SetLineColor(34);
		fContPlot3->SetLineWidth(3);
          
                fContPlot3->Draw("al");
                fContPlot2->Draw("l");
                fContPlot1->Draw("l");
		
                c2->Draw();	
	}

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

	f = fLikely->EvaluateLikelihood(fX,fM,fType,fNcand);

	t.Stop();
	AliDebug(2,Form("Real time spent to calculate function == %f \n", t.RealTime()));
	AliDebug(2,Form("CPU time spent to calculate function == %f \n", t.CpuTime()));
	AliDebug(2,Form("Actual value of the AliDielectronBtoJPSItoEleCDFfitFCN function == %f \n", f));

	return;
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetParamStartValues(Double_t inputparamvalues[45])
{
	for(Int_t index=0; index < 45; index++) fParamStartValues[index] = inputparamvalues[index];
}
//_______________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitHandler::SetResolutionConstants(Double_t* resolutionConst, Int_t type)
{
	//
	// Sets constants for the resolution function
	//
	fLikely->SetResolutionConstants(resolutionConst,type);

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

