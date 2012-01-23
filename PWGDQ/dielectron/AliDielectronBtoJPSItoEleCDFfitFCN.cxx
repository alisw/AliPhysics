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
#include <TFormula.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include <AliLog.h>

#include "AliDielectronBtoJPSItoEleCDFfitFCN.h"

//_________________________________________________________________________
//                        Class AliDielectronBtoJPSItoEleCDFfitFCN
//                   Definition of main function used in 
//                     unbinned log-likelihood fit for
//                 the channel B -> JPsi + X -> e+e- + X
//      
//                           Origin: C.Di Giglio
//       Contact: Carmelo.Digiglio@ba.infn.it , Giuseppe.Bruno@ba.infn.it
//_________________________________________________________________________

ClassImp(AliDielectronBtoJPSItoEleCDFfitFCN)

	//_________________________________________________________________________________________________
	AliDielectronBtoJPSItoEleCDFfitFCN::AliDielectronBtoJPSItoEleCDFfitFCN() :
		fFPlus(0.),
		fFMinus(0.),
		fFSym(0.),
		fintmMassSig(1.),
		fintmMassBkg(1.),
		fhCsiMC(0x0),
		fMassWndHigh(0.),
		fMassWndLow(0.),
		fCrystalBallParam(kFALSE)
{
	//
	// constructor
	//
	SetCrystalBallFunction(kFALSE);
	SetMassWndHigh(0.2);
	SetMassWndLow(0.5);
	for(Int_t iPar = 0; iPar < 20; iPar++) fParameters[iPar] = 0.;
	fParameters[9] = 1.;fParameters[11] = 1.;fParameters[12] = 1.;
	for(Int_t index=0; index<4; index++) fResolutionConstants[index] = 0.;
	AliInfo("Instance of AliDielectronBtoJPSItoEleCDFfitFCN-class created");
}
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCN::AliDielectronBtoJPSItoEleCDFfitFCN(const AliDielectronBtoJPSItoEleCDFfitFCN& source) :
	TNamed(source),
	fFPlus(source.fFPlus),
	fFMinus(source.fFMinus),
	fFSym(source.fFSym),
	fintmMassSig(source.fintmMassSig),
	fintmMassBkg(source.fintmMassBkg),
	fhCsiMC(source.fhCsiMC),
	fMassWndHigh(source.fMassWndHigh),
	fMassWndLow(source.fMassWndLow),
	fCrystalBallParam(source.fCrystalBallParam)
{
	//
	// Copy constructor
	//
	for(Int_t iPar = 0; iPar < 20; iPar++) fParameters[iPar] = source.fParameters[iPar];
	for(Int_t index=0; index<4; index++) fResolutionConstants[index] = source.fResolutionConstants[index];
}
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCN& AliDielectronBtoJPSItoEleCDFfitFCN::operator=(const AliDielectronBtoJPSItoEleCDFfitFCN& source) 
{
	//
	// Assignment operator
	//
	if(&source == this) return *this;
	fFPlus = source.fFPlus;
	fFMinus = source.fFMinus;
	fFSym = source.fFSym;
	fintmMassSig = source.fintmMassSig;
	fintmMassBkg = source.fintmMassBkg;
	fhCsiMC = source.fhCsiMC;
	fCrystalBallParam = source.fCrystalBallParam;

	for(Int_t iPar = 0; iPar < 20; iPar++) fParameters[iPar] = source.fParameters[iPar];
	for(Int_t index=0; index<4; index++) fResolutionConstants[index] = source.fResolutionConstants[index];

	return *this;
}  
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCN::~AliDielectronBtoJPSItoEleCDFfitFCN()
{
	//
	// Default destructor
	//

	delete fhCsiMC;
	for(Int_t iPar = 0; iPar < 20; iPar++) fParameters[iPar] = 0.;
	for(Int_t index=0; index<4; index++) fResolutionConstants[index] = 0.;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
		const Double_t* invariantmass, const Int_t ncand) const
{
	//
	// This function evaluates the Likelihood fnction
	// It returns the -Log(of the likelihood function)
	//
	Double_t f = 0.;
	Double_t ret = 0.;

	for(Int_t i=0; i < ncand; i++) {
		f = EvaluateCDFfuncNorm(pseudoproperdecaytime[i],invariantmass[i]);
		if(f <= 0.) continue;  
		ret += -1*TMath::Log(f);  
	}
	return ret;
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetAllParameters(const Double_t* parameters)
{ 
	//
	// Sets array of FCN parameters
	//
	for(Int_t index = 0; index < 20; index++) fParameters[index] = parameters[index];
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::ComputeMassIntegral() 
{ 
	//
	// this function compute the integral of the likelihood function 
	// (theoretical function) in order to normalize it to unity
	//
	Double_t npm = 20000.;
	Double_t stepm;
	Double_t mx=0.;
	stepm = (fMassWndHigh-fMassWndLow)/npm; 
	// compute integrals for  invariant mass terms        

	Double_t iMassSig;
	Double_t intmMassSig = 0.0;
	Double_t summMassSig = 0.0;
	for(iMassSig = 1.0;  iMassSig<= npm/2.; iMassSig++)  {
		mx = fMassWndLow + (iMassSig - .5)*stepm;
		summMassSig += EvaluateCDFInvMassSigDistr(mx);
		mx = fMassWndHigh - (iMassSig - .5)*stepm;
		summMassSig += EvaluateCDFInvMassSigDistr(mx);
	}
	intmMassSig = summMassSig*stepm;
	SetIntegralMassSig(intmMassSig);
	//
	Double_t iMassBkg;
	Double_t intmMassBkg = 0.0;
	Double_t summMassBkg = 0.0;
	for(iMassBkg = 1.0; iMassBkg <= npm/2.; iMassBkg++)  {
		mx = fMassWndLow + (iMassBkg - .5)*stepm;
		summMassBkg += EvaluateCDFInvMassBkgDistr(mx);
		mx = fMassWndHigh - (iMassBkg - .5)*stepm;
		summMassBkg += EvaluateCDFInvMassBkgDistr(mx);
	}
	intmMassBkg = summMassBkg*stepm;
	SetIntegralMassBkg(intmMassBkg);
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::PrintStatus()
{
	//
	//  Print the parameters of the fits 
	//
	printf("\n");
	// background param
	printf("actual value of fWeightRes------------------------------------->> | %f \n", GetResWeight());
	printf("actual value of fPos ------------------------------------------>> | %f \n", GetFPlus());
	printf("actual value of fNeg ------------------------------------------>> | %f \n", GetFMinus());
	printf("actual value of fSym ------------------------------------------>> | %f \n", GetFSym()); 
	printf("actual value of fOneOvLamPlus --------------------------------->> | %f \n", GetLamPlus());
	printf("actual value of fOneOvLamMinus -------------------------------->> | %f \n", GetLamMinus());
	printf("actual value of fOneOvLamSym ---------------------------------->> | %f \n", GetLamSym());
	printf("actual value of fFractionJpsiFromBeauty ----------------------->> | %f \n", GetFractionJpsiFromBeauty());
	printf("actual value of fFsig ----------------------------------------->> | %f \n", GetFsig());

	if(fCrystalBallParam){
		printf("actual value of fCrystalBallMmean ----------------------------->> | %f \n", GetCrystalBallMmean());
		printf("actual value of fCrystalBallNexp ------------------------------>> | %f \n", GetCrystalBallNexp());
		printf("actual value of fCrystalBallSigma ----------------------------->> | %f \n", GetCrystalBallSigma());
		printf("actual value of fCrystalBallAlpha ----------------------------->> | %f \n", GetCrystalBallAlpha());
		printf("actual value of fCrystalBallNorm  ----------------------------->> | %f \n", GetCrystalBallNorm());
	}else{
		printf("actual value of fMpv ------------------------------------------>> | %f \n", GetCrystalBallMmean());
		printf("actual value of fConstRovL ------------------------------------>> | %f \n", GetCrystalBallNexp());
		printf("actual value of fSigmaL --------------------------------------->> | %f \n", GetCrystalBallSigma());
		printf("actual value of fSigmaR --------------------------------------->> | %f \n", GetCrystalBallAlpha());
	}

	// back Mass func
	printf("actual value of normBkg ----------------------------------------->> | %f \n", GetBkgInvMassNorm());
	printf("actual value of meanBkg ----------------------------------------->> | %f \n", GetBkgInvMassMean());
	printf("actual value of slopeBkg ---------------------------------------->> | %f \n", GetBkgInvMassSlope());
	printf("actual value of constBkg ---------------------------------------->> | %f \n", GetBkgInvMassConst());
	// resolution func
	printf("actual value of norm1Gauss -------------------------------------->> | %f \n",GetNormGaus1ResFunc());
	printf("actual value of norm2Gauss -------------------------------------->> | %f \n",GetNormGaus2ResFunc());

	printf("\n");
	// integrals constants
	printf("Actual value of normalization integral for MassSig ---------------->> | %f \n", GetIntegralMassSig());
	printf("Actual value of normalization integral for MassBkg ---------------->> | %f \n", GetIntegralMassBkg());

	printf("\n");
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetResolutionConstants(Double_t* resolutionConst)
{
	//
	// Resolution function is parametrized as the sum of two gaussian
	//
	fResolutionConstants[0]  = resolutionConst[0]; // mean 1
	fResolutionConstants[1]  = resolutionConst[1]; // sigma 1
	fResolutionConstants[2]  = resolutionConst[2]; // mean 2
	fResolutionConstants[3]  = resolutionConst[3]; // sigma 2
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfunc(Double_t x, Double_t m) const 
{
	return fParameters[8]*EvaluateCDFfuncSignalPart(x,m) + (1. - fParameters[8])*EvaluateCDFfuncBkgPart(x,m);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncNorm(Double_t x, Double_t m) const
{
	return EvaluateCDFfunc(x,m);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncSignalPart(Double_t x, Double_t m) const 
{
	return EvaluateCDFDecayTimeSigDistr(x)*(EvaluateCDFInvMassSigDistr(m)/fintmMassSig); 
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistr(Double_t x) const
{
	//
	// Implementation of the Background part of the Likelyhood function
	// 

	Double_t retvalue = 0.;
	Double_t funBnorm = FunB(x);
	Double_t funPnorm = ResolutionFunc(x);
	retvalue = fParameters[7]*funBnorm + (1. - fParameters[7])*funPnorm;
	return retvalue;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassSigDistr(Double_t m) const
{ 
	//
	// Parametrization of signal part invariant mass distribution
	// It can be either Crystal Ball function or sum of two Landau
	//

	Double_t fitval = 0.;

	if(fCrystalBallParam){
		Double_t t = (m-fParameters[9])/fParameters[11]; ;
		if (fParameters[12] < 0) t = -t;

		Double_t absAlpha = TMath::Abs((Double_t)fParameters[12]);

		if (t >= -absAlpha) {
			return fParameters[13]*TMath::Exp(-0.5*t*t);
		}
		else {
			Double_t a =  TMath::Power(fParameters[10]/absAlpha,fParameters[10])* TMath::Exp(-0.5*absAlpha*absAlpha);
			Double_t b= fParameters[10]/absAlpha - absAlpha;
			fitval = (fParameters[13]*a/TMath::Power(b - t, fParameters[10]));
			return fitval;
		}
	}else{
		Double_t t=-1*m;
		Double_t tmpv=-1*fParameters[9];
		fitval=TMath::Sqrt(TMath::Landau(t,tmpv,fParameters[11]));
		fitval += fParameters[10]*(TMath::Landau(m,fParameters[9],fParameters[12]));
		return fitval;
	}
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunB(Double_t x) const  
{
	//  
	// Parameterisation of the fit function for the x part of the Background
	//
	Double_t np = 1000.0;
	Double_t sc = 10.;
	Double_t sigma3 = 1000.; // valore usato nella macro
	Double_t xprime;
	Double_t sum = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;
	xlow = x - sc * sigma3 ;
	xupp = x + sc * sigma3 ;
	step = (xupp-xlow) / np;
	Double_t csiMCxprime = 0.;
	Double_t resolutionxdiff = 0.;
	Double_t xdiff = 0.;

	for(i=1.0; i<=np; i++){
		xprime = xlow + (i-.5) * step;
		csiMCxprime = CsiMC(xprime);
		xdiff = xprime - x;
		resolutionxdiff = ResolutionFunc(xdiff); // normalized value
		sum += csiMCxprime * resolutionxdiff;
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunP(Double_t x) const 
{
	//
	//  Parameterisation of the Prompt part for the x distribution
	//
	return ResolutionFunc(x);
}


//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::CsiMC(Double_t x) const 
{
	//
	//  Distribution (template) of the x distribution for the x variable 
	//  for the J/psi coming from Beauty hadrons
	//
	Double_t returnvalue = 0.;

	if((fhCsiMC->FindBin(x) > 0) && (fhCsiMC->FindBin(x) < fhCsiMC->GetNbinsX()+1))  
		returnvalue = fhCsiMC->Interpolate(x);

	return returnvalue;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncBkgPart(Double_t x,Double_t m) const 
{
	//
	// Return the part of the likelihood function for the background hypothesis
	//
	return EvaluateCDFDecayTimeBkgDistr(x)*(EvaluateCDFInvMassBkgDistr(m)/fintmMassBkg);
}  

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistr(Double_t x) const 
{
	//
	// it returns the value of the probability to have a given x for the background 
	//

	Double_t ret = fParameters[0]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*ResolutionFunc(x) + fParameters[1]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgPos(x) +  fParameters[2]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgNeg(x) + fParameters[3]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgSym(x);
	return ret;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassBkgDistr(Double_t m) const 
{
	//
	// it returns the value of the probability to have a given mass for the background
	//
	Double_t value = 0.;
	value = fParameters[14]*TMath::Exp(-1*(m-fParameters[15])/fParameters[16]) + fParameters[17];  
	return value;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgPos(Double_t x) const 
{
	//
	// exponential with positive slopes for the background part (x)
	//

	Double_t np = 1000.0;
	Double_t sc = 10.;      
	Double_t sigma3 = 1000.; // valore usato nella macro 
	Double_t xprime;
	Double_t sum = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;
	xlow = x - sc * sigma3 ;
	xupp = x + sc * sigma3 ;
	step = (xupp-xlow) / np;

	for(i=1.0; i<=np/2; i++) {
		xprime = xlow + (i-.5) * step;
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x));}
		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x));}
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgNeg(Double_t x) const 
{
	//
	// exponential with negative slopes for the background part (x)
	//
	Double_t np = 1000.0;
	Double_t sc = 10.;      
	Double_t sigma3 = 1000.; 
	Double_t xprime;
	Double_t sum = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;
	xlow = x - sc * sigma3 ;
	xupp = x + sc * sigma3 ;
	step = (xupp-xlow) / np;

	for(i=1.0; i<=np/2; i++) {

		xprime = xlow + (i-.5) * step;
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x));}

		xprime = xupp - (i-.5) * step;
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x));}
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgSym(Double_t x) const 
{
	//
	// exponential with both positive and negative slopes for the background part (x)
	//
	Double_t np = 1000.0;
	Double_t sc = 10.;      
	Double_t sigma3 = 1000.; 
	Double_t xprime;
	Double_t sum1 = 0.0;
	Double_t sum2 = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;
	xlow = x - sc * sigma3 ;
	xupp = x + sc * sigma3 ;
	step = (xupp-xlow) / np;

	for(i=1.0; i<=np/2; i++) {

		xprime = xlow + (i-.5) * step;
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x));}
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x));}

		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x));} 
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x));}
	}

	return step*(sum1 + sum2) ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFunc(Double_t x) const 
{
	//
	// parametrization with 2 gaus
	//
	Double_t ret = 0.;
	Double_t mean1 = fResolutionConstants[0]; 
	Double_t mean2 = fResolutionConstants[2];
	Double_t norm1 = fParameters[18];
	Double_t sigma1 = fResolutionConstants[1];
	Double_t sigma2 = fResolutionConstants[3]; 
	Double_t norm2 = fParameters[19];

	ret = (norm1/(norm1+norm2))*((1/(sigma1*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean1)/sigma1)*((x-mean1)/sigma1)))+(norm2/(norm1+norm2))*((1/(sigma2*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean2)/sigma2)*((x-mean2)/sigma2))); 

	return ret;
}

//_________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetCsiMC(Double_t xmin, Double_t xmax) 
{
	// return the pointer to the templateMC function 
	TF1* templateMC = new TF1("MCtemplate",this,&AliDielectronBtoJPSItoEleCDFfitFCN::CsiMCfunc,xmin,xmax,0);
	templateMC->SetNpx(5000);
        return (TF1*)templateMC->Clone();
}

//__________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetResolutionFunc(Double_t xmin, Double_t xmax){
	// return the pointer to the resolution function
	TF1* resFunc = new TF1("resolutionFunc",this,&AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFuncf,xmin,xmax,0);
	resFunc->SetNpx(5000);
        return (TF1*)resFunc->Clone();
}

//___________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeBkgDistr(Double_t xmin, Double_t xmax){
	// return the pointer to the background x distribution function
	TF1 *backFunc = new TF1("backFunc",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrFunc,xmin,xmax,0);
	backFunc->SetNpx(5000);
        return (TF1*)backFunc->Clone();
}

//__________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeSigDistr(Double_t xmin, Double_t xmax){
	// return the pointer to the signal x distribution function
	TF1 *signFunc = new TF1("signalFunc",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistrFunc,xmin,xmax,0);
        signFunc->SetNpx(5000);
	return (TF1*)signFunc->Clone();
}
