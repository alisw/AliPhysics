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
#include "AliLog.h"
#include "AliDielectronBtoJPSItoEleCDFfitFCN.h"
#include "TFormula.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

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
        fWeightType[0] = 1.; fWeightType[1] = 1.; fWeightType[2] = 1.;
        for(Int_t iPar = 0; iPar < 45; iPar++) fParameters[iPar] = 0.;
        fParameters[9] = 1.;fParameters[11] = 1.;fParameters[12] = 1.;
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
	for(Int_t iPar = 0; iPar < 45; iPar++) fParameters[iPar] = source.fParameters[iPar];
        for(Int_t iW=0; iW<2; iW++) fWeightType[iW] = source.fWeightType[iW]; 
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

     	for(Int_t iPar = 0; iPar < 45; iPar++) fParameters[iPar] = source.fParameters[iPar];
	return *this;
}  
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCN::~AliDielectronBtoJPSItoEleCDFfitFCN()
{
	//
	// Default destructor
	//

	delete fhCsiMC;
	for(Int_t iPar = 0; iPar < 45; iPar++) fParameters[iPar] = 0.;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
		const Double_t* invariantmass, const Int_t *type, const Int_t ncand) const
{
	//
	// This function evaluates the Likelihood fnction
	// It returns the -Log(of the likelihood function)
	//
	Double_t f = 0.;
	Double_t ret = 0.;

	for(Int_t i=0; i < ncand; i++) {
		f = EvaluateCDFfuncNorm(pseudoproperdecaytime[i],invariantmass[i],type[i]);
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
	for(Int_t index = 0; index < 45; index++) fParameters[index] = parameters[index];
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
	if(intmMassBkg < 1.e-05) intmMassBkg = 1.;
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
	// resolution func (First-First)
	printf("actual value of norm1Gauss (FF) --------------------------------->> | %f \n", GetNormGaus1ResFunc(2));
	printf("actual value of norm2Gauss (FF) --------------------------------->> | %f \n", GetNormGaus2ResFunc(2));
        printf("actual value of fMean1Res (FF) ---------------------------------->> | %f \n", GetResMean1(2));
        printf("actual value of fSigma1Res (FF) --------------------------------->> | %f \n", GetResSigma1(2));
        printf("actual value of fMean2Res (FF) ---------------------------------->> | %f \n", GetResMean2(2));
        printf("actual value of fSigma2Res (FF) --------------------------------->> | %f \n", GetResSigma2(2));
        
        printf("actual value of alfaRes (FF) ------------------------------------>> | %f \n", GetResAlfa(2));
        printf("actual value of lambdaRes (FF) ---------------------------------->> | %f \n", GetResLambda(2)); 
        printf("actual value of normExpRes (FF) --------------------------------->> | %f \n", GetResNormExp(2)); 
       
        // resolution func (First-Second)
        printf("actual value of norm1Gauss (FS) --------------------------------->> | %f \n", GetNormGaus1ResFunc(1));
        printf("actual value of norm2Gauss (FS) --------------------------------->> | %f \n", GetNormGaus2ResFunc(1));
        printf("actual value of fMean1Res (FS) ---------------------------------->> | %f \n", GetResMean1(1));
        printf("actual value of fSigma1Res (FS) --------------------------------->> | %f \n", GetResSigma1(1));
        printf("actual value of fMean2Res (FS) ---------------------------------->> | %f \n", GetResMean2(1));
        printf("actual value of fSigma2Res (FS) --------------------------------->> | %f \n", GetResSigma2(1));
        
        printf("actual value of alfaRes (FS) ------------------------------------>> | %f \n", GetResAlfa(1));
        printf("actual value of lambdaRes (FS) ---------------------------------->> | %f \n", GetResLambda(1));    
        printf("actual value of normExpRes (FS) --------------------------------->> | %f \n", GetResNormExp(1));    
        
        // resolution func (Second-Second) 
        printf("actual value of norm1Gauss (SS) --------------------------------->> | %f \n", GetNormGaus1ResFunc(0));
        printf("actual value of norm2Gauss (SS) --------------------------------->> | %f \n", GetNormGaus2ResFunc(0));
        printf("actual value of fMean1Res (SS) ---------------------------------->> | %f \n", GetResMean1(0));
        printf("actual value of fSigma1Res (SS) --------------------------------->> | %f \n", GetResSigma1(0));
        printf("actual value of fMean2Res (SS) ---------------------------------->> | %f \n", GetResMean2(0));
        printf("actual value of fSigma2Res (SS) --------------------------------->> | %f \n", GetResSigma2(0));
        
        printf("actual value of alfaRes (SS) ------------------------------------>> | %f \n", GetResAlfa(0));
        printf("actual value of lambdaRes (SS) ---------------------------------->> | %f \n", GetResLambda(0));    
        printf("actual value of normExpRes (SS) --------------------------------->> | %f \n", GetResNormExp(0));    

        printf("\n");
	// integrals constants
	printf("Actual value of normalization integral for MassSig ---------------->> | %f \n", GetIntegralMassSig());
	printf("Actual value of normalization integral for MassBkg ---------------->> | %f \n", GetIntegralMassBkg());

	printf("\n");
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetResolutionConstants(const Double_t* resolutionConst, Int_t type) 
{
	//
	// Resolution function is parametrized as the sum of two gaussian
	//
        Int_t index = (2-type)*9; 
        fParameters[20+index]=resolutionConst[1]; //mean1
        fParameters[22+index]=resolutionConst[4]; //mean2
        fParameters[18+index]=resolutionConst[0]; //norm1
        fParameters[21+index]=resolutionConst[2]; //sigma1
        fParameters[23+index]=resolutionConst[5]; //sigma2
        fParameters[19+index]=resolutionConst[3]; //norm2

        //exp values
        fParameters[24+index]=resolutionConst[6]; //alfa
        fParameters[25+index]=resolutionConst[7]; //lambda
        fParameters[26+index]=resolutionConst[8]; //norm3
        return;
}

//________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfunc(Double_t x, Double_t m, Int_t type) const 
{
        // evaluate likelihood function
	return fParameters[8]*EvaluateCDFfuncSignalPart(x,m,type) + (1. - fParameters[8])*EvaluateCDFfuncBkgPart(x,m,type);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncNorm(Double_t x, Double_t m, Int_t type) const
{
        // evaluate likelihood function
	return EvaluateCDFfunc(x,m,type);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncSignalPart(Double_t x, Double_t m, Int_t type) const 
{
  // evaluate psproper signal	
  return EvaluateCDFDecayTimeSigDistr(x,type)*(EvaluateCDFInvMassSigDistr(m)/fintmMassSig); 
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistr(Double_t x, Int_t type) const
{
	//
	// Implementation of the Background part of the Likelyhood function
	// 
	Double_t retvalue = 0.;
	Double_t funBnorm = FunB(x,type);
	Double_t funPnorm = ResolutionFunc(x,type);
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
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunB(Double_t x, Int_t type) const  
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
		resolutionxdiff = ResolutionFunc(xdiff, type); // normalized value
		sum += csiMCxprime * resolutionxdiff;
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunP(Double_t x, Int_t type) const 
{
	//
	//  Parameterisation of the Prompt part for the x distribution
	//
	return ResolutionFunc(x,type);
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
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncBkgPart(Double_t x,Double_t m, Int_t type) const 
{
	//
	// Return the part of the likelihood function for the background hypothesis
	//
	return EvaluateCDFDecayTimeBkgDistr(x,type)*(EvaluateCDFInvMassBkgDistr(m)/fintmMassBkg);
}  

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistr(Double_t x, Int_t type) const 
{
	//
	// it returns the value of the probability to have a given x for the background 
	//

	Double_t ret = fParameters[0]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*ResolutionFunc(x,type) + fParameters[1]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgPos(x,type) +  fParameters[2]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgNeg(x,type) + fParameters[3]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3])*FunBkgSym(x,type);
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
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgPos(Double_t x,Int_t type) const 
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
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x,type));}
		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x,type));}
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgNeg(Double_t x, Int_t type) const 
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
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x,type));}

		xprime = xupp - (i-.5) * step;
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x,type));}
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgSym(Double_t x, Int_t type) const 
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
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x,type));}
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x,type));}

		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x,type));} 
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x,type));}
	}

	return step*(sum1 + sum2) ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFunc(Double_t x, Int_t type) const  
{
	//
	// parametrization with 2 gaus
	//
        Int_t index = (2-type)*9;
        Double_t mean1 = fParameters[20+index];
        Double_t mean2 = fParameters[22+index];
        Double_t norm1 = fParameters[18+index];
	Double_t sigma1 = fParameters[21+index];
	Double_t sigma2 = fParameters[23+index];
        Double_t norm2 = fParameters[19+index];
        //exp values
        Double_t alfa = fParameters[24+index];
        Double_t lambda = fParameters[25+index];
        Double_t norm3 = fParameters[26+index];
 
        Double_t ret = 0.; Double_t fitval = 0; 
        if(TMath::Abs(x)<=alfa) fitval = (lambda-1)/(2*alfa*lambda);
        else  fitval = ((lambda-1)/(2*alfa*lambda))*TMath::Power(alfa,lambda)*(TMath::Power(TMath::Abs(x),-1*lambda));

        ret = (norm1/(norm1+norm2+norm3))*((1/(sigma1*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean1)/sigma1)*((x-mean1)/sigma1))) + (norm2/(norm1+norm2+norm3))*((1/(sigma2*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean2)/sigma2)*((x-mean2)/sigma2))) + (norm3/(norm1+norm2+norm3))*fitval;

        return ret;
}     

//_________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetCsiMC(Double_t xmin, Double_t xmax, Double_t normalization) 
{
	// return the pointer to the templateMC function 
	TF1* templateMC = new TF1("MCtemplate",this,&AliDielectronBtoJPSItoEleCDFfitFCN::CsiMCfunc,xmin,xmax,1);
	templateMC->SetParameter(0,normalization);
        templateMC->SetNpx(5000);
        return (TF1*)templateMC->Clone();
}

//__________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetResolutionFunc(Double_t xmin, Double_t xmax, Double_t normalization, Double_t type){
	// return the pointer to the resolution function
	TF1* resFunc = new TF1(Form("resolutionFunc_%f",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFuncf,xmin,xmax,2);
        resFunc->SetParameter(0,normalization);
        resFunc->SetParameter(1,type);
	resFunc->SetNpx(5000);
        return (TF1*)resFunc->Clone();
}

//__________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetResolutionFuncAllTypes(Double_t xmin, Double_t xmax, Double_t normalization){
        // return the pointer to the resolution function
        TF1* resFunc = new TF1(Form("resolutionFunc"),this,&AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFuncAllTypes,xmin,xmax,1);
        resFunc->SetParameter(0,normalization);
        resFunc->SetNpx(5000);
        return (TF1*)resFunc->Clone();
}

//___________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeBkgDistr(Double_t xmin, Double_t xmax, Double_t normalization, Double_t type){
	// return the pointer to the background x distribution function
	TF1 *backFunc = new TF1(Form("backFunc_%f",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrFunc,xmin,xmax,2);
        backFunc->SetParameter(0,normalization);
        backFunc->SetParameter(1,type);
	backFunc->SetNpx(5000);
        return (TF1*)backFunc->Clone();
}

//___________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeBkgDistrAllTypes(Double_t xmin, Double_t xmax, Double_t normalization){
        // return the pointer to the background x distribution function
        TF1 *backFunc = new TF1(Form("backFunc"),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrFuncAllTypes,xmin,xmax,1);
        backFunc->SetParameter(0,normalization);
        backFunc->SetNpx(5000);
        return (TF1*)backFunc->Clone();
}

//__________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeSigDistr(Double_t xmin, Double_t xmax, Double_t normalization, Double_t type){
	// return the pointer to the signal x distribution function
	TF1 *signFunc = new TF1("signalFunc",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistrFunc,xmin,xmax,2);
        signFunc->SetParameter(0,normalization);
        signFunc->SetParameter(1,type); 
        signFunc->SetNpx(5000);
	return (TF1*)signFunc->Clone();
}

//____________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFInvMassBkgDistr(Double_t mMin, Double_t mMax, Double_t normalization){
  // return the pointer to the invariant mass distribution for the background 
  TF1 * invMassBkg = new TF1("invMassBkg",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassBkgDistrFunc,mMin,mMax,1);
  invMassBkg->SetParameter(0,normalization);
  invMassBkg->SetNpx(5000);
  return (TF1*)invMassBkg->Clone();
}


//____________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFInvMassSigDistr(Double_t mMin, Double_t mMax, Double_t normalization){
  // return the pointer to the invariant mass distribution for the signal
  TF1 * invMassSig = new TF1("invMassSig",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassSigDistrFunc,mMin,mMax,1);
  invMassSig->SetParameter(0,normalization);
  invMassSig->SetNpx(5000);
  return (TF1*)invMassSig->Clone();
}

//____________________________________________________________________________________________________
TF1 *AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFInvMassTotalDistr(Double_t mMin, Double_t mMax, Double_t normalization){
  // return the pointer to the invariant mass distribution
  TF1 * invMassTot = new TF1("invMassTot",this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassTotalDistr,mMin,mMax,1);
  invMassTot->SetParameter(0,normalization);
  invMassTot->SetNpx(5000);
  return (TF1*)invMassTot->Clone();
}

//____________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassTotalDistr(const Double_t* x, const Double_t *par) const
{
  // evaluate invariant mass total distribution
  Double_t value = 0;
  Double_t xx = x[0];
  value = ((EvaluateCDFInvMassSigDistr(xx)/fintmMassSig)*fParameters[8] + (1-fParameters[8])*(EvaluateCDFInvMassBkgDistr(xx)/fintmMassBkg))*par[0];
  return value;  
}

//____________________________________________________________________________________________________
TF1 *AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeTotalDistr(Double_t xMin, Double_t xMax, Double_t normalization,Double_t type){
 // return the pointer to the pseudoproper distribution for the background
 TF1 *decayTimeTot = new TF1(Form("decayTimeTot_%f",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistr,xMin,xMax,2);
 decayTimeTot->SetParameter(0,normalization);
 decayTimeTot->SetParameter(1,type);
 decayTimeTot->SetNpx(5000);
 return (TF1*)decayTimeTot->Clone();
}

//____________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistr(const Double_t* x, const Double_t *par) const
{
 // evaluate the total pseudoproper distribution for a given candidate's type. par[1] should be the candidate's type.
 Double_t value = 0;
 Double_t xx = x[0];
 value = (fParameters[8]*EvaluateCDFDecayTimeSigDistr(xx,(Int_t)par[1]) + (1-fParameters[8])*EvaluateCDFDecayTimeBkgDistr(xx,(Int_t)par[1]))*par[0];
 return value;
}

//____________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistrAllTypes(const Double_t* x, const Double_t *par) const
{
 // evaluate the total pseudoproper distribution considering all candidate's types
 Double_t value = 0;
 Double_t xx = x[0];

 value = (fParameters[8]*(fWeightType[2]*EvaluateCDFDecayTimeSigDistr(xx,2)+fWeightType[1]*EvaluateCDFDecayTimeSigDistr(xx,1)+fWeightType[0]*EvaluateCDFDecayTimeSigDistr(xx,0)))+((1-fParameters[8])*(fWeightType[2]*EvaluateCDFDecayTimeBkgDistr(xx,2) + fWeightType[1]*EvaluateCDFDecayTimeBkgDistr(xx,1)+fWeightType[0]*EvaluateCDFDecayTimeBkgDistr(xx,0))); 

 return value*par[0];
}

//____________________________________________________________________________________________________
TF1 *AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeTotalDistrAllTypes(Double_t xMin, Double_t xMax, Double_t normalization){
 // return the pointer to the pseudoproper function for the background considering all candidate's types
 TF1 *decayTimeTot = new TF1(Form("decayTimeTot"),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistrAllTypes,xMin,xMax,1);
 decayTimeTot->SetParameter(0,normalization);
 decayTimeTot->SetNpx(5000);
 return (TF1*)decayTimeTot->Clone();
}


//____________________________________________________________________________________________________
TF1 * AliDielectronBtoJPSItoEleCDFfitFCN::GetFunB(Double_t xmin, Double_t xmax, Double_t normalization, Double_t type){
 // return the pointer to the function that describe secondary jpsi pseudoproper distribution for a given candidate's type
 TF1* funb = new TF1(Form("secondaryJpsiConvolution_%f",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::FunBfunc,xmin,xmax,2);
 funb->SetParameter(0,normalization);
 funb->SetParameter(1,type);
 funb->SetNpx(5000);
 return (TF1*)funb->Clone();
 }

//____________________________________________________________________________________________________
TF1 * AliDielectronBtoJPSItoEleCDFfitFCN::GetFunBAllTypes(Double_t xmin, Double_t xmax, Double_t normalization){
// return the pointer to the function that describe secondary jpsi pseudoproper distribution for all types
 TF1* funb = new TF1(Form("secondaryJpsiConvolution"),this,&AliDielectronBtoJPSItoEleCDFfitFCN::FunBfuncAllTypes,xmin,xmax,1);
 funb->SetParameter(0,normalization);
 funb->SetNpx(5000);
 return (TF1*)funb->Clone();
 }


