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
#include <AliLog.h>
#include <TFormula.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFile.h>

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
                fShiftTemplate(0.),
		fMassWndHigh(0.),
		fMassWndLow(0.),
		fCrystalBallParam(kFALSE),
                fChangeResolution(1.),
                fChangeMass(1.),
		fWeights(0),
                fLoadFunctions(kFALSE),
                fMultivariate(kFALSE),
                fFunBSaved(0x0),
                fFunBkgSaved(0x0),
                fResParams(0x0),
                fBkgParams(0x0),
		fMassWindows(0x0),
		fPtWindows(0x0),
		fExponentialParam(kTRUE),
		fSignalBinForExtrapolation(0)
{
	//
	// constructor
	//
	SetMassWndHigh(0.2);
	SetMassWndLow(0.5);
        fWeightType[0] = 1.; fWeightType[1] = 1.; fWeightType[2] = 1.;
        for(Int_t iPar = 0; iPar < 49; iPar++) fParameters[iPar] = 0.;
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
        fShiftTemplate(source.fShiftTemplate),
        fMassWndHigh(source.fMassWndHigh),
        fMassWndLow(source.fMassWndLow),
        fCrystalBallParam(source.fCrystalBallParam),
        fChangeResolution(source.fChangeResolution),
        fChangeMass(source.fChangeMass),
        fWeights(source.fWeights),
        fLoadFunctions(source.fLoadFunctions),
        fMultivariate(source.fMultivariate),
        fFunBSaved(source.fFunBSaved),
        fFunBkgSaved(source.fFunBkgSaved),
        fResParams(source.fResParams),
        fBkgParams(source.fBkgParams),
        fMassWindows(source.fMassWindows),
        fPtWindows(source.fPtWindows),
        fExponentialParam(source.fExponentialParam),
        fSignalBinForExtrapolation(source.fSignalBinForExtrapolation)
{
	//
	// Copy constructor
	//
	for(Int_t iPar = 0; iPar < 49; iPar++) fParameters[iPar] = source.fParameters[iPar];
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
        fLoadFunctions = source.fLoadFunctions;
        fMultivariate = source.fMultivariate;
        fFunBSaved = source.fFunBSaved;
        fFunBkgSaved = source.fFunBkgSaved;
        fResParams = source.fResParams;
        fBkgParams = source.fBkgParams;
        fMassWindows = source.fMassWindows;
        fPtWindows = source.fPtWindows;
        fShiftTemplate = source.fShiftTemplate;
	fCrystalBallParam = source.fCrystalBallParam;
	fExponentialParam = source.fExponentialParam;

     	for(Int_t iPar = 0; iPar < 49; iPar++) fParameters[iPar] = source.fParameters[iPar];
	return *this;
}  
//_________________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCN::~AliDielectronBtoJPSItoEleCDFfitFCN()
{
	//
	// Default destructor
	//

	delete fhCsiMC;
        delete fFunBSaved;
        delete fFunBkgSaved;
	for(Int_t iPar = 0; iPar < 49; iPar++) fParameters[iPar] = 0.;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
		const Double_t* invariantmass, const Double_t *pt, const Int_t *type, Int_t ncand) const
{
	//
	// This function evaluates the Likelihood fnction
	// It returns the -Log(of the likelihood function)
	//
	Double_t f = 0.;
	Double_t ret = 0.;

	for(Int_t i=0; i < ncand; i++) {
		f = EvaluateCDFfuncNorm(pseudoproperdecaytime[i],invariantmass[i],pt[i],type[i]);
		if(f <= 0.) continue;   
                ret += -2.*TMath::Log(f);  
	}
        return ret;
}
//_________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetAllParameters(const Double_t* parameters)
{ 
	//
	// Sets array of FCN parameters
	//
	for(Int_t index = 0; index < 49; index++) fParameters[index] = parameters[index];
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
	printf("actual value of fSym1 ----------------------------------------->> | %f \n", GetFSym1()); 
	printf("actual value of fOneOvLamPlus --------------------------------->> | %f \n", GetLamPlus());
	printf("actual value of fOneOvLamMinus -------------------------------->> | %f \n", GetLamMinus());
	printf("actual value of fOneOvLamSym ---------------------------------->> | %f \n", GetLamSym());
	printf("actual value of fOneOvLamSym1 --------------------------------->> | %f \n", GetLamSym1());
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
        if(fExponentialParam){
		printf("actual value of normBkg ----------------------------------------->> | %f \n", GetBkgInvMassNorm());
		printf("actual value of meanBkg ----------------------------------------->> | %f \n", GetBkgInvMassMean());
		printf("actual value of slopeBkg ---------------------------------------->> | %f \n", GetBkgInvMassSlope());
		printf("actual value of constBkg ---------------------------------------->> | %f \n", GetBkgInvMassConst());
        }else{
		printf("actual value of m^{0} ------------------------------------------->> | %f \n", GetBkgInvMassNorm());
		printf("actual value of m^{1} ------------------------------------------->> | %f \n", GetBkgInvMassMean());
		printf("actual value of m^{2} ------------------------------------------->> | %f \n", GetBkgInvMassSlope());
		printf("actual value of m^{3} ------------------------------------------->> | %f \n", GetBkgInvMassConst());
		printf("actual value of m^{4} ------------------------------------------->> | %f \n", GetPolyn4());
		printf("actual value of m^{5} ------------------------------------------->> | %f \n", GetPolyn5());
        }

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
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfunc(Double_t x, Double_t m, Double_t pt, Int_t type) const 
{
        // evaluate likelihood function
	//printf("CDF func ---> x = %f m = %f pt = %f type = %d \n",x,m,pt,type);
	return fParameters[8]*EvaluateCDFfuncSignalPart(x,m, pt, type) + (1. - fParameters[8])*EvaluateCDFfuncBkgPart(x,m, pt, type);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncNorm(Double_t x, Double_t m, Double_t pt, Int_t type) const
{
        // evaluate likelihood function
	return EvaluateCDFfunc(x,m,pt, type);
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncSignalPart(Double_t x, Double_t m, Double_t pt, Int_t type) const 
{
  // evaluate psproper signal	
  return EvaluateCDFDecayTimeSigDistr(x,pt, type)*(EvaluateCDFInvMassSigDistr(m)/fintmMassSig); 
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistr(Double_t x, Double_t pt, Int_t type) const
{
	//
	// Implementation of the Background part of the Likelyhood function
	// 
	Double_t retvalue = 0.;
	Double_t funBnorm =  fMultivariate ? FunBsaved(x, pt, type) : FunB(x,pt, type)  ;
        Double_t funPnorm = ResolutionFunc(x, pt, type);
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
        // change inv Mass RMS fChangeMass 
	if(fCrystalBallParam){
		Double_t t = ((m-fParameters[9])/fChangeMass)/fParameters[11]; ;
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
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunB(Double_t x, Double_t pt, Int_t type) const  
{
	//  
	// Parameterisation of the fit function for the x part of the Background
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
	Double_t csiMCxprime = 0.;
	Double_t resolutionxdiff = 0.;
	Double_t xdiff = 0.;

	for(i=1.0; i<=np; i++){
		xprime = xlow + (i-.5) * step;
		csiMCxprime = CsiMC(xprime);
		xdiff = xprime - x;
		resolutionxdiff = ResolutionFunc(xdiff, pt, type); // normalized value
		sum += csiMCxprime * resolutionxdiff;
	}
     
	return step * sum ;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunP(Double_t x, Double_t pt, Int_t type) const 
{
	//
	//  Parameterisation of the Prompt part for the x distribution
	//
	return ResolutionFunc(x, pt, type);
}


//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::CsiMC(Double_t x) const 
{
	//
	//  Distribution (template) of the x distribution for the x variable 
	//  for the J/psi coming from Beauty hadrons
	//
	Double_t returnvalue = 0.;
       
	if((fhCsiMC->FindBin(x-fShiftTemplate) > 0) && (fhCsiMC->FindBin(x-fShiftTemplate) < fhCsiMC->GetNbinsX()+1))  
	returnvalue = fhCsiMC->Interpolate(x-fShiftTemplate);


	return returnvalue;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncBkgPart(Double_t x,Double_t m, Double_t pt, Int_t type) const 
{
	//
	// Return the part of the likelihood function for the background hypothesis
	//
         Double_t bkgValx = fMultivariate ? EvaluateCDFDecayTimeBkgDistrSaved(x,type,m,pt) : EvaluateCDFDecayTimeBkgDistr(x,type,m,pt);
        return bkgValx*(EvaluateCDFInvMassBkgDistr(m)/fintmMassBkg);
}
  

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistr(Double_t x, Int_t type, Double_t m, Double_t pt) const
{
        //
        // it returns the value of the probability to have a given x for the background 
        // in the pt, m , type correspondent range 
        //
	Double_t ret = 0.;
	if(fMultivariate) 
	ret = EvaluateCDFDecayTimeBkgDistrDifferential(x,type,m,pt);
	else ret = fParameters[0]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*ResolutionFunc(x, pt, type) + fParameters[1]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgPos(x, pt,type) +  fParameters[2]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgNeg(x,pt,type) + fParameters[3]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgSym(x, pt,type) + fParameters[46]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgSym1(x,pt,type);
	return ret;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassBkgDistr(Double_t m) const 
{
	//
	// it returns the value of the probability to have a given mass for the background
	//
	Double_t value = 0.;
	if(fExponentialParam) 
        value = fParameters[14]*TMath::Exp(-1*(m-fParameters[15])/fParameters[16]) + fParameters[17];  
        else value = fParameters[14] + fParameters[15]*m + fParameters[16]*m*m + fParameters[17]*m*m*m + fParameters[47]*m*m*m*m + fParameters[48]*m*m*m*m*m; 
        return value;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgPos(Double_t x, Double_t pt, Int_t type) const 
{
        //
	// exponential with positive slopes for the background part (x)
	//
	Double_t np, sc, sigma3;
	sc = 10.;
	if(fMultivariate){ np = 10000.0; sigma3 = 5000.;}
	else{ np = 1000.0; sigma3 = 1000.;}

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
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x,pt,type));}
		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum += fParameters[4] * TMath::Exp(-1*xprime*fParameters[4])*(ResolutionFunc(xprime-x,pt,type));}
	        }

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgNeg(Double_t x, Double_t pt, Int_t type) const 
{ 
	//
	// exponential with negative slopes for the background part (x)
	//
	Double_t np, sc, sigma3;
        sc = 10.;
        if(fMultivariate){ np = 10000.0;  sigma3 = 5000.;}
        else{ np = 1000.0; sigma3 = 1000.;}
	
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
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x,pt,type));}

		xprime = xupp - (i-.5) * step;
		if (xprime < 0) {sum += fParameters[5] * TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x,pt,type));}
	}

	return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgSym(Double_t x, Double_t pt, Int_t type) const 
{
	//
	// exponential with both positive and negative slopes for the background part (x)
	//
	Double_t np, sc, sigma3;
        sc = 10.;
        if(fMultivariate){ np = 10000.0; sigma3 = 5000.;}
        else{ np = 1000.0; sigma3 = 1000.;}

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
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x,pt,type));}
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x,pt,type));}

		xprime = xupp - (i-.5) * step;
		if (xprime > 0) {sum1 += 0.5 * fParameters[6]*TMath::Exp(-1*xprime*fParameters[6]) * (ResolutionFunc(xprime-x,pt,type));} 
		if (xprime < 0) {sum2 += 0.5 * fParameters[6]*TMath::Exp(xprime*fParameters[6]) * (ResolutionFunc(xprime-x,pt,type));}
	}

	return step*(sum1 + sum2) ;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBkgSym1(Double_t x, Double_t pt, Int_t type) const
{
	//
        // exponential with both positive and negative slopes for the background part (x)
        //
        Double_t np, sc, sigma3;
        sc = 10.;
        if(fMultivariate){ np = 10000.0; sigma3 = 5000.;}
        else{ np = 1000.0; sigma3 = 1000.;}

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
                if (xprime > 0) {sum1 += 0.5 * fParameters[45]*TMath::Exp(-1*xprime*fParameters[45]) * (ResolutionFunc(xprime-x,pt,type));}
                if (xprime < 0) {sum2 += 0.5 * fParameters[45]*TMath::Exp(xprime*fParameters[45]) * (ResolutionFunc(xprime-x,pt,type));}

                xprime = xupp - (i-.5) * step;
                if (xprime > 0) {sum1 += 0.5 * fParameters[45]*TMath::Exp(-1*xprime*fParameters[45]) * (ResolutionFunc(xprime-x,pt,type));} 
                if (xprime < 0) {sum2 += 0.5 * fParameters[45]*TMath::Exp(xprime*fParameters[45]) * (ResolutionFunc(xprime-x,pt,type));}
        }

        return step*(sum1 + sum2) ;
}


//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFunc(Double_t x, Double_t pt, Int_t type) const  
{
	//
	// parametrization with 2 gaus
	//
        x = x/fChangeResolution;
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
 
        if(fMultivariate) // set parameters from matrix
        {
          //pt dependence
          Int_t binPt = -1.;
          for(int j=0; j<fPtWindows->GetSize()-1; j++) {if(fPtWindows->At(j)<pt && pt<fPtWindows->At(j+1)) binPt = j;}
          norm1 = fResParams[binPt][type][0];
          mean1 = fResParams[binPt][type][1];
          sigma1 = fResParams[binPt][type][2];
          norm2 = fResParams[binPt][type][3];
          mean2 = fResParams[binPt][type][4];
          sigma2 = fResParams[binPt][type][5];
          alfa = fResParams[binPt][type][6];
          lambda = fResParams[binPt][type][7];
          norm3 = fResParams[binPt][type][8];
        }

        Double_t ret = 0.; Double_t fitval = 0; 
        if(TMath::Abs(x)<=alfa) fitval = (lambda-1)/(2*alfa*lambda);
        else  fitval = ((lambda-1)/(2*alfa*lambda))*TMath::Power(alfa,lambda)*(TMath::Power(TMath::Abs(x),-1*lambda));

        ret = (norm1/(norm1+norm2+norm3))*((1/(sigma1*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean1)/sigma1)*((x-mean1)/sigma1))) + (norm2/(norm1+norm2+norm3))*((1/(sigma2*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x-mean2)/sigma2)*((x-mean2)/sigma2))) + (norm3/(norm1+norm2+norm3))*fitval;

        return ret/fChangeResolution;
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
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetResolutionFunc(Double_t xmin, Double_t xmax, Double_t normalization, Double_t pt, Int_t type){
	// return the pointer to the resolution function
	TF1* resFunc = new TF1(Form("resolutionFunc_%d",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::ResolutionFuncf,xmin,xmax,3);
        resFunc->SetParameter(0,normalization);
        resFunc->SetParameter(1,pt);
        resFunc->SetParameter(2,(Double_t)type);
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
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeBkgDistr(Double_t xmin, Double_t xmax, Double_t normalization, Int_t type, Double_t mass, Double_t pt, Int_t npx){
	// return the pointer to the background x distribution function
	TF1 *backFunc = new TF1(Form("backFunc_%d_%1.2f_%1.2f",type,mass,pt),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrFunc,xmin,xmax,4);
        backFunc->SetParameter(0,normalization);
        backFunc->SetParameter(1,(Double_t)type);
        backFunc->SetParameter(2,mass);
        backFunc->SetParameter(3,pt);
	backFunc->SetNpx(npx);
        return (TF1*)backFunc->Clone();
}

//___________________________________________________________________________________________________
TF1* AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeBkgDistrAllTypes(Double_t xmin, Double_t xmax, Double_t normalization){
        // return the pointer to the background x distribution function
        TF1 *backFuncNew = new TF1(Form("backFunc_%f",normalization),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrFuncAllTypes,xmin,xmax,1);
        backFuncNew->SetParameter(0,normalization);
        backFuncNew->SetNpx(5000);
        return (TF1*)backFuncNew->Clone();
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
TF1 *AliDielectronBtoJPSItoEleCDFfitFCN::GetEvaluateCDFDecayTimeTotalDistr(Double_t xMin, Double_t xMax, Double_t normalization,Double_t pt, Int_t type){
 // return the pointer to the pseudoproper distribution for the background
 TF1 *decayTimeTot = new TF1(Form("decayTimeTot_%d",type),this,&AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistr,xMin,xMax,3);
 decayTimeTot->SetParameter(0,normalization);
 decayTimeTot->SetParameter(1,pt);
 decayTimeTot->SetParameter(2,(Double_t)type);
 decayTimeTot->SetNpx(5000);
 return (TF1*)decayTimeTot->Clone();
}

//____________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistr(const Double_t* x, const Double_t *par) const
{
 // evaluate the total pseudoproper distribution for a given candidate's type. par[1] should be the candidate's type.
 Double_t value = 0;
 Double_t xx = x[0];
 value = (fParameters[8]*EvaluateCDFDecayTimeSigDistr(xx,par[1],(Int_t)par[2]) + (1-fParameters[8])*EvaluateCDFDecayTimeBkgDistr(xx,(Int_t)par[2],par[1]))*par[0];
 return value;
}

//____________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeTotalDistrAllTypes(const Double_t* x, const Double_t *par) const
{
 // evaluate the total pseudoproper distribution considering all candidate's types
 Double_t value = 0;
 Double_t xx = x[0];

 value = (fParameters[8]*(fWeightType[2]*EvaluateCDFDecayTimeSigDistr(xx,200.,2)+fWeightType[1]*EvaluateCDFDecayTimeSigDistr(xx,200.,1)+fWeightType[0]*EvaluateCDFDecayTimeSigDistr(xx,200.,0)))+((1-fParameters[8])*(fWeightType[2]*EvaluateCDFDecayTimeBkgDistr(xx,2,3.09,200.) + fWeightType[1]*EvaluateCDFDecayTimeBkgDistr(xx,1,3.09,200.)+fWeightType[0]*EvaluateCDFDecayTimeBkgDistr(xx,0,3.09,200.))); 

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
TF1 * AliDielectronBtoJPSItoEleCDFfitFCN::GetFunB(Double_t xmin, Double_t xmax, Double_t normalization, Double_t pt, Int_t type, Int_t npx){
 // return the pointer to the function that describe secondary jpsi pseudoproper distribution for a given candidate's type
 TF1* funb = new TF1(Form("secondaryJpsiConvolution_%d_%1.2f",type,pt),this,&AliDielectronBtoJPSItoEleCDFfitFCN::FunBfunc,xmin,xmax,3);
 funb->SetParameter(0,normalization);
 funb->SetParameter(1,pt);
 funb->SetParameter(2,(Double_t)type);
 funb->SetNpx(npx);
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

//
// methods below are needed to perform the multivariate fit (pt,mass,type); this can be enabled 
// by the boolean fMultivariate; if functions to describe pseudoproper
// decay lenght in pt,m, type bins have been
// already computed, they can be loaded from the file 
// switching-on the option fLoadFunction (this to avoid the
// calculation of convolutions and speed-up the likelihood fit)
//

//____________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetBackgroundSpecificParameters(Int_t pt, Int_t mb, Int_t tp){
  //
  // methods to set specific background parameters in bins(pt, inv. mass, type)
  //
  for(int j=0; j<4;j++) fParameters[j]=fBkgParams[pt][mb][tp][j];
  for(int k=5; k<8;k++) fParameters[k-1]=fBkgParams[pt][mb][tp][k];
  fParameters[45] = fBkgParams[pt][mb][tp][8];
  fParameters[46] = fBkgParams[pt][mb][tp][4];
  return;
}

//_______________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::InitializeFunctions(Int_t ptSize, Int_t massSize){
  //
  // initialize pointers to save functions for the  multivariate fit
  //
  fFunBSaved = new TF1**[ptSize];
  for(int kpt=0; kpt<ptSize; kpt++) fFunBSaved[kpt] = new TF1*[3]; // type
  fFunBkgSaved = new TF1***[ptSize];
  for(int kpt=0; kpt<ptSize; kpt++){ fFunBkgSaved[kpt] = new TF1**[massSize];
  for(int ks = 0; ks<massSize; ks++) fFunBkgSaved[kpt][ks] = new TF1*[3];
  for(int kk=0; kk<3;kk++) {
     fFunBSaved[kpt][kk] = new TF1();
     for(int kk1=0; kk1<massSize;kk1++){ 
     fFunBkgSaved[kpt][kk1][kk] = new TF1();
     }
    }
  }

  // to extrapolate the function under the signal region
  fWeights = new Double_t**[massSize - 1]; // mass
  for(int km =0; km < (massSize - 1); km++) {fWeights[km] = new Double_t*[ptSize];
        for(int kpt =0; kpt<ptSize; kpt++) fWeights[km][kpt] = new Double_t[3];
  } // pt
 return;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::FunBsaved(Double_t x, Double_t pt, Int_t type) const
{
        //
        //   functions to describe non-prompt J/psi x distributions
        //
        Double_t returnvalue = 0.;
        Int_t binPt = -1;
        for(int j=0; j<fPtWindows->GetSize()-1; j++) {if(fPtWindows->At(j)<pt && pt<fPtWindows->At(j+1)) binPt = j;}
        returnvalue = fFunBSaved[binPt][type]->Eval(x);
        return returnvalue;
}

//________________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCN::SetFunctionsSaved(Int_t npxFunB, Int_t npxFunBkg, Double_t funBLimits, Double_t funBkgLimits, Int_t signalRegion){
 	//
 	// save functions for the multivariate fit 
 	//
 	if(!fMultivariate)
 	{AliInfo("Warning: fMultivariate is kFALSE! Functions are not saved! \n"); return;}
 	SetExtrapolationRegion(signalRegion);

                     for(int tp=0;tp<3;tp++)  // type
                      {
                        // pt
                        for(int pt=0; pt<fPtWindows->GetSize()-1;pt++){
                        if(fResParams) SetResolutionConstants(fResParams[pt][tp],tp);
                        SetFunBFunction(tp,pt,GetFunB(-1.*funBLimits,funBLimits,1.,(fPtWindows->At(pt) + (fPtWindows->At(pt+1)-fPtWindows->At(pt))/2.),tp,npxFunB));
                        }
                      }

 	AliInfo("+++++++  Pseudoproper-decay-length function for secondary J/psi saved  ++++++ \n");

 	if(!fLoadFunctions){
 	for(int ij = 0; ij<fMassWindows->GetSize()-1;ij++){
 	if(ij == signalRegion) continue;

    	Int_t mbin = (ij > signalRegion) ?  ij-1 : ij;
    	for(int tp=0;tp<3;tp++)  {
      	for(int pt =0; pt<fPtWindows->GetSize()-1; pt++){
         if(fBkgParams) SetBackgroundSpecificParameters(pt,mbin,tp);
         SetBkgFunction(ij, tp, pt, GetEvaluateCDFDecayTimeBkgDistr(-1.*funBkgLimits,funBkgLimits,1.,tp,(fMassWindows->At(ij) + (fMassWindows->At(ij+1)-fMassWindows->At(ij))/2.),(fPtWindows->At(pt) + (fPtWindows->At(pt+1)-fPtWindows->At(pt))/2.),npxFunBkg));

        	}

	}

	}
	AliInfo("+++++++  Pseudoproper-decay-length function for background saved  +++++++++++ \n");

 	} // loadFunctions
 	// evaluate under signal
 	for(int tp=0;tp<3;tp++)
       		{
        	for(int pt =0; pt<fPtWindows->GetSize()-1; pt++){
        	SetBkgFunction(signalRegion, tp, pt, GetEvaluateCDFDecayTimeBkgDistr(-1.*funBkgLimits,funBkgLimits,1.,tp,(fMassWindows->At(signalRegion) + (fMassWindows->At(signalRegion+1)-fMassWindows->At(signalRegion))/2.),(fPtWindows->At(pt) + (fPtWindows->At(pt+1)-fPtWindows->At(pt))/2.),npxFunBkg));
        	}
   
  	}
 	// save functions
 	TFile func("functions.root","RECREATE");
 	for(int kpt =0; kpt<fPtWindows->GetSize()-1; kpt++){
 	for(int ss=0; ss<3;ss++) {fFunBSaved[kpt][ss]->Write();
 	for(int kk=0; kk<fMassWindows->GetSize()-1; kk++) fFunBkgSaved[kpt][kk][ss]->Write();}}
 	return;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrDifferential(Double_t x, Int_t type, Double_t m, Double_t pt) const
{
        //
        // it returns the value of the probability to have a given x for the background 
        // in the pt, m , type correspondent range 
        //
        Int_t binPt = -1;
        for(int j=0; j<fPtWindows->GetSize()-1; j++)
        {if(fPtWindows->At(j)<pt && pt<fPtWindows->At(j+1)) binPt = j;}
        Bool_t isSignal = (fMassWindows->At(fSignalBinForExtrapolation)<m && m<fMassWindows->At(fSignalBinForExtrapolation+1));
        Double_t ret = 0.;
        if(!isSignal)
        ret = fParameters[0]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*ResolutionFunc(x, pt, type) + fParameters[1]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgPos(x, pt,type) +  fParameters[2]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgNeg(x,pt,type) + fParameters[3]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgSym(x, pt,type) + fParameters[46]/(fParameters[0]+fParameters[1]+fParameters[2]+fParameters[3]+fParameters[46])*FunBkgSym1(x,pt,type);
       else{
       for(int k=0; k<fMassWindows->GetSize()-2;k++) {
        Int_t mbin = (k > (fSignalBinForExtrapolation-1)) ?  k+1 : k;
        ret +=  fWeights[k][binPt][type]*EvaluateCDFDecayTimeBkgDistrSaved(x,type,(fMassWindows->At(mbin) + (fMassWindows->At(mbin+1)-fMassWindows->At(mbin))/2.),pt);}

       }
        return ret;
}

//_________________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistrSaved(Double_t x, Int_t type, Double_t m, Double_t pt) const
{
        //
        // it returns the value of the probability to have a given x for the background 
        // in the pt, m , type correspondent range 
        //
        Double_t returnvalue = 0.;
        Int_t binM = -1.;
        for(int j=0; j<fMassWindows->GetSize()-1; j++) {if(fMassWindows->At(j)<m && m<fMassWindows->At(j+1)) binM = j;}
        Int_t binPt = -1;
        for(int j=0; j<fPtWindows->GetSize()-1; j++) {if(fPtWindows->At(j)<pt && pt<fPtWindows->At(j+1)) binPt = j;}
        returnvalue = fFunBkgSaved[binPt][binM][type]->Eval(x);
        return returnvalue;
}

