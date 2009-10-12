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
#include "AliBtoJPSItoEleCDFfitFCN.h"

//_________________________________________________________________________
//                        Class AliBtoJPSItoEleCDFfitFCN
//                   Definition of main function used in 
//                     unbinned log-likelihood fit for
//                 the channel B -> JPsi + X -> e+e- + X
//      
//                           Origin: C.Di Giglio
//       Contact: Carmelo.Digiglio@ba.infn.it , Giuseppe.Bruno@ba.infn.it
//_________________________________________________________________________

ClassImp(AliBtoJPSItoEleCDFfitFCN)

//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN::AliBtoJPSItoEleCDFfitFCN() :
fFPlus(0.),
fFMinus(0.),
fFSym(0.),
fIntegral(1.),
fintxFunB(1.),
fintxDecayTimeBkgPos(1.),
fintxDecayTimeBkgNeg(1.),
fintxDecayTimeBkgSym(1.),
fintmMassSig(1.),
fintxRes(1.),
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
  for(Int_t iPar = 0; iPar < 16; iPar++) fParameters[iPar] = 0.;
  fParameters[9] = 1.;fParameters[11] = 1.;fParameters[12] = 1.;
  for(Int_t index=0; index<6; index++) fResolutionConstants[index] = 0.;
  AliInfo("Instance of AliBtoJPSItoEleCDFfitFCN-class created");
}
//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN::AliBtoJPSItoEleCDFfitFCN(const AliBtoJPSItoEleCDFfitFCN& source) :
TNamed(source),
fFPlus(source.fFPlus),
fFMinus(source.fFMinus),
fFSym(source.fFSym),
fIntegral(source.fIntegral),
fintxFunB(source.fintxFunB),
fintxDecayTimeBkgPos(source.fintxDecayTimeBkgPos),
fintxDecayTimeBkgNeg(source.fintxDecayTimeBkgNeg),
fintxDecayTimeBkgSym(source.fintxDecayTimeBkgSym),
fintmMassSig(source.fintmMassSig),
fintxRes(source.fintxRes),
fintmMassBkg(source.fintmMassBkg),
fhCsiMC(source.fhCsiMC),
fMassWndHigh(source.fMassWndHigh),
fMassWndLow(source.fMassWndLow),
fCrystalBallParam(source.fCrystalBallParam)
{
  //
  // Copy constructor
  //
  for(Int_t iPar = 0; iPar < 16; iPar++) fParameters[iPar] = source.fParameters[iPar];
  for(Int_t index=0; index<6; index++) fResolutionConstants[index] = source.fResolutionConstants[index];
}
//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN& AliBtoJPSItoEleCDFfitFCN::operator=(const AliBtoJPSItoEleCDFfitFCN& source) 
{
  //
  // Assignment operator
  //
  if(&source == this) return *this;
  fFPlus = source.fFPlus;
  fFMinus = source.fFMinus;
  fFSym = source.fFSym;
  fIntegral = source.fIntegral;
  fintxFunB = source.fintxFunB;
  fintxDecayTimeBkgPos = source.fintxDecayTimeBkgPos;
  fintxDecayTimeBkgNeg = source.fintxDecayTimeBkgNeg;
  fintxDecayTimeBkgSym = source.fintxDecayTimeBkgSym;
  fintmMassSig = source.fintmMassSig;
  fintxRes = source.fintxRes;
  fintmMassBkg = source.fintmMassBkg;
  fhCsiMC = source.fhCsiMC;
  fCrystalBallParam = source.fCrystalBallParam;

  for(Int_t iPar = 0; iPar < 16; iPar++) fParameters[iPar] = source.fParameters[iPar];
  for(Int_t index=0; index<6; index++) fResolutionConstants[index] = source.fResolutionConstants[index];

  return *this;
}  
//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN::~AliBtoJPSItoEleCDFfitFCN()
{
  //
  // Default destructor
  //
  
  delete fhCsiMC;
  for(Int_t iPar = 0; iPar < 16; iPar++) fParameters[iPar] = 0.;
  for(Int_t index=0; index<6; index++) fResolutionConstants[index] = 0.;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
           const Double_t* invariantmass, const Int_t ncand)
{
  //
  // This function evaluates the Likelihood fnction
  // It returns the -Log(of the likelihood function)
  //
  Double_t f = 0.;
  Double_t ret = 0.;

  for(Int_t i=0; i < ncand; i++) {
      f = EvaluateCDFfuncNorm(pseudoproperdecaytime[i],invariantmass[i]);
      if(f < 0.) {
        continue;  
      }
      ret+=-1.*TMath::Log(f);
      }
  return ret;
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::SetAllParameters(const Double_t* parameters)
{ 
  //
  // Sets array of FCN parameters
  //
  for(Int_t index = 0; index < 16; index++) fParameters[index] = parameters[index];
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::ComputeIntegral() 
{ 
  //
  // this function compute the integral of the likelihood function 
  // (theoretical function) in order to normalize it to unity
  //
  Double_t np = 100.0;                 //number of integration steps 
  Double_t npres = 200.0;              //number of integration steps for the resolution function
  Double_t npm = 200.;
  Double_t stepm;Double_t stepx;Double_t stepxres;       //integration step width in variable m,x 
  Double_t mx=0.;Double_t xprime=0.;
  Double_t xlow = -4000.; Double_t xup = 4000.;
  Double_t i; Double_t j;
  Double_t sumx = 0.0;Double_t intx = 0.0;Double_t intm = 0.0;
  stepm = (fMassWndHigh-fMassWndLow)/npm; 
  stepx = (xup-xlow)/np;
  stepxres = (xup-xlow)/npres;

// compute integrals for all the terms        

  Double_t iRes;
  Double_t intxRes = 0.0;
  Double_t sumxRes = 0.0;
       for(iRes = 1.0; iRes <= npres/2.; iRes++)  {
           xprime = xlow + (iRes - .5)*stepxres;
           sumxRes += ResolutionFunc(xprime);
           xprime = xup - (iRes - .5)*stepxres;
           sumxRes += ResolutionFunc(xprime);
           }
       intxRes = sumxRes*stepxres;
       SetIntegralRes(intxRes);

//
  Double_t iFunB;
  Double_t intxFunB = 0.0;
  Double_t sumxFunB = 0.0;
       for(iFunB = 1.0; iFunB <= np/2; iFunB++)  {
           xprime = xlow + (iFunB - .5)*stepx;
           sumxFunB += FunB(xprime);
           xprime = xup - (iFunB - .5)*stepx;
           sumxFunB += FunB(xprime);
           }
       intxFunB = sumxFunB*stepx;
       SetIntegralFunB(intxFunB);

//
  Double_t iDecayTimeBkgPos;
  Double_t intxDecayTimeBkgPos = 0.0;
  Double_t sumxDecayTimeBkgPos = 0.0;
       for(iDecayTimeBkgPos = 1.0; iDecayTimeBkgPos <= np/2; iDecayTimeBkgPos++)  {
           xprime = xlow + (iDecayTimeBkgPos - .5)*stepx;
           sumxDecayTimeBkgPos += FunBkgPos(xprime);
           xprime = xup - (iDecayTimeBkgPos - .5)*stepx;
           sumxDecayTimeBkgPos += FunBkgPos(xprime);
           }
       intxDecayTimeBkgPos = sumxDecayTimeBkgPos*stepx;
       SetIntegralBkgPos(intxDecayTimeBkgPos);

//
  Double_t iDecayTimeBkgNeg;
  Double_t intxDecayTimeBkgNeg = 0.0;
  Double_t sumxDecayTimeBkgNeg = 0.0;
       for(iDecayTimeBkgNeg = 1.0;  iDecayTimeBkgNeg<= np/2; iDecayTimeBkgNeg++)  {
           xprime = xlow + (iDecayTimeBkgNeg - .5)*stepx;
           sumxDecayTimeBkgNeg += FunBkgNeg(xprime);
           xprime = xup - (iDecayTimeBkgNeg - .5)*stepx;
           sumxDecayTimeBkgNeg += FunBkgNeg(xprime);
           }
       intxDecayTimeBkgNeg = sumxDecayTimeBkgNeg*stepx;
       SetIntegralBkgNeg(intxDecayTimeBkgNeg);
//
  Double_t iDecayTimeBkgSym;
  Double_t intxDecayTimeBkgSym = 0.0;
  Double_t sumxDecayTimeBkgSym = 0.0;
       for(iDecayTimeBkgSym = 1.0; intxDecayTimeBkgSym <= np/2; intxDecayTimeBkgSym++)  {
           xprime = xlow + (intxDecayTimeBkgSym - .5)*stepx;
           sumxDecayTimeBkgSym += FunBkgSym(xprime);
           xprime = xup - (intxDecayTimeBkgSym - .5)*stepx;
           sumxDecayTimeBkgSym += FunBkgSym(xprime);
           }
       intxDecayTimeBkgSym = sumxDecayTimeBkgSym*stepx;
       SetIntegralBkgSym(intxDecayTimeBkgSym);

//
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

//
// Compute integral of the whole distribution function
//
       for(i = 1.0; i <= np; i++)  {
           Double_t summ = 0.0;
           xprime = xlow + (i - .5)*stepx;
         for(j = 1.0; j <= npm/2; j++)  { 
             mx = fMassWndLow + (j - .5)*stepm;
             summ += EvaluateCDFfunc(xprime,mx);
             mx = fMassWndHigh - (j - .5)*stepm;
             summ += EvaluateCDFfunc(xprime,mx);
             }
           intm = summ*stepm; 
           sumx += intm; 
           }
       intx = sumx*stepx;
       SetIntegral(intx);
  
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::PrintStatus()
{
  //
  //  Print the parameters of the fits 
  //
  printf("\n");
  printf("actual value of fRadius---------------------------------------->> | %f \n", GetRadius());
  printf("actual value of fTheta ---------------------------------------->> | %f \n", GetTheta());
  printf("actual value of fPhi ------------------------------------------>> | %f \n", GetPhi());
  printf("actual value of fFPlus ---------------------------------------->> | %f \n", GetFPlus());
  printf("actual value of fFMinus --------------------------------------->> | %f \n", GetFMinus());
  printf("actual value of fFSym ----------------------------------------->> | %f \n", GetFSym());
  printf("actual value of fOneOvLamPlus --------------------------------->> | %f \n", GetLamPlus());
  printf("actual value of fOneOvLamMinus -------------------------------->> | %f \n", GetLamMinus());
  printf("actual value of fOneOvLamSym ---------------------------------->> | %f \n", GetLamSym());
  printf("actual value of fMassBkgSlope --------------------------------->> | %f \n", GetMassSlope());
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
  printf("actual value of fSigmaResol ----------------------------------->> | %f \n", GetSigmaResol());
  printf("actual value of fNResol --------------------------------------->> | %f \n", GetNResol());
  printf("\n");
  printf("Actual value of normalization integral for FunB ------------------->> | %f \n", GetIntegralFunB());
  printf("Actual value of normalization integral for BkgPos ----------------->> | %f \n", GetIntegralBkgPos());
  printf("Actual value of normalization integral for BkgNeg ----------------->> | %f \n", GetIntegralBkgNeg());
  printf("Actual value of normalization integral for BkgSym ----------------->> | %f \n", GetIntegralBkgSym());
  printf("Actual value of normalization integral for MassSig ---------------->> | %f \n", GetIntegralMassSig());
  printf("Actual value of normalization integral for MassBkg ---------------->> | %f \n", GetIntegralMassBkg());
  printf("Actual value of normalization integral for Resolution ------------->> | %f \n", GetIntegralRes());
  printf("Actual value of normalization integral for FCN -------------------->> | %f \n", GetIntegral());

  printf("\n");
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::SetResolutionConstants()
{
  //
  //  This method must be update: 
  //  for the time beeing the values are hard-wired. 
  //  Implementations have to be done to set the values from outside 
  //  (e.g. from a ConfigHF file) starting from an indipendent fit 
  //  of primary JPSI distribution.
  //

  fResolutionConstants[0]  = 8.; // mean sigma2/sigma1 
  fResolutionConstants[1]  = 0.1675; // mean Integral2/Integral1
  fResolutionConstants[2]   = 1374.; // sigma2
  fResolutionConstants[3]  = 0.001022; // N2
  fResolutionConstants[4]  = 686.6; // mu2
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFfunc(Double_t x, Double_t m) const 
{
  return fParameters[8]*EvaluateCDFfuncSignalPart(x,m) + (1. - fParameters[8])*EvaluateCDFfuncBkgPart(x,m);
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncNorm(Double_t x, Double_t m) const
{
  return EvaluateCDFfunc(x,m)/fIntegral;
} 
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncSignalPart(Double_t x, Double_t m) const 
{
  return EvaluateCDFDecayTimeSigDistr(x)*(EvaluateCDFInvMassSigDistr(m)/fintmMassSig); 
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistr(Double_t x) const
{
  //
  // Implementation of the Background part of the Likelyhood function
  // 

  Double_t retvalue = 0.;
  Double_t FunBnorm = FunB(x)/fintxFunB;
  Double_t FunPnorm = ResolutionFunc(x)/fintxRes;
  retvalue = fParameters[7]*FunBnorm + (1. - fParameters[7])*FunPnorm;
  return retvalue;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassSigDistr(Double_t m) const
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
Double_t AliBtoJPSItoEleCDFfitFCN::FunB(Double_t x) const 
{
  //  
  // Parameterisation of the fit function for the x part of the Background
  //

  Double_t np = 100.0;
  Double_t sc = 10.;
  Double_t sigma3 =  fResolutionConstants[2];
  Double_t xprime;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  xlow = x - sc * sigma3 ;
  xupp = x + sc * sigma3 ;
  step = (xupp-xlow) / np;
  Double_t CsiMCxprime = 0.;
  Double_t Resolutionxdiff = 0.;
  Double_t xdiff = 0.;

  for(i=1.0; i<=np; i++) {

      xprime = xlow + (i-.5) * step;
      CsiMCxprime = CsiMC(xprime);
      xdiff = xprime - x;
      Resolutionxdiff = ResolutionFunc(xdiff)/fintxRes; // normalized value
      sum += CsiMCxprime * Resolutionxdiff;

      }

  return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunP(Double_t x) const 
{
  //
  //  Parameterisation of the Prompt part for the x distribution
  //

  return ResolutionFunc(x)/fintxRes; // normalized value
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::CsiMC(Double_t x) const 
{
  //
  //  Distribution (template) of the x distribution for the x variable 
  //  for the J/psi coming from Beauty hadrons
  //

  Double_t returnvalue = 0.; 
  returnvalue = fhCsiMC->GetBinContent(fhCsiMC->FindBin(x));
  return returnvalue;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFfuncBkgPart(Double_t x,Double_t m) const 
{
  //
  // Return the part of the likelihood function for the background hypothesis
  //

  return EvaluateCDFDecayTimeBkgDistr(x)*(EvaluateCDFInvMassBkgDistr(m)/fintmMassBkg);
}  
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistr(Double_t x) const 
{
  //
  // it returns the value of the probability to have a given x for the background 
  //
 
  Double_t ret = (1. - TMath::Power(fParameters[0],2.))*(ResolutionFunc(x)/fintxRes)
                     + TMath::Power(fParameters[0]*TMath::Cos(fParameters[1]),2.)*
                                              (FunBkgPos(x)/fintxDecayTimeBkgPos)
                     + TMath::Power(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2]),2.)*
                                              (FunBkgNeg(x)/fintxDecayTimeBkgNeg)
                     + TMath::Power(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]),2.)*
                                              (FunBkgSym(x)/fintxDecayTimeBkgSym);
  return ret;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassBkgDistr(Double_t m) const 
{
  //
  // it returns the value of the probability to have a given mass for the background
  //

  return 1/(fMassWndHigh-fMassWndLow) + 
         fParameters[6] * m - 
         fParameters[6] * ((fMassWndHigh+fMassWndLow)/2);
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgPos(Double_t x) const 
{
  //
  // exponential with positive slopes for the background part (x)
  //

  Double_t np = 100.0;
  Double_t sc = 10.;      
  Double_t sigma3 = fResolutionConstants[2];
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
       if (xprime > 0) {sum += fParameters[3] * TMath::Exp(-1*xprime*fParameters[3]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum = 0.;}
 
      xprime = xupp - (i-.5) * step;
      if (xprime > 0) {sum += fParameters[3] * TMath::Exp(-1*xprime*fParameters[3]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum = 0.;}

      }

  return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgNeg(Double_t x) const 
{
  //
  // exponential with negative slopes for the background part (x)
  //

  Double_t np = 100.0;
  Double_t sc = 10.;      
  Double_t sigma3 =  fResolutionConstants[2];
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
      if (xprime < 0) {sum += fParameters[4] * TMath::Exp(xprime*fParameters[4]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum = 0.;}

      xprime = xupp - (i-.5) * step;
      if (xprime < 0) {sum += fParameters[4] * TMath::Exp(xprime*fParameters[4]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum = 0.;}
      }

  return step * sum ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgSym(Double_t x) const 
{
  //
  // exponential with both positive and negative slopes for the background part (x)
  //

  Double_t np = 100.0;
  Double_t sc = 10.;      
  Double_t sigma3 =  fResolutionConstants[2];
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
      if (xprime > 0) {sum1 += 0.5 * fParameters[5]*TMath::Exp(-1*xprime*fParameters[5]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum1 = 0.;}

      xprime = xupp - (i-.5) * step;
      if (xprime > 0) {sum1 += 0.5 * fParameters[5]*TMath::Exp(-1*xprime*fParameters[5]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum1 = 0.;}
      }

  for(i=1.0; i<=np/2; i++) {

      xprime = xlow + (i-.5) * step;
      if (xprime < 0) {sum2 += 0.5 * fParameters[5]*TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum2 = 0.;}

      xprime = xupp - (i-.5) * step;
      if (xprime < 0) {sum2 += 0.5 * fParameters[5]*TMath::Exp(xprime*fParameters[5]) * (ResolutionFunc(xprime-x)/fintxRes);} else {sum2 = 0.;}
      }

  return step*(sum1 + sum2) ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::ResolutionFunc(Double_t x) const 
{
  //
  // parametrization with 2 gaus
  //

  Double_t ret = 0.;
  Double_t x1 = x;
  Double_t x2 = x;
  //Double_t mean1 = 0.; 
  Double_t mean2 = fResolutionConstants[4];
  Double_t sigma1 = fParameters[14]; 
  Double_t sigma2 = fResolutionConstants[2];
  Double_t n1 = fParameters[15]; 
  Double_t n2 = fResolutionConstants[3];
  Double_t arg1 = x1/sigma1;
  Double_t arg2 = (x2-mean2)/sigma2;
  Double_t sqrt2Pi = TMath::Sqrt(2*TMath::Pi());

  ret = n2*((n1/n2)*TMath::Exp(-0.5*arg1*arg1) + TMath::Exp(-0.5*arg2*arg2));

  return ret/(sqrt2Pi*(n1*sigma1+n2*sigma2));//return value is normalized

}

