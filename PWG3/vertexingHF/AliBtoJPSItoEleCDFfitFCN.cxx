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
fIntegral(0.),
fhCsiMC(0x0),
fMassWndHigh(0.),
fMassWndLow(0.),
fCrystalBallParam(kFALSE)
{
  //
  // constructor
  //
  SetCrystalBallParam(kFALSE);
  SetMassWndHigh(0.2);
  SetMassWndLow(0.5);
  for(Int_t iPar = 0; iPar < 13; iPar++) fParameters[iPar] = 0.;
  fParameters[9] = 1.;fParameters[11] = 1.;fParameters[12] = 1.;
  for(Int_t index=0; index<7; index++) fResolutionConstants[index] = 0.;
  AliInfo("Instance of AliBtoJPSItoEleCDFfitFCN-class created");
}
//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN::AliBtoJPSItoEleCDFfitFCN(const AliBtoJPSItoEleCDFfitFCN& source) :
TNamed(source),
fFPlus(source.fFPlus),
fFMinus(source.fFMinus),
fFSym(source.fFSym),
fIntegral(source.fIntegral),
fhCsiMC(source.fhCsiMC),
fMassWndHigh(source.fMassWndHigh),
fMassWndLow(source.fMassWndLow),
fCrystalBallParam(source.fCrystalBallParam)
{
  //
  // Copy constructor
  //
  for(Int_t iPar = 0; iPar < 13; iPar++) fParameters[iPar] = source.fParameters[iPar];
  for(Int_t index=0; index<7; index++) fResolutionConstants[index] = source.fResolutionConstants[index];
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
  fhCsiMC = source.fhCsiMC;
  fCrystalBallParam = source.fCrystalBallParam;

  for(Int_t iPar = 0; iPar < 13; iPar++) fParameters[iPar] = source.fParameters[iPar];
  for(Int_t index=0; index<7; index++) fResolutionConstants[index] = source.fResolutionConstants[index];

  return *this;
}  
//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitFCN::~AliBtoJPSItoEleCDFfitFCN()
{
  //
 // Default destructor
  //
  
  delete fhCsiMC;
  for(Int_t iPar = 0; iPar < 13; iPar++) fParameters[iPar] = 0.;
  for(Int_t index=0; index<7; index++) fResolutionConstants[index] = 0.;
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
        //AliWarning("One negative contributors in the Log(Likely) ! ");
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
  for(Int_t index = 0; index < 13; index++) fParameters[index] = parameters[index];
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::ComputeIntegral() 
{ 
//
// this function compute the integral of the likelihood function 
// (theoretical function) in order to normalize it to unity
//
  Double_t np = 50.0;                  //number of integration steps 
  Double_t stepm;Double_t stepx;       //integration step width in variable m,x 
  Double_t mx;Double_t xx;
  Double_t xlow = -4000.; Double_t xup = 4000.;
  Double_t i; Double_t j;
  Double_t sumx = 0.0;Double_t intx = 0.0;Double_t intm = 0.0;
  stepm = (fMassWndHigh-fMassWndLow)/np; 
  stepx = (xup-xlow)/np;
                   
        for(i = 1.0; i <= np; i++)  {
            Double_t summ = 0.0;
            xx = xlow + (i - .5)*stepx;
          for(j = 1.0; j <= np/2; j++)  { 
              mx = fMassWndLow + (j - .5)*stepm;
              summ += EvaluateCDFfunc(xx,mx);
              mx = fMassWndHigh - (j - .5)*stepm;
              summ += EvaluateCDFfunc(xx,mx);
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
  }else{
   printf("actual value of fMpv ------------------------------------------>> | %f \n", GetCrystalBallMmean());
   printf("actual value of fConstRovL ------------------------------------>> | %f \n", GetCrystalBallNexp());
   printf("actual value of fSigmaL --------------------------------------->> | %f \n", GetCrystalBallSigma());
   printf("actual value of fSigmaR --------------------------------------->> | %f \n", GetCrystalBallAlpha());
  }
  printf("\n");
  printf("Actual value of normalization integral for FCN ---------------->> | %f \n", GetIntegral());
  printf("\n");
}
//_________________________________________________________________________________________________
void AliBtoJPSItoEleCDFfitFCN::SetResolutionConstants(Int_t BinNum)
{
//  This method must be update. 
//  For the time beeing the values are hard-wired. 
//  Implementations have to be done to set the values from outside (e.g. from a ConfigHF file)
//
  switch(BinNum){

   case(kallpt):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*535.9 ; fResolutionConstants[4] = 1.171*5535.9  ;//from fit integrated in pt
    fResolutionConstants[1]  = 0.0998*535.9; fResolutionConstants[3] = 0.1072   ; fResolutionConstants[5] = 0.04115    ; fResolutionConstants[6] = 1e-04;
    break;
   case(kptbin1):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*1087  ; fResolutionConstants[4] = 1.171*1087    ;//from fit integrated in pt
    fResolutionConstants[1]  = 0.04253*1087 ; fResolutionConstants[3] = 0.1482  ; fResolutionConstants[5] = 0.09778    ; fResolutionConstants[6] = 3.773e-04;
    break;
   case(kptbin2):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*661.5 ; fResolutionConstants[4] = 1.171*661.5   ;//from fit integrated in pt
    fResolutionConstants[1]  = 0.1*661.5    ; fResolutionConstants[3] = 0.2809  ; fResolutionConstants[5] =  0.09771   ; fResolutionConstants[6] = 1.916e-04;
    break;
   case(kptbin3):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1]  = 0.1578*502.8 ; fResolutionConstants[3] = 0.3547  ; fResolutionConstants[5] = 0.09896    ; fResolutionConstants[6] = 5.241e-04;
    break;
   case(kptbin4):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] = 0.2048*415.9  ; fResolutionConstants[3] = 0.4265  ; fResolutionConstants[5] = 0.09597    ; fResolutionConstants[6] = 6.469e-04;
    break;
   case(kptbin5):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] = 0.2219*379.7  ; fResolutionConstants[3] = 0.5414 ;  fResolutionConstants[5] = 0.07506    ; fResolutionConstants[6] = 7.465e-04;
    break;
   case(kptbin6):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] = 0.2481*307.   ; fResolutionConstants[3] = 0.8073 ;  fResolutionConstants[5] = 0.09664    ; fResolutionConstants[6] = 8.481e-04;
    break;
   case(kptbin7):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] = 0.262*283.5   ; fResolutionConstants[3] = 0.9639 ;  fResolutionConstants[5] = 0.07943    ; fResolutionConstants[6] = 6.873e-04;
    break;
   case(kptbin8):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] =  0.4514*204.8 ; fResolutionConstants[3] = 0.98  ;   fResolutionConstants[5] = 0.1192     ; fResolutionConstants[6] = 8.646e-04;
    break;
   case(kptbin9):
    fResolutionConstants[0]  = 0.326     ; fResolutionConstants[2] = 0.3622*502.8 ; fResolutionConstants[4] = 1.171*502.8   ;//from fit integrated in pt
    fResolutionConstants[1] =  0.525*181.   ; fResolutionConstants[3] = 0.99  ;   fResolutionConstants[5] = 0.1097     ; fResolutionConstants[6] = 9.637e-04;
    break;
   }
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
  return EvaluateCDFDecayTimeSigDistr(x)*EvaluateCDFInvMassSigDistr(m); 
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeSigDistr(Double_t x) const
{
//
// Implementation of the Background part of the Likelyhood function
// 
//
  Double_t retvalue = 0.;
  retvalue = fParameters[7]*FunB(x) + (1. - fParameters[7])*FunP(x);
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
   Double_t fitvalCB = 0.;
   Double_t normFactorCB = 1.;
   Double_t arg = (m - fParameters[9])/fParameters[11];
   Double_t numfactor = fParameters[10];
   Double_t denomfactor = numfactor - TMath::Abs(fParameters[12]) - arg;

   if(arg <= -1*TMath::Abs(fParameters[12])){
      Double_t exponent = fParameters[10]*TMath::Abs(fParameters[12]);
      Double_t numer = TMath::Exp(-0.5*fParameters[12]*fParameters[12])*TMath::Power(numfactor,exponent);
      Double_t denom = TMath::Power(denomfactor,exponent);
      fitvalCB += numer/denom;
      }
   if(arg > -1*TMath::Abs(fParameters[12])){
      fitvalCB += TMath::Exp(-0.5*arg*arg);
      }
   fitval = normFactorCB*fitvalCB;
   return fitval;
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
  Double_t np = 50.0;
  Double_t sc = 100.;
  Double_t sigma3 =  fResolutionConstants[2];
  Double_t xx;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  xlow = x - sc * sigma3 ;
  xupp = x + sc * sigma3 ;
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      sum += CsiMC(xx) * ResolutionFunc(xx-x);

      xx = xupp - (i-.5) * step;
      sum += CsiMC(xx) * ResolutionFunc(xx-x);
      }
  return (step * sum) ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunP(Double_t x) const 
//
//  Parameterisation of the Prompt part for the x distribution
//
{
  return ResolutionFunc(x);
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::CsiMC(Double_t x) const 
{
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
  return EvaluateCDFDecayTimeBkgDistr(x)*EvaluateCDFInvMassBkgDistr(m);
}  
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFDecayTimeBkgDistr(Double_t x) const 
{
// it returns the value of the probability to have a given x for the background 
//
//
  Double_t ret = (1 - fParameters[0]*fParameters[0])*ResolutionFunc(x);
  ret += FunBkgPos(x);
  ret += FunBkgNeg(x);
  ret += FunBkgSym(x);
  return ret;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::EvaluateCDFInvMassBkgDistr(Double_t m) const 
{
//
// it returns the value of the probability to have a given mass for the background
//
  return 1/(fMassWndHigh-fMassWndLow)+fParameters[6]*(m-(fMassWndHigh+fMassWndLow)/2);
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgPos(Double_t x) const 
{
//
// exponential with positive slopes for the background part (x)
//
  Double_t np = 50.0;
  Double_t sc = 100.;      
  Double_t sigma3 = fResolutionConstants[2];
  Double_t xx;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  xlow = x - sc * sigma3 ;
  xupp = x + sc * sigma3 ;
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      if (xx > 0) sum += (fParameters[0]*TMath::Cos(fParameters[1]))*(fParameters[0]*TMath::Cos(fParameters[1]))*fParameters[3]*TMath::Exp(-1*xx*fParameters[3]) * ResolutionFunc(xx-x);
 
      xx = xupp - (i-.5) * step;
      if (xx > 0) sum += fParameters[0]*TMath::Cos(fParameters[1])*(fParameters[0]*TMath::Cos(fParameters[1]))*fParameters[3]*TMath::Exp(-1*xx*fParameters[3]) * ResolutionFunc(xx-x);
      }
  return (step * sum) ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgNeg(Double_t x) const 
{
//
// exponential with negative slopes for the background part (x)
//
  Double_t np = 50.0;
  Double_t sc = 100.;      
  Double_t sigma3 =  fResolutionConstants[2];
  Double_t xx;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  xlow = x - sc * sigma3 ;
  xupp = x + sc * sigma3 ;
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      if (xx < 0) sum += (fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2]))*
                         fParameters[4]*TMath::Exp(xx*fParameters[4]) * ResolutionFunc(xx-x);
      xx = xupp - (i-.5) * step;
      if (xx < 0) sum += (fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2]))*
                         fParameters[4]*TMath::Exp(xx*fParameters[4]) * ResolutionFunc(xx-x);
      }
  return (step * sum) ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::FunBkgSym(Double_t x) const 
{
//
// exponential with both positive and negative slopes for the background part (x)
//
  Double_t np = 50.0;
  Double_t sc = 100.;      
  Double_t sigma3 =  fResolutionConstants[2];
  Double_t xx;
  Double_t sum1 = 0.0;
  Double_t sum2 = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  xlow = x - sc * sigma3 ;
  xupp = x + sc * sigma3 ;
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      if (xx > 0) sum1 += 0.5*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*
                              fParameters[5]*TMath::Exp(-1*xx*fParameters[5]) * ResolutionFunc(xx-x);

      xx = xupp - (i-.5) * step;
      if (xx > 0) sum1 += 0.5*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*
                              fParameters[5]*TMath::Exp(-1*xx*fParameters[5]) * ResolutionFunc(xx-x);
      }
  for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      if (xx < 0) sum2 += 0.5*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*
                              fParameters[5]*TMath::Exp(xx*fParameters[5]) * ResolutionFunc(xx-x);

      xx = xupp - (i-.5) * step;
      if (xx < 0) sum2 += 0.5*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*(fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2]))*
                              fParameters[5]*TMath::Exp(xx*fParameters[5]) * ResolutionFunc(xx-x);
      }
  return step*(sum1 + sum2) ;
}
//_________________________________________________________________________________________________
Double_t AliBtoJPSItoEleCDFfitFCN::ResolutionFunc(Double_t x) const 
{
  //
  //parametrization with 2 gaus + 1 exponential + 1 constant
  //
  Double_t arg=0;
  arg=x/fResolutionConstants[1];
  Double_t ret=TMath::Exp(-0.5*arg*arg);
  arg=x/fResolutionConstants[2];
  ret+=fResolutionConstants[3]*TMath::Exp(-0.5*arg*arg);
  arg=x/fResolutionConstants[4];
  if(x > 0)  { ret+=fResolutionConstants[5]*TMath::Exp(-arg); }
  if(x <= 0) { ret+=fResolutionConstants[5]*TMath::Exp(arg); }
  return fResolutionConstants[0]*(ret + fResolutionConstants[6]);
}
