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
#include "AliBtoJPSItoEleCDFfitHandler.h"
#include "AliBtoJPSItoEleCDFfitFCN.h"

//-------------------------------------------------------------------------
//                      Class AliBtoJPSItoEleCDFfitHandler
//            Class to perform unbinned log-likelihood fit
//      
//                        Origin: C. Di Giglio
//     Contact: carmelo.digiglio@ba.infn.it , giuseppe.bruno@ba.infn.it
//-------------------------------------------------------------------------

void CDFFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag);

ClassImp(AliBtoJPSItoEleCDFfitHandler)

//______________________________________________________________________________
void CDFFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // This function is called by minuit
  // The corresponding member method is called
  // using SetObjectFit/GetObjectFit methods of TMinuit
  AliBtoJPSItoEleCDFfitHandler* dummy = (AliBtoJPSItoEleCDFfitHandler *)TVirtualFitter::GetFitter()->GetObjectFit();
  dummy->CdfFCN(npar, gin, f, par, iflag);
}


//_________________________________________________________________________________________________
AliBtoJPSItoEleCDFfitHandler::AliBtoJPSItoEleCDFfitHandler():
fIsParamFixed(16),
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
AliBtoJPSItoEleCDFfitHandler::AliBtoJPSItoEleCDFfitHandler(Double_t* decaytime, 
  Double_t* invariantmass, Int_t ncand) :
fIsParamFixed(16),
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
  AliInfo("\n+++\n+++ Minimization object AliBtoJPSItoEleCDFfitHandler created\n+++\n");
  fLikely = new AliBtoJPSItoEleCDFfitFCN();
  AliInfo("\n+++\n+++ CDF fit function object AliBtoJPSItoEleCDFfitFCN created\n+++\n");
  AliInfo("Parameter 0  ----> fRadius");
  AliInfo("Parameter 1  ----> fTheta");
  AliInfo("Parameter 2  ----> fPhi");
  AliInfo("Parameter 3  ----> fOneOvLamPlus");
  AliInfo("Parameter 4  ----> fOneOvLamMinus");
  AliInfo("Parameter 5  ----> fOneOvLamSym");
  AliInfo("Parameter 6  ----> fMSlope");
  AliInfo("Parameter 7  ----> fB");
  AliInfo("Parameter 8  ----> fFsig");
  AliInfo("Parameter 9  ----> fMmean");
  AliInfo("Parameter 10 ----> fNexp");
  AliInfo("Parameter 11 ----> fSigma");
  AliInfo("Parameter 12 ----> fAlpha");
  AliInfo("Parameter 13 ----> fNorm");
  AliInfo("Parameter 14 ----> fSigmaResol");
  AliInfo("Parameter 15 ----> fNResol");
  AliInfo(Form("\n+++\n+++ Number of candidates ---> %d\n+++\n ", ncand));
}
//___________________________________________________________________________
AliBtoJPSItoEleCDFfitHandler& AliBtoJPSItoEleCDFfitHandler::operator=(const AliBtoJPSItoEleCDFfitHandler& c)
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
AliBtoJPSItoEleCDFfitHandler::AliBtoJPSItoEleCDFfitHandler(const AliBtoJPSItoEleCDFfitHandler& c) :
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
AliBtoJPSItoEleCDFfitHandler::~AliBtoJPSItoEleCDFfitHandler()
{
  //
  //destructor
  //
  delete fLikely;
}
//_______________________________________________________________________________________
Int_t AliBtoJPSItoEleCDFfitHandler::DoMinimization()
{
  //
  // performs the minimization
  //
  static TVirtualFitter *fitter = TVirtualFitter::Fitter(this,16);
  fitter->SetFCN(CDFFunction);

  fitter->SetParameter(0,"fRadius",fParamStartValues[0], 1.e-06, 0., 1.);
  fitter->SetParameter(1,"fTheta",fParamStartValues[1], 1.e-06, 0.,2*TMath::Pi());
  fitter->SetParameter(2,"fPhi",fParamStartValues[2], 1.e-06, 0.,2*TMath::Pi());
//  fitter->SetParameter(3,"fOneOvLamPlus",fParamStartValues[3], 1.e-10, 0., 5.e+01);
//  fitter->SetParameter(4,"fOneOvLamMinus",fParamStartValues[4], 1.e-10, 0., 5.e+01);
//  fitter->SetParameter(5,"fOneOvLamSym",fParamStartValues[5], 1.e-10, 0., 5.e+01);
  fitter->SetParameter(3,"fOneOvLamPlus",fParamStartValues[3], 1.e-10, 0.0000001, 5.e+01);
  fitter->SetParameter(4,"fOneOvLamMinus",fParamStartValues[4], 1.e-10, 0.00000001, 5.e+01);
  fitter->SetParameter(5,"fOneOvLamSym",fParamStartValues[5], 1.e-10, 0.00000001, 5.e+01);
  fitter->SetParameter(6,"fMSlope",fParamStartValues[6], 1.e-04, -2.5, 2.5);
  fitter->SetParameter(7,"fB",fParamStartValues[7], 1.e-08, 0., 1.);
  fitter->SetParameter(8,"fFsig",fParamStartValues[8], 1.e-08, 0., 1.);
  fitter->SetParameter(9,"fMmean",fParamStartValues[9], 1.e-08, 0., 1.e+04);
  fitter->SetParameter(10,"fNexp",fParamStartValues[10], 1.e-08, 0., 1.e+02);
  fitter->SetParameter(11,"fSigma",fParamStartValues[11], 1.e-08, 0., 1.e+04);
  fitter->SetParameter(12,"fAlpha",fParamStartValues[12], 1.e-08, 0., 1.e+04);
  fitter->SetParameter(13,"fNorm",fParamStartValues[13], 1.e-08, 0., 1.e+01);
  fitter->SetParameter(14,"fSigmaResol",fParamStartValues[14], 1.e-08, 0., 1.e+04);
  fitter->SetParameter(15,"fNResol",fParamStartValues[15], 1.e-08, 0., 1.e+05);

  for(UInt_t indexparam = 0; indexparam < 16; indexparam++){
     if(IsParamFixed(indexparam))fitter->FixParameter((Int_t)indexparam);
  }

  Double_t arglist[2]={10000,1.0}; 
  Int_t iret=fitter->ExecuteCommand("MIGRAD", arglist ,2);
  fitter->PrintResults(4,0);

  AliInfo("Minimization procedure finished\n");

  return iret;
}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::CdfFCN(Int_t & /* npar */, 
              Double_t * /* gin */,Double_t &f,Double_t *par,Int_t /* iflag */)
{
// 
// Definition of the FCN to be used by minuit
//
  fLikely->SetAllParameters(par);
  fLikely->ConvertFromSpherical();
  fLikely->ComputeIntegral();
  if(fPrintStatus)fLikely->PrintStatus();

  TStopwatch t;
  t.Start();

  f = fLikely->EvaluateLikelihood(fX,fM,fNcand);

  t.Stop();
  AliDebug(2,Form("Real time spent to calculate function == %f \n", t.RealTime()));
  AliDebug(2,Form("CPU time spent to calculate function == %f \n", t.CpuTime()));
  AliDebug(2,Form("Actual value of the AliBtoJPSItoEleCDFfitFCN function == %f \n", f));

  return;
}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::SetParamStartValues(Double_t inputparamvalues[16])
{
  for(Int_t index=0; index < 16; index++) fParamStartValues[index] = inputparamvalues[index];
}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::SetResolutionConstants()
{
  //
  // Sets constants for the resolution function
  //
  fLikely->SetResolutionConstants();

}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::SetCrystalBallFunction(Bool_t okCB) 
{
  //
  // Sets the CB as the parametrization for the signal invariant mass spectrum 
  // (otherwise Landau is chosen)
  //
  fLikely->SetCrystalBallFunction(okCB);
}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::SetMassWndHigh(Double_t limit) 
{ 
  //
  // Sets upper limit for the invariant mass window (under J/PSI mass peak)
  //
  fLikely->SetMassWndHigh(limit);
}
//_______________________________________________________________________________________
void AliBtoJPSItoEleCDFfitHandler::SetMassWndLow(Double_t limit)
{
  //
  // Sets lower limit for the invariant mass window (under J/PSI mass peak)
  //
  fLikely->SetMassWndLow(limit);
}

