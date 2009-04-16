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
  AliInfo("\n+++\n+++ Creating AliBtoJPSItoEleCDFfitFCN object\n+++\n");
  fLikely = new AliBtoJPSItoEleCDFfitFCN();
  fLikely->SetCrystalBallParam(kFALSE); //Landau selected; otherwise Crystal Ball is selected
  SetErrorDef(1.);
  AliInfo(Form("\n+++\n+++ Number of candidates ---> %d\n+++\n ", ncand));
}
//___________________________________________________________________________
AliBtoJPSItoEleCDFfitHandler& AliBtoJPSItoEleCDFfitHandler::operator=(const AliBtoJPSItoEleCDFfitHandler& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
      fUp     = c.fUp;
      fX      = c.fX;
      fM      = c.fM;
      fLikely = c.fLikely;
      fNcand = c.fNcand;
     }
  return *this;
}

//___________________________________________________________________________
AliBtoJPSItoEleCDFfitHandler::AliBtoJPSItoEleCDFfitHandler(const AliBtoJPSItoEleCDFfitHandler& c) :
TNamed(c),
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
//________________________________________________________________________
AliBtoJPSItoEleCDFfitHandler::~AliBtoJPSItoEleCDFfitHandler()
{
  //
  //destructor
  //
  delete fLikely;
}
//_________________________________________________________________________________________________
Int_t AliBtoJPSItoEleCDFfitHandler::DoMinimization()
{
  //
  // performs the minimization
  //
  static TVirtualFitter *fitter = TVirtualFitter::Fitter(this,13);
  fitter->SetFCN(CDFFunction);
  Double_t startingParamValues[13] =
         /* startfPlus
            startfMinus
            startfSym
            startfOneOvLamPlus
            startfOneOvLamMinus
            startfOneOvLamSym
            startfMSlope
            startfB
            startfFsig
            startfMmean
            startfNexp
            startfSigma
            startfAlpha */
            {5.00e-01,
             TMath::Pi()/4.,
             TMath::Pi()/4.,
             2.0964360e-03,
             4.8309180e-03,
             1.582530e-04,
             -1.5720e-02,
             0.1800e+00,
             0.7000e+00,
            3.0910e+00,
             1.0500e+00,
             1.4250e-02,
             6.758e-01};

  fitter->SetParameter(0,"fRadius",startingParamValues[0], 0.01, 0., 1.);
  fitter->SetParameter(1,"fTheta",startingParamValues[1], 0.001, 0., TMath::Pi()/2);
  fitter->SetParameter(2,"fPhi",startingParamValues[2], 0.001, 0., TMath::Pi()/2);
  fitter->SetParameter(3,"fOneOvLamPlus",startingParamValues[3], 0.0001, 0., 5.e-01);
  fitter->SetParameter(4,"fOneOvLamMinus",startingParamValues[4], 0.0001, 0., 5.e-01);
  fitter->SetParameter(5,"fOneOvLamSym",startingParamValues[5], 0.00001, 0., 5.e-01);
  fitter->SetParameter(6,"fMSlope",startingParamValues[6], 0.001, -1.e-00, 1.e+00);
  fitter->SetParameter(7,"fB",startingParamValues[7], 0.1, 0., 1.);
  fitter->SetParameter(8,"fFsig",startingParamValues[8], 0.1, 0., 1.);
  fitter->SetParameter(9,"fMmean",startingParamValues[9], 0.1, 0., 1.e+02);
  fitter->SetParameter(10,"fNexp",startingParamValues[10], 0.1, 0., 1.e+02);
  fitter->SetParameter(11,"fSigma",startingParamValues[11], 0.001, 0., 1.e+02);
  fitter->SetParameter(12,"fAlpha",startingParamValues[12], 0.01, 0., 1.e+02);
  fitter->FixParameter(9);

  Double_t arglist[2]={10000,0.5}; 
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
  //printf("\n+++\n+++\n+++\n");
  //fLikely->PrintStatus();

  TStopwatch t;
  t.Start();

  f = fLikely->EvaluateLikelihood(fX,fM,fNcand);

  t.Stop();
  AliDebug(2,Form("Real time spent to calculate function == %f \n", t.RealTime()));
  AliDebug(2,Form("CPU time spent to calculate function == %f \n", t.CpuTime()));
  AliDebug(2,Form("Actual value of the AliBtoJPSItoEleCDFfitFCN function == %f \n", f));

  return;
}
