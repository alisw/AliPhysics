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


//_________________________________________________________________________
// Class for the management by the Emc reconstruction.
////                  
//*-- Author   : Boris Polichtchouk (IHEP, Protvino) 6 Mar 2001
//
// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecEmcManager.h"

ClassImp(AliPHOSRecEmcManager) 

//____________________________________________________________________________

  AliPHOSRecEmcManager::AliPHOSRecEmcManager()
{
  // default ctor
//    fOneGamChisqCut = 3.;
  fOneGamChisqCut = 1.3; // bvp 31.08.2001 

  fOneGamInitialStep = 0.00005;
  fOneGamChisqMin = 1.;
  fOneGamStepMin = 0.0005;
  fOneGamNumOfIterations = 50;

  fTwoGamInitialStep = 0.00005;
  fTwoGamChisqMin = 1.;
  fTwoGamEmin = 0.1;
  fTwoGamStepMin = 0.00005;
  fTwoGamNumOfIterations = 50;  

//    fThr0 = 0.6;
  fThr0 = 0.;
//    fSqdCut = 3.;
//    fSqdCut = 0.5; // bvp 31.08.2001
  fSqdCut = 0.;

  SetTitle("Emc Reconstruction Manager");

}

AliPHOSRecEmcManager::~AliPHOSRecEmcManager(void) {}

Float_t AliPHOSRecEmcManager::Dispersion(Float_t ei) const
{
  //"Dispresion" of energy deposition in the cell.
  // eTot is the total shower energy, ai is the
  // calculated cell response,
  // ei is the measured cell response.

  return ei;
}

Float_t AliPHOSRecEmcManager::OneGamChi2(Float_t ai, Float_t ei, Float_t fi, Float_t& gi) const
{
  // Chi2 used in OneGam (one-gamma fitting).
  // gi is d(Chi2)/d(ai).

  fi = 0 ; 
  Float_t da = ai - ei;
  Float_t d = ei; // we assume that sigma(E) = sqrt(E)
  gi = 2.*(ai-ei)/d;

  return da*da/d;

}

Float_t AliPHOSRecEmcManager::TwoGamChi2(Float_t ai, Float_t ei, Float_t fi, Float_t& gi) const
{
  // calculates chi^2
  fi = 0 ; 
  Float_t da = ai - ei;
  Float_t d = ei; // we assume that sigma(E) = sqrt(E)
  gi = 2.*(ai-ei)/d;

  return da*da/d;

}

void AliPHOSRecEmcManager::AG(Float_t ei, Float_t xi, Float_t yi, Float_t& ai, Float_t& gxi, Float_t& gyi )
{
  //Calculates amplitude (ai) and gradients (gxi, gyi) of CPV pad response.
  //Integrated response (total "shower energy") is E, 
  //xi and yi are the distances along x and y from reference point 
  // to the pad center.
  //Shape of the shower is from PHOS TDR.


  Float_t r = TMath::Sqrt(xi*xi + yi*yi);
  Float_t r4    = r*r*r*r ;
  Float_t r295  = TMath::Power(r, 2.95) ;
  Float_t shape = ei*TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  ai = shape;


  //d(shape)/d(xi)
  gxi = (-(TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2)*
       ((-0.006077944*xi*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),
             0.4750000000000001))/
          TMath::Power(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475),2) - 
         (1.04*xi*(TMath::Power(xi,2) + TMath::Power(yi,2)))/
          TMath::Power(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2),2))) - 
    4*xi*(TMath::Power(xi,2) + TMath::Power(yi,2))*
     (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475)) + 
       1./(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2))))/
  TMath::Power(TMath::E(),TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2)*
    (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475)) + 
     1./(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2))));
  
  gxi = gxi*ei;

  //d(shape)/d(yi)
  gyi = (-(TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2)*
       ((-0.006077944*yi*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),
             0.4750000000000001))/
          TMath::Power(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475),2) - 
         (1.04*yi*(TMath::Power(xi,2) + TMath::Power(yi,2)))/
          TMath::Power(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2),2))) - 
    4*yi*(TMath::Power(xi,2) + TMath::Power(yi,2))*
     (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475)) + 
       1./(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2))))/
  TMath::Power(TMath::E(),TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2)*
    (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),1.475)) + 
     1./(2.32 + 0.26*TMath::Power(TMath::Power(xi,2) + TMath::Power(yi,2),2))));


  gyi = gyi*ei;

}







