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
//                  
//*-- Author   : Boris Polichtchouk (IHEP, Protvino) 6 Mar 2001
//
// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecEmcManager.h"
#include "AliPHOS.h"
#include "AliRun.h"

ClassImp(AliPHOSRecEmcManager) 

//____________________________________________________________________________

  AliPHOSRecEmcManager::AliPHOSRecEmcManager()
{

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

Float_t AliPHOSRecEmcManager::Dispersion(Float_t Etot, Float_t Ai, Float_t Ei) const
{
  //"Dispresion" of energy deposition in the cell.
  // Etot is the total shower energy, Ai is the
  // calculated cell response,
  // Ei is the measured cell response.

  return Ei;
}

Float_t AliPHOSRecEmcManager::OneGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi)
{
  // Chi2 used in OneGam (one-gamma fitting).
  // Gi is d(Chi2)/d(Ai).

  Float_t da = Ai - Ei;
  Float_t D = Ei; // we assume that sigma(E) = sqrt(E)
  Gi = 2.*(Ai-Ei)/D;

  return da*da/D;

}

Float_t AliPHOSRecEmcManager::TwoGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi) const
{

  Float_t da = Ai - Ei;
  Float_t D = Ei; // we assume that sigma(E) = sqrt(E)
  Gi = 2.*(Ai-Ei)/D;

  return da*da/D;

}

void AliPHOSRecEmcManager::AG(Float_t Ei, Float_t Xi, Float_t Yi, Float_t& Ai, Float_t& GXi, Float_t& GYi )
{
  //Calculates amplitude (Ai) and gradients (GXi, GYi) of CPV pad response.
  //Integrated response (total "shower energy") is E, 
  //Xi and Yi are the distances along x and y from reference point 
  // to the pad center.
  //Shape of the shower is from PHOS TDR.


  Float_t r = TMath::Sqrt(Xi*Xi + Yi*Yi);
  Float_t r4    = r*r*r*r ;
  Float_t r295  = TMath::Power(r, 2.95) ;
  Float_t shape = Ei*TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  Ai = shape;


  //d(shape)/d(Xi)
  GXi = (-(TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2)*
       ((-0.006077944*Xi*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),
             0.4750000000000001))/
          TMath::Power(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475),2) - 
         (1.04*Xi*(TMath::Power(Xi,2) + TMath::Power(Yi,2)))/
          TMath::Power(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2),2))) - 
    4*Xi*(TMath::Power(Xi,2) + TMath::Power(Yi,2))*
     (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475)) + 
       1./(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2))))/
  TMath::Power(TMath::E(),TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2)*
    (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475)) + 
     1./(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2))));
  
  GXi = GXi*Ei;

  //d(shape)/d(Yi)
  GYi = (-(TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2)*
       ((-0.006077944*Yi*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),
             0.4750000000000001))/
          TMath::Power(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475),2) - 
         (1.04*Yi*(TMath::Power(Xi,2) + TMath::Power(Yi,2)))/
          TMath::Power(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2),2))) - 
    4*Yi*(TMath::Power(Xi,2) + TMath::Power(Yi,2))*
     (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475)) + 
       1./(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2))))/
  TMath::Power(TMath::E(),TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2)*
    (0.0316/(1 + 0.0652*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),1.475)) + 
     1./(2.32 + 0.26*TMath::Power(TMath::Power(Xi,2) + TMath::Power(Yi,2),2))));


  GYi = GYi*Ei;

}







