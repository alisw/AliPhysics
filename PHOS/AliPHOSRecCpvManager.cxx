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
// Class for the management by the CPV reconstruction.
//                  
//*-- Author   : Boris Polichtchouk (IHEP, Protvino) 6 Mar 2001
//
// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecCpvManager.h"
#include "AliPHOS.h"
#include "AliRun.h"
#include "AliPHOSGetter.h"

ClassImp(AliPHOSRecCpvManager) 

//____________________________________________________________________________

  AliPHOSRecCpvManager::AliPHOSRecCpvManager()
{

  fOneGamChisqCut = 3.; // If Chi2/dof > fOneGamChisqCut, split point.

  fOneGamInitialStep = 0.00005;
  fOneGamChisqMin = 1.;
  fOneGamStepMin = 0.0005;
  fOneGamNumOfIterations = 50;

  fTwoGamInitialStep = 0.00005;
  fTwoGamChisqMin = 1.;
  fTwoGamEmin = 0.1;
  fTwoGamStepMin = 0.00005;
  fTwoGamNumOfIterations = 50;  

//    fThr0 = 0.0285; // Min. energy of rec. point. If E<fThr0, point deleted.
//    fSqdCut = 3.; // Min. distance (in cm) between two rec points.

  fThr0 = 0.; 
  fSqdCut = 0.; 
  
  SetTitle("Cpv Reconstruction Manager");
}

AliPHOSRecCpvManager::~AliPHOSRecCpvManager(void) {}

Float_t AliPHOSRecCpvManager::Dispersion(Float_t Etot, Float_t Ai, Float_t Ei)
{
  //"Dispresion" of energy deposition in the cell.
  // Etot is the total shower energy, Ai is the
  // calculated cell response,
  // Ei is the measured cell response.

  const Float_t Const = 1.5;
  return Const*Ai*(1.-Ai/Etot);
}

Float_t AliPHOSRecCpvManager::OneGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi)
{
  //"Chi2" for one cell.
  // Etot is the total "shower" energy, Ai is the
  // calculated cell response,
  // Ei is the measured cell response.

  const Float_t Const = 1.5;

  Float_t da = Ai - Ei;
  Float_t D = Const*Ai*(1.-Ai/Etot);

  Float_t dd = da/D;
  Gi = dd*(2.- dd*Const*(1.-2.*Ai/Etot));

  Info("OneGamChi2", " OneGamChi2 (Ai,Ei,Etot,&Gi,chi2) %f %f %f %f %f", Ai, Ei, Etot, Gi, da*da/D );

  return da*da/D;

}

Float_t AliPHOSRecCpvManager::TwoGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi)
{

  const Float_t Const = 1.5;

  Float_t da = Ai - Ei;
  Float_t D = Const*Ai*(1.-Ai/Etot);

  Float_t dd = da/D;
  Gi = dd*(2.- dd*Const*(1.-2.*Ai/Etot));

  return da*da/D;
 
}

void AliPHOSRecCpvManager::AG(Float_t Ei, Float_t Xi, Float_t Yi, Float_t& Ai, Float_t& GXi, Float_t& GYi )
{
  //Calculates amplitude (Ai) and gradients (GXi, GYi) of CPV pad response.
  //Integrated response (total "shower energy") is E, 
  //Xi and Yi are the distances along x and y from reference point 
  // to the pad center.

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  const AliPHOSGeometry* geom = gime->PHOSGeometry();
  Float_t CelZ = geom->GetPadSizeZ();
  Float_t CelY = geom->GetPadSizePhi();

//  //    Info("AG", "CelZ: %f CelY: %f", CelZ, CelY) ;

  Float_t dx = CelZ/2.;
  Float_t dy = CelY/2.;

//  //    Float_t x = Xi*CelZ;
//  //    Float_t y = Yi*CelZ;

  Float_t x = Xi*CelZ;
  Float_t y = Yi*CelY;

  Float_t E = Ei;

  Float_t A = Fcml(x+dx,y+dy) - Fcml(x+dx,y-dy) - Fcml(x-dx,y+dy) + Fcml(x-dx,y-dy);
  Ai = A*E;


  Float_t Gx = GradX(x+dx,y+dy) - GradX(x+dx,y-dy) - GradX(x-dx,y+dy) + GradX(x-dx,y-dy);
  GXi = Gx*E*E;

  Float_t Gy = GradY(x+dx,y+dy) - GradY(x+dx,y-dy) - GradY(x-dx,y+dy) + GradY(x-dx,y-dy);
  GYi = Gy*E*E;

}

Float_t AliPHOSRecCpvManager::Fcml(Float_t x, Float_t y)
{
  //Cumulative function

  const Float_t A = 1.0;
  const Float_t b = 0.70;

  Float_t Fff  = TMath::ATan(x*y/(  b*TMath::Sqrt(  (b*b) + x*x+y*y)))
    - TMath::ATan(x*y/(3*b*TMath::Sqrt((3*b)*(3*b) + x*x+y*y)))
    + TMath::ATan(x*y/(5*b*TMath::Sqrt((5*b)*(5*b) + x*x+y*y))) 
    - TMath::ATan(x*y/(7*b*TMath::Sqrt((7*b)*(7*b) + x*x+y*y)))
    + TMath::ATan(x*y/(9*b*TMath::Sqrt((9*b)*(9*b) + x*x+y*y))); 
  
  Float_t Fcml = A*Fff/6.2831853071796;
//    Info("Fcml", "Fcml: %f", Fcml) ;
  return Fcml;

}


Float_t AliPHOSRecCpvManager::GradX(Float_t x, Float_t y)
{

  const Float_t A = 1.0;
  const Float_t b = 0.70;

  Float_t skv      = b*b + x*x + y*y;

  Float_t Gradient = y*(1.-x*x/skv)*  b*TMath::Sqrt(skv)/( b*b*skv+x*x*y*y)
    - y*(1.-x*x/skv)*3*b*TMath::Sqrt(skv)/((3*b)*(3*b)*skv+x*x*y*y)
    + y*(1.-x*x/skv)*5*b*TMath::Sqrt(skv)/((5*b)*(5*b)*skv+x*x*y*y)
    - y*(1.-x*x/skv)*7*b*TMath::Sqrt(skv)/((7*b)*(7*b)*skv+x*x*y*y)
    + y*(1.-x*x/skv)*9*b*TMath::Sqrt(skv)/((9*b)*(9*b)*skv+x*x*y*y);
      
  Float_t Grad    = A*Gradient/6.2831853071796;
  return Grad;
}


Float_t AliPHOSRecCpvManager::GradY(Float_t x, Float_t y)
{
  
  const Float_t A = 1.0;
  const Float_t b = 0.70;
 
  Float_t skv      = b*b + x*x + y*y;
  Float_t Gradient = x*(1.-y*y/skv)*  b*TMath::Sqrt(skv)/( b*b*skv+x*x*y*y)
    - x*(1.-y*y/skv)*3*b*TMath::Sqrt(skv)/((3*b)*(3*b)*skv+x*x*y*y)
    + x*(1.-y*y/skv)*5*b*TMath::Sqrt(skv)/((5*b)*(5*b)*skv+x*x*y*y)
    - x*(1.-y*y/skv)*7*b*TMath::Sqrt(skv)/((7*b)*(7*b)*skv+x*x*y*y)
    + x*(1.-y*y/skv)*9*b*TMath::Sqrt(skv)/((9*b)*(9*b)*skv+x*x*y*y);
  
  Float_t Grad    = A*Gradient/6.2831853071796;
  return Grad;
}


