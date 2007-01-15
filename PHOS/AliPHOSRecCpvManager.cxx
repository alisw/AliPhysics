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

/* $Id$ */

//_________________________________________________________________________
// Class for the management by the CPV reconstruction.
////                  
//*-- Author   : Boris Polichtchouk (IHEP, Protvino) 6 Mar 2001
//
// --- ROOT system ---

#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecCpvManager.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSLoader.h"

ClassImp(AliPHOSRecCpvManager) 

//____________________________________________________________________________
AliPHOSRecCpvManager::AliPHOSRecCpvManager() :
  fOneGamChisqCut(3.),
  fOneGamInitialStep(0.00005),
  fOneGamChisqMin(1.),
  fOneGamStepMin(0.0005),
  fOneGamNumOfIterations(50),
  fTwoGamInitialStep(0.00005),
  fTwoGamChisqMin(1.),
  fTwoGamEmin(0.1),
  fTwoGamStepMin(0.00005),
  fTwoGamNumOfIterations(50),  
  fThr0(0.),
  fSqdCut(0.) 
{
  // Put a comment here
  SetTitle("Cpv Reconstruction Manager");
}

AliPHOSRecCpvManager::~AliPHOSRecCpvManager(void) 
{
  // Put a comment here
}

Float_t AliPHOSRecCpvManager::Dispersion(Float_t etot, Float_t ai) const
{
  //"Dispresion" of energy deposition in the cell.
  // etot is the total shower energy, ai is the
  // calculated cell response,
  // ei is the measured cell response.

  const Float_t kConst = 1.5;
  return kConst*ai*(1.-ai/etot);
}

Float_t AliPHOSRecCpvManager::OneGamChi2(Float_t ai, Float_t ei, Float_t etot, Float_t& Gi) const
{
  //"Chi2" for one cell.
  // etot is the total "shower" energy, ai is the
  // calculated cell response,
  // ei is the measured cell response.

  const Float_t kConst = 1.5;

  Float_t da = ai - ei;
  Float_t d = kConst*ai*(1.-ai/etot);

  Float_t dd = da/d;
  Gi = dd*(2.- dd*kConst*(1.-2.*ai/etot));

  Info("OneGamChi2", " OneGamChi2 (ai,ei,etot,&Gi,chi2) %f %f %f %f %f", ai, ei, etot, Gi, da*da/d );

  return da*da/d;

}

Float_t AliPHOSRecCpvManager::TwoGamChi2(Float_t ai, Float_t ei, Float_t etot, Float_t& gi) const 
{
  // Put a comment here

  const Float_t kConst = 1.5;

  Float_t da = ai - ei;
  Float_t d = kConst*ai*(1.-ai/etot);

  Float_t dd = da/d;
  gi = dd*(2.- dd*kConst*(1.-2.*ai/etot));

  return da*da/d;
 
}

void AliPHOSRecCpvManager::AG(Float_t ei, Float_t xi, Float_t yi, Float_t& ai, Float_t& gxi, Float_t& gyi )
{
  //Calculates amplitude (ai) and gradients (gxi, gyi) of CPV pad response.
  //Integrated response (total "shower energy") is e, 
  //xi and yi are the distances along x and y from reference point 
  // to the pad center.

  AliPHOSGeometry * geom = AliPHOSLoader::GetPHOSGeometry();

  Float_t celZ = geom->GetPadSizeZ();
  Float_t celY = geom->GetPadSizePhi();

//  //    Info("AG", "celZ: %f celY: %f", celZ, celY) ;

  Float_t dx = celZ/2.;
  Float_t dy = celY/2.;

//  //    Float_t x = xi*celZ;
//  //    Float_t y = yi*celZ;

  Float_t x = xi*celZ;
  Float_t y = yi*celY;

  Float_t e = ei;

  Float_t a = Fcml(x+dx,y+dy) - Fcml(x+dx,y-dy) - Fcml(x-dx,y+dy) + Fcml(x-dx,y-dy);
  ai = a*e;


  Float_t gx = GradX(x+dx,y+dy) - GradX(x+dx,y-dy) - GradX(x-dx,y+dy) + GradX(x-dx,y-dy);
  gxi = gx*e*e;

  Float_t gy = GradY(x+dx,y+dy) - GradY(x+dx,y-dy) - GradY(x-dx,y+dy) + GradY(x-dx,y-dy);
  gyi = gy*e*e;

}

Float_t AliPHOSRecCpvManager::Fcml(Float_t x, Float_t y)
{
  //Cumulative function

  const Float_t ka = 1.0;
  const Float_t kb = 0.70;

  Float_t fff  = TMath::ATan(x*y/(  kb*TMath::Sqrt(  (kb*kb) + x*x+y*y)))
    - TMath::ATan(x*y/(3*kb*TMath::Sqrt((3*kb)*(3*kb) + x*x+y*y)))
    + TMath::ATan(x*y/(5*kb*TMath::Sqrt((5*kb)*(5*kb) + x*x+y*y))) 
    - TMath::ATan(x*y/(7*kb*TMath::Sqrt((7*kb)*(7*kb) + x*x+y*y)))
    + TMath::ATan(x*y/(9*kb*TMath::Sqrt((9*kb)*(9*kb) + x*x+y*y))); 
  
  Float_t fcml = ka*fff/TMath::TwoPi();
//    Info("Fcml", "fcml: %f", fcml) ;
  return fcml;

}


Float_t AliPHOSRecCpvManager::GradX(Float_t x, Float_t y)
{
  // Put a comment here

  const Float_t ka = 1.0;
  const Float_t kb = 0.70;

  Float_t skv      = kb*kb + x*x + y*y;

  Float_t gradient = y*(1.-x*x/skv)*  kb*TMath::Sqrt(skv)/( kb*kb*skv+x*x*y*y)
    - y*(1.-x*x/skv)*3*kb*TMath::Sqrt(skv)/((3*kb)*(3*kb)*skv+x*x*y*y)
    + y*(1.-x*x/skv)*5*kb*TMath::Sqrt(skv)/((5*kb)*(5*kb)*skv+x*x*y*y)
    - y*(1.-x*x/skv)*7*kb*TMath::Sqrt(skv)/((7*kb)*(7*kb)*skv+x*x*y*y)
    + y*(1.-x*x/skv)*9*kb*TMath::Sqrt(skv)/((9*kb)*(9*kb)*skv+x*x*y*y);
      
  Float_t grad    = ka*gradient/TMath::TwoPi();
  return grad;
}


Float_t AliPHOSRecCpvManager::GradY(Float_t x, Float_t y)
{
  // Put a comment here
  
  const Float_t ka = 1.0;
  const Float_t kb = 0.70;
 
  Float_t skv      = kb*kb + x*x + y*y;
  Float_t gradient = x*(1.-y*y/skv)*  kb*TMath::Sqrt(skv)/( kb*kb*skv+x*x*y*y)
    - x*(1.-y*y/skv)*3*kb*TMath::Sqrt(skv)/((3*kb)*(3*kb)*skv+x*x*y*y)
    + x*(1.-y*y/skv)*5*kb*TMath::Sqrt(skv)/((5*kb)*(5*kb)*skv+x*x*y*y)
    - x*(1.-y*y/skv)*7*kb*TMath::Sqrt(skv)/((7*kb)*(7*kb)*skv+x*x*y*y)
    + x*(1.-y*y/skv)*9*kb*TMath::Sqrt(skv)/((9*kb)*(9*kb)*skv+x*x*y*y);
  
  Float_t grad    = ka*gradient/TMath::TwoPi();
  return grad;
}


