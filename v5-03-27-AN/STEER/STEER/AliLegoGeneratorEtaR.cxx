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

/* $Id $ */

//------------------------------------------------------------------------
//    Lego generator in Eta bins
//    Uses geantino rays to check the material distributions and detector's
//    geometry
//    Author: A.Morsch 
//------------------------------------------------------------------------

#include "TRandom.h"

#include "AliLegoGeneratorEtaR.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliLegoGeneratorEtaR)


//___________________________________________
void AliLegoGeneratorEtaR::Generate()
{
// Create a geantino with kinematics corresponding to the current bins
// Here: Coor1 =  phi 
//       Coor2 =  eta.
   
  //
  // Rootinos are 0
   const Int_t kMpart = 0;
   Float_t orig[3], pmom[3];
   Float_t t, cost, sint, cosp, sinp;
   fCurCoor1 = fCoor1Min + (fCoor1Max - fCoor1Min) * gRandom->Rndm();
   fCurCoor2 = fCoor2Min + (fCoor2Max - fCoor2Min) * gRandom->Rndm();

   Float_t phi   = fCurCoor1 * TMath::Pi() / 180.;
   Float_t theta = 2. * TMath::ATan(TMath::Exp(-fCurCoor2));

   
   cost      = TMath::Cos(theta);
   sint      = TMath::Sin(theta);
   cosp      = TMath::Cos(phi);
   sinp      = TMath::Sin(phi);
   
   pmom[0] = cosp * sint;
   pmom[1] = sinp * sint;
   pmom[2] = cost;
   
   // --- Where to start
   orig[0] = orig[1] = orig[2] = 0;
   Float_t dalicz = 3000;
   if (fRadMin > 0) {
       t = PropagateCylinder(orig, pmom, fRadMin, dalicz);
       orig[0] = pmom[0]*t;
       orig[1] = pmom[1]*t;
       orig[2] = pmom[2]*t;
       if (TMath::Abs(orig[2]) > fZMax) return;
   }
   
   Float_t polar[3]={0.,0.,0.};
   Int_t ntr;
   gAlice->GetMCApp()->PushTrack(1, -1, kMpart, pmom, orig, polar, 0, kPPrimary, ntr);
   
}
