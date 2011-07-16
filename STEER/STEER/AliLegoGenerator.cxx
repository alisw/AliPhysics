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

//------------------------------------------------------------------------
//        Generic Lego generator code
//    Uses geantino rays to check the material distributions and detector's
//    geometry
//    Author: A.Morsch
//------------------------------------------------------------------------

#include "AliLegoGenerator.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliLegoGenerator)

//_______________________________________________________________________
AliLegoGenerator::AliLegoGenerator():
  fRadMin(0),
  fRadMax(0),
  fZMax(0),
  fNCoor1(0),
  fNCoor2(0),
  fCoor1Min(0),
  fCoor1Max(0),
  fCoor2Min(0),
  fCoor2Max(0),
  fCoor1Bin(-1),
  fCoor2Bin(-1),
  fCurCoor1(0),
  fCurCoor2(0)
{
  //
  // Default Constructor
  //
  SetName("Lego");
}

//_______________________________________________________________________
AliLegoGenerator::AliLegoGenerator(Int_t nc1, Float_t c1min,
                                   Float_t c1max, Int_t nc2, 
                                   Float_t c2min, Float_t c2max,
                                   Float_t rmin, Float_t rmax, Float_t zmax):
  AliGenerator(0), 
  fRadMin(rmin),
  fRadMax(rmax),
  fZMax(zmax),
  fNCoor1(nc1),
  fNCoor2(nc2),
  fCoor1Min(0),
  fCoor1Max(0),
  fCoor2Min(0),
  fCoor2Max(0),
  fCoor1Bin(nc1),
  fCoor2Bin(-1),
  fCurCoor1(0),
  fCurCoor2(0)
{
  //
  // Standard generator for Lego rays
  //
  SetName("Lego");
  SetCoor1Range(nc1, c1min, c1max);
  SetCoor2Range(nc2, c2min, c2max);
}

//_______________________________________________________________________
void AliLegoGenerator::Generate()
{
  // Create a geantino with kinematics corresponding to the current bins
  // Here: Coor1 =  theta 
  //       Coor2 =  phi.
  
  //
  // Rootinos are 0
   const Int_t kMpart = 0;
   Float_t orig[3], pmom[3];
   Float_t t, cost, sint, cosp, sinp;
   if (fCoor1Bin==-1) fCoor1Bin=fNCoor1;
   // Prepare for next step
   if(fCoor1Bin>=fNCoor1-1)
     if(fCoor2Bin>=fNCoor2-1) {
       AliWarning("End of Lego Generation");
       return;
     } else { 
       fCoor2Bin++;
       AliDebug(1, Form("Generating rays in phi bin:%d",fCoor2Bin));
       fCoor1Bin=0;
     } else fCoor1Bin++;

   fCurCoor1 = (fCoor1Min+(fCoor1Bin+0.5)*(fCoor1Max-fCoor1Min)/fNCoor1);
   fCurCoor2 = (fCoor2Min+(fCoor2Bin+0.5)*(fCoor2Max-fCoor2Min)/fNCoor2);
   cost      = TMath::Cos(fCurCoor1 * TMath::Pi()/180.);
   sint      = TMath::Sin(fCurCoor1 * TMath::Pi()/180.);
   cosp      = TMath::Cos(fCurCoor2 * TMath::Pi()/180.);
   sinp      = TMath::Sin(fCurCoor2 * TMath::Pi()/180.);
   
   pmom[0] = cosp*sint;
   pmom[1] = sinp*sint;
   pmom[2] = cost;
   
   // --- Where to start
   orig[0] = orig[1] = orig[2] = 0;
   Float_t dalicz = 3000;
   if (fRadMin > 0) {
       t = PropagateCylinder(orig,pmom,fRadMin,dalicz);
       orig[0] = pmom[0]*t;
       orig[1] = pmom[1]*t;
       orig[2] = pmom[2]*t;
       if (TMath::Abs(orig[2]) > fZMax) return;
   }
   
   Float_t polar[3]={0.,0.,0.};
   Int_t ntr;
   gAlice->GetMCApp()->PushTrack(1, -1, kMpart, pmom, orig, polar, 0, kPPrimary, ntr);
   
}

//_______________________________________________________________________
Float_t AliLegoGenerator::PropagateCylinder(Float_t *x, Float_t *v, Float_t r, 
                                            Float_t z)
{
  //
  // Propagate to cylinder from inside
  //
   Double_t hnorm, sz, t, t1, t2, t3, sr;
   Double_t d[3];
   const Float_t kSmall  = 1e-8;
   const Float_t kSmall2 = kSmall*kSmall;

// ---> Find intesection with Z planes
   d[0]  = v[0];
   d[1]  = v[1];
   d[2]  = v[2];
   hnorm = TMath::Sqrt(1/(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]));
   d[0] *= hnorm;
   d[1] *= hnorm;
   d[2] *= hnorm;
   if (d[2] > kSmall)       sz = (z-x[2])/d[2];
   else if (d[2] < -kSmall) sz = -(z+x[2])/d[2];
   else                     sz = 1.e10;  // ---> Direction parallel to X-Y, no intersection

// ---> Intersection with cylinders
//      Intersection point (x,y,z)
//      (x,y,z) is on track :    x=X(1)+t*D(1)
//                               y=X(2)+t*D(2)
//                               z=X(3)+t*D(3)
//      (x,y,z) is on cylinder : x**2 + y**2 = R**2
//
//      (D(1)**2+D(2)**2)*t**2
//      +2.*(X(1)*D(1)+X(2)*D(2))*t
//      +X(1)**2+X(2)**2-R**2=0
// ---> Solve second degree equation
   t1 = d[0]*d[0] + d[1]*d[1];
   if (t1 <= kSmall2) {
      t = sz;  // ---> Track parallel to the z-axis, take distance to planes
   } else {
      t2 = x[0]*d[0] + x[1]*d[1];
      t3 = x[0]*x[0] + x[1]*x[1];
      // ---> It should be positive, but there may be numerical problems
      sr = (-t2 +TMath::Sqrt(TMath::Max(t2*t2-(t3-r*r)*t1,0.)))/t1;
      // ---> Find minimum distance between planes and cylinder
      t  = TMath::Min(sz,sr);
   }
   return t;
}


