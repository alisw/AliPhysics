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

/*
$Log$
Revision 1.7  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//                            
//  Utility class to evaluate the material budget from
//  a given radius to the surface of an arbitrary cylinder
//  along radial directions from the centre:
// 
//   - radiation length
//   - Interaction length
//   - g/cm2
// 
//  Geantinos are shot in the bins in the fNtheta bins in theta
//  and fNphi bins in phi with specified rectangular limits.
//  The statistics are accumulated per
//    fRadMin < r < fRadMax    and  <0 < z < fZMax
//
//  To activate this option, you can do:
//      Root > gAlice.RunLego();
//  or  Root > .x menu.C  then select menu item "RunLego"
//  Note that when calling gAlice->RunLego, an optional list
//  of arguments may be specified.
//
//Begin_Html
/*
<img src="picts/alilego.gif">
*/
//End_Html
//
//////////////////////////////////////////////////////////////

#include "TMath.h"
#include "AliLego.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliLego)


//___________________________________________
AliLego::AliLego()
{
   fHistRadl = 0;
   fHistAbso = 0;
   fHistGcm2 = 0;
   fHistReta = 0;
}

//___________________________________________
AliLego::AliLego(const char *name, const char *title) 
        : TNamed(name,title)
{
   fHistRadl = 0;
   fHistAbso = 0;
   fHistGcm2 = 0;
   fHistReta = 0;
}

//___________________________________________
AliLego::~AliLego()
{
   delete fHistRadl;
   delete fHistAbso;
   delete fHistGcm2;
   delete fHistReta;
}

//___________________________________________
void AliLego::GenerateKinematics()
{
// Create a geantino with kinematics corresponding to the current
// bins in theta and phi.
   
  //
  // Rootinos are 0
   const Int_t mpart = 0;
   Float_t orig[3], pmom[3];
   Float_t t, cost, sint, cosp, sinp;
   
// --- Set to 0 radiation length, absorption length and g/cm2 ---
   fTotRadl = 0;
   fTotAbso = 0;
   fTotGcm2 = 0;

   fCurTheta = (fThetaMin+(fThetaBin-0.5)*(fThetaMax-fThetaMin)/fNtheta)*kDegrad;
   fCurPhi   = (fPhiMin+(fPhiBin-0.5)*(fPhiMax-fPhiMin)/fNphi)*kDegrad;
   cost      = TMath::Cos(fCurTheta);
   sint      = TMath::Sin(fCurTheta);
   cosp      = TMath::Cos(fCurPhi);
   sinp      = TMath::Sin(fCurPhi);

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

// --- We do start here
   Float_t polar[3]={0.,0.,0.};
   Int_t ntr;
   gAlice->SetTrack(1, 0, mpart, pmom, orig, polar, 0, "LEGO ray", ntr);
}

//___________________________________________
void AliLego::Init(Int_t ntheta,Float_t themin,Float_t themax,
          Int_t nphi,Float_t phimin,Float_t phimax,Float_t rmin,Float_t rmax,
		  Float_t zmax)
{
// specify the angular limits and the size of the rectangular box
   fNtheta   = ntheta;
   fThetaMin = themin;
   fThetaMax = themax;
   fNphi     = nphi;
   fPhiMin   = phimin;
   fPhiMax   = phimax;
   fRadMin   = rmin;
   fRadMax   = rmax;
   fZMax     = zmax;
   Float_t etamin = -TMath::Log(TMath::Tan(TMath::Min((Double_t)fThetaMax*kDegrad/2,TMath::Pi()/2-1.e-10)));
   Float_t etamax = -TMath::Log(TMath::Tan(TMath::Max((Double_t)fThetaMin*kDegrad/2,              1.e-10)));

   fHistRadl = new TH2F("hradl","Radiation length map",    nphi,phimin,phimax,ntheta,themin,themax);
   fHistAbso = new TH2F("habso","Interaction length map",  nphi,phimin,phimax,ntheta,themin,themax);
   fHistGcm2 = new TH2F("hgcm2","g/cm2 length map",        nphi,phimin,phimax,ntheta,themin,themax);
   fHistReta = new TH2F("hetar","Radiation length vs. eta",nphi,phimin,phimax,ntheta,etamin,etamax);

}

//___________________________________________
Float_t AliLego::PropagateCylinder(Float_t *x, Float_t *v, Float_t r, Float_t z)
{
// Propagate to cylinder from inside

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
      sr = (t2 +TMath::Sqrt(TMath::Max(t2*t2-(t3-r*r)*t1,0.)))/t1;
      // ---> Find minimum distance between planes and cylinder
      t  = TMath::Min(sz,sr);
   }
   return t;
}

//___________________________________________
void AliLego::Run()
{
   // loop on phi,theta bins
   gMC->InitLego();
   Float_t thed, phid, eta;
   for (fPhiBin=1; fPhiBin<=fNphi; fPhiBin++) {
      printf("AliLego Generating rays in phi bin:%d\n",fPhiBin);
      for (fThetaBin=1; fThetaBin<=fNtheta; fThetaBin++) {
         gMC->Gtrigi();
         gMC->Gtrigc();
         GenerateKinematics();
         gMC->Gtreve_root();

         thed = fCurTheta*kRaddeg;
         phid = fCurPhi*kRaddeg;
	 eta  = -TMath::Log(TMath::Tan(TMath::Max(
                     TMath::Min((Double_t)fCurTheta/2,TMath::Pi()/2-1.e-10),1.e-10)));
         fHistRadl->Fill(phid,thed,fTotRadl);
         fHistAbso->Fill(phid,thed,fTotAbso);
         fHistGcm2->Fill(phid,thed,fTotGcm2);
	 fHistReta->Fill(phid,eta,fTotRadl);
         gAlice->FinishEvent();
      }
   }
   // store histograms in current Root file
   fHistRadl->Write();
   fHistAbso->Write();
   fHistGcm2->Write();
   fHistReta->Write();
}

//___________________________________________
void AliLego::StepManager()
{
// called from AliRun::Stepmanager from gustep.
// Accumulate the 3 parameters step by step
  
   Float_t t, tt;
   Float_t a,z,dens,radl,absl;
   Int_t i;
   
   Float_t step  = gMC->TrackStep();
       
   Float_t vect[3], dir[3];
   TLorentzVector pos, mom;

   gMC->TrackPosition(pos);  
   gMC->TrackMomentum(mom);
   gMC->CurrentMaterial(a,z,dens,radl,absl);
   
   if (z < 1) return;
   
// --- See if we have to stop now
   if (TMath::Abs(pos[2]) > fZMax  || 
       pos[0]*pos[0] +pos[1]*pos[1] > fRadMax*fRadMax) {
     gMC->StopEvent();
     return;
   }

// --- See how long we have to go
   for(i=0;i<3;++i) {
     vect[i]=pos[i];
     dir[i]=mom[i];
   }
   t  = PropagateCylinder(vect,dir,fRadMax,fZMax);
   tt = TMath::Min(step,t);

   fTotAbso += tt/absl;
   fTotRadl += tt/radl;
   fTotGcm2 += tt*dens;
}

