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
Revision 1.20  2000/11/30 07:12:48  alibrary
Introducing new Rndm and QA classes

Revision 1.19  2000/10/26 14:13:05  morsch
- Change from coordinates theta, phi to general coordinates Coor1 and Coor2.
- Lego generator instance can be passed in constructor.

Revision 1.18  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.17  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

Revision 1.16  2000/05/26 08:35:03  fca
Move the check on z after z has been retrieved

Revision 1.15  2000/05/16 13:10:40  fca
New method IsNewTrack and fix for a problem in Father-Daughter relations

Revision 1.14  2000/04/27 10:38:21  fca
Correct termination of Lego Run and introduce Lego getter in AliRun

Revision 1.13  2000/04/26 10:17:31  fca
Changes in Lego for G4 compatibility

Revision 1.12  2000/04/07 11:12:33  fca
G4 compatibility changes

Revision 1.11  2000/03/22 13:42:26  fca
SetGenerator does not replace an existing generator, ResetGenerator does

Revision 1.10  2000/02/23 16:25:22  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.9  1999/12/03 10:54:01  fca
Fix lego summary

Revision 1.8  1999/10/01 09:54:33  fca
Correct logics for Lego StepManager

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
#include "AliLegoGenerator.h"
#include "AliConst.h"
#include "AliMC.h"
#include "TH2.h"

ClassImp(AliLego)


//___________________________________________
AliLego::AliLego()
{
  //
  // Default constructor
  //
  fHistRadl = 0;
  fHistAbso = 0;
  fHistGcm2 = 0;
  fHistReta = 0;
}

//___________________________________________
AliLego::AliLego(const char *title, Int_t ntheta, Float_t thetamin, 
		 Float_t thetamax, Int_t nphi, Float_t phimin, Float_t phimax,
		 Float_t rmin, Float_t rmax, Float_t zmax)
  : TNamed("Lego Generator",title)
{
  //
  // specify the angular limits and the size of the rectangular box
  //
   fGener = new AliLegoGenerator(ntheta, thetamin, thetamax,
		       nphi, phimin, phimax, rmin, rmax, zmax);
   

   
   fHistRadl = new TH2F("hradl","Radiation length map",    
			ntheta,thetamin,thetamax,nphi,phimin,phimax);
   fHistAbso = new TH2F("habso","Interaction length map",  
			ntheta,thetamin,thetamax,nphi,phimin,phimax);
   fHistGcm2 = new TH2F("hgcm2","g/cm2 length map",        
			ntheta,thetamin,thetamax,nphi,phimin,phimax);
}

AliLego::AliLego(const char *title, AliLegoGenerator* generator)
  : TNamed("Lego Generator",title)
{
  //
  // specify the angular limits and the size of the rectangular box
  //
   fGener = generator;
   Float_t c1min, c1max, c2min, c2max;
   Int_t n1 = fGener->NCoor1();
   Int_t n2 = fGener->NCoor2();
   
   fGener->Coor1Range(c1min, c1max);
   fGener->Coor2Range(c2min, c2max);   

   fHistRadl = new TH2F("hradl","Radiation length map",    
			n2, c2min, c2max, n1, c1min, c1max);
   fHistAbso = new TH2F("habso","Interaction length map",  
			n2, c2min, c2max, n1, c1min, c1max);
   fHistGcm2 = new TH2F("hgcm2","g/cm2 length map",        
			n2, c2min, c2max, n1, c1min, c1max);
}

//___________________________________________
AliLego::~AliLego()
{
  //
  // Destructor
  //
  delete fHistRadl;
  delete fHistAbso;
  delete fHistGcm2;
  delete fGener;
}

//___________________________________________
void AliLego::BeginEvent()
{
  //
  // --- Set to 0 radiation length, absorption length and g/cm2 ---
  //
  fTotRadl = 0;
  fTotAbso = 0;
  fTotGcm2 = 0;
}

//___________________________________________
void AliLego::FinishEvent()
{
  //
  // Finish the event and update the histos
  //
  Double_t c1, c2;
  c1 = fGener->CurCoor1();
  c2 = fGener->CurCoor2();
  fHistRadl->Fill(c2,c1,fTotRadl);
  fHistAbso->Fill(c2,c1,fTotAbso);
  fHistGcm2->Fill(c2,c1,fTotGcm2);
}

//___________________________________________
void AliLego::FinishRun()
{
  //
  // Store histograms in current Root file
  //
  fHistRadl->Write();
  fHistAbso->Write();
  fHistGcm2->Write();

  // Delete histograms from memory
  fHistRadl->Delete(); fHistRadl=0;
  fHistAbso->Delete(); fHistAbso=0;
  fHistGcm2->Delete(); fHistGcm2=0;
}

//___________________________________________
void AliLego::Copy(AliLego &lego) const
{
  //
  // Copy *this onto lego -- not implemented
  //
  Fatal("Copy","Not implemented!\n");
}

//___________________________________________
void AliLego::StepManager()
{
// called from AliRun::Stepmanager from gustep.
// Accumulate the 3 parameters step by step
  
   static Float_t t;
   Float_t a,z,dens,radl,absl;
   Int_t i;
   
   Float_t step  = gMC->TrackStep();
       
   Float_t vect[3], dir[3];
   TLorentzVector pos, mom;

   gMC->CurrentMaterial(a,z,dens,radl,absl);
   
   if (z < 1) return;
    
   gMC->TrackPosition(pos);  
   gMC->TrackMomentum(mom);
// --- See if we have to stop now
   if (TMath::Abs(pos[2]) > fGener->ZMax()  || 
       pos[0]*pos[0] +pos[1]*pos[1] > fGener->RadMax()*fGener->RadMax()) {
     if (!gMC->IsNewTrack()) {
       // Not the first step, add past contribution
       fTotAbso += t/absl;
       fTotRadl += t/radl;
       fTotGcm2 += t*dens;
//       Int_t copy;
//       Int_t id = gMC->CurrentVolID(copy);
//       char* vol = gMC->VolName(id);
       
//       printf("\n %f %f %f %f %s  ", fTotRadl, vect[0], vect[1], vect[2], vol);

     }
     gMC->StopTrack();
     return;
   }

// --- See how long we have to go
   for(i=0;i<3;++i) {
     vect[i]=pos[i];
     dir[i]=mom[i];
   }

   t  = fGener->PropagateCylinder(vect,dir,fGener->RadMax(),fGener->ZMax());

   if(step) {
     fTotAbso += step/absl;
     fTotRadl += step/radl;
     fTotGcm2 += step*dens;
//     Int_t copy;
//     Int_t id = gMC->CurrentVolID(copy);
//     char* vol = gMC->VolName(id);
//     printf("\n %f %f %f %f %s %f ", fTotRadl, vect[0], vect[1], vect[2], vol, t);
   }
}





