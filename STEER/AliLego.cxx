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
Revision 1.29  2002/10/14 14:57:32  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.28.6.2  2002/07/24 10:08:13  alibrary
Updating VirtualMC

Revision 1.28.6.1  2002/06/10 14:43:06  hristov
Merged with v3-08-02

Revision 1.28  2001/10/21 18:38:44  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.27  2001/07/20 09:32:18  morsch
Protection against uncomplete backward stepping in dumping added.

Revision 1.26  2001/05/30 12:18:13  morsch
Fastidious printf commented.

Revision 1.25  2001/05/23 11:59:46  morsch
Use RemoveAt method instead of delete to remove objects from TClonesArray.

Revision 1.24  2001/05/21 08:39:13  morsch
Use fStepBack = 1 only in debug mode.

Revision 1.23  2001/05/20 10:10:39  morsch
- Debug output at the beginning of new event and end of run.
- Filter out boundary loops.

Revision 1.22  2001/05/11 13:22:40  morsch
If run with debug option (from gAlice) geantinos are sent back and volume sequence forward/backward is compared.
Can be very verbous in some cases.

Revision 1.21  2000/12/15 10:33:59  morsch
Invert coordinates to make meaningful zylindrical plots.

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

#include "TClonesArray.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

#include "AliConst.h"
#include "AliDebugVolume.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliMC.h"
#include "AliRun.h"

ClassImp(AliLego)

//_______________________________________________________________________
AliLego::AliLego():
  fGener(0),
  fTotRadl(0),
  fTotAbso(0),
  fTotGcm2(0),
  fHistRadl(0),
  fHistAbso(0),
  fHistGcm2(0),
  fHistReta(0),
  fVolumesFwd(0),
  fVolumesBwd(0),
  fStepBack(0),
  fStepsBackward(0),
  fStepsForward(0),
  fErrorCondition(0),
  fDebug(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliLego::AliLego(const AliLego &lego):
  TNamed(lego),
  fGener(0),
  fTotRadl(0),
  fTotAbso(0),
  fTotGcm2(0),
  fHistRadl(0),
  fHistAbso(0),
  fHistGcm2(0),
  fHistReta(0),
  fVolumesFwd(0),
  fVolumesBwd(0),
  fStepBack(0),
  fStepsBackward(0),
  fStepsForward(0),
  fErrorCondition(0),
  fDebug(0)
{
  //
  // Copy constructor
  //
  lego.Copy(*this);
}


//_______________________________________________________________________
AliLego::AliLego(const char *title, Int_t ntheta, Float_t thetamin, 
                 Float_t thetamax, Int_t nphi, Float_t phimin, Float_t phimax,
                 Float_t rmin, Float_t rmax, Float_t zmax):
  TNamed("Lego Generator",title),
  fGener(0),
  fTotRadl(0),
  fTotAbso(0),
  fTotGcm2(0),
  fHistRadl(0),
  fHistAbso(0),
  fHistGcm2(0),
  fHistReta(0),
  fVolumesFwd(0),
  fVolumesBwd(0),
  fStepBack(0),
  fStepsBackward(0),
  fStepsForward(0),
  fErrorCondition(0),
  fDebug(0)
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
//
   fVolumesFwd     = new TClonesArray("AliDebugVolume",1000);
   fVolumesBwd     = new TClonesArray("AliDebugVolume",1000);
   fDebug          = gAlice->GetDebug();
}

//_______________________________________________________________________
AliLego::AliLego(const char *title, AliLegoGenerator* generator):
  TNamed("Lego Generator",title),
  fGener(0),
  fTotRadl(0),
  fTotAbso(0),
  fTotGcm2(0),
  fHistRadl(0),
  fHistAbso(0),
  fHistGcm2(0),
  fHistReta(0),
  fVolumesFwd(0),
  fVolumesBwd(0),
  fStepBack(0),
  fStepsBackward(0),
  fStepsForward(0),
  fErrorCondition(0),
  fDebug(0)
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
  //
  //

  fVolumesFwd     = new TClonesArray("AliDebugVolume",1000);
  fVolumesBwd     = new TClonesArray("AliDebugVolume",1000);
  fDebug          = gAlice->GetDebug();
}

//_______________________________________________________________________
AliLego::~AliLego()
{
  //
  // Destructor
  //
  delete fHistRadl;
  delete fHistAbso;
  delete fHistGcm2;
  delete fGener;
  delete fVolumesFwd;
  delete fVolumesBwd;
}

//_______________________________________________________________________
void AliLego::BeginEvent()
{
  //
  // --- Set to 0 radiation length, absorption length and g/cm2 ---
  //
  fTotRadl = 0;
  fTotAbso = 0;
  fTotGcm2 = 0;

  if (fDebug) {
    if (fErrorCondition) DumpVolumes();
    fVolumesFwd->Delete();
    fVolumesBwd->Delete();
    fStepsForward    = 0;
    fStepsBackward   = 0;		  
    fErrorCondition  = 0;
    if (gAlice->CurrentTrack() == 0) fStepBack = 0;
  }
}

//_______________________________________________________________________
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

//_______________________________________________________________________
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
  //
  if (fErrorCondition) DumpVolumes();
}

//_______________________________________________________________________
void AliLego::Copy(AliLego&) const
{
  //
  // Copy *this onto lego -- not implemented
  //
  Fatal("Copy","Not implemented!\n");
}

//_______________________________________________________________________
void AliLego::StepManager() 
{
  //
  // called from AliRun::Stepmanager from gustep.
  // Accumulate the 3 parameters step by step
  //
  static Float_t t;
  Float_t a,z,dens,radl,absl;
  Int_t i, id, copy;
  const char* vol;
  static Float_t vect[3], dir[3];
  
  TString tmp1, tmp2;
  copy = 1;
  id  = gMC->CurrentVolID(copy);
  vol = gMC->VolName(id);
  Float_t step  = gMC->TrackStep();
  
  TLorentzVector pos, mom; 
  gMC->TrackPosition(pos);  
  gMC->TrackMomentum(mom);
  
  Int_t status = 0;
  if (gMC->IsTrackEntering()) status = 1;
  if (gMC->IsTrackExiting())  status = 2; 
  
  if (! fStepBack) {
    //      printf("\n volume %s %d", vol, status);      
    //
    // Normal Forward stepping
    //
    if (fDebug) {
      //	  printf("\n steps fwd %d %s %d %f", fStepsForward, vol, fErrorCondition, step);	  
      
      //
      // store volume if different from previous
      //
	  
	  TClonesArray &lvols = *fVolumesFwd;
	  if (fStepsForward > 0) {
        AliDebugVolume* tmp = dynamic_cast<AliDebugVolume*>((*fVolumesFwd)[fStepsForward-1]);
        if (tmp->IsVEqual(vol, copy) && gMC->IsTrackEntering()) {
		  fStepsForward -= 2;
		  fVolumesFwd->RemoveAt(fStepsForward);
		  fVolumesFwd->RemoveAt(fStepsForward+1);		  
        }
	  }
      
	  new(lvols[fStepsForward++]) 
        AliDebugVolume(vol,copy,step,pos[0], pos[1], pos[2], status);
      
    } // Debug
    //
    // Get current material properties
    
    gMC->CurrentMaterial(a,z,dens,radl,absl);
    
    if (z < 1) return;
    
    // --- See if we have to stop now
    if (TMath::Abs(pos[2]) > fGener->ZMax()  || 
        pos[0]*pos[0] +pos[1]*pos[1] > fGener->RadMax()*fGener->RadMax()) {
      if (!gMC->IsNewTrack()) {
        // Not the first step, add past contribution
        fTotAbso += t/absl;
        fTotRadl += t/radl;
        fTotGcm2 += t*dens;
        if (fDebug) {
          //
          //  generate "mirror" particle flying back
          //
          fStepsBackward = fStepsForward;
          
          Float_t pmom[3], orig[3];
          Float_t polar[3] = {0.,0.,0.};
          Int_t ntr;
          pmom[0] = -dir[0];
          pmom[1] = -dir[1];	   
          pmom[2] = -dir[2];
          orig[0] =  vect[0];
          orig[1] =  vect[1];	   
          orig[2] =  vect[2];
          
          gAlice->SetTrack(1, gAlice->CurrentTrack(), 
                           0, pmom, orig, polar, 0., kPNoProcess, ntr);
        } // debug
        
      } // not a new track !
      
      if (fDebug) fStepBack = 1;
      gMC->StopTrack();
      return;
    } // outside scoring region ?
    
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
    }
    
  } else {
      if (fDebug) {
        //
        // Geometry debugging
        // Fly back and compare volume sequence
        //
        TClonesArray &lvols = *fVolumesBwd;
        if (fStepsBackward < fStepsForward) {
	      AliDebugVolume* tmp = dynamic_cast<AliDebugVolume*>((*fVolumesBwd)[fStepsBackward]);
	      if (tmp->IsVEqual(vol, copy) && gMC->IsTrackEntering()) {
            fStepsBackward += 2;
            fVolumesBwd->RemoveAt(fStepsBackward-1);
            fVolumesBwd->RemoveAt(fStepsBackward-2);		  
	      }
        } 
        
        fStepsBackward--;
        //	  printf("\n steps %d %s %d", fStepsBackward, vol, fErrorCondition);	  
        if (fStepsBackward < 0) {
	      gMC->StopTrack();
	      fStepBack = 0;
	      return;
        }
        
        new(lvols[fStepsBackward]) AliDebugVolume(vol,copy,step,pos[0], pos[1], pos[2], status);
        
        AliDebugVolume* tmp = dynamic_cast<AliDebugVolume*>((*fVolumesFwd)[fStepsBackward]);
        if (! (tmp->IsVEqual(vol, copy)) && (!fErrorCondition)) 
          {
            printf("\n Problem at (x,y,z): %d %f %f %f, volumes: %s %s step: %f\n", 
                   fStepsBackward, pos[0], pos[1], pos[2], tmp->GetName(), vol, step);
            fErrorCondition = 1;
          } 
      } // Debug
  } // bwd/fwd
}

//_______________________________________________________________________
void AliLego::DumpVolumes()
{
  //
  // Dump volume sequence in case of error
  //
  printf("\n Dumping Volume Sequence:");
  printf("\n ==============================================");
  
  for (Int_t i = fStepsForward-1; i>=0; i--)
    {
      AliDebugVolume* tmp1 = dynamic_cast<AliDebugVolume*>((*fVolumesFwd)[i]);
      AliDebugVolume* tmp2 = dynamic_cast<AliDebugVolume*>((*fVolumesBwd)[i]);
      if (tmp1)
        printf("\n Volume Fwd: %3d: %5s (%3d) step: %12.5e (x,y,z) (%12.5e %12.5e %12.5e) status: %9s \n"
               , i, 
               tmp1->GetName(), tmp1->CopyNumber(), tmp1->Step(), 
               tmp1->X(), tmp1->Y(), tmp1->Z(), tmp1->Status());
      if (tmp2 && i>= fStepsBackward)
        printf("\n Volume Bwd: %3d: %5s (%3d) step: %12.5e (x,y,z) (%12.5e %12.5e %12.5e) status: %9s \n"
               , i, 
               tmp2->GetName(), tmp2->CopyNumber(), tmp2->Step(), 
               tmp2->X(), tmp2->Y(), tmp2->Z(), tmp2->Status());
      
      printf("\n ............................................................................\n");
    }
  printf("\n ==============================================\n");
}
