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
#include "TVirtualMC.h"

#include "AliLog.h"
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
  fDebug(0),
  fStopped(0)
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
  fDebug(0),
  fStopped(0)
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
  fDebug(0),
  fStopped(0)
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
   fDebug          = AliDebugLevel();
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
  fDebug(0),
  fStopped(0)
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
  fDebug          = AliDebugLevel();
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
  fStopped = 0;
  
  if (fDebug) {
    if (fErrorCondition) ToAliDebug(1, DumpVolumes());
    fVolumesFwd->Delete();
    fVolumesBwd->Delete();
    fStepsForward    = 0;
    fStepsBackward   = 0;		  
    fErrorCondition  = 0;
    if (gAlice->GetMCApp()->GetCurrentTrackNumber() == 0) fStepBack = 0;
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
  if (fErrorCondition) ToAliError(DumpVolumes());
}

//_______________________________________________________________________
void AliLego::Copy(TObject&) const
{
  //
  // Copy *this onto lego -- not implemented
  //
  AliFatal("Not implemented!");
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
//	    printf("\n steps fwd %d %s %d %f", fStepsForward, vol, fErrorCondition, step);	  
      
      //
      // store volume if different from previous
      //
	  
	    TClonesArray &lvols = *fVolumesFwd;
	    if (fStepsForward > 0) {
		AliDebugVolume* tmp = dynamic_cast<AliDebugVolume*>((*fVolumesFwd)[fStepsForward-1]);
		if (tmp && tmp->IsVEqual(vol, copy) && gMC->IsTrackEntering()) {
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
	if (TMath::Abs(pos[2]) > TMath::Abs(fGener->ZMax())  || 
	    pos[0]*pos[0] +pos[1]*pos[1] > fGener->RadMax()*fGener->RadMax()) {
	    if (!gMC->IsNewTrack()) {
		// Not the first step, add past contribution
		if (!fStopped) {
		    if (absl) fTotAbso += t/absl;
		    if (radl) fTotRadl += t/radl;
		    fTotGcm2 += t*dens;
		}
		
//		printf("We will stop now %5d %13.3f !\n", fStopped, t);
//		printf("%13.3f %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f %s %13.3f\n",
//		       pos[2], TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]), step, a, z, radl, absl, gMC->CurrentVolName(), fTotRadl);
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
		    
		    gAlice->GetMCApp()->PushTrack(1, gAlice->GetMCApp()->GetCurrentTrackNumber(), 
						  0, pmom, orig, polar, 0., kPNoProcess, ntr);
		} // debug
		
	    } // not a new track !
	    
	    if (fDebug) fStepBack = 1;
	    fStopped = kTRUE;
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
	    
	    if (absl) fTotAbso += step/absl;
	    if (radl) fTotRadl += step/radl;
	    fTotGcm2 += step*dens;
//	     printf("%13.3f %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f %s %13.3f\n",
//	     pos[2],  TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]), step, a, z, radl, absl, gMC->CurrentVolName(), fTotRadl);
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
		if (tmp && tmp->IsVEqual(vol, copy) && gMC->IsTrackEntering()) {
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
	    if (tmp && !(tmp->IsVEqual(vol, copy)) && (!fErrorCondition)) 
	    {
		AliWarning(Form("Problem at (x,y,z): %d %f %f %f, volumes: %s %s step: %f\n", 
				fStepsBackward, pos[0], pos[1], pos[2], tmp->GetName(), vol, step));
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
