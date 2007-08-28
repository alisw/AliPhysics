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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.7  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.6  2007/08/03 14:41:37  cvetan
 * Missing header files
 *
 * Revision 1.5  2007/08/03 13:52:16  kharlov
 * Working skeleton of matching the ESD tracks and ESD clusters (Iouri Belikov)
 *
 */

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include "AliPHOSTracker.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliESDEvent.h"
#include "AliPHOSQualAssDataMaker.h" 
#include "AliPHOSGetter.h"
#include "AliESDtrack.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"

//-------------------------------------------------------------------------
//                          PHOS tracker.
// Matches ESD tracks with the PHOS and makes the PID.  
//
//-------------------------------------------------------------------------

ClassImp(AliPHOSTracker)

Bool_t AliPHOSTracker::fgDebug = kFALSE ;  


// ***** Some geometrical constants (used in PropagateBack) 

const Double_t kR=460.+ 9;  // Radial coord. of the centre of EMC module (cm)

const Double_t kAlpha=20.*TMath::Pi()/180.;     // Segmentation angle (rad)
const Double_t kYmax=kR*TMath::Tan(0.5*kAlpha); // Maximal possible y-coord.(cm)
const Double_t kZmax=65.; // Approximately: the maximal possible z-coord.(cm)



//____________________________________________________________________________
AliPHOSTracker::AliPHOSTracker(): 
  AliTracker()
{
  //--------------------------------------------------------------------
  // The default constructor
  //--------------------------------------------------------------------
  for (Int_t i=0; i<5; i++) 
      fModules[i]=new TClonesArray("AliPHOSEmcRecPoint",777);
}

//____________________________________________________________________________
AliPHOSTracker::~AliPHOSTracker() 
{
  //--------------------------------------------------------------------
  // The destructor
  //--------------------------------------------------------------------
  for (Int_t i=0; i<5; i++) {
      (fModules[i])->Delete();
      delete fModules[i];
  }
}

//____________________________________________________________________________
Int_t AliPHOSTracker::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  // This function loads the PHOS clusters
  //--------------------------------------------------------------------
  TObjArray *arr=NULL;
  TBranch *branch=cTree->GetBranch("PHOSEmcRP");
  if (branch==0) {
    AliError("No branch with the EMC clusters found !");
    return 1;
  }
  branch->SetAddress(&arr);

  Int_t nclusters=0;
  Int_t nentr=(Int_t)branch->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!branch->GetEvent(i)) continue;
    Int_t ncl=arr->GetEntriesFast();
    while (ncl--) {
      AliPHOSEmcRecPoint *cl=(AliPHOSEmcRecPoint*)arr->UncheckedAt(ncl);

      Int_t m=cl->GetPHOSMod();
      if ((m<1)||(m>5)) {
         AliError("Wrong module index !");
         return 1;
      }

      // Here is how the alignment is treated
      if (!cl->Misalign()) AliWarning("Can't misalign this cluster !");

      cl->SetBit(14,kFALSE); // The clusters are not yet attached to any track

      TClonesArray &module=*fModules[m-1];
      Int_t idx=module.GetEntriesFast();
      new (module[idx]) AliPHOSEmcRecPoint(*cl); 

      nclusters++;

    }
  }  

  Info("LoadClusters","Number of loaded clusters: %d",nclusters);

  return 0;
}

//____________________________________________________________________________
Int_t AliPHOSTracker::PropagateBack(AliESDEvent *esd) {
  //--------------------------------------------------------------------
  // Called by AliReconstruction 
  // Performs the track matching with the PHOS modules
  // Makes the PID
  //--------------------------------------------------------------------

  Int_t nt=esd->GetNumberOfTracks();

  // *** Select and sort the ESD track in accordance with their quality
  Double_t *quality=new Double_t[nt];
  Int_t *index=new Int_t[nt];  
  for (Int_t i=0; i<nt; i++) {
     AliESDtrack *esdTrack=esd->GetTrack(i);
     quality[i] = esdTrack->GetSigmaY2() + esdTrack->GetSigmaZ2();
  }
  TMath::Sort(nt,quality,index,kFALSE);


  // *** Start the matching
  Double_t bz=GetBz(); 
  Int_t matched=0;
  for (Int_t i=0; i<nt; i++) {
     AliESDtrack *esdTrack=esd->GetTrack(index[i]);

     // Skip the tracks having "wrong" status (has to be checked/tuned)
     ULong_t status = esdTrack->GetStatus();
     if ((status & AliESDtrack::kTRDout)   == 0) continue;
     if ((status & AliESDtrack::kTRDrefit) == 1) continue;

     AliExternalTrackParam t(*esdTrack);

     Int_t isec=Int_t(t.GetAlpha()/kAlpha);
     Int_t imod=-isec-2; // PHOS module

     Double_t y;                       // Some tracks do not reach the PHOS
     if (!t.GetYAt(kR,bz,y)) continue; //    because of the bending

     Double_t z; t.GetZAt(kR,bz,z);   
     if (TMath::Abs(z) > kZmax) continue; // Some tracks miss the PHOS in Z
 
     Bool_t ok=kTRUE;
     while (TMath::Abs(y) > kYmax) {   // Find the matching module
        Double_t alp=t.GetAlpha();
        if (y > kYmax) {
	  if (!t.Rotate(alp+kAlpha)) {ok=kFALSE; break;}
          imod--;
        } else if (y < -kYmax) {
	  if (!t.Rotate(alp-kAlpha)) {ok=kFALSE; break;}
          imod++;
        }
        if (!t.GetYAt(kR,bz,y)) {ok=kFALSE; break;}
     }
     if (!ok) continue; // Track rotation failed
 
  
     if ((imod<0)||(imod>4)) continue; // Some tracks miss the PHOS in azimuth

     //t.CorrectForMaterial(...); // Correct for the TOF material, if needed
     t.PropagateTo(kR,bz);        // Propagate to the matching module


    // *** Search for the "best" cluster (can be improved)
     TClonesArray &cArray=*fModules[imod];
     Int_t ncl=cArray.GetEntriesFast();
     AliPHOSEmcRecPoint *bestCluster=0;            // The "best" cluster
     Double_t maxd2=400; // (cm^2)
     for (Int_t i=0; i<ncl; i++) {
       AliPHOSEmcRecPoint *c=(AliPHOSEmcRecPoint *)cArray.UncheckedAt(i);

       if (c->TestBit(14)) continue; // This clusters is "used"

       Double_t dy = t.GetY() - c->GetY(), dz = t.GetZ() - c->GetZ();
       Double_t d2 = dy*dy + dz*dz;
       if (d2 < maxd2) {
	  maxd2=d2;
          bestCluster=c;
       }
     }

     if (!bestCluster) continue;   // No reasonable matching found 

     bestCluster->SetBit(14,kTRUE); // This clusters is now attached to a track

     matched++;

     // *** Now, do the PID with the "bestCluster"
     // and add the corresponding info to the ESD track pointed by "esdTrack"  

     /*
     printf("%e %e %e %e\n",t.GetSign(), t.GetX() - bestCluster->GetX(),
	                                 t.GetY() - bestCluster->GetY(),
	                                 t.GetZ() - bestCluster->GetZ());
     */
  }
    
  Info("PropagateBack","Number of matched tracks: %d",matched);

  delete[] quality;
  delete[] index;

  return 0;
}

//____________________________________________________________________________
AliCluster *AliPHOSTracker::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  // Returns the pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t m=(index & 0xf0000000) >> 28;  // Module number
  Int_t i=(index & 0x0fffffff) >> 00;  // Index within the module
  
  return (AliCluster*)(fModules[m])->UncheckedAt(i);
}

//____________________________________________________________________________
void AliPHOSTracker::UnloadClusters() {
  //--------------------------------------------------------------------
  // This function unloads the PHOS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<5; i++) (fModules[i])->Delete();
}
