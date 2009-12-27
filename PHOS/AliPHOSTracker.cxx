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
 * Revision 1.9  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.8  2007/08/28 12:55:08  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
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
Int_t AliPHOSTracker::LoadClusters(TTree *) {
  //--------------------------------------------------------------------
  // This function loads the PHOS clusters
  //--------------------------------------------------------------------
  return 0 ; //At this stage we can not strore result 
             // the closest track and distance to it
             //We perform same task later in AliPHOSTrackSegmentMakerv1
}

//____________________________________________________________________________
Int_t AliPHOSTracker::PropagateBack(AliESDEvent *) {
  //--------------------------------------------------------------------
  // Called by AliReconstruction 
  // Performs the track matching with the PHOS modules
  // Makes the PID
  //--------------------------------------------------------------------

  return 0 ; //At this stage we can not strore result 
             // the closest track and distance to it
             //We perform same task later in AliPHOSTrackSegmentMakerv1
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
//  for (Int_t i=0; i<5; i++) (fModules[i])->Delete();
}
