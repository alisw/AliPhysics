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
// A brief description of the class
//*-- Author : Gines MARTINEZ  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSReconstructioner.h"
#include "AliPHOSClusterizer.h"

ClassImp(AliPHOSReconstructioner)


//____________________________________________________________________________
AliPHOSReconstructioner::AliPHOSReconstructioner() 
{
  // ctor
}        

//____________________________________________________________________________
AliPHOSReconstructioner::AliPHOSReconstructioner(AliPHOSClusterizer & Clusterizer, AliPHOSTrackSegmentMaker & Tracker)
{
  fClusterizer         = &Clusterizer ;
  fTrackSegmentMaker  = &Tracker ;
} 

//____________________________________________________________________________
AliPHOSReconstructioner::~AliPHOSReconstructioner() 
{
  // dtor
}  

//____________________________________________________________________________
 void AliPHOSReconstructioner:: Make(TClonesArray * dl, RecPointsList * emccl, RecPointsList * ppsdl, TrackSegmentsList * trsl)
{
  fClusterizer->MakeClusters(dl, emccl, ppsdl);

  fTrackSegmentMaker->MakeTrackSegments(dl,emccl,ppsdl,trsl) ;
}
