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
// Algorithm class to identify the type of particle from the PHOS TrackSegment alone 
//*-- Author : Y. Schutz SUBATECH
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

#include <iostream>

// --- AliRoot header files ---

#include "AliPHOSPIDv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"

ClassImp( AliPHOSPIDv1) 


//____________________________________________________________________________
 AliPHOSPIDv1::AliPHOSPIDv1() 
{
  // ctor

}

//____________________________________________________________________________
 AliPHOSPIDv1::~AliPHOSPIDv1()
{ 
  // dtor
}


//____________________________________________________________________________
void  AliPHOSPIDv1::GetParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl)
{
  // main function, does the job

  TIter next(trsl) ; 
  AliPHOSTrackSegment * tracksegment ; 
  Int_t index = 0 ; 

  while ( (tracksegment = (AliPHOSTrackSegment *)next()) ) {
    new( (*rpl)[index] ) AliPHOSRecParticle(tracksegment) ; 
    index++ ; 
  }
    
}

