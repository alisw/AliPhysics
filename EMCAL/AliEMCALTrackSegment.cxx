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
//_________________________________________________________________________
//  Track segment in EMCAL
//  Can be : 1 EmcRecPoint
//           1 EmcRecPoint + 1 PPSD
//           1 EmcRecPoint + 1 PPSD + 1 PPSD     
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)
//             Adapted from PHOS by Y. Schutz (SUBATECH)

// --- ROOT system ---
 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALTrackSegment.h"

ClassImp(AliEMCALTrackSegment)

//____________________________________________________________________________
AliEMCALTrackSegment::AliEMCALTrackSegment( AliEMCALRecPoint * eca)
{
  // ctor
  if( eca )   
    fECARecPoint =  eca->GetIndexInList() ;
  else 
    fECARecPoint = -1 ;
  fIndexInList = -1 ;
}

//____________________________________________________________________________
AliEMCALTrackSegment::AliEMCALTrackSegment( const AliEMCALTrackSegment & ts) 
  : TObject(ts)
{
  // Copy ctor

  ( (AliEMCALTrackSegment &)ts ).Copy(*this) ; 
}


//____________________________________________________________________________
void AliEMCALTrackSegment::Copy(TObject & obj) 
{
  // Copy of a track segment into another track segment

   TObject::Copy(obj) ;
   ( (AliEMCALTrackSegment &)obj ).fECARecPoint = fECARecPoint ; 
   ( (AliEMCALTrackSegment &)obj ).fIndexInList = fIndexInList ; 
}

//____________________________________________________________________________
void AliEMCALTrackSegment::Print(Option_t *) const
{
  // Print all information on this track Segment
  printf("Print: TrackSegment information:") ; 
  printf("--------AliEMCALTrackSegment-------- \n");
  printf("Stored at position %d\n", fIndexInList) ;
  if (fECARecPoint) 
    printf("EC RecPoint  #     %d\n", fECARecPoint) ;
  printf("------------------------------------ \n") ;  
}
