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
AliEMCALTrackSegment::AliEMCALTrackSegment( AliEMCALTowerRecPoint * eca, AliEMCALTowerRecPoint * pre, AliEMCALTowerRecPoint * hca)
{
  // ctor

  if( pre )   
    fPRERecPoint =  pre->GetIndexInList() ;
  else 
    fPRERecPoint = -1 ;

  if( eca )   
    fECARecPoint =  eca->GetIndexInList() ;
  else 
    fECARecPoint = -1 ;

  if( hca )   
    fHCARecPoint =  hca->GetIndexInList() ;
  else 
    fHCARecPoint = -1 ;

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
   ( (AliEMCALTrackSegment &)obj ).fPRERecPoint = fPRERecPoint ; 
   ( (AliEMCALTrackSegment &)obj ).fECARecPoint = fECARecPoint ; 
   ( (AliEMCALTrackSegment &)obj ).fHCARecPoint = fHCARecPoint ; 
   ( (AliEMCALTrackSegment &)obj ).fIndexInList = fIndexInList ; 
}

//____________________________________________________________________________
void AliEMCALTrackSegment::Print(Option_t *) const
{
  // Print all information on this track Segment
  
  
  Info("Print", "TrackSegment information:") ; 
  printf("--------AliEMCALTrackSegment-------- \n");
  printf("Stored at position %d\n", fIndexInList) ;
  if (fPRERecPoint) 
    printf("PRE RecPoint #     %d\n", fPRERecPoint) ;
  if (fECARecPoint) 
    printf("EC RecPoint  #     %d\n", fECARecPoint) ;
  if (fHCARecPoint) 
    printf("HC RecPoint  #     %d\n", fHCARecPoint) ;

  printf("------------------------------------ \n") ; 
  
}

//____________________________________________________________________________
void AliEMCALTrackSegment::SetPRERecPoint(AliEMCALRecPoint * pre) 
{
  // gives an id from its position in the list
  if( pre )  
    fPRERecPoint = pre->GetIndexInList() ;
 else 
    fPRERecPoint = -1 ;
}

//____________________________________________________________________________
void AliEMCALTrackSegment::SetHCARecPoint(AliEMCALRecPoint * hca) 
{
  // gives an id from its position in the list
  if( hca )  
    fHCARecPoint = hca->GetIndexInList() ;
 else 
    fHCARecPoint = -1 ;
}
