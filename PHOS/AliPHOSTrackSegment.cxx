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
//  Track segment in PHOS
//  Can be : 1 EmcRecPoint
//           1 EmcRecPoint + 1 PPSD
//           1 EmcRecPoint + 1 PPSD + 1 PPSD     
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---
 

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h" 

ClassImp(AliPHOSTrackSegment)

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( AliPHOSEmcRecPoint * emc , 
					  AliPHOSRecPoint * ppsdrp1,
					  AliPHOSRecPoint * ppsdrp2  ) 
{
  // ctor

  if( emc )   
    fEmcRecPoint =  emc->GetIndexInList() ;
  else 
    fEmcRecPoint = -1 ;

  if( ppsdrp1 )  
    fPpsdUpRecPoint = ppsdrp1->GetIndexInList() ;
 else 
    fPpsdUpRecPoint = -1 ;

  if( ppsdrp2  ) 
    fPpsdLowRecPoint = ppsdrp2->GetIndexInList() ;
  else 
    fPpsdLowRecPoint = -1 ;

  fIndexInList = -1 ;
}

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( const AliPHOSTrackSegment & ts) 
{
  // Copy ctor

  ( (AliPHOSTrackSegment &)ts ).Copy(*this) ; 
}


//____________________________________________________________________________
void AliPHOSTrackSegment::Copy(TObject & obj) 
{
  // Copy of a track segment into another track segment

   TObject::Copy(obj) ;
   ( (AliPHOSTrackSegment &)obj ).fEmcRecPoint     = fEmcRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdLowRecPoint = fPpsdLowRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fPpsdUpRecPoint  = fPpsdUpRecPoint ; 
   ( (AliPHOSTrackSegment &)obj ).fIndexInList     = fIndexInList ; 
}

//____________________________________________________________________________
void AliPHOSTrackSegment::Print(Option_t * opt)
{
  // Print all information on this track Segment
  

  cout << "--------AliPHOSTrackSegment-------- "<<endl ;
  cout << "Stored at position " << fIndexInList << endl ;
  cout << "Emc RecPoint #     " << fEmcRecPoint << endl ;
  if(fPpsdUpRecPoint >= 0)
    cout << "CPV RecPoint #      " << fPpsdUpRecPoint << endl ;
  else
    cout << "No CPV RecPoint " << endl ;

  if(fPpsdLowRecPoint >= 0)
    cout << "PPSD RecPoint #     " << fPpsdLowRecPoint << endl ;
  else
    cout << "No PPSD RecPoint " << endl ;
  
  cout << "------------------------------------ " << endl ; 
  
}
//____________________________________________________________________________
void AliPHOSTrackSegment::SetCpvRecPoint(AliPHOSRecPoint * PpsdUpRecPoint) 
{
  // gives an id from its position in the list
  if( PpsdUpRecPoint )  
    fPpsdUpRecPoint = PpsdUpRecPoint->GetIndexInList() ;
 else 
    fPpsdUpRecPoint = -1 ;
}

