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
// class of PHOS Sub Track
//*-- Author : Dmitri Peressounko RRC KI 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
 
#include "TVector3.h"

// --- Standard library ---

#include "iostream.h"

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h" 
#include "AliPHOSv0.h"

ClassImp(AliPHOSTrackSegment)

//____________________________________________________________________________
AliPHOSTrackSegment::AliPHOSTrackSegment( AliPHOSEmcRecPoint * emc , AliPHOSPpsdRecPoint * ppsdRP1,
                  AliPHOSPpsdRecPoint * ppsdRP2  ) 
{     
  if( emc )   
    fEmcRecPoint =  emc ;

  if( ppsdRP1 )  
    fPpsdUp = ppsdRP1 ;

  if( ppsdRP2  ) 
    fPpsdLow = ppsdRP2 ;

  fCutOnDispersion = 1.5 ; 
}

//____________________________________________________________________________
AliPHOSTrackSegment::~AliPHOSTrackSegment() // dtor
{
//    fEmcRecPoint.Delete() ;   Not Owners !!!
//    fPpsdUp.Delete() ;
//    fPpsdLow.Delete() ;
}

//____________________________________________________________________________
Float_t AliPHOSTrackSegment::GetDistanceInPHOSPlane()
{
 
  TVector3 vecEmc ;
  fEmcRecPoint->GetLocalPosition(vecEmc) ;

  TVector3 vecPpsd ;
  if( fPpsdLow->GetMultiplicity() )  
    fPpsdLow->GetLocalPosition(vecPpsd)  ; 
  else { 
    vecPpsd.SetX(10000.) ;
  } 
  vecEmc -= vecPpsd ;

  Float_t R = vecEmc.Mag();;

  return R ;
}

//____________________________________________________________________________
Bool_t AliPHOSTrackSegment::GetMomentumDirection( TVector3 & dir ) 
{   
  // True if determined
  Bool_t ifdeterm = kTRUE ;

  if( fPpsdLow ){
    TMatrix mdummy ;
    if( fPpsdUp->GetMultiplicity() ) { // draw line trough 2 points
      TVector3 posEmc ;
      fEmcRecPoint->GetGlobalPosition(posEmc,mdummy) ;
      TVector3 posPpsd ;
      fPpsdLow->GetGlobalPosition(posPpsd,mdummy) ; 
      dir = posEmc - posPpsd ;
      dir.SetMag(1.) ;
    }

    else { // draw line through 3 pionts
      TVector3 posEmc ;
      fEmcRecPoint->GetGlobalPosition(posEmc,mdummy) ;
      TVector3 posPpsdl ;
      fPpsdLow->GetGlobalPosition(posPpsdl,mdummy) ; 
      TVector3 posPpsdup ;
      fPpsdUp->GetGlobalPosition(posPpsdup,mdummy) ;
      posPpsdl = 0.5*( posPpsdup+posPpsdl ) ; 
      dir = posEmc - posPpsdl ;
      dir.SetMag(1.) ;
    }
  
  }
  else 
    ifdeterm = kFALSE ;
 
  return ifdeterm ;
}

//____________________________________________________________________________
Int_t AliPHOSTrackSegment::GetPartType()  
{  
  // Returns 0 - gamma
  //         1 - e+, e-
  //         2 - neutral hadron  
  //         3 - charged hadron

  Int_t PartType =0;                            
  if( fPpsdUp ){     // Neutral

    if( fPpsdLow ) // Neutral hadron  
      PartType = 2 ;   
    else                                // Gamma
      PartType = 0 ;               

  }

  else {             // Charged           

    if( fEmcRecPoint->GetDispersion() > fCutOnDispersion) 
      PartType = 3 ;
    else  
      PartType = 1 ;  
  
  }
  
  return   PartType ;                     

}

//____________________________________________________________________________
void AliPHOSTrackSegment::GetPosition( TVector3 & pos ) 
{  
  // Returns positions of hits
  TMatrix Dummy ;
  fEmcRecPoint->GetGlobalPosition(pos, Dummy) ;
}

//____________________________________________________________________________
void AliPHOSTrackSegment::Print()
{
  cout << "--------AliPHOSTrackSegment-------- "<<endl ;
  cout << "EMC Reconstructed Point: " << fEmcRecPoint << endl;
  
  TVector3 pos ;
  TMatrix Dummy ;  

  fEmcRecPoint->GetGlobalPosition( pos, Dummy ) ;
 
  cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << "      Energy " << fEmcRecPoint->GetTotalEnergy() << endl ;
  cout << "PPSD Low Reconstructed Point: " << endl;
  
  if(fPpsdLow){
    fPpsdLow->GetGlobalPosition( pos , Dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z() << endl ;
  }

  cout << "PPSD Up Reconstructed Point: " << endl;
  
  if(fPpsdUp ){
    fPpsdUp->GetGlobalPosition( pos, Dummy ) ;
    cout << "    position " << pos.X() << "   " << pos.Y() << "  " << pos.Z()  << endl ;
  }

}

