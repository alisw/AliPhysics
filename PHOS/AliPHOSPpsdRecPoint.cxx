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
// RecPoint in the PHOS PPSD: a list of AliPHOSDigit's
//*-- Author : Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TPad.h"
// --- Standard library ---

#include <iostream>
#include <cassert>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliRun.h"

ClassImp(AliPHOSPpsdRecPoint)

//____________________________________________________________________________
AliPHOSPpsdRecPoint::AliPHOSPpsdRecPoint(void)
{ 
  fMulDigit = 0 ;  
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;  
  fDelta = geom->GetCrystalSize(0) ; 
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
AliPHOSPpsdRecPoint::~AliPHOSPpsdRecPoint(void) // dtor
{ 
  //dtor  
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::AddDigit(AliDigitNew & digit, Float_t Energy)
{
  // adds a digit to the digits list
  // and accumulates the total amplitude and the multiplicity 
  
  
  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    int * tempo = new ( int[fMaxDigit] ) ; 
    Int_t index ; 
    
    for ( index = 0 ; index < fMulDigit ; index++ )
      tempo[index] = fDigitsList[index] ;
    
    delete []  fDigitsList ; 
    fDigitsList =  new ( int[fMaxDigit] ) ;
   
    for ( index = 0 ; index < fMulDigit ; index++ )
      fDigitsList[index] = tempo[index] ;
 
    delete [] tempo ;
  }

  fDigitsList[fMulDigit++]  =  (int) &digit  ; 
  fAmp += Energy ; 
}




//____________________________________________________________________________
Int_t AliPHOSPpsdRecPoint::Compare(TObject * obj)
{
  Int_t rv ; 
  
  AliPHOSPpsdRecPoint * clu = (AliPHOSPpsdRecPoint *)obj ; 
  
 
  Float_t x1 , z1 ;
  Float_t x2 , z2 ;
  
  Int_t PhosMod1 ;
  Int_t PhosMod2 ;
  
  Int_t up1 ;
  Int_t up2 ; 
  
  if(GetUp()) // upper layer
    up1 = 0 ; 
  else        // lower layer
    up1 = 1 ;       
  
  if(clu->GetUp()) // upper layer
    up2 = 0 ; 
  else            // lower layer
    up2 = 1 ;       

  TVector3 PosLoc ;
  this->GetLocalPosition(PosLoc) ;
  x1 = PosLoc.X() ;
  z1 = PosLoc.Z() ; 
  PhosMod1 = this->GetPHOSMod();  
  clu->GetLocalPosition(PosLoc) ;
  x2 = PosLoc.X() ;
  z2 = PosLoc.Z() ; 
  PhosMod2 = clu->GetPHOSMod();

  if(PhosMod1 == PhosMod2 ) {
 
    if(up1 == up2 ){
      Int_t rowdif = (Int_t)TMath::Ceil(x1/fDelta) - (Int_t) TMath::Ceil(x2/fDelta) ;

      if (rowdif> 0) 
	rv = -1 ;
      else if(rowdif < 0) 
	rv = 1 ;
      else if(z1>z2) 
	rv = -1 ;
      else 
	rv = 1 ; 
    }

    else {

      if(up1 < up2 ) // Upper level first (up = True or False, True > False)
	rv = 1 ;   
      else 
	rv = - 1 ;
    }
    
  } // if PhosMod1 == PhosMod2

  else {

    if(PhosMod1 < PhosMod2 ) 
      rv = -1 ;
    else 
      rv = 1 ;
  
}
  
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::GetLocalPosition(TVector3 &LPos){

  if( fLocPos.X() < 1000000.) { //allready evaluated
   LPos = fLocPos ;
   return ;
  }

  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;

  AliPHOSGeometry * PHOSGeom = (AliPHOSGeometry *) fGeom ;
  
  AliPHOSDigit * digit ;
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit]; 
 
    Float_t xi ;
    Float_t zi ;
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    x += xi ;
    z += zi ;
  }

  x   /= fMulDigit ;
  z   /= fMulDigit ;

  fLocPos.SetX(x)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z)  ;

  LPos = fLocPos ;
}

//____________________________________________________________________________
Bool_t AliPHOSPpsdRecPoint::GetUp() 
{
  Int_t relid[4] ;
  
  AliPHOSGeometry * PHOSGeom = (AliPHOSGeometry *) fGeom ;
  
  AliPHOSDigit *digit = (AliPHOSDigit *)fDigitsList[0] ; 
  
  PHOSGeom->AbsToRelNumbering(digit->GetId(),relid);
  Bool_t up ;

  if((Int_t)TMath::Ceil((Float_t)relid[1]/
			(PHOSGeom->GetNumberOfModulesPhi()*PHOSGeom->GetNumberOfModulesZ())-0.0001 ) > 1) 
    up = kFALSE ;
  else  
    up = kTRUE ;
  
  return up ;
}
//______________________________________________________________________________
void AliPHOSPpsdRecPoint::Paint(Option_t *)
{
//*-*-*-*-*-*-*-*-*-*-*Paint this ALiRecPoint as a TMarker  with its current attributes*-*-*-*-*-*-*
//*-*                  =============================================
   TVector3 pos(0.,0.,0.) ;
   GetLocalPosition(pos) ;
   Coord_t x = pos.X() ;
   Coord_t y = pos.Z() ;
   Color_t MarkerColor = 1 ;
   Size_t  MarkerSize = 1. ;
   Style_t MarkerStyle = 2 ;
   if (GetUp()) 
     MarkerStyle = 3 ;

   if (!gPad->IsBatch()) {
     gVirtualX->SetMarkerColor(MarkerColor);
     gVirtualX->SetMarkerSize (MarkerSize);
     gVirtualX->SetMarkerStyle(MarkerStyle);
   }
   gPad->SetAttMarkerPS(MarkerColor,MarkerStyle,MarkerSize);
   gPad->PaintPolyMarker(1,&x,&y,"");


}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::Print(Option_t * option) 
{
  cout << "AliPHOSPpsdRecPoint: " << endl ;
  
  AliPHOSDigit * digit ; 
  Int_t iDigit;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;

  Float_t xi ;
  Float_t zi ;
  Int_t relid[4] ; 

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) fDigitsList[iDigit];
    PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    PHOSGeom->RelPosInModule(relid, xi, zi);
    cout << " Id = " << digit->GetId() ;  
    cout << "  Phos mod = " << relid[0] ;  
    cout << "  PPSD mod = " << relid[1] ;  
    cout << "  x = " << xi ;  
    cout << "  z = " << zi ;  
    cout << "   Energy = " << digit->GetAmp() << endl ;
  }
  cout << "       Multiplicity    = " << fMulDigit  << endl ;
}


// //____________________________________________________________________________
// AliPHOSPpsdRecPoint& AliPHOSPpsdRecPoint::operator = (AliPHOSPpsdRecPoint Clu)
// {
//  int* DL = Clu.GetDigitsList() ; 
  
//   if(fDigitsList) 
//     delete fDigitsList ;

//   AliPHOSDigit * digit ;
 
//   Int_t iDigit;

//   for(iDigit=0; iDigit<fMulDigit; iDigit++) {
//     digit = (AliPHOSDigit *) DL[iDigit];
//     AddDigit(*digit) ;
//   }

//   fDelta = Clu.GetDelta() ; 
//   delete DL ; 

//   fAmp       = Clu.GetEnergy() ;
//   fGeom      = Clu.GetGeom() ;
//   TVector3 LocPos;
//   Clu.GetLocalPosition(LocPos) ;
//   fLocPos    = LocPos;
//   fMulDigit  = Clu.GetMultiplicity() ;
//   fMaxDigit  = Clu.GetMaximumMultiplicity() ;
//   fPHOSMod   = Clu.GetPHOSMod() ;
  

//   return *this ;
// }
