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
//  A RecPoint (cluster) in the PPSD 
//  A PPSD RecPoint ends up to be a single digit
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
//                
//*--  Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TPad.h"

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliRun.h"

ClassImp(AliPHOSPpsdRecPoint)

//____________________________________________________________________________
AliPHOSPpsdRecPoint::AliPHOSPpsdRecPoint(void)
{ 
  // ctor

  fMulDigit = 0 ;  
  fGeom = AliPHOSGeometry::GetInstance() ;  
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::AddDigit(AliPHOSDigit & digit, Float_t Energy)
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

  fDigitsList[fMulDigit++]  =  digit.GetIndexInList() ; 
  fAmp += Energy ; 
  EvalPHOSMod(&digit) ;
}

//____________________________________________________________________________
Int_t AliPHOSPpsdRecPoint::Compare(const TObject * obj) const
{
  // Compares according to the position
  Float_t delta = 1 ; //width of the "Sorting row"

  Int_t rv ; 
  
  if( (strcmp(obj->ClassName() , "AliPHOSPpsdRecPoint" )) == 0)  // PPSD Rec Point
    {
     AliPHOSPpsdRecPoint * clu = (AliPHOSPpsdRecPoint *)obj ; 

     Float_t x1 , z1 ;    //This rec point
     Float_t x2 , z2 ;    //
     
     Int_t phosmod1 ;
     Int_t phosmod2 ;
     
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

     TVector3 posloc ;
     GetLocalPosition(posloc) ;
     x1 = posloc.X() ;
     z1 = posloc.Z() ; 
     phosmod1 = GetPHOSMod();  
     clu->GetLocalPosition(posloc) ;
     x2 = posloc.X() ;
     z2 = posloc.Z() ; 
     phosmod2 = clu->GetPHOSMod();
     
     if(phosmod1 == phosmod2 ) {
       
       if(up1 == up2 ){
	 Int_t rowdif = (Int_t)TMath::Ceil(x1/delta) - (Int_t) TMath::Ceil(x2/delta) ;
	 
	 if (rowdif> 0) 
	   rv = 1 ;
	 else if(rowdif < 0) 
	   rv = -1 ;
	 else if(z1>z2) 
	   rv = -1 ;
	 else 
	   rv = 1 ; 
       }
       
       else {
	 
	 if(up1 < up2 ) // Upper level first (up = True or False, True > False)
	   rv = -1 ;   
	 else 
	   rv = 1 ;
       }
       
     } // if phosmod1 == phosmod2
     
     else {
       
       if(phosmod1 < phosmod2 ) 
	 rv = -1 ;
       else 
	 rv = 1 ;
       
     }
     
     return rv ;      
    }
  else
    {
      AliPHOSCpvRecPoint * clu  = (AliPHOSCpvRecPoint *) obj ;   
      if(GetPHOSMod()  < clu->GetPHOSMod() ) 
	rv = -1 ;
      else 
	rv = 1 ;
      return rv ;
    }

  
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::EvalAll(Float_t logWeight,TClonesArray * digits ){
  AliPHOSRecPoint::EvalAll(logWeight,digits) ;
  EvalLocalPosition(logWeight,digits) ;
  EvalUp(digits) ;
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::EvalLocalPosition(Float_t logWeight,TClonesArray * digits )
{
  // Calculates the local position in the PHOS-PPSD-module corrdinates
  
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;

  AliPHOSGeometry * phosgeom = (AliPHOSGeometry *) fGeom ;
  
  AliPHOSDigit * digit ;
  Int_t iDigit;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ; 
 
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    x += xi ;
    z += zi ;
  }

  x /= fMulDigit ;
  z /= fMulDigit ;

  fLocPos.SetX(x)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z)  ;

}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::EvalUp(TClonesArray * digits)
{
  // Are we in the uper PPSD module ?

  Int_t relid[4] ;
  
  AliPHOSGeometry * phosgeom = (AliPHOSGeometry *) fGeom ;
  
  
  AliPHOSDigit *digit = (AliPHOSDigit *) digits->At(fDigitsList[0]) ; 
  
  phosgeom->AbsToRelNumbering(digit->GetId(),relid);

  if((Int_t)TMath::Ceil((Float_t)relid[1]/
			(phosgeom->GetNumberOfModulesPhi()*phosgeom->GetNumberOfModulesZ())-0.0001 ) > 1) 
    fUp = kFALSE ;
  else  
    fUp = kTRUE ;
  
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
  Color_t markercolor = 1 ;
  Size_t  markersize  = 1. ;
  Style_t markerstyle = 2 ;
  if (GetUp()) 
    markerstyle = 3 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor);
    gVirtualX->SetMarkerSize (markersize);
    gVirtualX->SetMarkerStyle(markerstyle);
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize);
  gPad->PaintPolyMarker(1,&x,&y,"");
  
  
}

//____________________________________________________________________________
void AliPHOSPpsdRecPoint::Print(Option_t * option) 
{
  // Print the digits information 
  
  cout << "AliPHOSPpsdRecPoint: " << endl ;
  
  Int_t iDigit; 
  cout << " Digit{s} # " ; 
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    cout  << fDigitsList[iDigit] << "  " ;  
  cout << endl   ;  

  cout << "       Multiplicity    = " << fMulDigit  << endl ;
  cout << "       Stored at position " << fIndexInList << endl ; 
}


