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
// Geometry class for PHOS version SUBATECH
//*-- Author : Y. Schutz SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TVector3.h"
#include "TRotation.h" 

// --- Standard library ---

#include <iostream.h>
#include "assert.h"

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliConst.h"

ClassImp(AliPHOSGeometry)

  AliPHOSGeometry * AliPHOSGeometry::fGeom = 0 ;

//____________________________________________________________________________
AliPHOSGeometry::~AliPHOSGeometry(void)
{
  fRotMatrixArray->Delete() ; 
  delete fRotMatrixArray ; 
}

//____________________________________________________________________________
Bool_t AliPHOSGeometry::AbsToRelNumbering(const Int_t AbsId, Int_t * RelId)
{
  // RelId[0] = PHOS Module number 1:fNModules 
  // RelId[1] = 0 if PbW04
  //          = PPSD Module number 1:fNumberOfModulesPhi*fNumberOfModulesZ*2 (2->up and bottom level)
  // RelId[2] = Row number inside a PHOS or PPSD module
  // RelId[3] = Column number inside a PHOS or PPSD module

  Bool_t rv  = kTRUE ; 
  Float_t Id = AbsId ;

  Int_t PHOSModuleNumber = (Int_t)TMath:: Ceil( Id / ( GetNPhi() * GetNZ() ) ) ; 
  
  if ( PHOSModuleNumber >  GetNModules() ) { // its a PPSD pad

    Id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
    Float_t tempo = 2 *  GetNumberOfModulesPhi() * GetNumberOfModulesZ() *  GetNumberOfPadsPhi() * GetNumberOfPadsZ() ; 
    RelId[0] = (Int_t)TMath::Ceil( Id / tempo ) ; 
    Id -= ( RelId[0] - 1 ) * tempo ;
    RelId[1] = (Int_t)TMath::Ceil( Id / ( GetNumberOfPadsPhi() * GetNumberOfPadsZ() ) ) ; 
    Id -= ( RelId[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ() ;
    RelId[2] = (Int_t)TMath::Ceil( Id / GetNumberOfPadsPhi() ) ;
    RelId[3] = (Int_t) ( Id - ( RelId[2] - 1 )  * GetNumberOfPadsPhi() ) ; 
  } 
  else { // its a PW04 crystal

    RelId[0] = PHOSModuleNumber ;
    RelId[1] = 0 ;
    Id -= ( PHOSModuleNumber - 1 ) *  GetNPhi() * GetNZ() ; 
    RelId[2] = (Int_t)TMath::Ceil( Id / GetNPhi() ) ;
    RelId[3] = (Int_t)( Id - ( RelId[2] - 1 ) * GetNPhi() ) ; 
  } 
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos, TMatrix & gmat)
{

  AliPHOSRecPoint * tmpPHOS = (AliPHOSRecPoint *) RecPoint ;  
  TVector3 LocalPosition ;

  tmpPHOS->GetLocalPosition(gpos) ;


  if ( tmpPHOS->IsEmc() ) // it is a EMC crystal 
    {  gpos.SetY( -(GetIPtoOuterCoverDistance() + GetUpperPlateThickness() +
		    GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()) ) ;  

    }
  else
    { // it is a PPSD pad
      AliPHOSPpsdRecPoint * tmpPpsd = (AliPHOSPpsdRecPoint *) RecPoint ;
      if (tmpPpsd->GetUp() ) // it is an upper module
	{
	  gpos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() - 
		       GetLeadToMicro2Gap() - GetLeadConverterThickness() -  
		       GetMicro1ToLeadGap() - GetMicromegas1Thickness() / 2.0 )  ) ; 
	} 
      else // it is a lower module
	gpos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() / 2.0) ) ; 
    }  

  Float_t Phi           = GetPHOSAngle( tmpPHOS->GetPHOSMod()) ; 
  Double_t const RADDEG = 180.0 / kPI ;
  Float_t rPhi          = Phi / RADDEG ; 
  
  TRotation Rot ;
  Rot.RotateZ(-rPhi) ; // a rotation around Z by angle  
  
  TRotation dummy = Rot.Invert() ;  // to transform from original frame to rotate frame
  gpos.Transform(Rot) ; // rotate the baby 
}

//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos)
{
  AliPHOSRecPoint * tmpPHOS = (AliPHOSRecPoint *) RecPoint ;  
  TVector3 LocalPosition ;
  tmpPHOS->GetLocalPosition(gpos) ;


  if ( tmpPHOS->IsEmc() ) // it is a EMC crystal 
    {  gpos.SetY( -(GetIPtoOuterCoverDistance() + GetUpperPlateThickness() +
		    GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()) ) ;  
    }
  else
    { // it is a PPSD pad
      AliPHOSPpsdRecPoint * tmpPpsd = (AliPHOSPpsdRecPoint *) RecPoint ;
      if (tmpPpsd->GetUp() ) // it is an upper module
	{
	  gpos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() - 
		       GetLeadToMicro2Gap() - GetLeadConverterThickness() -  
		       GetMicro1ToLeadGap() - GetMicromegas1Thickness() / 2.0 )  ) ; 
	} 
      else // it is a lower module
	gpos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() / 2.0) ) ; 
    }  

  Float_t Phi           = GetPHOSAngle( tmpPHOS->GetPHOSMod()) ; 
  Double_t const RADDEG = 180.0 / kPI ;
  Float_t rPhi          = Phi / RADDEG ; 
  
  TRotation Rot ;
  Rot.RotateZ(-rPhi) ; // a rotation around Z by angle  
  
  TRotation dummy = Rot.Invert() ;  // to transform from original frame to rotate frame
  gpos.Transform(Rot) ; // rotate the baby 
}

//____________________________________________________________________________
void AliPHOSGeometry::Init(void)
{
  fRotMatrixArray = new TObjArray(fNModules) ; 

  cout << "PHOS geometry setup: parameters for option " << fName << " " << fTitle << endl ;
  if ( ((strcmp( fName, "default" )) == 0)  || ((strcmp( fName, "GPS2" )) == 0) ) {
    fInit     = kTRUE ; 
    this->InitPHOS() ; 
    this->InitPPSD() ;
    this->SetPHOSAngles() ; 
  }
 else {
   fInit = kFALSE ; 
   cout << "PHOS Geometry setup: option not defined " << fName << endl ; 
 }
}

//____________________________________________________________________________
void AliPHOSGeometry::InitPHOS(void)
{
     // PHOS 

  fNPhi     = 64 ; 
  fNZ       = 64 ; 
  fNModules =  5 ; 
  
  fPHOSAngle[0] = 0.0 ; // Module position angles are set in CreateGeometry()
  fPHOSAngle[1] = 0.0 ;
  fPHOSAngle[2] = 0.0 ;
  fPHOSAngle[3] = 0.0 ;
 
  fXtlSize[0] =  2.2 ;
  fXtlSize[1] = 18.0 ;
  fXtlSize[2] =  2.2 ;

  // all these numbers coming next are subject to changes

  fOuterBoxThickness[0] = 2.8 ;
  fOuterBoxThickness[1] = 5.0 ;      
  fOuterBoxThickness[2] = 5.0 ;
  
  fUpperPlateThickness  = 4.0 ;
  
  fSecondUpperPlateThickness = 5.0 ; 
  
  fCrystalSupportHeight   = 6.95 ; 
  fCrystalWrapThickness   = 0.01 ;
  fCrystalHolderThickness = 0.005 ;
  fModuleBoxThickness     = 2.0 ; 
  fIPtoOuterCoverDistance = 447.0 ;      
  fIPtoCrystalSurface     = 460.0 ;  
  
  fPinDiodeSize[0] = 1.0 ;    
  fPinDiodeSize[1] = 0.1 ;    
  fPinDiodeSize[2] = 1.0 ;    
  
  fUpperCoolingPlateThickness   = 0.06 ; 
  fSupportPlateThickness        = 10.0 ;
  fLowerThermoPlateThickness    =  3.0 ; 
  fLowerTextolitPlateThickness  =  1.0 ;
  fGapBetweenCrystals           = 0.03 ;
  
  fTextolitBoxThickness[0] = 1.5 ;  
  fTextolitBoxThickness[1] = 0.0 ;   
  fTextolitBoxThickness[2] = 3.0 ; 
  
  fAirThickness[0] =  1.56   ;
  fAirThickness[1] = 20.5175 ;  
  fAirThickness[2] =  2.48   ;  
  
  Float_t XtalModulePhiSize =  fNPhi * ( fXtlSize[0] + 2 * fGapBetweenCrystals ) ; 
  Float_t XtalModuleZSize   =  fNZ * ( fXtlSize[2] + 2 * fGapBetweenCrystals ) ;
  
  // The next dimensions are calculated from the above parameters
  
  fOuterBoxSize[0] =  XtalModulePhiSize + 2 * ( fAirThickness[0] + fModuleBoxThickness
						+ fTextolitBoxThickness[0] + fOuterBoxThickness[0] ) ; 
  fOuterBoxSize[1] = ( fXtlSize[1] + fCrystalSupportHeight + fCrystalWrapThickness + fCrystalHolderThickness )
    + 2 * (fAirThickness[1] +  fModuleBoxThickness + fTextolitBoxThickness[1] + fOuterBoxThickness[1] ) ;
  fOuterBoxSize[2] =  XtalModuleZSize +  2 * ( fAirThickness[2] + fModuleBoxThickness 
					       + fTextolitBoxThickness[2] + fOuterBoxThickness[2] ) ; 
  
  fTextolitBoxSize[0]  = fOuterBoxSize[0] - 2 * fOuterBoxThickness[0] ;
  fTextolitBoxSize[1]  = fOuterBoxSize[1] -  fOuterBoxThickness[1] - fUpperPlateThickness ;
  fTextolitBoxSize[2]  = fOuterBoxSize[2] - 2 * fOuterBoxThickness[2] ;
  
  fAirFilledBoxSize[0] =  fTextolitBoxSize[0] - 2 * fTextolitBoxThickness[0] ; 
  fAirFilledBoxSize[1] =  fTextolitBoxSize[1] - fSecondUpperPlateThickness ; 
  fAirFilledBoxSize[2] =  fTextolitBoxSize[2] - 2 * fTextolitBoxThickness[2] ; 
  
}

//____________________________________________________________________________
void AliPHOSGeometry::InitPPSD(void)
{
    // PPSD
    
  fAnodeThickness           = 0.0009 ; 
  fAvalancheGap             = 0.01 ; 
  fCathodeThickness         = 0.0009 ;
  fCompositeThickness       = 0.3 ; 
  fConversionGap            = 0.3 ; 
  fLeadConverterThickness   = 0.56 ; 
  fLeadToMicro2Gap          = 0.1 ; 
  fLidThickness             = 0.2 ; 
  fMicro1ToLeadGap          = 0.1 ; 
  fMicromegasWallThickness  = 0.6 ; 
  fNumberOfModulesPhi       = 4 ; 
  fNumberOfModulesZ         = 4 ; 
  fNumberOfPadsPhi          = 24 ; 
  fNumberOfPadsZ            = 24 ;   
  fPCThickness              = 0.1 ; 
  fPhiDisplacement          = 0.8 ;  
  fZDisplacement            = 0.8 ;  

  fMicromegas1Thickness   = fLidThickness + 2 * fCompositeThickness + fCathodeThickness + fPCThickness 
                              + fAnodeThickness + fConversionGap + fAvalancheGap ; 
  fMicromegas2Thickness   = fMicromegas1Thickness ; 


  fPPSDModuleSize[0] = 38.0 ; 
  fPPSDModuleSize[1] = fMicromegas1Thickness ; 
  fPPSDModuleSize[2] = 38.0 ; 
 
  fPPSDBoxSize[0] = fNumberOfModulesPhi * fPPSDModuleSize[0] + 2 * fPhiDisplacement ;  
  fPPSDBoxSize[1] = fMicromegas2Thickness + fMicromegas2Thickness + fLeadConverterThickness + fMicro1ToLeadGap + fLeadToMicro2Gap ;    
  fPPSDBoxSize[2] = fNumberOfModulesZ *  fPPSDModuleSize[2] + 2 * fZDisplacement ;

  fIPtoTopLidDistance     = fIPtoOuterCoverDistance -  fPPSDBoxSize[1] - 1. ;  
  
}

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance() 
{ 
  assert(fGeom!=0) ; 
  return (AliPHOSGeometry *) fGeom ; 
}

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance(const Text_t* name, const Text_t* title) 
{
  AliPHOSGeometry * rv = 0  ; 
  if ( fGeom == 0 ) {
    fGeom = new AliPHOSGeometry(name, title) ; 
    rv = (AliPHOSGeometry * ) fGeom ; 
  }
  else {
    if ( strcmp(fGeom->GetName(), name) != 0 ) {
      cout << "AliPHOSGeometry <E> : current geometry is " << fGeom->GetName() << endl
	   << "                      you cannot call     " << name << endl ; 
    }
    else
      rv = (AliPHOSGeometry *) fGeom ; 
  } 
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSGeometry::RelToAbsNumbering(const Int_t * RelId, Int_t &  AbsId)
{

  // AbsId = 1:fNModules * fNPhi * fNZ  -> PbWO4
  // AbsId = 1:fNModules * 2 * (fNumberOfModulesPhi * fNumberOfModulesZ) * fNumberOfPadsPhi * fNumberOfPadsZ -> PPSD

  Bool_t rv = kTRUE ; 
 
  if ( RelId[1] > 0 ) { // its a PPSD pad

    AbsId =    GetNPhi() * GetNZ() *  GetNModules()                          // the offset to separate emcal crystals from PPSD pads
      + ( RelId[0] - 1 ) * GetNumberOfModulesPhi() * GetNumberOfModulesZ()   // the pads offset of PHOS modules 
                         * GetNumberOfPadsPhi() * GetNumberOfPadsZ() * 2
      + ( RelId[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ()         // the pads offset of PPSD modules 
      + ( RelId[2] - 1 ) * GetNumberOfPadsPhi()                              // the pads offset of a PPSD row
      + RelId[3] ;                                                           // the column number
  } 
  else {
    if ( RelId[1] == 0 ) { // its a Phos crystal
      AbsId =  ( RelId[0] - 1 ) *  GetNPhi() * GetNZ() // the offset of PHOS modules
        + ( RelId[2] - 1 ) * GetNPhi()                 // the offset of a xtal row
        + RelId[3] ;                                   // the column number
    }
  }

  return rv ; 
}

//____________________________________________________________________________

void AliPHOSGeometry::RelPosInAlice(const Int_t Id, TVector3 & pos ) 
{
   if (Id > 0) { 

  Int_t RelId[4] ;
 
  AbsToRelNumbering(Id , RelId) ;

  Int_t PHOSModule = RelId[0] ; 

  
  if ( RelId[1] == 0 ) // it is a PbW04 crystal 
  {  pos.SetY( -(GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
      + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()) ) ;  
  }
  if ( RelId[1] > 0 ) { // its a PPSD pad
    if ( RelId[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() ) // its an bottom module
     {
       pos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() / 2.0) ) ;
     } 
    else // its an upper module
      pos.SetY(-( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() - GetLeadToMicro2Gap()
	-  GetLeadConverterThickness() -  GetMicro1ToLeadGap() - GetMicromegas1Thickness() / 2.0) ) ; 
  }

  Float_t x, z ; 
  RelPosInModule(RelId, x, z) ; 

  pos.SetX(x);
  pos.SetZ(z);


   Float_t Phi           = GetPHOSAngle( PHOSModule) ; 
   Double_t const RADDEG = 180.0 / kPI ;
   Float_t rPhi          = Phi / RADDEG ; 

   TRotation Rot ;
   Rot.RotateZ(-rPhi) ; // a rotation around Z by angle  
  
   TRotation dummy = Rot.Invert() ;  // to transform from original frame to rotate frame
  
   pos.Transform(Rot) ; // rotate the baby 
  }
  else {
 pos.SetX(0.);
 pos.SetY(0.);
 pos.SetZ(0.);
       }
} 

//____________________________________________________________________________
void AliPHOSGeometry::RelPosInModule(const Int_t * RelId, Float_t & x, Float_t & z) 
{
  Int_t PPSDModule  ; 
  Int_t Row        = RelId[2] ; //offset along z axiz
  Int_t Column     = RelId[3] ; //offset along x axiz

  Float_t PadSizeZ = GetPPSDModuleSize(2)/ GetNumberOfPadsZ();
  Float_t PadSizeX = GetPPSDModuleSize(0)/ GetNumberOfPadsPhi();

  if ( RelId[1] == 0 ) { // its a PbW04 crystal 
    x = -( GetNPhi()/2. - Row   + 0.5 ) *  GetCrystalSize(0) ; // position ox Xtal with respect
    z = -( GetNZ() /2. - Column + 0.5 ) *  GetCrystalSize(2) ; // of center of PHOS module  
   }  
   else  {    
    if ( RelId[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() )
       PPSDModule =  RelId[1]-GetNumberOfModulesPhi() *  GetNumberOfModulesZ(); 
    else PPSDModule =  RelId[1] ;
    Int_t ModRow = 1+(Int_t)TMath::Ceil( (Float_t)PPSDModule / GetNumberOfModulesPhi()-1. ) ; 
    Int_t ModCol = PPSDModule -  ( ModRow-1 ) * GetNumberOfModulesPhi() ;     
    Float_t x0 = (  GetNumberOfModulesPhi() / 2.  - ModRow  + 0.5 ) * GetPPSDModuleSize(0) ;
    Float_t z0 = (  GetNumberOfModulesZ() / 2.  - ModCol  + 0.5 ) * GetPPSDModuleSize(2)  ;     
    x = - ( GetNumberOfPadsPhi()/2. - Row - 0.5 ) * PadSizeX + x0 ; // position of pad  with respect
    z = - ( GetNumberOfPadsZ()/2.   - Column - 0.5 ) * PadSizeZ + z0 ; // of center of PHOS module  
         }
}

//____________________________________________________________________________
void AliPHOSGeometry:: SetPHOSAngles() 
{ 
  Double_t const RADDEG = 180.0 / kPI ;
  Float_t PPHI =  TMath::ATan( fOuterBoxSize[0]  / ( 2.0 * fIPtoOuterCoverDistance ) ) ;
  PPHI *= RADDEG ;
  
  for( Int_t i = 1; i <= fNModules ; i++ ) {
    Float_t angle = PPHI * 2 * ( i - fNModules / 2.0 - 0.5 ) ;
    fPHOSAngle[i-1] = -  angle ;
 } 
}

