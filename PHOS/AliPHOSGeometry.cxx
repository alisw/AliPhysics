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

#include <iostream>
#include <cassert>

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
Bool_t AliPHOSGeometry::AbsToRelNumbering(const Int_t AbsId, Int_t * relid)
{
  // relid[0] = PHOS Module number 1:fNModules 
  // relid[1] = 0 if PbW04
  //          = PPSD Module number 1:fNumberOfModulesPhi*fNumberOfModulesZ*2 (2->up and bottom level)
  // relid[2] = Row number inside a PHOS or PPSD module
  // relid[3] = Column number inside a PHOS or PPSD module

  Bool_t rv  = kTRUE ; 
  Float_t id = AbsId ;

  Int_t phosmodulenumber = (Int_t)TMath:: Ceil( id / ( GetNPhi() * GetNZ() ) ) ; 
  
  if ( phosmodulenumber >  GetNModules() ) { // its a PPSD pad

    id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
    Float_t tempo = 2 *  GetNumberOfModulesPhi() * GetNumberOfModulesZ() *  GetNumberOfPadsPhi() * GetNumberOfPadsZ() ; 
    relid[0] = (Int_t)TMath::Ceil( id / tempo ) ; 
    id -= ( relid[0] - 1 ) * tempo ;
    relid[1] = (Int_t)TMath::Ceil( id / ( GetNumberOfPadsPhi() * GetNumberOfPadsZ() ) ) ; 
    id -= ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ() ;
    relid[2] = (Int_t)TMath::Ceil( id / GetNumberOfPadsPhi() ) ;
    relid[3] = (Int_t) ( id - ( relid[2] - 1 )  * GetNumberOfPadsPhi() ) ; 
  } 
  else { // its a PW04 crystal

    relid[0] = phosmodulenumber ;
    relid[1] = 0 ;
    id -= ( phosmodulenumber - 1 ) *  GetNPhi() * GetNZ() ; 
    relid[2] = (Int_t)TMath::Ceil( id / GetNPhi() ) ;
    relid[3] = (Int_t)( id - ( relid[2] - 1 ) * GetNPhi() ) ; 
  } 
  return rv ; 
}
//____________________________________________________________________________  
void AliPHOSGeometry::EmcModuleCoverage(const Int_t mod, Double_t & tm, Double_t & tM, Double_t & pm, Double_t & pM, Option_t * opt) 
{
  // calculates the angular coverage in theta and phi of a EMC module

 Double_t conv ; 
  if ( opt == kRadian ) 
    conv = 1. ; 
  else if ( opt == kDegre )
    conv = 180. / TMath::Pi() ; 
  else {
    cout << "<I>  AliPHOSGeometry::EmcXtalCoverage : " << opt << " unknown option; result in radian " << endl ; 
    conv = 1. ;
      }

  Float_t phi =  GetPHOSAngle(mod) *  (TMath::Pi() / 180.)  ;  
  Float_t y0  =  GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
		  + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()  ;  
  
  Double_t angle = TMath::ATan( GetCrystalSize(0)*GetNPhi() / (2 * y0) ) ;
  phi = phi + 1.5 * TMath::Pi() ; // to follow the convention of the particle generator(PHOS is between 230 and 310 deg.)
  Double_t max  = phi - angle ;
  Double_t min   = phi + angle ;
  pM = TMath::Max(max, min) * conv ;
  pm = TMath::Min(max, min) * conv ; 
  
  angle =  TMath::ATan( GetCrystalSize(2)*GetNZ() / (2 * y0) ) ;
  max  = TMath::Pi() / 2.  + angle ; // to follow the convention of the particle generator(PHOS is at 90 deg.)
  min  = TMath::Pi() / 2.  - angle ;
  tM = TMath::Max(max, min) * conv ;
  tm = TMath::Min(max, min) * conv ; 
 
}

//____________________________________________________________________________  
void AliPHOSGeometry::EmcXtalCoverage(Double_t & theta, Double_t & phi, Option_t * opt) 
{
  // calculates the angular coverage in theta and phi of a single crystal in a EMC module

  Double_t conv ; 
  if ( opt == kRadian ) 
    conv = 1. ; 
  else if ( opt == kDegre )
    conv = 180. / TMath::Pi() ; 
  else {
    cout << "<I>  AliPHOSGeometry::EmcXtalCoverage : " << opt << " unknown option; result in radian " << endl ; 
    conv = 1. ;
      }

  Float_t y0   =  GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
    + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()  ;  
  theta = 2 * TMath::ATan( GetCrystalSize(2) / (2 * y0) ) * conv ;
  phi   = 2 * TMath::ATan( GetCrystalSize(0) / (2 * y0) ) * conv ;
}
 

//____________________________________________________________________________
void AliPHOSGeometry::ImpactOnEmc(const Double_t theta, const Double_t phi, Int_t & ModuleNumber, Double_t & z, Double_t & x) 
{
  // calculates the impact coordinates of a neutral particle  
  // emitted in direction theta and phi in ALICE

  // searches for the PHOS EMC module
  ModuleNumber = 0 ; 
  Double_t tm, tM, pm, pM ; 
  Int_t index = 1 ; 
  while ( ModuleNumber == 0 && index <= GetNModules() ) { 
    EmcModuleCoverage(index, tm, tM, pm, pM) ; 
    if ( (theta >= tm && theta <= tM) && (phi >= pm && phi <= pM ) ) 
      ModuleNumber = index ; 
    index++ ;    
  }
  if ( ModuleNumber != 0 ) {
    Float_t phi0 =  GetPHOSAngle(ModuleNumber) *  (TMath::Pi() / 180.) + 1.5 * TMath::Pi()  ;  
    Float_t y0  =  GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
      + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness()  ;   
    Double_t angle = phi - phi0; 
    x = y0 * TMath::Tan(angle) ; 
    angle = theta - TMath::Pi() / 2 ; 
    z = y0 * TMath::Tan(angle) ; 
  }
}

//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos, TMatrix & gmat)
{

  AliPHOSRecPoint * tmpPHOS = (AliPHOSRecPoint *) RecPoint ;  
  TVector3 localposition ;

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

  Float_t phi           = GetPHOSAngle( tmpPHOS->GetPHOSMod()) ; 
  Double_t const kRADDEG = 180.0 / kPI ;
  Float_t rphi          = phi / kRADDEG ; 
  
  TRotation rot ;
  rot.RotateZ(-rphi) ; // a rotation around Z by angle  
  
  TRotation dummy = rot.Invert() ;  // to transform from original frame to rotate frame
  gpos.Transform(rot) ; // rotate the baby 

}

//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos)
{
  AliPHOSRecPoint * tmpPHOS = (AliPHOSRecPoint *) RecPoint ;  
  TVector3 localposition ;
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

  Float_t phi           = GetPHOSAngle( tmpPHOS->GetPHOSMod()) ; 
  Double_t const kRADDEG = 180.0 / kPI ;
  Float_t rphi          = phi / kRADDEG ; 
  
  TRotation rot ;
  rot.RotateZ(-rphi) ; // a rotation around Z by angle  
  
  TRotation dummy = rot.Invert() ;  // to transform from original frame to rotate frame
  gpos.Transform(rot) ; // rotate the baby 
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
  
  fPinDiodeSize[0] = 1.71 ;   //Values given by Odd Harald feb 2000  
  fPinDiodeSize[1] = 0.0280 ; // 0.0280 is the depth of active layer in the silicon     
  fPinDiodeSize[2] = 1.61 ;    
  
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
  
  Float_t xtalModulePhiSize =  fNPhi * ( fXtlSize[0] + 2 * fGapBetweenCrystals ) ; 
  Float_t xtalModuleZSize   =  fNZ * ( fXtlSize[2] + 2 * fGapBetweenCrystals ) ;
  
  // The next dimensions are calculated from the above parameters
  
  fOuterBoxSize[0] =  xtalModulePhiSize + 2 * ( fAirThickness[0] + fModuleBoxThickness
						+ fTextolitBoxThickness[0] + fOuterBoxThickness[0] ) ; 
  fOuterBoxSize[1] = ( fXtlSize[1] + fCrystalSupportHeight + fCrystalWrapThickness + fCrystalHolderThickness )
    + 2 * (fAirThickness[1] +  fModuleBoxThickness + fTextolitBoxThickness[1] + fOuterBoxThickness[1] ) ;
  fOuterBoxSize[2] =  xtalModuleZSize +  2 * ( fAirThickness[2] + fModuleBoxThickness 
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
  fConversionGap            = 0.6 ; 
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
Bool_t AliPHOSGeometry::RelToAbsNumbering(const Int_t * relid, Int_t &  AbsId)
{

  // AbsId = 1:fNModules * fNPhi * fNZ  -> PbWO4
  // AbsId = 1:fNModules * 2 * (fNumberOfModulesPhi * fNumberOfModulesZ) * fNumberOfPadsPhi * fNumberOfPadsZ -> PPSD

  Bool_t rv = kTRUE ; 
 
  if ( relid[1] > 0 ) { // its a PPSD pad

    AbsId =    GetNPhi() * GetNZ() *  GetNModules()                          // the offset to separate emcal crystals from PPSD pads
      + ( relid[0] - 1 ) * GetNumberOfModulesPhi() * GetNumberOfModulesZ()   // the pads offset of PHOS modules 
                         * GetNumberOfPadsPhi() * GetNumberOfPadsZ() * 2
      + ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ()         // the pads offset of PPSD modules 
      + ( relid[2] - 1 ) * GetNumberOfPadsPhi()                              // the pads offset of a PPSD row
      + relid[3] ;                                                           // the column number
  } 
  else {
    if ( relid[1] == 0 ) { // its a Phos crystal
      AbsId =  ( relid[0] - 1 ) *  GetNPhi() * GetNZ() // the offset of PHOS modules
        + ( relid[2] - 1 ) * GetNPhi()                 // the offset of a xtal row
        + relid[3] ;                                   // the column number
    }
  }

  return rv ; 
}

//____________________________________________________________________________

void AliPHOSGeometry::RelPosInAlice(const Int_t id, TVector3 & pos ) 
{
   if (id > 0) { 

  Int_t relid[4] ;
 
  AbsToRelNumbering(id , relid) ;

  Int_t phosmodule = relid[0] ; 

  Float_t y0 = 0 ; 

  if ( relid[1] == 0 ) // it is a PbW04 crystal 
  {  y0 =  -(GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
      + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness())  ;  
  }
  if ( relid[1] > 0 ) { // its a PPSD pad
    if ( relid[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() ) // its an bottom module
     {
       y0 = -( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() / 2.0)  ;
     } 
    else // its an upper module
      y0 = -( GetIPtoOuterCoverDistance() - GetMicromegas2Thickness() - GetLeadToMicro2Gap()
	-  GetLeadConverterThickness() -  GetMicro1ToLeadGap() - GetMicromegas1Thickness() / 2.0) ; 
  }

  Float_t x, z ; 
  RelPosInModule(relid, x, z) ; 

  pos.SetX(x) ;
  pos.SetZ(z) ;
  pos.SetY( TMath::Sqrt(x*x + z*z + y0*y0) ) ; 



   Float_t phi           = GetPHOSAngle( phosmodule) ; 
   Double_t const kRADDEG = 180.0 / kPI ;
   Float_t rphi          = phi / kRADDEG ; 

   TRotation rot ;
   rot.RotateZ(-rphi) ; // a rotation around Z by angle  
  
   TRotation dummy = rot.Invert() ;  // to transform from original frame to rotate frame
  
   pos.Transform(rot) ; // rotate the baby 
  }
  else {
 pos.SetX(0.);
 pos.SetY(0.);
 pos.SetZ(0.);
       }
} 

//____________________________________________________________________________
void AliPHOSGeometry::RelPosInModule(const Int_t * relid, Float_t & x, Float_t & z) 
{
  Int_t ppsdmodule  ; 
  Int_t row        = relid[2] ; //offset along z axiz
  Int_t column     = relid[3] ; //offset along x axiz

  Float_t padsizeZ = GetPPSDModuleSize(2)/ GetNumberOfPadsZ();
  Float_t padsizeX = GetPPSDModuleSize(0)/ GetNumberOfPadsPhi();

  if ( relid[1] == 0 ) { // its a PbW04 crystal 
    x = -( GetNPhi()/2. - row   + 0.5 ) *  GetCrystalSize(0) ; // position ox Xtal with respect
    z = ( GetNZ() /2. - column + 0.5 ) *  GetCrystalSize(2) ; // of center of PHOS module  
   }  
   else  {    
    if ( relid[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() )
       ppsdmodule =  relid[1]-GetNumberOfModulesPhi() *  GetNumberOfModulesZ(); 
    else ppsdmodule =  relid[1] ;
    Int_t modrow = 1+(Int_t)TMath::Ceil( (Float_t)ppsdmodule / GetNumberOfModulesPhi()-1. ) ; 
    Int_t modcol = ppsdmodule -  ( modrow - 1 ) * GetNumberOfModulesPhi() ;     
    Float_t x0 = (  GetNumberOfModulesPhi() / 2.  - modrow  + 0.5 ) * GetPPSDModuleSize(0) ;
    Float_t z0 = (  GetNumberOfModulesZ() / 2.  - modcol  + 0.5 ) * GetPPSDModuleSize(2)  ;     
    x = - ( GetNumberOfPadsPhi()/2. - row - 0.5 ) * padsizeX + x0 ; // position of pad  with respect
    z = ( GetNumberOfPadsZ()/2.   - column - 0.5 ) * padsizeZ - z0 ; // of center of PHOS module  
         }
}

//____________________________________________________________________________
void AliPHOSGeometry::SetPHOSAngles() 
{ 
  Double_t const kRADDEG = 180.0 / kPI ;
  Float_t pphi =  TMath::ATan( fOuterBoxSize[0]  / ( 2.0 * fIPtoOuterCoverDistance ) ) ;
  pphi *= kRADDEG ;
  
  for( Int_t i = 1; i <= fNModules ; i++ ) {
    Float_t angle = pphi * 2 * ( i - fNModules / 2.0 - 0.5 ) ;
    fPHOSAngle[i-1] = -  angle ;
  } 
}

//____________________________________________________________________________
void AliPHOSGeometry::SetLeadConverterThickness(Float_t e) 
{
  // should ultimately disappear 
  cout << " AliPHOSGeometry WARNING : You have changed LeadConverterThickness from " 
       << fLeadConverterThickness << " to " << e << endl ;

  fLeadConverterThickness = e ; 
}
