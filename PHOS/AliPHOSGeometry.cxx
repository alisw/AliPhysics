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
// Geometry class  for PHOS : singleton  
// PHOS consists of the electromagnetic calorimeter (EMCA)
// and a charged particle veto either in the Subatech's version (PPSD)
// or in the IHEP's one (CPV).
// The EMCA/PPSD/CPV modules are parametrized so that any configuration
// can be easily implemented 
// The title is used to identify the version of CPV used.
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TVector3.h"
#include "TRotation.h" 

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliConst.h"

ClassImp(AliPHOSGeometry) ;

AliPHOSGeometry * AliPHOSGeometry::fgGeom = 0 ;
Bool_t            AliPHOSGeometry::fgInit = kFALSE ;

//____________________________________________________________________________
AliPHOSGeometry::~AliPHOSGeometry(void)
{
  // dtor

  if (fRotMatrixArray) fRotMatrixArray->Delete() ; 
  if (fRotMatrixArray) delete fRotMatrixArray ; 
  if (fPHOSAngle     ) delete fPHOSAngle ; 
}

//____________________________________________________________________________

void AliPHOSGeometry::Init(void)
{
  // Initializes the PHOS parameters

  if ( ((strcmp( fName, "default" )) == 0) || 
       ((strcmp( fName, "GPS2" ))    == 0) ||
       ((strcmp( fName, "IHEP" ))    == 0) ||
       ((strcmp( fName, "MIXT" ))    == 0) ) {
    fgInit     = kTRUE ; 

    fNModules     = 5;
    fNPPSDModules = 0;
    fAngle        = 20;

      fGeometryEMCA = new AliPHOSEMCAGeometry();
    if      ( ((strcmp( fName, "GPS2" ))  == 0) ) {
      fGeometryPPSD = new AliPHOSPPSDGeometry();
      fGeometryCPV  = 0;
      fNPPSDModules = fNModules;
    }
    else if ( ((strcmp( fName, "IHEP" ))  == 0) ) {
      fGeometryCPV  = new AliPHOSCPVGeometry ();
      fGeometryPPSD = 0;
      fNPPSDModules = 0;
    }
    else if ( ((strcmp( fName, "MIXT" ))  == 0) ) {
      fGeometryCPV  = new AliPHOSCPVGeometry ();
      fGeometryPPSD = new AliPHOSPPSDGeometry();
      fNPPSDModules = 1;
    }
      fGeometrySUPP = new AliPHOSSupportGeometry();

    fPHOSAngle = new Float_t[fNModules] ;
    Int_t index ;
    for ( index = 0; index < fNModules; index++ )
    fPHOSAngle[index] = 0.0 ; // Module position angles are set in CreateGeometry()

    this->SetPHOSAngles() ; 
    fRotMatrixArray = new TObjArray(fNModules) ; 
  }
  else {
    fgInit = kFALSE ; 
    cout << "PHOS Geometry setup: option not defined " << fName << endl ; 
  }
}

//____________________________________________________________________________
Float_t AliPHOSGeometry::GetCPVBoxSize(Int_t index)  const { 
    if      (strcmp(fName,"GPS2") ==0 ) 
      return fGeometryPPSD->GetCPVBoxSize(index);
    else if (strcmp(fName,"IHEP")==0) 
      return fGeometryCPV ->GetCPVBoxSize(index);
    else if (strcmp(fName,"MIXT")==0) 
      return TMath::Max(fGeometryCPV ->GetCPVBoxSize(index), fGeometryPPSD->GetCPVBoxSize(index));
    else                              
      return 0;
}  

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance() 
{ 
  // Returns the pointer of the unique instance
  return (AliPHOSGeometry *) fgGeom ; 
}

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance(const Text_t* name, const Text_t* title) 
{
  // Returns the pointer of the unique instance
  AliPHOSGeometry * rv = 0  ; 
  if ( fgGeom == 0 ) {
    if ( strcmp(name,"") == 0 ) 
      rv = 0 ;
    else {    
      fgGeom = new AliPHOSGeometry(name, title) ;
      if ( fgInit )
	rv = (AliPHOSGeometry * ) fgGeom ;
      else {
	rv = 0 ; 
	delete fgGeom ; 
	fgGeom = 0 ; 
      }
    }
  }
  else {
    if ( strcmp(fgGeom->GetName(), name) != 0 ) {
      cout << "AliPHOSGeometry <E> : current geometry is " << fgGeom->GetName() << endl
	   << "                      you cannot call     " << name << endl ; 
    }
    else
      rv = (AliPHOSGeometry *) fgGeom ; 
  } 
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSGeometry::SetPHOSAngles() 
{ 
  // Calculates the position in ALICE of the PHOS modules
  
  Double_t const kRADDEG = 180.0 / kPI ;
  Float_t pphi =  2 * TMath::ATan( GetOuterBoxSize(0)  / ( 2.0 * GetIPtoOuterCoverDistance() ) ) ;
  pphi *= kRADDEG ;
  if (pphi > fAngle) cout << "AliPHOSGeometry: PHOS modules overlap!\n";
  pphi = fAngle;
  
  for( Int_t i = 1; i <= fNModules ; i++ ) {
    Float_t angle = pphi * ( i - fNModules / 2.0 - 0.5 ) ;
    fPHOSAngle[i-1] = -  angle ;
  } 
}

//____________________________________________________________________________
Bool_t AliPHOSGeometry::AbsToRelNumbering(const Int_t AbsId, Int_t * relid)
{
  // Converts the absolute numbering into the following array/
  //  relid[0] = PHOS Module number 1:fNModules 
  //  relid[1] = 0 if PbW04
  //           = PPSD Module number 1:fNumberOfModulesPhi*fNumberOfModulesZ*2 (2->up and bottom level)
  //  relid[2] = Row number inside a PHOS or PPSD module
  //  relid[3] = Column number inside a PHOS or PPSD module

  Bool_t rv  = kTRUE ; 
  Float_t id = AbsId ;

  Int_t phosmodulenumber = (Int_t)TMath:: Ceil( id / ( GetNPhi() * GetNZ() ) ) ; 
  
  if ( phosmodulenumber >  GetNModules() ) { // it is a PPSD or CPV pad

    if      ( strcmp(fName,"GPS2") == 0 ) {
      id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
      Float_t tempo = 2 *  GetNumberOfModulesPhi() * GetNumberOfModulesZ() *  GetNumberOfPadsPhi() * GetNumberOfPadsZ() ; 
      relid[0] = (Int_t)TMath::Ceil( id / tempo ) ; 
      id -= ( relid[0] - 1 ) * tempo ;
      relid[1] = (Int_t)TMath::Ceil( id / ( GetNumberOfPadsPhi() * GetNumberOfPadsZ() ) ) ; 
      id -= ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ() ;
      relid[2] = (Int_t)TMath::Ceil( id / GetNumberOfPadsPhi() ) ;
      relid[3] = (Int_t) ( id - ( relid[2] - 1 )  * GetNumberOfPadsPhi() ) ; 
    }
    else if ( strcmp(fName,"IHEP") == 0 ) {
      id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
      Float_t nCPV  = GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ() ;
      relid[0] = (Int_t) TMath::Ceil( id / nCPV ) ;
      relid[1] = 1 ;
      id -= ( relid[0] - 1 ) * nCPV ; 
      relid[2] = (Int_t) TMath::Ceil( id / GetNumberOfCPVPadsZ() ) ;
      relid[3] = (Int_t) ( id - ( relid[2] - 1 ) * GetNumberOfCPVPadsZ() ) ; 
    }
    else if ( strcmp(fName,"MIXT") == 0 ) {
      id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
      Float_t nPPSD = 2 * GetNumberOfModulesPhi() * GetNumberOfModulesZ() *  GetNumberOfPadsPhi() * GetNumberOfPadsZ() ; 
      Float_t nCPV  = GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ() ;
      if (id <= nCPV*GetNCPVModules()) { // this pad belons to CPV
	relid[0] = (Int_t) TMath::Ceil( id / nCPV ) ;
	relid[1] = 1 ;
	id -= ( relid[0] - 1 ) * nCPV ; 
	relid[2] = (Int_t) TMath::Ceil( id / GetNumberOfCPVPadsZ() ) ;
	relid[3] = (Int_t) ( id - ( relid[2] - 1 ) * GetNumberOfCPVPadsZ() ) ; 
      }
      else {                             // this pad belons to PPSD
	id -= nCPV*GetNCPVModules();
	relid[0] = (Int_t)TMath::Ceil( id / nPPSD ); 
	id -= ( relid[0] - 1 ) * nPPSD ;
	relid[0] += GetNCPVModules();
	relid[1] = (Int_t)TMath::Ceil( id / ( GetNumberOfPadsPhi() * GetNumberOfPadsZ() ) ) ; 
	id -= ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ() ;
	relid[2] = (Int_t)TMath::Ceil( id / GetNumberOfPadsPhi() ) ;
	relid[3] = (Int_t) ( id - ( relid[2] - 1 )  * GetNumberOfPadsPhi() ) ; 
      }
    }
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
  if ( opt == Radian() ) 
    conv = 1. ; 
  else if ( opt == Degre() )
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
  if ( opt == Radian() ) 
    conv = 1. ; 
  else if ( opt == Degre() )
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
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos, TMatrix & gmat) const
{
  // Calculates the ALICE global coordinates of a RecPoint and the error matrix
 
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
void AliPHOSGeometry::GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos) const 
{
  // Calculates the ALICE global coordinates of a RecPoint 

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
Bool_t AliPHOSGeometry::RelToAbsNumbering(const Int_t * relid, Int_t &  AbsId)
{
  // Converts the relative numbering into the absolute numbering
  // EMCA crystals:
  //  AbsId = from 1 to fNModules * fNPhi * fNZ
  // PPSD gas cell:
  //  AbsId = from N(total EMCA crystals) + 1
  //          to NCPVModules * fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ +
  //          fNModules * 2 * (fNumberOfModulesPhi * fNumberOfModulesZ) * fNumberOfPadsPhi * fNumberOfPadsZ
  // CPV pad:
  //  AbsId = from N(total PHOS crystals) + 1
  //          to NCPVModules * fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ

  Bool_t rv = kTRUE ; 
 
  if      ( relid[1] > 0 && strcmp(fName,"GPS2")==0) { // it is a PPSD pad
    AbsId =    GetNPhi() * GetNZ() * GetNModules()                         // the offset to separate EMCA crystals from PPSD pads
      + ( relid[0] - 1 ) * GetNumberOfModulesPhi() * GetNumberOfModulesZ() // the pads offset of PPSD modules 
                         * GetNumberOfPadsPhi() * GetNumberOfPadsZ() * 2
      + ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ()       // the pads offset of PPSD modules 
      + ( relid[2] - 1 ) * GetNumberOfPadsPhi()                            // the pads offset of a PPSD row
      +   relid[3] ;                                                       // the column number
  } 

  else if ( relid[1] > 0 && strcmp(fName,"MIXT")==0) { // it is a PPSD pad
    AbsId =    GetNPhi() * GetNZ() * GetNModules()                         // the offset to separate EMCA crystals from PPSD pads
      + GetNCPVModules() * GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ() // the pads offset of CPV modules if any
      + ( relid[0] - 1 - GetNCPVModules())
                         * GetNumberOfModulesPhi() * GetNumberOfModulesZ() // the pads offset of PPSD modules 
                         * GetNumberOfPadsPhi() * GetNumberOfPadsZ() * 2
      + ( relid[1] - 1 ) * GetNumberOfPadsPhi() * GetNumberOfPadsZ()       // the pads offset of PPSD modules 
      + ( relid[2] - 1 ) * GetNumberOfPadsPhi()                            // the pads offset of a PPSD row
      +   relid[3] ;                                                       // the column number
  } 

  else if ( relid[1] ==  0 ) { // it is a Phos crystal
    AbsId =
        ( relid[0] - 1 ) * GetNPhi() * GetNZ()                               // the offset of PHOS modules
      + ( relid[2] - 1 ) * GetNPhi()                                         // the offset of a xtal row
      +   relid[3] ;                                                         // the column number
  }

  else if ( relid[1] == -1 ) { // it is a CPV pad
    AbsId =    GetNPhi() * GetNZ() *  GetNModules()                          // the offset to separate EMCA crystals from CPV pads
      + ( relid[0] - 1 ) * GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ()         // the pads offset of PHOS modules 
      + ( relid[2] - 1 ) * GetNumberOfCPVPadsZ()                                // the pads offset of a CPV row
      +   relid[3] ;                                                         // the column number
  }
  
  return rv ; 
}

//____________________________________________________________________________

void AliPHOSGeometry::RelPosInAlice(const Int_t id, TVector3 & pos ) 
{
  // Converts the absolute numbering into the global ALICE coordinates
  // It works only for the GPS2 geometry
  
  if (id > 0 && strcmp(fName,"GPS2")==0) { 
    
    Int_t relid[4] ;
    
    AbsToRelNumbering(id , relid) ;
    
    Int_t phosmodule = relid[0] ; 
    
    Float_t y0 = 0 ; 
    
    if ( relid[1] == 0 ) { // it is a PbW04 crystal 
      y0 =  -(GetIPtoOuterCoverDistance() + GetUpperPlateThickness()
	      + GetSecondUpperPlateThickness() + GetUpperCoolingPlateThickness())  ;  
    }
    if ( relid[1] > 0 ) { // its a PPSD pad
      if ( relid[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() ) { // its an bottom module
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
  // Converts the relative numbering into the local PHOS-module (x, z) coordinates
  // Note: sign of z differs from that in the previous version (Yu.Kharlov, 12 Oct 2000)

  Bool_t padOfCPV  = (strcmp(fName,"IHEP")==0) ||
                    ((strcmp(fName,"MIXT")==0) && relid[0]<=GetNCPVModules()) ;
  Bool_t padOfPPSD = (strcmp(fName,"GPS2")==0) ||
                    ((strcmp(fName,"MIXT")==0) && relid[0]> GetNCPVModules()) ;
  
  Int_t ppsdmodule  ; 
  Float_t x0,z0;
  Int_t row        = relid[2] ; //offset along x axiz
  Int_t column     = relid[3] ; //offset along z axiz

  Float_t padsizeZ = 0;
  Float_t padsizeX = 0;
  Int_t   nOfPadsPhi = 0;
  Int_t   nOfPadsZ   = 0;
  if      ( padOfPPSD ) {
    padsizeZ   = GetPPSDModuleSize(2) / GetNumberOfPadsZ();
    padsizeX   = GetPPSDModuleSize(0) / GetNumberOfPadsPhi();
    nOfPadsPhi = GetNumberOfPadsPhi();
    nOfPadsZ   = GetNumberOfPadsZ();
  }
  else if ( padOfCPV  ) {
    padsizeZ   = GetPadSizeZ();
    padsizeX   = GetPadSizePhi();
    nOfPadsPhi = GetNumberOfCPVPadsPhi();
    nOfPadsZ   = GetNumberOfCPVPadsZ();
  }
  
  if ( relid[1] == 0 ) { // its a PbW04 crystal 
    x = - ( GetNPhi()/2. - row    + 0.5 ) *  GetCrystalSize(0) ; // position ox Xtal with respect
    z =   ( GetNZ()  /2. - column + 0.5 ) *  GetCrystalSize(2) ; // of center of PHOS module  
  }  
  else  {    
    if ( padOfPPSD ) {
      if ( relid[1] >  GetNumberOfModulesPhi() *  GetNumberOfModulesZ() )
	ppsdmodule =  relid[1]-GetNumberOfModulesPhi() *  GetNumberOfModulesZ(); 
      else
	ppsdmodule =  relid[1] ;
      Int_t modrow = 1+(Int_t)TMath::Ceil( (Float_t)ppsdmodule / GetNumberOfModulesPhi()-1. ) ; 
      Int_t modcol = ppsdmodule -  ( modrow - 1 ) * GetNumberOfModulesPhi() ;     
      x0 = (  GetNumberOfModulesPhi() / 2.  - modrow  + 0.5 ) * GetPPSDModuleSize(0) ;
      z0 = (  GetNumberOfModulesZ()   / 2.  - modcol  + 0.5 ) * GetPPSDModuleSize(2)  ;     
    } else {
      x0 = 0;
      z0 = 0;
    }
    x = - ( nOfPadsPhi/2. - row    - 0.5 ) * padsizeX + x0 ; // position of pad  with respect
    z =   ( nOfPadsZ  /2. - column - 0.5 ) * padsizeZ - z0 ; // of center of PHOS module  
  }
}
