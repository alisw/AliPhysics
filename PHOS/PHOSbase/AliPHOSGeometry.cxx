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
// -- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC "KI" & SUBATECH)

// --- ROOT system ---

#include "TVector3.h"
#include "TRotation.h" 
#include "TParticle.h"
#include <TGeoManager.h>
#include <TGeoMatrix.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEMCAGeometry.h" 
#include "AliPHOSRecPoint.h"

ClassImp(AliPHOSGeometry)

// these initialisations are needed for a singleton
AliPHOSGeometry  * AliPHOSGeometry::fgGeom = 0 ;
Bool_t             AliPHOSGeometry::fgInit = kFALSE ;

//____________________________________________________________________________
AliPHOSGeometry::AliPHOSGeometry() : 
                    AliPHOSGeoUtils(),
	            fAngle(0.f),
	            fPHOSAngle(0),
	            fIPtoUpperCPVsurface(0),
		    fCrystalShift(0),
		    fCryCellShift(0),
	            fRotMatrixArray(0)
{
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
    fgGeom          = 0 ;

    fPHOSParams[0] = 0.;
    fPHOSParams[1] = 0.;
    fPHOSParams[2] = 0.;
    fPHOSParams[3] = 0.;
}  

//____________________________________________________________________________
AliPHOSGeometry::AliPHOSGeometry(const AliPHOSGeometry & rhs)
		    : AliPHOSGeoUtils(rhs),
		      fAngle(rhs.fAngle),
		      fPHOSAngle(0),
		      fIPtoUpperCPVsurface(rhs.fIPtoUpperCPVsurface),
		      fCrystalShift(rhs.fCrystalShift),
		      fCryCellShift(rhs.fCryCellShift),
		      fRotMatrixArray(0)
{
  Fatal("cpy ctor", "not implemented") ; 
}

//____________________________________________________________________________
AliPHOSGeometry::AliPHOSGeometry(const Text_t* name, const Text_t* title) 
	          : AliPHOSGeoUtils(name, title),
	            fAngle(0.f),
	            fPHOSAngle(0),
	            fIPtoUpperCPVsurface(0),
		    fCrystalShift(0),
		    fCryCellShift(0),
	            fRotMatrixArray(0)
{ 
  // ctor only for internal usage (singleton)
  Init() ; 
  fgGeom = this;
}

//____________________________________________________________________________
AliPHOSGeometry::~AliPHOSGeometry(void)
{
  // dtor

  if (fRotMatrixArray) fRotMatrixArray->Delete() ; 
  if (fRotMatrixArray) delete fRotMatrixArray ; 
  if (fPHOSAngle     ) delete[] fPHOSAngle ; 
}

//____________________________________________________________________________
void AliPHOSGeometry::Init(void)
{
  // Initializes the PHOS parameters :
  //  IHEP is the Protvino CPV (cathode pad chambers)
  
  fgInit     = kTRUE ; 

  fAngle        = 20;

  
  fPHOSAngle = new Float_t[fNModules] ;
  
  const Float_t * emcParams = fGeometryEMCA->GetEMCParams() ;
  
  fPHOSParams[0] =  TMath::Max((Double_t)fGeometryCPV->GetCPVBoxSize(0)/2., 
 			       (Double_t)(emcParams[0] - (emcParams[1]-emcParams[0])*
					  fGeometryCPV->GetCPVBoxSize(1)/2/emcParams[3]));
  fPHOSParams[1] = emcParams[1] ;
  fPHOSParams[2] = TMath::Max((Double_t)emcParams[2], (Double_t)fGeometryCPV->GetCPVBoxSize(2)/2.);
  fPHOSParams[3] = emcParams[3] + fGeometryCPV->GetCPVBoxSize(1)/2. ;
  
  fIPtoUpperCPVsurface = fGeometryEMCA->GetIPtoOuterCoverDistance() - fGeometryCPV->GetCPVBoxSize(1) ;

  //calculate offset to crystal surface
  const Float_t * inthermo = fGeometryEMCA->GetInnerThermoHalfSize() ;
  const Float_t * strip = fGeometryEMCA->GetStripHalfSize() ;
  const Float_t * splate = fGeometryEMCA->GetSupportPlateHalfSize();
  const Float_t * crystal = fGeometryEMCA->GetCrystalHalfSize() ;
  const Float_t * pin = fGeometryEMCA->GetAPDHalfSize() ;
  const Float_t * preamp = fGeometryEMCA->GetPreampHalfSize() ;
  fCrystalShift=-inthermo[1]+strip[1]+splate[1]+crystal[1]-fGeometryEMCA->GetAirGapLed()/2.+pin[1]+preamp[1] ;
  fCryCellShift=crystal[1]-(fGeometryEMCA->GetAirGapLed()-2*pin[1]-2*preamp[1])/2;
 
  Int_t index ;
  for ( index = 0; index < fNModules; index++ )
    fPHOSAngle[index] = 0.0 ; // Module position angles are set in CreateGeometry()
  
  fRotMatrixArray = new TObjArray(fNModules) ; 

  // Geometry parameters are calculated

  SetPHOSAngles();
  Double_t const kRADDEG = 180.0 / TMath::Pi() ;
  Float_t r = GetIPtoOuterCoverDistance() + fPHOSParams[3] - GetCPVBoxSize(1) ;
  for (Int_t iModule=0; iModule<fNModules; iModule++) {
    fModuleCenter[iModule][0] = r * TMath::Sin(fPHOSAngle[iModule] / kRADDEG );
    fModuleCenter[iModule][1] =-r * TMath::Cos(fPHOSAngle[iModule] / kRADDEG );
    fModuleCenter[iModule][2] = 0.;
    
    fModuleAngle[iModule][0][0] =  90;
    fModuleAngle[iModule][0][1] =   fPHOSAngle[iModule];
    fModuleAngle[iModule][1][0] =   0;
    fModuleAngle[iModule][1][1] =   0;
    fModuleAngle[iModule][2][0] =  90;
    fModuleAngle[iModule][2][1] = 270 + fPHOSAngle[iModule];
  }

}

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance() 
{ 
  // Returns the pointer of the unique instance; singleton specific
  
  return static_cast<AliPHOSGeometry *>( fgGeom ) ; 
}

//____________________________________________________________________________
AliPHOSGeometry *  AliPHOSGeometry::GetInstance(const Text_t* name, const Text_t* title) 
{
  // Returns the pointer of the unique instance
  // Creates it with the specified options (name, title) if it does not exist yet

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
    if ( strcmp(fgGeom->GetName(), name) != 0 ) 
      ::Error("GetInstance", "Current geometry is %s. You cannot call %s", 
		      fgGeom->GetName(), name) ; 
    else
      rv = (AliPHOSGeometry *) fgGeom ; 
  } 
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSGeometry::SetPHOSAngles() 
{ 
  // Calculates the position of the PHOS modules in ALICE global coordinate system
  // in ideal geometry
  
  Double_t const kRADDEG = 180.0 / TMath::Pi() ;
  Float_t pphi =  2 * TMath::ATan( GetOuterBoxSize(0)  / ( 2.0 * GetIPtoUpperCPVsurface() ) ) ;
  pphi *= kRADDEG ;
  if (pphi > fAngle){ 
    AliError(Form("PHOS modules overlap!\n pphi = %f fAngle = %f", 
		  pphi, fAngle));

  }
  pphi = fAngle;
  
  for( Int_t i = 1; i <= fNModules ; i++ ) {
    Float_t angle = pphi * ( i - fNModules / 2.0 - 0.5 ) ;
    fPHOSAngle[i-1] = -  angle ;
  } 
}
//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* , TVector3 & ) const
{
  AliFatal(Form("Please use GetGlobalPHOS(recPoint,gpos) instead of GetGlobal!"));
}

//____________________________________________________________________________
void AliPHOSGeometry::GetGlobalPHOS(const AliPHOSRecPoint* recPoint, TVector3 & gpos) const
{
  // Calculates the coordinates of a RecPoint and the error matrix in the ALICE global coordinate system
 
  const AliPHOSRecPoint * tmpPHOS = recPoint ;  
  TVector3 localposition ;

  tmpPHOS->GetLocalPosition(gpos) ;

  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
  //construct module name
  TGeoHMatrix *m = 0x0;
  char path[100] ; 
  Double_t dy ;
  if(tmpPHOS->IsEmc()){
    snprintf(path,100,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",tmpPHOS->GetPHOSMod()) ;
    if (!gGeoManager->CheckPath(path)){
      snprintf(path,100,"/ALIC_1/PHOC_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",tmpPHOS->GetPHOSMod()) ;
      if (!gGeoManager->CheckPath(path)){
        snprintf(path,100,"/ALIC_1/PHOH_%d/PEMH_1/PCLH_1/PIOH_1/PCOH_1/PAGH_1/PTIH_1",tmpPHOS->GetPHOSMod()) ;
        if(!gGeoManager->CheckPath(path)){     
          AliFatal("Geo manager can not find path \n");
	}
      }
    }
    gGeoManager->cd(path) ;
    m = gGeoManager->GetCurrentMatrix();
    dy=fCrystalShift ;
  }
  else{
    snprintf(path,100,"/ALIC_1/PHOC_%d/PCPV_1",tmpPHOS->GetPHOSMod());
    if (!gGeoManager->CheckPath(path)){
      snprintf(path,100,"/ALIC_1/PHOH_%d/PCPV_1",tmpPHOS->GetPHOSMod());
      if (!gGeoManager->CheckPath(path))
        AliFatal(Form("Geo manager can not find path /ALIC_1/PHOC(H)_%d/PCPV_1 \n",tmpPHOS->GetPHOSMod()));
    }
    gGeoManager->cd(path) ;
    m = gGeoManager->GetCurrentMatrix();
    dy= GetCPVBoxSize(1)/2. ; //center of CPV module 
  }
  Double_t pos[3]={gpos.X(),gpos.Y()-dy,gpos.Z()} ;
  if(tmpPHOS->IsEmc())
    pos[2]=-pos[2] ; //Opposite z directions in EMC matrix and local frame!!!
  Double_t posC[3] = {};
  //now apply possible shifts and rotations
  if (m){
     m->LocalToMaster(pos,posC);
  }
  else{
    AliFatal("Geo matrixes are not loaded \n") ;
  }
  gpos.SetXYZ(posC[0],posC[1],posC[2]) ;

}
//____________________________________________________________________________

void AliPHOSGeometry::GetModuleCenter(TVector3& center, 
				      const char *det,
				      Int_t module) const
{
  // Returns a position of the center of the CPV or EMC module
  // in ideal (not misaligned) geometry
  Float_t rDet = 0.;
  if      (strcmp(det,"CPV") == 0) rDet  = GetIPtoCPVDistance   ();
  else if (strcmp(det,"EMC") == 0) rDet  = GetIPtoCrystalSurface();
  else 
    AliFatal(Form("Wrong detector name %s",det));

  Float_t angle = GetPHOSAngle(module); // (40,20,0,-20,-40) degrees
  angle *= TMath::Pi()/180;
  angle += 3*TMath::Pi()/2.;
  center.SetXYZ(rDet*TMath::Cos(angle), rDet*TMath::Sin(angle), 0.);
}

