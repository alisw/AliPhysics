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
	            fNModules(0),
	            fAngle(0.f),
	            fPHOSAngle(0),
	            fIPtoUpperCPVsurface(0),
	            fRotMatrixArray(0),
	            fGeometryEMCA(0),
	            fGeometryCPV(0),
	            fGeometrySUPP(0)
{
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
    fgGeom          = 0 ;
}  

//____________________________________________________________________________
AliPHOSGeometry::AliPHOSGeometry(const AliPHOSGeometry & rhs)
		    : AliGeometry(rhs),
		      fNModules(rhs.fNModules),
		      fAngle(rhs.fAngle),
		      fPHOSAngle(0),
		      fIPtoUpperCPVsurface(rhs.fIPtoUpperCPVsurface),
		      fRotMatrixArray(0),
		      fGeometryEMCA(0),
		      fGeometryCPV(0),
		      fGeometrySUPP(0)
{
  Fatal("cpy ctor", "not implemented") ; 
}

//____________________________________________________________________________
AliPHOSGeometry::AliPHOSGeometry(const Text_t* name, const Text_t* title) 
	          : AliGeometry(name, title),
	            fNModules(0),
	            fAngle(0.f),
	            fPHOSAngle(0),
	            fIPtoUpperCPVsurface(0),
	            fRotMatrixArray(0),
	            fGeometryEMCA(0),
	            fGeometryCPV(0),
	            fGeometrySUPP(0)
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
  
  TString test(GetName()) ; 
  if (test != "IHEP" && test != "noCPV") {
    AliFatal(Form("%s is not a known geometry (choose among IHEP)", 
		  test.Data() )) ; 
  }

  fgInit     = kTRUE ; 

  fNModules     = 5;
  fAngle        = 20;

  fGeometryEMCA = new AliPHOSEMCAGeometry();
  
  fGeometryCPV  = new AliPHOSCPVGeometry ();
  
  fGeometrySUPP = new AliPHOSSupportGeometry();
  
  fPHOSAngle = new Float_t[fNModules] ;
  
  Float_t * emcParams = fGeometryEMCA->GetEMCParams() ;
  
  fPHOSParams[0] =  TMath::Max((Double_t)fGeometryCPV->GetCPVBoxSize(0)/2., 
 			       (Double_t)(emcParams[0] - (emcParams[1]-emcParams[0])*
					  fGeometryCPV->GetCPVBoxSize(1)/2/emcParams[3]));
  fPHOSParams[1] = emcParams[1] ;
  fPHOSParams[2] = TMath::Max((Double_t)emcParams[2], (Double_t)fGeometryCPV->GetCPVBoxSize(2)/2.);
  fPHOSParams[3] = emcParams[3] + fGeometryCPV->GetCPVBoxSize(1)/2. ;
  
  fIPtoUpperCPVsurface = fGeometryEMCA->GetIPtoOuterCoverDistance() - fGeometryCPV->GetCPVBoxSize(1) ;

  //calculate offset to crystal surface
  Float_t * inthermo = fGeometryEMCA->GetInnerThermoHalfSize() ;
  Float_t * strip = fGeometryEMCA->GetStripHalfSize() ;
  Float_t* splate = fGeometryEMCA->GetSupportPlateHalfSize();
  Float_t * crystal = fGeometryEMCA->GetCrystalHalfSize() ;
  Float_t * pin = fGeometryEMCA->GetAPDHalfSize() ;
  Float_t * preamp = fGeometryEMCA->GetPreampHalfSize() ;
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
Bool_t AliPHOSGeometry::AbsToRelNumbering(Int_t absId, Int_t * relid) const
{
  // Converts the absolute numbering into the following array
  //  relid[0] = PHOS Module number 1:fNModules 
  //  relid[1] = 0 if PbW04
  //           = -1 if CPV
  //  relid[2] = Row number inside a PHOS module
  //  relid[3] = Column number inside a PHOS module

  Bool_t rv  = kTRUE ; 
  Float_t id = absId ;

  Int_t phosmodulenumber = (Int_t)TMath:: Ceil( id / GetNCristalsInModule() ) ; 
  
  if ( phosmodulenumber >  GetNModules() ) { // it is a CPV pad
    
    id -=  GetNPhi() * GetNZ() *  GetNModules() ; 
    Float_t nCPV  = GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ() ;
    relid[0] = (Int_t) TMath::Ceil( id / nCPV ) ;
    relid[1] = -1 ;
    id -= ( relid[0] - 1 ) * nCPV ; 
    relid[2] = (Int_t) TMath::Ceil( id / GetNumberOfCPVPadsZ() ) ;
    relid[3] = (Int_t) ( id - ( relid[2] - 1 ) * GetNumberOfCPVPadsZ() ) ; 
  } 
  else { // it is a PW04 crystal

    relid[0] = phosmodulenumber ;
    relid[1] = 0 ;
    id -= ( phosmodulenumber - 1 ) *  GetNPhi() * GetNZ() ; 
    relid[2] = (Int_t)TMath::Ceil( id / GetNZ() )  ;
    relid[3] = (Int_t)( id - ( relid[2] - 1 ) * GetNZ() ) ; 
  } 
  return rv ; 
}
//____________________________________________________________________________
void AliPHOSGeometry::GetGlobal(const AliRecPoint* recPoint, TVector3 & gpos) const
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
  char path[100] ; 
  Double_t dy ;
  if(tmpPHOS->IsEmc()){
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",tmpPHOS->GetPHOSMod()) ;
//    sprintf(path,"/ALIC_1/PHOS_%d",tmpPHOS->GetPHOSMod()) ;
    dy=fCrystalShift ;
  }
  else{
    sprintf(path,"/ALIC_1/PHOS_%d/PCPV_1",tmpPHOS->GetPHOSMod()) ;
    dy= GetCPVBoxSize(1)/2. ; //center of CPV module 
  }
  Double_t pos[3]={gpos.X(),gpos.Y()-dy,gpos.Z()} ;
  if(tmpPHOS->IsEmc())
    pos[2]=-pos[2] ; //Opposite z directions in EMC matrix and local frame!!!
  Double_t posC[3];
  //now apply possible shifts and rotations
  if (!gGeoManager->cd(path)){
    AliFatal("Geo manager can not find path \n");
  }
  TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
  if (m){
     m->LocalToMaster(pos,posC);
  }
  else{
    AliFatal("Geo matrixes are not loaded \n") ;
  }
  gpos.SetXYZ(posC[0],posC[1],posC[2]) ;

}
//____________________________________________________________________________
void AliPHOSGeometry::ImpactOnEmc(Double_t * vtx, Double_t theta, Double_t phi, 
                                  Int_t & moduleNumber, Double_t & z, Double_t & x) const
{
  // calculates the impact coordinates on PHOS of a neutral particle  
  // emitted in the vertex vtx[3] with direction theta and phi in the ALICE global coordinate system
  TVector3 p(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) ;
  TVector3 v(vtx[0],vtx[1],vtx[2]) ;

  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
 
  for(Int_t imod=1; imod<=GetNModules() ; imod++){
    //create vector from (0,0,0) to center of crystal surface of imod module
    Double_t tmp[3]={0.,-fCrystalShift,0.} ;
 
    char path[100] ;
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",imod) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    Double_t posG[3]={0.,0.,0.} ;
    if (m) m->LocalToMaster(tmp,posG);
    TVector3 n(posG[0],posG[1],posG[2]) ; 
    Double_t direction=n.Dot(p) ;
    if(direction<=0.)
      continue ; //momentum directed FROM module
    Double_t fr = (n.Mag2()-n.Dot(v))/direction ;  
    //Calculate direction in module plain
    n-=v+fr*p ;
    n*=-1. ;
    Float_t * sz = fGeometryEMCA->GetInnerThermoHalfSize() ; //Wery close to the zise of the Xtl set
    if(TMath::Abs(TMath::Abs(n.Z())<sz[2]) && n.Pt()<sz[0]){
      moduleNumber = imod ;
      z=n.Z() ;
      x=TMath::Sign(n.Pt(),n.X()) ;
      //no need to return to local system since we calcilated distance from module center
      //and tilts can not be significant.
      return ;
    }
  }
  //Not in acceptance
  x=0; z=0 ;
  moduleNumber=0 ;

}

//____________________________________________________________________________
Bool_t  AliPHOSGeometry::Impact(const TParticle * particle) const 
{
  // Tells if a particle enters PHOS
  Bool_t in=kFALSE;
  Int_t moduleNumber=0;
  Double_t vtx[3]={particle->Vx(),particle->Vy(),particle->Vz()} ;
  Double_t z,x;
  ImpactOnEmc(vtx,particle->Theta(),particle->Phi(),moduleNumber,z,x);
  if(moduleNumber!=0) 
    in=kTRUE;
  return in;
}

//____________________________________________________________________________
Bool_t AliPHOSGeometry::RelToAbsNumbering(const Int_t * relid, Int_t &  absId) const
{
  // Converts the relative numbering into the absolute numbering
  // EMCA crystals:
  //  absId = from 1 to fNModules * fNPhi * fNZ
  // CPV pad:
  //  absId = from N(total PHOS crystals) + 1
  //          to NCPVModules * fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ

  Bool_t rv = kTRUE ; 
  
  if ( relid[1] ==  0 ) {                            // it is a Phos crystal
    absId =
      ( relid[0] - 1 ) * GetNPhi() * GetNZ()         // the offset of PHOS modules
      + ( relid[2] - 1 ) * GetNZ()                   // the offset along phi
      +   relid[3] ;                                 // the offset along z
  }
  else { // it is a CPV pad
    absId =    GetNPhi() * GetNZ() *  GetNModules()         // the offset to separate EMCA crystals from CPV pads
      + ( relid[0] - 1 ) * GetNumberOfCPVPadsPhi() * GetNumberOfCPVPadsZ()   // the pads offset of PHOS modules 
      + ( relid[2] - 1 ) * GetNumberOfCPVPadsZ()                             // the pads offset of a CPV row
      +   relid[3] ;                                                         // the column number
  }
  
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSGeometry::RelPosInAlice(Int_t id, TVector3 & pos ) const
{
  // Converts the absolute numbering into the global ALICE coordinate system
  
  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
    
  Int_t relid[4] ;
    
  AbsToRelNumbering(id , relid) ;
    
  //construct module name
  char path[100] ;
  if(relid[1]==0){ //this is EMC
 
    Double_t ps[3]= {0.0,-fCryCellShift,0.}; //Position incide the crystal 
    Double_t psC[3]={0.0,0.0,0.}; //Global position
 
    //Shift and possibly apply misalignment corrections
    Int_t nCellsXInStrip=fGeometryEMCA->GetNCellsXInStrip() ;
    Int_t nCellsZInStrip=fGeometryEMCA->GetNCellsZInStrip() ;
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/nCellsXInStrip))*fGeometryEMCA->GetNStripZ()-
                (Int_t) TMath::Ceil((Double_t)relid[3]/nCellsZInStrip) ;
    Int_t cellraw= relid[3]%nCellsZInStrip ;
    if(cellraw==0)cellraw=nCellsZInStrip ;
    Int_t cell= ((relid[2]-1)%nCellsXInStrip)*nCellsZInStrip + cellraw ;
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d/PCEL_%d",
            relid[0],strip,cell) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(ps,psC);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
    pos.SetXYZ(psC[0],psC[1],psC[2]) ; 
  }
  else{
    //first calculate position with respect to CPV plain
    Int_t row        = relid[2] ; //offset along x axis
    Int_t column     = relid[3] ; //offset along z axis
    Double_t ps[3]= {0.0,GetCPVBoxSize(1)/2.,0.}; //Position on top of CPV
    Double_t psC[3]={0.0,0.0,0.}; //Global position
    pos[0] = - ( GetNumberOfCPVPadsPhi()/2. - row    - 0.5 ) * GetPadSizePhi()  ; // position of pad  with respect
    pos[2] = - ( GetNumberOfCPVPadsZ()  /2. - column - 0.5 ) * GetPadSizeZ()  ; // of center of PHOS module
 
    //now apply possible shifts and rotations
    sprintf(path,"/ALIC_1/PHOS_%d/PCPV_1",relid[0]) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(ps,psC);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
    pos.SetXYZ(psC[0],psC[1],-psC[2]) ; 
  }
} 

//____________________________________________________________________________
void AliPHOSGeometry::RelPosToAbsId(Int_t module, Double_t x, Double_t z, Int_t & absId) const
{
  // converts local PHOS-module (x, z) coordinates to absId 

  //find Global position
  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
  Double_t posL[3]={x,-fCrystalShift,-z} ; //Only for EMC!!!
  Double_t posG[3] ;
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",module) ;
  if (!gGeoManager->cd(path)){
    AliFatal("Geo manager can not find path \n");
  }
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS){
     mPHOS->LocalToMaster(posL,posG);
  }
  else{
    AliFatal("Geo matrixes are not loaded \n") ;
  }

  Int_t relid[4] ;
  gGeoManager->FindNode(posG[0],posG[1],posG[2]) ;
  //Check that path contains PSTR and extract strip number
  TString cpath(gGeoManager->GetPath()) ;
  Int_t indx = cpath.Index("PCEL") ;
  if(indx==-1){ //for the few events when particle hits between srips use ideal geometry
    relid[0] = module ;
    relid[1] = 0 ;
    relid[2] = static_cast<Int_t>(TMath::Ceil( x/ GetCellStep() + GetNPhi() / 2.) );
    relid[3] = static_cast<Int_t>(TMath::Ceil(-z/ GetCellStep() + GetNZ()   / 2.) ) ;
    if(relid[2]<1)relid[2]=1 ;
    if(relid[3]<1)relid[3]=1 ;
    if(relid[2]>GetNPhi())relid[2]=GetNPhi() ;
    if(relid[3]>GetNZ())relid[3]=GetNZ() ;
    RelToAbsNumbering(relid,absId) ;
  }
  else{
    Int_t indx2 = cpath.Index("/",indx) ;
    if(indx2==-1)
      indx2=cpath.Length() ;
    TString cell=cpath(indx+5,indx2-indx-5) ;
    Int_t icell=cell.Atoi() ;
    indx = cpath.Index("PSTR") ;
    indx2 = cpath.Index("/",indx) ;
    TString strip=cpath(indx+5,indx2-indx-5) ;
    Int_t iStrip = strip.Atoi() ; 

    Int_t row = fGeometryEMCA->GetNStripZ() - (iStrip - 1) % (fGeometryEMCA->GetNStripZ()) ;
    Int_t col = (Int_t) TMath::Ceil((Double_t) iStrip/(fGeometryEMCA->GetNStripZ())) -1 ;
 
    // Absid for 8x2-strips. Looks nice :)
    absId = (module-1)*GetNCristalsInModule() +
                  row * 2 + (col*fGeometryEMCA->GetNCellsXInStrip() + (icell - 1) / 2)*GetNZ() - (icell & 1 ? 1 : 0);
 
  }
 
}

//____________________________________________________________________________
void AliPHOSGeometry::RelPosInModule(const Int_t * relid, Float_t & x, Float_t & z) const 
{
  // Converts the relative numbering into the local PHOS-module (x, z) coordinates
  // Note: sign of z differs from that in the previous version (Yu.Kharlov, 12 Oct 2000)
  

  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
  //construct module name
  char path[100] ;
  if(relid[1]==0){ //this is PHOS

//   Calculations using ideal geometry (obsolete)
//    x = - ( GetNPhi()/2. - relid[2]    + 0.5 ) *  GetCellStep() ; // position of Xtal with respect
//    z = - ( GetNZ()  /2. - relid[3] + 0.5 ) *  GetCellStep() ; // of center of PHOS module  

    Double_t pos[3]= {0.0,-fCryCellShift,0.}; //Position incide the crystal 
    Double_t posC[3]={0.0,0.0,0.}; //Global position

    //Shift and possibly apply misalignment corrections
    Int_t nCellsXInStrip=fGeometryEMCA->GetNCellsXInStrip() ;
    Int_t nCellsZInStrip=fGeometryEMCA->GetNCellsZInStrip() ;
//    Int_t strip=1+(relid[3]-1)/fGeometryEMCA->GetNCellsZInStrip()+((relid[2]-1)/nCellsInStrip)*fGeometryEMCA->GetNStripZ() ;
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/nCellsXInStrip))*fGeometryEMCA->GetNStripZ()-
                (Int_t) TMath::Ceil((Double_t)relid[3]/nCellsZInStrip) ;
    Int_t cellraw= relid[3]%nCellsZInStrip ;
    if(cellraw==0)cellraw=nCellsZInStrip ;
    Int_t cell= ((relid[2]-1)%nCellsXInStrip)*nCellsZInStrip + cellraw ; 
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d/PCEL_%d",
            relid[0],strip,cell) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(pos,posC);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
    //    printf("Local: x=%f, y=%f, z=%f \n",pos[0],pos[1],pos[2]) ;
    //    printf("   gl: x=%f, y=%f, z=%f \n",posC[0],posC[1],posC[2]) ;
    //Return to PHOS local system  
    Double_t posL[3]={posC[0],posC[1],posC[2]};
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",relid[0]) ;
    //    sprintf(path,"/ALIC_1/PHOS_%d",relid[0]) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
    if (mPHOS) mPHOS->MasterToLocal(posC,posL);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
//printf("RelPosInMod: posL=[%f,%f,%f]\n",posL[0],posL[1],posL[2]) ;
//printf("old: x=%f, z=%f \n",x,z);
    x=posL[0] ;
    z=-posL[2];
    return ;
  }
  else{//CPV
    //first calculate position with respect to CPV plain 
    Int_t row        = relid[2] ; //offset along x axis
    Int_t column     = relid[3] ; //offset along z axis
    Double_t pos[3]= {0.0,0.0,0.}; //Position incide the CPV printed circuit
    Double_t posC[3]={0.0,0.0,0.}; //Global position
    //    x = - ( GetNumberOfCPVPadsPhi()/2. - row    - 0.5 ) * GetPadSizePhi()  ; // position of pad  with respect
    //    z = - ( GetNumberOfCPVPadsZ()  /2. - column - 0.5 ) * GetPadSizeZ()  ; // of center of PHOS module
    pos[0] = - ( GetNumberOfCPVPadsPhi()/2. - row    - 0.5 ) * GetPadSizePhi()  ; // position of pad  with respect
    pos[2] = - ( GetNumberOfCPVPadsZ()  /2. - column - 0.5 ) * GetPadSizeZ()  ; // of center of PHOS module

    //now apply possible shifts and rotations
    sprintf(path,"/ALIC_1/PHOS_%d/PCPV_1",relid[0]) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(pos,posC);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
    //Return to PHOS local system
    Double_t posL[3]={0.,0.,0.,} ;
    sprintf(path,"/ALIC_1/PHOS_%d",relid[0]) ;
    if (!gGeoManager->cd(path)){
      AliFatal("Geo manager can not find path \n");
    }
    TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
    if (mPHOS) mPHOS->MasterToLocal(posC,posL);
    else{
      AliFatal("Geo matrixes are not loaded \n") ;
    }
    x=posL[0] ;
    z=posL[1];
    return ;
 
  }
  
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

//____________________________________________________________________________

void AliPHOSGeometry::Global2Local(TVector3& localPosition,
				   const TVector3& globalPosition,
				   Int_t module) const
{
  // Transforms a global position of the rec.point to the local coordinate system
  //Return to PHOS local system
  Double_t posG[3]={globalPosition.X(),globalPosition.Y(),globalPosition.Z()} ;
  Double_t posL[3]={0.,0.,0.} ;
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",module) ;
//  sprintf(path,"/ALIC_1/PHOS_%d",module) ;
  if (!gGeoManager->cd(path)){
    AliFatal("Geo manager can not find path \n");
  }
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS) mPHOS->MasterToLocal(posG,posL);
  else{
    AliFatal("Geo matrixes are not loaded \n") ;
  }
  localPosition.SetXYZ(posL[0],posL[1]+fCrystalShift,-posL[2]) ;  
 
/*
  Float_t angle = GetPHOSAngle(module); // (40,20,0,-20,-40) degrees
  angle *= TMath::Pi()/180;
  angle += 3*TMath::Pi()/2.;
  localPosition = globalPosition;
  localPosition.RotateZ(-angle);
*/
}
//____________________________________________________________________________
void AliPHOSGeometry::Local2Global(Int_t mod, Float_t x, Float_t z,
				   TVector3& globalPosition) const 
{
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",mod) ;
//  sprintf(path,"/ALIC_1/PHOS_%d",mod) ;
  if (!gGeoManager->cd(path)){
    AliFatal("Geo manager can not find path \n");
  }
  Double_t posL[3]={x,-fCrystalShift,-z} ; //Only for EMC!!!
  Double_t posG[3] ;
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS){
     mPHOS->LocalToMaster(posL,posG);
  }    
  else{
    AliFatal("Geo matrixes are not loaded \n") ;
  }
  globalPosition.SetXYZ(posG[0],posG[1],posG[2]) ;
}
//____________________________________________________________________________
void AliPHOSGeometry::GetIncidentVector(const TVector3 &vtx, Int_t module, Float_t x,Float_t z, TVector3 &vInc) const {
  //Calculates vector pointing from vertex to current poisition in module local frame
  //Note that PHOS local system and ALICE global have opposite z directions

  Global2Local(vInc,vtx,module) ; 
  vInc.SetXYZ(vInc.X()+x,vInc.Y(),vInc.Z()+z) ;
}
