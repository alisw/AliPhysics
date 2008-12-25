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

/* $Id: AliPHOSGeometry.cxx 25590 2008-05-06 07:09:11Z prsnko $ */

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
#include "TParticle.h"
#include <TGeoManager.h>
#include <TGeoMatrix.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSEMCAGeometry.h"
#include "AliPHOSCPVGeometry.h"
#include "AliPHOSSupportGeometry.h"
#include "AliPHOSGeoUtils.h"

ClassImp(AliPHOSGeoUtils)

//____________________________________________________________________________
AliPHOSGeoUtils::AliPHOSGeoUtils():
  fGeometryEMCA(0x0),fGeometryCPV(0x0),fGeometrySUPP(0x0),
  fNModules(0),fNCristalsInModule(0),fNPhi(0),fNZ(0),
  fNumberOfCPVPadsPhi(0),fNumberOfCPVPadsZ(0),
  fNCellsXInStrip(0),fNCellsZInStrip(0),fNStripZ(0),
  fCrystalShift(0.),fCryCellShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.)
 
{
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
}  

//____________________________________________________________________________
AliPHOSGeoUtils::AliPHOSGeoUtils(const AliPHOSGeoUtils & rhs)
  : TNamed(rhs),
  fGeometryEMCA(0x0),fGeometryCPV(0x0),fGeometrySUPP(0x0),
  fNModules(0),fNCristalsInModule(0),fNPhi(0),fNZ(0),
  fNumberOfCPVPadsPhi(0),fNumberOfCPVPadsZ(0),
  fNCellsXInStrip(0),fNCellsZInStrip(0),fNStripZ(0),
  fCrystalShift(0.),fCryCellShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.)
{
  Fatal("cpy ctor", "not implemented") ; 
}

//____________________________________________________________________________
AliPHOSGeoUtils::AliPHOSGeoUtils(const Text_t* name, const Text_t* title) 
    : TNamed(name, title),
  fGeometryEMCA(0x0),fGeometryCPV(0x0),fGeometrySUPP(0x0),
  fNModules(0),fNCristalsInModule(0),fNPhi(0),fNZ(0),
  fNumberOfCPVPadsPhi(0),fNumberOfCPVPadsZ(0),
  fNCellsXInStrip(0),fNCellsZInStrip(0),fNStripZ(0),
  fCrystalShift(0.),fCryCellShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.)
{ 
  // ctor only for normal usage 

  fGeometryEMCA = new AliPHOSEMCAGeometry() ;
  fGeometryCPV  = new AliPHOSCPVGeometry() ;
  fGeometrySUPP = new AliPHOSSupportGeometry() ;

  fNModules     = 5;
  fNPhi  = fGeometryEMCA->GetNPhi() ;
  fNZ    = fGeometryEMCA->GetNZ() ;
  fNCristalsInModule = fNPhi*fNZ ;
  fNCellsXInStrip= fGeometryEMCA->GetNCellsXInStrip() ;
  fNCellsZInStrip= fGeometryEMCA->GetNCellsZInStrip() ;
  fNStripZ = fGeometryEMCA->GetNStripZ() ;
  fXtlArrSize[0]=fGeometryEMCA->GetInnerThermoHalfSize()[0] ; //Wery close to the zise of the Xtl set
  fXtlArrSize[1]=fGeometryEMCA->GetInnerThermoHalfSize()[1] ; //Wery close to the zise of the Xtl set
  fXtlArrSize[2]=fGeometryEMCA->GetInnerThermoHalfSize()[2] ; //Wery close to the zise of the Xtl set

  //calculate offset to crystal surface
  Float_t * inthermo = fGeometryEMCA->GetInnerThermoHalfSize() ;
  Float_t * strip = fGeometryEMCA->GetStripHalfSize() ;
  Float_t* splate = fGeometryEMCA->GetSupportPlateHalfSize();
  Float_t * crystal = fGeometryEMCA->GetCrystalHalfSize() ;
  Float_t * pin = fGeometryEMCA->GetAPDHalfSize() ;
  Float_t * preamp = fGeometryEMCA->GetPreampHalfSize() ;
  fCrystalShift=-inthermo[1]+strip[1]+splate[1]+crystal[1]-fGeometryEMCA->GetAirGapLed()/2.+pin[1]+preamp[1] ;
  fCryCellShift=crystal[1]-(fGeometryEMCA->GetAirGapLed()-2*pin[1]-2*preamp[1])/2;
  fCellStep = 2.*fGeometryEMCA->GetAirCellHalfSize()[0] ;


  fNumberOfCPVPadsPhi = fGeometryCPV->GetNumberOfCPVPadsPhi() ;
  fNumberOfCPVPadsZ   = fGeometryCPV->GetNumberOfCPVPadsZ() ;
  fPadSizePhi = fGeometryCPV->GetCPVPadSizePhi() ;
  fPadSizeZ   = fGeometryCPV->GetCPVPadSizeZ() ; 
  fCPVBoxSizeY= fGeometryCPV->GetCPVBoxSize(1) ;
}

//____________________________________________________________________________
AliPHOSGeoUtils & AliPHOSGeoUtils::operator = (const AliPHOSGeoUtils  & /*rvalue*/) { 

  Fatal("assignment operator", "not implemented") ; 
    return *this ;
}

//____________________________________________________________________________
AliPHOSGeoUtils::~AliPHOSGeoUtils(void)
{
  // dtor
  if(fGeometryEMCA){
    delete fGeometryEMCA; fGeometryEMCA = 0 ;
  }
  if(fGeometryCPV){
    delete fGeometryCPV; fGeometryCPV=0 ;
  }
  if(fGeometrySUPP){
    delete fGeometrySUPP ; fGeometrySUPP=0 ;
  }

}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::AbsToRelNumbering(Int_t absId, Int_t * relid) const
{
  // Converts the absolute numbering into the following array
  //  relid[0] = PHOS Module number 1:fNModules 
  //  relid[1] = 0 if PbW04
  //           = -1 if CPV
  //  relid[2] = Row number inside a PHOS module
  //  relid[3] = Column number inside a PHOS module

  Float_t id = absId ;

  Int_t phosmodulenumber = (Int_t)TMath:: Ceil( id / fNCristalsInModule ) ; 
  
  if ( phosmodulenumber >  fNModules ) { // it is a CPV pad
    
    id -=  fNPhi * fNZ *  fNModules ; 
    Float_t nCPV  = fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ ;
    relid[0] = (Int_t) TMath::Ceil( id / nCPV ) ;
    relid[1] = -1 ;
    id -= ( relid[0] - 1 ) * nCPV ; 
    relid[2] = (Int_t) TMath::Ceil( id / fNumberOfCPVPadsZ ) ;
    relid[3] = (Int_t) ( id - ( relid[2] - 1 ) * fNumberOfCPVPadsZ ) ; 
  } 
  else { // it is a PW04 crystal

    relid[0] = phosmodulenumber ;
    relid[1] = 0 ;
    id -= ( phosmodulenumber - 1 ) *  fNPhi * fNZ ; 
    relid[2] = (Int_t)TMath::Ceil( id / fNZ )  ;
    relid[3] = (Int_t)( id - ( relid[2] - 1 ) * fNZ ) ; 
  } 
  return kTRUE ; 
}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::RelToAbsNumbering(const Int_t * relid, Int_t &  absId) const
{
  // Converts the relative numbering into the absolute numbering
  // EMCA crystals:
  //  absId = from 1 to fNModules * fNPhi * fNZ
  // CPV pad:
  //  absId = from N(total PHOS crystals) + 1
  //          to NCPVModules * fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ

  if ( relid[1] ==  0 ) {                            // it is a Phos crystal
    absId =
      ( relid[0] - 1 ) * fNPhi * fNZ         // the offset of PHOS modules
      + ( relid[2] - 1 ) * fNZ                   // the offset along phi
      +   relid[3] ;                                 // the offset along z
  }
  else { // it is a CPV pad
    absId =    fNPhi * fNZ *  fNModules         // the offset to separate EMCA crystals from CPV pads
      + ( relid[0] - 1 ) * fNumberOfCPVPadsPhi * fNumberOfCPVPadsZ   // the pads offset of PHOS modules 
      + ( relid[2] - 1 ) * fNumberOfCPVPadsZ                            // the pads offset of a CPV row
      +   relid[3] ;                                                         // the column number
  }
  
  return kTRUE ; 
}

//____________________________________________________________________________
void AliPHOSGeoUtils::RelPosInModule(const Int_t * relid, Float_t & x, Float_t & z) const 
{
  // Converts the relative numbering into the local PHOS-module (x, z) coordinates

  if (!gGeoManager){
    printf("Geo manager not initialized\n");
    abort() ;
  }
  //construct module name
  char path[100] ;
  if(relid[1]==0){ //this is PHOS

    Double_t pos[3]= {0.0,-fCryCellShift,0.}; //Position incide the crystal 
    Double_t posC[3]={0.0,0.0,0.}; //Global position

    //Shift and possibly apply misalignment corrections
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/fNCellsXInStrip))*fNStripZ-
                (Int_t) TMath::Ceil((Double_t)relid[3]/fNCellsZInStrip) ;
    Int_t cellraw= relid[3]%fNCellsZInStrip ;
    if(cellraw==0)cellraw=fNCellsZInStrip ;
    Int_t cell= ((relid[2]-1)%fNCellsXInStrip)*fNCellsZInStrip + cellraw ; 
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d/PCEL_%d",
            relid[0],strip,cell) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort() ;
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(pos,posC);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
    //Return to PHOS local system  
    Double_t posL[3]={posC[0],posC[1],posC[2]};
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",relid[0]) ;
    //    sprintf(path,"/ALIC_1/PHOS_%d",relid[0]) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort();
    }
    TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
    if (mPHOS) mPHOS->MasterToLocal(posC,posL);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
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
    pos[0] = - ( fNumberOfCPVPadsPhi/2. - row    - 0.5 ) * fPadSizePhi  ; // position of pad  with respect
    pos[2] = - ( fNumberOfCPVPadsZ  /2. - column - 0.5 ) * fPadSizeZ  ; // of center of PHOS module

    //now apply possible shifts and rotations
    sprintf(path,"/ALIC_1/PHOS_%d/PCPV_1",relid[0]) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort() ;
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(pos,posC);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
    //Return to PHOS local system
    Double_t posL[3]={0.,0.,0.,} ;
    sprintf(path,"/ALIC_1/PHOS_%d",relid[0]) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort() ;
    }
    TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
    if (mPHOS) mPHOS->MasterToLocal(posC,posL);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
    x=posL[0] ;
    z=posL[1];
    return ;
 
  }
  
}
//____________________________________________________________________________
void AliPHOSGeoUtils::RelPosToAbsId(Int_t module, Double_t x, Double_t z, Int_t & absId) const
{
  // converts local PHOS-module (x, z) coordinates to absId 

  //find Global position
  if (!gGeoManager){
    printf("Geo manager not initialized\n");
    abort() ;
  }
  Double_t posL[3]={x,-fCrystalShift,-z} ; //Only for EMC!!!
  Double_t posG[3] ;
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",module) ;
  if (!gGeoManager->cd(path)){
    printf("Geo manager can not find path \n");
    abort() ;
  }
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS){
     mPHOS->LocalToMaster(posL,posG);
  }
  else{
    printf("Geo matrixes are not loaded \n") ;
    abort() ;
  }

  Int_t relid[4] ;
  gGeoManager->FindNode(posG[0],posG[1],posG[2]) ;
  //Check that path contains PSTR and extract strip number
  TString cpath(gGeoManager->GetPath()) ;
  Int_t indx = cpath.Index("PCEL") ;
  if(indx==-1){ //for the few events when particle hits between srips use ideal geometry
    relid[0] = module ;
    relid[1] = 0 ;
    relid[2] = static_cast<Int_t>(TMath::Ceil( x/ fCellStep + fNPhi / 2.) );
    relid[3] = static_cast<Int_t>(TMath::Ceil(-z/ fCellStep + fNZ   / 2.) ) ;
    if(relid[2]<1)relid[2]=1 ;
    if(relid[3]<1)relid[3]=1 ;
    if(relid[2]>fNPhi)relid[2]=fNPhi ;
    if(relid[3]>fNZ)relid[3]=fNZ ;
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

    Int_t row = fNStripZ - (iStrip - 1) % (fNStripZ) ;
    Int_t col = (Int_t) TMath::Ceil((Double_t) iStrip/(fNStripZ)) -1 ;
 
    // Absid for 8x2-strips. Looks nice :)
    absId = (module-1)*fNCristalsInModule +
                  row * 2 + (col*fNCellsXInStrip + (icell - 1) / 2)*fNZ - (icell & 1 ? 1 : 0);
 
  }
 
}

//____________________________________________________________________________
void AliPHOSGeoUtils::RelPosToRelId(Int_t module, Double_t x, Double_t z, Int_t * relId) const
{
  //Evaluates RelId of the crystall with given coordinates

  Int_t absId ;
  RelPosToAbsId(module, x,z,absId) ;
  AbsToRelNumbering(absId,relId) ;
}

//____________________________________________________________________________
void AliPHOSGeoUtils::RelPosInAlice(Int_t id, TVector3 & pos ) const
{
  // Converts the absolute numbering into the global ALICE coordinate system
  
  if (!gGeoManager){
    printf("Geo manager not initialized\n");
    abort();
  }
    
  Int_t relid[4] ;
    
  AbsToRelNumbering(id , relid) ;
    
  //construct module name
  char path[100] ;
  if(relid[1]==0){ //this is EMC
 
    Double_t ps[3]= {0.0,-fCryCellShift,0.}; //Position incide the crystal 
    Double_t psC[3]={0.0,0.0,0.}; //Global position
 
    //Shift and possibly apply misalignment corrections
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/fNCellsXInStrip))*fNStripZ-
                (Int_t) TMath::Ceil((Double_t)relid[3]/fNCellsZInStrip) ;
    Int_t cellraw= relid[3]%fNCellsZInStrip ;
    if(cellraw==0)cellraw=fNCellsZInStrip ;
    Int_t cell= ((relid[2]-1)%fNCellsXInStrip)*fNCellsZInStrip + cellraw ;
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d/PCEL_%d",
            relid[0],strip,cell) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort() ;
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(ps,psC);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
    pos.SetXYZ(psC[0],psC[1],psC[2]) ; 
  }
  else{
    //first calculate position with respect to CPV plain
    Int_t row        = relid[2] ; //offset along x axis
    Int_t column     = relid[3] ; //offset along z axis
    Double_t ps[3]= {0.0,fCPVBoxSizeY/2.,0.}; //Position on top of CPV
    Double_t psC[3]={0.0,0.0,0.}; //Global position
    pos[0] = - ( fNumberOfCPVPadsPhi/2. - row    - 0.5 ) * fPadSizePhi  ; // position of pad  with respect
    pos[2] = - ( fNumberOfCPVPadsZ  /2. - column - 0.5 ) * fPadSizeZ  ; // of center of PHOS module
 
    //now apply possible shifts and rotations
    sprintf(path,"/ALIC_1/PHOS_%d/PCPV_1",relid[0]) ;
    if (!gGeoManager->cd(path)){
      printf("Geo manager can not find path \n");
      abort();
    }
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
    if (m) m->LocalToMaster(ps,psC);
    else{
      printf("Geo matrixes are not loaded \n") ;
      abort() ;
    }
    pos.SetXYZ(psC[0],psC[1],-psC[2]) ; 
  }
} 

//____________________________________________________________________________
void AliPHOSGeoUtils::Local2Global(Int_t mod, Float_t x, Float_t z,
				   TVector3& globalPosition) const 
{
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",mod) ;
  if (!gGeoManager->cd(path)){
    printf("Geo manager can not find path \n");
    abort() ;
  }
  Double_t posL[3]={x,-fCrystalShift,-z} ; //Only for EMC!!!
  Double_t posG[3] ;
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS){
     mPHOS->LocalToMaster(posL,posG);
  }    
  else{
    printf("Geo matrixes are not loaded \n") ;
    abort() ;
  }
  globalPosition.SetXYZ(posG[0],posG[1],posG[2]) ;
}
//____________________________________________________________________________
void AliPHOSGeoUtils::Global2Local(TVector3& localPosition,
				   const TVector3& globalPosition,
				   Int_t module) const
{
  // Transforms a global position to the local coordinate system
  // of the module 
  //Return to PHOS local system
  Double_t posG[3]={globalPosition.X(),globalPosition.Y(),globalPosition.Z()} ;
  Double_t posL[3]={0.,0.,0.} ;
  char path[100] ;
  sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",module) ;
  if (!gGeoManager->cd(path)){
    printf("Geo manager can not find path \n");
    abort() ;
  }
  TGeoHMatrix *mPHOS = gGeoManager->GetCurrentMatrix();
  if (mPHOS) mPHOS->MasterToLocal(posG,posL);
  else{
    printf("Geo matrixes are not loaded \n") ;
    abort() ;
  }
  localPosition.SetXYZ(posL[0],posL[1]+fCrystalShift,-posL[2]) ;  
 
}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::GlobalPos2RelId(TVector3 & global, Int_t * relId){
  //Converts position in global ALICE coordinates to relId 
  //returns false if x,z coordinates are beyond PHOS
  //distande to PHOS surface is NOT calculated 
  TVector3 loc ;
  for(Int_t mod=1; mod<fNModules; mod++){
    Global2Local(loc,global,mod) ;
    //If in Acceptance
    if((TMath::Abs(loc.Z())<fXtlArrSize[2]) && (TMath::Abs(loc.X())<fXtlArrSize[0])){
       RelPosToRelId(mod,loc.X(),loc.Z(),relId);
       return kTRUE ;
    }
  }
  return kFALSE ; 

}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::ImpactOnEmc(const TParticle * particle,
       Int_t & moduleNumber, Double_t & z, Double_t & x) const
{
  // Tells if a particle enters PHOS and evaluates hit position
  Double_t vtx[3]={particle->Vx(),particle->Vy(),particle->Vz()} ;
  return ImpactOnEmc(vtx,particle->Theta(),particle->Phi(),moduleNumber,z,x);
}
 
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::ImpactOnEmc(const Double_t * vtx, Double_t theta, Double_t phi, 
                                  Int_t & moduleNumber, Double_t & z, Double_t & x) const
{
  // calculates the impact coordinates on PHOS of a neutral particle  
  // emitted in the vertex vtx[3] with direction vec(p) in the ALICE global coordinate system
  TVector3 p(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) ;
  return ImpactOnEmc(vtx,p,moduleNumber,z,x) ;

}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::ImpactOnEmc(const Double_t * vtx, const TVector3 &p,
                                  Int_t & moduleNumber, Double_t & z, Double_t & x) const
{
  // calculates the impact coordinates on PHOS of a neutral particle  
  // emitted in the vertex vtx[3] with direction theta and phi in the ALICE global coordinate system
  TVector3 v(vtx[0],vtx[1],vtx[2]) ;

  if (!gGeoManager){
    printf("Geo manager not initialized\n");
    abort() ;
    return kFALSE ;
  }
 
  for(Int_t imod=1; imod<=fNModules ; imod++){
    //create vector from (0,0,0) to center of crystal surface of imod module
    Double_t tmp[3]={0.,-fCrystalShift,0.} ;
 
    char path[100] ;
    sprintf(path,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",imod) ;
    if (!gGeoManager->cd(path)){ //Module does not present
      continue ;
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
    if(TMath::Abs(TMath::Abs(n.Z())<fXtlArrSize[2]) && n.Pt()<fXtlArrSize[0]){
      moduleNumber = imod ;
      z=n.Z() ;
      x=TMath::Sign(n.Pt(),n.X()) ;
      //no need to return to local system since we calcilated distance from module center
      //and tilts can not be significant.
      return kTRUE ;
    }
  }
  //Not in acceptance
  x=0; z=0 ;
  moduleNumber=0 ;
  return kFALSE ;

}
//____________________________________________________________________________
void AliPHOSGeoUtils::GetIncidentVector(const TVector3 &vtx, Int_t module, Float_t x,Float_t z, TVector3 &vInc) const {
  //Calculates vector pointing from vertex to current poisition in module local frame
  //Note that PHOS local system and ALICE global have opposite z directions

  Global2Local(vInc,vtx,module) ; 
  vInc.SetXYZ(vInc.X()+x,vInc.Y(),vInc.Z()+z) ;
}

