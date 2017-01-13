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
// Geometry class  for PHOS 
// PHOS consists of the electromagnetic calorimeter (EMCA)
// and a charged particle veto (CPV)
// The EMCA/CPV modules are parametrized so that any configuration
// can be easily implemented 
// The title is used to identify the version of CPV used.
//                  
// -- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC "KI" & SUBATECH)

// --- ROOT system ---

#include "TClonesArray.h"
#include "TVector3.h"
#include "TParticle.h"
#include <TGeoManager.h>
#include <TGeoMatrix.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
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
  fCrystalShift(0.),fCryCellShift(0.),fCryStripShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.),fMisalArray(0x0)
 
{
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
  
  fXtlArrSize[0]=0.;   
  fXtlArrSize[1]=0.;                                                                           
  fXtlArrSize[2]=0.; 
  
  for(Int_t mod=0; mod<5; mod++){
    fEMCMatrix[mod]=0 ;
    for(Int_t istrip=0; istrip<224; istrip++)
      fStripMatrix[mod][istrip]=0 ;
    fCPVMatrix[mod]=0;
    fPHOSMatrix[mod]=0 ;
  }

}  

//____________________________________________________________________________
AliPHOSGeoUtils::AliPHOSGeoUtils(const AliPHOSGeoUtils & rhs)
  : TNamed(rhs),
  fGeometryEMCA(0x0),fGeometryCPV(0x0),fGeometrySUPP(0x0),
  fNModules(0),fNCristalsInModule(0),fNPhi(0),fNZ(0),
  fNumberOfCPVPadsPhi(0),fNumberOfCPVPadsZ(0),
  fNCellsXInStrip(0),fNCellsZInStrip(0),fNStripZ(0),
  fCrystalShift(0.),fCryCellShift(0.),fCryStripShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.),fMisalArray(0x0)
{
  Fatal("cpy ctor", "not implemented") ; 
  for(Int_t mod=0; mod<5; mod++){
    fEMCMatrix[mod]=0 ;
    for(Int_t istrip=0; istrip<224; istrip++)
      fStripMatrix[mod][istrip]=0 ;
    fCPVMatrix[mod]=0;
    fPHOSMatrix[mod]=0 ;
  }
}

//____________________________________________________________________________
AliPHOSGeoUtils::AliPHOSGeoUtils(const Text_t* name, const Text_t* title) 
    : TNamed(name, title),
  fGeometryEMCA(0x0),fGeometryCPV(0x0),fGeometrySUPP(0x0),
  fNModules(0),fNCristalsInModule(0),fNPhi(0),fNZ(0),
  fNumberOfCPVPadsPhi(0),fNumberOfCPVPadsZ(0),
  fNCellsXInStrip(0),fNCellsZInStrip(0),fNStripZ(0),
  fCrystalShift(0.),fCryCellShift(0.),fCryStripShift(0.),fCellStep(0.),
  fPadSizePhi(0.),fPadSizeZ(0.),fCPVBoxSizeY(0.),fMisalArray(0x0)
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
  const Float_t * inthermo = fGeometryEMCA->GetInnerThermoHalfSize() ;
  const Float_t * strip    = fGeometryEMCA->GetStripHalfSize() ;
  const Float_t * splate   = fGeometryEMCA->GetSupportPlateHalfSize();
  const Float_t * crystal  = fGeometryEMCA->GetCrystalHalfSize() ;
  const Float_t * pin      = fGeometryEMCA->GetAPDHalfSize() ;
  const Float_t * preamp   = fGeometryEMCA->GetPreampHalfSize() ;
  fCrystalShift=-inthermo[1]+strip[1]+splate[1]+crystal[1]-fGeometryEMCA->GetAirGapLed()/2.+pin[1]+preamp[1] ;
  fCryCellShift=crystal[1]-(fGeometryEMCA->GetAirGapLed()-2*pin[1]-2*preamp[1])/2;
  fCryStripShift=fCryCellShift+splate[1] ;
  fCellStep = 2.*fGeometryEMCA->GetAirCellHalfSize()[0] ;

  fNumberOfCPVPadsPhi = fGeometryCPV->GetNumberOfCPVPadsPhi() ;
  fNumberOfCPVPadsZ   = fGeometryCPV->GetNumberOfCPVPadsZ() ;
  fPadSizePhi = fGeometryCPV->GetCPVPadSizePhi() ;
  fPadSizeZ   = fGeometryCPV->GetCPVPadSizeZ() ; 
  fCPVBoxSizeY= fGeometryCPV->GetCPVBoxSize(1) ;

  for(Int_t mod=0; mod<5; mod++){
    fEMCMatrix[mod]=0 ;
    for(Int_t istrip=0; istrip<224; istrip++)
      fStripMatrix[mod][istrip]=0 ;
    fCPVMatrix[mod]=0;
    fPHOSMatrix[mod]=0 ;
  }
 
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
  if(fMisalArray){
    delete fMisalArray; fMisalArray=0 ;
  }

  for(Int_t mod=0; mod<5; mod++){
     if(fPHOSMatrix[mod]){
      delete fPHOSMatrix[mod];
      fPHOSMatrix[mod]=0x0 ;
    }
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

  if(relid[1]==0){ //this is PHOS

    Double_t pos[3]= {0.0,-fCryCellShift,0.}; //Position incide the crystal 
    Double_t posC[3]={0.0,0.0,0.}; //Global position

    //Shift and possibly apply misalignment corrections
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/fNCellsXInStrip))*fNStripZ-
                (Int_t) TMath::Ceil((Double_t)relid[3]/fNCellsZInStrip) ;
    pos[0]=((relid[2]-1)%fNCellsXInStrip-fNCellsXInStrip/2+0.5)*fCellStep ;
    pos[2]=(-(relid[3]-1)%fNCellsZInStrip+fNCellsZInStrip/2-0.5)*fCellStep ;

    Int_t mod = relid[0] ;
    const TGeoHMatrix * m2 = GetMatrixForStrip(mod, strip) ;
    if(m2)
      m2->LocalToMaster(pos,posC);
    else{ //shold not happen!
      AliError(Form("Can not find matrix for mod=%d, strip=%d",mod, strip)) ;
      //posC contains fixed valued to identify problem in analysis
    }  
    //Return to PHOS local system  
    Double_t posL2[3]={posC[0],posC[1],posC[2]};
    const TGeoHMatrix *mPHOS2 = GetMatrixForModule(mod) ;
    if(mPHOS2){
      mPHOS2->MasterToLocal(posC,posL2);
      x=posL2[0] ;
      z=-posL2[2];
      return ;
    }
    else{
      AliError(Form("Can not find matrix for mod=%d",mod)) ;
      //Return wrong fixed value to notice in analysis 
      x=0. ;
      z=0.;
      return ;
    }
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
    const TGeoHMatrix *m = GetMatrixForCPV(relid[0]) ;
    if(m)
      m->LocalToMaster(pos,posC);
    else{
      AliError(Form("Can not find CPV matrix for mod=%d",relid[0])) ;
      //posC contains fixed valued to identify problem in analysis
    }  
    //Return to PHOS local system
    Double_t posL[3]={0.,0.,0.,} ;
    const TGeoHMatrix *mPHOS = GetMatrixForPHOS(relid[0]) ;
    if(mPHOS){
      mPHOS->MasterToLocal(posC,posL);
      x=posL[0] ;
      z=posL[1];
      return ;
    }
    else{
      AliError(Form("Can not find (CPV) matrix for mod=%d",relid[0])) ;
      //Return wrong fixed value to notice in analysis 
      x=0. ;
      z=0.;
      return ;
    }
  }
  
}
//____________________________________________________________________________
void AliPHOSGeoUtils::GetCrystalsEdges(Int_t mod, Float_t & xmin, Float_t &zmin, Float_t &xmax, Float_t &zmax){
  //return coordinated of crystal matrix edges in local frame
  Int_t relid[4]={mod,0,1,1} ;
  relid[0]=mod ;
  //check if this is 1/2 of the module 
  Bool_t halfMod=kFALSE ;
  if(gGeoManager){
    if (gGeoManager->CheckPath(Form("/ALIC_1/PHOH_%d",mod))){  
      halfMod=kTRUE ;
    }
  }
  else{ //hardcoded!
     halfMod=(mod==4) ;
  }
  if(halfMod)
    relid[2]=33 ;
  RelPosInModule(relid,xmin,zmin) ; //coordinate of the corner cell
  relid[2]=64 ;
  relid[3]=56 ;
  RelPosInModule(relid,xmax,zmax) ; //coordinate of the corner cell
 
}
//____________________________________________________________________________
void AliPHOSGeoUtils::RelPosToAbsId(Int_t module, Double_t x, Double_t z, Int_t & absId) const
{
  // converts local PHOS-module (x, z) coordinates to absId 

  //Calculate AbsId using ideal geometry. Should be sufficient for primary particles calculation
  //(the only place where this method used currently)
  Int_t relid[4]={module,0,1,1} ;
  relid[2] = static_cast<Int_t>(TMath::Ceil( x/ fCellStep + fNPhi / 2.) );
  relid[3] = fNZ+1-static_cast<Int_t>(TMath::Ceil(-z/ fCellStep + fNZ   / 2.) ) ;
  if(relid[2]<1)relid[2]=1 ;
  if(relid[3]<1)relid[3]=1 ;
  if(relid[2]>fNPhi)relid[2]=fNPhi ;
  if(relid[3]>fNZ)relid[3]=fNZ ;
  RelToAbsNumbering(relid,absId) ;

/*
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
*/
 
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
    AliFatal("Geo manager not initialized\n");
  }
    
  Int_t relid[4] ;
    
  AbsToRelNumbering(id , relid) ;
    
  //construct module name
  if(relid[1]==0){ //this is EMC
 
    Double_t ps[3]= {0.0,-fCryStripShift,0.}; //Position incide the crystal
    Double_t psC[3]={0.0,0.0,0.}; //Global position
 
    //Shift and possibly apply misalignment corrections
    Int_t strip=1+((Int_t) TMath::Ceil((Double_t)relid[2]/fNCellsXInStrip))*fNStripZ-
                (Int_t) TMath::Ceil((Double_t)relid[3]/fNCellsZInStrip) ;
    ps[0]=((relid[2]-1)%fNCellsXInStrip-fNCellsXInStrip/2+0.5)*fCellStep ;
    ps[2]=(-(relid[3]-1)%fNCellsZInStrip+fNCellsZInStrip/2-0.5)*fCellStep ;
 
    Int_t mod = relid[0] ;
    const TGeoHMatrix * m2 = GetMatrixForStrip(mod, strip) ;
    if(m2){
      m2->LocalToMaster(ps,psC);
      pos.SetXYZ(psC[0],psC[1],psC[2]) ; 
    }
    else{
      AliError(Form("Can not find matrix for mod=%d, strip=%d",mod,strip)) ;
      //Return wrong fixed value to notice in analysis 
      pos.SetXYZ(0.,0.,0.) ;
    }
  }
  else{
    //first calculate position with respect to CPV plane
    Int_t row        = relid[2] ; //offset along x axis
    Int_t column     = relid[3] ; //offset along z axis
    Double_t ps[3]= {0.0,fCPVBoxSizeY/2.,0.}; //Position on top of CPV
    Double_t psC[3]={0.0,0.0,0.}; //Global position
    pos[0] = - ( fNumberOfCPVPadsPhi/2. - row    - 0.5 ) * fPadSizePhi  ; // position of pad  with respect
    pos[2] = - ( fNumberOfCPVPadsZ  /2. - column - 0.5 ) * fPadSizeZ  ; // of center of PHOS module
 
    //now apply possible shifts and rotations
    const TGeoHMatrix *m = GetMatrixForCPV(relid[0]) ;
    if(m){
      m->LocalToMaster(ps,psC);
      pos.SetXYZ(psC[0],psC[1],-psC[2]) ; 
    }
    else{
      AliError(Form("Can not find (CPV) matrix for mod=%d",relid[0])) ;
      //Return wrong fixed value to notice in analysis 
      pos.SetXYZ(0.,0.,0.) ;
      
    }
  }
} 

//____________________________________________________________________________
void AliPHOSGeoUtils::Local2Global(Int_t mod, Float_t x, Float_t z,
				   TVector3& globalPosition) const 
{
  Double_t posL[3]={x,-fCrystalShift,-z} ; //Only for EMC!!!
  Double_t posG[3] ;
  const TGeoHMatrix *mPHOS = GetMatrixForModule(mod) ;
  mPHOS->LocalToMaster(posL,posG);
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
  const TGeoHMatrix *mPHOS = GetMatrixForModule(module) ;
  if(mPHOS){
    mPHOS->MasterToLocal(posG,posL);
    localPosition.SetXYZ(posL[0],posL[1]+fCrystalShift,-posL[2]) ;  
  }
  else{
    localPosition.SetXYZ(999.,999.,999.) ; //module does not exist in given configuration
  }
 
}
//____________________________________________________________________________
Bool_t AliPHOSGeoUtils::GlobalPos2RelId(TVector3 & global, Int_t * relId){
  //Converts position in global ALICE coordinates to relId 
  //returns false if x,z coordinates are beyond PHOS
  //distande to PHOS surface is NOT calculated 
  TVector3 loc ;
  for(Int_t mod=1; mod<=fNModules; mod++){
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

  for(Int_t imod=1; imod<=fNModules ; imod++){
    //create vector from (0,0,0) to center of crystal surface of imod module
    Double_t tmp[3]={0.,-fCrystalShift,0.} ;
 
    const TGeoHMatrix *m = GetMatrixForModule(imod) ;
    if(!m) //module does not exist in given configuration
      continue ; 
    Double_t posG[3]={0.,0.,0.} ;
    m->LocalToMaster(tmp,posG);
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
//____________________________________________________________________________
const TGeoHMatrix * AliPHOSGeoUtils::GetMatrixForModule(Int_t mod)const {
  //Provides shift-rotation matrix for module mod

  //If GeoManager exists, take matrixes from it
  if(gGeoManager){
    char path[255] ;
    snprintf(path,255,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",mod) ;
    //    sprintf(path,"/ALIC_1/PHOS_%d",relid[0]) ;
    if (!gGeoManager->CheckPath(path)){ //Module with CPV
      snprintf(path,255,"/ALIC_1/PHOC_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",mod) ;
      if (!gGeoManager->CheckPath(path)){
        //Try half-mod name      
        snprintf(path,255,"/ALIC_1/PHOH_%d/PEMH_1/PCLH_1/PIOH_1/PCOH_1/PAGH_1/PTIH_1",mod) ;
        if (!gGeoManager->CheckPath(path)){
//          AliWarning(Form("Geo manager can not find path %s \n",path));
          return 0;
	}
      }
    }
    gGeoManager->cd(path) ;
    return gGeoManager->GetCurrentMatrix();
  }
  if(fEMCMatrix[mod-1]){
    return fEMCMatrix[mod-1] ;
  }
  else{
 //   AliWarning("Can not find PHOS misalignment matrixes\n") ;
 //   AliWarning("Either import TGeoManager from geometry.root or \n");
 //   AliWarning("read stored matrixes from AliESD Header: \n") ;
 //   AliWarning("AliPHOSGeoUtils::SetMisalMatrixes(header->GetPHOSMisalMatrix()) \n") ; 
    return 0 ;
  }
  return 0 ;
}
//____________________________________________________________________________
const TGeoHMatrix * AliPHOSGeoUtils::GetMatrixForStrip(Int_t mod, Int_t strip)const {
  //Provides shift-rotation matrix for strip unit of the module mod

  //If GeoManager exists, take matrixes from it
  if(gGeoManager){
    char path[255] ;
    snprintf(path,255,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d",mod,strip) ;
    if (!gGeoManager->CheckPath(path)){ //Test module with CPV
      snprintf(path,255,"/ALIC_1/PHOC_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1/PSTR_%d",mod,strip) ;
      if (!gGeoManager->CheckPath(path)){
        //Look for half-module path
        snprintf(path,255,"/ALIC_1/PHOH_%d/PEMH_1/PCLH_1/PIOH_1/PCOH_1/PAGH_1/PTIH_1/PSTR_%d",mod,strip) ;
        if (!gGeoManager->CheckPath(path)){    
//          AliWarning(Form("Geo manager can not find path %s \n",path));
          return 0 ;
	}
      }
    }
    gGeoManager->cd(path) ;
    return gGeoManager->GetCurrentMatrix();
  } 
  if(fStripMatrix[mod-1][strip-1]){
    return fStripMatrix[mod-1][strip-1] ;
  }
  else{
    AliWarning("Can not find PHOS misalignment matrixes\n") ;
    AliWarning("Either import TGeoManager from geometry.root or \n");
    AliWarning("read stored matrixes from AliESD Header: \n") ; 
    AliWarning("AliPHOSGeoUtils::SetMisalMatrixes(header->GetPHOSMisalMatrix()) \n") ;
    return 0 ;
  } 
  return 0 ;
}
//____________________________________________________________________________
const TGeoHMatrix * AliPHOSGeoUtils::GetMatrixForCPV(Int_t mod)const {
  //Provides shift-rotation matrix for CPV of the module mod

  //If GeoManager exists, take matrixes from it
  if(gGeoManager){ 
    char path[255] ;
    //now apply possible shifts and rotations
    snprintf(path,255,"/ALIC_1/PHOC_%d/PCPV_1",mod) ;
    if (!gGeoManager->CheckPath(path)){
      snprintf(path,255,"/ALIC_1/PHOH_%d/PCPV_1",mod) ;
      if (!gGeoManager->CheckPath(path)){
//        AliWarning(Form("Geo manager can not find path %s \n",path));
        return 0 ;
      }
    }
    gGeoManager->cd(path) ;
    return gGeoManager->GetCurrentMatrix();
  }
  if(fCPVMatrix[mod-1]){
    return fCPVMatrix[mod-1] ;
  }
  else{
    AliWarning("Can not find PHOS misalignment matrixes\n") ;
    AliWarning("Either import TGeoManager from geometry.root or \n");
    AliWarning("read stored matrixes from AliESD Header: \n") ;  
    AliWarning("AliPHOSGeoUtils::SetMisalMatrixes(header->GetPHOSMisalMatrix()) \n") ;
    return 0 ;
  }
  return 0 ;
} 
//____________________________________________________________________________
const TGeoHMatrix * AliPHOSGeoUtils::GetMatrixForPHOS(Int_t mod)const {
  //Provides shift-rotation matrix for PHOS (EMC+CPV) 

  //If manually set matrises exist, use them rather than TGeoManager (analysis case)  
  if(fPHOSMatrix[mod-1]){
    return fPHOSMatrix[mod-1] ;
  }    
    
  //If GeoManager exists, take matrixes from it
  if(gGeoManager){

    char path[255] ;
    snprintf(path,255,"/ALIC_1/PHOS_%d",mod) ;
    if (!gGeoManager->CheckPath(path)){ //Module with CPV
      snprintf(path,255,"/ALIC_1/PHOC_%d",mod) ;
      if (!gGeoManager->CheckPath(path)){ //1/2 module
        snprintf(path,255,"/ALIC_1/PHOH_%d",mod) ;
        if (!gGeoManager->CheckPath(path)){
//          AliWarning(Form("Geo manager can not find path %s \n",path));
          return 0 ;
	}
      }
    }
    gGeoManager->cd(path) ;
    return gGeoManager->GetCurrentMatrix();
  }
  else{
    AliWarning("Can not find PHOS misalignment matrixes\n") ;
    AliWarning("Either import TGeoManager from geometry.root or \n");
    AliWarning("read stored matrixes from AliESD Header:  \n") ;   
    AliWarning("AliPHOSGeoUtils::SetMisalMatrixes(header->GetPHOSMisalMatrix()) \n") ;
    return 0 ;
  }
  return 0 ;
}
//____________________________________________________________________________
void AliPHOSGeoUtils::SetMisalMatrix(const TGeoHMatrix * m, Int_t mod){
  //Fills pointers to geo matrixes
 
  if(fPHOSMatrix[mod]){ //have been set already. Can not be changed any more
    return ;
  }
  if(m==NULL) //Matrix for non-existing modules? Remain zero, no need to re-set
    return ;
  fPHOSMatrix[mod]= new TGeoHMatrix(*m) ;
  
  //Calculate maxtrices for PTII
  if(!fMisalArray)
    fMisalArray = new TClonesArray("TGeoHMatrix",1120+10) ;
  Int_t nr = fMisalArray->GetEntriesFast() ;
  Double_t rotEMC[9]={1.,0.,0.,0.,0.,-1.,0.,1.,0.} ;
  const Float_t * inthermo = fGeometryEMCA->GetInnerThermoHalfSize() ;
  const Float_t * strip    = fGeometryEMCA->GetStripHalfSize() ;
  const Float_t * covparams = fGeometryEMCA->GetAlCoverParams() ;
  const Float_t * warmcov = fGeometryEMCA->GetWarmAlCoverHalfSize() ;
  Float_t z = fGeometryCPV->GetCPVBoxSize(1) / 2. - warmcov[2] + covparams[3]-inthermo[1] ;
  Double_t locTII[3]={0.,0.,z} ; 
  Double_t globTII[3] ;

  if (fEMCMatrix[mod] == NULL)
    fEMCMatrix[mod] = new((*fMisalArray)[nr])TGeoHMatrix() ;
  nr++ ;
  fEMCMatrix[mod]->SetRotation(rotEMC) ;
  fEMCMatrix[mod]->MultiplyLeft(fPHOSMatrix[mod]) ;
  fPHOSMatrix[mod]->LocalToMaster(locTII,globTII) ;
  fEMCMatrix[mod]->SetTranslation(globTII) ;
 
  //Now calculate ideal matrixes for strip misalignment.
  //For the moment we can not store them in ESDHeader

  Double_t loc[3]={0.,inthermo[1] - strip[1],0.} ; 
  Double_t glob[3] ;

  Int_t istrip=0 ;
  for(Int_t irow = 0; irow < fGeometryEMCA->GetNStripX(); irow ++){
    loc[0] = (2*irow + 1 - fGeometryEMCA->GetNStripX())* strip[0] ;
    for(Int_t icol = 0; icol < fGeometryEMCA->GetNStripZ(); icol ++){
      loc[2] = (2*icol + 1 - fGeometryEMCA->GetNStripZ()) * strip[2] ;
      fEMCMatrix[mod]->LocalToMaster(loc,glob) ;
      if (fStripMatrix[mod][istrip] == NULL)
	fStripMatrix[mod][istrip] = new((*fMisalArray)[nr])TGeoHMatrix(*(fEMCMatrix[mod])) ; //Use same rotation as PHOS module
      nr++ ;
      fStripMatrix[mod][istrip]->SetTranslation(glob) ;
      istrip++;
    }
  }
 
  //Now calculate CPV matrixes
  const Float_t * emcParams = fGeometryEMCA->GetEMCParams() ;
  Double_t globCPV[3] ;
  Double_t locCPV[3]={0.,0.,- emcParams[3]} ;
  Double_t rot[9]={1.,0.,0.,0.,0.,1.,0.,-1.,0.} ;

  if (fCPVMatrix[mod] == NULL)
    fCPVMatrix[mod] = new((*fMisalArray)[nr])TGeoHMatrix() ;
  nr++ ;
  fCPVMatrix[mod]->SetRotation(rot) ;
  fCPVMatrix[mod]->MultiplyLeft(fPHOSMatrix[mod]) ;
  fCPVMatrix[mod]->ReflectY(kFALSE) ;
  fPHOSMatrix[mod]->LocalToMaster(locCPV,globCPV) ;
  fCPVMatrix[mod]->SetTranslation(globCPV) ;

}
//____________________________________________________________________________ 
void AliPHOSGeoUtils::TestSurvey(Int_t module, const Float_t *point, TVector3 &globaPos) const { 
  //method used in PHOS alignment check
  //Input is module number and point is Photogrammetry reference point wrt top right crystal
  //output- point coordinates in ALICE global system
  
  Double_t x0=31.5*fCellStep;      //Number of crystals 
  Double_t z0=26.5*fCellStep;  //from module center
  Double_t posL[3]={-x0+point[0],-fCrystalShift-point[1],z0-point[2]} ; 
  Double_t posG[3] ;
  const TGeoHMatrix *mPHOS = GetMatrixForModule(module) ;
  if(mPHOS){
    mPHOS->LocalToMaster(posL,posG);
    globaPos.SetXYZ(posG[0],posG[1],posG[2]) ;
  }
  else{
     AliError(Form("Can not find matrix for mod=%d",module)) ;
      //Return wrong fixed value to notice in analysis 
      globaPos.SetXYZ(0.,0.,0.) ;
  }
}
