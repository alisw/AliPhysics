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
//  AliAODCaloCluster extension for PHOS to recalculate cluster 
//  parameters in case of recalibration.
//*--
//*-- Author: Dmitri Peressounko (RRC KI)


// --- ROOT system ---
#include "TVector3.h"
#include "TMath.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h" 
#include "AliPHOSGeometry.h" 
#include "AliPHOSPIDv1.h" 
#include "AliPHOSReconstructor.h" 
#include "AliPHOSAodCluster.h" 
#include "AliPHOSCalibData.h"
#include "AliAODCaloCells.h"

ClassImp(AliPHOSAodCluster)

//____________________________________________________________________________
AliPHOSAodCluster::AliPHOSAodCluster() : 
  AliAODCaloCluster(),fRecalibrated(0)
{
  // ctor
}
//____________________________________________________________________________
AliPHOSAodCluster::AliPHOSAodCluster(const AliAODCaloCluster & clu) : 
  AliAODCaloCluster(clu),fRecalibrated(0)
{
  // cpy ctor
}

//____________________________________________________________________________
AliPHOSAodCluster::~AliPHOSAodCluster()
{
  // dtor
}
//____________________________________________________________________________
void AliPHOSAodCluster::Recalibrate(AliPHOSCalibData * calibData,AliAODCaloCells *phsCells){
  //If not done yet, apply recalibration coefficients to energies list
  //NOTE that after recalibration fCellsAmpFraction contains not FRACTION but FULL energy 
  
  if(fRecalibrated)
   return ;
  
  if(!calibData)
    return ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;
  if(!phosgeom)
    AliFatal("AliPHOSGeometry was not contructed\n") ;
	
  Double32_t * cellsAmpFraction = GetCellsAmplitudeFraction(); 
	
  for(Int_t i=0; i < GetNCells(); i++){
    Int_t relId[4];
    phosgeom->AbsToRelNumbering(GetCellAbsId(i),relId) ;
    Int_t   module = relId[0];
    Int_t   column = relId[3];
    Int_t   row    = relId[2];
    Double_t energy = phsCells->GetCellAmplitude(GetCellAbsId(i)) ;
    cellsAmpFraction[i]*=energy*calibData->GetADCchannelEmc(module,column,row);
  }
	
  SetCellsAmplitudeFraction(cellsAmpFraction);
  fRecalibrated=kTRUE; 
}
//____________________________________________________________________________
void  AliPHOSAodCluster::EvalAll(Float_t logWeight, TVector3 &vtx){
    //If recalibrated - recalculate all cluster parameters
  if(!fRecalibrated)
    return ;

  EvalEnergy() ; //Energy should be evaluated first
  EvalCoord(logWeight, vtx) ;
  
}
//____________________________________________________________________________
void AliPHOSAodCluster::EvalEnergy(){
  if(!fRecalibrated) // no need to recalibrate
    return ;
    
  Float_t energy=0. ;
  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
    energy+=GetCellAmplitudeFraction(iDigit) ;
  }
  
  SetE(energy);
  //Correct for nonlinearity later
   
}
//____________________________________________________________________________
void AliPHOSAodCluster::EnergyCorrection(){
  //apply nonlinearity correction same as in AliPHOSPIDv1.
  SetE(AliPHOSReconstructor::CorrectNonlinearity(E())) ;
}
//____________________________________________________________________________
void AliPHOSAodCluster::EvalPID(AliPHOSPIDv1 * /*pid*/){           
  //re-evaluate identification parameters
//  pid->CalculatePID(E(),GetDispersion(),GetEmcCpvDistance(),GetTOF(),GetPID()) ;  
//  pid->CalculatePID(E(),GetDispersion(),GetM20(),GetM02(),GetEmcCpvDistance(),GetTOF(),GetPID()) ;
}
//____________________________________________________________________________
void AliPHOSAodCluster::EvalCoord(Float_t logWeight, TVector3 &vtx)
{
  // Calculates new center of gravity in the local PHOS-module coordinates 
  // and tranfers into global ALICE coordinates
  // Calculates Dispersion and main axis
  if(!fRecalibrated) // no need to recalibrate
    return ;
 
  Float_t wtot = 0. ;
  Int_t relid[4] ;
  Int_t phosMod=0 ;
  Float_t xMean = 0. ;
  Float_t zMean = 0. ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;
  if(!phosgeom)
    AliFatal("AliPHOSGeometry was not contructed\n") ;

  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(GetCellAbsId(iDigit), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    phosMod=relid[0] ;
    Double_t ei=GetCellAmplitudeFraction(iDigit) ;
    if (E() > 0 && ei > 0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log(ei/E()) ) ;
      xMean+= xi * w ;
      zMean+= zi * w ;
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", ei, E()));
  }
  if (wtot>0) {
    xMean /= wtot ;
    zMean /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));


// Calculates the dispersion and second momenta
  Double_t d=0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  for(Int_t iDigit=0; iDigit < GetNCells(); iDigit++) {
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(GetCellAbsId(iDigit), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    Double_t ei=GetCellAmplitudeFraction(iDigit) ;
    if (E() > 0 && ei > 0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log(ei/E()) ) ;
      d += w*((xi-xMean)*(xi-xMean) + (zi-zMean)*(zi-zMean) ) ; 
      dxx  += w * xi * xi ;
      dzz  += w * zi * zi ;
      dxz  += w * xi * zi ; 
   }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", ei, E()));
  }
  
  if (wtot>0) {
    d /= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xMean * xMean ;
    dzz -= zMean * zMean ;
    dxz -= xMean * zMean ;
    SetM02(0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )) ;
    SetM20(0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ));
  }
  else{
    AliError(Form("Wrong weight %f\n", wtot));
    d=0. ;
    SetM20(0.) ;
    SetM02(0.) ;
  }

  if (d>=0)
    SetDispersion(TMath::Sqrt(d)) ;
  else    
    SetDispersion(0) ;


  // Correction for the depth of the shower starting point (TDR p 127)  
  Float_t para = 0.925 ; 
  Float_t parb = 6.52 ; 

  TVector3 vInc ;
  phosgeom->GetIncidentVector(vtx,phosMod,xMean,zMean,vInc) ;

  Float_t depthx = 0.; 
  Float_t depthz = 0.;
  if (E()>0&&vInc.Y()!=0.) {
    depthx = ( para * TMath::Log(E()) + parb ) * vInc.X()/TMath::Abs(vInc.Y()) ;
    depthz = ( para * TMath::Log(E()) + parb ) * vInc.Z()/TMath::Abs(vInc.Y()) ;
  }
  else 
    AliError(Form("Wrong amplitude %f\n", E()));

  xMean-= depthx  ;
  zMean-= depthz  ;

  //Go to the global system
  TVector3 gps ;
  phosgeom->Local2Global(phosMod, xMean, zMean, gps) ;
  SetPositionAt(gps[0],0) ;
  SetPositionAt(gps[1],1) ;  
  SetPositionAt(gps[2],2) ;
}
