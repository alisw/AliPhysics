#ifndef ALIPHOSESDCLUSTER_H
#define ALIPHOSESDCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  AliESDCaloCluster version for PHOS (used for recalibration)
//           
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSPIDv1 ;
class AliPHOSCalibData ;
class AliESDCaloCells ;

#include "AliESDCaloCluster.h"

class AliPHOSEsdCluster : public AliESDCaloCluster  {

public:

  AliPHOSEsdCluster() ;
  AliPHOSEsdCluster(const AliESDCaloCluster & clu) ; 
 
  virtual ~AliPHOSEsdCluster() ;  

  void  EvalAll(Float_t logWeight, TVector3 &vtx) ; //re-calculate all cluster parameters
  void  Recalibrate(AliPHOSCalibData * calibData,AliESDCaloCells *phsCells) ; //Apply recalibration to this cluster
  void  EnergyCorrection() ;  //Apply non-linearity correction
  void  EvalPID(AliPHOSPIDv1 * pid) ;           //re-evaluate identification parameters

protected:
 
  void EvalCoord(Float_t logWeight, TVector3 &vtx) ;  //calculate coordinate-related parameters (position, dispersion)
  void EvalEnergy() ; //re-calculate energy of the cluster

  Bool_t fRecalibrated ;  //Have this cluster been recalibrated
    
  ClassDef(AliPHOSEsdCluster,3)  // (PHOS ESD cluster)

};

#endif // AliPHOSESDCLUSTER_H
