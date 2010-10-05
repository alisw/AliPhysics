#ifndef ALIPHOSAODCLUSTER_H
#define ALIPHOSAODCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  AliAodCaloCluster version for PHOS (used for recalibration)
//           
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSPIDv1 ;
class AliPHOSCalibData ;
class AliAODCaloCells ;

#include "AliAODCaloCluster.h"

class AliPHOSAodCluster : public AliAODCaloCluster  {

public:

  AliPHOSAodCluster() ;
  AliPHOSAodCluster(const AliAODCaloCluster & clu) ; 
 
  virtual ~AliPHOSAodCluster() ;  

  void  EvalAll(Float_t logWeight, TVector3 &vtx) ; //re-calculate all cluster parameters
  void  Recalibrate(AliPHOSCalibData * calibData,AliAODCaloCells *phsCells) ; //Apply recalibration to this cluster
  void  EnergyCorrection() ;  //Apply non-linearity correction
  void  EvalPID(AliPHOSPIDv1 * pid) ;           //re-evaluate identification parameters

protected:
 
  void EvalCoord(Float_t logWeight, TVector3 &vtx) ;  //calculate coordinate-related parameters (position, dispersion)
  void EvalEnergy() ; //re-calculate energy of the cluster

  Bool_t fRecalibrated ;  //Have this cluster been recalibrated
    
  ClassDef(AliPHOSAodCluster,1)  // (PHOS AOD cluster)

};

#endif // AliPHOSAODCLUSTER_H
