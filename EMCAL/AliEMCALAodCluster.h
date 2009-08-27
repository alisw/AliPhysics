#ifndef ALIEMCALAODCLUSTER_H
#define ALIEMCALAODCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  AliAodCaloCluster version for EMCAL (used for recalibration)
//  Copy-paste from methods in AliEMCALRecPoint.       
//
//*-- Author: Dmitri Peressounko (RRC KI) for PHOS
//*-- Adapted for EMCAL: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TVector3;
// --- Standard library ---

// --- AliRoot header files ---
//class AliEMCALPID ;
class AliEMCALCalibData ;
class AliAODCaloCells ;
#include "AliAODCaloCluster.h"

class AliEMCALAodCluster : public AliAODCaloCluster  {

public:

  AliEMCALAodCluster() ;
  AliEMCALAodCluster(const AliAODCaloCluster & clu) ; 
 
  virtual ~AliEMCALAodCluster() ;  

  void  EvalAll(Float_t logWeight, TString geoname) ; //re-calculate all cluster parameters
  void  Recalibrate(AliEMCALCalibData * calibData, AliAODCaloCells *phsCells, TString geoname) ; //Apply recalibration to this cluster
//  void  EnergyCorrection(AliEMCALPID * pid) ;  //Apply non-linearity correction
  void  EvalPID() ;           //re-evaluate identification parameters
	
protected:
	
  Double_t TmaxInCm(const Double_t e , const Int_t key) const ; //Cluster max depth used in EvalPositionAndShowerShape
  void EvalPositionAndShowerShape(Float_t logWeight, TString geoname) ;  //calculate coordinate-related parameters (position, dispersion)
  void EvalEnergy() ; //re-calculate energy of the cluster

  Bool_t fRecalibrated ;  //Has this cluster been recalibrated?
	
  ClassDef(AliEMCALAodCluster,1)  // (EMCAL AOD cluster)

};

#endif // AliEMCALAODCLUSTER_H
