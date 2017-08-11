#ifndef ALICALOTRIGGERSIMULATOR_H
#define ALICALOTRIGGERSIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */
 
//_________________________________________________________________________
// Class to fill two-photon invariant mass hisograms
// to be used to extract pi0 raw yield.
//
//-- Author: Dmitri Peressounko (RRC "KI")
// This class contains all (minimal) necessary information about photon to 
// calculate invarint mass distr for pi0
// and for tagging and isolation analysis

class AliVCluster;

#include "TLorentzVector.h"
#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloTrigger.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEmcCalibData.h"
#include "AliPHOSCalibData.h"
#include "AliCDBStorage.h"
#include "TRandom.h"

class AliCaloTriggerSimulator :public TLorentzVector{
  
 public:
  
  AliCaloTriggerSimulator() ;
  AliCaloTriggerSimulator(AliAODEvent* fEvent) ;
  ~AliCaloTriggerSimulator(){} 
 
  AliAODCaloTrigger* CreateTriggerMap();
    void SetThreshold(Double_t thre){fThreshold = thre;}
  void SetDB(AliPHOSGeometry* fPHOSGeo,AliPHOSEmcCalibData* fCalibDataEmc, AliPHOSCalibData* fPHOSCalibData, AliCDBStorage *fCDBstorage);

 private:
  AliCaloTriggerSimulator(const AliCaloTriggerSimulator&); // not implemented
  AliCaloTriggerSimulator& operator=(const AliCaloTriggerSimulator&);
  
  void SetEvent(AliAODEvent* fEvent);
  
  Int_t GetTRUNum(Int_t cellX, Int_t cellZ);

  AliAODEvent         *fAOD;
  AliAODCaloCluster   *fCluster;
  AliAODCaloCells     *fCells;
  AliPHOSGeometry     *fPHOSGeo;
  AliPHOSEmcCalibData *fCalibDataEmc;
  AliPHOSCalibData    *fPHOSCalibData;
  AliCDBStorage       *fCDBstorage;

  Double_t             fThreshold;

  ClassDef(AliCaloTriggerSimulator,7);
};

#endif // #ifdef ALICALOPHOTON_H

  
