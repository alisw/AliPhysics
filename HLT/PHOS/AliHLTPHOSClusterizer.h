/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
/** @file   AliHLTPHOSClusterizer.h
    @author Ãystein Djuvsland
    @date   
    @brief  A temporary clusterizer for PHOS
*/


#ifndef ALIHLTPHOSCLUSTERIZER_H
#define ALIHLTPHOSCLUSTERIZER_H

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliPHOSGeometry.h"
#include "TClonesArray.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGetter.h"
#include "AliPHOSRecoParamEmc.h"

class AliHLTPHOSClusterizer : public AliHLTPHOSBase
{
  
public:
  
  AliHLTPHOSClusterizer();    
  
  virtual ~AliHLTPHOSClusterizer();
  
  AliHLTPHOSClusterizer(const AliHLTPHOSClusterizer &);
  AliHLTPHOSClusterizer & operator = (const AliHLTPHOSClusterizer &) {return *this;}

  void SetRecPointContainer(AliHLTPHOSRecPointContainerStruct *RecPointContainerPtr)
  { fRecPointContainerPtr = RecPointContainerPtr; }

  void SetRecoParameters(AliPHOSRecoParamEmc*);

  void SetEmcClusteringThreshold(Float_t threshold) { fEmcClusteringThreshold = threshold; }
  void SetEmcMinEnergyThreshold(Float_t threshold) { fEmcMinEnergyThreshold = threshold; }
  void SetEmcTimeGate(Float_t gate) { fEmcTimeGate = gate; }
  void SetLogWeight(Float_t weight) { fLogWeight = weight; }  
    
  void SetOfflineMode(AliPHOSGetter*); 
  
  virtual Int_t ClusterizeEvent();
  virtual Int_t GetEvent(Int_t);
  
  Int_t GetNEvents();

  virtual void ScanForNeighbourDigits(Int_t, AliHLTPHOSRecPointDataStruct*);
  virtual Int_t AreNeighbours(AliHLTPHOSDigitDataStruct*, AliHLTPHOSDigitDataStruct*);
  virtual void CalculateCenterOfGravity();

private:
  
  Float_t fEmcClusteringThreshold;
  Float_t fEmcMinEnergyThreshold;
  Float_t fEmcTimeGate;
  Float_t fLogWeight;
  Int_t fDigitsInCluster;

  Bool_t fOnlineMode;
 
  TClonesArray *fDigitArrayPtr; 
  TObjArray *fEmcRecPointsPtr;
  AliPHOSDigit *fDigitPtr;

  AliHLTPHOSDigitContainerDataStruct *fDigitContainerPtr;
  AliHLTPHOSRecPointContainerStruct *fRecPointContainerPtr;
  AliPHOSGeometry *fPHOSGeometry;
  
  AliPHOSGetter *fGetterPtr;
  
  ClassDef(AliHLTPHOSClusterizer, 1);
};

#endif
