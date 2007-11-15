/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSCLUSTERIZER_H
#define ALIHLTPHOSCLUSTERIZER_H

#include "AliHLTPHOSBase.h"
#include "AliPHOSGetter.h"

#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"

#include "AliPHOSGeometry.h"

class TClonesArray;
class AliPHOSDigit;
class AliPHOSRecoParamEmc;

class AliHLTPHOSClusterizer : public AliHLTPHOSBase
{
  
public:
  
  AliHLTPHOSClusterizer();    
  
  virtual ~AliHLTPHOSClusterizer();
  
  AliHLTPHOSClusterizer(const AliHLTPHOSClusterizer &);
  AliHLTPHOSClusterizer & operator = (const AliHLTPHOSClusterizer &) {return *this;}

  void SetRecPointContainer(AliHLTPHOSRecPointContainerStruct *RecPointContainerPtr)
  { fRecPointContainerPtr = RecPointContainerPtr; }

  void SetRecoParameters(AliPHOSRecoParamEmc* recoPars);

  void SetEmcClusteringThreshold(Float_t threshold) { fEmcClusteringThreshold = threshold; }
  void SetEmcMinEnergyThreshold(Float_t threshold) { fEmcMinEnergyThreshold = threshold; }
  void SetEmcTimeGate(Float_t gate) { fEmcTimeGate = gate; }
  void SetLogWeight(Float_t weight) { fLogWeight = weight; }  
    
  void SetOfflineMode(AliPHOSGetter* getter); 
  
  virtual Int_t ClusterizeEvent();
  virtual Int_t GetEvent(Int_t evtNr);
  
  Int_t GetNEvents();

  virtual void ScanForNeighbourDigits(Int_t digIndex, AliHLTPHOSRecPointDataStruct* recPoint);
  virtual Int_t AreNeighbours(AliHLTPHOSDigitDataStruct* d1, AliHLTPHOSDigitDataStruct* d2);
  virtual void CalculateCenterOfGravity();

private:
  
  Float_t fEmcClusteringThreshold;  //comment
  Float_t fEmcMinEnergyThreshold; //comment
  Float_t fEmcTimeGate; //comment
  Float_t fLogWeight; //comment
  Int_t fDigitsInCluster; //comment

  Bool_t fOnlineMode; //comment
  
  TClonesArray *fDigitArrayPtr; //comment 
  TObjArray *fEmcRecPointsPtr; //comment
  AliPHOSDigit *fDigitPtr; //comment

  AliHLTPHOSDigitContainerDataStruct *fDigitContainerPtr; //comment
  AliHLTPHOSRecPointContainerStruct *fRecPointContainerPtr; //comment
  AliPHOSGeometry *fPHOSGeometry; //comment
  
  AliPHOSGetter *fGetterPtr; //comment
  
  ClassDef(AliHLTPHOSClusterizer, 1);
};

#endif
