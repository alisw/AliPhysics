/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTPHOSClusterizer.h
    @author Ãystein Djuvsland
    @date   
    @brief  A temporary clusterizer for PHOS
*/

#ifndef ALIHLTPHOSCLUSTERIZER_H
#define ALIHLTPHOSCLUSTERIZER_H

//#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;

struct AliHLTPHOSClusterDataStruct;
struct AliHLTPHOSRecPointDataStruct;
struct AliHLTPHOSValidCellDataStruct;
struct AliHLTPHOSRecPointListDataStruct;
struct AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSClusterizer
{
  
 public: 
  
  AliHLTPHOSClusterizer();
  virtual ~AliHLTPHOSClusterizer();
  AliHLTPHOSClusterizer(const AliHLTPHOSClusterizer &);
  AliHLTPHOSClusterizer & operator = (const AliHLTPHOSClusterizer &) {return *this;}
   
  void    SetThreshold(Float_t threshold) {fThreshold = threshold;}
  void    SetClusterThreshold(Float_t clusterThreshold) {fClusterThreshold = clusterThreshold;}
  
  void    SetHighGainFactor(Float_t highGain) {fHighGainFactor = highGain;}
  void    SetLowGainFactor(Float_t lowGain) {fLowGainFactor = lowGain;}
  void    SetArraySize(Int_t size) 
  { 
    fArraySize = size;
    fMultiplicity = fArraySize * fArraySize;
  }
  
  Float_t GetHighGainFactor() {return fHighGainFactor;}
  Float_t GetLowGainFactor() {return fLowGainFactor;}
  
  Int_t   BuildCellEnergyArray(AliHLTPHOSRcuCellEnergyDataStruct *structPtr, AliHLTPHOSRecPointListDataStruct* recPointList);
  Int_t   CreateRecPointStructArray(AliHLTPHOSRecPointDataStruct* rectStructsPtr, AliHLTPHOSRecPointListDataStruct* list, Int_t nPoints);
  Int_t   CalculateCenterOfGravity(AliHLTPHOSRecPointDataStruct* recPointPtr);
  Int_t   ClusterizeStruct(AliHLTPHOSRecPointDataStruct* recArrayPtr, AliHLTPHOSClusterDataStruct* clusterArrayPtr);
  Int_t   ResetCellEnergyArray();

  
 private:

  AliHLTUInt8_t fPHOSModule;                                     /**<Number of the PHOSModule*/
  Float_t fEnergyArray[N_COLUMNS_MOD][N_ROWS_MOD];               /**<2D array of cell energies*/
  Float_t fThreshold;                                            /**<Energy threshold*/
  Float_t fClusterThreshold;                                     /**<Cluster threshold*/
  Float_t fHighGainFactor;                                       /**<High gain factor*/
  Float_t fLowGainFactor;                                        /**<Low gain factor*/
  Int_t   fArraySize;                                            /**<Size of the array which the energies are summed*/
  Int_t   fMultiplicity;                                         /**<Number of crystals the energies are summed for*/

  ClassDef(AliHLTPHOSClusterizer, 1);
};

#endif
