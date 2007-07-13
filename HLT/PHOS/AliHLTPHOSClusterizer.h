/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
/** @file   AliHLTPHOSClusterizer.h
    @author Ãystein Djuvsland
    @date   
    @brief  A temporary clusterizer for PHOS
*/


#ifndef ALIHLTPHOSCLUSTERIZER_H
#define ALIHLTPHOSCLUSTERIZER_H

#include "AliHLTPHOSCommonDefs.h"
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
   
  void    SetThreshold(float threshold) {fThreshold = threshold;}
  void    SetClusterThreshold(float clusterThreshold) {fClusterThreshold = clusterThreshold;}
  
  void    SetHighGainFactor(float highGain) {fHighGainFactor = highGain;}
  void    SetLowGainFactor(float lowGain) {fLowGainFactor = lowGain;}
  void    SetArraySize(int size) 
  { 
    fArraySize = size;
    fMultiplicity = fArraySize * fArraySize;
  }
  float GetThreshold() {return fThreshold;}
  float GetClusterThreshold() {return fClusterThreshold;}
  float GetHighGainFactor() {return fHighGainFactor;}
  float GetLowGainFactor() {return fLowGainFactor;}
  float GetArraySize() {return fArraySize;}
  float GetMultiplicity() {return fMultiplicity;}
  
  int   BuildCellEnergyArray(AliHLTPHOSRcuCellEnergyDataStruct *structPtr, AliHLTPHOSRecPointListDataStruct* recPointList);
  int   CreateRecPointStructArray(AliHLTPHOSRecPointDataStruct* rectStructsPtr, AliHLTPHOSRecPointListDataStruct* list, int nPoints);
  int   CalculateCenterOfGravity(AliHLTPHOSRecPointDataStruct* recPointPtr);
  int   CalculateMoments(AliHLTPHOSRecPointDataStruct* recPointPtr, Bool_t axisOnly);
  int   ClusterizeStruct(AliHLTPHOSRecPointDataStruct* recArrayPtr, AliHLTPHOSClusterDataStruct* clusterArrayPtr);
  int   ResetCellEnergyArray();

  
 private:

  AliHLTUInt8_t fPHOSModule;                                     /**<Number of the PHOSModule*/
  Float_t fEnergyArray[N_XCOLUMNS_MOD][N_ZROWS_MOD];             /**<2D array of cell energies*/
  Float_t fThreshold;                                            /**<Energy threshold*/
  Float_t fClusterThreshold;                                     /**<Cluster threshold*/
  Float_t fHighGainFactor;                                       /**<High gain factor*/
  Float_t fLowGainFactor;                                        /**<Low gain factor*/
  Int_t   fArraySize;                                            /**<Size of the array which the energies are summed*/
  Int_t   fMultiplicity;                                         /**<Number of crystals the energies are summed for*/

  ClassDef(AliHLTPHOSClusterizer, 1);
};

#endif
