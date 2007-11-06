// 1
// 2
// 3
// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTPHOSDEFINITIONS_H
#define ALIHLTPHOSDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTDataTypes.h"

class AliHLTPHOSDefinitions
{
public:
  static const AliHLTComponentDataType fgkCellEnergyDataType;    /**<Reconstructed cell/crystal energies*/
  static const AliHLTComponentDataType fgkDDLPackedRawDataType;  /**<DDL raw data on the RCU data format*/
  static const AliHLTComponentDataType fgkCellEnergyHistogramDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellAverageEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellAccumulatedEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellTimingHistogramDataType;  /**<Histogram for per cell/gain time distribution*/      
  static const AliHLTComponentDataType fgkCellTimingAverageDataType;  /**<Histogram for per cell/gain time distribution*/  
  static const AliHLTComponentDataType fgkCellChannelDataDataType;  /**<Time dependent signal from the readout channels*/  
  static const AliHLTComponentDataType fgkAliHLTClusterDataType;  //Cluster data type
  static const AliHLTComponentDataType fgkAliHLTHistDataType;     //hist data type
  static const AliHLTComponentDataType fgkAliHLTSpectrumDataType; //spectrum data type
  static const AliHLTComponentDataType fgkAliHLTDigitDataType; //Digit data type
  static const AliHLTComponentDataType fgkAliHLTRootTreeDataType; //Root tree type
  static const AliHLTComponentDataType fgkAliHLTBaselineDataType; //Baseline type
  static const AliHLTComponentDataType fgkAliHLTMIPDataType; //"MIP" data type
  static const AliHLTComponentDataType fgkAliHLTNoiseMapDataType; //Noise map data type
  static const AliHLTComponentDataType fgkAliHLTSandboxDataType; //General data type
};

#endif
