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
  static const AliHLTComponentDataType fgkPhosHistDataType;    /**Fourier transform of time dependent signals*/
  static const AliHLTComponentDataType fgkFourierTransform;    /**Fourier transform of time dependent signals*/
  static const AliHLTComponentDataType fgkChannelDataType;    /**<Reconstructed channels*/
  static const AliHLTComponentDataType fgkCellEnergyDataType;    /**<Reconstructed cell/crystal energies*/
  static const AliHLTComponentDataType fgkDDLPackedRawDataType;  /**<DDL raw data on the RCU data format*/
  static const AliHLTComponentDataType fgkCellEnergyHistogramDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellAverageEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellAccumulatedEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
  static const AliHLTComponentDataType fgkCellTimingHistogramDataType;  /**<Histogram for per cell/gain time distribution*/      
  static const AliHLTComponentDataType fgkCellTimingAverageDataType;  /**<Histogram for per cell/gain time distribution*/  
  static const AliHLTComponentDataType fgkCellChannelDataDataType;  /**<Time dependent signal from the readout channels*/  
  static const AliHLTComponentDataType fgkClusterDataType;  //Cluster data type
  static const AliHLTComponentDataType fgkRecPointDataType; //RecPoint data type
  static const AliHLTComponentDataType fgkHistDataType;     //hist data type
  static const AliHLTComponentDataType fgkSpectrumDataType; //spectrum data type
  static const AliHLTComponentDataType fgkDigitDataType; //Digit data type
  static const AliHLTComponentDataType fgkRootTreeDataType; //Root tree type
  static const AliHLTComponentDataType fgkBaselineDataType; //Baseline type
  static const AliHLTComponentDataType fgkMIPDataType; //"MIP" data type
  static const AliHLTComponentDataType fgkNoiseMapDataType; //Noise map data type
  static const AliHLTComponentDataType fgkSandboxDataType; //General data type
  static const AliHLTComponentDataType fgkEmcCalibDataType; //Calibration data type
  static const AliHLTComponentDataType fgkCaloClusterDataType; //Calo cluster data type
};

#endif
