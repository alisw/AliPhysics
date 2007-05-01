// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTPHOSDEFINITIONS_H
#define ALIHLTPHOSDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTPHOSDefinitions
 */

#include "AliHLTDataTypes.h"
#include "Rtypes.h"

class AliHLTPHOSDefinitions
    {
    public:
      static const AliHLTComponentDataType gkCellEnergyDataType;    /**<Reconstructed cell/crystal energies*/
      static const AliHLTComponentDataType gkDDLPackedRawDataType;  /**<DDL raw data on the RCU data format*/
      static const AliHLTComponentDataType gkCellEnergyHistogramDataType;  /**<Histogram for per cell/gain energy distribution*/
      static const AliHLTComponentDataType gkCellAverageEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
      static const AliHLTComponentDataType gkCellAccumulatedEnergyDataType;  /**<Histogram for per cell/gain energy distribution*/
      static const AliHLTComponentDataType gkCellTimingHistogramDataType;  /**<Histogram for per cell/gain time distribution*/      
      static const AliHLTComponentDataType gkCellTimingAverageDataType;  /**<Histogram for per cell/gain time distribution*/  
      static const AliHLTComponentDataType gkCellChannelDataDataType;  /**<Time dependent signal from the readout channels*/  
    };

#endif
