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
    };

#endif
