//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHANDLEREQUID_H
#define ALIHLTOUTHANDLEREQUID_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTHandlerEquId.h
    @author Matthias Richter
    @date   
    @brief  HLTOUT handler returning equipment id from data type and spec.
*/

#include "AliHLTOUTHandler.h"

/**
 * @class AliHLTOUTHandlerEquId
 * Base class for DDL raw data redirection handlers.
 *
 * HLTOUT handlers of this type are used for the replacement of detector
 * reconstruction input by the corresponding data from the HLT output.
 * From the data type and specification of an HLT output block the
 * corresponding equipment id of the original detector streams is determined.
 *
 * The class introduces a layer in the HLTOUT handler hierarchy in order to
 * collect all handlers which can be used by the AliRawReaderHLT for
 * redirection of HLTOUT data blocks and replacement of original detector
 * data. The common denominator of all those handlers is the data format.
 * The data block itself or the decoded data produced by the handler have
 * exactly the detector raw data format. Thus, the data streams can be
 * transparently replaced in the AliRoot reconstruction.
 *
 * The handler might produce decoded data from the data block in order
 * the get the right data format. The AliRawReaderHLT will try to fetch
 * those data by calling AliHLTOUTHandler::GetProcessedData(). If no data
 * is provided, the input block itself is treated as the data to redirect.
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTOUTHandlerEquId : public AliHLTOUTHandler {
 public:
  /** standard constructor */
  AliHLTOUTHandlerEquId();
  /** standard destructor */
  virtual ~AliHLTOUTHandlerEquId();

  /**
   * Process a data block.
   * The handler retrieves the data and it's properties and derives the
   * equipment id from it. The default behavior returns the specification as
   * equipment id and does not touch the data itself.
   * @return equipment id the block should be used for.
   */
  virtual int ProcessData(AliHLTOUT* pData);

 private:
  /** copy constructor prohibited */
  AliHLTOUTHandlerEquId(const AliHLTOUTHandlerEquId&);
  /** assignment operator prohibited */
  AliHLTOUTHandlerEquId& operator=(const AliHLTOUTHandlerEquId&);

  ClassDef(AliHLTOUTHandlerEquId, 0)
};
#endif
