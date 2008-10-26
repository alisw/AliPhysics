//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHANDLERDETECTORDDL_H
#define ALIHLTOUTHANDLERDETECTORDDL_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTHandlerDetectorDDL.h
    @author Matthias Richter
    @date   2008-09-09
    @brief  Default HLTOUT handler returning equipment id from data type and
            bit pattern in spec.
*/

#include "AliHLTOUTHandlerEquId.h"

/**
 * @class AliHLTOUTHandlerDetectorDDL
 * A default handler class for DDL raw data redirection handlers.
 *
 * This class implements an AliHLTOUTHandlerEquId which extracts the
 * equipment Id from the bit pattern in the specification. All detectors
 * with up to 32 DDL links follow this convention. The bit no in the
 * data specification word corresponds to the DDL number within the
 * sub-detector.
 *
 * DDL offsets for sub-detectors can be fetched by means of AliDAQ.
 * The class must be initialized with the detector identification and the
 * data type it should be used for. 
 * @note The detector identification is according to AliDAQ. E.g. for ITS
 * and MUON there are ITSSPD, ITSSDD, ITSSSD, and MUONTRK and MUONTRG
 * respectively.
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTOUTHandlerDetectorDDL : public AliHLTOUTHandlerEquId {
 public:
  /** constructor 
   * the class is initialized with the detector identification and the
   * data type it should be used for. Note: the detector identification
   * is according to AliDAQ. E.g. for ITS and MUON there are ITSSPD,
   * ITSSDD, ITSSSD, and MUONTRK and MUONTRG respectively.
   */
  AliHLTOUTHandlerDetectorDDL(const char* detector, AliHLTComponentDataType dt);
  /** standard destructor */
  virtual ~AliHLTOUTHandlerDetectorDDL();

  /**
   * Process a data block.
   * Derives the eqipment ID from the DDL offset of the detector and
   * the DDL no within the detector which corresponds to a bit in the
   * data specification. Only one bit is allowed to be set.
   * @return equipment id the block should be used for.
   */
  virtual int ProcessData(AliHLTOUT* pData);

 private:
  /** standard constructor prohibited */
  AliHLTOUTHandlerDetectorDDL();
  /** copy constructor prohibited */
  AliHLTOUTHandlerDetectorDDL(const AliHLTOUTHandlerDetectorDDL&);
  /** assignment operator prohibited */
  AliHLTOUTHandlerDetectorDDL& operator=(const AliHLTOUTHandlerDetectorDDL&);

  int fDDLOffset; //!transient
  int fNumberOfDDLs; //!transient
  AliHLTComponentDataType fDt; //!transient

  ClassDef(AliHLTOUTHandlerDetectorDDL, 0)
};
#endif
