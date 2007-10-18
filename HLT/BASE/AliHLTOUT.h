//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTOUT_H
#define ALIHLTOUT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUT.h
    @author Matthias Richter
    @date   
    @brief  The control class for HLTOUT data.

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "AliHLTLogging.h"

/**
 * @class AliHLTOUT
 * The control class for HLTOUT data.
 * The output of the HLT, either from the HLTOUT nodes or simulated output,
 * is transferred and stored in the HOMER format. The AliHLTOUT class 
 * implements scanning of the HOMER data for all HLTOUT DDL links and
 * abstracts access to the complete HLTOUT data.
 * 
 */
class AliHLTOUT : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTOUT();
  /** standard destructor */
  virtual ~AliHLTOUT();

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTOUT(const AliHLTOUT&);
  /** assignment operator prohibited */
  AliHLTOUT& operator=(const AliHLTOUT&);

  ClassDef(AliHLTOUT, 0)
};
#endif
