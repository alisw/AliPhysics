//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTOUTHANDLER_H
#define ALIHLTOUTHANDLER_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTHandler.h
    @author Matthias Richter
    @date   
    @brief  Base class declaration of HLTOUT handlers

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "AliHLTLogging.h"

class AliHLTOUT;

/**
 * @class AliHLTOUTHandler
 * Base class declaration of HLT output handlers.
 * The library implementation of the AliHLTModuleAgent allows to generate
 * handlers for data blocks of the HLT output. This can be the output of
 * the real HLT coming from the HLTOUT nodes, or simulated HLT output.   <br>
 * \em Note: The created instance of AliHLTOUTHandler is
 * deleted by the framework.
 */
class AliHLTOUTHandler : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTOUTHandler();
  /** standard destructor */
  virtual ~AliHLTOUTHandler();

 protected:

  /**
   * Process the data.
   */
  virtual int ProcessData(AliHLTOUT* pData) = 0;

 private:
  /** copy constructor prohibited */
  AliHLTOUTHandler(const AliHLTOUTHandler&);
  /** assignment operator prohibited */
  AliHLTOUTHandler& operator=(const AliHLTOUTHandler&);

  ClassDef(AliHLTOUTHandler, 0)
};
#endif
