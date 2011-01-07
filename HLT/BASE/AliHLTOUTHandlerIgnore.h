//-*- Mode: C++ -*-
// $Id: $

#ifndef ALIHLTOUTHANDLERIGNORE_H
#define ALIHLTOUTHANDLERIGNORE_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTHandlerIgnore.h
    @author Artur Szostak <artursz@iafrica.com>
    @date   7 Jan 2011
    @brief  HLT output handler for ignoring data blocks.

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "AliHLTOUTHandler.h"

/**
 * @class AliHLTOUTHandlerIgnore
 * HLT output handler used to ignore data block types completely.
 * It will not inspect the data at all.
 */
class AliHLTOUTHandlerIgnore : public AliHLTOUTHandler
{
public:
  
  /// Default constructor
  AliHLTOUTHandlerIgnore() : AliHLTOUTHandler() {}
  
  /// Default destructor
  virtual ~AliHLTOUTHandlerIgnore() {}

  /**
   * Process the data will simply ignore the input data.
   * @param pData  instance of the AliHLTOUT data
   * @return always returns zero.
   */
  virtual int ProcessData(AliHLTOUT* data);

private:

  /// Do not allow copying of this class
  AliHLTOUTHandlerIgnore(const AliHLTOUTHandlerIgnore&);
  /// Do not allow copying of this class
  AliHLTOUTHandlerIgnore& operator = (const AliHLTOUTHandlerIgnore&);

  ClassDef(AliHLTOUTHandlerIgnore, 0)
};

#endif // ALIHLTOUTHANDLERIGNORE_H
