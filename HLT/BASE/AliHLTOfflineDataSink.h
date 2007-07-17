// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOFFLINEDATASINK_H
#define ALIHLTOFFLINEDATASINK_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOfflineDataSink.h
    @author Matthias Richter
    @date   
    @brief  AliRoot data sink component base class.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTDataSink.h"
#include "AliHLTOfflineInterface.h"

/******************************************************************************/

/**
 * @class AliHLTOfflineDataSink
 * The class implements a AliRoot data sink component base class. 
 * The child class must implement the functions:
 * - @ref DoInit (optional)
 * - @ref DoDeinit (optional)
 * - @ref DumpEvent
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref Spawn
 * - @ref FillESD
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTOfflineDataSink 
: public AliHLTDataSink, public AliHLTOfflineInterface 
{
 public:
  /** standard constructor */
  AliHLTOfflineDataSink();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTOfflineDataSink(const AliHLTOfflineDataSink&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTOfflineDataSink& operator=(const AliHLTOfflineDataSink&);
  /** destructor */
  virtual ~AliHLTOfflineDataSink();

 private:
  ClassDef(AliHLTOfflineDataSink, 1);
};

#endif
