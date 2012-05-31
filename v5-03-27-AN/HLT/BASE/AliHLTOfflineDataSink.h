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
 * @note This class is very likely to be deprecated. According to the new
 * reconstruction scheme, the esd is no longer filled by components in the
 * reconstruction chain, but added to the HLTOUT data. The HLTOUT is
 * processed during AliReconstruction at the end of the HLT event processing,
 * literally during the FillESD method of the AliRoot reconstruction
 * interface. The HLT module must implement HLTOUT handlers and provide
 * those through the module agent.
 *
 * @ingroup alihlt_system
 */
class AliHLTOfflineDataSink 
: public AliHLTDataSink, public AliHLTOfflineInterface 
{
 public:
  /** standard constructor */
  AliHLTOfflineDataSink();
  /** destructor */
  virtual ~AliHLTOfflineDataSink();

 private:
  /** copy constructor prohibited */
  AliHLTOfflineDataSink(const AliHLTOfflineDataSink&);
  /** assignment operator prohibited */
  AliHLTOfflineDataSink& operator=(const AliHLTOfflineDataSink&);

  ClassDef(AliHLTOfflineDataSink, 1);
};

#endif
