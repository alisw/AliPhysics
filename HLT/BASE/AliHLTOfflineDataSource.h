// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOFFLINEDATASOURCE_H
#define ALIHLTOFFLINEDATASOURCE_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOfflineDataSource.h
    @author Matthias Richter
    @date   
    @brief  AliRoot data sink component base class.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTDataSource.h"
#include "AliHLTOfflineInterface.h"

/******************************************************************************/

/**
 * @class AliHLTOfflineDataSource
 * The class implements a AliRoot data source component base class. 
 * The child class must implement the functions:
 * - @ref DoInit (optional)
 * - @ref DoDeinit (optional)
 * - @ref GetEvent
 * - @ref GetComponentID
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTOfflineDataSource 
: public AliHLTDataSource, public AliHLTOfflineInterface {
 public:
  /** standard constructor */
  AliHLTOfflineDataSource();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTOfflineDataSource(const AliHLTOfflineDataSource&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTOfflineDataSource& operator=(const AliHLTOfflineDataSource&);
  /** destructor */
  virtual ~AliHLTOfflineDataSource();

  /**
   * Default implementation as sources do not have a real FillESD method.
   */
  int FillESD(int /*eventNo*/, AliRunLoader* /*runLoader*/, AliESDEvent* /*esd*/) {
    return 0;
  }

 private:
 
  ClassDef(AliHLTOfflineDataSource, 1);
};

#endif
