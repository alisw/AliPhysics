// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOFFLINEDATASOURCE_H
#define ALIHLTOFFLINEDATASOURCE_H
///* This file is property of and copyright by the                          * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///* See cxx source for full Copyright notice                               *

/// @file   AliHLTOfflineDataSource.h
/// @author Matthias Richter
/// @date   
/// @brief  AliRoot data sink component base class.
///

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
  /** destructor */
  virtual ~AliHLTOfflineDataSource();

  /**
   * Default implementation as sources do not have a real FillESD method.
   */
  int FillESD(int /*eventNo*/, AliRunLoader* /*runLoader*/, AliESDEvent* /*esd*/) {
    return 0;
  }

 private:
  /** copy constructor prohibited */
  AliHLTOfflineDataSource(const AliHLTOfflineDataSource&);
  /** assignment operator prohibited */
  AliHLTOfflineDataSource& operator=(const AliHLTOfflineDataSource&);
 
  ClassDef(AliHLTOfflineDataSource, 0);
};

#endif
