//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVENTSTATISTICS_H
#define ALIHLTEVENTSTATISTICS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEventStatistics.h
    @author Jochen Thaeder
    @date   
    @brief  Base class for event statistics, for all detectors
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

/**
 * @defgroup alihlt_run_statistics Event and run statistics for the HLT
 * This section describes the event and run statistics as well as the 
 * event and run summary handling for the HLT chain.
 */

/**
 * @defgroup alihlt_trigger Trigger components for the HLT.
 * This section describes the handling of different triggers of the HLT.
 * @ingroup alihlt_modules
 */

#include "TObject.h"
#include "TString.h"

#include "AliHLTDataTypes.h"

/**
 * @class  AliHLTEventStatistics
 * @brief  Base class for event statistics, for all detectors
 *
 * The event statistic classes hold information about certain characteristica 
 * of the processed events. They are devided into 3 parts. A base class 
 * @see AliHLTEventStatistics for general Information, detector specific
 * classes like @see AliHLTTPCEventStatistics for the TPC and a summary class
 * @see AliHLTEventStatisticsSummary which can hold several detector classes.
 *
 * This is the base class.
 *
 * Currently implemented detecor classes <br>
 * * @see AliHLTTPCEventStatistics <br>
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTEventStatistics : public TObject {
  
public:
  
  /** constructor */
  AliHLTEventStatistics();
  /** destructor */
  virtual ~AliHLTEventStatistics();

  /** Get detector name
   *  @return name of detector
   */
  TString GetDetectorName()                              { return fDetectorName; }

  /** Set Total number of tracks 
   *  @param s  number of tracks
   */
  void SetDetectorName( TString s )                      { fDetectorName = s; }

private:

  /** copy constructor prohibited */
  AliHLTEventStatistics (const AliHLTEventStatistics&);

  /** assignment operator prohibited */
  AliHLTEventStatistics& operator= (const AliHLTEventStatistics&);

  /** Detector Name */
  TString fDetectorName;                       // see above

  ClassDef(AliHLTEventStatistics, 0);

};
#endif

