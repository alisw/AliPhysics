//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTRUNSTATISTICS_H
#define ALIHLTRUNSTATISTICS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRunStatistics.h
    @author Jochen Thaeder
    @date   
    @brief  Base class for run statistics, for all detectors
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "TObject.h"
#include "TString.h"

#include "AliHLTDataTypes.h"

/**
 * @class  AliHLTRunStatistics
 * @brief  Base class for run statistics, for all detectors
 *
 * The run statistic classes hold information / histograms about certain 
 * characteristica of the processed events. They are devided into 3 parts. 
 * A base class  @see AliHLTRunStatistics for general Information, detector 
 * specific classes like @see AliHLTTPCRunStatistics for the TPC and a summary 
 * class @see AliHLTRunStatisticsSummary which can hold several detector classes.
 *
 * This is the base class.
 *
 * Currently implemented detecor classes<BR/>
 * * @see AliHLTTPCRunStatistics<BR/>
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTRunStatistics : public TObject {
  
public:
  
  /** constructor */
  AliHLTRunStatistics();
  /** destructor */
  virtual ~AliHLTRunStatistics();

  /** Get detector name
   *  @return name of detector
   */
  TString GetDetectorName()              { return fDetectorName; }

  /** Set Total number of tracks 
   *  @param s name of detector
   */
  void SetDetectorName( TString s )      { fDetectorName = s; }

  // -- event parameters ------------------------

  /** Set Number of events 
   *  @param i number of events
   */
  void SetNEvents( ULong_t i )           { fNEvents = i; }

  /** Add events */
  void AddNEvents()                      { fNEvents++; }

  /** Get number of events
   *  @return number of events
   */
  ULong_t GetNEvents()                   { return fNEvents; }

private:

  /** copy constructor prohibited */
  AliHLTRunStatistics (const AliHLTRunStatistics&);

  /** assignment operator prohibited */
  AliHLTRunStatistics& operator= (const AliHLTRunStatistics&);

  /** Detector Name */
  TString fDetectorName;                         // see above

  /** Number of events */
  ULong_t fNEvents;                              // see above

  ClassDef(AliHLTRunStatistics, 0);

};
#endif

