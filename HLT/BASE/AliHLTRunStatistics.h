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

#include "TNamed.h"
#include "TObjArray.h"
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
 * Currently implemented detecor classes <br>
 * * @see AliHLTTPCRunStatistics <br>
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTRunStatistics : public TNamed {
  
public:
  
  /** constructor */
  AliHLTRunStatistics();
  /** copy constructor */
  AliHLTRunStatistics (const AliHLTRunStatistics&);
  /** assignment operator */
  AliHLTRunStatistics& operator= (const AliHLTRunStatistics&);
  /** destructor */
  virtual ~AliHLTRunStatistics();

  /// Get detector name
  TString GetDetectorName() const        { return GetName(); }

  /// Set detector name
  void SetDetectorName( TString s )      { SetName(s); }
  /// Set detector name
  void SetDetectorName(const char* name) { SetName(name); }

  /// inherited from TObject
  virtual void Print(Option_t* option) const;

  /// inherited from TObject, copy to the target object
  virtual void Copy(TObject &object) const;

  /// Inherited from TObject. Create a new clone.
  virtual TObject *Clone(const char *newname="") const;

  /// Inherited from TObject. Clear content.
  virtual void  Clear(Option_t * option ="");

  // -- event parameters ------------------------

  /// Set Number of events 
  void SetNEvents( ULong_t i )           { fNEvents = i; }

  /// Increment event count
  void IncrementNEvents()                { fNEvents++; }

  /// Get number of events
  ULong_t GetNEvents() const             { return fNEvents; }

  /// Add clone of object
  int Add(const TObject* pObject);

  /// Find an object
  virtual TObject *FindObject(const char *name) const {
    return fMyObjects.FindObject(name); }

private:

  /** Number of events */
  ULong_t fNEvents;                              // see above

  /// array of statistics objects owned by the array
  TObjArray fMyObjects;                          // see above

  ClassDef(AliHLTRunStatistics, 1);

};
#endif
