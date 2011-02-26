//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTRUNSUMMARY_H
#define ALIHLTRUNSUMMARY_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRunSummary.h
    @author Jochen Thaeder
    @date  
    @brief  Summary class for a run, merges all detectors 
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "TObject.h"
#include "TObjArray.h"

#include "AliHLTDataTypes.h"

/**
 * @class  AliHLTRunSummary
 * @brief  Summary class for a run, merges all detectors 
 *
 * The run statistic classes hold information / histograms about certain 
 * characteristica of the processed events. They are devided into 3 parts. 
 * A base class  @see AliHLTRunStatistics for general Information, detector 
 * specific classes like @see AliHLTTPCRunStatistics for the TPC and a summary 
 * class @see AliHLTRunSummary which can hold several detector classes.
 * 
 * This is the summary class.
 * 
 * See base class @see AliHLTRunStatistics for further information.
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTRunSummary : public TObject {
  
public:
  
  /** constructor */
  AliHLTRunSummary();
  /** destructor */
  virtual ~AliHLTRunSummary();

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
  ULong_t GetNEvents() const             { return fNEvents; }
 
  /** Set number of rejected events 
   *  @param i number of events rejected
   */
  void SetNEventsRejected( ULong_t i )   { fNEventsRejected = i; }

  /** Add rejected events */
  void AddNEventsRejected()              { fNEventsRejected++; }

  /** Get number of reejcted events 
   *  @return number of events rejected
   */
  ULong_t GetNEventsRejected() const     { return fNEventsRejected; }

  // -- run parameters ------------------------

  /** Set Run Number 
   *  @param i run number
   */
  void SetRunNumber( AliHLTUInt32_t i )           { fRunNumber = i; }

  /** Get Run Number 
   *  @return run number
   */
  AliHLTUInt32_t GetRunNumber() const             { return  fRunNumber;}

  /** Set Run Type 
   *  @param i run type
   */
  void SetRunType( AliHLTUInt32_t i )             { fRunType = i; }

  /** Get Run Type s
   *  @return run type
   */
  AliHLTUInt32_t GetRunType() const               { return fRunType; }

  // -- trigger parameters ------------------------

  /** Set ocurrance of trigger classe with index ndx
   *  @param ndx index of Trigger Class   
   *  @return kTRUE on success
   */
  Bool_t AddTriggerClass( Int_t ndx );

  /** Get ocurrance of trigger classes 
   *  @return ptr to array of trigger classes
   */
  const ULong_t* GetTriggerClasses() const { return fTriggerClasses; }
  ULong_t* GetTriggerClasses()             { return fTriggerClasses; }
  
  // -- detector parameters ------------------------

  /** Detector run statistics classes 
   *  @return ptr to Detector arry
   */
  TObjArray* GetDetectorArray ()         { return fDetectorArray; }

  /** Rest the Detector array, all elements are removed */
  void ResetDetectorArray()              { if ( fDetectorArray ) fDetectorArray->Clear(); }

  /** Add another detector run statistics class to detector arry 
   *  @param obj to TObject, which should be added
   *  @return kTRUE on success
   */
  Bool_t AddDetector( TObject *obj );

private:

  /** copy constructor prohibited */
  AliHLTRunSummary (const AliHLTRunSummary&);

  /** assignment operator prohibited */
  AliHLTRunSummary& operator= (const AliHLTRunSummary&);

  /** Number of events */
  ULong_t fNEvents;                              // see above

  /** Number of rejected events */
  ULong_t fNEventsRejected;                      // see above

  /** Run Number */   
  AliHLTUInt32_t fRunNumber;                     // see above

  /** Run Type */
  AliHLTUInt32_t fRunType;                       // see above

  /** Ocurrance of trigger classes */
  ULong_t fTriggerClasses[gkNCTPTriggerClasses]; // see above

  /** Detector run statistics classes */
  TObjArray* fDetectorArray;                     //! transient

  ClassDef(AliHLTRunSummary, 1);  
};
#endif
