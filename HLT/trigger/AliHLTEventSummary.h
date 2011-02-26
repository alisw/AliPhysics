//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVENTSUMMARY_H
#define ALIHLTEVENTSUMMARY_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEventSummary.h
    @author Jochen Thaeder
    @date   
    @brief  Summary class for run statistics, merges all detectors
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
#include "AliHLTDataTypes.h"

/**
 * @class  AliHLTEventSummary
 * @brief  Summary  class for event statistics, merges all detectors
 *
 * The event statistic classes hold information about certain characteristica 
 * of the processed events. They are devided into 3 parts. A base class 
 * @see AliHLTEventStatistics for general Information, detector specific
 * classes like @see AliHLTTPCEventStatistics for the TPC and a summary class
 * @see AliHLTEventSummary which can hold several detector classes.
 *
 * This is the summary class.
 *
 * See base class @see AliHLTEventStatistics for further information.
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTEventSummary : public TObject {
  
public:
  
  /** constructor */
  AliHLTEventSummary();
  /** destructor */
  virtual ~AliHLTEventSummary();

  // -- event reject  ------------------------

  /** Reject this event */
  void RejectEvent()                              { fRejected = kFALSE;} 

  /** Check event is accepted
   *  @return kTRUE if accepted, kFALSE if rejected
   */
  Bool_t IsAccepted() const                       { return fRejected;} 

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

  /** Set ocurrance of trigger classe with index i
   *  @param u index of Trigger Class   
   *  @return kTRUE on success
   */
  void SetTriggerClass( AliHLTUInt64_t u )        { fTriggerClass = u; } 

  /** Get ocurrance of trigger classes 
   *  @return ptr to array of trigger classes
   */
  AliHLTUInt64_t GetTriggerClasses() const        { return fTriggerClass; }
  
  // -- detector parameters ------------------------

  /** Detector run statistics classes 
   *  @return ptr to Detector arry
   */
  const TObjArray* GetDetectorArray () const      { return fDetectorArray; }
  TObjArray* GetDetectorArray ()                  { return fDetectorArray; }

  /** Rest the Detector array, all elements are removed */
  void ResetDetectorArray()                       { if ( fDetectorArray ) fDetectorArray->Clear(); }

  /** Add another detector event statistics class to detector arry 
   *  @param obj to TObject, which should be added
   *  @return kTRUE on success
   */
  Bool_t AddDetector( TObject *obj );

private:

  /** copy constructor prohibited */
  AliHLTEventSummary (const AliHLTEventSummary&);

  /** assignment operator prohibited */
  AliHLTEventSummary& operator= (const AliHLTEventSummary&);

  /** Event rejected 
   *  - kFALSE is rejected
   *  - kTRUE is default
   */
  Bool_t fRejected;                              // see above

  /** Run Number */   
  AliHLTUInt32_t fRunNumber;                     // see above

  /** Run Type */
  AliHLTUInt32_t fRunType;                       // see above

  /** Trigger class */
  AliHLTUInt64_t fTriggerClass;                  // see above

  /** Detector run statistics classes */
  TObjArray* fDetectorArray;                     //! transient

  ClassDef(AliHLTEventSummary, 1);

};
#endif
