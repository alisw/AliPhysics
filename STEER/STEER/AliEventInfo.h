#ifndef ALIEVENTINFO_H
#define ALIEVENTINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                          Class AliEventInfo                              //
//   Container class for all the information related to                     //
//   event type, trigger mask and trigger clusters.                         //
//   It is used together with AliRunInfo in order to provide the detector's //
//   AliRecoParam object with                                               //
//   the necessary information so that it can decide which instance of      //
//   AliDetectorRecoParam objects to use in reconstruction one particular   //
//   event.                                                                 //
//                                                                          //
//   cvetan.cheshkov@cern.ch 12/06/2008                                     //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjString.h>

class AliEventInfo : public TObject {

 public:
  AliEventInfo();
  AliEventInfo(UInt_t evType,
	       const char *classes,
	       ULong64_t mask,
	       const char *cluster,
	       const char *decision);
  virtual ~AliEventInfo() {}

  void SetEventType(UInt_t evType) { fEventType = evType; }
  void SetTriggerClasses(const char *classes) { fTriggerClasses = classes; }
  void SetTriggerMask(ULong64_t mask) { fTriggerMask = mask; }
  void SetTriggerCluster(const char *cluster) { fTriggerCluster = cluster; }
  void SetHLTDecision(const char *decision) { fHLTDecision = decision; }

  virtual void Print(Option_t */*option=""*/) const { Dump(); }

  UInt_t      GetEventType() const { return fEventType; }
  const char *GetTriggerClasses() const { return fTriggerClasses.Data(); }
  ULong64_t   GetTriggerMask() const { return fTriggerMask; }
  const char *GetTriggerCluster() const { return fTriggerCluster.Data(); }
  const char *GetHLTDecision() const { return fHLTDecision.Data(); }

  AliEventInfo(const AliEventInfo &evInfo);
  AliEventInfo& operator= (const AliEventInfo& evInfo);

  void Reset();

 private:

  UInt_t      fEventType;      // event type as defined by DAQ (start_of_*,calibration,physics etc) (per event)
  TString     fTriggerClasses; // list of fired trigger classes (per event)
  ULong64_t   fTriggerMask;    // trigger mask as received from DAQ or CTP raw-data payload (per event)
  TString     fTriggerCluster; // list of detectors that have been read out (per event)
  TString     fHLTDecision;    // string describing the HLT decision (per event)

  ClassDef(AliEventInfo,3)     // Event info class
};

#endif
