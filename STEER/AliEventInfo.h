#ifndef ALIEVENTINFO_H
#define ALIEVENTINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                          Class AliEventInfo                              //
//   Container class for all the information related to LHCstate, run and   //
//   event types, trigger mask and trigger clusters.                        //
//   It is used in order to provide the detector's AliRecoParam objects with//
//   the necessary information so that they can decide which instance of    //
//   AliDetectorRecoParam to use in reconstruction one particular event.    //
//                                                                          //
//   cvetan.cheshkov@cern.ch 12/06/2008                                     //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjString.h>

class AliEventInfo : public TObject {

 public:
  AliEventInfo();
  AliEventInfo(const char *lhcState,
	       const char *beamType,
	       const char *runType,
	       const char *activeDetectors);
  virtual ~AliEventInfo() {}

  void SetEventType(UInt_t evType) { fEventType = evType; }
  void SetTriggerClasses(const char *classes) { fTriggerClasses.SetString(classes); }
  void SetTriggerMask(ULong64_t mask) { fTriggerMask = mask; }
  void SetTriggerCluster(const char *cluster) { fTriggerCluster.SetString(cluster); }
  void SetHLTDecision(const char *decision) { fHLTDecision.SetString(decision); }

  virtual void Print(Option_t */*option=""*/) const { Dump(); }

  const char *GetLHCState() const { return fLHCState.GetString().Data(); }
  const char *GetBeamType() const { return fBeamType.GetString().Data(); }
  const char *GetRunType() const { return fRunType.GetString().Data(); }
  const char *GetActiveDetectors() const { return fActiveDetectors.GetString().Data(); }
  UInt_t      GetEventType() const { return fEventType; }
  const char *GetTriggerClasses() const { return fTriggerClasses.GetString().Data(); }
  ULong64_t   GetTriggerMask() const { return fTriggerMask; }
  const char *GetTriggerCluster() const { return fTriggerCluster.GetString().Data(); }
  const char *GetHLTDecision() const { return fHLTDecision.GetString().Data(); }

  AliEventInfo(const AliEventInfo &evInfo);
  AliEventInfo& operator= (const AliEventInfo& evInfo);

  void Reset();
 private:

  TObjString  fLHCState;       // state of the machine as provided by DCS and DAQ log-book (per run)
  TObjString  fBeamType;       // beam type for Alice
  TObjString  fRunType;        // run type accoring to ECS (per run)
  TObjString  fActiveDetectors;// list of active detectors (per run)
  UInt_t      fEventType;      // event type as defined by DAQ (start_of_*,calibration,physics etc) (per event)
  TObjString  fTriggerClasses; // list of fired trigger classes (per event)
  ULong64_t   fTriggerMask;    // trigger mask as received from DAQ or CTP raw-data payload (per event)
  TObjString  fTriggerCluster; // list of detectors that have been read out (per event)
  TObjString  fHLTDecision;    // string describing the HLT decision (per event)

  ClassDef(AliEventInfo,2)     // Event info class
};

#endif
