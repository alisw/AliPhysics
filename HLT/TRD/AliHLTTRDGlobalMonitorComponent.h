//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRDGLOBALMONITORCOMPONENT_H
#define ALIHLTTRDGLOBALMONITORCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTTRDGlobalMonitorComponent.h
/// @author Felix Rettig, Stefan Kirsch
/// @date   2011-08-02
/// @brief  A processing component for TRD tracking/trigger data on CN-level
/// @ingroup alihlt_trd_components

#include "AliHLTProcessor.h"

class TObjArray;
class TH1I;
class TH2I;
class TH1F;

class AliHLTTRDGlobalMonitorComponent : public AliHLTProcessor {
public:
  AliHLTTRDGlobalMonitorComponent();
  virtual ~AliHLTTRDGlobalMonitorComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTTRDGlobalMonitorComponent(const AliHLTTRDGlobalMonitorComponent&);
  /** assignment operator prohibited */
  AliHLTTRDGlobalMonitorComponent& operator=(const AliHLTTRDGlobalMonitorComponent&);

  int Configure(const char* arguments);

  /** data **/
  TObjArray* fHistArray;
  TH1I *fHistTrackletY;   // tracklet y-positions from all stacks
  TH1I *fHistTrackletDy;  // tracklet deflections from all stacks
  TH1I *fHistTrackletZ;   // tracklet z-positions from all stacks
  TH1I *fHistTrackletPID; // tracklet PID values from all stacks
  TH2I *fHistTrackletYDy; // tracklet deflection vs. position from all stacks
  TH2I *fHistTrackletHC;  // tracklet numbers by half-chamber from all stacks
  TH2I *fHistTrackletBadY; // tracklet numbers with invalid y-position by stack from all stacks
  TH2I *fHistTrackletBadPID; // tracklet numbers with invalid PID value by stack from all stacks
  TH1F *fHistFirstTrackletTime; // arrival time of the first tracklet of each link
  TH1F *fHistLastTrackletTime; // arrival time of the last tracklet of each link
  TH1F *fHistTmuTime;     // tracking done time for each TMU
  TH1F *fHistSmuTime;     // tracking done time for each SMU
  TH1F *fHistTrackPt;     // transverse momentum of GTU tracks from all stacks
  TH1I *fHistTrackPID;    // PID of GTU tracks from all stacks
  TH2I *fHistTrackStack;  // GTU track numbers by stack from all stacks
  TH1I *fHistTrackletsTrack; // tracklets per GTU track from all stacks
  TH1I *fHistTrackletsTrackHpt; // tracklets per high-pt GTU track from all stacks
  TH2I *fHistTriggerContribs; // sector-level trigger contributions from all stacks
  ClassDef(AliHLTTRDGlobalMonitorComponent, 0)
};
#endif
