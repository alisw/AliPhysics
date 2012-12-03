//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRDMONITORCOMPONENT_H
#define ALIHLTTRDMONITORCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTTRDMonitorComponent.h
/// @author Felix Rettig, Stefan Kirsch
/// @date   2012-08-16
/// @brief  The TRD monitoring component
/// @ingroup alihlt_trd_components

#include "AliHLTProcessor.h"

class TObjArray;
class TH1I;
class TH2I;
class TH2F;
class AliTRDonlineTrackingDataContainer;

class AliHLTTRDMonitorComponent : public AliHLTProcessor {
public:
  AliHLTTRDMonitorComponent();
  virtual ~AliHLTTRDMonitorComponent();

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
  int ScanConfigurationArgument(int argc, const char** argv);
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTTRDMonitorComponent(const AliHLTTRDMonitorComponent&);
  /** assignment operator prohibited */
  AliHLTTRDMonitorComponent& operator=(const AliHLTTRDMonitorComponent&);

  int Configure(const char* arguments);

  void DbgLog(const char* prefix, ... );

  int PrepareTRDData();
  void DumpTrackingData();
  int ProcessTRDData();

  static const AliHLTEventID_t fgkInvalidEventId = 0xffffffffffffffffllu;
  static const unsigned int fkTRDChambers = 540;             //! number of chambers in TRD
  static const unsigned int fkTRDStacks = 90;                //! number of stacks in TRD
  static const unsigned int fkTRDStacksPerSector = 5;        //! number of stacks per sector in TRD
  static const unsigned int fkTRDSectors = 18;               //! number of sectors in TRD
  static const unsigned int fkMaxRefTracksPerStack = 1000;   //! maximum number of ref tracks per stack

  Double_t fTrackHighPtThreshold;           //! high-pt track pt threshold
  Bool_t fHistoMode;                        //! histogramming mode, 0: single event, 1: accumulative (debugging)
  Bool_t fTrackingDataDebugOutput;          //! switch on/off tracking data text dump

  UShort_t fDebugLevel;                                      //! debug level, 0: debug off
  Bool_t fWriteHistos;                                       //! switch on/off histogram writing

  AliHLTEventID_t fEventId;                                  //! hlt internal event id
  AliTRDonlineTrackingDataContainer* fTrackingData;          //! container for TRD tracking data

  TObjArray* fHistArray;
  TH1I* fHistEventTypes;                                     //! counting of event types
  TH1I* fHistTrackletY;                                      //! tracklet y-positions from all stacks
  TH1I* fHistTrackletDy;                                     //! tracklet deflections from all stacks
  TH2I* fHistTrackletYDy;                                    //! tracklet deflections vs. y-positions
  TH1I* fHistTrackletZ;                                      //! tracklet z-positions from all stacks
  TH1I* fHistTrackletPID;                                    //! tracklet PID values from all stacks
  TH2F* fHistTrackletsHCId;                                  //! number of tracklets per half-chamber
  TH1I* fHistTrackPt;                                        //! transverse momentum of GTU tracks from all stacks
  TH1I* fHistTrackPID;                                       //! PID of GTU tracks from all stacks
  TH1I* fHistTrackLayers;                                    //! contributing layers per GTU track
  TH1I* fHistTrackLayersHighPt;                              //! contributing layer per high-pt GTU track
  TH2F* fHistTracksStack;                                    //! GTU tracks per stack
  TH2I* fHistTrackletTimingStack;                            //! tracklet arrival timing by stack
  TH2I* fHistTrackingTiming;                                 //! tracking timing
  TH2I* fHistTriggerContribs;                                //! trigger contributions by sector

  ClassDef(AliHLTTRDMonitorComponent, 0)
};

#endif
