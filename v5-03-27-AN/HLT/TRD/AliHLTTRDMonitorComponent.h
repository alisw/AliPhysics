//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRDMONITORCOMPONENT_H
#define ALIHLTTRDMONITORCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTTRDMonitorComponent.h
/// @author Felix Rettig, Stefan Kirsch
/// @date   2011-08-02
/// @brief  A FEP-level processing component for TRD tracking/trigger data
/// @ingroup alihlt_trd_components

#include "AliHLTProcessor.h"

class AliRawReaderMemory;
class TTree;
class AliTRDdigitsManager;
class AliTRDrawStream;
class TH1F;
class TH1I;
class TH2I;
class TH2F;

/**
 * @class AliHLTTRDMonitorComponent
 * Component fetches raw data input objects in DDL format and extracts tracklets and GTU tracks.
 *  It also instantiates a RawReader in order to be used with some reconstruction.
 *
 * More information and examples can be found here (relative to $ALICE_ROOT):
 * 
 * -- HLT/BASE/AliHLTComponent.h/.cxx,  HLT/BASE/AliHLTProcessor.h/.cxx
 *    Interface definition and description
 * -- HLT/SampleLib: example implementations of components
 *
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b AliHLTTRDMonitorComponent <br>
 * Library: \b libAliHLTTRD.so     <br>
 * Input Data Types: @ref kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTRD <br>
 * Output Data Types: @ref kAliHLTTrackDataTypeID|kAliHLTDataOriginTRD <br>
 *
 * <h2>Mandatory arguments:</h2>
 * none     
 *
 * <h2>Optional arguments:</h2>
 * none
 *
 * <h2>Configuration:</h2>
 * none
 *
 * <h2>Default CDB entries:</h2>
 * none
 *
 * <h2>Performance:</h2>
 * minmal
 *
 * <h2>Memory consumption:</h2>
 * don't know yet
 *
 * <h2>Output size:</h2>
 * not very much
 *
 * @ingroup The component has no output data.
 */
class AliHLTTRDMonitorComponent : public AliHLTProcessor {
public:
  AliHLTTRDMonitorComponent();
  virtual ~AliHLTTRDMonitorComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  void GetOCDBObjectDescription( TMap* const targetMap);

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int ScanConfigurationArgument(int argc, const char** argv);
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTTRDMonitorComponent(const AliHLTTRDMonitorComponent&);
  /** assignment operator prohibited */
  AliHLTTRDMonitorComponent& operator=(const AliHLTTRDMonitorComponent&);

  // trd specific data
  TClonesArray* fTrackletArray;
  TClonesArray *fGtuTrackArray;

  // rawreader instance
  AliRawReaderMemory* fRawReaderMem;
  AliTRDdigitsManager *fDigitsManagerTrd;
  AliTRDrawStream*    fRawReaderTrd; 

  // FEE statistics data
  TObjArray* fHistArray;
  TH1I *fHistTrackletY;
  TH1I *fHistTrackletDy;
  TH1I *fHistTrackletZ;
  TH1I *fHistTrackletPID;
  TH2I *fHistTrackletYDy;
  TH2I *fHistTrackletHC;
  TH2I *fHistTrackletBadY;
  TH2I *fHistTrackletBadPID;
  TH1F *fHistFirstTrackletTime;
  TH1F *fHistLastTrackletTime;
  // GTU statistics data
  TH1F *fHistTmuTime;
  TH1F *fHistSmuTime;
  TH1F *fHistTrackPt;
  TH2I *fHistTrackStack;
  TH1I *fHistTrackletsTrack;
  TH1I *fHistTrackletsTrackHpt;
  TH1I *fHistTrackPID;
  TH2I *fHistTriggerContribs;

  ClassDef(AliHLTTRDMonitorComponent, 0)
};
#endif
