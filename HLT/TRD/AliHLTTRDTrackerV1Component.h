//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRDTRACKERV1COMPONENT_H
#define ALIHLTTRDTRACKERV1COMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDTrackerV1Component.h
    @author Theodor Rascanu
    @date   
    @brief  Declaration of a TRDTracker component. */

#include "AliHLTProcessor.h"
#include "AliESDEvent.h"

class TFile;
class TGeoManager;
class AliCDBManager;
class AliMagF;
class AliTRDtrackerV1;
class AliTRDrecoParam;
class AliTRDReconstructor;
class TClonesArray;
class TBuffer;

/**
 * @class AliHLTTRDTrackerV1Component
 * @brief A TRDTrackerV1 HLT processing component. 
 *
 * Uses the second generation TRD tracker AliTRDtrackerV1
 */

class AliHLTTRDTrackerV1Component : public AliHLTProcessor
{
protected:

  class AliHLTTRDESDEvent : public AliESDEvent
  {
    // Here we use our own slim version of AliESDEvent 
    // by overwriting AliESDEvent::CreateStdContent
  public:
    void CreateStdContent();

    // streamer function is prohibited, class for internal use only
    void Streamer(TBuffer &R__b);
  };

public:
  AliHLTTRDTrackerV1Component();
  virtual ~AliHLTTRDTrackerV1Component();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();
  void GetOCDBObjectDescription(TMap* const targetMap);

protected:
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, 
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponent_TriggerData& /*trigData*/, 
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, 
	       vector<AliHLTComponent_BlockData>& outputBlocks );
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);

  using AliHLTProcessor::DoEvent;
  
  int Configure(const char* arguments);
  int SetParams();

protected:
  /** copy constructor prohibited */
  AliHLTTRDTrackerV1Component(const AliHLTTRDTrackerV1Component&);
  /** assignment operator prohibited */
  AliHLTTRDTrackerV1Component& operator=(const AliHLTTRDTrackerV1Component&);

  // The size of the output data produced, as a percentage of the input data's size.
  // Can be greater than 100 (%)
  unsigned fOutputPercentage; // Output volume in percentage of the input
	
  AliTRDtrackerV1 *fTracker;//! Offline-pure/HLT tracker V1
  AliTRDrecoParam *fRecoParam; //! Offline reco params
  AliTRDReconstructor* fReconstructor;
  AliHLTTRDESDEvent* fESD;

  TClonesArray* fClusterArray;

  Int_t fRecoParamType;       // default will be the low flux
  Int_t fNtimeBins;           // number of time bins for the tracker to use
  Int_t fPIDmethod;           // 0=LikelyHood(LH) 1=NeuronalNetwork(NN) 2=TruncatedMean(TM)
  TString fgeometryFileName;
  Bool_t fHLTflag;            // use HLT flag in reconstructor
  Bool_t fOutputV1Tracks;     // output TRD tracks, needed for calibration and ESDfriends
  Bool_t fHighLevelOutput;    // do we what to have high level output (only for debuging)
  Bool_t fEmulateHLTTracks;   // for debugging data containers
  Bool_t fImproveTracklets;   // improve tracklets?

  ClassDef(AliHLTTRDTrackerV1Component, 4)

};
#endif
