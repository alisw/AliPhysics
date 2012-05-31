//-*- Mode: C++ -*-
// $Id: AliHLTTRDCalibHistoComponent.h 40269 2010-04-08 22:08:53Z richterm $

#ifndef ALIHLTTRDCALIBHISTOCOMPONENT_H
#define ALIHLTTRDCALIBHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTRDCalibHistoComponent.h
//  @author 
//  @date   
//  @brief  Declaration of a TRDCalibration component. 
// 


#include "AliHLTProcessor.h"
class AliCDBManager;
class AliTRDCalibraFillHisto;
class TClonesArray;

class AliHLTTRDCalibHistoComponent : public AliHLTProcessor
{
public:
  AliHLTTRDCalibHistoComponent();
  virtual ~AliHLTTRDCalibHistoComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();
	
protected:
	
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
	
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
	       AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;

  int Configure(const char* arguments);
  int SetParams();
  int TakeHistos(Int_t nTimeBins);

private:
  /** copy constructor prohibited */
  AliHLTTRDCalibHistoComponent(const AliHLTTRDCalibHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTRDCalibHistoComponent& operator=(const AliHLTTRDCalibHistoComponent&);
  void FormOutput();

  AliHLTUInt32_t fOutputSize;    // output size
  AliHLTUInt32_t fSpec;         // accumulated specification
  TClonesArray* fTracksArray;    // array containing the input
  TObjArray* fOutArray;          // array containing the output
  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;  // TRD calibration object	
  Bool_t fSavedTimeBins;         // already saved the number of time bins?
  TObjArray *fTrgStrings;        // name of trigger classes to accept or reject
  Int_t  fAccRejTrg;             // do we actually accept or reject the trigger strings?
  Int_t fMinClusters;           // minimal number of clusters/tracklet accepted to fill histos
  Int_t fMinTracklets;          // minimal number of tracklets/track accepted to fill histos
  Bool_t fTakeAllEvents;         // take all events, disregarding the triggers
  
  ClassDef(AliHLTTRDCalibHistoComponent, 2)

};
#endif

