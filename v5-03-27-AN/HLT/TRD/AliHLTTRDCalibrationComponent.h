//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCALIBRATIONCOMPONENT_H
#define ALIHLTTRDCALIBRATIONCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTRDCalibrationComponent.h
//  @author 
//  @date   
//  @brief  Declaration of a TRDCalibration component. 
// 


#include "AliHLTCalibrationProcessor.h"
class AliCDBManager;
class AliTRDCalibraFillHisto;
class TClonesArray;

/**
 * @class AliHLTTRDCalibrationComponent
 * @brief A TRDCalibration HLT processing component. 
 *
 * - @ref InitCalibration (optional)
 * - @ref ScanArgument (optional)
 * - @ref DeinitCalibration (optional)
 * - @ref ProcessCalibration
 * - @ref ShipDataToFXS
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn
 * @ingroup alihlt_tutorial
 */
class AliHLTTRDCalibrationComponent : public AliHLTCalibrationProcessor
{
public:
  AliHLTTRDCalibrationComponent();
  virtual ~AliHLTTRDCalibrationComponent();

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
	
  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;
	
  virtual Int_t InitCalibration();
  virtual Int_t ScanArgument(int argc, const char** argv);
  virtual Int_t DeinitCalibration();
  virtual Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData,
				   const AliHLTComponent_BlockData* blocks,
				   AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
				   AliHLTUInt32_t& size,
				   vector<AliHLTComponent_BlockData>& outputBlocks);
  /* 	virtual Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData); */
  virtual Int_t ShipDataToFXS(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  virtual Int_t  EORCalibration();
	
  using AliHLTCalibrationProcessor::ProcessCalibration;
  //using AliHLTCalibrationProcessor::ShipDataToFXS;

private:
  /** copy constructor prohibited */
  AliHLTTRDCalibrationComponent(const AliHLTTRDCalibrationComponent&);
  /** assignment operator prohibited */
  AliHLTTRDCalibrationComponent& operator=(const AliHLTTRDCalibrationComponent&);
  void FormOutput(Int_t param);

  AliHLTUInt32_t fOutputSize;    // output size
  TClonesArray* fTracksArray;    // array containing the input
  TObjArray* fOutArray;          // array containing the output
  TObjArray* fAfterRunArray;     // array with after run processing output 
  TObjArray* fDisplayArray;      // array with online display histos
  Bool_t fSavedTimeBins;         // already saved the number of time bins?
  TObjArray *fTrgStrings;        // name of trigger classes to accept or reject
  Int_t  fAccRejTrg;             // do we actually accept or reject the trigger strings?
  Int_t fMinClusters;           // minimal number of clusters/tracklet accepted to fill histos
  Int_t fMinTracklets;          // minimal number of tracklets/track accepted to fill histos
  Bool_t fTakeAllEvents;         // take all events, disregarding the triggers
  
  ClassDef(AliHLTTRDCalibrationComponent, 2)

};
#endif

