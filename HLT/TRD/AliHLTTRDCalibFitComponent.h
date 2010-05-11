//-*- Mode: C++ -*-
// $Id: AliHLTTRDCalibFitComponent.h 40269 2010-04-08 22:08:53Z richterm $

#ifndef ALIHLTTRDCALIBFITCOMPONENT_H
#define ALIHLTTRDCALIBFITCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTRDCalibFitComponent.h
//  @author 
//  @date   
//  @brief  Declaration of a TRDCalibration component. 
// 


#include "AliHLTCalibrationProcessor.h"
class AliCDBManager;
class AliTRDCalibraFillHisto;
class TClonesArray;

/**
 * @class AliHLTTRDCalibFitComponent
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
class AliHLTTRDCalibFitComponent : public AliHLTCalibrationProcessor
{
public:
  AliHLTTRDCalibFitComponent();
  virtual ~AliHLTTRDCalibFitComponent();

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
	
  virtual Int_t InitCalibration();
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
  AliHLTTRDCalibFitComponent(const AliHLTTRDCalibFitComponent&);
  /** assignment operator prohibited */
  AliHLTTRDCalibFitComponent& operator=(const AliHLTTRDCalibFitComponent&);

  AliHLTUInt32_t fOutputSize;    // output size
  TObjArray* fOutArray;          // array containing the output
  TObjArray* fAfterRunArray;     // array with after run processing output 
  Bool_t fIncSM[18];             // array for telling which super module was already added
  Int_t fNoOfSM;                 // number of known SM
  Int_t fNoOfIncSM;              // number of SM already added
  
  ClassDef(AliHLTTRDCalibFitComponent, 2)

};
#endif

