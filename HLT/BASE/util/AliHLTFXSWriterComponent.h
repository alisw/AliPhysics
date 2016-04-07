// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTFXSWRITERCOMPONENT_H
#define ALIHLTFXSWRITERCOMPONENT_H
//* This file is property of and copyright by the ALICE                    * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTCalibrationProcessor.h"
#include "TString.h"

class AliHLTFXSWriterComponent : public AliHLTCalibrationProcessor
{
 public:
  AliHLTFXSWriterComponent();
  virtual ~AliHLTFXSWriterComponent();

  virtual const char* GetComponentID() {return "FXSWriter";};
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long&, double&);
  virtual void GetOCDBObjectDescription( TMap* const targetArray);

  virtual AliHLTComponent* Spawn() {return new AliHLTFXSWriterComponent;}

 protected:
  virtual int InitCalibration();
  int ScanArgument( int argc, const char** argv ) {
    int result=ScanConfigurationArgument(argc, argv); return result>0?result-1:result;
  }
  virtual int DeinitCalibration();
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  virtual int ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  virtual int ScanConfigurationArgument(int argc, const char** argv);

private:
  AliHLTFXSWriterComponent(const AliHLTFXSWriterComponent&);
  AliHLTFXSWriterComponent& operator=(const AliHLTFXSWriterComponent&);
  
  TString fFXSName;
  TString fFXSDetector;
  AliHLTComponentDataType fDataType;
  bool fRootObject; //Deserialize the input and store as root file

  ClassDef(AliHLTFXSWriterComponent, 0)
};
#endif
