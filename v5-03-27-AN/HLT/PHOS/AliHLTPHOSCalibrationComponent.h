//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSCALIBRATIONCOMPONENT_H
#define ALIHLTPHOSCALIBRATIONCOMPONENT_H

#include "AliHLTCalibrationProcessor.h"

class AliHLTPHOSEmcCalibData;

class AliHLTPHOSCalibrationComponent: public AliHLTCalibrationProcessor
//class AliHLTPHOSCalibrationComponent:  public AliHLTPHOSBase, public AliHLTProcessor
{
public:

  /** constructor */
  AliHLTPHOSCalibrationComponent();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTPHOSCalibrationComponent(const AliHLTPHOSCalibrationComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTPHOSCalibrationComponent& operator=(const AliHLTPHOSCalibrationComponent&);
  /** destructor */
  virtual ~AliHLTPHOSCalibrationComponent();
      
  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

protected:
      
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
      
  /** Initialize the calibration component. */
  Int_t InitCalibration();

  /** Scan commandline arguments of the calibration component. */
  Int_t ScanArgument( Int_t argc, const char** argv );

  /** DeInitialize the calibration component. */
  Int_t DeinitCalibration();

  /** Process the data in the calibration component. */
  //Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
  virtual Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData,
			  const AliHLTComponent_BlockData* blocks,
			  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
			  AliHLTUInt32_t& size,
			  vector<AliHLTComponent_BlockData>& outputBlocks);
 
  /** Ship the data to the FXS at end of run or eventmodulo. */
  Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

private:

  AliHLTPHOSEmcCalibData *fEmcCalibData;

};

#endif
