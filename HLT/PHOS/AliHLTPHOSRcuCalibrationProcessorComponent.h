//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUCALIBRATIONPROCESSORCOMPONENT_H
#define ALIHLTPHOSRCUCALIBRATIONPROCESSORCOMPONENT_H

#include "AliHLTCalibrationProcessor.h"

class AliHLTPHOSRcuCalibrationProcessor;
class AliHLTPHOSSharedMemoryInterface;
class TObjArray;

class AliHLTPHOSRcuCalibrationProcessorComponent: public AliHLTCalibrationProcessor
//class AliHLTPHOSRcuCalibrationProcessorComponent:  public AliHLTPHOSBase, public AliHLTProcessor
{
public:

  /** constructor */
  AliHLTPHOSRcuCalibrationProcessorComponent();

 
  /** destructor */
  virtual ~AliHLTPHOSRcuCalibrationProcessorComponent();
      
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
  using  AliHLTCalibrationProcessor::ProcessCalibration;
  Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  /** Ship the data to the FXS at end of run or eventmodulo. */
  using AliHLTCalibrationProcessor::ShipDataToFXS; 
  Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

private:
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTPHOSRcuCalibrationProcessorComponent(const AliHLTPHOSRcuCalibrationProcessorComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTPHOSRcuCalibrationProcessorComponent& operator=(const AliHLTPHOSRcuCalibrationProcessorComponent&);
  TObjArray* fCalibDataPtr;                                  //! transient
  AliHLTPHOSRcuCalibrationProcessor* fRcuCalibProcessorPtr;   /**<Pointer to a phos histoproducer object*/
  AliHLTPHOSSharedMemoryInterface *fShmPtr; // Interface to read altro channel data from shared memory

};

#endif
