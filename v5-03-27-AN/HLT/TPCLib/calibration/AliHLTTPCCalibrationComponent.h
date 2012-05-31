//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCCALIBRATIONCOMPONENT_H
#define ALIHLTTPCCALIBRATIONCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibrationComponent.h
    @author Kalliopi Kanaki
    @date   
    @brief  
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTCalibrationProcessor.h"
class AliTPCcalibTime;
class AliTPCcalibTimeGain;
class AliHLTTPCAnalysisTaskcalib;
class AliESDEvent;
class TObjArray;

/**
 * @class AliHLTTPCCalibrationComponent
 * 
 * Interface of the offline algorithm (AliTPCcalibTimeGain) for updating
 * the TPC gain as a function of time.
 * 
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCCalibrationComponent : public AliHLTCalibrationProcessor
    {
    public:
      /** constructor */
      AliHLTTPCCalibrationComponent();
      /** destructor */
      virtual ~AliHLTTPCCalibrationComponent();

      const char* GetComponentID();
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
      AliHLTComponentDataType GetOutputDataType();
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
      AliHLTComponent* Spawn();

    protected:

      using AliHLTCalibrationProcessor::ProcessCalibration;
      using AliHLTCalibrationProcessor::ShipDataToFXS;
           
      /** Initialize the calibration component. */
      Int_t InitCalibration();

      /** Scan commandline arguments of the calibration component. */
      Int_t ScanArgument( Int_t argc, const char** argv );

      /** DeInitialize the calibration component. */
      Int_t DeinitCalibration();

      /** Process the data in the calibration component. */
      Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

      /** Ship the data to the FXS at end of run or eventmodulo. */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    private:
      /** copy constructor prohibited */
      AliHLTTPCCalibrationComponent(const AliHLTTPCCalibrationComponent&);
      /** assignment operator prohibited */
      AliHLTTPCCalibrationComponent& operator=(const AliHLTTPCCalibrationComponent&);

      AliHLTTPCAnalysisTaskcalib *fCalibTask;
      AliTPCcalibTime            *fCalibTime;
      AliTPCcalibTimeGain        *fCalibTimeGain;
     
      AliESDEvent                *fESDEvent;  //!transient
      TObjArray                  *fSeedArray; //!transient
      
      AliHLTUInt8_t  fMinPartition;   // see above
      AliHLTUInt8_t  fMaxPartition;   // see above
      AliHLTUInt8_t  fMinSlice;       // see above
      AliHLTUInt8_t  fMaxSlice;       // see above
      AliHLTUInt32_t fSpecification;  // see above

      /** Analyze calibration data before shipping to FXS */
      Bool_t fEnableAnalysis;  // see above

      ClassDef(AliHLTTPCCalibrationComponent, 0)
    };
#endif
