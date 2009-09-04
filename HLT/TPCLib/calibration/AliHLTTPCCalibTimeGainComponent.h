//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCCALIBTIMEGAINCOMPONENT_H
#define ALIHLTTPCCALIBTIMEGAINCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibTimeGainComponent.h
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  A calibration component for interfacing the offline calculation of TPC gain correction
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTCalibrationProcessor.h"

class AliTPCcalibTimeGain;
//class AliExternalTrackParam;
class AliESDEvent;
class TObjArray;

/**
 * @class AliHLTTPCCalibTimeGainComponent
 * 
 * Interface of the offline algorithm (AliTPCcalibTimeGain) for updating
 * the TPC gain as a function of time.
 * 
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCCalibTimeGainComponent : public AliHLTCalibrationProcessor
    {
    public:
      /** constructor */
      AliHLTTPCCalibTimeGainComponent();
      /** destructor */
      virtual ~AliHLTTPCCalibTimeGainComponent();

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
      AliHLTTPCCalibTimeGainComponent(const AliHLTTPCCalibTimeGainComponent&);
      /** assignment operator prohibited */
      AliHLTTPCCalibTimeGainComponent& operator=(const AliHLTTPCCalibTimeGainComponent&);

      AliTPCcalibTimeGain *fCalibTimeGain; //!transient
      AliESDEvent         *fESDEvent;      //!transient
      TObjArray           *fSeedArray;     //!transient
      
      AliHLTUInt8_t  fMinPartition;   // see above
      AliHLTUInt8_t  fMaxPartition;   // see above
      AliHLTUInt8_t  fMinSlice;       // see above
      AliHLTUInt8_t  fMaxSlice;       // see above
      AliHLTUInt32_t fSpecification;  // see above

      /** Analyze calibration data before shipping to FXS */
      Bool_t fEnableAnalysis;  // see above

      ClassDef(AliHLTTPCCalibTimeGainComponent, 0)
    };
#endif
