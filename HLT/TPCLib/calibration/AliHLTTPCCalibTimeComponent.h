//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCCALIBTIMECOMPONENT_H
#define ALIHLTTPCCALIBTIMECOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibTimeComponent.h
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  A calibration component for interfacing the offline calculation of TPC drift velocity correction
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTCalibrationProcessor.h"

class AliTPCcalibTime;
class AliExternalTrackParam;
class AliESDEvent;
class TObjArray;

/**
 * @class AliHLTTPCCalibTimeComponent
 * 
 * Interface of the offline algorithm (AliTPCcalibTime) for correcting the 
 * drift velocity for changes of p and T as a function of
 * time.
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCCalibTimeComponent : public AliHLTCalibrationProcessor
    {
    public:
      /** constructor */
      AliHLTTPCCalibTimeComponent();
      /** destructor */
      virtual ~AliHLTTPCCalibTimeComponent();
      
      const char* GetComponentID();
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
      AliHLTComponentDataType GetOutputDataType();
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
      AliHLTComponent* Spawn();

    protected:

      using AliHLTCalibrationProcessor::ProcessCalibration;
      using AliHLTCalibrationProcessor::ShipDataToFXS;
      
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
      Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

      /** Ship the data to the FXS at end of run or eventmodulo. */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    private:
      /** copy constructor prohibited */
      AliHLTTPCCalibTimeComponent(const AliHLTTPCCalibTimeComponent&);
      /** assignment operator prohibited */
      AliHLTTPCCalibTimeComponent& operator=(const AliHLTTPCCalibTimeComponent&);

      AliTPCcalibTime *fCalibTime; //!transient
      AliESDEvent     *fESDEvent;  //!transient
      TObjArray       *fSeedArray; //!transient
      
      AliHLTUInt8_t  fMinPartition;  // see above
      AliHLTUInt8_t  fMaxPartition;  // see above
      AliHLTUInt8_t  fMinSlice;      // see above
      AliHLTUInt8_t  fMaxSlice;      // see above
      AliHLTUInt32_t fSpecification; // see above

      /** Analyze calibration data before shipping to FXS */
      Bool_t fEnableAnalysis;  // see above

      ClassDef(AliHLTTPCCalibTimeComponent, 0)
    };
#endif
