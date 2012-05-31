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
class AliTPCcalibCalib;
class AliESDEvent;
class AliESDtrack;
class AliESDfriend;
class TObjArray;
class AliTPCclusterMI;

/**
 * @class AliHLTTPCCalibTimeComponent
 * 
 * Interface of the offline algorithm (AliTPCcalibTime) for estimating the 
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
      
      /** the name of the component registered in macros and configurations */
      const char* GetComponentID();
      
      /** the data types the component subscribes to in memory */
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list); 
      
      /** the data type the component produces */
      AliHLTComponentDataType GetOutputDataType();
      
      /** the size of the output data buffer */
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
      AliHLTComponent* Spawn();

    protected:

      using AliHLTCalibrationProcessor::ProcessCalibration;
      using AliHLTCalibrationProcessor::ShipDataToFXS;
      
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing
      // capabilities of the component. 
      
      /** Initialize the calibration component */
      Int_t InitCalibration();

      /** Scan commandline arguments of the calibration component */
      Int_t ScanConfigurationArgument( Int_t argc, const char** argv );

      /** Clean up memory at the end of the run */
      Int_t DeinitCalibration();

      /** Process the data in the calibration component, called once per event */
      Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
      
      /** inherited from AliHLTComponent: handle re-configuration event */
      int Reconfigure(const char* cdbEntry, const char* chainId);

      /** Ship the data to the FXS at end of run or event modulo (the first by default, the latter to be implemented if necessary). */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    private:
      
      /** copy constructor prohibited */
      AliHLTTPCCalibTimeComponent(const AliHLTTPCCalibTimeComponent&);
      
      /** assignment operator prohibited */
      AliHLTTPCCalibTimeComponent& operator=(const AliHLTTPCCalibTimeComponent&);

      AliTPCcalibTime  *fCalibTime; //!transient
      AliTPCcalibCalib *fCal;       //!transient
      AliESDEvent      *fESDevent;  //!transient
      AliESDtrack      *fESDtrack;  //!transient
      AliESDfriend     *fESDfriend; //!transient
      TObjArray        *fSeedArray; //!transient
      AliHLTUInt32_t    fOutputSize;// output size
      
      static const Int_t fkNPartition = 36*6; // number of partitions in TPC
      AliTPCclusterMI   *fPartitionClusters[fkNPartition];  //! arrays of cluster data for each TPC partition
      Int_t              fNPartitionClusters[fkNPartition]; //! number of clusters for each TPC partition

      /** the default configuration entry for this component */
      static const char* fgkOCDBEntry; //!transient

      ClassDef(AliHLTTPCCalibTimeComponent, 5)
    };
#endif
