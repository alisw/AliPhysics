//-*- Mode: C++ -*-
#ifndef ALIHLTTPCCALIBPEDESTALCOMPONENT_H
#define ALIHLTTPCCALIBPEDESTALCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibPedestalComponent.h
    @author Jochen Thaeder
    @date   
    @brief  A pedestal calibration component for the TPC.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTCalibrationProcessor.h"
#include "AliHLTTPCDefinitions.h"

class AliTPCRawStream;
class AliRawReaderMemory;
class AliTPCCalibPedestal;

/**
 * @class AliHLTTPCCalibPedestalComponent
 * 
 * This class is the calibration component for the AliTPCCalibPedestal class 
 * used for pedestal calibration of the TPC. 
 * 
 * It inherits from the AliHLTCalibrationProcessor and uses the high-level 
 * interface. The output is the class AliTPCCalibPedestal as a TObject.
 *
 * The component has the following component arguments:
 *   -rcuformat <old/new>  : Wether to use old or new rcuformat ( default is new )
 *   -enableanalysis       : Wether to enable analyis before shipping data to FXS
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCalibPedestalComponent : public AliHLTCalibrationProcessor
    {
    public:
      /** constructor */
      AliHLTTPCCalibPedestalComponent();
      /** destructor */
      virtual ~AliHLTTPCCalibPedestalComponent();
      
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
      Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

      /** Ship the data to the FXS at end of run or eventmodulo. */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    private:
      /** copy constructor prohibited */
      AliHLTTPCCalibPedestalComponent(const AliHLTTPCCalibPedestalComponent&);
      /** assignment operator prohibited */
      AliHLTTPCCalibPedestalComponent& operator=(const AliHLTTPCCalibPedestalComponent&);

      /** The reader object for reading from memory */
      AliRawReaderMemory* fRawReader;                                              //!transient

      /** The reader object for reading TPC raw data */  
      AliTPCRawStream* fRawStream;                                                 //!transient

      /** Pedestal Calibration class */
      AliTPCCalibPedestal * fCalibPedestal;                                        //!transient
      
      /** Wether to use old RCU format */
      Bool_t fRCUFormat;                                                           // see above

      /** Minimum patch specifcation for this component */
      AliHLTUInt8_t fMinPatch;                                                     // see above

      /** Minimum patch specifcation for this component */
      AliHLTUInt8_t fMaxPatch;                                                     // see above

      /** The Specification for this component */
      AliHLTUInt32_t fSpecification;                                               // see above

      /** Analysze calibration data before shipping to FXS */
      Bool_t fEnableAnalysis;                                                      // see above

      ClassDef(AliHLTTPCCalibPedestalComponent, 1)

    };
#endif
