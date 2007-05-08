

#ifndef ALIHLTTPCCALIBSIGNAlCOMPONENT_H
#define ALIHLTTPCCALIBSIGNALCOMPONENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibSignalComponent.h
    @author Jochen Thaeder, Sorina Popescu
    @date   
    @brief HLT Pulser calibration component for the TPC.
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCDefinitions.h"

class AliTPCRawStream;
class AliRawReaderMemory;
class AliTPCCalibSignal;

/**
 * @class AliHLTTPCCalibSignalComponent
 * 
 * This class is the component for the AliTPCCalibSignal class used for 
 * pedestal calibration of the TPC. It uses the high-level interface and
 * the output is the TObject of AliTPCCalibSignal.
 *
 * The component has the following component arguments:
 * -    The RCU format:  rcuformat  
 *      which can be either 
 *        - old ( used in the TPC Commissioning )
 *        - new
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCalibSignalComponent : public AliHLTProcessor
    {
    public:
      /**
       * constructor 
       */
      AliHLTTPCCalibSignalComponent();
      /** not a valid copy constructor, defined according to effective C++ style */
      AliHLTTPCCalibSignalComponent(const AliHLTTPCCalibSignalComponent&);
      /** not a valid assignment op, but defined according to effective C++ style */
      AliHLTTPCCalibSignalComponent& operator=(const AliHLTTPCCalibSignalComponent&);
      /** destructor */
      virtual ~AliHLTTPCCalibSignalComponent();
      
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
      
      Int_t DoInit( int argc, const char** argv );
      Int_t DoDeinit();
      Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
      
    private:

      /** the reader object for reading from memory */
      AliRawReaderMemory* fRawReader;                                              //!transient

      /** the reader object for reading TPC raw data */  
      AliTPCRawStream* fRawStream;                                                 //!transient

      /** Signal Calibration class */
      AliTPCCalibSignal * fCalibSignal;                                            //!transient
      
      /** if use old RCU format */
      Bool_t fRCUFormat;                                                           // see description

      ClassDef(AliHLTTPCCalibSignalComponent, 0)

    };
#endif
