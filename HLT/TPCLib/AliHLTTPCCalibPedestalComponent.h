
#ifndef ALIHLTTPCCALIBPEDESTALCOMPONENT_H
#define ALIHLTTPCCALIBPEDESTALCOMPONENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCCalibPedestalComponent.h
    @author Jochen Thaeder
    @date   
    @brief  A pedestal calibration component for the TPC.
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCDefinitions.h"

class AliTPCRawStream;
class AliRawReaderMemory;
class AliTPCCalibPedestal;

/**
 * @class AliHLTTPCCalibPedestalComponent
 * 
 * This class is the component for the AliTPCCalibPedestal class used for 
 * pedestal calibration of the TPC. It uses the high-level interface and
 * the output is the TObject of AliTPCCalibPedestal.
 *
 * The component has the following component arguments:
 * -    The RCU format:  rcuformat  
 *      which can be either 
 *        - old ( used in the TPC Commissioning )
 *        - new
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCalibPedestalComponent : public AliHLTProcessor
    {
    public:
      /**
       * constructor 
       */
      AliHLTTPCCalibPedestalComponent();
      /** not a valid copy constructor, defined according to effective C++ style */
      AliHLTTPCCalibPedestalComponent(const AliHLTTPCCalibPedestalComponent&);
      /** not a valid assignment op, but defined according to effective C++ style */
      AliHLTTPCCalibPedestalComponent& operator=(const AliHLTTPCCalibPedestalComponent&);
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
      
      Int_t DoInit( int argc, const char** argv );
      Int_t DoDeinit();
      Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
   
    private:

      /** the reader object for reading from memory */
      AliRawReaderMemory* fRawReader;                                              //!transient

      /** the reader object for reading TPC raw data */  
      AliTPCRawStream* fRawStream;                                                 //!transient

      /** Pedestal Calibration class */
      AliTPCCalibPedestal * fCalibPedestal;                                        //!transient
      
      /** if use old RCU format */
      Bool_t fRCUFormat;                                                           // see above

      ClassDef(AliHLTTPCCalibPedestalComponent, 0)

    };
#endif
