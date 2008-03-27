// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCNOISEMAPCOMPONENT_H
#define ALIHLTTPCNOISEMAPCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCNoiseMapComponent.h
    @author Kalliopi Kanaki
    @date   
    @brief  Component for Noise Map
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

//forward declarations
class AliHLTTPCDigitReader;
class TH2;

/**
 * @class AliHLTTPCNoiseMapComponent
 * Implementation of the component to fill histograms with TPC noise by request.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * It reads the data pad by pad and fills histograms. The output is unpacked and 
 * sent to the clulsterfinder.
 * 
 * The component has the following component arguments:
 * - adc-threshold   ADC count threshold for zero suppression.
 *
 * - rms-threshold   RMS threshold for zero suppression.
 *          
 * - first-timebin   The first timebin for zero suppression
 *
 * - last-timebin    The last timebin for zero suppression
 *
 * - occupancy-limit Minimum number of timebins with signal
 *
 * - sort-pads Flag to switch on pad sorting(needed by the SORTED clusterfinder)
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCNoiseMapComponent : public AliHLTProcessor {
    
   public:
   
   /** standard constructor */    
   AliHLTTPCNoiseMapComponent();           
   /** destructor */
   virtual ~AliHLTTPCNoiseMapComponent();


      // Public functions to implement AliHLTComponent's interface.
      // These functions are required for the registration process
      
      /** interface function, see @ref AliHLTComponent for description */
      const char* GetComponentID();							     
      /** interface function, see @ref AliHLTComponent for description */
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
      /** interface function, see @ref AliHLTComponent for description */
      AliHLTComponentDataType GetOutputDataType();					     
      /** interface function, see @ref AliHLTComponent for description */
      int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   
      /** interface function, see @ref AliHLTComponent for description */
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
      /** interface function, see @ref AliHLTComponent for description */
      AliHLTComponent* Spawn(); 							   
      /** function for acting on the saving and cleaning histograms, after they are filled */
      void SaveAndResetHistograms();
  
   protected:
	
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing capabilities of the component. 

      int DoInit( int argc, const char** argv );
      int DoDeinit();
      int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		   AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		   AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
      int Reconfigure(const char* cdbEntry, const char* chainId);

      using AliHLTProcessor::DoEvent;

   private:
   
      int Configure(const char* arguments);
          
      /** copy constructor prohibited */
      AliHLTTPCNoiseMapComponent(const AliHLTTPCNoiseMapComponent&);

      /** assignment operator prohibited */
      AliHLTTPCNoiseMapComponent& operator=(const AliHLTTPCNoiseMapComponent&);

      /** the reader object for data decoding */
      AliHLTUInt32_t fSpecification;  //!transient
      //AliHLTUInt8_t fMinPartition;    //!transient
      //AliHLTUInt8_t fMaxPartition;    //!transient

      Bool_t fNoiseMap;    //!transient
      Bool_t fIsPacked;    //!transient   
      Bool_t fIsUnpacked;  //!transient
      
      Int_t  fCurrentPartition; //!transient
      Int_t  fCurrentRow;       //!transient
      Int_t  rowOffset;         //!transient
      
      TH2 *fHistSideC; //!transient    

      ClassDef(AliHLTTPCNoiseMapComponent, 0)
    };

#endif
