#ifndef ALIHLTTPCNOISEMAPCOMPONENT_H
#define ALIHLTTPCNOISEMAPCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

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
#include "TH2.h"

class AliHLTTPCDigitReader;

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
        
   AliHLTTPCNoiseMapComponent();          //constructor	   
   virtual ~AliHLTTPCNoiseMapComponent(); //destructor


      // Public functions to implement AliHLTComponent's interface.
      // These functions are required for the registration process
      
      const char* GetComponentID();							   //interface function, see @ref AliHLTComponent for description  
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			   //interface function, see @ref AliHLTComponent for description  
      AliHLTComponentDataType GetOutputDataType();					   //interface function, see @ref AliHLTComponent for description  
      int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   //interface function, see @ref AliHLTComponent for description
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); //interface function, see @ref AliHLTComponent for description
      AliHLTComponent* Spawn(); 							   //interface function, see @ref AliHLTComponent for description
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
      AliHLTTPCDigitReader* fDigitReader;  //!transient
      AliHLTUInt32_t fSpecification;       //!transient
      AliHLTUInt8_t fMinPatch;             //!transient
      AliHLTUInt8_t fMaxPatch;             //!transient

      UInt_t fFirstTimeBin;
      UInt_t fLastTimeBin;
      UInt_t fNSigmaThreshold;
      UInt_t fSignalThreshold;
      UInt_t fMinimumNumberOfSignals;
      UInt_t fOldRCUFormat;
      Bool_t fSortPads;
      Bool_t fNoiseMap;
      Bool_t fIsPacked;    
      Bool_t fIsUnpacked;
      Int_t  fCurrentPatch;
      Int_t  fCurrentRow;
      
      TH2F *hPatch;     

      ClassDef(AliHLTTPCNoiseMapComponent, 0)  
    };

#endif
