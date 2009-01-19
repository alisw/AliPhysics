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
    @brief  Component for plotting TPC data and applying noise map
*/

#include "AliHLTProcessor.h"

//forward declarations
class AliHLTTPCDigitReader;
class TH1;
class TH2;
class AliTPCCalPad;

/**
 * @class AliHLTTPCNoiseMapComponent
 * 
 * Implementation of the component to fill histograms with TPC output and read
 * the noise map from HCDB by request.
 * 
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * It reads the data pad by pad and fills histograms per partition
 * 
 * The component has the following component arguments:
 * 
 * -read-noisemap  Reads the noise map from the HCDB (and plots it in a histogram)
 *
 * It returns an AliTPCCalPad object that can be used 
 * for other plots (e.g. noise subtraction in histograms) 
 *
 * It makes use of the class AliHLTTPCNoise (see header for usage).
 *
 * -reset-histograms Resets histograms
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
      
      /** interface function, see AliHLTComponent for description */
      const char* GetComponentID();							     
      /** interface function, see AliHLTComponent for description */
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
      /** interface function, see AliHLTComponent for description */
      AliHLTComponentDataType GetOutputDataType();					     
      /** interface function, see AliHLTComponent for description */
      int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   
      /** interface function, see AliHLTComponent for description */
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
      /** interface function, see AliHLTComponent for description */
      AliHLTComponent* Spawn(); 							   
      /** function for acting on the saving and cleaning histograms, after they are filled */
      void MakeHistosPublic();
  
   protected:
	
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing capabilities of the component. 

      int DoInit( int argc, const char** argv );
      int DoDeinit();
      int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
      int Reconfigure(const char* cdbEntry, const char* chainId);

      using AliHLTProcessor::DoEvent;

   private:
   
      int Configure(const char* arguments);
          
      /** copy constructor prohibited */
      AliHLTTPCNoiseMapComponent(const AliHLTTPCNoiseMapComponent&);

      /** assignment operator prohibited */
      AliHLTTPCNoiseMapComponent& operator=(const AliHLTTPCNoiseMapComponent&);
      
      void InitializeHistograms(UInt_t minSlice, UInt_t maxSlice, UInt_t minPartition, UInt_t maxPartition);
      void ResetHistograms();
      
      AliHLTUInt32_t fSpecification;  //!transient

      Bool_t fReadNoiseMap;    //!transient
      Bool_t fResetHistograms; //!transient      
      Bool_t fInitHist;        //!transient

      Int_t fCurrentRow; //!transient
      
      TH1 *fHistSignal;     //!transient 
     
      TH2 *fHistSideAMaxSignal;  //!transient 
      TH2 *fHistSideATotSignal;  //!transient 
      TH2 *fHistSideAPadRMS;     //!transient 
     
      TH2 *fHistSideCMaxSignal;  //!transient 
      TH2 *fHistSideCTotSignal;  //!transient 
      TH2 *fHistSideCPadRMS;	 //!transient 
     
      TH2 *fHistCDBMap;     //!transient 
           
      ClassDef(AliHLTTPCNoiseMapComponent, 4)
    };

#endif
