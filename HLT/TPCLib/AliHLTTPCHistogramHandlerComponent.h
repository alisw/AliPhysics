// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCHISTOGRAMHANDLERCOMPONENT_H
#define ALIHLTTPCHISTOGRAMHANDLERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCHistogramHandlerComponent.h
    @author Kalliopi Kanaki
    @date   
    @brief  Component for acting upon histograms
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class TH1;
class TH2;

/**
 * @class AliHLTTPCHistogramHandlerComponent
 * Implementation of the component to read histograms from other
 * components and add, divide etc.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 *  
 * The component has the following component arguments:
 *
 * -sum-noise-histograms Loops over the output of TPCNoiseMap and adds the histograms
 *
 * It loops over histogram input and sums up the TPC histograms per side (at the moment).
 * 
 * @ingroup alihlt_tpc
 */
class AliHLTTPCHistogramHandlerComponent : public AliHLTProcessor {
    
   public:
   
   /** standard constructor */    
   AliHLTTPCHistogramHandlerComponent();           
   /** destructor */
   virtual ~AliHLTTPCHistogramHandlerComponent();

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
      AliHLTTPCHistogramHandlerComponent(const AliHLTTPCHistogramHandlerComponent&);

      /** assignment operator prohibited */
      AliHLTTPCHistogramHandlerComponent& operator=(const AliHLTTPCHistogramHandlerComponent&);

      /** the reader object for data decoding */
      AliHLTUInt32_t fSpecification;  //!transient
      
      
      Bool_t fNoiseHistograms;   //!transient
      Bool_t fKryptonHistograms; //!transient
 
      AliHLTUInt32_t fSpecificationTPCA; //!transient
      AliHLTUInt32_t fSpecificationTPCC; //!transient
      
      Int_t fSlice;  //!transient
      
      TH1 *fHistTH1Tmp;    //!transient           
      TH2 *fHistTH2Tmp;    //!transient
      TH2 *fHistTPCSideA;  //!transient	
      TH2 *fHistTPCSideC;  //!transient  

            
      ClassDef(AliHLTTPCHistogramHandlerComponent, 0)
    };

#endif
