// -*- Mode: C++ -*-

#ifndef ALIHLTTPCCLUSTERTRANSFORMATIONMERGERCOMPONENT_H
#define ALIHLTTPCCLUSTERTRANSFORMATIONMERGERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCClusterTransformationMergerComponent.h
    @author Sergey Gorbunov
    @date   
    @brief
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class AliHLTTPCClusterTransformation;


/**
 * @class AliHLTTPCClusterTransformationMergerComponent
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCClusterTransformationMergerComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTTPCClusterTransformationMergerComponent();           
  /** destructor */
  virtual ~AliHLTTPCClusterTransformationMergerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
      
  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();							     
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();					     
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn(); 							   
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);
  
protected:
	
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing capabilities of the component. 

  int DoInit( int argc, const char** argv );
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  int DoDeinit();
  
  int Reconfigure(const char* cdbEntry, const char* chainId);  

  int ScanConfigurationArgument(int argc, const char** argv);

  using AliHLTProcessor::DoEvent;
  
private:
   
  int fCumulative;
  int fTotalInputs;
  TObject* fObj;
          
  /** copy constructor prohibited */
  AliHLTTPCClusterTransformationMergerComponent(const AliHLTTPCClusterTransformationMergerComponent&);

  /** assignment operator prohibited */
  AliHLTTPCClusterTransformationMergerComponent& operator=(const AliHLTTPCClusterTransformationMergerComponent&);

  ClassDef(AliHLTTPCClusterTransformationMergerComponent, 0)
};

#endif
