// -*- Mode: C++ -*-

#ifndef ALIHLTTPCCLUSTERTRANSFORMATIONPREPARECOMPONENT_H
#define ALIHLTTPCCLUSTERTRANSFORMATIONPREPARECOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCClusterTransformationPrepareComponent.h
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
#include "AliHLTAsyncMemberProcessor.h"

class AliHLTTPCClusterTransformation;
class AliHLTTPCFastTransformObject;
class AliCDBEntry;

/**
 * @class AliHLTTPCClusterTransformationPrepareComponent
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCClusterTransformationPrepareComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTTPCClusterTransformationPrepareComponent();           
  /** destructor */
  virtual ~AliHLTTPCClusterTransformationPrepareComponent();

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
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int DoDeinit();
  
  int Reconfigure(const char* cdbEntry, const char* chainId);  

  int ScanConfigurationArgument(int argc, const char** argv);

  using AliHLTProcessor::DoEvent;
  
private:
   
          
  /** copy constructor prohibited */
  AliHLTTPCClusterTransformationPrepareComponent(const AliHLTTPCClusterTransformationPrepareComponent&);

  /** assignment operator prohibited */
  AliHLTTPCClusterTransformationPrepareComponent& operator=(const AliHLTTPCClusterTransformationPrepareComponent&);

  Int_t fMinInitSec;	//Min sector for parallel initialization
  Int_t fMaxInitSec;	//Max sector for parallel initialization
  bool fNoInitialObject; //Do not create an initial object at doinit

  static const char* fgkOCDBEntryClusterTransformation; //!transient
  static AliHLTTPCClusterTransformation fgTransform; //!transient
  static Bool_t fgTimeInitialisedFromEvent; //!transient
  
  AliHLTTPCFastTransformObject* fTmpFastTransformObject;
  AliHLTTPCFastTransformObject* GenerateFastTransformObject();
  void* AsyncGenerateFastTransformObject(void*);
  
  AliHLTAsyncMemberProcessor<AliHLTTPCClusterTransformationPrepareComponent> fAsyncProcessor;
  int fAsyncProcessorQueueDepth;
  
  AliCDBEntry* fNewCalibObject;

  ClassDef(AliHLTTPCClusterTransformationPrepareComponent, 0)
};

#endif
