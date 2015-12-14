// -*- Mode: C++ -*-

#ifndef ALIHLTROOTOBJECTMERGERCOMPONENT_H
#define ALIHLTROOTOBJECTMERGERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTRootObjectMergerComponent.h
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

class TObject;
class TList;

/**
 * @class AliHLTRootObjectMergerComponent
 * @ingroup alihlt_tpc_components
 */

class AliHLTRootObjectMergerComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTRootObjectMergerComponent();           
  /** destructor */
  virtual ~AliHLTRootObjectMergerComponent();

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
  void* cleanup(void*);
  
  int Reconfigure(const char* cdbEntry, const char* chainId);  

  int ScanConfigurationArgument(int argc, const char** argv);

  using AliHLTProcessor::DoEvent;
  
private:
  int BuildMergeList(TObject*& returnObj, TList*& mergeList, TObject* obj);
  void ClearBuffers(void* buffer, bool isMergeObjectStruct = false);
	
  struct MergeObjectStruct
  {
	  TObject* fObject;
	  TList* fList;
  };
   
  int fCumulative;
  int fTotalInputs;
  TObject* fObj;
  int fQueueDepth;
  int fAsyncProcess;
  
  bool fDataTypeSet;
  AliHLTComponentDataType fDataType;
  
  void* MergeObjects(void*);
          
  /** copy constructor prohibited */
  AliHLTRootObjectMergerComponent(const AliHLTRootObjectMergerComponent&);

  /** assignment operator prohibited */
  AliHLTRootObjectMergerComponent& operator=(const AliHLTRootObjectMergerComponent&);
  
  AliHLTAsyncMemberProcessor<AliHLTRootObjectMergerComponent> fAsyncProcessor; //Processor for asynchronous processing

  ClassDef(AliHLTRootObjectMergerComponent, 0)
};

#endif
