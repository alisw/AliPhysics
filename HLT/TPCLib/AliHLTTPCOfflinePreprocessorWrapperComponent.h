// -*- Mode: C++ -*-

#ifndef ALIHLTTPCOFFLINEPREPROCESSORWRAPPERCOMPONENT_H
#define ALIHLTTPCOFFLINEPREPROCESSORWRAPPERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCOfflinePreprocessorWrapperComponent.h
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

class TObjArray;
class AliCDBEntry;

/**
 * @class AliHLTTPCOfflinePreprocessorWrapperComponent
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCOfflinePreprocessorWrapperComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTTPCOfflinePreprocessorWrapperComponent();           
  /** destructor */
  virtual ~AliHLTTPCOfflinePreprocessorWrapperComponent();

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
  AliHLTTPCOfflinePreprocessorWrapperComponent(const AliHLTTPCOfflinePreprocessorWrapperComponent&);

  /** assignment operator prohibited */
  AliHLTTPCOfflinePreprocessorWrapperComponent& operator=(const AliHLTTPCOfflinePreprocessorWrapperComponent&);
  
  TObjArray* GetDataContainer(TObject* obj);

  AliCDBEntry* RunPreprocessor(TObjArray* dataContainer);
  void* AsyncRunPreprocessor(void*);
  
  AliHLTAsyncMemberProcessor<AliHLTTPCOfflinePreprocessorWrapperComponent> fAsyncProcessor;
  int fAsyncProcessorQueueDepth;
  int fAsyncProcess;
  
  ClassDef(AliHLTTPCOfflinePreprocessorWrapperComponent, 0)
};

#endif
