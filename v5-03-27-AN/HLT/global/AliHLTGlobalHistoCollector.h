// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALHISTOCOLLECTOR_H
#define ALIHLTGLOBALHISTOCOLLECTOR_H

//* This file is property of and copyright by the ALICE HLT Project        * 
  //* ALICE Experiment at CERN, All rights reserved.                         *
    //* See cxx source for full Copyright notice                               *

    /** @file   AliHLTGlobalHistoCollector.h
	@author Kalliopi Kanaki, Kenneth Aamodt
	@date   
	@brief  Component for acting upon histograms
    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponentBenchmark.h"
#include <vector>

class TH1;
class TH2;

/**
 * @class AliHLTGlobalHistoCollector
 * Implementation of the component to read histograms from other
 * components and add, divide etc.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 *  
 * The component has the following component arguments:
 *
 * -sum-noise-histograms Loops over the output of TPCNoiseMap and sums the partition histograms
 *  They are sorted per TPC side.
 *
 * -sum-krypton-histograms Loops over the output of the krypton CF and sums the histograms
 * (it will become obsolete, when the next option does all the work)
 *
 * -use-general It will become the standard general option for summing histograms
 *
 * -ignore-specification It ignores the last part of the histogram name, if it has 
 * the form "_Slice_%.2d%.2d_Partition_%.2d%.2d, minSlice, maxSlice, minPartition, maxPartition".
 * It keeps the first part of the hist name and uses it to name the summed histogram.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTGlobalHistoCollector : public AliHLTProcessor {
    
public:
  struct AliHLTGlobalHCInstance
  {
    TObject *fObject;
    AliHLTUInt32_t fHLTSpecification;
  };

  struct AliHLTGlobalHCCollection
  {
  public:
    AliHLTGlobalHCCollection():fMergedObject(0),fHLTDataType(kAliHLTVoidDataType),fInstances(),fNeedToMerge(0){}
    AliHLTGlobalHCCollection( const AliHLTGlobalHCCollection &x):fMergedObject(x.fMergedObject),fHLTDataType(x.fHLTDataType),fInstances(x.fInstances),fNeedToMerge(x.fNeedToMerge){}
    AliHLTGlobalHCCollection &operator=( const AliHLTGlobalHCCollection &x){
      if( &x!=this ){
	fMergedObject = x.fMergedObject;
	fHLTDataType = x.fHLTDataType;
	fInstances = x.fInstances;   
	fNeedToMerge = x.fNeedToMerge;
      }
      return *this;
    }

   ~AliHLTGlobalHCCollection(){}
    
    TObject *fMergedObject;
    AliHLTComponentDataType fHLTDataType;
    std::vector<AliHLTGlobalHCInstance> fInstances;
    bool fNeedToMerge;
  };

  /** standard constructor */    
  AliHLTGlobalHistoCollector();           
  /** destructor */
  virtual ~AliHLTGlobalHistoCollector();

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
  AliHLTGlobalHistoCollector(const AliHLTGlobalHistoCollector&);

  /** assignment operator prohibited */
  AliHLTGlobalHistoCollector& operator=(const AliHLTGlobalHistoCollector&);

  void Clear(); // reset the store

  AliHLTUInt32_t fUID;// uID of the component

  std::vector<AliHLTGlobalHCCollection> fStore;
  AliHLTComponentBenchmark fBenchmark;// benchmark

};

#endif
