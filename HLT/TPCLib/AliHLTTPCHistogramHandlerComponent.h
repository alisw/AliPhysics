// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCHISTOGRAMHANDLERCOMPONENT_H
#define ALIHLTTPCHISTOGRAMHANDLERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
  //* ALICE Experiment at CERN, All rights reserved.                         *
    //* See cxx source for full Copyright notice                               *

    /** @file   AliHLTTPCHistogramHandlerComponent.h
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
#include <vector>

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
class AliHLTTPCHistogramHandlerComponent : public AliHLTProcessor {
    
public:
  struct AliHLTHistogramData
  {
    TH1 *fHistogram;
    UInt_t fMinSlice;
    UInt_t fMaxSlice;
    UInt_t fMinPartition;
    UInt_t fMaxPartition;
  };
  typedef struct AliHLTHistogramData AliHLTHistogramData; //!

  /** standard constructor */    
  AliHLTTPCHistogramHandlerComponent();           
  /** destructor */
  virtual ~AliHLTTPCHistogramHandlerComponent();

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
  AliHLTTPCHistogramHandlerComponent(const AliHLTTPCHistogramHandlerComponent&);

  /** assignment operator prohibited */
  AliHLTTPCHistogramHandlerComponent& operator=(const AliHLTTPCHistogramHandlerComponent&);

  /** the reader object for data decoding */
  AliHLTUInt32_t fSpecification;  //!transient
            
  Bool_t fNoiseHistograms;     //!transient
  Bool_t fKryptonHistograms;   //!transient
  Bool_t fUseGeneral;          //!transient
  Bool_t fIgnoreSpecification; //!transient 
      
  Int_t fSlice;  //!transient
      
  TH1 *fHistTH1Tmp;                //!transient  
  TH1 *fTotalClusterChargeIROCAll; //!transient
  TH1 *fTotalClusterChargeOROCAll; //!transient
  TH1 *fQMaxPartitionAll;          //!transient
  TH1 *fPlotQmaxROCAll;            //!transient
  TH1 *fNumberOfClusters;          //!transient
            
  TH2 *fHistTH2Tmp;       //!transient
  TH2 *fHistTPCSideAmax;  //!transient 
  TH2 *fHistTPCSideCmax;  //!transient  
  TH2 *fHistTPCSideAtot;  //!transient 
  TH2 *fHistTPCSideCtot;  //!transient  
  TH2 *fHistTPCSideArms;  //!transient 
  TH2 *fHistTPCSideCrms;  //!transient  

  vector<AliHLTHistogramData> fHistogramData;

  
  ClassDef(AliHLTTPCHistogramHandlerComponent, 4)
};

#endif
