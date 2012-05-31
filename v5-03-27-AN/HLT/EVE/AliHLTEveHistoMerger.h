// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTEVEHISTOMERGER_H
#define ALIHLTEVEHISTOMERGER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
  //* ALICE Experiment at CERN, All rights reserved.                         *
    //* See cxx source for full Copyright notice                               *

    /** @file   AliHLTEveHistoMerger.h
	@author Kalliopi Kanaki, Kenneth Aamodt
	@date   
	@brief  Component for acting upon histograms
    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTDataTypes.h"
#include <vector>

class TH1;
class TH2;
class TObject;

/**
 * @class AliHLTEveHistoMerger
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

class AliHLTEveHistoMerger  {
    
public:
  struct AliHLTGlobalHCInstance
  {
    TObject *fObject;
    AliHLTUInt32_t fHLTSpecification;
  };

  struct AliHLTGlobalHCCollection
  {
  public:
    AliHLTGlobalHCCollection():fMergedObject(0),fInstances(),fNeedToMerge(0){}
    AliHLTGlobalHCCollection( const AliHLTGlobalHCCollection &x):fMergedObject(x.fMergedObject),fInstances(x.fInstances),fNeedToMerge(x.fNeedToMerge){}
    AliHLTGlobalHCCollection &operator=( const AliHLTGlobalHCCollection &x){
      if (this==&x) return *this;
      fMergedObject = x.fMergedObject;
      fInstances = x.fInstances;   
      fNeedToMerge = x.fNeedToMerge;
      return *this;
    }

   ~AliHLTGlobalHCCollection(){}
    
    TObject *fMergedObject;
    std::vector<AliHLTGlobalHCInstance> fInstances;
    bool fNeedToMerge;
  };

  /** standard constructor */    
  AliHLTEveHistoMerger();           
  /** destructor */
  virtual ~AliHLTEveHistoMerger();

  TObject* Process( const TObject* evtData, AliHLTUInt32_t spec);
  
private:
   
  /** copy constructor prohibited */
  AliHLTEveHistoMerger(const AliHLTEveHistoMerger&);

  /** assignment operator prohibited */
  AliHLTEveHistoMerger& operator=(const AliHLTEveHistoMerger&);

  void Clear(); // reset the store


  std::vector<AliHLTGlobalHCCollection> fStore;

};

#endif
