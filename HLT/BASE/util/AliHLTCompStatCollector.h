// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTCOMPSTATCOLLECTOR_H
#define ALIHLTCOMPSTATCOLLECTOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTCompStatCollector.h
    @author Matthias Richter
    @date   
    @brief  Collector component for the component statistics information.
*/

#include "AliHLTProcessor.h"

class TStopwatch;
class TH1F;
class TH2F;
class TH2C;
class TTree;
class TFolder;

/**
 * @class AliHLTCompStatCollector
 * Collector component for the statistics entries produced by the
 * AliHLTComponent base class.
 *
 * <h2>General properties:</h2>
 * This components collects all data blocks of types
 * ::kAliHLTDataTypeComponentStatistics and ::kAliHLTDataTypeComponentTable
 * which can be produced by the AliHLTComponent base class for every component
 * and event. Component statistics entries are data blocks of
 * ::AliHLTComponentStatistics arrays containing a couple of informations
 * about each component. The information is extracted and stored into a
 * TTree. The component table entries (AliHLTComponentTable structs) are sent
 * on SOR and EOR events and are arranged in a TFolder hierarchy.
 *
 * Component ID: \b StatisticsCollector                                 <br>
 * Library: \b libAliHLTUtil.so					        <br>
 * Input Data Types: kAliHLTDataTypeComponentStatistics                 <br>
 * Output Data Types: kAliHLTDataTypeHistogram, kAliHLTDataTypeTNtuple  <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Output size:</h2>
 *
 * @ingroup alihlt_util_components
 */
class AliHLTCompStatCollector : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTCompStatCollector();
  /** destructor */
  virtual ~AliHLTCompStatCollector();

  const char* GetComponentID() {return "StatisticsCollector";};
  AliHLTComponent* Spawn() {return new AliHLTCompStatCollector;}
  void GetInputDataTypes( vector<AliHLTComponentDataType>& );
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

 protected:
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  
  using AliHLTProcessor::DoEvent;

 private:
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTCompStatCollector(const AliHLTCompStatCollector&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTCompStatCollector& operator=(const AliHLTCompStatCollector&);

  /**
   * Reset all filling variables and lists.
   */
  void ResetFillingVariables();

  /**
   * Fill the lists from the component statistics block.
   */
  int FillVariablesSorted(void* ptr, int size);

  /** delete all internal objects */
  void ClearAll();

  /**
   * Remove entries from the parent list if they occur further down in the
   * hierarchy.
   */
  int RemoveRecurrence(TFolder* pRoot) const;

  /** event cycle timer */
  TStopwatch* fpTimer; //!transient

  /** top folder */
  TFolder* fpFolder; //!transient

  /** statistics tree */
  TTree* fpStatTree; //!transient

  /** branch filling variable */
  Float_t fCycleTime; //!transient
  /** branch filling variable */
  Int_t fNofSets; //!transient
  /** array size */
  UInt_t fArraySize; //!transient
  /** current position */
  UInt_t fPosition; //!transient
  /** branch filling variable */
  UInt_t* fpLevelArray; //!transient
  /** branch filling variable */
  UInt_t* fpSpecArray; //!transient
  /** branch filling variable */
  UInt_t* fpBlockNoArray; //!transient
  /** branch filling variable */
  UInt_t* fpIdArray; //!transient
  /** branch filling variable */
  UInt_t* fpTimeArray; //!transient
  /** branch filling variable */
  UInt_t* fpCTimeArray; //!transient
  /** branch filling variable */
  UInt_t* fpInputBlockCountArray; //!transient
  /** branch filling variable */
  UInt_t* fpTotalInputSizeArray; //!transient
  /** branch filling variable */
  UInt_t* fpOutputBlockCountArray; //!transient
  /** branch filling variable */
  UInt_t* fpTotalOutputSizeArray; //!transient

  /** const base of GetOutputSize, updated on error in DoEvent */
  int fSizeEstimator; //! transient

  ClassDef(AliHLTCompStatCollector, 2)
};
#endif
