//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTGLOBALFLATESDCONVERTERCOMPONENT_H
#define ALIHLTGLOBALFLATESDCONVERTERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTGlobalFlatEsdConverterComponent.h
//  @author Matthias Richter
//  @date   
//  @brief  Global ESD converter component.
//  @note

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
#include <vector>


/**
 * @class AliHLTGlobalFlatEsdConverterComponent
 * Global collector for information designated for the HLT ESD.
 *
 * componentid: \b GlobalEsdConverter <br>
 * componentlibrary: \b libAliHLTGlobal.so <br>
 * Arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -notree                                                          <br>
 *      write ESD directly to output (::kAliHLTDataTypeESDObject)
 *      this has been made the default behavior in Sep 2008.
 * \li -tree                                                            <br>
 *      write ESD directly to TTree and to output (::kAliHLTDataTypeESDTree)
 * \li -skipobject=name1,name2,...                                   <br>
 *      comma separated list of ESD object names to be skipped, default is
 *      AliESDZDC,AliESDFMD,Cascades,Kinks,AliRawDataErrorLogs,AliESDACORDE
 *      leave blank to disable the option
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTGlobalFlatEsdConverterComponent : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTGlobalFlatEsdConverterComponent();
  /** destructor */
  virtual ~AliHLTGlobalFlatEsdConverterComponent();

  // interface methods of base class
  const char* GetComponentID() {return "GlobalFlatEsdConverter";};
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& list);

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn() {return new AliHLTGlobalFlatEsdConverterComponent;}

 protected:
  // interface methods of base class
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks );

  using AliHLTProcessor::DoEvent;


 private:
  /** copy constructor prohibited */
  AliHLTGlobalFlatEsdConverterComponent(const AliHLTGlobalFlatEsdConverterComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalFlatEsdConverterComponent& operator=(const AliHLTGlobalFlatEsdConverterComponent&);

  /**
   * (Re)Configure from the CDB
   * Loads the following objects:
   * - HLT/ConfigHLT/SolenoidBz
   */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);

  /// verbosity level
  int fVerbosity; //!transient

protected:

  static const Int_t fkNPartition = 36*6;           // number of patches in TPC

  AliHLTComponentBenchmark fBenchmark; // benchmark
  
  Bool_t fProduceFriend; // should it produce the flat ESD friend

  ClassDef(AliHLTGlobalFlatEsdConverterComponent, 0)
};
#endif
