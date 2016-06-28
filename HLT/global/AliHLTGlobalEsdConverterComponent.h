//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALESDCONVERTERCOMPONENT_H
#define ALIHLTGLOBALESDCONVERTERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTGlobalEsdConverterComponent.h
//  @author Matthias Richter
//  @date   
//  @brief  Global ESD converter component.
//  @note

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
#include <vector>

// forward declarations
class AliESDEvent;
class AliESDfriend;
class TTree;
struct AliHLTTracksData;
class AliTPCclusterMI;

/**
 * @class AliHLTGlobalEsdConverterComponent
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
 * @ingroup alihlt_global_components
 */
class AliHLTGlobalEsdConverterComponent : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTGlobalEsdConverterComponent();
  /** destructor */
  virtual ~AliHLTGlobalEsdConverterComponent();

  // interface methods of base class
  const char* GetComponentID() {return "GlobalEsdConverter";};
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn() {return new AliHLTGlobalEsdConverterComponent;}
  void SetMakeFriends(Bool_t make) {fMakeFriends=make;}

 protected:
  // interface methods of base class
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData,
	       AliHLTComponentTriggerData& trigData);

  using AliHLTProcessor::DoEvent;

  /**
   * Process the input data blocks.
   * @param pTree    tree to be filled
   * @param pESD     ESD to be filled
   * @return neg. error code if failed
   */
  int ProcessBlocks(TTree* pTree, AliESDEvent* pESD, AliESDfriend *pESDfriend);

 // void FillBenchmarkHistos(Double_t *statistics, TString *names);
 private:
  /** copy constructor prohibited */
  AliHLTGlobalEsdConverterComponent(const AliHLTGlobalEsdConverterComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalEsdConverterComponent& operator=(const AliHLTGlobalEsdConverterComponent&);

  /**
   * (Re)Configure from the CDB
   * Loads the following objects:
   * - HLT/ConfigHLT/SolenoidBz
   */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /// write object to TTree or directly
  int fWriteTree; //!transient

  /// verbosity level
  int fVerbosity; //!transient

protected:

  static const Int_t fkNPartition = 36*6;           // number of patches in TPC

  /// the ESD
  AliESDEvent* fESD; //! transient value

  /// the ESD friend
  AliESDfriend* fESDfriend; //! transient value

  /// solenoid b field
  Double_t fSolenoidBz; //! transient
  Int_t fScaleDownTracks; //!

  Bool_t fMakeFriends; // flag to create friends
  AliTPCclusterMI   *fPartitionClusters[fkNPartition];  //! arrays of cluster data for each TPC partition
  Int_t              fNPartitionClusters[fkNPartition]; //! number of clusters for each TPC partition
 
  AliHLTComponentBenchmark fBenchmark; // benchmark

  ClassDef(AliHLTGlobalEsdConverterComponent, 0)
};
#endif
