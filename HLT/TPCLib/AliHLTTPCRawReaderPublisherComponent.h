//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCRAWREADERPUBLISHERCOMPONENT_H
#define ALIHLTTPCRAWREADERPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCRawReaderPublisherComponent.h
/// @author Matthias Richter
/// @date   2011-11-18
/// @brief  Specific publisher for TPC raw data from the AliRawReader
///         

#include "AliHLTRawReaderPublisherComponent.h"
#include <map>

/**
 * @class AliHLTTPCRawReaderPublisherComponent
 * This component uses the functionality of AliHLTRawReaderPublisherComponent
 * and overloads IsSelected and GetSpecificationFromEquipmentId. Blocks are
 * only generated if the corresponding partition is missing in HLTOUT.
 *
 * It is used in an emulation chain which produces all compressed cluster
 * blocks which are missing in HLTOUT. If TPC reconstruction requires HLT
 * clusters, the emulator is automatically executed and the compressed
 * data produced if raw data is available.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCRawReaderPublisher      <br>
 * Library: \b libAliHLTTPC.so     <br>
 * Input Data Types:  <br>
 * Output Data Types: <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->

 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Output size:</h2>
 *
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCRawReaderPublisherComponent : public AliHLTRawReaderPublisherComponent {
public:
  /// standard constructor
  AliHLTTPCRawReaderPublisherComponent();
  /// destructor
  ~AliHLTTPCRawReaderPublisherComponent();

  /// inherited from AliHLTComponent: id of the component
  virtual const char* GetComponentID();

  /// inherited from AliHLTComponent: spawn function.
  virtual AliHLTComponent* Spawn();

protected:
  /// inherited from AliHLTDataSource: get one event
  int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks );

  /// inherited from AliHLTComponent: initialize
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent: cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: argument scan
  int ScanConfigurationArgument(int argc, const char** argv);

  /// check the HLTOUT for availability of compressed data blocks
  int InitMapFromHLTOUT(std::map<AliHLTUInt32_t, bool>& hltoutmap);

  /// inherited from AliHLTRawReaderPublisherComponent: get specification
  virtual int GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t &specification) const;

  /// inherited from AliHLTRawReaderPublisherComponent: check if a block is selected or not
  virtual bool IsSelected(int equipmentId) const;

private:
  AliHLTTPCRawReaderPublisherComponent(const AliHLTTPCRawReaderPublisherComponent&);
  AliHLTTPCRawReaderPublisherComponent& operator=(const AliHLTTPCRawReaderPublisherComponent&);

  bool* fArraySelected; //! transient

  ClassDef(AliHLTTPCRawReaderPublisherComponent, 0)
};

#endif //ALIHLTTPCRAWREADERPUBLISHERCOMPONENT_H
