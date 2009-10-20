// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTROOTSCHEMAEVOLUTIONCOMPONENT_H
#define ALIHLTROOTSCHEMAEVOLUTIONCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTRootSchemaEvolutionComponent.h
    @author Matthias Richter
    @date   2009-10-18
    @brief  Handler component for ROOT schema evolution of streamed objects
*/

#include "AliHLTProcessor.h"

class TObjArray;

/**
 * @class AliHLTRootSchemaEvolutionComponent
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ROOTFileWriter                                      <br>
 * Library: \b libAliHLTUtil.so						<br>
 * Input Data Types: ::kAliHLTAnyDataType				<br>
 * Output Data Types: none						<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -fxs<=off> <br>
 *      push streamer info to FXS, fetched by Shuttle and store in the entry
 *      HLT/Calib/StreamerInfo
 *      default off
 * \li -hltout<=[all,first,eor,off]> <br>
 *      push streamer info to output, the streamer info is stored in the
 *      events in all, the first, and/or the EOR.
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * Depending on the mode.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRootSchemaEvolutionComponent : public AliHLTProcessor
{
 public:
  /// standard constructor
  AliHLTRootSchemaEvolutionComponent();
  /// destructor
  virtual ~AliHLTRootSchemaEvolutionComponent();

  /// inherited from AliHLTComponent: return id of the component.
  virtual const char* GetComponentID() {return "ROOTSchemaEvolutionComponent";};
  /// inherited from AliHLTComponent: input data types
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList&);
  /// inherited from AliHLTComponent: output data types
  virtual AliHLTComponentDataType GetOutputDataType();
  /// inherited from AliHLTComponent: output data size
  virtual void GetOutputDataSize(unsigned long&, double&);

  /// inherited from AliHLTComponent: spawn function, create an instance.
  virtual AliHLTComponent* Spawn() {return new AliHLTRootSchemaEvolutionComponent;}

  enum {
    /// push streamer info to the HLTOUT for the first event
    kHLTOUTatFirstEvent   = 0x1,
    /// push streamer info to the HLTOUT for all events
    kHLTOUTatAllEvents    = 0x2,
    /// push streamer info to the HLTOUT at EOR, this has no relevance
    /// for reconstruction as it is too late and just in one raw file,
    /// but it allows archival at the end of the run
    kHLTOUTatEOR          = 0x4,
    /// push streamer info to FXS
    kFXS                  = 0x100,
  };

  /// Update the array of known streamer infos from a list of infos
  /// Checks whether the provided infos are already there in the present version
  /// and adds if it is a new info. 
  int UpdateStreamerInfos(const TList* list, TObjArray* infos) const;

 protected:
  /// inherited from AliHLTComponent: custom initialization
  int DoInit( int argc, const char** argv );
  /// inherited from AliHLTComponent: cleanup
  int DoDeinit();

  /// inherited from AliHLTProcessor processing
  virtual int DoEvent( const AliHLTComponentEventData& evtData,
		       AliHLTComponentTriggerData& trigData );
  
  using AliHLTProcessor::DoEvent;

  /**
   * Inherited from AliHLTComponent
   * Scan one argument and adjacent parameters.
   * @return number of scanned parameters, neg. error code if failed
   */
  virtual int ScanConfigurationArgument(int argc, const char** argv);

  void SetBits(AliHLTUInt32_t b) {fFlags|=b;}
  void ClearBits(AliHLTUInt32_t b) {fFlags&=~b;}
  bool TestBits(AliHLTUInt32_t b) {return (fFlags&b) != 0;}

private:
  /** copy constructor prohibited */
  AliHLTRootSchemaEvolutionComponent(const AliHLTRootSchemaEvolutionComponent&);
  /** assignment operator prohibited */
  AliHLTRootSchemaEvolutionComponent& operator=(const AliHLTRootSchemaEvolutionComponent&);

  AliHLTUInt32_t fFlags; //! property flags

  TObjArray* fpStreamerInfos; //! array of streamer infos

  AliHLTUInt32_t fFXSPrescaler; //! prescalar for the publishing to FXS

  static const char* fgkConfigurationObject; //! configuration object

  ClassDef(AliHLTRootSchemaEvolutionComponent, 1) // ROOT schema evolution component
};
#endif
