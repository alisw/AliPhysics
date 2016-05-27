// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTROOTSCHEMAEVOLUTIONCOMPONENT_H
#define ALIHLTROOTSCHEMAEVOLUTIONCOMPONENT_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTRootSchemaEvolutionComponent.h
/// @author Matthias Richter
/// @date   2009-10-18
/// @brief  Handler component for ROOT schema evolution of streamed objects
///

#include "AliHLTCalibrationProcessor.h"
#include "TString.h"
#include <vector>

class TObjArray;
class TObject;
class TStopwatch;
class AliHLTMessage;

using std::vector;

/**
 * @class AliHLTRootSchemaEvolutionComponent
 * Collects streamer info for all input objects and produces the corresponding
 * calibration object for reconstruction of HLT. The component runs with a
 * configurable rate constraint and skips the processing of known data blocks
 * for the sake of performance. New data blocks are always processed and added
 * to the list.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ROOTSchemaEvolutionComponent                        <br>
 * Library: \b libAliHLTUtil.so						<br>
 * Input Data Types: ::kAliHLTAnyDataType				<br>
 * Output Data Types: none						<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -fxs<=[n,off]> <br>
 *      push streamer info to FXS, fetched by Shuttle and store in the entry
 *      HLT/Calib/StreamerInfo                                          <br>
 *      if a scalar greather then 0 is specified the calibration object is
 *      pushed during the event processing with the specified scale down<br>
 *      always pushed on EOR, default on
 * \li -hltout<=[all,first,eor,off]> <br>
 *      push streamer info to output, the streamer info is stored in the
 *      events in all, the first, and/or the EOR.
 * \li -file=filename <br>
 *      write to file at EOR
 * \li -rate=hz <br>
 *      required processing rate in Hz, default 2000Hz
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * TODO: update performance requirements for unpacking ESDs and creating the
 * streamer info
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * Depending on the mode.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRootSchemaEvolutionComponent : public AliHLTCalibrationProcessor
{
 public:
  /// standard constructor
  AliHLTRootSchemaEvolutionComponent();
  /// destructor
  virtual ~AliHLTRootSchemaEvolutionComponent();

  /// inherited from AliHLTComponent: return id of the component.
  virtual const char* GetComponentID() {return "ROOTSchemaEvolutionComponent";};
  /// inherited from AliHLTComponent: input data types
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
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
  int UpdateStreamerInfos(const TObjArray* list, TObjArray* infos) const;

  class AliHLTDataBlockItem
  {
  public:
    AliHLTDataBlockItem(AliHLTComponentDataType dt=kAliHLTVoidDataType,
			AliHLTUInt32_t spec=kAliHLTVoidDataSpec);
    ~AliHLTDataBlockItem();

    /// extract data block to root object, and update performance parameters
    /// object needs to be deleted externally
    TObject* Extract(const AliHLTComponentBlockData* bd);

    /// stream object and update performance parameters
    int Stream(const TObject* obj, AliHLTMessage& msg);

    bool IsObject() const {return fIsObject;}
    bool operator==(const AliHLTDataBlockItem& i) const {return fDt==i.fDt && fSpecification==i.fSpecification;}
    bool operator==(AliHLTComponentDataType dt) const {return fDt==dt;}
    bool operator==(AliHLTUInt32_t spec) const {return fSpecification==spec;}
    operator const AliHLTComponentDataType&() const {return fDt;}
    AliHLTUInt32_t GetSpecification() const {return fSpecification;}
    
    /// average extraction time in usec
    AliHLTUInt32_t GetExtractionTime() const {return fNofExtractions>0?fExtractionTimeUsec/fNofExtractions:0;}
    /// average streaming time in usec
    AliHLTUInt32_t GetStreamingTime() const {return fNofStreamings>0?fStreamingTimeUsec/fNofStreamings:0;}
    /// average total time in usec
    AliHLTUInt32_t GetTotalTime() const {return GetExtractionTime() + GetStreamingTime();}

    /// print status
    void Print(const char* option) const;

  private:
    /// data type of the block
    AliHLTComponentDataType fDt; //! transient
    /// specification of the block
    AliHLTUInt32_t fSpecification; //! transient
    /// flag for TObject
    bool fIsObject; //! transient

    /// number of extractions
    AliHLTUInt32_t fNofExtractions; //! transient
    /// object extraction time in usec
    AliHLTUInt32_t fExtractionTimeUsec; //! transient
    /// timestamp of last extraction in usec
    AliHLTUInt32_t fLastExtraction; //! transient
    /// number of streamings
    AliHLTUInt32_t fNofStreamings; //! transient
    /// object streaming time in usec
    AliHLTUInt32_t fStreamingTimeUsec; //! transient
    /// timestamp of last streaming in usec
    AliHLTUInt32_t fLastStreaming; // !transient
  };

  /// find item in the list
  AliHLTDataBlockItem* FindItem(AliHLTComponentDataType dt,
				AliHLTUInt32_t spec);

 protected:
  /// inherited from AliHLTCalibrationProcessor: custom initialization
  int InitCalibration();
  /// inherited from AliHLTCalibrationProcessor: custom argument scan
  /// the AliHLTCalibrationProcessor so far does not use the base class
  /// methods for argument scan.
  int ScanArgument( int argc, const char** argv ) {
    int result=ScanConfigurationArgument(argc, argv); return result>0?result-1:result;
  }
  /// inherited from AliHLTCalibrationProcessor: cleanup
  int DeinitCalibration();

  /// inherited from AliHLTCalibrationProcessor processing
  virtual int ProcessCalibration( const AliHLTComponentEventData& evtData,
				  AliHLTComponentTriggerData& trigData );
  
  using AliHLTCalibrationProcessor::ProcessCalibration;

  /// inherited from AliHLTCalibrationProcessor processing
  int ShipDataToFXS( const AliHLTComponentEventData& evtData,
		     AliHLTComponentTriggerData& trigData);

  using AliHLTCalibrationProcessor::ShipDataToFXS;

  /**
   * Inherited from AliHLTComponent
   * Scan one argument and adjacent parameters.
   * @return number of scanned parameters, neg. error code if failed
   */
  virtual int ScanConfigurationArgument(int argc, const char** argv);

  void SetBits(AliHLTUInt32_t b) {fPropertyFlags|=b;}
  void ClearBits(AliHLTUInt32_t b) {fPropertyFlags&=~b;}
  bool TestBits(AliHLTUInt32_t b) const {return (fPropertyFlags&b) != 0;}
  int WriteToFile(const char* filename, const TObjArray* infos) const;

private:
  /** copy constructor prohibited */
  AliHLTRootSchemaEvolutionComponent(const AliHLTRootSchemaEvolutionComponent&);
  /** assignment operator prohibited */
  AliHLTRootSchemaEvolutionComponent& operator=(const AliHLTRootSchemaEvolutionComponent&);

  vector<AliHLTDataBlockItem> fList; //! list of block properties

  AliHLTUInt32_t fPropertyFlags; //! property flags

  TObjArray* fpStreamerInfos; //! array of streamer infos
  TStopwatch* fpEventTimer; //! stopwatch for event processing
  TStopwatch* fpCycleTimer; //! stopwatch for event cycle
  AliHLTUInt32_t fMaxEventTime; //! required maximum processing time in usec

  AliHLTUInt32_t fFXSPrescaler; //! prescalar for the publishing to FXS

  TString fFileName; //! file name for dump at EOR

  static const char* fgkConfigurationObject; //! configuration object
  static const AliHLTUInt32_t fgkTimeScale; //! timescale base

  ClassDef(AliHLTRootSchemaEvolutionComponent, 0) // ROOT schema evolution component
};
#endif
