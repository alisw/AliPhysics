//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTBENCHEXTERNALTRACKCOMPONENT_H
#define ALIHLTBENCHEXTERNALTRACKCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTBenchExternalTrackComponent.h
    @author Matthias Richter
    @date   2008-10-30
    @brief  Benchmark component for AliExternalTrackParam transportation.
*/

// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"

class AliExternalTrackParam;
class AliHLTExternalTrackParam;
class TClonesArray;
class TObjArray;
class TRandom;
class TList;

/**
 * @class AliHLTBenchExternalTrackComponent
 * A component publishing arrays of randomly filled AliExternalTrackParam
 * structures in different forms. Currently supported is TClonesArray,
 * TObjArray, and C Array.
 * 
 * The component can both publish and receive data blocks. If data blocks
 * are available on the input, the array is restored and forwarded in the
 * specified publishing mode with specified compression level.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b BenchmarkAliExternalTrackParam <br>
 * Library: \b libAliHLTBenchmark.so     <br>
 * Input Data Types: ::kAliHLTDataTypeTrack, ::kAliHLTDataTypeTObjArray <br>
 * Output Data Types: ::kAliHLTDataTypeTrack, ::kAliHLTDataTypeTObjArray <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * By default, publishing of the array is switched off. It has to be enabled
 * by one of the three options
 * \li -tobjarray <br>
 *      publish TObjArray
 * \li -tclonesarray <br>
 *      publish TClonesArray
 * \li -carray <br>
 *      publish as C-array
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -maxsize      <i> max array size  </i> <br>
 * \li -minsize      <i> min array size  </i> <br>
 *
 * \li -rangemodulo  <i> count  </i> <br>
 *      number of events after which the range is changed
 * \li -rangeoffset  <i> offset  </i> <br>
 *      offset is added to min and max range after the specified number of
 *      events. Idially, this offset is negative, causing decreasing range
 *      over run time.
 * \li -rangefactor  <i> factor  </i> <br>
 *      min and max range are multiplied by the factor (float) after the
 *      specified number of events.
 *
 * \li -verbosity      <i> level  </i> <br>
 *      different levels of verbosity: 0 (default) is silent
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has no default CDB entries.
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 *
 * By using the range alterator options, the range can be changed after a
 * certain number of events.
 *
 * @ingroup alihlt_benchmark_components
 */
class AliHLTBenchExternalTrackComponent : public AliHLTProcessor {
 public:
  /** default constructor */
  AliHLTBenchExternalTrackComponent();
  /** destructor */
  virtual ~AliHLTBenchExternalTrackComponent();

  // interface functions: property getters
  const char* GetComponentID();
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn();

  /**
   * Fill members of the track object with random float values.
   * @param track     reference to AliExternalTrackParam object
   * @param fillCov   number of elements of the covar members to be filled
   *                  max 15
   */
  static int FillRandom(AliExternalTrackParam* track, int fillCov=15);

  /**
   * Serialize an object or clones array of AliExternalTrackParam objects
   * to a buffer using the definition AliHLTExternalTrackParam.
   * @param pArray    the array
   * @param buffer    target buffer
   * @param size      buffer size
   * @return bytes filled to the buffer, neg. error code if failed
   */
  static int SerializeToStruct(const TObjArray* pArray, AliHLTUInt8_t* buffer, unsigned int size);

  /**
   * Translate AliHLTExternalTrackParam struct into objects.
   * The TObjArray can also be a TClonesArray and must have enough capacity.
   */
  static int ReadFromStruct(TObjArray* pTgtArray, AliHLTExternalTrackParam* pArray, unsigned int arraySize);

  /**
   * Calculate a crc checksum for the object array
   */
  static AliHLTUInt32_t CalcChecksum(const TObjArray* pArray);

  /**
   * Compare two arrays of AliExternalTrackParamElements
   */
  static bool Compare(const TObjArray* array1, const TObjArray* array2);

  /**
   * Find object in the registry
   */
  static TObject* FindObject(AliHLTUInt32_t id);

  enum {
    kPublishingOff = 0,
    ktobjarray = 1,
    ktclonesarray,
    kcarray,
  };

 protected:
  // interface functions: processing
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks );

  /**
   * Register object and create reference id
   * The registry is used for consistency checks in a single threaded system.
   * The dump component can access the original object and compare the two.
   */
  AliHLTUInt32_t Register(TObject* pObject);

  /**
   * Remove object from registry
   */
  int Unregister(TObject* pObject);
 
 private:
  /** copy constructor prohibited */
  AliHLTBenchExternalTrackComponent(const AliHLTBenchExternalTrackComponent&);
  /** assignment operator prohibited */
  AliHLTBenchExternalTrackComponent& operator=(const AliHLTBenchExternalTrackComponent&);

  /** verbosity */
  int fVerbosity; //!transient

  /** disable global object registry, no check pf received objects */
  bool fDisableRegistry; //! transient

  /** publishing mode: off, ktclonesarray, ktobjarray, kcarray */
  int fMode; //! transient

  /** maximum array size */
  int fMaxSize; //! transient

  /** minimum array size */
  int fMinSize; //! transient

  /** event count for the next automatic change of range */
  int fEventModulo;

  /** offset is added to range after fEventModulo events */
  int fRangeOffset;

  /** range is multiplied after fEventModulo events */
  float fRangeMultiplicator;

  /** the internal array */
  TClonesArray* fpTcArray; //! transient

  /** TObjArray this publishing mode was used */
  TObjArray* fpTObjArray; //! transient

  /** random number generator */
  TRandom* fpDice; //! transient

  /** registry of object arrays */
  static TList* fgpRegistry; //! transient

  ClassDef(AliHLTBenchExternalTrackComponent, 0);
};

#endif
