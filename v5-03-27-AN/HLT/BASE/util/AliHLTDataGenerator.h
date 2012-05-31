// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDATAGENERATOR_H
#define ALIHLTDATAGENERATOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTDataGenerator.h
    @author Matthias Richter
    @date   
    @brief  An HLT file publishing (data source) component.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTProcessor.h"

class TRandom;

/**
 * @class AliHLTDataGenerator
 * An HLT data source component to produce random data.
 *
 * The component produces fake data blocks according to the size range
 * or the total input size. For the former, the initial size can be set
 * by the \em -size, \em -minsize arguments, and the range respectivly.
 * Tn this mode, the size can be decremented (\em -decrement) after n
 * events (\em -modulo).
 *
 * When producing fake blocks with respect to the total input data,
 * \em -offset and \em -multiplier can be used to set ratio and offset.
 *
 * Component ID: \b DataGenerator <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments:<br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype 'CLUSTERS' 'TPC ' </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 * \li -size      <i> size </i> <br>
 *      initial size of the data to be produced
 * \li -range      <i> n </i> <br>
 *      range of the data to be produced [size,size+range]
 * \li -minsize      <i> size </i> <br>
 *      the minimum size of the data to be produced
 * \li -maxsize      <i> size </i> <br>
 *      the maximum size of the data to be produced, default = minsize
 *
 * \li -decrement <i> m </i> <br>
 *      subtract after \em modulo events
 * \li -modulo <i> n </i> <br>
 *      size manipulated after \em n events
 *
 * \li -offset <i> m </i> <br>
 *      output size = offset + factor * (total size of input blocks)
 * \li -multiplier <i> factor </i> <br>
 *      output size = offset + factor * (total size of input blocks)
 *
 * The component produces data blocks of random content and random size in the
 * range of [\em minsize , \em maxsize ]. The size arguments can contain 'k' or
 * 'M' to indicate kByte or MByte.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTDataGenerator : public AliHLTProcessor  {
 public:
  /** standard constructor */
  AliHLTDataGenerator();
  /** destructor */
  virtual ~AliHLTDataGenerator();

  const char* GetComponentID();
  void GetInputDataTypes( AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(vector<AliHLTComponentDataType>& tgtList);
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

 protected:
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       AliHLTComponentBlockDataList& outputBlocks );

  using AliHLTProcessor::DoEvent;

  /**
   * Scan one argument and adjacent parameters.
   * Can be overloaded by child classes in order to add additional arguments
   * beyond the standard arguments of the file publisher. The method is called
   * whenever a non-standard argument is recognized.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual int ScanArgument(int argc, const char** argv);

 protected:
  /**
   * Scan a size argument.
   * The argument is expected to be an integer, which can be suffixed by 'k'
   * or 'M' in order to indicate the base, kByte or MByte.
   * @param size      target to store the size
   * @param arg       the argument to scan
   */
  int ScanSizeArgument(AliHLTUInt32_t &size, const char* arg);

  /**
   * Scan a float argument.
   * @param value     target to store the size
   * @param arg       the argument to scan
   */
  int ScanFloatArgument(float &value, const char* arg);

 private:
  /** prohibit copy constructor */
  AliHLTDataGenerator(const AliHLTDataGenerator&);
  /** prohibit assignment operator */
  AliHLTDataGenerator& operator=(const AliHLTDataGenerator&);

  /** data type */
  AliHLTComponentDataType fDataType;                                //! transient

  /** specification */
  AliHLTUInt32_t fSpecification;                                    //! transient

  // mode 1: just fake data independent of the input data

  /** the original size size */
  AliHLTUInt32_t fSize;                                             //! transient
  /** range: [size, size+range] */
  AliHLTUInt32_t fRange;                                                     //! transient
  /** the manipulated size */
  AliHLTUInt32_t fCurrSize;                                         //! transient
  /** divisor: each modulo event ignoring the input data size) */
  AliHLTUInt32_t fDivisor;                                          //! transient
  /** decrement: each modulo event ignoring the input data size */
  AliHLTUInt32_t fDecrement;                                        //! transient
  /** modulo for size manipulation */
  AliHLTUInt32_t fModulo;                                           //! transient

  // mode 2: generate data depending on input data size

  /** offset (generation of data from input data size) */
  AliHLTUInt32_t fOffset;                                           //! transient
  /** multiplier (generation of data from input data size) */
  float fMultiplier;                                                //! transient

  /** random number generator */
  TRandom* fpDice;                                                  //! transient

  ClassDef(AliHLTDataGenerator, 0)
};
#endif
