// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTDATAGENERATOR_H
#define ALIHLTDATAGENERATOR_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataGenerator.h
    @author Matthias Richter
    @date   
    @brief  An HLT file publishing (data source) component.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTDataSource.h"

/**
 * @class AliHLTDataGenerator
 * An HLT data source component to produce random data.
 *
 * Component ID: \b DataGenerator <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype 'CLUSTERS' 'TPC ' </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 * \li -minsize      <i> size </i> <br>
 *      the minimum size of the data to be produced
 * \li -maxsize      <i> size </i> <br>
 *      the maximum size of the data to be produced, default = minsize
 *
 * Optional arguments:<br>
 * \li -disisor <i> m </i> <br>
 *      a disisor to shrink the size after \em modulo events
 * \li -offset <i> m </i> <br>
 *      an offset to subtract from the size after \em modulo events
 * \li -modulo <i> n </i> <br>
 *      size manipulated by the disisor or subtractor after \em n events
 *
 * The component produces data blocks of random content and random size in the
 * range of [\em minsize , \em maxsize ]. The size arguments can contain 'k' or
 * 'M' to indicate kByte or MByte.
 * @ingroup alihlt_component
 */
class AliHLTDataGenerator : public AliHLTDataSource  {
 public:
  /** standard constructor */
  AliHLTDataGenerator();
  /** destructor */
  virtual ~AliHLTDataGenerator();

  const char* GetComponentID();
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(vector<AliHLTComponentDataType>& tgtList);
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Deinit method.
   */
  int DoDeinit();

  /**
   * Data processing method for the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   * @return
   */
  int GetEvent( const AliHLTComponentEventData& evtData,
		        AliHLTComponentTriggerData& trigData,
		        AliHLTUInt8_t* outputPtr, 
		        AliHLTUInt32_t& size,
		        vector<AliHLTComponentBlockData>& outputBlocks );

  using AliHLTDataSource::GetEvent;

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

 private:
  /** prohibit copy constructor */
  AliHLTDataGenerator(const AliHLTDataGenerator&);
  /** prohibit assignment operator */
  AliHLTDataGenerator& operator=(const AliHLTDataGenerator&);

  /** data type */
  AliHLTComponentDataType fDataType;                                //! transient

  /** specification */
  AliHLTUInt32_t fSpecification;                                    //! transient

  /** the original size size */
  AliHLTUInt32_t fSize;                                             //! transient

  /** the manipulated size */
  AliHLTUInt32_t fCurrSize;                                         //! transient

  /** range */
  AliHLTUInt32_t fRange;                                            //! transient

  /** divisor */
  AliHLTUInt32_t fDivisor;                                          //! transient

  /** subtractor */
  AliHLTUInt32_t fSubtractor;                                       //! transient

  /** modulo for size manipulation */
  AliHLTUInt32_t fModulo;                                           //! transient

  ClassDef(AliHLTDataGenerator, 0)
};
#endif
