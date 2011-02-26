//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITDUMPCOMPONENT_H
#define ALIHLTTPCDIGITDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCDigitDumpComponent.h
    @author Matthias Richter
    @date   
    @brief  Special file writer converting TPC digit input to ASCII.
*/

#include "AliHLTFileWriter.h"

class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCDigitDumpComponent
 * A converter for digit data of the TPC input to ASCII output.
 * Data blocks of type ::kAliHLTDataTypeDDLRaw and origin 'TPC ' is docoded
 * written in readable ASCII format to a file.
 * 
 * The component supports different types of readers in order to
 * choose different data formats (raw/digits) and reading modes.
 *
 * Component ID: \b TPCDigitDump <br>
 * Library: \b libAliHLTTPC
 *
 * See AliHLTFileWriter for arguments, further specific options
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -digitreader    <i> reader   </i> <br>
 *      type of the digit reader: <i>unpacked, packed, raw, decoder</i>
 *      default 'decoder' 
 * \li -rcutrailersize    <i> size   </i> <br>
 *      size of the RCU trailer in 32bit words (default 2), if digitreader
 *      'decoder' is used, the trailer size is determined automatically
 * \li -unsorted <br>
 *      unsorted mode of digit readers (default mode)
 * \li -sorted <br>
 *      sorted mode of digit readers (default is unsorted)
 * \li -bulk <br>
 *      bulk read mode: NextChannel/Bunch
 * \li -stream <br>
 *      stream read mode: Next
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCDigitDumpComponent : public AliHLTFileWriter {
 public:
  /** default constructor */
  AliHLTTPCDigitDumpComponent();
  /** destructor */
  virtual ~AliHLTTPCDigitDumpComponent();

  // interface functions: property getters
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  virtual AliHLTComponent* Spawn();

 protected:
  // interface functions: processing
  int InitWriter();
  int CloseWriter();
  int DumpEvent( const AliHLTComponentEventData& evtData,
  		 const AliHLTComponentBlockData* blocks, 
  		 AliHLTComponentTriggerData& trigData );
  using AliHLTDataSink::DumpEvent;

  int ScanArgument(int argc, const char** argv);

 private:
  /** copy constructor prohibited */
  AliHLTTPCDigitDumpComponent(const AliHLTTPCDigitDumpComponent&);
  /** assignment operator prohibited */
  AliHLTTPCDigitDumpComponent& operator=(const AliHLTTPCDigitDumpComponent&);

  /**
   * Print slice/partition/row/pad header if changed.
   */
  int PrintHeaders(int slice, int &iPrintedSlice,
		   int part, int &iPrintedPart,
		   AliHLTTPCDigitReader* pReader,
		   int &iPrintedRow, int &iPrintedPad,
		   ofstream &dump) const;

  enum {
    kDigitReaderInvalid,
    kDigitReaderUnpacked,
    kDigitReader32Bit
  };

  /** the digit reader to use */
  short fDigitReaderType; //!transient

  /** size of the RCU trailer in 32bit words */
  short fRcuTrailerSize; //! transient

  /** unsorted/sorted mode of digit readers */
  bool fUnsorted; //!transient

  /** bulk read mode */
  bool fbBulkMode; //!transient

  /** the digit reader */
  AliHLTTPCDigitReader* fpReader; //!transient

  /** flag if 32 bit format is used */
  Bool_t f32BitFormat;

  ClassDef(AliHLTTPCDigitDumpComponent, 0);
};

#endif
