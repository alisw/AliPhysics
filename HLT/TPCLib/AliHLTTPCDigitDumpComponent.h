//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTTPCDIGITDUMPCOMPONENT_H
#define ALIHLTTPCDIGITDUMPCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDigitDumpComponent.h
    @author Matthias Richter
    @date   
    @brief  Special file writer converting TPC digit input to ASCII.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTFileWriter.h"

/**
 * @class AliHLTTPCDigitDumpComponent
 * A converter for digit data of the TPC input to ASCII output.
 * Data blocks of type @ref kAliHLTDataTypeDDLRaw and origin 'TPC ' is docoded
 * written in readable ASCII format to a file.
 * 
 * Component ID: \b TPCDigitDump <br>
 * Library: \b libAliHLTTPC
 *
 * See AliHLTFileWriter for arguments, further specific options
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -rawreadermode  <i> mode   </i> <br>
 *      data mode of the AliHLTTPCDigitReaderRaw</i> 
 * \li -digitreader    <i> reader   </i> <br>
 *      type of the digit reader: <i>unpacked, packed, raw</i> 
 *
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

  enum {
    kDigitReaderInvalid,
    kDigitReaderUnpacked,
    kDigitReaderPacked,
    kDigitReaderRaw
  };

  /** the mode for the DigitReaderRaw */
  unsigned fRawreaderMode; //!transient

  /** the digit reader to use */
  short fDigitReaderType; //!transient

  ClassDef(AliHLTTPCDigitDumpComponent, 0);
};

#endif
