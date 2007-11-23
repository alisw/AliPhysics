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
 * Data is written to file.
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
  virtual int InitWriter();
  virtual int CloseWriter();
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData );
  virtual int ScanArgument(int argc, const char** argv);

 
 private:
  /** copy constructor prohibited */
  AliHLTTPCDigitDumpComponent(const AliHLTTPCDigitDumpComponent&);
  /** assignment operator prohibited */
  AliHLTTPCDigitDumpComponent& operator=(const AliHLTTPCDigitDumpComponent&);

  ClassDef(AliHLTTPCDigitDumpComponent, 0);
};

#endif
