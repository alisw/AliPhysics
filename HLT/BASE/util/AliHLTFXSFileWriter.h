// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTFXSFILEWRITER_H
#define ALIHLTFXSFILEWRITER_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTFXSFileWriter.h
/// @author Timo Breitner
/// @date   
/// @brief  An HLT FXS file dump (data sink) component.
///

#include "AliHLTDataSink.h"
#include <TString.h>

class AliHLTFXSFileWriter : public AliHLTDataSink  {
 public:
  /** standard constructor */
  AliHLTFXSFileWriter();
  /** destructor */
  virtual ~AliHLTFXSFileWriter();

  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( AliHLTComponentDataTypeList& list);
  virtual AliHLTComponent* Spawn();

  static const Int_t gkCalibObjectMaxVersion;

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Data processing method for the component.
   * The function can be overloaded by other file writer components.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& trigData );

  /**
   * Get the target directory
   */
  TString GetDirectory() {return fDirectory;}


 private:
  /** copy constructor prohibited */
  AliHLTFXSFileWriter(const AliHLTFXSFileWriter&);
  /** assignment operator prohibited */
  AliHLTFXSFileWriter& operator=(const AliHLTFXSFileWriter&);

  /** target directory */
  TString    fDirectory;                                           // see above
  /** base name of the event sub directories */
  
  ClassDef(AliHLTFXSFileWriter, 0)
};
#endif
