// @(#) $Id$

#ifndef ALIHLTFILEWRITER_H
#define ALIHLTFILEWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTFileWriter.h
    @author Matthias Richter
    @date   
    @brief  An HLT file dump (data sink) component.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTDataSink.h"
#include <TList.h>

/**
 * @class AliHLTFileWriter
 * An HLT data sink component which writes data to file(s)
 *
 * @ingroup alihlt_component
 */
class AliHLTFileWriter : public AliHLTDataSink  {
 public:
  /** standard constructor */
  AliHLTFileWriter();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTFileWriter(const AliHLTFileWriter&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTFileWriter& operator=(const AliHLTFileWriter&);
  /** destructor */
  virtual ~AliHLTFileWriter();

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
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
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  int DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData );

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

 private:
  /**
   * Build file name from eventID data type and the specified directory and basename.
   * @param eventID [in]   the ID of the event
   * @param blockID [in]   the ID of the current block
   * @param dataType [in]  the data type of the data block
   * @param filename [out] string to receive the file name
   */
  int BuildFileName(const AliHLTEventID_t eventID, const int blockID, const AliHLTComponentDataType& dataType, TString& filename);

  /** the basename of the output file */
  TString    fBaseName;
  /** target directory */
  TString    fDirectory;
  /** enumeration format string */
  TString    fEnumeration;
  /**
   * flag to indicate whether to write each incoming block to separate files
   * or all blocks of one event to one file.
   */
  Int_t      fbSeparate;

  ClassDef(AliHLTFileWriter, 0)
};
#endif
