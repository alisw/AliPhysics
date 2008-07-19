// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTROOTFILEWRITERCOMPONENT_H
#define ALIHLTROOTFILEWRITERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTRootFileWriterComponent.h
    @author Matthias Richter
    @date   
    @brief  Base class for writer components to store data in a ROOT file
*/

#include "AliHLTFileWriter.h"

class TFile;

/**
 * @class AliHLTRootFileWriterComponent
 * The RootFileWriter provides a stand alone component to write incoming
 * TObject like structures into a Root file. Furthermore it provides a
 * base class for customized writers.
 * By default, the \em -concatenate-blocks option of the AliHLTFileWriter
 * is set. If you want to accumulate all events in the same file set the
 * -concatenate-events option as well. In that case the \em -overwrite
 * option might be a good choice in order to avoid multiple keys for the
 * same object in the root file.
 *
 * All non-root object data blocks are just ignored.
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
 * See AliHLTFileWriter for full list of arguments.
 * \li -overwrite <br>
 *      write objects with the TObject::kOverwrite flag and avoid multiple
 *      keys in the file
 *
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
 * No data published (AliHLTDataSink).
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRootFileWriterComponent : public AliHLTFileWriter
{
 public:
  /** standard constructor */
  AliHLTRootFileWriterComponent();
  /** destructor */
  virtual ~AliHLTRootFileWriterComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  virtual const char* GetComponentID() {return "ROOTFileWriter";};

  /**
   * Spawn function.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn() {return new AliHLTRootFileWriterComponent;}

 protected:
  // interface functions
  int InitWriter();
  int CloseWriter();

  /**
   * Data processing method for the component.
   * The function can be overloaded by specific ROOT file writer
   * components. The RootFileWriter processes only TObject like data
   * structures of the input blocks and uses the
   * @ref alihltcomponent-high-level-interface. Despite of that it implements
   * the lox-level DumpEvent method in order to allow child classes to use the
   * low-level method.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData );
  
  using AliHLTFileWriter::DumpEvent;

  /**
   * Scan one argument and adjacent parameters.
   * \b IMPORTANT: if  overloaded by child class, call this function
   * as the default from the cutomized switch, e.g.
   * <pre>
   * </pre>
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual int ScanArgument(int argc, const char** argv);

  /**
   * Write ROOT object to current file.
   * @param eventID    ID of the current event
   * @param pOb        pointer to ROOT object
   * @return neg. error code if failed
   */
  int WriteObject(const AliHLTEventID_t eventID, const TObject *pOb);

  /**
   * Open a ROOT file.
   * The function calls @ref AliHLTFileWriter::BuildFileName in order to
   * create a file name and opens it as a root file.
   * @param eventID    ID of the current event
   * @param blockID    ID of the current block
   * @param option     option as specified in TFile
   * @return pointer to TFile object, the called has to clean-up the object after use.
   */
  TFile* OpenFile(const AliHLTEventID_t eventID, const int blockID=-1, const char* option="recreate");

  /** the event ID associated with the current file */
  AliHLTEventID_t fEventID; // see above

  /** the name of the current file */
  TFile* fCurrentFile; //! transient value

  /** options for the TObject::Write function */
  Int_t fOptions; //!transient

private:
  /** copy constructor prohibited */
  AliHLTRootFileWriterComponent(const AliHLTRootFileWriterComponent&);
  /** assignment operator prohibited */
  AliHLTRootFileWriterComponent& operator=(const AliHLTRootFileWriterComponent&);

  ClassDef(AliHLTRootFileWriterComponent, 1) // ROOT file writer component
};
#endif
