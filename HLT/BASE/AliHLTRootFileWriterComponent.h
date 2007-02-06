// @(#) $Id$

#ifndef ALIHLTROOTFILEWRITERCOMPONENT_H
#define ALIHLTROOTFILEWRITERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRootFileWriterComponent.h
    @author Matthias Richter
    @date   
    @brief  Base class for writer components to store data in a ROOT file

                                                                          */
#include "AliHLTFileWriter.h"
#include "TObject.h" 

class TFile;

/**
 * @class AliHLTRootFileWriterComponent
 * @see AliHLTFileWriter for parameters
 */
class AliHLTRootFileWriterComponent : public AliHLTFileWriter
{
 public:
  /** standard constructor */
  AliHLTRootFileWriterComponent();
  /** destructor */
  ~AliHLTRootFileWriterComponent();

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
  /**
   * Data processing method for the component.
   * The function can be overloaded by specific ROOT file writer
   * components.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData );

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
  int WriteObject(const AliHLTEventID_t eventID, TObject *pOb);

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
  AliHLTEventID_t fEventID;

  /** the name of the current file */
  TFile* fCurrentFile; //! transient value

  ClassDef(AliHLTRootFileWriterComponent, 0)
};
#endif
