// @(#) $Id$

#ifndef ALIHLTTPCESDWRITERCOMPONENT_H
#define ALIHLTTPCESDWRITERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCEsdWriterComponent.h
    @author Matthias Richter
    @date   
    @brief  Writer component to store tracks of the HLT TPC conformal
            mapping tracker in the AliESD format

                                                                          */
#include "AliHLTRootFileWriterComponent.h"

// forward declarations
class TTree;
class AliESD;
class AliHLTTPCTrackArray;

/**
 * @class AliHLTRootFileWriterComponent
 * @see AliHLTFileWriter for parameters
 */
class AliHLTTPCEsdWriterComponent : public AliHLTRootFileWriterComponent
{
 public:
  /** standard constructor */
  AliHLTTPCEsdWriterComponent();
  /** destructor */
  ~AliHLTTPCEsdWriterComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "TPCEsdWriter";};

  /**
   * Spawn function.
   * @return new class instance
   */
  AliHLTComponent* Spawn() {return new AliHLTTPCEsdWriterComponent;}

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
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  int ScanArgument(int argc, const char** argv);

 private:
  /**
   * Init the writer.
   * The DoInit function is not available for this child class. InitWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  int InitWriter();

  /**
   * Init the writer.
   * The DoDeinit function is not available for this child class. CloseWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  int CloseWriter();

  /**
   * Covert tracks to AliTPCtracks (AliKalmanTracks) and add them to ESD.
   * @param pTracks  array of tracks
   * @param pESD     pointer to ESD
   * @return neg. error code if failed
   */
  int Tracks2ESD(AliHLTTPCTrackArray* pTracks, AliESD* pESD);

  /** the ESD tree */
  TTree* fTree; //! transient value

  /** the ESD */
  AliESD* fESD; //! transient value

  ClassDef(AliHLTTPCEsdWriterComponent, 0)
};
#endif
