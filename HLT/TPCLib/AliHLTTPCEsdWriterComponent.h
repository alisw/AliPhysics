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
class AliESDEvent;
class AliHLTTPCTrackArray;

/**
 * @class AliHLTTPCEsdWriterComponent
 * This class translates incoming track segments structures from the TPC
 * conformal mapping tracker (datatype TRAKSEGS/TPC) or tracks in global 
 * coordinates from the AliHLTTPCGlobalMergerComponent (TRACKS/TPC) into
 * the ESD fromat and writes it to a ROOT file. In case of TRAKSEGS, the
 * component can only process data block containing data of one slice, but
 * it can read an unlimeted number of data blocks.
 *
 * componentid: TPCEsdWriter <br>
 * componentlibrary libAliHLTTPC.so <br>
 * Arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -datafile     <i> filename   </i> <br>
 *      file name base
 * \li -directory    <i> directory  </i> <br>
 *      target directory
 *
 * <pre>
 * Example usage (HLT configuration file):
 *         \<Proc ID="EsdWriter" type="prc">
 *             \<Cmd>AliRootWrapperSubscriber -eventmodulo 1
 *                 -componentid TPCEsdWriter
 *                 -componentlibrary libAliHLTTPC.so
 *                 -componentargs "-datafile AliESDs.root"
 *            \</Cmd>
 *
 *            \<Parent>TR0-SC\</Parent>
 *            \<Node>master\</Node>
 *            \<Shm blocksize="1k" blockcount="1" type="sysv"/>
 *        \</Proc>
 * </pre>
 *
 * @see AliHLTFileWriter and AliHLTRootFileWriterComponent for more parameters
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
  /** copy constructor prohibited */
  AliHLTTPCEsdWriterComponent(const AliHLTTPCEsdWriterComponent&);
  /** assignment operator prohibited */
  AliHLTTPCEsdWriterComponent& operator=(const AliHLTTPCEsdWriterComponent&);

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
  int Tracks2ESD(AliHLTTPCTrackArray* pTracks, AliESDEvent* pESD);

  /** the ESD tree */
  TTree* fTree; //! transient value

  /** the ESD */
  AliESDEvent* fESD; //! transient value

  ClassDef(AliHLTTPCEsdWriterComponent, 1)
};
#endif
