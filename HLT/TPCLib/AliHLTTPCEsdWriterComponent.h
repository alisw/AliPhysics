// @(#) $Id$

#ifndef ALIHLTTPCESDWRITERCOMPONENT_H
#define ALIHLTTPCESDWRITERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCEsdWriterComponent.h
    @author Matthias Richter
    @date   
    @brief  Writer component to store tracks of the HLT TPC conformal
            mapping tracker in the AliESD format
*/

#include "AliHLTRootFileWriterComponent.h"
#include "AliHLTProcessor.h"

// forward declarations
class TTree;
class AliESDEvent;
class AliHLTTPCTrackArray;
#include "AliHLTTPCClusterFinder.h"

/**
 * @class AliHLTTPCEsdWriterComponent
 * Translation of internal TPC track data to the ESD format.
 * This class translates incoming track segments structures from the TPC
 * conformal mapping tracker (datatype TRAKSEGS/TPC) or tracks in global 
 * coordinates from the AliHLTTPCGlobalMergerComponent (TRACKS/TPC) into
 * the ESD format.
 *
 * The \em TPCEsdWriter writes it to a ROOT file, the \em TPCEsdConverter
 * to a TTree and sends the whole object to the output stream with data
 * type ::kAliHLTDataTypeESDTree and origin TPC.
 *
 * In case of TRAKSEGS, the component can only process data block containing
 * data of one slice, but it can read an unlimeted number of data blocks.
 *
 * componentid: \b TPCEsdWriter <br>
 * componentid: \b TPCEsdConverter <br>
 * componentlibrary: \b libAliHLTTPC.so <br>
 * Arguments TPCEsdWriter: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -datafile     <i> filename   </i> <br>
 *      file name base
 * \li -directory    <i> directory  </i> <br>
 *      target directory
 *
 * Arguments TPCEsdConverter: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -notree                                                          <br>
 *      write ESD directly to output (::kAliHLTDataTypeESDObject)
 *      this has been made the default behavior in Sep 2008.
 * \li -tree                                                            <br>
 *      write ESD directly to TTree and to output (::kAliHLTDataTypeESDTree)
 * \li -solenoidBz                                                      <br>
 *      magnetic field like -solenoidBz 5.0
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
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCEsdWriterComponent : public AliHLTLogging
{
 public:
  /** standard constructor */
  AliHLTTPCEsdWriterComponent();
  /** destructor */
  ~AliHLTTPCEsdWriterComponent();

  /**
   * class AliHLTTPCEsdWriterComponent::AliWriter
   * The writer component of the AliHLTTPCEsdWriterComponent.
   */
  class AliWriter : public AliHLTRootFileWriterComponent
  {
  public:
  /** standard constructor */
  AliWriter();
  /** destructor */
  virtual ~AliWriter();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "TPCEsdWriter";};

  void GetInputDataTypes(AliHLTComponentDataTypeList& list);

  /**
   * Spawn function.
   * @return new class instance
   */
  AliHLTComponent* Spawn() {return new AliHLTTPCEsdWriterComponent::AliWriter;}

 protected:
  /**
   * Data processing method for the component.
   * The function can be overloaded by specific ROOT file writer
   * components.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );

  using AliHLTRootFileWriterComponent::DumpEvent;

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
  AliWriter(const AliWriter&);
  /** assignment operator prohibited */
  AliWriter& operator=(const AliWriter&);

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

  /** the ESD tree */
  TTree* fTree; //! transient value

  /** the ESD */
  AliESDEvent* fESD; //! transient value

  /** pointer to the basic ESD conversion methods */
  AliHLTTPCEsdWriterComponent* fBase; //! transient value
  };

  /**
   * class AliHLTTPCEsdWriterComponent::AliConverter
   * The converter component of the AliHLTTPCEsdWriterComponent.
   * 
   */
  class AliConverter : public AliHLTProcessor
  {
  public:
    /** standard constructor */
    AliConverter();
    /** destructor */
    virtual ~AliConverter();

    // interface methods of base class
    const char* GetComponentID() {return "TPCEsdConverter";};
    void GetInputDataTypes(AliHLTComponentDataTypeList& list);
    AliHLTComponentDataType GetOutputDataType();
    void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
    AliHLTComponent* Spawn() {return new AliHLTTPCEsdWriterComponent::AliConverter;}

  protected:
    // interface methods of base class
    int DoInit(int argc, const char** argv);
    int DoDeinit();
    int DoEvent(const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks );

    using AliHLTProcessor::DoEvent;

  private:
    /** copy constructor prohibited */
    AliConverter(const AliConverter&);
    /** assignment operator prohibited */
    AliConverter& operator=(const AliConverter&);

    /** the ESD */
    AliESDEvent* fESD; //! transient value

    /** pointer to the basic ESD conversion methods */
    AliHLTTPCEsdWriterComponent* fBase; //! transient value

    /** write object to TTree or directly */
    int fWriteTree; //!transient
  };

 protected:
  /**
   * Process the input data blocks.
   * @param pTree    tree to be filled
   * @param pESD     ESD to be filled
   * @param blocks   data block descriptor array
   * @param nBlocks  size of the array
   * @param pMinSlice [OUT] receives the minimum slice no
   * @param pMaxSlice [OUT] receives the maximum slice no
   * @return neg. error code if failed
   */
  int ProcessBlocks(TTree* pTree, AliESDEvent* pESD, const AliHLTComponentBlockData* blocks,
		    int nBlocks, int* pMinSlice=NULL, int* pMaxSlice=NULL);

  /**
   * Covert tracks to AliTPCtracks (AliKalmanTracks) and add them to ESD.
   * @param pTracks  array of tracks
   * @param pESD     pointer to ESD
   * @return neg. error code if failed
   */
  int Tracks2ESD(AliHLTTPCTrackArray* pTracks, AliESDEvent* pESD);

 private:
  /** copy constructor prohibited */
  AliHLTTPCEsdWriterComponent(const AliHLTTPCEsdWriterComponent&);
  /** assignment operator prohibited */
  AliHLTTPCEsdWriterComponent& operator=(const AliHLTTPCEsdWriterComponent&);

  /**
   * (Re)Configure from the CDB
   * Loads the following objects:
   * - HLT/ConfigHLT/SolenoidBz
   */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */  
  int Configure(const char* arguments);

  /** solenoid b field */
  Double_t fSolenoidBz; //! transient

  /** array of pointers to cluster MC labels **/

  AliHLTTPCClusterFinder::ClusterMCInfo *fClusterLabels[36*6]; //! cluster MC labels for each TPC patch
  Int_t fNClusterLabels[36*6]; //! Number of MC labels, for check of consistensy
  
  /** flag for calculating MC labels for tracks **/

  Bool_t fDoMCLabels;

  ClassDef(AliHLTTPCEsdWriterComponent, 3)
};
#endif
