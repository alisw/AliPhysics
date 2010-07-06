//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTReconstructor.h
//  @author Matthias Richter
//  @date   
//  @brief  Binding class for HLT reconstruction in AliRoot
//          Implements bot the interface to run HLT chains embedded into
//          AliReconstruction and the unpacking and treatment of HLTOUT

#include "AliReconstructor.h"

class AliHLTSystem;
class AliRawReader;
class AliESDEvent;
class AliHLTOUT;
class AliHLTEsdManager;
class AliHLTPluginBase;
/**
 * @defgroup alihlt_aliroot_reconstruction AliRoot reconstruction.
 *
 * @section alihlt_aliroot_reconstruction_intro General Remarks
 * Like all other ALICE detectors, HLT utilizes the AliReconstruction interface
 * to implement a plugin for the AliRoot reconstruction. The reconstructor can be
 * used to
 * -# run HLT analysis chains in the AliRoot reconstruction <br>
 *    This option is mainly intended for the <em>development and debugging cycle</em>. 
 *    HLT chains can be defined by means of AliHLTConfiguration and can be run either
 *    stand-alone or embedded into the AliReconstruction cycle.
 * -# run the default analysis chains <br>
 *    HLT modules can define default analysis chains to be run during AliRoot
 *    reconstruction.
 * -# handle the HLTOUT data<br>
 *    The HLT output stream contains multiple data blocks produced by the various
 *    components of the HLT chain. Each block might need different and even
 *    detector specific processing, like e.g. the processing of ESD objects or the
 *    handling of compressed data.
 *
 * @section alihlt_aliroot_reconstruction_steering Steering
 * The AliHLTReconstructor provides the main interface for the reconstruction. 
 * The handling of the HLTOUT data is described in AliRawReaderHLT.
 *
 * @section alihlt_aliroot_reconstruction_examples Examples
 * @subsection alihlt_aliroot_reconstruction_examples_reco Run chains
 * - @ref tut_reconstruction
 *
 * @subsection alihlt_aliroot_reconstruction_examples_hltout Handle HLTOUT
 * - @ref tut_alirawreaderhlt
 *
 * @ingroup alihlt_system
 */

/**
 * @class AliHLTReconstructor
 * AliRoot event reconstruction plug-in for the HLT.
 * The AliHLTReconstructor holds an instance of the @ref AliHLTSystem
 * steering class. The actual reconstruction depends on the loaded component
 * libraries. Each library must implement a module agent (@ref AliHLTModuleAgent)
 * in order to provide information on the supported features and the
 * configurations to be run.
 *
 * The AliHLTReconstructor provides both the functionality to run customized
 * analysis chains and the treatment of the HLTOUT data.
 *
 * @section sec_alihltreconstructor_options Options
 * The default component libraries which are loaded through the initialization
 * are determined by the @ref AliHLTSystem::fgkHLTDefaultLibs array. The library
 * loading can be overridden by an option to the AliHLTReconstructor through the
 * <tt>SetOption</tt> method of <tt>AliReconstruction</tt>, e.g.
 * <pre>
 * AliReconstruction rec;
 * rec.SetOption("HLT", "libAliHLTSample.so");
 * </pre>
 * will only load <tt>libAliHLTSample.so</tt>
 * 
 * Optional arguments:<br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li loglevel=<i>level</i><br>
 *     level can be a hex number encoding the @ref AliHLTComponentLogSeverity
 * \li alilog=off <br>
 *     disables the logging of HLT log messages through <tt>AliLog</tt> <br>
 * \li ignore-hltout <br>
 *     ignore data from the HLTOUT data links
 * \li ignore-ctp <br>
 *     ignore CTP trigger setup
 * \li esdmanager=<option> <br>
 *     options passed to the AliHLTEsdManager
 *
 * For further information on the AliRoot reconstruction refer to the AliRoot
 * documentation, namely <tt>AliReconstruction</tt>.
 *
 * @section sec_alihltreconstructor_chains Custom reconstruction chains
 * In order to run an HLT chain during the AliRoot reconstruction, a chain
 * configuration must be defined. This can be done by
 * - specifying a configuration macro defining a configuration macro and
 *   the name of the chain as parameters
 *   <pre>
 *   rec.SetOption("HLT", "config=[macro.C] chains=[name]")
 *   </pre>
 * - providing the configuration and the name by the module agent.
 *   AliHLTModuleAgent and the functions AliHLTModuleAgent::CreateConfigurations
 *   and AliHLTModuleAgent::GetReconstructionChains
 *
 * @section sec_alihltreconstructor_hltout Treatment of HLTOUT data.
 * The HLTOUT data is a collation of output blocks produced by the various
 * components running in an HLT chain. Typically its the output of the last
 * component(s) in the chain, or components specifically connected to the HLT
 * output.
 *
 * The treatment of the HLTOUT data blocks is implemented in handlers of type
 * AliHLTOUTHandler. The AliHLTModuleAgent of the module  creates the appropriate
 * handler for a data block.
 * The data type of the individual blocks is set by the producer component and
 * specifies the nature of the data processing. There are 5 overall groups:
 * - output is in ESD format:
 *      @ref sec_alihltreconstructor_hltout_esd
 * - data describes DDL raw format:
 *      @ref sec_alihltreconstructor_hltout_rawreader
 * - pre-analyzed data to be fed into the normal reconstruction:
 *      @ref sec_alihltreconstructor_hltout_rawstream
 * - data is fed into an analysis chain:
 *      @ref sec_alihltreconstructor_hltout_chain
 * - detector specific handler:
 *      @ref sec_alihltreconstructor_hltout_proprietary
 *
 * @subsection sec_alihltreconstructor_hltout_esd ESD HLTOUT data
 * The framework implements a standard handling of ESD data
 * blocks of type ::kAliHLTDataTypeESDObject {ALIESDV0:ANY} and 
 * ::kAliHLTDataTypeESDTree {ESD_TREE:ANY}. Please note that the V0 refers to
 * a foreseen version number, not the AliESDV0 class. \em ANY can be
 * any detector origin. Each ESD block contains the data of only one event,
 * the ESDs are merged by the AliHLTEsdManager into the hltEsd. Optionally,
 * ESD contributions are written to files following the naming scheme
 * AliHLT\em DET ESDs.root. This debugging feature can be enabled by option
 * esdmanager='-writelocal -directory=dir'. The specification of target
 * directory is optional.
 *
 * The module agent can provide a handler for multiple ESD data blocks, e.g.
 * for merging within one event prior to the writing. Instead of the individual
 * ESDs the one provided by the handler is passed to the AliHLTEsdManager. The
 * handler is of type \link AliHLTModuleAgent::AliHLTOUTHandlerType kEsd \endlink.
 * 
 * A specific handler AliHLTOUTHandlerEsdBranch allows to simply merge a
 * streamed Root object into the ESD. The class can be used as-is by just
 * specifying the data type and specification of the relevant data block and
 * the branch name. Alternatively, a child class can implement conversion of
 * binary data to a Root object or merging of several data blocks into one
 * object to be added to the ESD. \b Note: in order to create the branch at
 * the beginning of reconstruction the hltEsd layout needs to be adjusted.
 *
 * @subsection sec_alihltreconstructor_hltout_rawreader DDL raw HLTOUT data
 * The HLT can perform selective readout and produces a reduced amount of data
 * in the original raw ddl format. In order to feed this data from the HLTOUT
 * DDL links into the normal reconstruction, a handler of type 
 * \link AliHLTModuleAgent::AliHLTOUTHandlerType kRawReader \endlink must be
 * implemented and provided by the
 * module agent. The handler has to derive the original equipment id from the
 * data type and specification of the block. The offline reconstruction does
 * not need to be changed or adapted at all. See AliRawReaderHLT for details.
 *
 * @subsection sec_alihltreconstructor_hltout_rawstream Preprocessed Raw HLTOUT data
 * Handlers of type \link AliHLTModuleAgent::AliHLTOUTHandlerType kRawStream \endlink
 * are foreseen though at the time of writing (May 08) the
 * concept is not fixed. Advanced data compression algorithms can produce a
 * raw data format which is not convertible into the raw DDL data, e.g. lossy
 * compression techniques storing clusters parametrized regarding to tracks. A
 * specific RawStream is needed here since the data is detector specific and the
 * first stage of the offline reconstruction might need some adaptions.
 *
 * @subsection sec_alihltreconstructor_hltout_chain HLTOUT data fed into a chain
 * Handlers of type \link AliHLTModuleAgent::AliHLTOUTHandlerType kChain \endlink
 * can execute a normal HLT chain and thus process HLTOUT data blocks by normal
 * HLT components just as if the components were running on-line. The base class
 * is provided by AliHLTOUTHandlerChain and can be used as it is just specifying
 * the chain to be run.<br>
 * Example:
 *
 * @subsection sec_alihltreconstructor_hltout_proprietary Proprietary HLTOUT data
 * This is a handler of proprietary detector data. Handlers of type 
 * \link AliHLTModuleAgent::AliHLTOUTHandlerType kProprietary \endlink
 * do not have any standard output to the framework. Data can be processed and
 * stored to files.
 *
 * @section sec_alihltreconstructor_helper Tools and helper functions
 * Some helper functions of the AliHLTReconstruction can be used in stand-alone
 * mode. Remember to Init() the reconstructor.
 * <pre>
 * {
 * gSystem->Load("libHLTrec.so");
 * AliHLTReconstructor rec;
 * rec.Init();
 * // do something
 * }
 * </pre>
 * @subsection sec_alihltreconstructor_hltout_standalone Stand-alone HLTOUT processing
 * - HLTOUT processing from a digit file:
 * <pre>
 *  void ProcessHLTOUT(const char*, AliESDEvent*) const;
 * </pre>
 * - HLTOUT processing from an AliRawReader
 * <pre>
 *  void ProcessHLTOUT(AliRawReader*, AliESDEvent*) const;
 * </pre>
 *
 * @ingroup alihlt_aliroot_reconstruction
 * @section sec_alihltreconstructor_members Class members
 */
class AliHLTReconstructor: public AliReconstructor {
public:
  /** standard constructor */
  AliHLTReconstructor();
  /** constructor */
  AliHLTReconstructor(const char* options);
  /** destructor */
  virtual ~AliHLTReconstructor();

  /** init the reconstructor */
  void Init();

  /** init the reconstructor */
  void Init(const char* options);

  /**
   * Init streamer infos for the relevent classes to be extracted from HLT raw
   * data payload. Reads the info from HLT/Calib/StreamerInfo
   */
  int InitStreamerInfos();

  /**
   * Init streamer infos for the relevant classes to be extracted from HLT raw
   * data payload.
   */
  int InitStreamerInfos(TObjArray* pSchemas) const;

  /**
   * Build the CTP_TRIGGER_CLASSES string from CTP trigger configuration
   */
  int BuildCTPTriggerClassString(TString& triggerclasses) const;

  /**
   * This Reconstructor function is not applicable for the AliHLTReconstructor
   * as it gets a detector specific digits tree. But HLT processes all detectors.
   * Furthermore it's purely simulated data. <br>
   * The function forwards to the default bahavior of AliReconstructor but gives
   * a warning if there were options set, i.e. the user runs customized
   * reconstruction.
   *
   * @note HLT reconstruction on simulated data is processed at the end of
   * simulation. <br>
   */
  void Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  /**
   * Reconstruction from RAW data.
   * The rawReader holds data for all detectors and this version of Reconstruct
   * is thus applicable for the HLT. The clustersTree is just ignored.
   */
  void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  /**
   * This function treats the simulated HLTOUT data.
   * Opens a handler for simulated HLTOUT data and forwards to ::ProcessHLTOUT.
   */
  void FillESD(TTree* digitsTree, TTree* clustersTree, AliESDEvent* esd) const;

  /**
   * Process the raw HLTOUT data and fill ESD.
   * Opens a handler for raw HLTOUT data and forwards to ::ProcessHLTOUT.
   */
  void FillESD(AliRawReader* rawReader, TTree* clustersTree, AliESDEvent* esd) const;

  /**
   * Process HLTOUT data and fill ESD.
   * This is the final treatment of the HLTOUT data, either simulated or real.
   * HLTOUT data is stored in HOMER format, the AliHLTOUT object provides the interface
   * to the individual data blocks.
   *
   * During reconstruction (::Reconstruct), module or user defined chains can be
   * processed and may add additional data to the HLTOUT object. This data is then
   * treated in the same way.
   */
  void ProcessHLTOUT(AliHLTOUT* pHLTOUT, AliESDEvent* esd, bool bVerbose=false) const;

  /**
   * Process HLTOUT data.
   * Open digit file and process the HLTOUT digit data.
   * This function is mostly intended for debugging purposes and stand-alone
   * processing of the output from the simulation. Loops over all events.
   * @param digitFile        path of the digit file
   * @param pEsd             optional ESD to be filled
   */
  void ProcessHLTOUT(const char* digitFile="HLT.Digits.root", AliESDEvent* pEsd=NULL) const;

  /**
   * Process HLTOUT data.
   * Process the HLTOUT from the raw reader.
   * This function is mostly intended for debugging purposes and stand-alone
   * processing of simulated or real raw data. 
   * \em Note: Loops over all events, i.e. the current event of the the raw
   * reader will change. Not to be called inside the normal AliRoot processsing.
   * @param pRawReader       raw reader instance
   * @param pEsd             optional ESD to be filled
   */
  void ProcessHLTOUT(AliRawReader* pRawReader, AliESDEvent* pEsd=NULL) const;

  /**
   * Print a short info about the HLTOUT content.
   */
  void PrintHLTOUTContent(AliHLTOUT* pHLTOUT) const;

  enum {
    /// ignore the blocks from the HLTOUT payload
    kAliHLTReconstructorIgnoreHLTOUT = 0x1,
    kAliHLTReconstructorIgnoreCTP = 0x2,
    kAliHLTReconstructorLastFlag
  };

private:
  /** copy constructor prohibited */
  AliHLTReconstructor(const AliHLTReconstructor& src);
  /** assignment operator prohibited */
  AliHLTReconstructor& operator=(const AliHLTReconstructor& src);

  /** ESD manger instance for this reconstruction */
  AliHLTEsdManager* fpEsdManager; //!transient

  /** base class for AliRoot HLT plugins */
  AliHLTPluginBase* fpPluginBase;                                     //!transient

  UInt_t fFlags; //! transient

  static const char* fgkCalibStreamerInfoEntry; //! OCDB path

  ClassDef(AliHLTReconstructor, 8)   // class for the HLT reconstruction

};

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
