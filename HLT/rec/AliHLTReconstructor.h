// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTReconstructor.h
    @author Matthias Richter
    @date   
    @brief  Binding class for HLT simulation in AliRoot
*/

#include "AliReconstructor.h"
#include "AliHLTReconstructorBase.h"

class AliHLTSystem;
class AliRawReader;
class AliESDEvent;
class AliHLTOUT;
class AliHLTEsdManager;

/**
 * @class AliHLTReconstructor
 * AliHLTReconstructor AliRoot event reconstruction plugin for the HLT.
 * The AliHLTReconstructor holds an instance of the @ref AliHLTSystem
 * steering class. The actual reconstruction depends on the loaded component
 * libraries. Each library must implement a module agent (@ref AliHLTModuleAgent)
 * in order to provide information on the supported features and the
 * configurations to be run.
 *
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
 *
 * For further information on the AliRoot reconstruction refer to the AliRoot
 * documentation, namely <tt>AliReconstruction</tt>.
 */
class AliHLTReconstructor: public AliReconstructor, public AliHLTReconstructorBase {
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
  void ProcessHLTOUT(AliHLTOUT* pHLTOUT, AliESDEvent* esd) const;

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

private:
  /** copy constructor prohibited */
  AliHLTReconstructor(const AliHLTReconstructor& src);
  /** assignment operator prohibited */
  AliHLTReconstructor& operator=(const AliHLTReconstructor& src);

  /** function pointer: processing of HLTOUT data */
  void* fFctProcessHLTOUT; //!transient

  /** ESD manger instance for this reconstruction */
  AliHLTEsdManager* fpEsdManager; //!transient

  ClassDef(AliHLTReconstructor, 5)   // class for the HLT reconstruction

};

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
