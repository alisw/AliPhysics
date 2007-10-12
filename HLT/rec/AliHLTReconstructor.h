// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTReconstructor.h
    @author Matthias Richter
    @date   
    @brief  Binding class for HLT simulation in AliRoot

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */

#include "AliReconstructor.h"

class AliHLTSystem;
class AliRawReader;
class AliESDEvent;

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
 * are determined by the @ref kHLTDefaultLibs array. The library loading can
 * be overridden by an option to the AliHLTReconstructor through the
 * <tt>SetOption</tt> method of <tt>AliReconstruction</tt>, e.g.
 * <pre>
 * AliReconstruction rec;
 * rec.SetOption("HLT", "libAliHLTSample.so");
 * </pre>
 * will only load <tt>libAliHLTSample.so</tt>
 * 
 * Optional arguments:<br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li loglevel=<i>level</i><br>
 *     level can be a hex number encoding the @ref AliHLTComponentLogSeverity
 * \li alilog=off <br>
 *     disables the logging of HLT log messages through <tt>AliLog</tt> <br>
 *
 * For further information on the AliRoot reconstruction refer to the AliRoot
 * documentation, namely <tt>AliReconstruction</tt>.
 */
class AliHLTReconstructor: public AliReconstructor {
public:
  AliHLTReconstructor();
  /** destructor */
  virtual ~AliHLTReconstructor();

  /** init the reconstructor */
  void Init();

  /**
   * This Reconstructor function is not applicable for the AliHLTReconstructor
   * as it gets a detector specific digits tree. But HLT processes all detectors.
   * Furthermore it's purely simulated data. <br>
   * The function forwards to the default bahavior of AliReconstructor but gives
   * a warning if there were options set, i.e. the user runs customized
   * reconstruction.
   */
  void Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  /**
   * Reconstruction from RAW data.
   * The rawReader holds data for all detectors and this version of Reconstruct
   * is thus applicable for the HLT. The clustersTree is just ignored.
   */
  void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  /**
   * This function is purely for simulated data and not applicable for HLT.
   * HLT reconstruction on simulated data is processed at the end of
   * simulation. <br>
   * The function forwards to the default bahavior of AliReconstructor but gives
   * a warning if there were options set, i.e. the user runs customized
   * reconstruction.
   */
  void FillESD(TTree* digitsTree, TTree* clustersTree, AliESDEvent* esd) const;

  /**
   * Fill the ESD from RAW data.
   * This is the main entry for HLT reconstruction of RAW data. It performs both
   * the analysis by the defined chains and the filling of the ESD.
   */
  void FillESD(AliRawReader* rawReader, TTree* clustersTree, AliESDEvent* esd) const;

private:
  /** copy constructor prohibited */
  AliHLTReconstructor(const AliHLTReconstructor& src);
  /** assignment operator prohibited */
  AliHLTReconstructor& operator=(const AliHLTReconstructor& src);

  AliHLTSystem* fpSystem; //! HLT steering object

  ClassDef(AliHLTReconstructor, 3)   // class for the HLT reconstruction
};

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
