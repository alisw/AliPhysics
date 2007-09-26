// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliReconstructor.h"

class AliHLTSystem;

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
  AliHLTReconstructor(Bool_t doTracker, Bool_t doHough);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTReconstructor(const AliHLTReconstructor& src);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTReconstructor& operator=(const AliHLTReconstructor& src);
  /** destructor */
  virtual ~AliHLTReconstructor();

  /** init the reconstructor */
  void Init();

  /** create a tracker */
  // Deprecated and must be removed.
  //  AliTracker*  CreateTracker() const;

  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const{
    AliReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(rawReader,clustersTree);
  }

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESDEvent* esd) const {
    AliReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESDEvent* esd) const {
    AliReconstructor::FillESD(rawReader,clustersTree,esd);
  }
  void SetDoBench(Bool_t b){fDoBench=b;}
  void SetDoCleanup(Bool_t b){fDoCleanUp=b;}
  
  // Deprecated and must be removed.
//  virtual void         FillDHLTRecPoint(AliRawReader* rawReader, Int_t nofEvent, Int_t dcCut) const;

private:
/*   void ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const; */
/*   void ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const; */
/*   void FillESDforConformalMapping(AliESDEvent* esd,Int_t iEvent) const; */
/*   void FillESDforHoughTransform(AliESDEvent* esd,Int_t iEvent) const; */

  Bool_t fDoHough;   //do the hough transform
  Bool_t fDoTracker; //do the standard conformal tracker
  Bool_t fDoBench;   //store the benchmark results
  Bool_t fDoCleanUp; //delete tmp tracking files

  AliHLTSystem* fpSystem; //! HLT steering object

  ClassDef(AliHLTReconstructor, 2)   // class for the HLT reconstruction
};

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
