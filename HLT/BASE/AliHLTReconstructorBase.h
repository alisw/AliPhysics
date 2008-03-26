// $Id$

#ifndef ALIHLTRECONSTRUCTORBASE_H
#define ALIHLTRECONSTRUCTORBASE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTReconstructorBase.h
    @author Matthias Richter
    @date   
    @brief  Base class for HLT reconstruction classes.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"

class AliHLTSystem;

/**
 * @class AliHLTReconstructorBase
 * Base class for HLT reconstruction classes. AliHLTReconstructor and
 * AliRawReaderHLT both use the global AliHLTSystem instance. This
 * base class hosts the global instance.
 */
class AliHLTReconstructorBase {
 public:
  AliHLTReconstructorBase();
  /** destructor */
  virtual ~AliHLTReconstructorBase();

  /**
   * Init the global AliHLTSystem instance.
   */
  static void InitInstance();

  /**
   * Get the global AliHLTSystem instance.
   */
  static AliHLTSystem* GetInstance();

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTReconstructorBase(const AliHLTReconstructorBase& src);
  /** assignment operator prohibited */
  AliHLTReconstructorBase& operator=(const AliHLTReconstructorBase& src);

  static AliHLTSystem* fpSystem; //! HLT steering object

  static int fNofInstances;

  ClassDef(AliHLTReconstructorBase, 0)   // base class for the HLT reconstruction
};
#endif
