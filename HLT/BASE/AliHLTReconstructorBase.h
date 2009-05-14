// $Id$

#ifndef ALIHLTRECONSTRUCTORBASE_H
#define ALIHLTRECONSTRUCTORBASE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTReconstructorBase.h
    @author Matthias Richter
    @date   
    @brief  AliHLTPluginBase child for backward compatibility.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPluginBase.h"

/**
 * @class AliHLTReconstructorBase
 * This class was the former name of the base class for HLT
 * reconstruction instances. It has evolved to a more general
 * 'plugin' base class and the name has been changed to
 * AliHLTPluginBase
 */
class AliHLTReconstructorBase : public AliHLTPluginBase {
 public:
  AliHLTReconstructorBase();
  ~AliHLTReconstructorBase();

 protected:

 private:
  ClassDef(AliHLTReconstructorBase, 0)   // AliHLTPluginBase child for backward compatibility
};
#endif
