//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDAQINTERFACEIMPLEMENTATION_H
#define ALIHLTDAQINTERFACEIMPLEMENTATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTDAQInterfaceImplementation.h
    @author Matthias Richter
    @date   
    @brief  Implementation of the AliHLTDAQInterfaceImplementation
*/

#include "AliHLTDAQ.h"

/**
 * @class AliHLTDAQInterfaceImplementation
 * Implementation of the AliHLTDAQVirtualInterface
 *
 * For the sake of library (in)dependencies, AliDAQ can not be used directly in
 * libHLTbase as this would  introduce dependencies to AliRoot libraries.
 * The AliHLTDAQVirtualInterface provides a virtual interface to AliDAQ with
 * the implementation in libHLTrec.so.
 * See AliHLTDAQVirtualInterface for usage.
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTDAQInterfaceImplementation : public AliHLTDAQ {
 public:
  /** constructor */
  AliHLTDAQInterfaceImplementation();
  /** destructor */
  virtual ~AliHLTDAQInterfaceImplementation();

  Int_t       VirtNumberOfDetectors();

  Int_t       VirtHLTId();
  Int_t       VirtDetectorID(const char *detectorName);
  const char *VirtDetectorName(Int_t detectorID);

  Int_t       VirtDdlIDOffset(const char *detectorName);
  Int_t       VirtDdlIDOffset(Int_t detectorID);

  const char *VirtDetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex);
  Int_t       VirtDetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex);

  Int_t       VirtDdlID(const char *detectorName, Int_t ddlIndex);
  Int_t       VirtDdlID(Int_t detectorID, Int_t ddlIndex);
  const char *VirtDdlFileName(const char *detectorName, Int_t ddlIndex);
  const char *VirtDdlFileName(Int_t detectorID, Int_t ddlIndex);

  Int_t       VirtNumberOfDdls(const char *detectorName);
  Int_t       VirtNumberOfDdls(Int_t detectorID);

  const char *VirtListOfTriggeredDetectors(UInt_t detectorPattern);
  UInt_t      VirtDetectorPattern(const char *detectorList);

  const char *VirtOfflineModuleName(const char *detectorName);
  const char *VirtOfflineModuleName(Int_t detectorID);

  const char *VirtOnlineName(const char *detectorName);
  const char *VirtOnlineName(Int_t detectorID);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTDAQInterfaceImplementation(const AliHLTDAQInterfaceImplementation&);
  /** assignment operator prohibited */
  AliHLTDAQInterfaceImplementation& operator=(const AliHLTDAQInterfaceImplementation&);

  ClassDef(AliHLTDAQInterfaceImplementation, 0)
};

#endif
