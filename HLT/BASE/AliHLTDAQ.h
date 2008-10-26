//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDAQ_H
#define ALIHLTDAQ_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTDAQ.h
    @author Matthias Richter
    @date   24.10.2008
    @brief  Virtual Interface to the AliDAQ class.
*/

#include "Rtypes.h"

/**
 * Virtual interface to the AliDAQ class.
 * In order to keep the libHLTbase free of AliRoot dependencies, the
 * implementation has been separated from libHLTbase.
 * Implementation in libHLTrec.
 */
class AliHLTDAQ {
 public:
  AliHLTDAQ();
  virtual ~AliHLTDAQ();
  static  Int_t       DetectorID(const char *detectorName);
  static  const char *DetectorName(Int_t detectorID);

  static  Int_t       DdlIDOffset(const char *detectorName);
  static  Int_t       DdlIDOffset(Int_t detectorID);

  static  const char *DetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex);
  static  Int_t       DetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex);

  static  Int_t       DdlID(const char *detectorName, Int_t ddlIndex);
  static  Int_t       DdlID(Int_t detectorID, Int_t ddlIndex);
  static  const char *DdlFileName(const char *detectorName, Int_t ddlIndex);
  static  const char *DdlFileName(Int_t detectorID, Int_t ddlIndex);

  static  Int_t       NumberOfDdls(const char *detectorName);
  static  Int_t       NumberOfDdls(Int_t detectorID);

  static AliHLTDAQ* GetInstance();

 private:
  virtual  Int_t       VirtDetectorID(const char *detectorName)=0;
  virtual  const char *VirtDetectorName(Int_t detectorID)=0;

  virtual  Int_t       VirtDdlIDOffset(const char *detectorName)=0;
  virtual  Int_t       VirtDdlIDOffset(Int_t detectorID)=0;

  virtual  const char *VirtDetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex)=0;
  virtual  Int_t       VirtDetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex)=0;

  virtual  Int_t       VirtDdlID(const char *detectorName, Int_t ddlIndex)=0;
  virtual  Int_t       VirtDdlID(Int_t detectorID, Int_t ddlIndex)=0;
  virtual  const char *VirtDdlFileName(const char *detectorName, Int_t ddlIndex)=0;
  virtual  const char *VirtDdlFileName(Int_t detectorID, Int_t ddlIndex)=0;

  virtual  Int_t       VirtNumberOfDdls(const char *detectorName)=0;
  virtual  Int_t       VirtNumberOfDdls(Int_t detectorID)=0;

  /** global instance */
  static AliHLTDAQ* fgpInstance; //!

  /** the name of the actual implementation */
  static const char* fgkImplName; //!

  /** the library of the implementation */
  static const char* fgkImplLibrary; //!
};

#endif //AliHLTDAQ
