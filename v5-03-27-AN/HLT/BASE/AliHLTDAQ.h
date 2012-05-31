//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDAQ_H
#define ALIHLTDAQ_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDAQ.h
/// @author Matthias Richter
/// @date   24.10.2008
/// @brief  Virtual Interface to the AliDAQ class.
///

#include "Rtypes.h"
#include <string>

/**
 * Virtual interface to the AliDAQ class.
 * In order to keep the libHLTbase free of AliRoot dependencies, the
 * implementation has been separated from libHLTbase.
 * Implementation in libHLTrec.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_system
 */
class AliHLTDAQ {
 public:
  AliHLTDAQ();
  virtual ~AliHLTDAQ();
  static  Int_t       NumberOfDetectors();

  static  Int_t       DetectorID(const char *detectorName);
  static  const char *DetectorName(Int_t detectorID);
  // Note: use specific number instead of kAliHLTComponentDataTypefOriginSize to avoid including AliHLTDataTypes.h
  static  Int_t       DetectorIDFromHLTOrigin(const char dataorigin[4]);
  static  const char *DetectorName(const char dataorigin[4]);

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

  static const char *ListOfTriggeredDetectors(UInt_t detectorPattern);
  static UInt_t      DetectorPattern(const char *detectorList);

  static const char *OfflineModuleName(const char *detectorName);
  static const char *OfflineModuleName(Int_t detectorID);

  static const char *OnlineName(const char *detectorName);
  static const char *OnlineName(Int_t detectorID);

  static std::string HLTOrigin(const char *detectorName);
  static std::string HLTOrigin(Int_t detectorID);
  
  static std::string HLTSpecificationFromDdlID(Int_t ddlID);
  // Note: use specific number instead of kAliHLTComponentDataTypefOriginSize to avoid including AliHLTDataTypes.h
  static Int_t       DdlIDFromHLTBlockData(const char dataorigin[4], UInt_t specification);

  static AliHLTDAQ* GetInstance();

 private:
  virtual  Int_t       VirtNumberOfDetectors()=0;

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

  virtual  const char *VirtListOfTriggeredDetectors(UInt_t detectorPattern)=0;
  virtual  UInt_t      VirtDetectorPattern(const char *detectorList)=0;

  virtual  const char *VirtOfflineModuleName(const char *detectorName)=0;
  virtual  const char *VirtOfflineModuleName(Int_t detectorID)=0;

  virtual  const char *VirtOnlineName(const char *detectorName)=0;
  virtual  const char *VirtOnlineName(Int_t detectorID)=0;

  /** global instance */
  static AliHLTDAQ* fgpInstance; //!

  /** the name of the actual implementation */
  static const char* fgkImplName; //!

  /** the library of the implementation */
  static const char* fgkImplLibrary; //!

  /// mapping between HLT data origin and AliDAQ detector number
  static const char* fgkOriginMapping[]; //!
};

#endif //AliHLTDAQ
