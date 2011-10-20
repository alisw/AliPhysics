#ifndef ALIDAQ_H
#define ALIDAQ_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The AliDAQ class is responsible for handling all the information about   //
// Data Acquisition configuration. It defines the detector indexing,        //
// the number of DDLs and LDCs per detector.                                //
// The number of LDCs per detector is used only in the simulation in order  //
// to define the configuration of the dateStream application. Therefore the //
// numbers in the corresponding array can be changed without affecting the  //
// rest of the aliroot code.                                                //
// The equipment ID (DDL ID) is an integer (32-bit) number defined as:      //
// Equipment ID = (detectorID << 8) + DDLIndex                              //
// where the detectorID is given by fgkDetectorName array and DDLIndex is   //
// the index of the corresponding DDL inside the detector partition.        //
// Due to DAQ/HLT limitations, the ddl indexes should be consequtive, or    //
// at least without big gaps in between.                                    //
// The sub-detector code use only this class in the simulation and reading  //
// of the raw data.                                                         //
//                                                                          //
// cvetan.cheshkov@cern.ch  2006/06/09                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliDAQ: public TObject {
 public:

  AliDAQ() {};
  AliDAQ(const AliDAQ& source);
  AliDAQ& operator = (const AliDAQ& source);
  virtual ~AliDAQ() {};

  static Int_t       DetectorID(const char *detectorName);
  static const char *DetectorName(Int_t detectorID);

  static Int_t       DdlIDOffset(const char *detectorName);
  static Int_t       DdlIDOffset(Int_t detectorID);

  static const char *DetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex);
  static Int_t       DetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex);

  static Int_t       DdlID(const char *detectorName, Int_t ddlIndex);
  static Int_t       DdlID(Int_t detectorID, Int_t ddlIndex);
  static const char *DdlFileName(const char *detectorName, Int_t ddlIndex);
  static const char *DdlFileName(Int_t detectorID, Int_t ddlIndex);

  static Int_t       NumberOfDdls(const char *detectorName);
  static Int_t       NumberOfDdls(Int_t detectorID);

  static Float_t     NumberOfLdcs(const char *detectorName);
  static Float_t     NumberOfLdcs(Int_t detectorID);

  static void        PrintConfig();

  static const char *ListOfTriggeredDetectors(UInt_t detectorPattern);
  static UInt_t      DetectorPattern(const char *detectorList);
  static UInt_t      DetectorPatternOffline(const char *detectorList);

  static const char *OfflineModuleName(const char *detectorName);
  static const char *OfflineModuleName(Int_t detectorID);

  static const char *OnlineName(const char *detectorName);
  static const char *OnlineName(Int_t detectorID);

  enum {
// #ifdef MFT_UPGRADE
//     kNDetectors = 22,    // Number of detectors
// #else
//     kNDetectors = 21,    // Number of detectors
// #endif
    kNDetectors = 22,    // Number of detectors   // AU
    kHLTId = 30          // HLT detector index
  };
  
  enum DetectorBits {kSPD = 0x0001, kSDD = 0x0002, kSSD = 0x0004, kITS = 0x0007, 
		     kTPC = 0x0008, kTRD = 0x0010, kTOF = 0x0020, kHMPID = 0x0040, 
		     kPHOS = 0x0080, kCPV = 0x0100, kPMD = 0x0200, kMUONTRK = 0x0400,
		     kMUONTRG = 0x0800, kMUON = 0x0c00, kFMD = 0x1000, kT0 = 0x2000, kVZERO = 0x4000,
		     kZDC = 0x8000, kACORDE = 0x10000, kTRG = 0x20000, kEMCAL = 0x40000,
		     kDAQTEST = 0x80000, kHLT = 0x40000000
// #ifdef MFT_UPGRADE
// 		     ,kMFT = 0x80000000
// #endif
		     ,kMFT = 0x80000000     // AU
  };
  
private:
  
  static const char *fgkDetectorName[kNDetectors]; // Detector names
  static Int_t       fgkNumberOfDdls[kNDetectors]; // Number of DDLs per detector
  static Float_t     fgkNumberOfLdcs[kNDetectors]; // Number of LDCs per detector (not fixed - used only for the raw data simulation)
  static const char* fgkOfflineModuleName[kNDetectors]; // Names of the offline modules corresponding to the detectors
  static const char* fgkOnlineName[kNDetectors]; // Online (DAQ/ECS) detector names

  ClassDef(AliDAQ, 3)   // ALICE DAQ Configuration class
};

#endif
