#ifndef ALIANALYSISTASKMUONCDBCONNECT_H
#define ALIANALYSISTASKMUONCDBCONNECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup muondep
/// \class AliAnalysisTaskMuonCDBConnect
/// \brief task to setup OCDB path and load OCDB data
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMuonCDBConnect : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonCDBConnect();
  AliAnalysisTaskMuonCDBConnect(const char *name);
  virtual ~AliAnalysisTaskMuonCDBConnect() {}
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  /// Set location of the specific OCDB storage for MUON alignment
  void SetAlignStorage(const char* ocdbPath) { fAlignStorage = ocdbPath; }
  
  /// Set location of the specific OCDB storage for MUON recoParam
  void SetRecoParamStorage(const char* ocdbPath) { fRecoParamStorage = ocdbPath; }
  
  /// Set the flag to load the magnetic field
  void LoadMagField() { fLoadMagField = kTRUE; }
  
  /// Set the flag to load the geometry
  void LoadGeometry() { fLoadGeometry = kTRUE; }
  
  /// Set the flag to load the mapping (only the segmentation if specified)
  void LoadMapping(Bool_t segmentationOnly = kFALSE) { fLoadMapping = kTRUE; fSegOnly = segmentationOnly; }
  
  virtual void   UserCreateOutputObjects() {}
  virtual void   UserExec(Option_t *) {}
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *) {}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonCDBConnect(const AliAnalysisTaskMuonCDBConnect& rhs);
  /// Not implemented
  AliAnalysisTaskMuonCDBConnect& operator = (const AliAnalysisTaskMuonCDBConnect& rhs);
  
private:
  
  TString fDefaultStorage;   ///< location of the default OCDB storage
  TString fAlignStorage;     ///< location of the specific OCDB storage for MUON alignment
  TString fRecoParamStorage; ///< location of the specific OCDB storage for MUON recoParam
  Bool_t  fLoadMagField;     ///< flag to load the magnetic field
  Bool_t  fLoadGeometry;     ///< flag to load the geometry
  Bool_t  fLoadMapping;      ///< flag to load the mapping
  Bool_t  fSegOnly;          ///< flag to load only the segmentation
  Bool_t  fOCDBSet;          //!< flag to setup the OCDB only once
  
  ClassDef(AliAnalysisTaskMuonCDBConnect, 1);
};

#endif

