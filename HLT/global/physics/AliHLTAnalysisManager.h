#ifndef ALIHLTANALYSISMANAGER_H
#define ALIHLTANALYSISMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch, 2014-12-15

//==============================================================================
//   AliHLTAnalysysManager 
//   implementation of the AliAnlysisManager for the HLT
//
//==============================================================================

#include "AliAnalysisManager.h"

class TList;

class AliHLTAnalysisManager : public AliAnalysisManager {

public:
  AliHLTAnalysisManager(const char* name="HLTanalysis", const char* title="");
  virtual ~AliHLTAnalysisManager();
  Bool_t InitAnalysis();
  Int_t WriteAnalysisToFile();
  Int_t ResetOutputData();

private:
  AliHLTAnalysisManager& operator=(const AliHLTAnalysisManager& that);
  AliHLTAnalysisManager(const AliHLTAnalysisManager& that);
  
  ClassDef(AliHLTAnalysisManager, 1)  // Analysis manager class
};

#endif

