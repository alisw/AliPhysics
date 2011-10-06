/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Author : A. Mastroserio
//-----------------------------------------------------------------------

#ifndef ALIANALYSISTASKSPD_H
#define ALIANALYSISTASKSPD_H

#include "AliAnalysisTaskSE.h"

class TString;
class TList;

class AliITSsegmentationSPD;

class AliAnalysisTaskSPD : public AliAnalysisTaskSE {
 public:


  AliAnalysisTaskSPD();
  AliAnalysisTaskSPD(const Char_t* name);
  AliAnalysisTaskSPD& operator= (const AliAnalysisTaskSPD& c);
  AliAnalysisTaskSPD(const AliAnalysisTaskSPD& c);
  virtual ~AliAnalysisTaskSPD();

  // ANALYSIS FRAMEWORK 
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);

  void     SetOCDBInfo(UInt_t runNb, const char *location) {fRunNb=runNb; fOCDBLocation=location;}
  void     LoadGeometryFromOCDB(); 
  
  void     SetHeavyIonMode() {fHI=kTRUE;}
  void     SetTestMode() {fTest=kTRUE;} 

 protected:
  AliITSsegmentationSPD *fSegSPD;  
  TList          *fOutput   ;  // user histograms list
  UInt_t fRunNb;               // run number
  TString fOCDBLocation;       // ocdb path
  Bool_t fHI;                  // changes to the histo limits 
  Bool_t fTest;                // ocdb settings 

  ClassDef(AliAnalysisTaskSPD,2);
};

#endif
