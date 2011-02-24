#ifndef ALIANALYSISTASKCHECKV0TENDER_H
#define ALIANALYSISTASKCHECKV0TENDER_H

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

/* $Id$ */ 

//
// Task for PID QA
// Using AliHFEpidQA and AliHFEMCpidQA
// More information can be found in the source file
//
#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1F;
class TList;
class AliHFEcollection;

class AliAnalysisTaskCheckV0tender : public AliAnalysisTaskSE{
 public:
   enum{ // Reconstructed V0
      kRecoGamma = 0,
      kRecoK0 = 1,
      kRecoLambda = 2,
      kRecoALambda = 3
      
      };
  enum{ // Identified Daughter particles
    kRecoElectron = 0,
      kRecoPionK0 = 1,
      kRecoPionL = 2,
      kRecoProton = 3
      };

  AliAnalysisTaskCheckV0tender();
  AliAnalysisTaskCheckV0tender(const Char_t *name);
  ~AliAnalysisTaskCheckV0tender();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
 
  
  void ProcessV0s();
  void ProcessDaughters();
  void ProcessV0sMC();
  void ProcessDaughtersMC();

  Int_t GetTenderPidV0(AliESDv0 * const v0);
  Int_t GetTenderPidDaughter(AliESDtrack * const track);

 private:
  AliAnalysisTaskCheckV0tender(const AliAnalysisTaskCheckV0tender &ref);
  AliAnalysisTaskCheckV0tender &operator=(const AliAnalysisTaskCheckV0tender &ref);

  Float_t MassV0(AliESDv0 * const v0, Int_t id);
  Bool_t CheckSigns(AliESDv0 * const v0);
  TList *fOutput;            //! Container for output histos
  AliHFEcollection *fColl;   //! collection of Data output
  AliHFEcollection *fCollMC; //! collection of MC output
  TH1 *fEvents;              //! Number of Events
  ClassDef(AliAnalysisTaskCheckV0tender, 1)

};

#endif
