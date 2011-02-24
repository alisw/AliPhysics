#ifndef ALIANALYSISTASKCHECKV0TENDERII_H
#define ALIANALYSISTASKCHECKV0TENDERII_H

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
// Task for perfomance studies of V0 selection code
// More information can be found in the source file
//
#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1F;
class TList;
class AliHFEcollection;
class AliESDv0KineCuts;
class AliKFVertex;

class AliAnalysisTaskCheckV0tenderII : public AliAnalysisTaskSE{
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

  AliAnalysisTaskCheckV0tenderII();
  AliAnalysisTaskCheckV0tenderII(const Char_t *name);
  ~AliAnalysisTaskCheckV0tenderII();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
 
  
  void ProcessV0s();
  void ProcessDaughters(AliESDv0 * const v0);
  void ProcessV0sMC(AliESDv0 * const v0);
  void ProcessDaughtersMC(AliESDv0 * const v0);
  void ProcessBackground(AliESDv0 * const v0);

 private:
  AliAnalysisTaskCheckV0tenderII(const AliAnalysisTaskCheckV0tenderII &ref);
  AliAnalysisTaskCheckV0tenderII &operator=(const AliAnalysisTaskCheckV0tenderII &ref);

  Float_t MassV0(AliESDv0 * const v0, Int_t id);
  Bool_t  CheckSigns(AliESDv0 * const v0);

  Int_t   PDGtoPIDv0(Int_t pdgV0) const;
  Int_t   PDGtoPID(Int_t pdg) const;  
  
  void    ResetPDGcodes();
  
  TList              *fOutput;        //! Container for output histos
  AliHFEcollection   *fColl;          //! collection of Data output
  AliHFEcollection   *fCollMC;        //! collection of MC output
  AliESDv0KineCuts   *fV0cuts;        //! standalone V0 selection class
  AliKFVertex        *fPrimaryVertex; //! primary vertex of the current event

  Int_t               fpdgV0;  // PDG code of the reconstructed V0
  Int_t               fpdgP;   // PDG code of the positive daughter (sign corrected)
  Int_t               fpdgN;   // PDG code of the negative daughter (sign coreccted)
  
  TH1 *fEvents;              //! Number of Events

  ClassDef(AliAnalysisTaskCheckV0tenderII, 1)

};


#endif
