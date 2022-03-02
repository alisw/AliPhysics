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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskMCPredictionsEE_H
#define AliAnalysisTaskMCPredictionsEE_H

class TTree;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliMCEvent;
class TDatabasePDG;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCPredictionsEE : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskMCPredictionsEE();
  AliAnalysisTaskMCPredictionsEE(const char *name);
  virtual ~AliAnalysisTaskMCPredictionsEE();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  //Check MC type
  Bool_t IsHijing()  const;
  Bool_t IsDPMJet()  const;
  Bool_t IsEPOSLHC() const;
    
  void SetSelectINELgtZERO ( Bool_t lOpt ) { fkSelectINELgtZERO = lOpt; }

  void SetEtaThreshold(float val) { fEtaThreshold = val; }
  float GetEtaThreshold() const { return fEtaThreshold; }
  void SetEnergyThreshold(float val) { fEnergyThreshold = val; }
  float GetEnergyThreshold() const { return fEnergyThreshold; }
  void SetEtaBarrel(float val) { fEtaBarrel = val; }
  float GetEtaBarrel() const { return fEtaBarrel; }

  //---------------------------------------------------------------------------------------
  void loopMC(AliMCEvent *mcEvent);

private:
  // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
  // your data member object is created on the worker nodes and streaming is not needed.
  // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TTree *fTree; //! tree with event infos
  TDatabasePDG *fDB = nullptr; //!

  Bool_t fkSelectINELgtZERO;

  float fEtaThreshold = 8.;
  float fEnergyThreshold = 0.;
  float fEtaBarrel = 0.5;

  // tree variables
  Int_t fNMPI = 0;           //! number of multiparton interaction
  Int_t fMC_NPart = 0;       //! number of participants
  Int_t fMC_NColl = 0;       //! number of collisions
  Float_t fMC_b = 0.;        //! impact parameter
  int   fInelGT0 = 0;        //!

  Int_t fNchEta = 0;          //! multiplicity in central region
  Int_t fNLambdaEta = 0;      //!
  Int_t fNXiEta = 0;          //!
  Int_t fNAntiXiEta = 0;      //!
  Int_t fNOmegaEta = 0;       //!
  Int_t fNPiEta = 0;          //!
  Int_t fNPi0Eta = 0;         //!
  Int_t fNKchEta = 0;         //!
  Int_t fNK0Eta = 0;          //!
  float fEffEnergy = 0.;      //!
    
  AliAnalysisTaskMCPredictionsEE(const AliAnalysisTaskMCPredictionsEE&);            // not implemented
  AliAnalysisTaskMCPredictionsEE& operator=(const AliAnalysisTaskMCPredictionsEE&); // not implemented
  
  ClassDef(AliAnalysisTaskMCPredictionsEE, 1);
};

#endif

