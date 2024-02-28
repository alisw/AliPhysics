/**
 *
 *  \class AliAnalysisTaskMCGenSpeedOfSound
 *
 *  Antonio Ortiz (ICN-UNAM), antonio.ortiz@nucleares.unam.mx
 *  First version: 	April 23, 2020
 *
 */
#ifndef ALIANALYSISTASKMCGENSPEEDOFSOUND_H
#define ALIANALYSISTASKMCGENSPEEDOFSOUND_H

// ROOT includes

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TList.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TSystem.h>

// AliRoot includes

#include <AliAnalysisTaskSE.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliInputEventHandler.h"

class TList;
class TH1D;
class TH1I;
class TH2D;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskMCGenSpeedOfSound : public AliAnalysisTaskSE {  //

 public:
  AliAnalysisTaskMCGenSpeedOfSound();
  AliAnalysisTaskMCGenSpeedOfSound(const char* name);

  virtual ~AliAnalysisTaskMCGenSpeedOfSound();

  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

 private:
  void GetMultipliciy();
  int GetPidCode(Int_t pdgCode) const;

 protected:
  bool IsMCEventSelected(TObject* obj);

  AliMCEvent* fMcEvent;              //!<! MC event
  AliInputEventHandler* fMcHandler;  //!<! MCEventHandler

  AliStack* fStack;
  TH1D* hNch08;
  TH1D* hNchV0M;
  TH1D* hNchEtaNeg;
  TH1D* hNchEtaPos;
  TH1D* hNchTPCEtaGap;
  TH1D* hNchSPDEtaGap;
  TH2D* hPtvsNch08;
  TH2D* hPtvsV0M;
  TH2D* hPtvsNchEtaPos;
  TH2D* hPtvsNchEtaNeg;
  TH2D* hPtvsNchTPCEtaGap;
  TH2D* hPtvsNchSPDEtaGap;

  TList* fListOfObjects;  //!<! Output list of objects

  AliAnalysisTaskMCGenSpeedOfSound(
      const AliAnalysisTaskMCGenSpeedOfSound&);  // not implemented
  AliAnalysisTaskMCGenSpeedOfSound& operator=(
      const AliAnalysisTaskMCGenSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskMCGenSpeedOfSound,
           1);  // Analysis task for LF spectra analysis
};

#endif
