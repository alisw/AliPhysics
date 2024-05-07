/**
 *
 *  \class AliAnalysisTaskMCGenSpeedOfSound
 *
 *  Omar Vazqquez (UH),
 *  *  First version: 	May 7, 2024
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
  bool IsGoodVertex() const;

  AliMCEvent* fMcEvent;              //!<! MC event
  AliInputEventHandler* fMcHandler;  //!<! MCEventHandler
  AliStack* fStack;
  TH1D* hNchFull;
  TH1D* hNchV0M;
  TH1D* hTPCEtaGap;
  TH1D* hSPDEtaGap;
  TH1D* hNchEtaGap;
  TH1D* hSPDFull;
  TH1D* hSPDEtaAdj;
  TH1D* hSPDEtaGapW;
  TH1D* hTPCFull;
  TH1D* hEtFull;
  TH1D* hEtEtaGap;
  TH2D* hPtvsNchFull;
  TH2D* hPtvsV0M;
  TH2D* hPtvsTPCEtaGap;
  TH2D* hPtvsSPDEtaGap;
  TH2D* hPtvsNchEtaGap;
  TH2D* hPtvsSPDFull;
  TH2D* hPtvsSPDEtaAdj;
  TH2D* hPtvsSPDEtaGapW;
  TH2D* hPtvsTPCFull;
  TH2D* hPtvsEtFull;
  TH2D* hPtvsEtEtaGap;
  TList* fListOfObjects;  //!<! Output list of objects

  AliAnalysisTaskMCGenSpeedOfSound(
      const AliAnalysisTaskMCGenSpeedOfSound&);  // not implemented
  AliAnalysisTaskMCGenSpeedOfSound& operator=(
      const AliAnalysisTaskMCGenSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskMCGenSpeedOfSound,
           2);  // Analysis task for LF spectra analysis
};

#endif
