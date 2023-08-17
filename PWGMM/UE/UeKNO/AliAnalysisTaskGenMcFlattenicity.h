/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef ALIANALYSISTASKGENMCFLATTENICITY_H
#define ALIANALYSISTASKGENMCFLATTENICITY_H
// ROOT includes

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TSystem.h>

// AliRoot includes

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliInputEventHandler.h"
#include <AliAnalysisTaskSE.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class TList;
class TH1F;
class TH1I;
class TH2D;
class TH3D;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskGenMcFlattenicity : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskGenMcFlattenicity();
  AliAnalysisTaskGenMcFlattenicity(const char *name);

  virtual ~AliAnalysisTaskGenMcFlattenicity();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void SetPtMin(Double_t val) { fPtMin = val; }
  virtual void SetEtaMax(Double_t val) { fEtaCut = val; }
  virtual void SetIsPP(Bool_t val) { fIsPP = val; }
  virtual void SetGenerator(int val) { fGenerator = val; }

private:
  void GetGenLeadingObject();
  void GetGenUEObservables();
  void GetFlattenicity();
  void MakeAnalysis();
  int GetPidCode(int pdgCode);

protected:
  Bool_t IsMCEventSelected(TObject *obj);
  AliMCEvent *fMC;                  //!<! MC event
  AliInputEventHandler *fMcHandler; //!<! MCEventHandler
  AliStack *fMCStack;

  int fGenerator;
  double fEtaCut;
  bool fIsPP;
  double fPtMin;
  bool fIsINEL0;
  int fNchv0;
  int fNchv0a;
  int fNchv0c;
  double fFlattV0;
  TH1D *hFlatt;
  TH2D *hflatVsNchV0;
  TH2D *hMultVsFlat;
  TH3D *hMultVsFlatVsPt[4];
  TH2D *hFlatVsPt[4];
  TList *fOutputList; //!<! Output list of objects

  AliAnalysisTaskGenMcFlattenicity(const AliAnalysisTaskGenMcFlattenicity &);
  AliAnalysisTaskGenMcFlattenicity &
  operator=(const AliAnalysisTaskGenMcFlattenicity &);

  ClassDef(AliAnalysisTaskGenMcFlattenicity, 1);
};

#endif
