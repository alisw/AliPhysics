#ifndef ALIANALYSISTASKSEMULTDQDEV_H
#define ALIANALYSISTASKSEMULTDQDEV_H

#include "AliAnalysisTaskSE.h"

class TList;
class TClonesArray;
class AliPicoDQheader;

class AliMuonPairCuts;
class AliMuonTrackCuts;

class AliAnalysisTaskSEMultDQdev : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskSEMultDQdev();
  AliAnalysisTaskSEMultDQdev(const char *s);
  virtual ~AliAnalysisTaskSEMultDQdev();
//=============================================================================

//virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option="");
  virtual void NotifyRun();
//=============================================================================

 private :

  AliAnalysisTaskSEMultDQdev &operator= (const AliAnalysisTaskSEMultDQdev &);
  AliAnalysisTaskSEMultDQdev(const AliAnalysisTaskSEMultDQdev &);
//=============================================================================

  Bool_t IsNotEventSelected();
//=============================================================================

  AliPicoDQheader *fHeader;   //!
  TClonesArray    *fDQClArr;  //!

  AliMuonPairCuts  *fCutsDimu; //!

  TList *fListOutputs; //!
//=============================================================================

  ClassDef(AliAnalysisTaskSEMultDQdev, 2);
};

#endif
