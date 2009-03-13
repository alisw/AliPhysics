#ifndef ALITPCTASKQA_H
#define ALITPCTASKQA_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliTPCseed.h>
#include <TString.h>
class AliMCEventHadnler;
class TParticle;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent;
class AliMCEvent;
class THnSparse;
class TObjArray;

class AliTPCtaskQA : public AliAnalysisTask {
public:
  AliTPCtaskQA();
  AliTPCtaskQA(const char *name);
  AliTPCtaskQA(const AliTPCtaskQA& info);
  virtual ~AliTPCtaskQA();  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  void           Init();
  //
  void           ProcessMCInfo();
  //
  THnSparse * GetTPCqa(){return fTPCqa;}
  static      AliTPCtaskQA* ReadFromFile(const char *fname="OutputPID.root");
  //
  static void BinLogX(TAxis *axis);
protected:
  AliTPCtaskQA& operator=(const AliTPCtaskQA& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;          //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  //
  //
  // 
  TObjArray  *fList; //TList output object
  THnSparse * fTPCqa;         //raw tpc QA
  ClassDef(AliTPCtaskQA, 1); // Analysis task base class for tracks
};

#endif
