
#ifndef ALIPHOTONISOLATION_H
#define ALIPHOTONISOLATION_H

#include "AliAnalysisTaskSE.h"
#include <vector>
#include <map>
#include <utility>

using namespace std;

//inherits properties of AliAnalysisTaskSE
class AliPhotonIsolation : public AliAnalysisTaskSE {

 public:
  AliPhotonIsolation(const char *name="PhotonIsolation_0", Int_t photonType=0);
  //Uncopyable & operator=(const Uncopyable&);

  virtual ~AliPhotonIsolation();                            //virtual destructor
  void UserCreateOutputObjects();

  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  TList* GetPhotonIsolationHistograms() {return fListHistos;}

  Bool_t GetIsolation(Int_t clusterID, Float_t R, Float_t isoPt);

 private:

  AliPhotonIsolation (const AliPhotonIsolation&); // not implemented
  AliPhotonIsolation & operator=(const AliPhotonIsolation&); // not implemented

  void SetV0ReaderName(TString name)    {fV0ReaderName = name; return;}
  void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting; return;}

  // private methods
  void Initialize();
  void ProcessEvent(AliVEvent *event);

  // debug methods
  void DebugIsolation();

  // basic variables/objects
  Int_t                 fPhotonType;             //
  TString               fV0ReaderName;           // Name of V0Reader
  TString               fCorrTaskSetting;        // Name of Corr Task Setting
  map<Int_t,Float_t> fMapClustertoPtR1;    // Map cluster ID to pTsum in cone R=0.1
  map<Int_t,Float_t> fMapClustertoPtR2;    // Map cluster ID to pTsum in cone R=0.2
  map<Int_t,Float_t> fMapClustertoPtR3;    // Map cluster ID to pTsum in cone R=0.3
  map<Int_t,Float_t> fMapClustertoPtR4;    // Map cluster ID to pTsum in cone R=0.4

  //histos
  TList*                fListHistos;             // list with histogram(s)
  TH1F*                 fHistTest;
  TH1F*                 fHistClusterEnergy;
  TH1F*                 fHistIso;

  ClassDef(AliPhotonIsolation,1)
    };

#endif
