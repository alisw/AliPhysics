#ifndef ALIANALYSISTASKEMCALJETV0CF_H
#define ALIANALYSISTASKEMCALJETV0CF_H

#include <TH1D.h>
#include "AliAnalysisTaskEmcalJet.h"
//=============================================================================

class TClonesArray;
class TVector3;

class AliAODEvent;
class AliESDEvent;
class AliCentrality;

class AliParticleContainer;
class AliClusterContainer;
class AliJetContainer;
//=============================================================================

class AliAnalysisTaskEmcalJetV0CF : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalJetV0CF();
  AliAnalysisTaskEmcalJetV0CF(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalJetV0CF();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t *opt);

  void SetKaCutNS(Double_t d) { fKaCutNS = d; }
  void SetLaCutNS(Double_t d) { fLaCutNS = d; }
  void SetHistoKshortInvM(TH1D const *h) { fHistoKshortInvM = new TH1D(*h); }
  void SetHistoLambdaInvM(TH1D const *h) { fHistoLambdaInvM = new TH1D(*h); }
  void SetHistoAntiLaInvM(TH1D const *h) { fHistoAntiLaInvM = new TH1D(*h); }
  void SetV0EtaRange(Double_t dMin, Double_t dMax) { fV0CutMinEta = dMin, fV0CutMaxEta = dMax; }
//=============================================================================

 protected :

  virtual void   ExecOnce();
  virtual Bool_t FillGeneralHistograms();
  virtual Bool_t FillHistograms();
  virtual Bool_t IsEventSelected();
  virtual Bool_t RetrieveEventObjects();
  virtual Bool_t Run();
//=============================================================================

 private :

  AliAnalysisTaskEmcalJetV0CF(const AliAnalysisTaskEmcalJetV0CF &);
  AliAnalysisTaskEmcalJetV0CF& operator=(const AliAnalysisTaskEmcalJetV0CF &);

  Bool_t FillRecoInfo();
  Bool_t FillKineInfo();
  void CreateUserOutputHistograms();

  Bool_t IsV0InJet(TVector3 vV0, Double_t dJetPtMin);

  Double_t fKaCutNS; //
  Double_t fLaCutNS; //

  Double_t fV0CutMinEta; //
  Double_t fV0CutMaxEta; //

  AliAODEvent   *fEventAOD;  //!
  AliESDEvent   *fEventESD;  //!
  AliCentrality *fCentInfo;  //!

  AliJetContainer      *fJetsContRD;          //!
  AliParticleContainer *fTracksContRD;        //!
  AliClusterContainer  *fCaloClustersContRD;  //!

  AliJetContainer      *fJetsContMC;    //!
  AliParticleContainer *fTracksContMC;  //!

  TClonesArray *fV0s;  //!

  TH1D *fHistoKshortInvM;  //!
  TH1D *fHistoLambdaInvM;  //!
  TH1D *fHistoAntiLaInvM;  //!

  TList *fListUserOutputs; //!

  ClassDef(AliAnalysisTaskEmcalJetV0CF, 2);
};

#endif
