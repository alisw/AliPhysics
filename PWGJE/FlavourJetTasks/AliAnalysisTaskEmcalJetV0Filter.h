#ifndef ALIANALYSISTASKEMCALJETV0FILTER_H
#define ALIANALYSISTASKEMCALJETV0FILTER_H

#include "AliAnalysisTaskEmcalJet.h"
//=============================================================================

class TString;
class TClonesArray;

class AliParticleContainer;
class AliClusterContainer;
class AliJetContainer;

class AliPicoHeaderCJ;
//=============================================================================

class AliAnalysisTaskEmcalJetV0Filter : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalJetV0Filter();
  AliAnalysisTaskEmcalJetV0Filter(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalJetV0Filter();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t *opt);

  void SetNameJetRD02(TString s) { fNameJetRD02 = s; }
  void SetNameJetRD03(TString s) { fNameJetRD03 = s; }
  void SetNameJetRD04(TString s) { fNameJetRD04 = s; }

  void SetNameJetMC02(TString s) { fNameJetMC02 = s; }
  void SetNameJetMC03(TString s) { fNameJetMC03 = s; }
  void SetNameJetMC04(TString s) { fNameJetMC04 = s; }

  void SetIsAnaPicoV0(Bool_t b)  { fIsAnaPicoV0 = b; }
  void SetAnaPicoV0MC(Bool_t b)  { fAnaPicoV0MC = b; }
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

  AliAnalysisTaskEmcalJetV0Filter(const AliAnalysisTaskEmcalJetV0Filter &);
  AliAnalysisTaskEmcalJetV0Filter& operator=(const AliAnalysisTaskEmcalJetV0Filter &);

  void CreateUserOutputHistograms();

  TString fNameJetRD02;  //
  TString fNameJetRD03;  //
  TString fNameJetRD04;  //

  TString fNameJetMC02;  //
  TString fNameJetMC03;  //
  TString fNameJetMC04;  //

  Bool_t fIsAnaPicoV0;  //
  Bool_t fAnaPicoV0MC;  //


  AliJetContainer      *fJetsContRD02;          //!
  AliParticleContainer *fTracksContRD02;        //!
  AliClusterContainer  *fCaloClustersContRD02;  //!

  AliJetContainer      *fJetsContRD03;          //!
  AliParticleContainer *fTracksContRD03;        //!
  AliClusterContainer  *fCaloClustersContRD03;  //!

  AliJetContainer      *fJetsContRD04;          //!
  AliParticleContainer *fTracksContRD04;        //!
  AliClusterContainer  *fCaloClustersContRD04;  //!

  AliJetContainer      *fJetsContMC02;    //!
  AliParticleContainer *fTracksContMC02;  //!

  AliJetContainer      *fJetsContMC03;    //!
  AliParticleContainer *fTracksContMC03;  //!

  AliJetContainer      *fJetsContMC04;    //!
  AliParticleContainer *fTracksContMC04;  //!

  TClonesArray *fV0s; //!


  AliPicoHeaderCJ *fPicoHeaderCJ; //!

  TClonesArray *fPicoJetsClArrRD02;  //!
  TClonesArray *fPicoJetsClArrRD03;  //!
  TClonesArray *fPicoJetsClArrRD04;  //!

  TClonesArray *fPicoJetsClArrMC02;  //!
  TClonesArray *fPicoJetsClArrMC03;  //!
  TClonesArray *fPicoJetsClArrMC04;  //!

  TClonesArray *fPicoV0sClArr;  //!

  TList *fListUserOutputs;  //!

  ClassDef(AliAnalysisTaskEmcalJetV0Filter, 2);
};

#endif
