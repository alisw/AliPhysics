#ifndef ALIANALYSISTASKEMCALTMPSPARSEMAKER_H
#define ALIANALYSISTASKEMCALTMPSPARSEMAKER_H

#include <TString.h>

#include "AliAnalysisTaskEmcalJet.h"
//=============================================================================

class THnSparse;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
//=============================================================================

class AliAnalysisTaskEmcalTmpSparseMaker : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalTmpSparseMaker();
  AliAnalysisTaskEmcalTmpSparseMaker(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalTmpSparseMaker();
//=============================================================================

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t *opt);

  TString GetNameJet() const { return fNameJet; };
  void SetNameJet(const TString s) { fNameJet = s; }
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

  AliAnalysisTaskEmcalTmpSparseMaker(const AliAnalysisTaskEmcalTmpSparseMaker &);
  AliAnalysisTaskEmcalTmpSparseMaker &operator=(const AliAnalysisTaskEmcalTmpSparseMaker &);
//=============================================================================

  void MakeSparseEveH();
  void MakeSparseTrks();
  void MakeSparseClus();
  void MakeSparseJets();

  void FillSparseEveH();
  void FillSparseTrks();
  void FillSparseClus();
  void FillSparseJets();

  Double_t CalcAysPlane();
  Double_t CalcRelPhiEP(Double_t dPhi);
//=============================================================================

  Double_t fMtCh; //
  Double_t fMtEm; //
  Double_t fRhoV; //
  Double_t fAPhi; //
//=============================================================================

  TString fNameJet; //

  AliJetContainer      *fContJets; //!
  AliParticleContainer *fContTrks; //!
  AliClusterContainer  *fContClus; //!

  THnSparse *fHnsEveH; //!
  THnSparse *fHnsTrks; //!
  THnSparse *fHnsClus; //!
  THnSparse *fHnsJets; //!

  TList *fListOutputEvH; //!
  TList *fListOutputTrk; //!
  TList *fListOutputClu; //!
  TList *fListOutputJet; //!
//=============================================================================

  ClassDef(AliAnalysisTaskEmcalTmpSparseMaker, 1)
};

#endif
