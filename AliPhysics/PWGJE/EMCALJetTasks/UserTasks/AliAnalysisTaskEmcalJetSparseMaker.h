#ifndef ALIANALYSISTASKEMCALJETSPARSEMAKER_H
#define ALIANALYSISTASKEMCALJETSPARSEMAKER_H

#include <TString.h>

#include "AliAnalysisTaskEmcalJet.h"
//=============================================================================

class THnSparse;

class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
//=============================================================================

class AliAnalysisTaskEmcalJetSparseMaker : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalJetSparseMaker();
  AliAnalysisTaskEmcalJetSparseMaker(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalJetSparseMaker();
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

  AliAnalysisTaskEmcalJetSparseMaker(const AliAnalysisTaskEmcalJetSparseMaker &);
  AliAnalysisTaskEmcalJetSparseMaker &operator=(const AliAnalysisTaskEmcalJetSparseMaker &);
//=============================================================================

  void MakeSparse();
  void FillSparse();

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
  THnSparse *fHnsJets; //!

  TList *fListOutputEvH; //!
  TList *fListOutputJet; //!
//=============================================================================

  ClassDef(AliAnalysisTaskEmcalJetSparseMaker, 1)
};

#endif
