#ifndef AliAnalysisTaskEmcalJetHMEC_H
#define AliAnalysisTaskEmcalJetHMEC_H

// $Id$

//class TClonesArray;
//class TList;
class TH1;
class TH2;
class TH3;
class THnSparse;
class AliEmcalJet;
class AliEventPoolManager;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetHMEC : public AliAnalysisTaskEmcalJet {
 public:
  enum jetBias_t {
    kDisableBias = 10000
  };

  AliAnalysisTaskEmcalJetHMEC();
  AliAnalysisTaskEmcalJetHMEC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHMEC() {}

  // TODO: Cleanup and organize arguments
  // TODO: Move functions to protected unless good reason not to be!
  // TODO: Consider moving to THistManager
  virtual void            UserCreateOutputObjects();
  //  virtual void            UserExec(Option_t *option);
  virtual void            Terminate(Option_t *);
  virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
  virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  virtual void            SetTrkBias(Double_t b)                   { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)                  { fClusBias   = b; }  //require a cluster with pt > b in jet

  virtual void            SetEventMixing(Int_t yesno)              { fDoEventMixing=yesno;}
  virtual void            SetMixingTracks(Int_t tracks)            { fMixingTracks = tracks; }

  // event trigger/mixed selection - setters
  virtual void            SetTrigType(UInt_t te)       { fTriggerEventType = te; }
  virtual void            SetMixType(UInt_t me)        { fMixingEventType = me; }
  virtual void            SetNMixedTracks(Int_t nmt)   { fNMIXtracks = nmt; }
  virtual void            SetNMixedEvents(Int_t nme)   { fNMIXevents = nme; }

  // switch to cut out some unneeded sparse axis
  void                    SetDoLessSparseAxes(Bool_t dlsa) { fDoLessSparseAxes = dlsa; }
  void                    SetDoWiderTrackBin(Bool_t wtrbin) { fDoWiderTrackBin = wtrbin; }

  virtual void            SetCentBinSize(Bool_t centbins) { fCentBinSize = centbins; }
  // set efficiency correction
  void                    SetDoEffCorr(Int_t effcorr)          { fDoEffCorrection = effcorr; }
  virtual void            SetEffCorrFunc(Double_t efffunc)     { fEffFunctionCorrection = efffunc; }

  // Set embedding correction
  void                    SetEmbeddingCorrectionHist(TH2F * embeddingCorrectionHist) { fEmbeddingCorrectionHist = embeddingCorrectionHist; }

 protected:
  enum binArrayLimits_t {
    kMaxJetPtBins = 5,
    kMaxTrackPtBins = 7,
    kMaxCentralityBins = 6,
    kMaxEtaBins = 3
  };

  void                   ExecOnce();
  Bool_t                 Run();
  void                   InitializeArraysToZero();
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetTrackPtBin(Double_t pt) const;
  virtual Int_t          GetJetPtBin(Double_t pt) const;
  virtual Double_t       EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function

  void GetDeltaEtaDeltaPhiDeltaR(AliVParticle * particleOne, AliVParticle * particleTwo, Double_t & deltaPhi, Double_t & deltaEta, Double_t & deltaR);
  Bool_t BiasedJet(AliEmcalJet * jet);

  // Fill methods which allow for the embedding correction
  void                   FillHist(TH1 * hist, Double_t fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   accessSetOfYBinValues(TH2F * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor = -1.0);

  Double_t               fTrkBias;
  Double_t               fClusBias;
  Int_t                  fDoEventMixing;           // flag to do evt mixing
  Int_t                  fMixingTracks;            // size of track buffer for event mixing
  Int_t                  fNMIXtracks;              // threshold to use event pool # tracks
  Int_t                  fNMIXevents;              // threshold to use event pool # events
  TObjArray*             CloneAndReduceTrackList();

  // event selection types
  UInt_t         fTriggerEventType;
  UInt_t         fMixingEventType;

  // efficiency correction
  Int_t          fDoEffCorrection;
  Double_t       fEffFunctionCorrection;

  // Embedding correction
  TH2F                   *fEmbeddingCorrectionHist;

  Bool_t         fDoLessSparseAxes;
  Bool_t         fDoWiderTrackBin;

  UInt_t         fCentBinSize;

  AliEventPoolManager   *fPoolMgr; //!
  TH1                   *fHistTrackPt; //! Pt spectrum
  TH2                   *fHistJetEtaPhi;//!
  TH2                   *fHistTrackEtaPhi[7];//!
  TH2                   *fHistJetHEtaPhi;//!

  TH3                   *fHistClusEtaPhiEn; //!

  TH1                   *fHistJetPt[6]; //!
  TH1                   *fHistJetPtBias[6];//!
  TH1                   *fHistLeadJetPtBias[6];//!
  TH2                   *fHistJetH[6][5][3];//!
  TH2                   *fHistJetHBias[6][5][3];//!
  THnSparse             *fhnMixedEvents;      //!mixed events matrix
  THnSparse             *fhnJH;      //!Fg events matrix
  TH3                   *fHistJHPsi; //! Psi angle distribution

 private:

  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetHMEC, 10); 
};
#endif
