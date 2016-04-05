#ifndef AliAnalysisTaskEmcalJetHMEC_H
#define AliAnalysisTaskEmcalJetHMEC_H

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


  // Jet bias - setters
  virtual void            SetTrackBias(Double_t b)                   { fTrackBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusterBias(Double_t b)                  { fClusterBias = b; }  //require a cluster with pt > b in jet
  // Event trigger/mixed selection - setters
  virtual void            SetTriggerType(UInt_t te)                  { fTriggerType = te; }
  virtual void            SetMixedEventTriggerType(UInt_t me)        { fMixingEventType = me; }
  // Mixed events
  virtual void            SetEventMixing(Bool_t enable)              { fDoEventMixing = enable;}
  virtual void            SetNumberOfMixingTracks(Int_t tracks)      { fNMixingTracks = tracks; }
  virtual void            SetMinNTracksForMixedEvents(Int_t nmt)           { fMinNTracksMixedEvents = nmt; }
  virtual void            SetMinNEventsForMixedEvents(Int_t nme)           { fMinNEventsMixedEvents = nme; }
  virtual void            SetNCentBinsMixedEvent(Bool_t centbins)    { fNCentBinsMixedEvent = centbins; }
  // Switch to cut out some unneeded sparse axis
  void                    SetDoLessSparseAxes(Bool_t dlsa)           { fDoLessSparseAxes = dlsa; }
  void                    SetDoWiderTrackBin(Bool_t wtrbin)          { fDoWiderTrackBin = wtrbin; }
  // set efficiency correction
  void                    SetDoEffCorr(Int_t effcorr)                { fDoEffCorrection = effcorr; }
  virtual void            SetEffCorrFunc(Double_t efffunc)           { fEffFunctionCorrection = efffunc; }
  // Set embedding correction
  void                    SetEmbeddingCorrectionHist(TH2F * embeddingCorrectionHist) { fEmbeddingCorrectionHist = embeddingCorrectionHist; }

  virtual void            UserCreateOutputObjects();
  virtual void            Terminate(Option_t *);

 protected:

  // NOTE: This is not an ideal way to resolve the size of histogram initialization.
  //       Will be resolved when we move fully to the THnSparse
  enum binArrayLimits_t {
    kMaxJetPtBins = 5,
    kMaxTrackPtBins = 7,
    kMaxCentralityBins = 6,
    kMixedEventMulitplictyBins = 8,
    kMaxEtaBins = 3
  };

  // EMCal framework functions
  void                   ExecOnce();
  Bool_t                 Run();

  // Reduce event mixing memory usage
  TObjArray*             CloneAndReduceTrackList();
  // Histogram helper functions
  virtual THnSparse*     NewTHnSparseF(const char* name, UInt_t entries);
  virtual void           GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  // Binning helper functions
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetTrackPtBin(Double_t pt) const;
  virtual Int_t          GetJetPtBin(Double_t pt) const;
  // Helper functions
  void                   InitializeArraysToZero();
  void                   GetDeltaEtaDeltaPhiDeltaR(AliVParticle * particleOne, AliVParticle * particleTwo, Double_t & deltaPhi, Double_t & deltaEta, Double_t & deltaR);
  // Test for biased jet
  Bool_t                 BiasedJet(AliEmcalJet * jet);
  // Corrections
  virtual Double_t       EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
  // Fill methods which allow for the embedding correction
  void                   FillHist(TH1 * hist, Double_t fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   accessSetOfYBinValues(TH2F * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor = -1.0);
  
  // Jet bias
  Double_t               fTrackBias;               ///< Jet track bias
  Double_t               fClusterBias;             ///< Jet cluster bias
  // Event Mixing
  Bool_t                 fDoEventMixing;           ///< flag to do evt mixing
  Int_t                  fNMixingTracks;           ///< size of track buffer for event mixing
  Int_t                  fMinNTracksMixedEvents;   ///< threshold to use event pool # tracks
  Int_t                  fMinNEventsMixedEvents;   ///< threshold to use event pool # events
  UInt_t                 fNCentBinsMixedEvent;     ///< N cent bins for the event mixing pool
  AliEventPoolManager   *fPoolMgr;                 //!<! Event pool manager
  // Event selection types
  UInt_t                 fTriggerType;             ///<
  UInt_t                 fMixingEventType;         ///<
  // Efficiency correction
  Int_t                  fDoEffCorrection;         ///<
  Double_t               fEffFunctionCorrection;   ///<
  // Embedding correction
  TH2F                  *fEmbeddingCorrectionHist; //!<!
  // Histogram binning variables
  Bool_t                 fDoLessSparseAxes;        ///<
  Bool_t                 fDoWiderTrackBin;         ///<

  // TODO: Consider moving to THistManager
  TH1                   *fHistTrackPt;             //!<! Pt spectrum
  TH2                   *fHistJetEtaPhi;           //!<!
  TH2                   *fHistTrackEtaPhi[7];      //!<!
  TH2                   *fHistJetHEtaPhi;          //!<!

  TH3                   *fHistClusEtaPhiEn;        //!<!

  TH1                   *fHistJetPt[6];            //!<!
  TH1                   *fHistJetPtBias[6];        //!<!
  TH1                   *fHistLeadJetPtBias[6];    //!<!
  TH2                   *fHistJetH[6][5][3];       //!<!
  TH2                   *fHistJetHBias[6][5][3];   //!<!
  TH3                   *fHistJHPsi;               //!<! Psi angle distribution
  THnSparse             *fhnMixedEvents;           //!<! Mixed events THnSparse
  THnSparse             *fhnJH;                    //!<! JetH THnSparse

 private:

  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetHMEC, 11); 
};
#endif
