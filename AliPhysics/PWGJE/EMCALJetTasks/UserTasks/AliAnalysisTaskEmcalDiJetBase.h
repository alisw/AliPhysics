#ifndef ALIANALYSISTASKEMCALDIJETBASE_H
#define ALIANALYSISTASKEMCALDIJETBASE_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalDiJetBase : public AliAnalysisTaskEmcalJet {
 public:
  enum JetFullChargedMatchingType {
    kFraction   = 0,     // match full and charged jets with largest shared charged pt fraction
    kGeo        = 1,     // match full and charged jets geometrically
    kNoMatching = 3      // include autocorrelation in dijet correlation
  };

  enum JetCorrelationType {
    kCorrelateAll = 0,   // correlate all jets with all jets in event
    kCorrelateTwo = 1,   // correlate all jets with leading jet in opposite hemisphere
    kCorrelateLS  = 2    // correlate leading and subleading jet
  };

  AliAnalysisTaskEmcalDiJetBase();
  AliAnalysisTaskEmcalDiJetBase(const char *name);
  virtual ~AliAnalysisTaskEmcalDiJetBase();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  Bool_t                      SelectEvent();              //decides if event is used for analysis

  //Setters
  void SetDebug(Int_t d)                                        { fDebug = d;}

  void SetJetCorrelationType(JetCorrelationType c)              { fJetCorrelationType = c; }

  void SetFullChargedMatchingType(JetFullChargedMatchingType m) { fJetFullChargedMatchingType = m; }

  void SetTriggerClass(const char *n)       { fTriggerClass = n; }

  void SetContainerFull(Int_t c)            { fContainerFull      = c;}
  void SetContainerCharged(Int_t c)         { fContainerCharged   = c;}
  void SetContainerFullMC(Int_t c)          { fContainerFullMC    = c;}
  void SetContainerChargedMC(Int_t c)       { fContainerChargedMC = c;} 

  void SetRhoType(Int_t i)                  { fRhoType = i;}

  void SetDoChargedCharged(Bool_t b)        { fDoChargedCharged = b;}
  void SetDoFullCharged(Bool_t b)           { fDoFullCharged    = b;}
  void SetDoFullFull(Bool_t b)              { fDoFullFull       = b;}

  void SetPtMinTriggerJet(Double_t ptmin)   { fPtMinTriggerJet  = ptmin;}
  void SetDoPtBias(Bool_t b)                { fDoPtBias         = b;}

  void SetMinSharedFraction(Double_t f)     { fMinFractionShared = f;}

  void ResetMatchFlag()                     { fMatchingDone = kFALSE; }

  //Getters
  Double_t GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2);
  Double_t GetDeltaPhi(Double_t phi1,Double_t phi2);
  Double_t GetDeltaR(const AliEmcalJet* jet1, const AliEmcalJet* jet2) const;

  Double_t GetZ(const AliVParticle *trk, const AliEmcalJet *jet)       const;
  Double_t GetZ(Double_t trkPx, Double_t trkPy, Double_t trkPz, Double_t jetPx, Double_t jetPy, Double_t jetPz) const;

  AliEmcalJet* GetLeadingJetOppositeHemisphere(Int_t type, Int_t typea, const AliEmcalJet *jetTrig);
  AliEmcalJet* GetSecondLeadingJetOppositeHemisphere(Int_t type, Int_t typea, const AliEmcalJet *jetTrig);

 protected:
  virtual Bool_t                      RetrieveEventObjects();

  Bool_t                      IsSameJet(Int_t jt, Int_t ja, Int_t type, Bool_t isMC = kFALSE);
  Double_t                    GetJetPt(const AliEmcalJet *jet, Int_t type);

  void                        MatchJetsGeo(Int_t cFull, Int_t cCharged,
					   Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t type = 0);
  Double_t                    GetFractionSharedPt(const AliEmcalJet *jetFull, const AliEmcalJet *jetCharged) const;

  void                        SetChargedFractionIndex();
  void                        SetChargedFractionIndexMC();

  Bool_t                     fDebug;                      // debug level
  JetCorrelationType         fJetCorrelationType;         // type of correlation between jets
  JetFullChargedMatchingType fJetFullChargedMatchingType; //matching type between full and charged jets to be used
  TString                    fTriggerClass;               // trigger class to analyze EJ1 or EJ2    

  Int_t             fContainerCharged;          //  number of container with charged jets DET
  Int_t             fContainerFull;             //  number of container with full jets DET
  Int_t             fContainerChargedMC;        //  number of container with charged jets MC
  Int_t             fContainerFullMC;           //  number of container with full jets MC

  Int_t             fRhoType;                   //  rho type
  Double_t          fRhoChVal;                  //  charged rho value
  Double_t          fRhoFullVal;                // scaled charged rho value

  Bool_t            fDoChargedCharged;          //  do charged-charged ana
  Bool_t            fDoFullCharged;             //  do full-charged ana
  Bool_t            fDoFullFull;                //  do full-full ana

  Double_t          fPtMinTriggerJet;           //  minimum pT of trigger jet
  Bool_t            fDoPtBias;                  //  pT trigger jet > pT assoc jet
  Double_t          fMinFractionShared;         //  minimum fraction charged pT

  Bool_t            fMatchingDone;              // flag to indicate if matching is done or not
  TArrayI           faFullFracIndex;            // index of charged jet with largest shared charged fraction - detector level
  TArrayI           faFullFracIndexMC;          // index of charged jet with largest shared charged fraction - particle level

  TH1F             *fhNEvents;                  //! Histo number of events
  TH1              *fHistTrialsSelEvents;       //!trials from pyxsec.root only for selected events

 private:
  AliAnalysisTaskEmcalDiJetBase(const AliAnalysisTaskEmcalDiJetBase&);            // not implemented
  AliAnalysisTaskEmcalDiJetBase &operator=(const AliAnalysisTaskEmcalDiJetBase&); // not implemented

  ClassDef(AliAnalysisTaskEmcalDiJetBase, 6) // dijet base task
};
#endif
