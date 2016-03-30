#ifndef AliAnalysisTaskEmcalJetHMEC_H
#define AliAnalysisTaskEmcalJetHMEC_H

// $Id$

class TClonesArray;
class TList;
class TH1;
class TH2;
class TH3;
class THnSparse;
class AliEmcalJet;
class AliESDEvent;
class AliAODEvent;
class AliEventPoolManager;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetHMEC : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetHMEC();
  AliAnalysisTaskEmcalJetHMEC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHMEC() {}
  
  virtual void            UserCreateOutputObjects();
  virtual Double_t        RelativePhi(Double_t mphi, Double_t vphi);
//  virtual void            UserExec(Option_t *option);
  virtual void            Terminate(Option_t *);
  virtual Int_t           AcceptthisJet(AliEmcalJet *jet);
  virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
  virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  /*virtual void            SetTracksName(const char *n)             {fTracksName=n;}
  virtual void            SetJetsName(const char *jn)              {fJetsName=jn;}
  virtual void            SetCaloClustersName(const char *cn)      {fCaloClustersName=cn;}*/

  //virtual void            SetAreaCut(Double_t a)                   { fAreacut    = a; }
  virtual void            SetTrkBias(Double_t b)                   { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)                  { fClusBias   = b; }  //require a cluster with pt > b in jet

  virtual void            SetTrkEta(Double_t e)                    { fTrkEta   = e; }  //eta range of the associated tracks

  /*virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }*/
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
  void					 ExecOnce();
  Bool_t			     Run();
  //virtual Int_t          GetCentBin(Double_t cent) const;
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetpTjetBin(Double_t pt) const;
  virtual Double_t       EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function

  // Fill methods which allow for the embedding correction
  void                   FillHist(TH1 * hist, Double_t fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight = 1.0, Bool_t noCorrection = kTRUE);
  void                   accessSetOfYBinValues(TH2F * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor = -1.0);

  /*TString                fTracksName;              // name of tracks collection
  TString                fJetsName;                // name of Jet collection
  TString                fCaloClustersName;        // name of Calo Cluster collection
  */
  /*Double_t               fPhimin;                  // phi min of jet
  Double_t               fPhimax;                  // phi max of jet
  Double_t               fEtamin;                  // eta min of jet
  Double_t               fEtamax;                  // eta max of jet
  Double_t               fAreacut;                 // area cut of jet
  */
  Double_t               fTrkBias;
  Double_t               fClusBias;
  Double_t               fTrkEta;                  // eta min/max of tracks
  Int_t                  fDoEventMixing;           // flag to do evt mixing
  Int_t  		 fMixingTracks;		       // size of track buffer for event mixing
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

  //AliESDEvent           *fESD;    //! ESD object
  //AliAODEvent			*fAOD;    //! AOD object
  AliEventPoolManager   *fPoolMgr; //!
  TH1                   *fHistTrackPt; //! Pt spectrum
  //TH1                   *fHistCentrality;//!
  TH2                   *fHistJetEtaPhi;//!
  TH2                   *fHistTrackEtaPhi[7];//!
  TH2                   *fHistJetHEtaPhi;//!

  TH3                   *fHistClusEtaPhiEn; //!

  TH1                   *fHistJetPt[6]; //!
  TH1                   *fHistJetPtBias[6];//!
  TH1                   *fHistLeadJetPtBias[6];//!
  TH1                   *fHistJetPtTT[6];//!
  TH2                   *fHistJetH[6][5][3];//!
  TH2                   *fHistJetHBias[6][5][3];//!
  TH2                   *fHistJetHTT[6][5][3];//!
  THnSparse             *fhnMixedEvents;      //!mixed events matrix
  THnSparse             *fhnJH;      //!Fg events matrix
  TH3                   *fHistJHPsi; //! Psi angle distribution

 private:
   
  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHMEC, 10); 
};
#endif
