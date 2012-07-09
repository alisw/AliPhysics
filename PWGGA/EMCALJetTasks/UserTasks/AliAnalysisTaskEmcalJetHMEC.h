#ifndef AliAnalysisTaskEmcalJetHMEC_H
#define AliAnalysisTaskEmcalJetHMEC_H

// $Id$

class TList;
class TH1;
class TH2;
class THnSparse;
class AliESDEvent;
class AliEventPoolManager;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcalJetHMEC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEmcalJetHMEC();
  AliAnalysisTaskEmcalJetHMEC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHMEC() {}
  
  virtual void            UserCreateOutputObjects();
  virtual Double_t        RelativePhi(Double_t mphi, Double_t vphi);
  virtual void            UserExec(Option_t *option);
  virtual void            Terminate(Option_t *);
  virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
  virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  virtual void            SetTracksName(const char *n)             {fTracksName=n;}
  virtual void            SetJetsName(const char *jn)              {fJetsName=jn;}

  virtual void            SetAreaCut(Double_t a)                   { fAreacut    = a; }
  virtual void            SetTrkBias(Double_t b)                   { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)                  { fClusBias   = b; }  //require a cluster with pt > b in jet

  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
  virtual void            SetEventMixing(Int_t yesno)              { fDoEventMixing=yesno;}
  virtual void            SetMixingTracks(Int_t tracks)            { fMixingTracks = tracks; }





 protected:
  virtual Int_t          GetCentBin(Double_t cent) const;
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetpTjetBin(Double_t pt) const;

  TString                fTracksName;              //name of tracks collection
  TString                fJetsName;                //name of Jet collection
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut
  Double_t               fTrkBias;
  Double_t               fClusBias;
  Int_t                  fDoEventMixing;
  Int_t  		 fMixingTracks;		// size of track buffer for event mixing
  TObjArray*             CloneAndReduceTrackList(TObjArray* tracks);

  AliESDEvent           *fESD;    //! ESD object
  AliEventPoolManager   *fPoolMgr;
  TList                 *fOutputList; //! Output list
  TH1                   *fHistTrackPt; //! Pt spectrum
  TH1                   *fHistCentrality;
  TH2                   *fHistJetEtaPhi;
  TH2                   *fHistTrackEtaPhi;
  TH2                   *fHistJetHEtaPhi;

  TH1                   *fHistJetPt[6];
  TH1                   *fHistJetPtBias[6];
  TH1                   *fHistJetPtTT[6];
  TH2                   *fHistJetH[6][5][3];
  TH2                   *fHistJetHBias[6][5][3];
  TH2                   *fHistJetHTT[6][5][3];
  THnSparse             *fhnMixedEvents;      //!mixed events matrix

 private:
   
  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHMEC, 6); 
};
#endif
