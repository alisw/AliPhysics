#ifndef AliAnalysisTaskEmcalJetHMEC_cxx
#define AliAnalysisTaskEmcalJetHMEC_cxx


class TList;
class TH1;
class TH2;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"
//#include "/project/projectdirs/alice/tschuste/AliceSoftware/aliroot/train_jet/ANALYSIS/AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcalJetHMEC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEmcalJetHMEC();
  AliAnalysisTaskEmcalJetHMEC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHMEC() {}
  
  virtual void      UserCreateOutputObjects();
  virtual Double_t   RelativePhi(Double_t mphi,Double_t vphi);
  virtual void      UserExec(Option_t *option);
  virtual void      Terminate(Option_t *);
  virtual void      SetTracksName(const char *n) {fTracksName=n;}
  virtual void      SetJetsName(const char *jn) {fJetsName=jn;}

  virtual void           SetAreaCut(Double_t a)                   { fAreacut    = a; }
  virtual void           SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void           SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }

 protected:
  virtual Int_t          GetCentBin(Double_t cent) const;
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetpTjetBin(Double_t pt) const;

 private:
  TString      fTracksName;  //name of tracks collection
  TString      fJetsName;  //name of Jet collection
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1        *fHistTrackPt; //! Pt spectrum
  TH1         *fHistCentrality;
  TH2         *fHistJetEtaPhi;
  TH2         *fHistTrackEtaPhi;
  TH2         *fHistJetHEtaPhi;
  TH1         *fHistJetPt[6];
  TH1         *fHistJetPtBias[6];
  TH1         *fHistJetPtTT[6];
  TH2         *fHistJetH[6][5][3];
  TH2         *fHistJetHBias[6][5][3];
  TH2         *fHistJetHTT[6][5][3];

   
  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHMEC, 2); 
};

#endif
