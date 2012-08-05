#ifndef AliAnalysisTaskEmcalJetSpectra_h
#define AliAnalysisTaskEmcalJetSpectra_h

// $Id$

class TList;
class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcalJetSpectra : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEmcalJetSpectra();
  AliAnalysisTaskEmcalJetSpectra(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectra() {}
  
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);
  virtual void           Terminate(Option_t *);

  virtual void           SetAreaCut(Double_t a)                   { fAreacut    = a; }
  virtual void           SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void           SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
  virtual void           SetJetsName(const char *n)               { fJetsName   = n; }
  virtual void           SetRhos1Name(const char *n)              { fRhos1Name  = n; }
  virtual void           SetRhos2Name(const char *n)              { fRhos2Name  = n; }
  virtual void           SetRhos3Name(const char *n)              { fRhos3Name  = n; }
  virtual void           SetTracksName(const char *n)             { fTracksName = n; }
  
 protected:
  virtual Int_t          GetCentBin(Double_t cent) const;
   
 private:
  TString                fTracksName;              // name of track collection
  TString                fJetsName;                // name of jet collection
  TString                fClustersName;            // name of clusters collection
  TString                fRhos1Name;               // name of Rho1 array output
  TString                fRhos2Name;               // name of Rho2 array output
  TString                fRhos3Name;               // name of Rho2 array output
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut

  AliESDEvent           *fESD;                     //!esd event
  TList                 *fOutputList;              //!output list
  TH1F                  *fHistCentrality;          //!centrality
  TH2F                  *fHistDeltaRho12vsCent;    //!delta rho1 and rho2 vs centrality
  TH2F                  *fHistDeltaRho13vsCent;    //!delta rho1 and rho3 vs centrality
  TH2F                  *fHistDeltaRho23vsCent;    //!delta rho2 and rho3 vs centrality
  TH2F                  *fHistDeltaJetPt12vsCent;  //!delta jet pt1 and pt2 vs centrality 
  TH2F                  *fHistDeltaJetPt13vsCent;  //!delta jet pt1 and pt3 vs centrality 
  TH2F                  *fHistDeltaJetPt23vsCent;  //!delta jet pt2 and pt3 vs centrality 
  TH2F                  *fHistRho1vsCent;          //!rho1 vs centrality
  TH2F                  *fHistRho2vsCent;          //!rho2 vs centrality
  TH2F                  *fHistRho3vsCent;          //!rho3 vs centrality
  TH2F                  *fHistNEFvsPt[6][4];       //!neutral energy fraction vs pt
  TH2F                  *fHistZvsPt[6][4];         //!z all vs pt
  TH2F                  *fHistZchvsPt[6][4];       //!z charged vs pt
  TH2F                  *fHistZemvsPt[6][4];       //!z neutral vs pt
  TH1F                  *fHistJetPt[6][4];         //!jet pt
  TH1F                  *fHistJetPt3[6][4];        //!jet pt>3
  TH1F                  *fHistJetPt5[6][4];        //!jet pt>5
  TH1F                  *fHistJetPt7[6][4];        //!jet pt>7
  TH1F                  *fHistJetPt9[6][4];        //!jet pt>9
  TH2F                  *fHistNconsvsPt[6][4];     //!constituents vs pt
  TH1F                  *fHistRawJetPt[6];         //!raw jet pt
  TH2F                  *fHistAreavsRawPt[6];      //!area vs raw pt

  AliAnalysisTaskEmcalJetSpectra(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  AliAnalysisTaskEmcalJetSpectra& operator=(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetSpectra, 3); // Emcal jet spectra task
};
#endif
