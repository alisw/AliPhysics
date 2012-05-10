#ifndef AliAnalysisTaskEmcalJetSpectra_cxx
#define AliAnalysisTaskEmcalJetSpectra_cxx


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
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetTracksName(const char *n)            { fTracksName = n; }
  virtual void   SetJetsName(const char *n)              { fJetsName = n;   }
  virtual void   SetRhos1Name(const char *rn)             { fRhos1Name=rn;    }
  virtual void   SetRhos2Name(const char *rn)             { fRhos2Name=rn;    }
  virtual void   SetJetPhi(Double_t pmin, Double_t pmax) { phimin = pmin; phimax = pmax; }
  virtual void   SetJetEta(Double_t emin, Double_t emax) { etamin = emin; etamax = emax; }
  virtual void SetAreaCut(Double_t a)                    { areacut = a;}
  
 protected:
  virtual Int_t GetCentBin(Double_t cent) const;
 
  
 private:
  TString                fTracksName;             // name of track collection
  TString                fJetsName;             // name of jet collection
  TString                fRhos1Name;             // name of Rho1 array output
  TString                fRhos2Name;             // name of Rho2 array output
  TString                fClustersName;             // name of clusters collection
  TClonesArray *fArrRhos1;
  TClonesArray *fArrRhos2;

  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistCentrality;
  TH1F        *fHistJetArea; 
  TH1F        *fHistJetMaxPt;
  TH1F        *fHistJetZ;
  TH1F        *fHistJetNEF;
  TH2F        *fHistJetPtvsCent;  
  TH2F        *fHistJetPtM3vsCent;  
  TH2F        *fHistLeadingJetPtvsCent;  
  TH2F        *fHistLeadingJetPtM3vsCent;  
  TH2F        *fHistJetAreavsCent;  
  TH2F        *fHistJetMaxPtvsCent;  
  TH2F        *fHistJetZvsCent;  
  TH2F        *fHistJetNEFvsCent;
  TH2F        *fHistNjetvsCent;
  
  TH2F        *fHistJetPtvsNtrack;
  TH2F        *fHistJetAreavsNtrack;
  TH2F        *fHistJetMaxPtvsNtrack;
  TH2F        *fHistJetZvsNtrack;
  TH2F        *fHistJetNEFvsNtrack;
  TH2F        *fHistNjetvsNtrack;

  TH2F       *fHistNEFvsPt[6][3];
  TH2F       *fHistZvsPt[6][3];
  TH2F       *fHistZchvsPt[6][3];
  TH2F       *fHistZemvsPt[6][3];
  TH2F       *fHistAreavsPt[6][3];
  TH1F       *fHistJetPt[6][3];
  TH2F       *fHistNconsvsPt[6][3];
  TH1F       *fHistRawJetPt[6];
  TH2F       *fHistAreavsRawPt[6];
  TH2F       *fHistRho1vsCent;
  TH2F       *fHistRho2vsCent;
  TH2F       *fHistRho3vsCent;

  TH2F       *fHistDeltaRho12vsCent;
  TH2F       *fHistDeltaRho13vsCent;
  TH2F       *fHistDeltaRho23vsCent;

  Double_t phimin;
  Double_t phimax;
  Double_t etamin;
  Double_t etamax;
  Double_t areacut;
  Int_t cSwitch;

   AliAnalysisTaskEmcalJetSpectra(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  AliAnalysisTaskEmcalJetSpectra& operator=(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetSpectra, 1); // example of analysis
};

#endif
