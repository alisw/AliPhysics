#ifndef ALIANALYSISTASKRHO_cxx
#define ALIANALYSISTASKRHO_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;
class TString;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskRho : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskRho();
  AliAnalysisTaskRho(const char *name);
  virtual ~AliAnalysisTaskRho() {}
  
  virtual void          UserCreateOutputObjects();
  virtual void          UserExec(Option_t*);
  virtual void          Terminate(Option_t*);

  void                  SetTracksName(const char *n)                          { fTracksName  = n    ; }
  void                  SetJetsName(const char *n)                            { fJetsName    = n    ; }
  void                  SetRhosName(const char *rn)                           { fRhosName    = rn   ; }
  void                  SetJetPhi(Double_t pmin = 1.8, Double_t pmax = 2.74)  { fPhiMin      = pmin ; fPhiMax = pmax  ; }
  void                  SetJetEta(Double_t emin = -0.3, Double_t emax = 0.3)  { fEtaMin      = emin ; fEtaMax = emax  ; }
  void                  SetAreaCut(Double_t a = 0.0)                          { fAreaCut     = a    ; }
  
 protected:
  virtual void          Sort(vector<Double_t>& v)            ;
  virtual Double_t      GetMedian(vector<Double_t> v, int c) ;
  virtual Double_t      GetScaleFactor(Double_t cent)        ;
  virtual Double_t      GetRhoFactor(Double_t cent)          ;

 private: 
  TString                fTracksName;                    // name of track collection
  TString                fJetsName;                      // name of jet collection
  TString                fRhosName;                      // name of Rho array output
  TString                fClustersName;                  // name of clusters collection

  TList                 *fOutputList;                    //! Output list
  TH1F                  *fHistCentrality;
  TH1F                  *fHistJetPt;
  TH1F                  *fHistJetArea; 
  TH2F                  *fHistRhovsCent;
  TH2F                  *fHistDeltaRhovsCent;  
  TH2F                  *fHistDeltaRhoScalevsCent;  
  TH2F                  *fHistJetPtvsCent;  
  TH2F                  *fHistJetAreavsCent;  
  TH2F                  *fHistNjetvsCent;

  TH2F                 *fHistRhovsNtrack;
  TH2F                 *fHistDeltaRhovsNtrack;   
  TH2F                 *fHistDeltaRhoScalevsNtrack;   
  TH2F                 *fHistJetPtvsNtrack;
  TH2F                 *fHistJetAreavsNtrack;
  TH2F                 *fHistNjetvsNtrack;

  TClonesArray         *fArrRhos;

  Double_t              fPhiMin;
  Double_t              fPhiMax;
  Double_t              fEtaMin;
  Double_t              fEtaMax;
  Double_t              fAreaCut;
  Int_t                 fCswitch;

  AliAnalysisTaskRho(const AliAnalysisTaskRho&);             // not implemented
  AliAnalysisTaskRho& operator=(const AliAnalysisTaskRho&);  // not implemented
  
  ClassDef(AliAnalysisTaskRho, 1); // Rho task
};

#endif
