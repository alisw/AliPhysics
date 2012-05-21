#ifndef ALIANALYSISTASKRHO_cxx
#define ALIANALYSISTASKRHO_cxx

// $Id: $

class TList;
class TH1F;
class TH2F;
class TClonesArray;
class TString;
class TF1;

#include <TParameter.h>

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRho : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRho();
  AliAnalysisTaskRho(const char *name);
  AliAnalysisTaskRho(const char *name, Bool_t histo);
  virtual ~AliAnalysisTaskRho() {}
  
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t*);
  virtual void           Terminate(Option_t*);

  void                   SetAreaCut(Double_t a = 0.0)                          { fAreaCut       = a    ; }
  void                   SetJetEta(Double_t emin, Double_t emax)               { fEtaMin        = emin ; fEtaMax = emax  ; }
  void                   SetJetPhi(Double_t pmin, Double_t pmax)               { fPhiMin        = pmin ; fPhiMax = pmax  ; }
  void                   SetJetsName(const char *n)                            { fJetsName      = n    ; }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction = sf   ; }
  void                   SetTracksName(const char *n)                          { fTracksName    = n    ; }
  void                   SetExcludeLeadJets(UInt_t n)                          { fNExclLeadJets = n    ; }
  
 protected:
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fTracksName;                    // name of track collection
  TString                fJetsName;                      // name of jet collection
  TString                fClustersName;                  // name of clusters collection
  TString                fRhoScaledName;                 // name of scaled rho object
  Double_t               fPhiMin;                        // minimum phi
  Double_t               fPhiMax;                        // maximum phi
  Double_t               fEtaMin;                        // minimum eta
  Double_t               fEtaMax;                        // maximum eta
  Double_t               fAreaCut;                       // cut on jet area
  UInt_t                 fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation
  TF1                   *fScaleFunction;                 // pre-computed scale factor as a function of centrality
  Bool_t                 fCreateHisto;                   // whether or not create histograms
  TList                 *fOutputList;                    //!output list
  TH1F                  *fHistCentrality;                //!centrality distribution
  TH1F                  *fHistJetPt;                     //!jet pt distribution
  TH1F                  *fHistJetArea;                   //!jet area
  TH2F                  *fHistRhovsCent;                 //!rho vs. centrality
  TH2F                  *fHistDeltaRhovsCent;            //!delta rho vs. centrality
  TH2F                  *fHistDeltaRhoScalevsCent;       //!delta rhoscaled vs. centrality
  TH2F                  *fHistJetPtvsCent;               //!jet pt vs. centrality
  TH2F                  *fHistJetAreavsCent;             //!jet area vs. centrality
  TH2F                  *fHistNjetvsCent;                //!no. of jets vs. centrality
 
  TH2F                  *fHistRhovsNtrack;               //!rho vs. no. of tracks
  TH2F                  *fHistDeltaRhovsNtrack;          //!delta rho vs. no. of tracks
  TH2F                  *fHistDeltaRhoScalevsNtrack;     //!delta rho scaled vs. no. of tracks
  TH2F                  *fHistJetPtvsNtrack;             //!jet pt vs. no. of tracks
  TH2F                  *fHistJetAreavsNtrack;           //!jet area vs. no. of tracks
  TH2F                  *fHistNjetvsNtrack;              //!no. of jets vs. no. of tracks
  TParameter<Double_t>  *fRhoScaled;                     //!per event scaled rho

  AliAnalysisTaskRho(const AliAnalysisTaskRho&);             // not implemented
  AliAnalysisTaskRho& operator=(const AliAnalysisTaskRho&);  // not implemented
  
  ClassDef(AliAnalysisTaskRho, 3); // Rho task
};
#endif
