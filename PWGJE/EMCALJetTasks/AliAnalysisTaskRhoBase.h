#ifndef ALIANALYSISTASKRHOBASE_H
#define ALIANALYSISTASKRHOBASE_H

// $Id$

class TString;
class TF1;
class TH1F;
class TH2F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskRhoBase : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskRhoBase();
  AliAnalysisTaskRhoBase(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoBase() {}

  void                   UserCreateOutputObjects();

  void                   SetRhoName(const char *name)                          { fRhoName              = name ; 
                                                                                 fRhoScaledName        = Form("%s_Scaled",name) ; }
  void                   SetCompareRhoName(const char *name)                   { fCompareRhoName       = name ;                   }
  void                   SetCompareRhoScaledName(const char *name)             { fCompareRhoScaledName = name ;                   }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf   ;                   }
  void                   SetRhoFunction(TF1* rf)                               { fRhoFunction          = rf   ;                   }
  void                   SetInEventSigmaRho(Double_t s)                        { fInEventSigmaRho      = s    ;                   }
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a    ;                   }

 protected:
  void                   ExecOnce();
  Bool_t                 Run();
  Bool_t                 FillHistograms();

  virtual Double_t       GetRhoFactor(Double_t cent);
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fRhoScaledName;                 // name of scaled rho object
  TString                fCompareRhoName;                // name of rho object to compare
  TString                fCompareRhoScaledName;          // name of scaled rho object to compare
  TF1                   *fRhoFunction;                   // pre-computed rho as a function of centrality
  TF1                   *fScaleFunction;                 // pre-computed scale factor as a function of centrality
  Double_t               fInEventSigmaRho;               // in-event sigma rho
  Bool_t                 fAttachToEvent;                 // whether or not attach rho to the event objects list

  AliRhoParameter       *fRhoScaled;                     //!scaled rho object
  AliRhoParameter       *fCompareRho;                    //!rho object to compare
  AliRhoParameter       *fCompareRhoScaled;              //!scaled rho object to compare

  TH2F                  *fHistJetPtvsCent;               //!jet pt vs. centrality
  TH2F                  *fHistJetAreavsCent;             //!jet area vs. centrality
  TH2F                  *fHistJetRhovsCent;              //!jet pt/area vs. centrality
  TH2F                  *fHistNjetvsCent;                //!no. of jets vs. centrality
  TH2F                  *fHistJetPtvsNtrack;             //!jet pt vs. no. of tracks
  TH2F                  *fHistJetAreavsNtrack;           //!jet area vs. no. of tracks
  TH2F                  *fHistNjetvsNtrack;              //!no. of jets vs. no. of tracks
  TH2F                  *fHistNjUEoverNjVsNj[12];        //!ratio no. of jets below rho*A+sigma_rho over. no. of jets vs. no. of jets
  TH2F                  *fHistJetNconstVsPt[4];          //!jet no. of constituents vs. pt
  TH2F                  *fHistJetRhovsEta[4];               //!rho vs. eta
  TH2F                  *fHistRhovsCent;                 //!rho vs. centrality
  TH2F                  *fHistRhoScaledvsCent;           //!rhoscaled vs. centrality
  TH2F                  *fHistDeltaRhovsCent;            //!delta rho vs. centrality
  TH2F                  *fHistDeltaRhoScalevsCent;       //!delta rhoscaled vs. centrality

  TH2F                  *fHistRhovsNtrack;               //!rho vs. no. of tracks
  TH2F                  *fHistRhoScaledvsNtrack;         //!rhoscaled vs. no. of tracks
  TH2F                  *fHistDeltaRhovsNtrack;          //!delta rho vs. no. of tracks
  TH2F                  *fHistDeltaRhoScalevsNtrack;     //!delta rho scaled vs. no. of tracks
 
  TH2F                  *fHistRhovsNcluster;             //!rho vs. no. of clusters
  TH2F                  *fHistRhoScaledvsNcluster;       //!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoBase(const AliAnalysisTaskRhoBase&);             // not implemented
  AliAnalysisTaskRhoBase& operator=(const AliAnalysisTaskRhoBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoBase, 9); // Rho base task
};
#endif
