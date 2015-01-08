#ifndef ALIANALYSISTASKRHOMASSBASE_H
#define ALIANALYSISTASKRHOMASSBASE_H

// $Id$

class TString;
class TF1;
class TH1F;
class TH2F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskRhoMassBase : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskRhoMassBase();
  AliAnalysisTaskRhoMassBase(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoMassBase() {}

  void                   UserCreateOutputObjects();

  void                   SetOutRhoMassName(const char *name)                   { fOutRhoMassName           = name ; 
                                                                                 fOutRhoMassScaledName     = Form("%s_Scaled",name) ; }
  void                   SetCompareRhoMassName(const char *name)               { fCompareRhoMassName       = name ;                   }
  void                   SetCompareRhoMassScaledName(const char *name)         { fCompareRhoMassScaledName = name ;                   }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf   ;                       }
  void                   SetRhoMassFunction(TF1* rf)                           { fRhoMassFunction      = rf   ;                       }
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a    ;                       }

  const TString&         GetOutRhoMassName() const                             { return fOutRhoMassName;                              }
  const TString&         GetOutRhoMassScaledName() const                       { return fOutRhoMassScaledName;                        } 

 protected:
  void                   ExecOnce();
  Bool_t                 Run();
  Bool_t                 FillHistograms();

  virtual Double_t       GetRhoMassFactor(Double_t cent);
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fOutRhoMassName;                // name of output rho mass object
  TString                fOutRhoMassScaledName;          // name of output scaled rho mass object
  TString                fCompareRhoMassName;            // name of rho mass object to compare
  TString                fCompareRhoMassScaledName;      // name of scaled rho mass object to compare
  TF1                   *fRhoMassFunction;               // pre-computed rho mass as a function of centrality
  TF1                   *fScaleFunction;                 // pre-computed scale factor as a function of centrality
  Bool_t                 fAttachToEvent;                 // whether or not attach rho mass to the event objects list

  AliRhoParameter       *fOutRhoMass;                    //!output rho object
  AliRhoParameter       *fOutRhoMassScaled;              //!output scaled rho object
  AliRhoParameter       *fCompareRhoMass;                //!rho object to compare
  AliRhoParameter       *fCompareRhoMassScaled;          //!scaled rho object to compare

  TH2F                  *fHistJetMassvsCent;             //!jet mass vs. centrality
  TH2F                  *fHistRhoMassvsCent;             //!rho mass vs. centrality
  TH2F                  *fHistRhoMassScaledvsCent;       //!rho mass scaled vs. centrality
  TH2F                  *fHistDeltaRhoMassvsCent;        //!delta rho mass vs. centrality
  TH2F                  *fHistDeltaRhoMassScalevsCent;   //!delta rho mass scaled vs. centrality

  TH2F                  *fHistRhoMassvsNtrack;           //!rho mass vs. no. of tracks
  TH2F                  *fHistRhoMassScaledvsNtrack;     //!rho mass scaled vs. no. of tracks
  TH2F                  *fHistDeltaRhoMassvsNtrack;      //!delta rho mass vs. no. of tracks
  TH2F                  *fHistDeltaRhoMassScalevsNtrack; //!delta rho mass scaled vs. no. of tracks
 
  TH2F                  *fHistRhoMassvsNcluster;         //!rho mass vs. no. of clusters
  TH2F                  *fHistRhoMassScaledvsNcluster;   //!rho mass scaled vs. no. of clusters

  TH2F                  *fHistGammaVsNtrack;             //!Gamma(<E>/<M>) vs Ntrack

  AliAnalysisTaskRhoMassBase(const AliAnalysisTaskRhoMassBase&);             // not implemented
  AliAnalysisTaskRhoMassBase& operator=(const AliAnalysisTaskRhoMassBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMassBase, 2); // Rho mass base task
};
#endif
