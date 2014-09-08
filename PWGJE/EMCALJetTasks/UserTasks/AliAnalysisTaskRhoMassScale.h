#ifndef ALIANALYSISTASKRHOMASSSCALE_H
#define ALIANALYSISTASKRHOMASSSCALE_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskRhoMassScale : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskRhoMassScale();
  AliAnalysisTaskRhoMassScale(const char *name);
  virtual ~AliAnalysisTaskRhoMassScale();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerNeutral(Int_t c)                { fContainerNeutral     = c   ; }
  void SetJetContainerCharged(Int_t c)                { fContainerCharged     = c   ; }

  void SetRhoMNeutralName(const char *n)              { fRhoMNeutralName = n ; }
  void SetRhoMChargedEmcalName(const char *n)         { fRhoMChargedEmcalName = n ; }
  void SetRhoMCharged2xEmcalName(const char *n)       { fRhoMCharged2xEmcalName = n ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

 private:
  Int_t                               fContainerNeutral;              // particle level jets
  Int_t                               fContainerCharged;              // detector level jets
  TString                             fRhoMNeutralName;               // Name of neutral rho mass object
  TString                             fRhoMChargedEmcalName;          // Name of charged rho mass object in EMCal acceptance
  TString                             fRhoMCharged2xEmcalName;        // Name of charged rho mass object in two times EMCal acceptance
  AliRhoParameter                    *fRhoMNeutral;                   //!neutral rho_m
  AliRhoParameter                    *fRhoMChargedEmcal;              //!charged rho_m in EMCal acceptance
  AliRhoParameter                    *fRhoMCharged2xEmcal;            //!charged rho_m in two times EMCal acceptance

  TH2                                *fHistScaleEmcalvsCent;          //!scale factor 1xEmcal vs centrality
  TH2                                *fHistScale2EmcalvsCent;         //!scale factor 2xEmcal vs centrality
  TH2                                *fHistDeltaScale2EmcalvsCent;    //!difference between scale factors vs centrality

  TH2                                *fHistScaleEmcalvsMult;          //!scale factor 1xEmcal vs track multiplicity
  TH2                                *fHistScale2EmcalvsMult;         //!scale factor 2xEmcal vs track multiplicity
  TH2                                *fHistDeltaScale2EmcalvsMult;    //!difference between scale factors vs track multiplicity

  AliAnalysisTaskRhoMassScale(const AliAnalysisTaskRhoMassScale&);            // not implemented
  AliAnalysisTaskRhoMassScale &operator=(const AliAnalysisTaskRhoMassScale&); // not implemented

  ClassDef(AliAnalysisTaskRhoMassScale, 1)
};
#endif

