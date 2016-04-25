#ifndef ALIEMCALJETUTILITYGENSUBTRACTOR_H
#define ALIEMCALJETUTILITYGENSUBTRACTOR_H

#include <TNamed.h>

#include "AliEmcalJetUtility.h"
#include "AliFJWrapper.h"

class AliEmcalJetTask;
class AliEmcalJet;
class AliFJWrapper;
class AliRhoParameter;

class AliEmcalJetUtilityGenSubtractor : public AliEmcalJetUtility
{
 public:

  AliEmcalJetUtilityGenSubtractor();
  AliEmcalJetUtilityGenSubtractor(const char* name);
  AliEmcalJetUtilityGenSubtractor(const AliEmcalJetUtilityGenSubtractor &jet);
  AliEmcalJetUtilityGenSubtractor& operator=(const AliEmcalJetUtilityGenSubtractor &jet);
  ~AliEmcalJetUtilityGenSubtractor() {;}

  void Init();
  void Prepare(AliFJWrapper& fjw);
  void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw);
  void Terminate(AliFJWrapper& fjw);

  void                   SetRhoName(const char *n)                                        { fRhoName      = n                      ; }
  void                   SetRhomName(const char *n)                                       { fRhomName     = n                      ; }

  void                   SetGenericSubtractionJetMass(Bool_t b)                           { fDoGenericSubtractionJetMass        = b; }
  void                   SetGenericSubtractionGR(Bool_t b, Double_t rmax = 2., 
                                                 Double_t dr = 0.04, Double_t ptmin = 0.) { fDoGenericSubtractionGR             = b; fRMax = rmax; fDRStep = dr; fPtMinGR = ptmin;}
  void                   SetGenericSubtractionExtraJetShapes(Bool_t b)                    { fDoGenericSubtractionExtraJetShapes = b; }
  void                   SetGenericSubtractionNsubjettiness(Bool_t b)                     { fDoGenericSubtractionNsubjettiness = b; }
  void                   SetUseExternalBkg(Bool_t b)                                      { fUseExternalBkg                     = b; }

 protected:

  Bool_t                 fDoGenericSubtractionJetMass;        // calculate generic subtraction
  Bool_t                 fDoGenericSubtractionGR;             // calculate generic subtraction for angular structure function GR
  Bool_t                 fDoGenericSubtractionExtraJetShapes; // calculate generic subtraction for other jet shapes like radialmoment,pTD etc
  Bool_t                 fDoGenericSubtractionNsubjettiness; // calculate generic subtraction for 1subjettiness, 2subjettiness and the opening Angle between subjets
  Bool_t                 fUseExternalBkg;                     // use external background for generic subtractor
  TString                fRhoName;                            // name of rho
  TString                fRhomName;                           // name of rhom
  Double_t               fRho;                                // pT background density
  Double_t               fRhom;                               // mT background density
  Double_t               fRMax;                               // R max for GR calculation
  Double_t               fDRStep;                             // step width for GR calculation
  Double_t               fPtMinGR;                            // min pT for GR calculation

  AliRhoParameter       *fRhoParam;                           //!event rho
  AliRhoParameter       *fRhomParam;                          //!event rhom

  ClassDef(AliEmcalJetUtilityGenSubtractor, 1) // Emcal jet utility that implements generic subtractors form the fastjet contrib
};
#endif
