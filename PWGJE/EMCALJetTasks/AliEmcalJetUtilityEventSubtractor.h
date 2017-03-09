#ifndef ALIEMCALJETUTILITYEVENTSUBTRACTOR_H
#define ALIEMCALJETUTILITYEVENTSUBTRACTOR_H

#include <TNamed.h>

#include "AliEmcalJetUtility.h"
#include "AliFJWrapper.h"

class TClonesArray;
class AliEmcalJetTask;
class AliEmcalJet;
class AliFJWrapper;
class AliRhoParameter;

class AliEmcalJetUtilityEventSubtractor : public AliEmcalJetUtility
{
 public:

  AliEmcalJetUtilityEventSubtractor();
  AliEmcalJetUtilityEventSubtractor(const char* name);
  AliEmcalJetUtilityEventSubtractor(const AliEmcalJetUtilityEventSubtractor &jet);
  AliEmcalJetUtilityEventSubtractor& operator=(const AliEmcalJetUtilityEventSubtractor &jet);
  ~AliEmcalJetUtilityEventSubtractor() {;}

  void                   SetRhoName(const char *n)           { fRhoName      = n         ; }
  void                   SetRhomName(const char *n)          { fRhomName     = n         ; }
  void                   SetUseExternalBkg(Bool_t b)         { fUseExternalBkg   = b     ; }
  void                   SetMaxDelR(Double_t r)                { fMaxDelR      = r         ; }

  void                   SetJetsSubName(const char *n)       { fJetsSubName      = n     ; }
  void                   SetParticlesSubName(const char *n)  { fParticlesSubName = n     ; }

  void Init();
  void InitEvent(AliFJWrapper& fjw);
  void Prepare(AliFJWrapper& fjw);
  void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw);
  void Terminate(AliFJWrapper& fjw);

 protected:

  TString                fJetsSubName;                        // name of subtracted jet collection
  TString                fParticlesSubName;                   // name of subtracted particle collection
  Bool_t                 fUseExternalBkg;                     // use external background for generic subtractor
  TString                fRhoName;                            // name of rho
  TString                fRhomName;                           // name of rhom
  Double_t               fRho;                                // pT background density
  Double_t               fRhom;                               // mT background density
  Double_t               fMaxDelR;

  TClonesArray          *fJetsSub;                            //!subtracted jet collection
  TClonesArray          *fParticlesSub;                       //!subtracted particle collection
  AliRhoParameter       *fRhoParam;                           //!event rho
  AliRhoParameter       *fRhomParam;                          //!event rhom

  ClassDef(AliEmcalJetUtilityEventSubtractor, 1) // Emcal jet utility that implements the constituent subtractor form the fastjet contrib
};
#endif
