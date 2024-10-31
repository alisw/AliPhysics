#ifndef ALIEMCALJETUTILITYSOFTDROP_H
#define ALIEMCALJETUTILITYSOFTDROP_H

#include <TNamed.h>

#include "AliEmcalJetUtility.h"
#include "AliFJWrapper.h"

class AliEmcalJetTask;
class AliEmcalJet;
class AliFJWrapper;
class AliRhoParameter;
class TClonesArray;

class AliEmcalJetUtilitySoftDrop : public AliEmcalJetUtility
{
 public:

  AliEmcalJetUtilitySoftDrop();
  AliEmcalJetUtilitySoftDrop(const char* name);
  AliEmcalJetUtilitySoftDrop(const AliEmcalJetUtilitySoftDrop &jet);
  AliEmcalJetUtilitySoftDrop& operator=(const AliEmcalJetUtilitySoftDrop &jet);
  ~AliEmcalJetUtilitySoftDrop() {;}

  void                   SetRhoName(const char *n)           { fRhoName      = n         ; }
  void                   SetRhomName(const char *n)          { fRhomName     = n         ; }
  void                   SetUseExternalBkg(Bool_t b)         { fUseExternalBkg   = b     ; }

  void                   SetGroomedJetsName(const char *n)       { fGroomedJetsName      = n     ; }
  void                   SetGroomedJetParticlesName(const char *n)  { fGroomedJetParticlesName = n     ; }

  void Init();
  void InitEvent(AliFJWrapper& fjw);
  void Prepare(AliFJWrapper& fjw);
  void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw);
  void Terminate(AliFJWrapper& fjw);
  void SetZCut(Double_t zCut) {fZCut = zCut;}
  void SetBeta(Double_t beta) {fBeta = beta;}
  void SetRecursiveDepth(Int_t n) {fRecursiveDepth = n;}

 protected:

  TString                fGroomedJetsName;                   ///< name of groomed jet collection
  TString                fGroomedJetParticlesName;          ///< name of groomed particle collection
  Bool_t                 fUseExternalBkg;                   ///< use external background for generic subtractor
  TString                fRhoName;                            ///< name of rho
  TString                fRhomName;                           ///< name of rhom
  Double_t               fRho;                                ///< pT background density
  Double_t               fRhom;                               ///< mT background density
    // condition to stop the grooming (rejection of soft splitting) z > fZCut theta^fBeta
  Double_t               fZCut;                              //< fZCut = 0.1                
  Double_t               fBeta;                              //< fBeta = 0
  Int_t                  fRecursiveDepth;                    //< 0: No Soft Drop; 1: Apply once; -1: Apply until the condition stops the grooming

  TClonesArray          *fGroomedJets;                        //!<! groomed jet collection
  TClonesArray          *fGroomedJetParticles;                       //!<! groomed particle collection
  AliRhoParameter       *fRhoParam;                           //!<! event rho
  AliRhoParameter       *fRhomParam;                          //!<! event rhom

  ClassDef(AliEmcalJetUtilitySoftDrop, 2) // Emcal jet utility that implements the constituent subtractor form the fastjet contrib
};
#endif
