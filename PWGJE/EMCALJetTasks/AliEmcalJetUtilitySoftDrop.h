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
  void Prepare(AliFJWrapper& fjw);
  void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw);
  void Terminate(AliFJWrapper& fjw);

 protected:

  TString                fGroomedJetsName;                   ///< name of groomed jet collection
  TString                fGroomedJetParticlesName;          ///< name of groomed particle collection
  Bool_t                 fUseExternalBkg;                   ///< use external background for generic subtractor
  TString                fRhoName;                            ///< name of rho
  TString                fRhomName;                           ///< name of rhom
  Double_t               fRho;                                ///< pT background density
  Double_t               fRhom;                               ///< mT background density

  TClonesArray          *fGroomedJets;                        //!<! groomed jet collection
  TClonesArray          *fGroomedJetParticles;                       //!<! groomed particle collection
  AliRhoParameter       *fRhoParam;                           //!<! event rho
  AliRhoParameter       *fRhomParam;                          //!<! event rhom

  ClassDef(AliEmcalJetUtilitySoftDrop, 1) // Emcal jet utility that implements the constituent subtractor form the fastjet contrib
};
#endif
