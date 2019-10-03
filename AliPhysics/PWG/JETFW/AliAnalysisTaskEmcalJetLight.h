#ifndef ALIANALYSISTASKEMCALJETLIGHT_H
#define ALIANALYSISTASKEMCALJETLIGHT_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TList;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliLocalRhoParameter;
class AliVCluster;
class AliVParticle;

#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalLight.h"

/**
 * @class AliAnalysisTaskEmcalJetLight
 * @brief Base task in the EMCAL jet framework (lighter version of AliAnalysisTaskEmcalJet)
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class is the base class for Analysis Tasks using the
 * core EMCAL jet framework. User tasks can choose to inherit from it
 * as an alternative to AliAnalysisTaskEmcalJet.
 * This class derives from AliAnalysisTaskEmcalLight and adds
 * some additional features, useful for jet analysis. The key feature
 * is the handling of the jet containers.
 */
class AliAnalysisTaskEmcalJetLight : public AliAnalysisTaskEmcalLight {
 public:
  typedef AliJetContainer::EJetType_t EJetType_t;
  typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
  typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

  AliAnalysisTaskEmcalJetLight();
  AliAnalysisTaskEmcalJetLight(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskEmcalJetLight();

  AliJetContainer            *AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
      UInt_t accType, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag = "Jet");
  AliJetContainer            *AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
      UInt_t accType, std::string partContName, std::string clusContName, TString tag = "Jet");

  void                        AdoptJetContainer(AliJetContainer* cont)           { fJetCollArray[cont->GetName()] = cont; }
  void                        RemoveJetContainer(std::string name)               { fJetCollArray.erase(name);}
  AliJetContainer            *GetJetContainer(std::string name)                                        const;

 protected:
  void                        ExecOnce()                                                                    ;
  Bool_t                      RetrieveEventObjects()                                                        ;

  std::map<std::string,
    AliJetContainer*>         fJetCollArray;               ///< jet collection array

 private:
  AliAnalysisTaskEmcalJetLight(const AliAnalysisTaskEmcalJetLight&);            // not implemented
  AliAnalysisTaskEmcalJetLight &operator=(const AliAnalysisTaskEmcalJetLight&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetLight, 2);
  /// \endcond
};
#endif
