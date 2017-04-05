/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliAnalysisTaskEmcalJetLight.h"

#include <TClonesArray.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliLocalRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetLight);
/// \endcond

/**
 * Default constructor (should only be used by ROOT I/O).
 */
AliAnalysisTaskEmcalJetLight::AliAnalysisTaskEmcalJetLight() :
  AliAnalysisTaskEmcalLight(),
  fJetCollArray()
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * Note: This constructor also handles the general histograms. In
 * case the second parameter is true, then general histograms (see
 * UserCreateOutputObjects and FillHistograms) are created and filled
 * by the task, and a container is provided handling the user histograms.
 * @param[in] name Name of the task
 * @param[in] histo If true then general histograms are filled by the task
 */
AliAnalysisTaskEmcalJetLight::AliAnalysisTaskEmcalJetLight(const char *name, Bool_t histo) :
  AliAnalysisTaskEmcalLight(name, histo),
  fJetCollArray()
{
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetLight::~AliAnalysisTaskEmcalJetLight()
{
  for (auto cont_it : fJetCollArray) delete cont_it.second;
}

/**
 * Perform steps needed to initialize the analysis.
 * This function relies on the presence of an input
 * event (ESD or AOD event). Consequently it is called
 * internally by UserExec for the first event.
 *
 * This function connects all containers attached to
 * this task to the corresponding arrays in the
 * input event. Furthermore it initializes the geometry.
 */
void AliAnalysisTaskEmcalJetLight::ExecOnce()
{
  AliAnalysisTaskEmcalLight::ExecOnce();

  //Load all requested jet branches - each container knows name already
  if (fJetCollArray.size() == 0) {
    AliWarning("There are no jet collections");
    return;
  }

  for (auto cont_it : fJetCollArray) {
    AliJetContainer *cont = cont_it.second;
    cont->SetRunNumber(InputEvent()->GetRunNumber());
    cont->SetArray(InputEvent());
    cont->LoadRho(InputEvent());
  }
}

/**
 * Retrieve objects from event. This operation needs to be performed
 * for every event.
 * @return kTRUE if successful, kFALSE otherwise
 */
Bool_t AliAnalysisTaskEmcalJetLight::RetrieveEventObjects()
{
  if (!AliAnalysisTaskEmcalLight::RetrieveEventObjects()) return kFALSE;

  for (auto cont_it : fJetCollArray) cont_it.second->NextEvent();

  return kTRUE;
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] jetType One of the AliJetContainer::EJetType_t enumeration values (charged, full, neutral)
 * @param[in] jetAlgo One of the AliJetContainer::EJetAlgo_t enumeration values (anti-kt, kt, ...)
 * @param[in] recoScheme One of the AliJetContainer::ERecoScheme_t enumeration values (pt-scheme, ...)
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @param[in] accType One of the AliEmcalJet::JetAcceptanceType enumeration values (TPC, EMCAL, user, ...)
 * @param[in] partContName Name of the particle container of the objects used to generate the jets
 * @param[in] clusContName Name of the cluster container of the objects used to generate the jets
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
    UInt_t accType, std::string partContName, std::string clusContName, TString tag)
{
  AliParticleContainer* partCont = nullptr;
  AliClusterContainer* clusCont =  nullptr;
  auto partContSearch = fParticleCollArray.find(partContName);
  if (partContSearch != fParticleCollArray.end()) partCont = partContSearch->second;
  auto clusContSearch = fClusterCollArray.find(partContName);
  if (clusContSearch != fClusterCollArray.end()) clusCont = clusContSearch->second;
  if (!partCont && !clusCont) {
    AliError(Form("Could not find neither particle nor cluster container with names '%s' and '%s'", partContName.c_str(), clusContName.c_str()));
    return nullptr;
  }

  return AddJetContainer(jetType, jetAlgo, recoScheme, radius, accType, partCont, clusCont, tag);
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] jetType One of the AliJetContainer::EJetType_t enumeration values (charged, full, neutral)
 * @param[in] jetAlgo One of the AliJetContainer::EJetAlgo_t enumeration values (anti-kt, kt, ...)
 * @param[in] recoScheme One of the AliJetContainer::ERecoScheme_t enumeration values (pt-scheme, ...)
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @param[in] accType One of the AliEmcalJet::JetAcceptanceType enumeration values (TPC, EMCAL, user, ...)
 * @param[in] partCont Particle container of the objects used to generate the jets
 * @param[in] clusCont Cluster container of the objects used to generate the jets
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, UInt_t accType,
    AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
{
  AliJetContainer *cont = new AliJetContainer(jetType, jetAlgo, recoScheme, radius, partCont, clusCont, tag);
  cont->SetJetAcceptanceType(accType);
  AdoptJetContainer(cont);

  return cont;
}

/**
 * Find jet container attached to this task according to its name
 * @param[in] name Name of the jet container
 * @return Jet container found under the given name
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::GetJetContainer(std::string name) const
{
  std::map<std::string, AliJetContainer*>::const_iterator cont_it = fJetCollArray.find(name);
  if (cont_it != fJetCollArray.end()) return cont_it->second;
  else return nullptr;
}
