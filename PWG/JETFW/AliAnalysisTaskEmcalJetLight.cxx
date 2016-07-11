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
  fJetCollArray.SetOwner(kTRUE);
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
  fJetCollArray.SetOwner(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetLight::~AliAnalysisTaskEmcalJetLight()
{
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
  if(fJetCollArray.GetEntriesFast()==0) {
    AliWarning("There are no jet collections");
    return;
  }

  for(Int_t i =0; i<fJetCollArray.GetEntriesFast(); i++) {
    AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
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
  if (!AliAnalysisTaskEmcalLight::RetrieveEventObjects())
    return kFALSE;

  AliEmcalContainer* cont = 0;

  TIter nextJetColl(&fJetCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextJetColl()))) cont->NextEvent();

  return kTRUE;
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] jetType One of the AliJetContainer::EJetType_t enumeration values (charged, full, neutral)
 * @param[in] jetAlgo One of the AliJetContainer::EJetAlgo_t enumeration values (anti-kt, kt, ...)
 * @param[in] recoScheme One of the AliJetContainer::ERecoScheme_t enumeration values (pt-scheme, ...)
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @param[in] accType One of the AliJetContainer::JetAcceptanceType enumeration values (TPC, EMCAL, user, ...)
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
    UInt_t accType, TString tag)
{
  AliParticleContainer* partCont = GetParticleContainer(0);
  AliClusterContainer* clusCont = GetClusterContainer(0);

  return AddJetContainer(jetType, jetAlgo, recoScheme, radius, accType, partCont, clusCont, tag);
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] jetType One of the AliJetContainer::EJetType_t enumeration values (charged, full, neutral)
 * @param[in] jetAlgo One of the AliJetContainer::EJetAlgo_t enumeration values (anti-kt, kt, ...)
 * @param[in] recoScheme One of the AliJetContainer::ERecoScheme_t enumeration values (pt-scheme, ...)
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @param[in] accType One of the AliJetContainer::JetAcceptanceType enumeration values (TPC, EMCAL, user, ...)
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
  fJetCollArray.Add(cont);

  return cont;
}

/**
 * Get \f$ i^{th} \f$ jet container attached to this task
 * @param[in] i Index of the jet container
 * @return Jet container found for the given index (NULL if no jet container exists for that index)
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::GetJetContainer(Int_t i) const
{
  if (i < 0 || i >= fJetCollArray.GetEntriesFast()) return 0;
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
  return cont;
}

/**
 * Find jet container attached to this task according to its name
 * @param[in] name Name of the jet container
 * @return Jet container found under the given name
 */
AliJetContainer* AliAnalysisTaskEmcalJetLight::GetJetContainer(const char* name) const{
  // Get the jet container with name
  
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.FindObject(name));
  return cont;
}
