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

#include "AliAnalysisTaskEmcalJet.h"

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
ClassImp(AliAnalysisTaskEmcalJet);
/// \endcond

/**
 * Default constructor (should only be used by ROOT I/O).
 */
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalJet"),
  fRhoName(),
  fLocalRhoName(),
  fJetCollArray(),
  fJets(0),
  fRho(0),
  fLocalRho(0),
  fRhoVal(0)
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
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fRhoName(),
  fLocalRhoName(),
  fJetCollArray(),
  fJets(0),
  fRho(0),
  fLocalRho(0),
  fRhoVal(0)
{
  fJetCollArray.SetOwner(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJet::~AliAnalysisTaskEmcalJet()
{
}

/**
 * Get the rho object from the event object list. The rho object contains
 * the value of the average event background.
 * @param name Name of the object used to search in the event object list
 * @return Pointer to the rho object
 */
AliRhoParameter *AliAnalysisTaskEmcalJet::GetRhoFromEvent(const char *name)
{
  AliRhoParameter *rho = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(sname));
    if (!rho) {
      AliWarning(Form("%s: Could not retrieve rho with name %s!", GetName(), name)); 
      return 0;
    }
  }
  return rho;
}

/**
 * Get the local rho object from the event object list. The local rho object contains
 * information about the average event background as a function of the event geometry.
 * @param name Name of the object used to search in the event object list
 * @return Pointer to the local rho object
 */
AliLocalRhoParameter *AliAnalysisTaskEmcalJet::GetLocalRhoFromEvent(const char *name)
{
  // Get local rho from event.
  AliLocalRhoParameter *rho = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    rho = dynamic_cast<AliLocalRhoParameter*>(InputEvent()->FindListObject(sname));
    if (!rho) {
      AliWarning(Form("%s: Could not retrieve local rho with name %s!", GetName(), name));
      return 0;
    }
  }
  return rho;
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
void AliAnalysisTaskEmcalJet::ExecOnce()
{
  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fRhoName.IsNull() && !fRho) { // get rho from the event
    fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (!fLocalRhoName.IsNull() && !fLocalRho) {
    fLocalRho = dynamic_cast<AliLocalRhoParameter*>(InputEvent()->FindListObject(fLocalRhoName));
    if (!fLocalRho) {
      AliError(Form("%s: Could not retrieve local rho %s!", GetName(), fLocalRhoName.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

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

  //Get Jets, cuts and rho for first jet container
  AliJetContainer *cont = GetJetContainer(0);
  
  if (!cont->GetArrayName().IsNull()) {
    fJets = cont->GetArray();
    if(!fJets && fJetCollArray.GetEntriesFast()>0) {
      AliError(Form("%s: Could not retrieve first jet branch!", GetName()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (!fRho) { // if rho name is not provided, tries to use the rho object of the first jet branch
    fRhoName = cont->GetRhoName();
    fRho = cont->GetRhoParameter();
  }
}

/**
 * Checks whether a cluster is among the constituents of a jet.
 * @param jet Pointer to an AliEmcalJet object
 * @param iclus Index of the cluster to look for
 * @param sorted If the constituents are sorted by index it will speed up computation
 * @return kTRUE if the cluster is among the constituents, kFALSE otherwise
 */
Bool_t AliAnalysisTaskEmcalJet::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfClusters(); ++i) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

/**
 * Checks whether a track is among the constituents of a jet.
 * @param jet Pointer to an AliEmcalJet object
 * @param itrack Index of the track to look for
 * @param sorted If the constituents are sorted by index it will speed up computation
 * @return kTRUE if the track is among the constituents, kFALSE otherwise
 */
Bool_t AliAnalysisTaskEmcalJet::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); ++i) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

/**
 * Retrieve objects from event. This operation needs to be performed
 * for every event.
 * @return kTRUE if successful, kFALSE otherwise
 */
Bool_t AliAnalysisTaskEmcalJet::RetrieveEventObjects()
{
  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  if (fRho) fRhoVal = fRho->GetVal();

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
 * @param[in] accType One of the AliJetContainer::JetAcceptanceType enumeration values (kTPC, kEMCAL, kDCAL, ...),
 * or a combination using bitwise OR: For example, (kEMCAL | kDCAL) will select all jets in either EMCal or DCal.
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJet::AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
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
 * @param[in] accType One of the AliJetContainer::JetAcceptanceType enumeration values (kTPC, kEMCAL, kDCAL, ...),
 * or a combination using bitwise OR: For example, (kEMCAL | kDCAL) will select all jets in either EMCal or DCal.
 * @param[in] partCont Particle container of the objects used to generate the jets
 * @param[in] clusCont Cluster container of the objects used to generate the jets
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJet::AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, UInt_t accType,
    AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
{
  AliJetContainer *cont = new AliJetContainer(jetType, jetAlgo, recoScheme, radius, partCont, clusCont, tag);
  cont->SetJetAcceptanceType(accType);
  fJetCollArray.Add(cont);

  return cont;
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] n Name of the jet branch
 * @param[in] accType One of the AliJetContainer::JetAcceptanceType enumeration values (kTPC, kEMCAL, kDCAL, ...),
 * or a combination using bitwise OR: For example, (kEMCAL | kDCAL) will select all jets in either EMCal or DCal.
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJet::AddJetContainer(const char *n, UInt_t accType, Float_t jetRadius)
{
  if (TString(n).IsNull()) return 0;

  AliJetContainer *cont = new AliJetContainer(n);
  cont->SetJetRadius(jetRadius);
  cont->SetJetAcceptanceType(accType);
  fJetCollArray.Add(cont);

  return cont;
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * @param[in] n Name of the jet branch
 * @param[in] defaultCutType String that correspond to a possible acceptance cut type
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskEmcalJet::AddJetContainer(const char *n, TString defaultCutType, Float_t jetRadius) {

  // Add particle container
  // will be called in AddTask macro

  if(TString(n).IsNull()) return 0;

  UInt_t acc = AliJetContainer::kUser;

  defaultCutType.ToUpper();

  if (defaultCutType.IsNull() || defaultCutType.EqualTo("USER")) {
    acc = AliJetContainer::kUser;
  }
  else if(defaultCutType.EqualTo("TPC")) {
    acc = AliJetContainer::kTPC;
  }
  else if(defaultCutType.EqualTo("TPCFID")) {
    acc = AliJetContainer::kTPCfid;
  }
  else if(defaultCutType.EqualTo("EMCAL")) {
    acc = AliJetContainer::kEMCAL;
  }
  else if(defaultCutType.EqualTo("EMCALFID")) {
    acc = AliJetContainer::kEMCALfid;
  }
  else if(defaultCutType.EqualTo("DCAL")) {
    acc = AliJetContainer::kDCAL;
  }
  else if(defaultCutType.EqualTo("DCALFID")) {
    acc = AliJetContainer::kDCALfid;
  }
  else {
    AliWarning(Form("%s: default cut type %s not recognized. Not setting cuts.",GetName(),defaultCutType.Data()));
  }

  return AddJetContainer(n, acc, jetRadius);
}

/**
 * Get \f$ i^{th} \f$ jet container attached to this task
 * @param[in] i Index of the jet container
 * @return Jet container found for the given index (NULL if no jet container exists for that index)
 */
AliJetContainer* AliAnalysisTaskEmcalJet::GetJetContainer(Int_t i) const
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
AliJetContainer* AliAnalysisTaskEmcalJet::GetJetContainer(const char* name) const
{
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.FindObject(name));
  return cont;
}

void AliAnalysisTaskEmcalJet::SetJetAcceptanceType(UInt_t t, Int_t c) 
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) {
    cont->SetJetAcceptanceType(t);
  }
  else {
    AliError(Form("%s in SetJetAcceptanceType(...): container %d not found!",GetName(),c));
  }
}

void AliAnalysisTaskEmcalJet::SetJetAcceptanceType(TString cutType, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (!cont) {
    AliError(Form("%s in SetJetAcceptanceType(...): container %d not found",GetName(),c));
    return;
  }

  cutType.ToUpper();

  if(!cutType.IsNull() && !cutType.EqualTo("USER")) {
    if(cutType.EqualTo("TPC"))
     cont->SetJetAcceptanceType(AliJetContainer::kTPC);
    else if(cutType.EqualTo("EMCAL"))
      cont->SetJetAcceptanceType(AliJetContainer::kEMCAL);
    else
      AliWarning(Form("%s: default cut type %s not recognized. Not setting cuts.",GetName(),cutType.Data()));
  } else
    cont->SetJetAcceptanceType(AliJetContainer::kUser);
}

void AliAnalysisTaskEmcalJet::SetRhoName(const char *n, Int_t c)
{
  if (c >= 0) {
    AliJetContainer *cont = GetJetContainer(c);
    if (cont) cont->SetRhoName(n);
    else AliError(Form("%s in SetRhoName(...): container %d not found",GetName(),c));
  }
  else {
    fRhoName = n;
  }
}

void AliAnalysisTaskEmcalJet::SetJetEtaLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetEtaLimits(min,max);
  else AliError(Form("%s in SetJetEtaLimits(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetJetPhiLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetPhiLimits(min,max);
  else AliError(Form("%s in SetJetPhiLimits(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetJetAreaCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetAreaCut(cut);
  else AliError(Form("%s in SetJetAreaCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetPercAreaCut(Float_t p, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPercAreaCut(p);
  else AliError(Form("%s in SetPercAreaCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetZLeadingCut(Float_t zemc, Float_t zch, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetZLeadingCut(zemc,zch);
  else AliError(Form("%s in SetZLeadingCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetNEFCut(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetNEFCut(min,max);
  else AliError(Form("%s in SetNEFCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetAreaEmcCut(Double_t a, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetAreaEmcCut(a);
  else AliError(Form("%s in SetAreaEmcCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetJetPtCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetPtCut(cut);
  else AliError(Form("%s in SetJetPtCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetJetRadius(Float_t r, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetRadius(r);
  else AliError(Form("%s in SetJetRadius(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetMaxClusterPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetMaxClusterPt(cut);
  else AliError(Form("%s in SetMaxClusterPt(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetMaxTrackPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetMaxTrackPt(cut);
  else AliError(Form("%s in SetMaxTrackPt(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetPtBiasJetClus(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPtBiasJetClus(cut);
  else AliError(Form("%s in SetPtBiasJetClus(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetPtBiasJetTrack(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPtBiasJetTrack(cut);
  else AliError(Form("%s in SetPtBiasJetTrack(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetLeadingHadronType(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetLeadingHadronType(t);
  else AliError(Form("%s in SetLeadingHadronType(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetNLeadingJets(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetNLeadingJets(t);
  else AliError(Form("%s in SetNLeadingJets(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetJetBitMap(UInt_t m, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetBitMap(m);
  else AliError(Form("%s in SetJetBitMap(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcalJet::SetIsParticleLevel(Bool_t b, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetIsParticleLevel(b);
  else AliError(Form("%s in SetIsParticleLevel(...): container %d not found",GetName(),c));
}

const TString& AliAnalysisTaskEmcalJet::GetRhoName(Int_t c) const
{
  if (c >= 0) {
    AliJetContainer *cont = GetJetContainer(c);
    if (cont) return cont->GetRhoName();
    else { AliError(Form("%s in GetRhoName(...): container %d not found. Returning fRhoName...",GetName(),c)); return fRhoName; }
  }
  else {
    return fRhoName;
  }
}

TClonesArray* AliAnalysisTaskEmcalJet::GetJetArray(Int_t i) const
{
  AliJetContainer *cont = GetJetContainer(i);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();
}

Double_t AliAnalysisTaskEmcalJet::GetJetRadius(Int_t i) const
{
  AliJetContainer *cont = GetJetContainer(i);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }

  return cont->GetJetRadius();
}

AliEmcalJet* AliAnalysisTaskEmcalJet::GetJetFromArray(Int_t j, Int_t c) const
{
  AliJetContainer *cont = GetJetContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }
  AliEmcalJet *jet = cont->GetJet(j);

  return jet;
}

AliEmcalJet* AliAnalysisTaskEmcalJet::GetAcceptJetFromArray(Int_t j, Int_t c) const
{
  AliJetContainer *cont = GetJetContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }
  AliEmcalJet *jet = cont->GetAcceptJet(j);

  return jet;
}

Int_t AliAnalysisTaskEmcalJet::GetNJets(Int_t i) const
{
  AliJetContainer *cont = GetJetContainer(i);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNJets();

}

Double_t AliAnalysisTaskEmcalJet::GetRhoVal(Int_t i) const
{
  AliJetContainer *cont = GetJetContainer(i);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetRhoVal();
}

Double_t AliAnalysisTaskEmcalJet::GetLeadingHadronPt(AliEmcalJet *jet, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->GetLeadingHadronPt(jet);
}

Bool_t AliAnalysisTaskEmcalJet::AcceptJet(AliEmcalJet *jet, Int_t c)
{
  if (!jet)
    return kFALSE;

  AliJetContainer *cont = GetJetContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  UInt_t rejectionReason = 0;
  return cont->AcceptJet(jet, rejectionReason);
}
