// $Id$
//
// Emcal jet analysis base task.
//
// Author: S.Aiola, M. Verweij

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
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskEmcalJet)

//________________________________________________________________________
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
  // Default constructor.

  fJetCollArray.SetOwner(kTRUE);
}

//________________________________________________________________________
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
  // Standard constructor.

  fJetCollArray.SetOwner(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::~AliAnalysisTaskEmcalJet()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::AcceptBiasJet(AliEmcalJet *jet, Int_t c)
{ 
  // Accept jet with a bias.

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptBiasJet(jet);
}

//________________________________________________________________________
Float_t* AliAnalysisTaskEmcalJet::GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const
{
  Float_t *bins = new Float_t[n+1];

  Float_t binWidth = (max-min)/n;
  bins[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    bins[i] = bins[i-1]+binWidth;
  }

  return bins;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJet::GetLeadingHadronPt(AliEmcalJet *jet, Int_t c)
{

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->GetLeadingHadronPt(jet);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::AcceptJet(AliEmcalJet *jet, Int_t c) 
{   
  // Return true if jet is accepted.
  if (!jet)
    return kFALSE;

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptJet(jet);
}

//________________________________________________________________________
AliRhoParameter *AliAnalysisTaskEmcalJet::GetRhoFromEvent(const char *name)
{
  // Get rho from event.

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

//________________________________________________________________________
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


//________________________________________________________________________
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
    cont->SetEMCALGeometry();
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

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  // Return true if cluster is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfClusters(); ++i) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  // Return true if track is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfTracks(); ++i) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  if (fRho)
    fRhoVal = fRho->GetVal();

  return kTRUE;
}

//________________________________________________________________________
AliJetContainer* AliAnalysisTaskEmcalJet::AddJetContainer(const char *n, TString defaultCutType, Float_t jetRadius) {

  // Add particle container
  // will be called in AddTask macro

  TString tmp = TString(n);
  if(tmp.IsNull()) return 0;

  AliJetContainer *cont = 0x0;
  cont = new AliJetContainer();
  cont->SetArrayName(n);
  cont->SetJetRadius(jetRadius);

  defaultCutType.ToUpper();

  if(!defaultCutType.IsNull() && !defaultCutType.EqualTo("USER")) {
    if(defaultCutType.EqualTo("TPC"))
      cont->SetJetAcceptanceType(AliJetContainer::kTPC);
    else if(defaultCutType.EqualTo("EMCAL"))
      cont->SetJetAcceptanceType(AliJetContainer::kEMCAL);
    else
      AliWarning(Form("%s: default cut type %s not recognized. Not setting cuts.",GetName(),defaultCutType.Data()));
  } else
    cont->SetJetAcceptanceType(AliJetContainer::kUser);
 
  fJetCollArray.Add(cont);

  return cont;
}

//________________________________________________________________________
AliJetContainer* AliAnalysisTaskEmcalJet::GetJetContainer(Int_t i) const{
  // Get i^th jet container

  if(i<0 || i>=fJetCollArray.GetEntriesFast()) return 0;
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
  return cont;
}

//________________________________________________________________________
AliJetContainer* AliAnalysisTaskEmcalJet::GetJetContainer(const char* name) const{
  // Get the jet container with name
  
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.FindObject(name));
  return cont;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetAcceptanceType(UInt_t t, Int_t c) 
{
  // Set acceptance cuts
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) {
    cont->SetJetAcceptanceType((AliJetContainer::JetAcceptanceType)t);
  }
  else {
    AliError(Form("%s in SetJetAcceptanceType(...): container %d not found!",GetName(),c));
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetAcceptanceType(TString cutType, Int_t c) {
  //set acceptance cuts

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

//________________________________________________________________________
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

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetEtaLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetEtaLimits(min,max);
  else AliError(Form("%s in SetJetEtaLimits(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetPhiLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetPhiLimits(min,max);
  else AliError(Form("%s in SetJetPhiLimits(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetAreaCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetAreaCut(cut);
  else AliError(Form("%s in SetJetAreaCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetPercAreaCut(Float_t p, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPercAreaCut(p);
  else AliError(Form("%s in SetPercAreaCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetZLeadingCut(Float_t zemc, Float_t zch, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetZLeadingCut(zemc,zch);
  else AliError(Form("%s in SetZLeadingCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetNEFCut(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetNEFCut(min,max);
  else AliError(Form("%s in SetNEFCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetAreaEmcCut(Double_t a, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetAreaEmcCut(a);
  else AliError(Form("%s in SetAreaEmcCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetPtCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetPtCut(cut);
  else AliError(Form("%s in SetJetPtCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetRadius(Float_t r, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetRadius(r);
  else AliError(Form("%s in SetJetRadius(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetMaxClusterPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetMaxClusterPt(cut);
  else AliError(Form("%s in SetMaxClusterPt(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetMaxTrackPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetMaxTrackPt(cut);
  else AliError(Form("%s in SetMaxTrackPt(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetPtBiasJetClus(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPtBiasJetClus(cut);
  else AliError(Form("%s in SetPtBiasJetClus(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetPtBiasJetTrack(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetPtBiasJetTrack(cut);
  else AliError(Form("%s in SetPtBiasJetTrack(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetLeadingHadronType(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetLeadingHadronType(t);
  else AliError(Form("%s in SetLeadingHadronType(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetNLeadingJets(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetNLeadingJets(t);
  else AliError(Form("%s in SetNLeadingJets(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetJetBitMap(UInt_t m, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetJetBitMap(m);
  else AliError(Form("%s in SetJetBitMap(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::SetIsParticleLevel(Bool_t b, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  if (cont) cont->SetIsParticleLevel(b);
  else AliError(Form("%s in SetIsParticleLevel(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
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

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcalJet::GetJetArray(Int_t i) const {
  // Get i^th TClonesArray with AliEmcalJet

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJet::GetJetRadius(Int_t i) const {
  // Get jet radius from jet container i

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }

  return cont->GetJetRadius();
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalJet::GetJetFromArray(Int_t j, Int_t c) const {
  // Get jet j if accepted from  container c
  // If jet not accepted return 0

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }
  AliEmcalJet *jet = cont->GetJet(j);

  return jet;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalJet::GetAcceptJetFromArray(Int_t j, Int_t c) const {
  // Get jet j if accepted from  container c
  // If jet not accepted return 0

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }
  AliEmcalJet *jet = cont->GetAcceptJet(j);

  return jet;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJet::GetNJets(Int_t i) const {
  // Get number of entries in jet array i

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNJets();

}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJet::GetRhoVal(Int_t i) const {
  // Get rho corresponding to jet array i

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetRhoVal();
}
