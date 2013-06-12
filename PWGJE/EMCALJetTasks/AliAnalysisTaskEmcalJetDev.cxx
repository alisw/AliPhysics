// $Id$
//
// Emcal jet analysis base task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcalJetDev.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskEmcalJetDev)

//________________________________________________________________________
AliAnalysisTaskEmcalJetDev::AliAnalysisTaskEmcalJetDev() : 
  AliAnalysisTaskEmcalDev("AliAnalysisTaskEmcalJetDev"),
  fJetsName(),
  fRhoName(),
  fJets(0),
  fRho(0),
  fRhoVal(0),
  fJetCollArray()
{
  // Default constructor.

  fJetCollArray.SetOwner(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetDev::AliAnalysisTaskEmcalJetDev(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcalDev(name, histo),
  fJetsName(),
  fRhoName(),
  fJets(0),
  fRho(0),
  fRhoVal(0),
  fJetCollArray()
{
  // Standard constructor.

  fJetCollArray.SetOwner(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetDev::~AliAnalysisTaskEmcalJetDev()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetDev::AcceptBiasJet(AliEmcalJet *jet, Int_t c)
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
Float_t* AliAnalysisTaskEmcalJetDev::GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const
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
Double_t AliAnalysisTaskEmcalJetDev::GetLeadingHadronPt(AliEmcalJet *jet, Int_t c)
{

  AliJetContainer *cont = GetJetContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->GetLeadingHadronPt(jet);

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetDev::AcceptJet(AliEmcalJet *jet, Int_t c) 
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
AliRhoParameter *AliAnalysisTaskEmcalJetDev::GetRhoFromEvent(const char *name)
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
void AliAnalysisTaskEmcalJetDev::ExecOnce()
{
  // Init the analysis.

  AliAnalysisTaskEmcalDev::ExecOnce();

  //Load all requested track branches - each container knows name already
  if(fJetCollArray.GetEntriesFast()==0) {
    AliWarning("There are no jet collections");
    return;
  }

  for(Int_t i =0; i<fJetCollArray.GetEntriesFast(); i++) {
    AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
    cont->SetJetArray(InputEvent());
    cont->SetEMCALGeometry();
    cont->LoadRho(InputEvent());
  }

  //Get Jets, cuts and rho for first jet container
  AliJetContainer *cont = GetJetContainer(0);
  // fJetsName = cont->GetArrayName();
  if (fAnaType == kTPC) {
    cont->SetJetAcceptanceType(AliJetContainer::kTPC);
    cont->SetJetEtaPhiTPC();
  }
  else if (fAnaType == kEMCAL) {
    cont->SetJetAcceptanceType(AliJetContainer::kEMCAL);
    cont->SetJetEtaPhiEMCAL();
  }
  
  fJets = GetJetArray(0);
  if(!fJets && fJetCollArray.GetEntriesFast()>0) {
    AliError(Form("%s: Could not retrieve first jet branch!", GetName()));
    fInitialized = kFALSE;
    return;
  }

  fRhoName = cont->GetRhoName();
  if(!fRhoName.IsNull()) {
    fRho = cont->GetRhoParameter();
    if(!fRho) {
      AliError(Form("%s: Could not retrieve rho of first jet branch!", GetName()));
      fInitialized = kFALSE;
      return;
    }
  }

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetDev::GetSortedArray(Int_t indexes[], TClonesArray *array, Double_t rho, Int_t c)
{
  // Get the leading jets.

  static Float_t pt[9999] = {0};

  if (!array)
    return 0;

  const Int_t n = array->GetEntriesFast();

  if (n < 1)
    return kFALSE;
  
  if (array->GetClass()->GetBaseClass("AliEmcalJet")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliEmcalJet* jet = static_cast<AliEmcalJet*>(array->At(i));
      
      if (!jet) {
	AliError(Form("Could not receive jet %d", i));
	continue;
      }
      
      if (!AcceptJet(jet,c))
	continue;
      
      pt[i] = jet->Pt() - rho * jet->Area();
    }
  }

  else if (array->GetClass()->GetBaseClass("AliVTrack")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliVTrack* track = static_cast<AliVTrack*>(array->At(i));
      
      if (!track) {
	AliError(Form("Could not receive track %d", i));
	continue;
      }  
      
      if (!AcceptTrack(track, c))
	continue;
      
      pt[i] = track->Pt();
    }
  }

  else if (array->GetClass()->GetBaseClass("AliVCluster")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliVCluster* cluster = static_cast<AliVCluster*>(array->At(i));
      
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", i));
	continue;
      }  
      
      if (!AcceptCluster(cluster, c))
	continue;

      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
      
      pt[i] = nPart.Pt();
    }
  }

  TMath::Sort(n, pt, indexes);

  if (pt[indexes[0]] == -FLT_MAX) 
    return 0;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetDev::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
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
Bool_t AliAnalysisTaskEmcalJetDev::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
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
Bool_t AliAnalysisTaskEmcalJetDev::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!AliAnalysisTaskEmcalDev::RetrieveEventObjects())
    return kFALSE;

  if (fRho)
    fRhoVal = fRho->GetVal();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::AddJetContainer(const char *n, TString defaultCutType) {

  // Add particle container
  // will be called in AddTask macro

  AliJetContainer *cont = 0x0;
  cont = new AliJetContainer();
  cont->SetArrayName(n);

  if(!defaultCutType.IsNull()) {
    if(defaultCutType.EqualTo("TPC"))
     cont->SetJetAcceptanceType(AliJetContainer::kTPC);
    else if(defaultCutType.EqualTo("EMCAL"))
      cont->SetJetAcceptanceType(AliJetContainer::kEMCAL);
    else
      AliWarning(Form("%s: default cut type %s not recognized. Not setting cuts.",GetName(),defaultCutType.Data()));
  } else
    cont->SetJetAcceptanceType(AliJetContainer::kUser);
 
  fJetCollArray.Add(cont);

}

//________________________________________________________________________
AliJetContainer* AliAnalysisTaskEmcalJetDev::GetJetContainer(Int_t i) const{
  // Get i^th jet container

  if(i<0 || i>fJetCollArray.GetEntriesFast()) return 0;
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
  return cont;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetRhoName(const char *n, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetRhoName(n);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetEtaLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetEtaLimits(min,max);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetPhiLimits(Float_t min, Float_t max, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetPhiLimits(min,max);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetAreaCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetAreaCut(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetPercAreaCut(Float_t p, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetPercAreaCut(p);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetAreaEmcCut(Double_t a, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetAreaEmcCut(a);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetPtCut(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetPtCut(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetRadius(Float_t r, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetRadius(r);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetMaxClusterPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetMaxClusterPt(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetMaxTrackPt(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetMaxTrackPt(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetPtBiasJetClus(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetPtBiasJetClus(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetPtBiasJetTrack(Float_t cut, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetPtBiasJetTrack(cut);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetLeadingHadronType(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetLeadingHadronType(t);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetNLeadingJets(Int_t t, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetNLeadingJets(t);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetDev::SetJetBitMap(UInt_t m, Int_t c)
{
  AliJetContainer *cont = GetJetContainer(c);
  cont->SetJetBitMap(m);
}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcalJetDev::GetJetArray(Int_t i) const {
  // Get i^th TClonesArray with AliEmcalJet

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();

}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalJetDev::GetAcceptJetFromArray(Int_t j, Int_t c) const {
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
Int_t AliAnalysisTaskEmcalJetDev::GetNJets(Int_t i) const {
  // Get number of entries in jet array i

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNJets();

}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetDev::GetRhoVal(Int_t i) const {
  // Get rho corresponding to jet array i

  AliJetContainer *cont = GetJetContainer(i);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetRhoVal();

}



