/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliTLorentzVector.h"

#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliClusterContainer)
/// \endcond

/**
 * Default constructor.
 */
AliClusterContainer::AliClusterContainer():
  AliEmcalContainer(),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1)
{
  fClassName = "AliVCluster";

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}

/**
 * Standard constructor.
 * @param name Name of the array connected to this container
 */
AliClusterContainer::AliClusterContainer(const char *name):
  AliEmcalContainer(name),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1)
{
  fClassName = "AliVCluster";

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}

/**
 * Get the leading cluster; use e if "e" is contained in opt (otherwise et)
 * @param opt
 * @return
 */
AliVCluster* AliClusterContainer::GetLeadingCluster(const char* opt)
{
  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliVCluster *clusterMax = GetNextAcceptCluster();
  AliVCluster *cluster = 0;

  if (option.Contains("e")) {
    while ((cluster = GetNextAcceptCluster())) {
      if (cluster->E() > clusterMax->E()) clusterMax = cluster;
    }
  }
  else {
    Double_t et = 0;
    Double_t etmax = 0;
    while ((cluster = GetNextAcceptCluster())) {
      TLorentzVector mom;
      cluster->GetMomentum(mom,const_cast<Double_t*>(fVertex));
      et = mom.Et();
      if (et > etmax) { 
	clusterMax = cluster;
	etmax = et;
      }
    }
  }

  fCurrentID = tempID;

  return clusterMax;
}

/**
 * Get i^th cluster in array
 * @param i
 * @return
 */
AliVCluster* AliClusterContainer::GetCluster(Int_t i) const 
{
  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVCluster *vp = static_cast<AliVCluster*>(fClArray->At(i));
  return vp;

}

/**
 * Return pointer to cluster if cluster is accepted
 * @param i
 * @return
 */
AliVCluster* AliClusterContainer::GetAcceptCluster(Int_t i) 
{
  AliVCluster *vc = GetCluster(i);
  if (!vc) return 0;

  if (AcceptCluster(vc))
    return vc;
  else {
    AliDebug(2,"Cluster not accepted.");
    return 0;
  }
}

/**
 * Get particle with label lab in array
 * @param lab
 * @return
 */
AliVCluster* AliClusterContainer::GetClusterWithLabel(Int_t lab) const 
{
  Int_t i = GetIndexFromLabel(lab);
  return GetCluster(i);
}

/**
 * Get particle with label lab in array
 * @param lab
 * @return
 */
AliVCluster* AliClusterContainer::GetAcceptClusterWithLabel(Int_t lab)  
{
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptCluster(i);
}

/**
 * Get next accepted cluster
 * @return
 */
AliVCluster* AliClusterContainer::GetNextAcceptCluster() 
{
  const Int_t n = GetNEntries();
  AliVCluster *c = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    c = GetAcceptCluster(fCurrentID);
  } while (!c);

  return c;
}

/**
 * Get next cluster
 * @return
 */
AliVCluster* AliClusterContainer::GetNextCluster() 
{
  const Int_t n = GetNEntries();
  AliVCluster *c = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    c = GetCluster(fCurrentID);
  } while (!c);

  return c;
}

Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, const AliVCluster* vc, Double_t mass)
{
  if (mass < 0) mass = 0;

  Double_t energy = 0;

  if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
    energy = vc->GetUserDefEnergy((AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
  }
  else {
    energy = vc->E();
  }

  Double_t p = TMath::Sqrt(energy*energy - mass*mass);

  Float_t pos[3];
  vc->GetPosition(pos);

  pos[0]-=fVertex[0];
  pos[1]-=fVertex[1];
  pos[2]-=fVertex[2];

  Double_t r = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]) ;

  if (r > 1e-12) {
    mom.SetPxPyPzE( p*pos[0]/r,  p*pos[1]/r,  p*pos[2]/r, energy) ;
  }
  else {
    AliInfo("Null cluster radius, momentum calculation not possible");
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, const AliVCluster* vc)
{
  if (fMassHypothesis > 0) return GetMomentum(mom, vc, fMassHypothesis);

  if (vc) {
    if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
      vc->GetMomentum(mom, fVertex, (AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
    }
    else {
      vc->GetMomentum(mom, fVertex);
    }
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
  return kTRUE;
}

/**
 * Get momentum of the i^th particle in array
 * @param mom
 * @param i
 * @return
 */
Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  AliVCluster *vc = GetCluster(i);
  return GetMomentum(mom, vc);
}

/**
 * Get momentum of the next particle in array
 * @param mom
 * @return
 */
Bool_t AliClusterContainer::GetNextMomentum(TLorentzVector &mom)
{
  AliVCluster *vc = GetNextCluster();
  return GetMomentum(mom, vc);
}

/**
 * Get momentum of the i^th particle in array
 * @param mom
 * @param i
 * @return
 */
Bool_t AliClusterContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  AliVCluster *vc = GetAcceptCluster(i);
  return GetMomentum(mom, vc);
}

/**
 * Get momentum of the next accepted particle in array
 * @param mom
 * @return
 */
Bool_t AliClusterContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  AliVCluster *vc = GetNextAcceptCluster();
  return GetMomentum(mom, vc);
}

Bool_t AliClusterContainer::AcceptCluster(Int_t i)
{
  Bool_t r = ApplyClusterCuts(GetCluster(i));
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom);
}

Bool_t AliClusterContainer::AcceptCluster(const AliVCluster* clus)
{
  Bool_t r = ApplyClusterCuts(clus);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, clus);

  return ApplyKinematicCuts(mom);
}

/**
 * Return true if cluster is accepted.
 * @param clus
 * @return
 */
Bool_t AliClusterContainer::ApplyClusterCuts(const AliVCluster* clus)
{
  fRejectionReason = 0;

  if (!clus) {
    fRejectionReason |= kNullObject;
    return kFALSE;
  }
      
  if (!clus->IsEMCAL()) {
    fRejectionReason |= kIsEMCalCut;
    return kFALSE;
  }

  if (clus->TestBits(fBitMap) != (Int_t)fBitMap) {
    fRejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(clus->GetLabel()) > fMinMCLabel) {
    fRejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(clus->GetLabel()) < fMaxMCLabel) {
    fRejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (clus->GetTOF() > fClusTimeCutUp || clus->GetTOF() < fClusTimeCutLow) {
    fRejectionReason |= kTimeCut;
    return kFALSE;
  }

  if (fExoticCut && clus->GetIsExotic()) {
    fRejectionReason |= kExoticCut;
    return kFALSE;
  }
  
  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    if (clus->GetUserDefEnergy((VCluUserDefEnergy_t)i) < fUserDefEnergyCut[i]) {
      fRejectionReason |= kEnergyCut;
      return kFALSE;
    }
  }
  
  return kTRUE;
}

/**
 * Get number of accepted particles
 * @return
 */
Int_t AliClusterContainer::GetNAcceptedClusters()
{
  Int_t nClus = 0;
  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliVCluster *clus = GetNextAcceptCluster();
  if(clus) nClus = 1;
  while (GetNextAcceptCluster())
    nClus++;

  fCurrentID = tempID;

  return nClus;
}

Double_t AliClusterContainer::GetClusUserDefEnergyCut(Int_t t) const
{
  if (t >= 0 && t <= AliVCluster::kLastUserDefEnergy){
    return fUserDefEnergyCut[t];
  }
  else {
    return fMinE;
  }
}

void AliClusterContainer::SetClusUserDefEnergyCut(Int_t t, Double_t cut)
{
  if (t >= 0 && t <= AliVCluster::kLastUserDefEnergy){
    fUserDefEnergyCut[t] = cut;
  }
  else {
    fMinE = cut;
  }
}


/**
 * Set the class name
 * @param clname
 */
void AliClusterContainer::SetClassName(const char *clname)
{
  TClass cls(clname);
  if (cls.InheritsFrom("AliVCluster")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliClusterContainer, it must inherits from AliVCluster!",clname));
}

const char* AliClusterContainer::GetTitle() const
{
  static TString clusterString;

  Double_t Ecut = GetClusUserDefEnergyCut(GetDefaultClusterEnergy());

  if (Ecut == 0) {
    clusterString = TString::Format("%s_E0000", GetArrayName().Data());
  }
  else if (Ecut < 1.0) {
    clusterString = TString::Format("%s_E0%3.0f", GetArrayName().Data(), Ecut*1000.0);
  }
  else {
    clusterString = TString::Format("%s_E%4.0f", GetArrayName().Data(), Ecut*1000.0);
  }

  return clusterString.Data();
}
