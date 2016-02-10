//
// Container with name, TClonesArray and cuts for calo clusters
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "TLorentzVector.h"

#include "AliClusterContainer.h"

ClassImp(AliClusterContainer)

//________________________________________________________________________
AliClusterContainer::AliClusterContainer():
  AliEmcalContainer(),
  fClusPtCut(0.),
  fClusECut(0.15),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1)
{
  // Default constructor.

  fClassName = "AliVCluster";

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}

//________________________________________________________________________
AliClusterContainer::AliClusterContainer(const char *name):
  AliEmcalContainer(name),
  fClusPtCut(0.),
  fClusECut(0.15),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1)
{
  // Standard constructor.

  fClassName = "AliVCluster";

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetLeadingCluster(const char* opt)
{
  // Get the leading cluster; use e if "e" is contained in opt (otherwise et)

  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;

  AliVCluster *clusterMax = GetNextAcceptCluster(0);
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

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetCluster(Int_t i) const 
{
  //Get i^th cluster in array

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVCluster *vp = static_cast<AliVCluster*>(fClArray->At(i));
  return vp;

}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetAcceptCluster(Int_t i) 
{
  //Return pointer to cluster if cluster is accepted

  AliVCluster *vc = GetCluster(i);
  if (!vc) return 0;

  if (AcceptCluster(vc))
    return vc;
  else {
    AliDebug(2,"Cluster not accepted.");
    return 0;
  }
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetClusterWithLabel(Int_t lab) const 
{
  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetCluster(i);
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetAcceptClusterWithLabel(Int_t lab)  
{
  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptCluster(i);
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetNextAcceptCluster(Int_t i) 
{
  //Get next accepted cluster; if i >= 0 (re)start counter from i; return 0 if no accepted cluster could be found

  if (i>=0) fCurrentID = i;

  const Int_t n = GetNEntries();
  AliVCluster *c = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    c = GetAcceptCluster(fCurrentID);
  } while (!c);

  return c;
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetNextCluster(Int_t i) 
{
  //Get next cluster; if i >= 0 (re)start counter from i; return 0 if no cluster could be found

  if (i>=0) fCurrentID = i;

  const Int_t n = GetNEntries();
  AliVCluster *c = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    c = GetCluster(fCurrentID);
  } while (!c);

  return c;
}

//________________________________________________________________________
Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVCluster *vc = GetCluster(i);
  if (vc) {
    if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
      vc->GetMomentum(mom, fVertex, (AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
    }
    else {
      vc->GetMomentum(mom, fVertex);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliClusterContainer::GetNextMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVCluster *vc = GetNextCluster(i);
  if (vc) {
    if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
      vc->GetMomentum(mom, fVertex, (AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
    }
    else {
      vc->GetMomentum(mom, fVertex);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliClusterContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVCluster *vc = GetAcceptCluster(i);
  if (vc) {
    if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
      vc->GetMomentum(mom, fVertex, (AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
    }
    else {
      vc->GetMomentum(mom, fVertex);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliClusterContainer::GetNextAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVCluster *vc = GetNextAcceptCluster(i);
  if (vc) {
    if (fDefaultClusterEnergy >= 0 &&  fDefaultClusterEnergy <= AliVCluster::kLastUserDefEnergy) {
      vc->GetMomentum(mom, fVertex, (AliVCluster::VCluUserDefEnergy_t)fDefaultClusterEnergy);
    }
    else {
      vc->GetMomentum(mom, fVertex);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliClusterContainer::AcceptCluster(AliVCluster *clus)
{
  // Return true if cluster is accepted.

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

  if (clus->E()<fClusECut) {
    fRejectionReason |= kEnergyCut;
    return kFALSE;
  }

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fClusPtCut) {
    fRejectionReason |= kPtCut;
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

//________________________________________________________________________
Int_t AliClusterContainer::GetNAcceptedClusters()
{
  // Get number of accepted particles

  Int_t nClus = 0;

  AliVCluster *clus = GetNextAcceptCluster(0);
  if(clus) nClus = 1;
  while (GetNextAcceptCluster())
    nClus++;

  return nClus;
}

//________________________________________________________________________
Double_t AliClusterContainer::GetClusUserDefEnergyCut(Int_t t) const
{
  if (t >= 0 && t <= AliVCluster::kLastUserDefEnergy){
    return fUserDefEnergyCut[t];
  }
  else {
    return fClusECut;
  }
}

//________________________________________________________________________
void AliClusterContainer::SetClusUserDefEnergyCut(Int_t t, Double_t cut)
{
  if (t >= 0 && t <= AliVCluster::kLastUserDefEnergy){
    fUserDefEnergyCut[t] = cut;
  }
  else {
    fClusECut = cut;
  }
}


//________________________________________________________________________
void AliClusterContainer::SetClassName(const char *clname)
{
  // Set the class name

  TClass cls(clname);
  if (cls.InheritsFrom("AliVCluster")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliClusterContainer, it must inherits from AliVCluster!",clname));
}

//________________________________________________________________________
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
