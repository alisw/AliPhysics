/************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <vector>
#include <TClonesArray.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliTLorentzVector.h"

#include "AliClusterContainer.h"

ClassImp(AliClusterContainer);

// string to enum map for use with the %YAML config
const std::map <std::string, AliVCluster::VCluUserDefEnergy_t> AliClusterContainer::fgkClusterEnergyTypeMap = {
  {"kNonLinCorr", AliVCluster::kNonLinCorr },
  {"kHadCorr", AliVCluster::kHadCorr },
  {"kUserDefEnergy1", AliVCluster::kUserDefEnergy1 },
  {"kUserDefEnergy2", AliVCluster::kUserDefEnergy2 }
};

// Properly instantiate the object
AliEmcalContainerIndexMap <TClonesArray, AliVCluster> AliClusterContainer::fgEmcalContainerIndexMap;

AliClusterContainer::AliClusterContainer():
  AliEmcalContainer(),
  fEMCALCells(nullptr),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1),
  fIncludePHOS(kFALSE),
  fIncludePHOSonly(kFALSE),
  fPhosMinNcells(0),
  fPhosMinM02(0),
  fEmcalMinNcells(0),
  fEmcalMinM02(-1.),
  fEmcalMaxM02(DBL_MAX),
  fEmcalMaxM02CutEnergy(DBL_MAX),
  fMaxFracEnergyLeadingCell(DBL_MAX),
  fSelectOneCellCluster(kTRUE),
  fEmcalMaxNCellEffCutEnergy(4.)
{
  fBaseClassName = "AliVCluster";
  SetClassName("AliVCluster");

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}

AliClusterContainer::AliClusterContainer(const char *name):
  AliEmcalContainer(name),
  fEMCALCells(nullptr),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fExoticCut(kTRUE),
  fDefaultClusterEnergy(-1),
  fIncludePHOS(kFALSE),
  fIncludePHOSonly(kFALSE),
  fPhosMinNcells(0),
  fPhosMinM02(0),
  fEmcalMinNcells(0),
  fEmcalMinM02(-1.),
  fEmcalMaxM02(DBL_MAX),
  fEmcalMaxM02CutEnergy(DBL_MAX),
  fMaxFracEnergyLeadingCell(DBL_MAX),
  fSelectOneCellCluster(kTRUE),
  fEmcalMaxNCellEffCutEnergy(4.)
{
  fBaseClassName = "AliVCluster";
  SetClassName("AliVCluster");

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    fUserDefEnergyCut[i] = 0.;
  }
}


AliVCluster* AliClusterContainer::GetLeadingCluster(const char* opt)
{
  TString option(opt);
  option.ToLower();

  double (AliTLorentzVector::*momentum)() const = 0;

  if (option.Contains("e")) {
    momentum = &AliTLorentzVector::E;
  }
  else {
    momentum = &AliTLorentzVector::Et;
  }

  AliClusterIterableMomentumContainer::momentum_object_pair clusterMax;

  for (auto cluster : accepted_momentum()) {
    if ((clusterMax.first.*momentum)() < (cluster.first.*momentum)()) {
      clusterMax = cluster;
    }
  }

  return clusterMax.second;
}

AliVCluster* AliClusterContainer::GetCluster(Int_t i) const 
{
  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVCluster *vp = static_cast<AliVCluster*>(fClArray->At(i));
  return vp;

}

AliVCluster* AliClusterContainer::GetAcceptCluster(Int_t i) const
{
  AliVCluster *vc = GetCluster(i);
  if (!vc) return 0;

  UInt_t rejectionReason = 0;
  if (AcceptCluster(vc, rejectionReason))
    return vc;
  else {
    AliDebug(2,"Cluster not accepted.");
    return 0;
  }
}

AliVCluster* AliClusterContainer::GetClusterWithLabel(Int_t lab) const 
{
  Int_t i = GetIndexFromLabel(lab);
  return GetCluster(i);
}

AliVCluster* AliClusterContainer::GetAcceptClusterWithLabel(Int_t lab) const
{
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptCluster(i);
}

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

Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, const AliVCluster* vc, Double_t mass) const
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

Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, const AliVCluster* vc) const
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

Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  AliVCluster *vc = GetCluster(i);
  return GetMomentum(mom, vc);
}

Bool_t AliClusterContainer::GetNextMomentum(TLorentzVector &mom)
{
  AliVCluster *vc = GetNextCluster();
  return GetMomentum(mom, vc);
}

Bool_t AliClusterContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  AliVCluster *vc = GetAcceptCluster(i);
  return GetMomentum(mom, vc);
}

Bool_t AliClusterContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  AliVCluster *vc = GetNextAcceptCluster();
  return GetMomentum(mom, vc);
}

Bool_t AliClusterContainer::AcceptCluster(Int_t i, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyClusterCuts(GetCluster(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliClusterContainer::AcceptCluster(const AliVCluster* clus, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyClusterCuts(clus, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, clus);

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliClusterContainer::ApplyClusterCuts(const AliVCluster* clus, UInt_t &rejectionReason) const
{
  if (!clus) {
    rejectionReason |= kNullObject;
    return kFALSE;
  }

  Bool_t bInAcceptance = clus->IsEMCAL();
  if (fIncludePHOS) {
    bInAcceptance = clus->IsEMCAL() || (clus->GetType() == AliVCluster::kPHOSNeutral);
  }
  if (fIncludePHOSonly) {
    bInAcceptance = (clus->GetType() == AliVCluster::kPHOSNeutral);
  }
  if (!bInAcceptance) {
      rejectionReason |= kIsEMCalCut;
      return kFALSE;
  }

  if (clus->TestBits(fBitMap) != (Int_t)fBitMap) {
    rejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(clus->GetLabel()) < fMinMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(clus->GetLabel()) > fMaxMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (clus->GetTOF() > fClusTimeCutUp || clus->GetTOF() < fClusTimeCutLow) {
    rejectionReason |= kTimeCut;
    return kFALSE;
  }
  // GetPassedSpecialNCell disabling implicit cut on the number of cells
  // for flagged 1-cell clusters
  if (fExoticCut && clus->GetIsExotic() && !(fSelectOneCellCluster && GetPassedSpecialNCell(clus))) {
    rejectionReason |= kExoticCut;
    return kFALSE;
  }

  if (clus->IsEMCAL()) {
    if(clus->GetNCells() < fEmcalMinNcells && !(fSelectOneCellCluster && GetPassedSpecialNCell(clus))) {
      // explicit Ncell cut
      rejectionReason |= kExoticCut;
      return kFALSE;
    }
    if (clus->GetNonLinCorrEnergy() < fEmcalMaxM02CutEnergy) {
      if((clus->GetM02() < fEmcalMinM02 || clus->GetM02() > fEmcalMaxM02) && !(fSelectOneCellCluster && GetPassedSpecialNCell(clus))) {
        rejectionReason |= kExoticCut; // Not really true, but there is a lack of leftover bits
        return kFALSE;
      }
    }

  }

  if(fEMCALCells) {
    double ecellmax = 0;
    for(int icell = 0; icell < clus->GetNCells(); icell++){
      double celltmp = fEMCALCells->GetCellAmplitude(clus->GetCellAbsId(icell));
      if(celltmp > ecellmax) ecellmax = celltmp;
    }
    if(ecellmax/clus->E() > fMaxFracEnergyLeadingCell && !(fSelectOneCellCluster && GetPassedSpecialNCell(clus))) {
      rejectionReason |= kExoticCut;
      return kFALSE;
    }
  }

  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    if (clus->GetUserDefEnergy((VCluUserDefEnergy_t)i) < fUserDefEnergyCut[i]) {
      rejectionReason |= kEnergyCut;
      return kFALSE;
    }
  }

  if (fIncludePHOS || fIncludePHOSonly) {
    if (clus->GetType() == AliVCluster::kPHOSNeutral) {
      if (clus->GetNCells() < fPhosMinNcells) {
        rejectionReason |= kExoticCut;
        return kFALSE;
      }

      if (clus->GetM02() < fPhosMinM02) {
        rejectionReason |= kExoticCut;
        return kFALSE;
      }
    }
  }

  return kTRUE;
}

Bool_t AliClusterContainer::GetPassedSpecialNCell(const AliVCluster* clus) const
{
  if(!clus->IsEMCAL()) return kFALSE;
  if(clus->GetNCells() > 1) return kFALSE;
  if(clus->Chi2() == 1 && clus->E() < fEmcalMaxNCellEffCutEnergy) return kTRUE;
  return kFALSE;
}

Int_t AliClusterContainer::GetNAcceptedClusters() const
{
  UInt_t rejectionReason = 0;
  Int_t nClus = 0;
  for(int iclust = 0; iclust < this->fClArray->GetEntries(); ++iclust){
    AliVCluster *clust = this->GetCluster(iclust);
    if(this->AcceptCluster(clust, rejectionReason)) nClus++;
  }
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

void AliClusterContainer::SetArray(const AliVEvent * event)
{
  AliEmcalContainer::SetArray(event);

  // Connect EMCAL cells
  fEMCALCells = event->GetEMCALCells();

  // Register TClonesArray in index map
  fgEmcalContainerIndexMap.RegisterArray(GetArray());
}

const AliClusterIterableContainer AliClusterContainer::all() const {
  return AliClusterIterableContainer(this, false);
}

const AliClusterIterableContainer AliClusterContainer::accepted() const {
  return AliClusterIterableContainer(this, true);
}

const AliClusterIterableMomentumContainer AliClusterContainer::all_momentum() const {
  return AliClusterIterableMomentumContainer(this, false);
}

const AliClusterIterableMomentumContainer AliClusterContainer::accepted_momentum() const {
  return AliClusterIterableMomentumContainer(this, true);
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

TString AliClusterContainer::GetDefaultArrayName(const AliVEvent * const ev) const {
  if(ev->IsA() == AliAODEvent::Class()) return "caloClusters";
  else if(ev->IsA() == AliESDEvent::Class()) return "CaloClusters";
  else return "";
}


/******************************************
 * Unit tests                             *
 ******************************************/

int TestClusterContainerIterator(const AliClusterContainer *const cont, int iteratorType, bool verbose){
  std::vector<AliVCluster *> reference, variation;
  AliVCluster *test = NULL;
  for(int iclust = 0; iclust < cont->GetNClusters(); iclust++){
    test = cont->GetCluster(iclust);
    if(!iteratorType){
      UInt_t rejectionReason = 0;
      if(!cont->AcceptCluster(test, rejectionReason)) continue;
    }
    reference.push_back(test);
  }

  if(!iteratorType){
    // test accept iterator
    for(auto cluster : cont->accepted()){
      variation.push_back(cluster);
    }
  } else {
    // test all iterator
    for(auto cluster : cont->all()){
      variation.push_back(cluster);
    }
  }

  if(verbose){
    // Printing cluster addresses:
    if(reference.size() < 30){
      std::cout << "Clusters in reference container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVCluster *>::iterator refit = reference.begin(); refit != reference.end(); ++refit){
        std::cout << "Address: " << *refit << std::endl;
      }
    }
    if(variation.size() < 30){
      std::cout << "Clusters in test container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVCluster *>::iterator varit = variation.begin(); varit != variation.end(); ++varit){
        std::cout << "Address: " << *varit << std::endl;
      }
    }
  }

  int testresult = 0;
  // compare distributions: all clusters in one vector needs to be found in the other vector and vice versa
  bool failure = false;
  // first test: cleck if all clusters are found by the iterator
  for(std::vector<AliVCluster *>::iterator clit = reference.begin(); clit != reference.end(); ++clit){
    if(std::find(variation.begin(), variation.end(), *clit) == variation.end()) {
      if(verbose)
        std::cout << "Could not find cluster with address " << *clit << " in test container" << std::endl;
      failure = true;
      break;
    }
  }
  if(!failure){
    // second test: check if there are no extra clusters found
    for(std::vector<AliVCluster *>::iterator clit = variation.begin(); clit != variation.end(); ++clit){
      if(std::find(reference.begin(), reference.end(), *clit) == reference.end()) {
        if(verbose)
          std::cout << "Could not find cluster with address " << *clit << " in reference container" << std::endl;
        failure = true;
        break;
      }
    }
    if(failure) testresult = 2;
  } else {
    testresult = 1;
  }

  if(verbose){
    std::cout << "Unit test cluster container, iterator type " << iteratorType << std::endl;
    std::cout << "Number of expected clusters: " << reference.size() << ", number of found clusters: " << variation.size() << std::endl;
    std::cout << "Test result: " << testresult << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
  }
  return testresult;
}
