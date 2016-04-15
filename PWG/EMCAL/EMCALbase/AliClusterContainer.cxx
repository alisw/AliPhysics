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
#include <algorithm>
#include <iostream>
#include <vector>
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
  fBaseClassName = "AliVCluster";
  SetClassName("AliVCluster");

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
  fBaseClassName = "AliVCluster";
  SetClassName("AliVCluster");

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
AliVCluster* AliClusterContainer::GetAcceptClusterWithLabel(Int_t lab) const
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

/**
 * Get momentum of the i^th particle in array
 * @param mom
 * @param i
 * @return
 */
Bool_t AliClusterContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
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
Bool_t AliClusterContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
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

/**
 * Return true if cluster is accepted.
 * @param clus
 * @return
 */
Bool_t AliClusterContainer::ApplyClusterCuts(const AliVCluster* clus, UInt_t &rejectionReason) const
{
  if (!clus) {
    rejectionReason |= kNullObject;
    return kFALSE;
  }
      
  if (!clus->IsEMCAL()) {
    rejectionReason |= kIsEMCalCut;
    return kFALSE;
  }

  if (clus->TestBits(fBitMap) != (Int_t)fBitMap) {
    rejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(clus->GetLabel()) > fMinMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(clus->GetLabel()) < fMaxMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (clus->GetTOF() > fClusTimeCutUp || clus->GetTOF() < fClusTimeCutLow) {
    rejectionReason |= kTimeCut;
    return kFALSE;
  }

  if (fExoticCut && clus->GetIsExotic()) {
    rejectionReason |= kExoticCut;
    return kFALSE;
  }
  
  for (Int_t i = 0; i <= AliVCluster::kLastUserDefEnergy; i++) {
    if (clus->GetUserDefEnergy((VCluUserDefEnergy_t)i) < fUserDefEnergyCut[i]) {
      rejectionReason |= kEnergyCut;
      return kFALSE;
    }
  }
  
  return kTRUE;
}

/**
 * Get number of accepted particles
 * @return
 */
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

/**
 * Create standard (forward) iterator over all clusters in the
 * container, for the beginning of the iteration.
 * @return Standard (forward) beginning iterator
 */
AliClusterContainer::all_iterator AliClusterContainer::begin() const{
  return all_iterator(this, 0, true);
}

/**
 * Create standard (forward) iterator over all clusters in the
 * container, marking the end of the iteration.
 * @return Standard (forward) end iterator
 */
AliClusterContainer::all_iterator AliClusterContainer::end() const{
  return all_iterator(this, this->GetNClusters(), true);
}

/**
 * Create standard (forward) iterator over all clusters in the
 * container, for the beginning of the iteration.
 * @return Reverse (backward) beginning operator
 */
AliClusterContainer::all_iterator AliClusterContainer::rbegin() const{
  return all_iterator(this, this->GetNClusters()-1, false);
}

/**
 * Create reverse (backward) iterator over all clusters in the
 * container, marking the end of the iteration.
 * @return Reverse (backward) end iterator
 */
AliClusterContainer::all_iterator AliClusterContainer::rend() const{
  return all_iterator(this, -1, false);
}

/**
 * Create standard (forward) iterator over accepted clusters in the
 * container, for the beginning of the iteration.
 * @return Standard (forward) beginning iterator
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_begin() const{
  return accept_iterator(this, 0, true);
}

/**
 * Create standard (forward) iterator over accepted clusters in the
 * container, marking the end of the iteration.
 * @return Standard (forward) end iterator
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_end() const{
  return accept_iterator(this, this->GetNAcceptedClusters(), true);
}

/**
 * Create standard (forward) iterator over accepted clusters in the
 * container, for the beginning of the iteration.
 * @return Reverse (backward) beginning operator
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_rbegin() const{
  return accept_iterator(this, this->GetNAcceptedClusters()-1, false);
}

/**
 * Create reverse (backward) iterator over accepted clusters in the
 * container, marking the end of the iteration.
 * @return Reverse (backward) end iterator
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_rend() const{
  return accept_iterator(this, -1, false);
}


///////////////////////////////////////////////////////////////
///
/// Content of accept_iterator
///
///////////////////////////////////////////////////////////////

/**
 * Constructor, inializing the iterator with
 * - Cluster container to iterate over
 * - Starting position for the iteration
 * - Direction of the iteration (default: forward)
 * This constructor also builds a lookup table of
 * accepted clusters.
 *
 * Note: Iterators should be created by the object
 * it iterates over and not by hand.
 *
 * @param[in] cont Container to be iterated over (not modified)
 * @param[in] start Start position for the iteration
 * @param[in] forward If true iterator is a forward iterator
 */
AliClusterContainer::accept_iterator::accept_iterator(
    const AliClusterContainer *cont,
    int start,
    bool forward):
    fkContainer(cont),
    fAcceptIndices(),
    fCurrentPos(start),
    fForward(forward)
{
  fAcceptIndices.Set(fkContainer->GetNAcceptedClusters());
  UInt_t rejection = 0;
  Int_t mappos = 0;
  for(int i = 0; i < fkContainer->GetNEntries(); i++){
    if(fkContainer->AcceptCluster(i, rejection)){
      fAcceptIndices[mappos++] = i;
    }
  }
}

/**
 * Copy constructor. For the underlying container the
 * pointer is copied.
 * @param[in] other Reference for the copy
 */
AliClusterContainer::accept_iterator::accept_iterator(const accept_iterator &other):
    fkContainer(other.fkContainer),
    fAcceptIndices(other.fAcceptIndices),
    fCurrentPos(other.fCurrentPos),
    fForward(other.fForward)
{
}

/**
 * Assignmet operator, copying properties of the
 * reference iterator
 * @param[in] other Reference for the assignment
 * @return This object after the assingment
 */
AliClusterContainer::accept_iterator &AliClusterContainer::accept_iterator::operator=(const AliClusterContainer::accept_iterator &other){
  if(this != &other){
    fkContainer = other.fkContainer;
    fAcceptIndices = other.fAcceptIndices;
    fCurrentPos = other.fCurrentPos;
    fForward = other.fForward;
  }
  return *this;
}

/**
 * Comparison operator. The comparison is performed
 * based on the current position of the iterator.
 * @param[in] other Object to compare to.
 * @return If true the positions of the two iterators match
 */
bool AliClusterContainer::accept_iterator::operator!=(const AliClusterContainer::accept_iterator &other) const{
  return fCurrentPos != other.fCurrentPos;
}

/**
 * Incrementation prefix operator. Incrementing/Decrementing
 * the postion of the iterator in the array based on the direction.
 * @return State of the iterator (reference) after iteration.
 */
AliClusterContainer::accept_iterator &AliClusterContainer::accept_iterator::operator++(){
  if(fForward){
    fCurrentPos++;
  } else {
    fCurrentPos--;
  }
  return *this;
}

/**
 * Incrementation postfix operator. Incrementing/Decrementing
 * the postion of the iterator in the array based on the direction.
 * @param[in] Not used, only defining the operator as postfix operator
 * @return State of the iterator (copy) before iteration
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_iterator::operator++(int){
  AliClusterContainer::accept_iterator result(*this);
  this->operator++();
  return result;
}

/**
 * Decrementation prefix operator. Decrementing/Incrementing
 * the postion of the iterator in the array based on the direction.
 * @return State of the iterator (reference) after iteration.
 */
AliClusterContainer::accept_iterator &AliClusterContainer::accept_iterator::operator--(){
  if(fForward){
    fCurrentPos--;
  } else {
    fCurrentPos++;
  }
  return *this;
}

/**
 * Decrementation postfix operator. Decrementing/Incrementing
 * the postion of the iterator in the array based on the direction.
 * @param[in] Not used, only defining the operator as postfix operator
 * @return State of the iterator (copy) before iteration
 */
AliClusterContainer::accept_iterator AliClusterContainer::accept_iterator::operator--(int){
  AliClusterContainer::accept_iterator result(*this);
  this->operator++();
  return result;
}

/**
 * Dereferencing operator. Providing access to the accepted cluster
 * at the current position.
 * @return Accepted cluster at the given position (NULL if iterator is out of range)
 */
AliVCluster *AliClusterContainer::accept_iterator::operator*() const{
  if(fCurrentPos < 0 || fCurrentPos >= fAcceptIndices.GetSize()){
    return NULL;
  }
  return fkContainer->GetCluster(fAcceptIndices[fCurrentPos]);
}

///////////////////////////////////////////////////////////////
///
/// Content of all_iterator
///
///////////////////////////////////////////////////////////////

/**
 * Constructor, inializing the iterator with
 * - Cluster container to iterate over
 * - Starting position for the iteration
 * - Direction of the iteration (default: forward)
 *
 * Note: Iterators should be created by the object
 * it iterates over and not by hand.
 *
 * @param[in] cont Container to be iterated over (not modified)
 * @param[in] startpos Start position for the iteration
 * @param[in] forward If true iterator is a forward iterator
 */
AliClusterContainer::all_iterator::all_iterator(
    const AliClusterContainer *cont,
    int startpos,
    bool forward):
    fkContainer(cont),
    fCurrentPos(startpos),
    fForward(forward)
{
}

/**
 * Copy constructor. For the underlying container the
 * pointer is copied.
 * @param[in] other Reference for the copy
 */
AliClusterContainer::all_iterator::all_iterator(const all_iterator &other):
    fkContainer(other.fkContainer),
    fCurrentPos(other.fCurrentPos),
    fForward(other.fForward)
{
}

/**
 * Assignmet operator, copying properties of the
 * reference iterator
 * @param[in] other Reference for the assignment
 * @return This object after the assingment
 */
AliClusterContainer::all_iterator &AliClusterContainer::all_iterator::operator=(const AliClusterContainer::all_iterator &other){
  if(&other != this){
    fkContainer = other.fkContainer;
    fCurrentPos = other.fCurrentPos;
    fForward = other.fForward;
  }
  return *this;
}

/**
 * Comparison operator. The comparison is performed
 * based on the current position of the iterator.
 * @param[in] other Object to compare to.
 * @return If true the positions of the two iterators match
 */
bool AliClusterContainer::all_iterator::operator!=(const AliClusterContainer::all_iterator &other) const{
  return other.fCurrentPos != fCurrentPos;
}

/**
 * Incrementation prefix operator. Incrementing/Decrementing
 * the postion of the iterator in the array based on the direction.
 * @return State of the iterator (reference) after iteration.
 */
AliClusterContainer::all_iterator &AliClusterContainer::all_iterator::operator++(){
  if(fForward){
    fCurrentPos++;
  } else {
    fCurrentPos--;
  }
  return *this;
}

/**
 * Incrementation postfix operator. Incrementing/Decrementing
 * the postion of the iterator in the array based on the direction.
 * @param[in] Not used, only defining the operator as postfix operator
 * @return State of the iterator (copy) before iteration
 */
AliClusterContainer::all_iterator AliClusterContainer::all_iterator::operator++(int){
  all_iterator tmp(*this);
  operator++();
  return tmp;
}

/**
 * Decrementation prefix operator. Decrementing/Incrementing
 * the postion of the iterator in the array based on the direction.
 * @return State of the iterator (reference) after iteration.
 */
AliClusterContainer::all_iterator &AliClusterContainer::all_iterator::operator--(){
  if(fForward){
    fCurrentPos--;
  } else {
    fCurrentPos++;
  }
  return *this;
}

/**
 * Decrementation postfix operator. Decrementing/Incrementing
 * the postion of the iterator in the array based on the direction.
 * @param[in] Not used, only defining the operator as postfix operator
 * @return State of the iterator (copy) before iteration
 */
AliClusterContainer::all_iterator AliClusterContainer::all_iterator::operator--(int){
  all_iterator tmp(*this);
  operator--();
  return tmp;
}

/**
 * Dereferencing operator. Accessing cluster at the current
 * position of the iterator in the cluster container
 * @return
 */
AliVCluster *AliClusterContainer::all_iterator::operator*() const{
  if(fCurrentPos >= 0 && fCurrentPos < fkContainer->GetNClusters())
    return fkContainer->GetCluster(fCurrentPos);
  return NULL;
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
    for(AliClusterContainer::accept_iterator iter = cont->accept_begin(); iter != cont->accept_end(); ++iter){
      variation.push_back(*iter);
    }
  } else {
    // test all iterator
    for(AliClusterContainer::all_iterator iter = cont->begin(); iter != cont->end(); ++iter){
      variation.push_back(*iter);
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
