/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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
#include "AliEmcalClusterJetConstituent.h"

/// \cond CLASSIMP
ClassImp(PWG::JETFW::AliEmcalClusterJetConstituent)
/// \endcond

namespace PWG {
namespace JETFW {

AliEmcalClusterJetConstituent::AliEmcalClusterJetConstituent() :
    AliEmcalJetConstituent(),
    fkCaloCluster(nullptr),
    fDefaultEnergyDefinition(AliVCluster::kUserDefEnergy1),
    fPVec()
{
}


AliEmcalClusterJetConstituent::AliEmcalClusterJetConstituent(const AliVCluster *const clust, AliVCluster::VCluUserDefEnergy_t energydef, Double_t *pvec):
    AliEmcalJetConstituent(),
    fkCaloCluster(clust),
    fDefaultEnergyDefinition(energydef),
    fPVec(pvec)
{
}

AliEmcalClusterJetConstituent::AliEmcalClusterJetConstituent(const AliVCluster *const clust):
    AliEmcalJetConstituent(),
    fkCaloCluster(clust),
    fDefaultEnergyDefinition(AliVCluster::kUserDefEnergy1),
    fPVec()
{
}

AliEmcalClusterJetConstituent::AliEmcalClusterJetConstituent(const AliEmcalClusterJetConstituent &other) :
    AliEmcalJetConstituent(other),
    fkCaloCluster(other.fkCaloCluster),
    fDefaultEnergyDefinition(other.fDefaultEnergyDefinition),
    fPVec(other.fPVec)
{
}

AliEmcalClusterJetConstituent &AliEmcalClusterJetConstituent::operator=(const AliEmcalClusterJetConstituent &other){
  AliEmcalJetConstituent::operator=(other);

  if(this != &other){
    fkCaloCluster = other.fkCaloCluster;
    fDefaultEnergyDefinition = other.fDefaultEnergyDefinition;
    fPVec = other.fPVec;
  }
  return *this;
}

bool AliEmcalClusterJetConstituent::operator==(const AliEmcalClusterJetConstituent &rhs) const{
  if(!fkCaloCluster || !rhs.fkCaloCluster) return false;
  return fkCaloCluster->GetUniqueID() == rhs.fkCaloCluster->GetUniqueID();
}


AliEmcalClusterJetConstituent::~AliEmcalClusterJetConstituent() {

}

double AliEmcalClusterJetConstituent::Px() const{
  return fPVec.X();
}

double AliEmcalClusterJetConstituent::Py() const{
  return fPVec.Y();
}

double AliEmcalClusterJetConstituent::Pz() const {
  return fPVec.Z();
}

double AliEmcalClusterJetConstituent::Pt() const {
  return fPVec.Pt();
}

double AliEmcalClusterJetConstituent::E() const {
  return fkCaloCluster ? fkCaloCluster->GetUserDefEnergy(fDefaultEnergyDefinition): 0;
}

double AliEmcalClusterJetConstituent::Eta() const {
  return fPVec.Eta();
}

double AliEmcalClusterJetConstituent::Phi() const {
  return fPVec.Phi();
}

}
}
