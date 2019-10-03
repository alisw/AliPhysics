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
#include "AliEmcalParticleJetConstituent.h"
#include "AliVParticle.h"

/// \cond CLASSIMP
ClassImp(PWG::JETFW::AliEmcalParticleJetConstituent)
/// \endcond

namespace PWG {
namespace JETFW {

AliEmcalParticleJetConstituent::AliEmcalParticleJetConstituent() :
    AliEmcalJetConstituent(),
    fkParticle(nullptr)
{

}

AliEmcalParticleJetConstituent::AliEmcalParticleJetConstituent(const AliVParticle *const part) :
    AliEmcalJetConstituent(),
    fkParticle(part)
{

}

AliEmcalParticleJetConstituent::AliEmcalParticleJetConstituent(const AliEmcalParticleJetConstituent &other) :
    AliEmcalJetConstituent(other),
    fkParticle(other.fkParticle)
{

}

AliEmcalParticleJetConstituent &AliEmcalParticleJetConstituent::operator=(const AliEmcalParticleJetConstituent &other){
  AliEmcalJetConstituent::operator=(other);
  if(this != &other){
    fkParticle = other.fkParticle;
  }
  return *this;
}

bool AliEmcalParticleJetConstituent::operator ==(const AliEmcalParticleJetConstituent &rhs) const {
  if(!fkParticle || !rhs.fkParticle) return false;
  return rhs.fkParticle->GetUniqueID() == fkParticle->GetUniqueID();
}

bool AliEmcalParticleJetConstituent::operator >(const AliEmcalParticleJetConstituent &rhs) const {
  if(!fkParticle || !rhs.fkParticle) return false;
  return rhs.fkParticle->Pt() < fkParticle->Pt();
}

bool AliEmcalParticleJetConstituent::operator <(const AliEmcalParticleJetConstituent &rhs) const {
  if(!fkParticle || !rhs.fkParticle) return false;
  return rhs.fkParticle->Pt() > fkParticle->Pt();
}


AliEmcalParticleJetConstituent::~AliEmcalParticleJetConstituent() {

}

double AliEmcalParticleJetConstituent::Px() const {
  return fkParticle->Px();
}

double AliEmcalParticleJetConstituent::Py() const {
  return fkParticle->Py();
}

double AliEmcalParticleJetConstituent::Pz() const {
  return fkParticle->Pz();
}

double AliEmcalParticleJetConstituent::Pt() const {
  return fkParticle->Pt();
}

double AliEmcalParticleJetConstituent::E() const {
  return fkParticle->E();
}

double AliEmcalParticleJetConstituent::Eta() const {
  return fkParticle->Eta();
}

double AliEmcalParticleJetConstituent::Phi() const {
  return fkParticle->Phi();
}

}
}
