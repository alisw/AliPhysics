/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef __ALIEMCALMCPARTONINFO_H__
#define __ALIEMCALMCPARTONINFO_H__

#include <TNamed.h>
#include <TList.h>
#include "AliTLorentzVector.h"

class AliVParticle;

namespace PWG {

namespace EMCAL {

class AliEmcalPartonData : public TObject {
  public:
    AliEmcalPartonData(Int_t pdg, Int_t positionInStack, Double_t px, Double_t py, Double_t pz, Double_t e) :
      TObject(),
      fPdg(pdg),
      fPositionInStack(positionInStack),
      fMomentumVec(px, py, pz, e)
    {
    }
    virtual ~AliEmcalPartonData() { }

    void SetPdg(int pdg) { fPdg = pdg; }
    void SetPositionInStack(int position) { fPositionInStack = position; }
    void SetMomentum(Double_t px, Double_t py, Double_t pz, Double_t e) { fMomentumVec.SetPxPyPzE(px, py, pz, e); }

    int GetPdg() const { return fPdg; }
    int GetPositionInStack() const { return fPdg; }
    const AliTLorentzVector &GetMomentum() const { return fMomentumVec; }

  private:
    Int_t   fPdg;                             ///< PDG code
    Int_t   fPositionInStack;                 ///< Position of the particle in stack
    AliTLorentzVector fMomentumVec;           ///< Momentum vector
};


class AliEmcalMCPartonInfo : public TNamed {
public:
  AliEmcalMCPartonInfo();
  AliEmcalMCPartonInfo(const char *name);
  virtual ~AliEmcalMCPartonInfo() {}
  virtual void Clear(Option_t *option="") { Reset(); }

  void AddDirectParton(const AliVParticle *part, Int_t positionInStack);

  const TList &GetListOfDirectPartons() const { return fPartons; }
  AliEmcalPartonData *GetHardestParton() const { return fHardestParton; }

  void Reset();

private:
  TList               fPartons;           //!<! List of partons directly coming from beam particles
  AliEmcalPartonData *fHardestParton;     //!<! Hardest direct parton

  AliEmcalMCPartonInfo(const AliEmcalMCPartonInfo &);
  AliEmcalMCPartonInfo &operator=(const AliEmcalMCPartonInfo &);

  ClassDef(AliEmcalMCPartonInfo, 1);
};

}

}

#endif