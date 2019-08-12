/************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                 *
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
#include "AliAnalysisTaskRhoAverage.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "TLorentzVector.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"

ClassImp(AliAnalysisTaskRhoAverage)

AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRhoAverage"),
  fRhoType(0),
  fNExclLeadPart(0),
  fUseMedian(kFALSE),
  fTotalArea(1)
{
}

AliAnalysisTaskRhoAverage::AliAnalysisTaskRhoAverage(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fRhoType(0),
  fNExclLeadPart(0),
  fUseMedian(kFALSE),
  fTotalArea(1)
{
}

Bool_t AliAnalysisTaskRhoAverage::Run() 
{
  const Int_t NMAX = 9999;
  static Double_t rhovec[NMAX];
  Int_t NpartAcc = 0;

  Int_t   maxPartIds[] = {0, 0};
  Float_t maxPartPts[] = {0, 0};

  // push all jets within selected acceptance into stack

  AliParticleContainer* tracks = GetParticleContainer(0);
  AliClusterContainer* clusters = GetClusterContainer(0);

  if (fNExclLeadPart > 0) {

    if (tracks && (fRhoType == 0 || fRhoType == 1)) {

      AliVParticle *track = 0;
      tracks->ResetCurrentID();
      while ((track = tracks->GetNextAcceptParticle())) {

        if (track->Pt() > maxPartPts[0]) {
          maxPartPts[1] = maxPartPts[0];
          maxPartIds[1] = maxPartIds[0];
          maxPartPts[0] = track->Pt();
          maxPartIds[0] = tracks->GetCurrentID()+1;
        }
        else if (track->Pt() > maxPartPts[1]) {
          maxPartPts[1] = track->Pt();
          maxPartIds[1] = tracks->GetCurrentID()+1;
        }
      }
    }

    if (clusters && (fRhoType == 0 || fRhoType == 2)) {

      AliVCluster *cluster = 0;
      clusters->ResetCurrentID();
      while ((cluster = clusters->GetNextAcceptCluster())) {
        TLorentzVector nPart;
        clusters->GetMomentum(nPart, clusters->GetCurrentID());

        if (nPart.Pt() > maxPartPts[0]) {
          maxPartPts[1] = maxPartPts[0];
          maxPartIds[1] = maxPartIds[0];
          maxPartPts[0] = nPart.Pt();
          maxPartIds[0] = -clusters->GetCurrentID()-1;
        }
        else if (nPart.Pt() > maxPartPts[1]) {
          maxPartPts[1] = nPart.Pt();
          maxPartIds[1] = -clusters->GetCurrentID()-1;
        }
      }
    }

    if (fNExclLeadPart < 2) {
      maxPartIds[1] = 0;
      maxPartPts[1] = 0;
    }
  }

  if (tracks && (fRhoType == 0 || fRhoType == 1)) {
    AliVParticle *track = 0;
    tracks->ResetCurrentID();
    while ((track = tracks->GetNextAcceptParticle()) && NpartAcc < NMAX) {

      // exlcuding lead particles
      if (tracks->GetCurrentID() == maxPartIds[0]-1 || tracks->GetCurrentID() == maxPartIds[1]-1)
        continue;

      rhovec[NpartAcc] = track->Pt();
      ++NpartAcc;
    }
  }

  if (clusters && (fRhoType == 0 || fRhoType == 2)) {

    AliVCluster *cluster = 0;
    clusters->ResetCurrentID();
    while ((cluster = clusters->GetNextAcceptCluster()) && NpartAcc < NMAX) {
      // exlcuding lead particles
      if (clusters->GetCurrentID() == -maxPartIds[0]-1 || clusters->GetCurrentID() == -maxPartIds[1]-1)
        continue;

      TLorentzVector nPart;
      clusters->GetMomentum(nPart, clusters->GetCurrentID());

      rhovec[NpartAcc] = nPart.Pt();
      ++NpartAcc;
    }
  }

  if (NpartAcc == NMAX) {
    AliError(Form("%s: NpartAcc >= %d", GetName(), NMAX));
  }

  Double_t rho = 0;

  if (NpartAcc > 0) {
    if (fUseMedian)
      rho = TMath::Median(NpartAcc, rhovec);
    else
      rho = TMath::Mean(NpartAcc, rhovec);

    rho *= NpartAcc / fTotalArea;
  }

  fOutRho->SetVal(rho);

  if (fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fOutRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}

void AliAnalysisTaskRhoAverage::ExecOnce() 
{
  AliAnalysisTaskRhoBase::ExecOnce();

  AliParticleContainer *partCont = GetParticleContainer(0);
  if (!partCont) {
    AliError(Form("%s: No particle container found! Assuming area = 1...",GetName()));
    fTotalArea = 1;
    return;
  }

  Float_t maxEta = partCont->GetParticleEtaMax();
  Float_t minEta = partCont->GetParticleEtaMin();
  Float_t maxPhi = partCont->GetParticlePhiMax();
  Float_t minPhi = partCont->GetParticlePhiMin();
  
  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;
  
  fTotalArea = (maxEta - minEta) * (maxPhi - minPhi);

  if (fTotalArea < 1e-6) {
    AliError(Form("%s: Area = %f < 1e-6, assuming area = 1", GetName(), fTotalArea));
    fTotalArea = 1;
  }
}
