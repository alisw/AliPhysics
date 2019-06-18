/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TH1F.h>
#include <TH2F.h>
#include <TFormula.h>
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliCSTrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"

/// \file AliCSTrackCuts.cxx
/// \brief Implementation of event cuts class within the correlation studies analysis

const char *AliCSTrackCuts::fgkCutsNames[AliCSTrackCuts::kNCuts] = {
    "Not constrained",
    "#eta",
    "Track type (ITS cls,TPC cls,DCA)",
    "p_{t}",
};


/// Default constructor for serialization
AliCSTrackCuts::AliCSTrackCuts() :
    AliCSTrackCutsBase(),
    fConstrain(kFALSE),
    fBaseSystem(kUnknownBase),
    fEtaCut(100.0),
    fMinPt(0.0),
    fMaxPt(9999.0),
    fEtaShift(0.0),
    fESDTrackCuts(NULL),
    fAODFilterBits(0),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhPtVsDCAxy{NULL},
    fhPtVsDCAz{NULL},
    fhPtVsTPCCls{NULL},
    fhPtVsTPCRowOverFindCls{NULL},
    fhEtaVsPhi{NULL},
    fhPtVsEta{NULL}
{
}

/// Constructor
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSTrackCuts::AliCSTrackCuts(const char *name, const char *title) :
    AliCSTrackCutsBase(kNCuts,kNCutsParameters,name,title),
    fConstrain(kFALSE),
    fBaseSystem(kUnknownBase),
    fEtaCut(100.0),
    fMinPt(0.0),
    fMaxPt(9999.0),
    fEtaShift(0.0),
    fESDTrackCuts(NULL),
    fAODFilterBits(0),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhPtVsDCAxy{NULL},
    fhPtVsDCAz{NULL},
    fhPtVsTPCCls{NULL},
    fhPtVsTPCRowOverFindCls{NULL},
    fhEtaVsPhi{NULL},
    fhPtVsEta{NULL}
{
}

/// Destructor
/// We don't own anything, everything we allocate is owned
/// by the output list
AliCSTrackCuts::~AliCSTrackCuts()
{

}

/// Processes a potential change in the run number
///
/// Checks if the current period under analysis has changes and if so
/// updates the needed members
void AliCSTrackCuts::NotifyRun() {

  /* checks the change in the analysis period */
  if (AliCSTrackCutsBase::GetGlobalPeriod() != fDataPeriod) {

    fDataPeriod = AliCSTrackCutsBase::GetGlobalPeriod();

    /* Configure the cuts accordingly to the data period */
    this->SetActualTypeOfTrackCuts();
    this->SetActualITSClustersCut();
    this->SetActualTPCClustersCut();
    this->SetActualDCACut();

    /* and now we ask for histogram allocation */
    DefineHistograms();
  }
}

/// Processes the start of a new event
///
/// Does nothin for the time being
void AliCSTrackCuts::NotifyEvent() {

}

/// Check whether the passed track is accepted by the different cuts
/// \param trk the track to analyze whether it is accepted or not
/// \return kTRUE if the track  is accepted, kFALSE otherwise
///
/// An internal copy of the track is used in case it were needed to constrain it if the cut
/// so decide it. The (constrained) copy is then used for further analysis or cut
/// decisions. The #fConstrain flag is kept to inform the user she needs to constrain the track.
Bool_t AliCSTrackCuts::IsTrackAccepted(AliVTrack *trk) {

  AliVTrack *ttrk = NULL; /* it will be created once the type is known */
  fConstrain = kFALSE; /* it will be updated with the track type */

  /* just to be sure */
  if (trk == NULL) return kFALSE;

  /* if MC analysis is ongoing get rid of "ghost" tracks */
  if (fgIsMC) if (trk->GetLabel() < 0) return kFALSE;

  /* for the time being */
  Bool_t accepted = kTRUE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask.ResetAllBits();

  /* is the constrained parameter info available */
  if (trk->IsA() == AliESDtrack::Class()) {
    if (!trk->GetConstrainedParam()) {
      fCutsActivatedMask.SetBitNumber(kNotConstrainedCut);
      accepted = kFALSE;
    }
  }

  /* Track type cut: TPC clusters, ITS clusters and DCA cuts. We will need to change this */
  /* should we need more granularity in the output information */
  /* The track will be constrained if the cut so decide it then the (potentially) constrained */
  /* is returned in ttrk for further use */
  if (!AcceptTrackType(trk,ttrk)) {
    fCutsActivatedMask.SetBitNumber(kTrackTypeCuts);
    accepted = kFALSE;
  }

  /* update the constrained potentially constrained track */
  if (!fConstrain) {
    /* not constrained so use the same track further on */
    ttrk = trk;
  }

  /* eta cut */
  if (fCutsEnabledMask.TestBitNumber(kEtaCut)) {
    if ((ttrk->Eta() < (- fEtaCut + fEtaShift)) ||
        ((fEtaCut + fEtaShift) < ttrk->Eta())) {
      fCutsActivatedMask.SetBitNumber(kEtaCut);
      accepted = kFALSE;
    }
  }

  /* Pt cut */
  if (fCutsEnabledMask.TestBitNumber(kPtCut)) {
    if ((ttrk->Pt() < fMinPt) || (fMaxPt < ttrk->Pt())) {
      fCutsActivatedMask.SetBitNumber(kPtCut);
      accepted = kFALSE;
    }
  }

  if (fQALevel > kQALevelNone) {
    /* let's fill the histograms */
    fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n tracks")));
    if (!accepted)
      fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n cut tracks")));

    for (Int_t i=0; i<kNCuts; i++) {
      if (fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i]) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", i, fgkCutsNames[i]));

      if (fCutsActivatedMask.TestBitNumber(i))
        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i])));

      if (fQALevel > kQALevelLight) {
        for (Int_t j=i; j<kNCuts; j++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[j]) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", j, fgkCutsNames[j]));

          if (fCutsActivatedMask.TestBitNumber(i) && fCutsActivatedMask.TestBitNumber(j)) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fgkCutsNames[i]));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fgkCutsNames[j]));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }

    /* get some values needed for histograms filling */
    Float_t dca[2], bCov[3];
    if (trk->IsA() == AliESDtrack::Class()) {
      ttrk->GetImpactParameters(dca,bCov);
    }
    else {
      /* TODO: if the track is constrained this needs to be considered */
      /* the selected tracks get a DCA of 0.0. Implement GetDCA from   */
      /* AliDielectronVarManager but it will require access to the     */
      /* event object                                                  */
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);
      Float_t pos[3];
      if (aodt->GetPosition(pos)) {
        dca[0] = pos[0];
        dca[1] = pos[1];
      }
      else {
        dca[0] = 0.0;
        dca[1] = 0.0;
      }
    }

    Float_t nCrossedRowsTPC = ttrk->GetTPCCrossedRows();
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (ttrk->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / ttrk->GetTPCNclsF();
    }

    for (Int_t i = 0; i < 2; i++) {

      fhPtVsDCAxy[i]->Fill(dca[0],ttrk->Pt());
      fhPtVsDCAz[i]->Fill(dca[1],ttrk->Pt());
      fhPtVsTPCCls[i]->Fill(ttrk->GetTPCNcls(),ttrk->Pt());
      fhPtVsTPCRowOverFindCls[i]->Fill(ratioCrossedRowsOverFindableClustersTPC,ttrk->Pt());
      fhPtVsEta[i]->Fill(ttrk->Eta(),ttrk->Pt());

      if (fQALevel > kQALevelLight) {
        fhEtaVsPhi[i]->Fill(ttrk->Phi()*180.0/TMath::Pi(),ttrk->Eta());
      }

      /* don't fill after if event not accepted */
      if (!accepted) break;
    }
  }
  /* delete the potentially constrained track */
  if (fConstrain) {
    /* constrained so use delete it */
    delete ttrk;
  }

  return accepted;
}


/// Check whether the true track associated to the passed track is accepted by the kinematic cuts
/// \param trk the track to analyze whether its associated true track is accepted or not
/// \return kTRUE if the associated true track  is accepted, kFALSE otherwise
///
Bool_t AliCSTrackCuts::IsTrueTrackAccepted(AliVTrack *trk) {

  /* reject ghost tracks */
  if (trk->GetLabel() < 0) return kFALSE;

  return IsTrueTrackAccepted(trk->GetLabel());
}


/// Check whether the passed true track is accepted by the kinematic and PID cuts
/// \param itrk the index of true track to analyze whether it is accepted or not
/// \return kTRUE if the track  is accepted, kFALSE otherwise
///
Bool_t AliCSTrackCuts::IsTrueTrackAccepted(Int_t itrk) {

  Double_t eta = 0.0;
  Double_t pt = 0.0;

  if (fgIsESD) {
    /* we stay only with primary tracks */
    if (!IsPhysicalPrimary(itrk))
      return kFALSE;

    /* get the associated particle */
    AliVParticle *particle = fgMCHandler->MCEvent()->GetTrack(itrk);

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    eta = particle->Eta();
    pt = particle->Pt();
  }
  else {
    /* get the associated particle */
    AliAODMCParticle *particle = (AliAODMCParticle *) fgMCArray->At(itrk);

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    /* we stay only with primary tracks */
    if (!particle->IsPhysicalPrimary())
      return kFALSE;

    eta = particle->Eta();
    pt = particle->Pt();
  }


  /* eta cut */
  if (fCutsEnabledMask.TestBitNumber(kEtaCut)) {
    if ((eta < (- fEtaCut + fEtaShift)) ||
        ((fEtaCut + fEtaShift) < eta)) {
      return kFALSE;
    }
  }

  /* Pt cut */
  if (fCutsEnabledMask.TestBitNumber(kPtCut)) {
    if ((pt < fMinPt) || (fMaxPt < pt)) {
      return kFALSE;
    }
  }

  return kTRUE;
}

/// Check whether the passed track is accepted by the track type cuts
/// \param trk the track to analyze whether it is accepted or not
/// \param ttrk the potential constrained track. NULL if not constrained
/// \return kTRUE if the track  is accepted, kFALSE otherwise
///
/// If the selected cut requires the track being constrained, the passed track is
/// cloned, its clone constrained, its reference returned in *ttrk* parameter and the
/// #fConstrain flag raised.
Bool_t AliCSTrackCuts::AcceptTrackType(AliVTrack *trk, AliVTrack *&ttrk) {
  Bool_t accepted = kTRUE;

  if (trk->IsA() == AliESDtrack::Class()) {
    if (!fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack *>(trk))) {
      accepted = kFALSE;
    }
    else {
      /* the track is accepted, check special handling depending on the cut */
      switch (fParameters[kTrackTypeCutParam]) {
      case 5: /* special handling for TPC only tracks to constrain to the SPD vertex */
        AliFatal("Functionality not yet supported");
        break;
      case 7: /* special handling for global hybrid constrained */
          /* we have to select complementary hybrid tracks but this is period dependent */
        AliFatal("Functionality not yet supported");
        switch(fBaseSystem) {
        case k2010based:
          /* to check if HG1 or GC we just require ITS:none which complements case 6 */
          if (trk->GetStatus() & AliESDtrack::kITSupg) {
            /* upgraded ITS with seven layers */
            if (trk->HasPointOnITSLayer(0) || (trk->HasPointOnITSLayer(1) && trk->HasPointOnITSLayer(2))) {
              accepted = kFALSE;
            }
            AliError("We should not be here within 2010 based periods");
          }
          else {
            if (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1)) {
              accepted = kFALSE;
            }
          }
          break;
        case k2011based:
          break;
        default:
          accepted = kFALSE;
        }
        /* so, if we reached here with the track accepted we have to constrain it */
        if (accepted) {
          /* constrain the track */
          if (trk->GetConstrainedParam()) {
            AliESDtrack *tmptrk = new AliESDtrack(*((AliESDtrack*) trk));
            tmptrk->Set(trk->GetConstrainedParam()->GetX(),
                trk->GetConstrainedParam()->GetAlpha(),
                trk->GetConstrainedParam()->GetParameter(),
                trk->GetConstrainedParam()->GetCovariance());
            ttrk = tmptrk;
            fConstrain = kTRUE;
          }
          else
            accepted = kFALSE;
        }
        break;
      default:
        ;
      }
    }
  }
  else if (trk->IsA() == AliAODTrack::Class()) {
    /* TODO: for the time being only the pure equivalence of FB with the */
    /* type of track selected is used. This means that DCA cut, ITS      */
    /* clusters cut and TPC clusters cut are not considered apart of     */
    /* what is already implicit in the FB. As an additional comment to   */
    /* consider, if the cut is loose than the FB, the FB is useless      */
    /* to discard, tests with the track information available
    Int_t a,bb,c,d;
    Int_t clsITS[6];
    Float_t b[2];
    Float_t bCov[3];

    trk->GetTPCClusterInfo(a,bb,c,d);
    trk->GetTPCCrossedRows();
    trk->GetTPCNcls();
    trk->GetTPCNclsF();
    trk->GetTPCSharedMapPtr();
    trk->GetImpactParameters(b,bCov);
    trk->GetITSclusters(clsITS);

    end tests */

    AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);

    if(!aodt->TestFilterBit(fAODFilterBits)) {
      accepted = kFALSE;
    }
    else {
      /* the track is accepted, check special handling depending on the cut */
      switch (fParameters[kTrackTypeCutParam]) {
      case 5: /* special handling for TPC only tracks to constrain to the SPD vertex */
        /* do nothing, the track should have already been constrained */
        break;
      case 7: /* special handling for global hybrid constrained */
          /* we have to select complementary hybrid tracks but this is period dependent */
        AliFatal("Functionality not yet supported");
        switch(fBaseSystem) {
        case k2010based:
          /* to check if HG1 or GC we just require ITS:none which complements case 6 */
          if (trk->GetStatus() & AliESDtrack::kITSupg) {
            /* upgraded ITS with seven layers */
            if (trk->HasPointOnITSLayer(0) || (trk->HasPointOnITSLayer(1) && trk->HasPointOnITSLayer(2))) {
              accepted = kFALSE;
            }
            AliError("We should not be here within 2010 based periods");
          }
          else {
            if (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1)) {
              accepted = kFALSE;
            }
          }
          break;
        case k2011based:
          break;
        default:
          accepted = kFALSE;
        }
        /* so, if we reached here with the track accepted we have to constrain it */
        /* TODO: this is wrong. If we are in AOD the tracks should already be constrained */
        if (accepted) {
          /* constrain the track */
          if (trk->GetConstrainedParam()) {
            AliFatal("Functionality not yet supported");
            fConstrain = kTRUE;
          }
          else
            accepted = kFALSE;
        }
        break;
      default:
        ;
      }
    }
  }
  else {
    AliError("Neither ESD nor AOD track!!!");
    accepted = kFALSE;
  }
  return accepted;
}


/// Check whether the index passed corresponds to a MC physical primary track
///
/// Basically this is needed in case of AMPT fast MC which does not discard
/// weak decays when asked for physical primaries.
///
/// (NOTE: learned from AliAnalysisTaskPhiCorrelations)
/// \param itrk the index of the true track
/// \return kTRUE if the track particle is MC physical primary kFALSE otherwise
Bool_t AliCSTrackCuts::IsPhysicalPrimary(Int_t itrk) {

  if (fgIsESD) {
    if (!fgIsMConlyTruth) {
      return fgMCHandler->MCEvent()->IsPhysicalPrimary(itrk);
    }
    else {
      /* taken from AliAnalysisTaskPhiCorrelations */
      // Exclude weak decay products (if not done by IsPhysicalPrimary)
      // In order to prevent analyzing daughters from weak decays
      // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it

      const Int_t kNWeakParticles = 7;
      const Int_t kWeakParticles[kNWeakParticles] = {
          3322, 3312, 3222, // Xi0 Xi+- Sigma-+
          3122, 3112, // Lambda0 Sigma+-
          130, 310 // K_L0 K_S0
      };

      AliMCEvent* mcEvent = fgMCHandler->MCEvent();

      if (mcEvent != NULL) {
        /* if it is rejected by the own event we have finished */
        if (!mcEvent->IsPhysicalPrimary(itrk))
          return kFALSE;

        AliVParticle* particle = mcEvent->GetTrack(itrk);

        if (particle != NULL) {
          Int_t motherix = particle->GetMother();

          if (motherix < 0) {
            return kTRUE;
          }
          else {
            AliVParticle* motherparticle = mcEvent->GetTrack(motherix);

            if (motherparticle != NULL) {
              Int_t pdgcode = TMath::Abs(motherparticle->PdgCode());

              for (Int_t j=0; j != kNWeakParticles; ++j) {
                if (kWeakParticles[j] == pdgcode) {
                  return kFALSE;
                }
              }
              return kTRUE;
            }
            else {
              return kTRUE;
            }
          }
        }
        else {
          return kFALSE;
        }
      }
      else {
        return kFALSE;
      }
    }
  }
  else {
    /* get the associated particle */
    AliAODMCParticle *particle = (AliAODMCParticle *) fgMCArray->At(TMath::Abs(itrk));

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    /* we stay only with primary tracks */
    return particle->IsPhysicalPrimary();
  }
}


/// Check whether the true particle associated to a reconstructed track is primary
/// \param trk the reconstructed track
/// \return kTRUE if the associated particle is primary kFALSE otherwise
Bool_t AliCSTrackCuts::IsTruePrimary(AliVTrack *trk) {

  if (fgIsESD) {
    return IsPhysicalPrimary(TMath::Abs(trk->GetLabel()));
  }
  else {
    /* get the associated particle */
    AliAODMCParticle *particle = (AliAODMCParticle *) fgMCArray->At(TMath::Abs(trk->GetLabel()));

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    /* we stay only with primary tracks */
    return particle->IsPhysicalPrimary();
  }
}


/// Sets the individual value for the cut parameter ID
/// \param paramID the ID of the cut parameter of interest
/// \param value the value for the cut parameter
/// \return kTRUE if the cut parameter value was accepted
Bool_t AliCSTrackCuts::SetCutAndParams(Int_t paramID, Int_t value) {

  switch (cutsParametersIds(paramID)) {
  case kTrackTypeCutParam:
    if (SetTypeOfTrackCut(value)) {
      fParameters[kTrackTypeCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kEtaCutParam:
    if (SetEtaCut(value)) {
      fParameters[kEtaCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSClsCutParam:
    if (SetITSClustersCut(value) ) {
      fParameters[kITSClsCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCClsCutParam:
    if (SetTPCClustersCut(value)) {
      fParameters[kTPCClsCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kDCACutParam:
    if (SetDCACut(value)) {
      fParameters[kDCACutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kPtCutParam:
    if (SetPtCut(value)) {
      fParameters[kPtCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  default:
    AliError(Form("Cut param id %d out of supported range", paramID));
    return kFALSE;
  }
}


/// Print the whole cut iformatio for the cut ID
/// \param paramID the ID of the cut of interest
void AliCSTrackCuts::PrintCutWithParams(Int_t paramID) const {

  switch (cutsParametersIds(paramID)) {
  case kTrackTypeCutParam:
    printf("  Type of track cut: ");
    this->PrintTypeOfTrackCut();
    break;
  case kEtaCutParam:
    printf("  Acceptance cut: ");
    if (fCutsEnabledMask.TestBitNumber(kEtaCut))
      printf("|eta| < %3.1f\n",fEtaCut);
    else
      printf("none\n");
    break;
  case kITSClsCutParam:
    printf("  ITS clusters cut: ");
    PrintITSClustersCut();
    break;
  case kTPCClsCutParam:
    printf("  TPC clusters cut: ");
    PrintTPCClustersCut();
    break;
  case kDCACutParam:
    printf("  DCA cut: ");
    PrintDCACut();
    break;
  case kPtCutParam:
    printf("  Pt cut: ");
    if (fCutsEnabledMask.TestBitNumber(kPtCut))
      printf("%3.1f < pT %s", fMinPt, (fMaxPt < 9990) ? Form("%3.1f\n", fMaxPt) : "\n");
    else
      printf("none\n");
    break;
  default:
    AliError(Form("Cut param id %d out of supported range", paramID));
  }
}

/// Configures the type of track cut
/// \param ttype The type of track the cut shall select
/// | code | type of track |
/// |:--:|:--|
/// | 0 | default standard cut for the considered data period |
/// | 1 | FB1. TPC only tracks |
/// | 2 | FB32. Global primaries |
/// | 3 | FB64. Global primaries complementary to code 2 ones |
/// | 4 | FB16. Global primaries plus secondaries |
/// | 5 | FB128. TPC only tracks to constrain to the SPD vertex |
/// | 6 | FB256. Global hybrid primaries plus secondaries |
/// | 7 | FB512. To constrain global hybrid primaries plus secondaries complementary to code 6 |
/// | 9 | FBxFFFF. no cut on the type of track, i.e. accept all of them |
/// \return kTRUE if proper and supported \f$ \eta \f$  cut code
///
Bool_t AliCSTrackCuts::SetTypeOfTrackCut(Int_t ttype) {
  switch (ttype) {
  case 0:
    break;
  case 1:
    break;
  case 2:
    break;
  case 3:
    break;
  case 4:
    break;
  case 5:
    break;
  case 6:
    break;
  case 7:
    break;
  case 9:
    break;
  default:
    AliError(Form("Track type cut code %d not supported", ttype));
    return kFALSE;
  }
  return kTRUE;
}

/// Prints the type of track cut
void AliCSTrackCuts::PrintTypeOfTrackCut() const {
  switch(fParameters[kTrackTypeCutParam]){
    case 0:
      printf("default standard cut for the considered data period\n");
      break;
    case 1:
      printf("FB1. TPC only tracks\n");
      break;
    case 2:
      printf("FB32. Global primaries. Tight DCA cut. A hit on SPD:ANY\n");
      break;
    case 3:
      printf("FB64. Global primaries. Tight DCA cut. No hit on SPD\n");
      break;
    case 4:
      printf("FB16. Global primaries plus secondaries. Loose DCA cut. A hit on SPD:ANY\n");
      break;
    case 5:
      printf("FB128. TPC only tracks to constrain to the SPD vertex\n");
      break;
    case 6:
      printf("FB256. Global hybrid primaries plus secondaries. Loose DCA cut. A hit on SPD:ANY. TPC shared clusters\n");
      break;
    case 7:
      printf("FB512. Global hybrid primaries plus secondaries to constrain. Loose DCA cut. No hit on SPD or ITS refit not required. TPC shared clusters\n");
      break;
    case 9:
      printf("No track type cut. Accept all\n");
      break;
    default:
      AliError(Form("Wrong track type cut code %d stored", fParameters[kTrackTypeCutParam]));
  }
}

/// Set the standard track cuts set according to the data period
void AliCSTrackCuts::SetActualTypeOfTrackCuts() {

  if (fESDTrackCuts != NULL)
    delete fESDTrackCuts;
  fAODFilterBits = 0;

  TString system = "";
  TString period = "";
  TString basename = "";
  fBaseSystem = kUnknownBase;


  switch (GetGlobalAnchorPeriod()) {
  case kLHC10bg:
    fBaseSystem = k2010based;
    basename = "2010";
    system = "p-p";
    period = "2010bg";
    break;
  case kLHC11a:
  case kLHC11b:
  case kLHC11cg:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "p-p";
    period = "2011ag";
    break;
  case kLHC12:
  case kLHC13g:
  case kLHC15fm:
  case kLHC15n:
  case kLHC16k:
  case kLHC16l:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "p-p";
    period = "various";
    break;
  case kLHC13bc:
  case kLHC13de:
  case kLHC13f:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "p-Pb";
    period = "2013bf";
    break;
  case kLHC10h:
    fBaseSystem = k2010based;
    basename = "2010";
    system = "Pb-Pb";
    period = "2010h";
    break;
  case kLHC11h:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2011h";
    break;
  case kLHC15oLIR:
  case kLHC15oHIR:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2015o";
    break;
  case kLHC17n:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "Xe-Xe";
    period = "2017n";
    break;
  case kLHC18q:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2018q";
    break;
  case kLHC18r:
    fBaseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2018r";
    break;
  default:
    fESDTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fAODFilterBits = 1;
    fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard TPC only");
    AliError("SYSTEM: no system set from data files. Standard TPC only track cuts");
    return;
  }

  TString tracktype = "";
  switch (fParameters[kTrackTypeCutParam]) {
  case 0:
    /* printf("default standard cut for the considered data period\n"); */
    switch (fBaseSystem) {
    case k2010based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard 2010 primaries");
      fAODFilterBits = 32;
      tracktype = "FB32. Standard ITS+TPC 2010. Tight DCA. ITS:ANY. Number of clusters";
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard 2011 primaries");
      fAODFilterBits = 32;
      tracktype = "FB32. Standard ITS+TPC 2011. Tight DCA. ITS:ANY. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 1:
    /* printf("TPC only tracks\n"); */
    fESDTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard TPC only");
    fAODFilterBits = 1;
    tracktype = "FB1. Standard TPC only";
    break;
  case 2:
    /* printf("Global primaries. Tight DCA cut. A hit on SPD:ANY\n"); */
    switch (fBaseSystem) {
    case k2010based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2010 primaries");
      fAODFilterBits = 32;
      tracktype = "FB32. Global primaries ITS+TPC 2010. Tight DCA. SPD:ANY. Number of clusters";
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2011 primaries");
      fAODFilterBits = 32;
      tracktype = "FB32. Global primaries ITS+TPC 2011. Tight DCA. SPD:ANY. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 3:
    /* printf("Global primaries. Tight DCA cut. No hit on SPD\n"); */
    switch (fBaseSystem) {
    case k2010based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global SPD:NONE 2010 primaries");
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
      fAODFilterBits = 64;
      tracktype = "FB64. Global primaries ITS+TPC 2010. Tight DCA. SPD:NONE, SDD:FIRST. Number of clusters";
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global SPD:NONE 2011 primaries");
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
      fAODFilterBits = 64;
      tracktype = "FB64. Global primaries ITS+TPC 2011. Tight DCA. SPD:NONE, SDD:FIRST. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 4:
    /* printf("Global primaries plus secondaries. Loose DCA cut. A hit on SPD:ANY\n"); */
    switch (fBaseSystem) {
    case k2010based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2010 primaries and secondaries");
      fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
      fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
      fESDTrackCuts->SetDCAToVertex2D(kTRUE);
      fAODFilterBits = 16;
      tracktype = "FB16. Global primaries plus secondaries ITS+TPC 2010. Loose DCA. ITS:ANY. Number of clusters";
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2011 primaries and secondaries");
      fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
      fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
      fESDTrackCuts->SetDCAToVertex2D(kTRUE);
      fAODFilterBits = 16;
      tracktype = "FB16. Global primaries plus secondaries ITS+TPC 2011. Loose DCA. ITS:ANY. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 5:
    /* printf("TPC only tracks to constrain to the SPD vertex\n"); */
    switch (fBaseSystem) {
    case k2010based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard TPC only, to constrain to SPD vertex");
      fESDTrackCuts->SetMinNClustersTPC(70);
      fAODFilterBits = 128;
      tracktype = "FB128. Standard TPC only, to constrain to SPD vertex";
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Standard TPC only, to constrain to SPD vertex");
      fAODFilterBits = 128;
      tracktype = "FB128. Standard TPC only, to constrain to SPD vertex";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 6:
    /* printf("Global hybrid primaries plus secondaries. Loose DCA cut. A hit on SPD:ANY. TPC shared clusters\n"); */
    switch (fBaseSystem) {
    case k2010based: {
        fESDTrackCuts = new AliESDtrackCuts(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2010 hybrid tracks");

        TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
        fESDTrackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
        fESDTrackCuts->SetMinNClustersTPC(70);
        fESDTrackCuts->SetMaxChi2PerClusterTPC(4);
        fESDTrackCuts->SetRequireTPCStandAlone(kTRUE);
        fESDTrackCuts->SetAcceptKinkDaughters(kFALSE);
        fESDTrackCuts->SetRequireTPCRefit(kTRUE);
        fESDTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
        // ITS
        fESDTrackCuts->SetRequireITSRefit(kTRUE);
        fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
        //accept secondaries
        fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
        fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
        fESDTrackCuts->SetDCAToVertex2D(kTRUE);
        //reject fakes
        fESDTrackCuts->SetMaxChi2PerClusterITS(36);
        fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

        fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);

        fESDTrackCuts->SetEtaRange(-0.9,0.9);
        fESDTrackCuts->SetPtRange(0.15);

        fAODFilterBits = 256;
        tracktype = "Global hybrid 2010 track. Loose DCA. ITS:ANY. Number of clusters";
      }
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2011 hybrid tracks");
      fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
      fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
      fESDTrackCuts->SetDCAToVertex2D(kTRUE);
      fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      fESDTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
      fAODFilterBits = 256;
      tracktype = "Global hybrid 2011 tracks. Loose DCA. ITS:ANY. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 7:
    /* printf("Global hybrid primaries plus secondaries to constrain. Loose DCA cut. No hit on SPD or ITS refit not required. TPC shared clusters\n"); */
    switch (fBaseSystem) {
    case k2010based: {
        fESDTrackCuts = new AliESDtrackCuts(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2010 hybrid tracks to constrain");

        TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
        fESDTrackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
        fESDTrackCuts->SetMinNClustersTPC(70);
        fESDTrackCuts->SetMaxChi2PerClusterTPC(4);
        fESDTrackCuts->SetRequireTPCStandAlone(kTRUE);
        fESDTrackCuts->SetAcceptKinkDaughters(kFALSE);
        fESDTrackCuts->SetRequireTPCRefit(kTRUE);
        fESDTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
        // ITS
        /* the full handling requires some additional processing at cut time */
        fESDTrackCuts->SetRequireITSRefit(kFALSE);
        fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
        //accept secondaries
        fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
        fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
        fESDTrackCuts->SetDCAToVertex2D(kTRUE);
        //reject fakes
        fESDTrackCuts->SetMaxChi2PerClusterITS(36);
        fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

        fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);

        fESDTrackCuts->SetPtRange(0.15);

        fAODFilterBits = 512;
        tracktype = "FB 512. Global hybrid 2010 track. Loose DCA. ITS:OFF no ITS refit required. Requires processing. Number of clusters";
      }
      break;
    case k2011based:
      fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      fESDTrackCuts->SetNameTitle(Form("TrackType_%s",GetCutsString()),"Type of tracks: Global 2011 hybrid tracks to constrain");
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
      fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
      fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
      fESDTrackCuts->SetDCAToVertex2D(kTRUE);
      fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      fESDTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
      fAODFilterBits = 512;
      tracktype = "FB 512. Global hybrid 2011 tracks. Loose DCA. ITS:none. Number of rows";
      break;
    default:
      AliFatal("Internal inconsistency between data period and base period for track cuts selection. ABORTING!!!");
      return;
    }
    break;
  case 9:
    fESDTrackCuts = new AliESDtrackCuts(Form("TrackType_%s", GetCutsString()), "Type of tracks: all");
    fAODFilterBits = 0xFFFF;
    printf("No track type cut. Accept all\n");
    break;
  default:
    AliError(Form("Stored track type cut code %d not supported", fParameters[kTrackTypeCutParam]));
  }

  /* report the selected cut */
  AliInfo("=============== Track type cut ===========================");
  AliInfo(Form("SYSTEM: %s; PERIOD: %s", system.Data(), period.Data()));
  AliInfo(Form("TRACK CUT: %s", tracktype.Data()));
  AliInfo(Form("BASED ON: %s standards cuts", basename.Data()));
  AliInfo("=============== Track type cut end =======================");
}

/// Configures the \f$ \eta \f$ cut
/// \param etacode the \f$ \eta \f$ cut code
/// | code | \f$ \eta \f$ cut|
/// |:--:|:--|
/// | 0 | no \f$ \eta \f$ cut|
/// | 1 | \f$ -1.4 < \eta < 1.4\f$ |
/// | 2 | \f$ -1.2 < \eta < 1.2\f$ |
/// | 3 | \f$ -1.0 < \eta < 1.0\f$ |
/// | 4 | \f$ -0.9 < \eta < 0.9\f$ |
/// | 5 | \f$ -0.8 < \eta < 0.8\f$ |
/// | 6 | \f$ -0.7 < \eta < 0.7\f$ |
/// | 7 | \f$ -0.6 < \eta < 0.6\f$ |
/// | 8 | \f$ -0.5 < \eta < 0.5\f$ |
/// | 9 | \f$ -0.4 < \eta < 0.4\f$ |
/// \return kTRUE if proper and supported \f$ \eta \f$  cut code
///
Bool_t AliCSTrackCuts::SetEtaCut(Int_t etacode){
  // Set eta Cut
  switch(etacode){
    case 0:
      fCutsEnabledMask.ResetBitNumber(kEtaCut);
      fEtaCut = 100.;
      break;
    case 1:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 1.4;
      break;
    case 2:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 1.2;
      break;
    case 3:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 1.0;
      break;
    case 4:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.9;
      break;
    case 5:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.8;
      break;
    case 6:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.75;
      break;
    case 7:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.6;
      break;
    case 8:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.5;
      break;
    case 9:
      fCutsEnabledMask.SetBitNumber(kEtaCut);
      fEtaCut = 0.4;
      break;
    default:
      AliError(Form("Eta cut code %d not supported", etacode));
      return kFALSE;
  }
  return kTRUE;
}

/// Configures the ITS clusters cut
/// \param ITSclscode the ITS clusters cut code
/// | code | hits in ITS | hits in SPD |
/// |:--:|:--| :--|
/// | 0 | as standard for selected tracks in data period | same |
/// | 1 | n/a | at least one on first layer |
/// | 2 | n/a | at least one on any layer |
/// | 3 | at least 4 | at least one on first layer |
/// | 4 | at least 3 | at least one on any layer |
/// | 5 | at least 4 | at least one on any layer |
/// | 6 | at least 5 | at least one on any layer |
/// | 7 | n/a | none |
/// | 8 | SDD:FIRST | none |
/// | 9 | no ITS required | no ITS required |
/// \return kTRUE if proper and supported ITS clusters cut code
///
Bool_t AliCSTrackCuts::SetITSClustersCut(Int_t ITSclscode){

  switch(ITSclscode){
    case 0: // as standard for selected tracks in data period
      break;
    case 1: //1 hit first layer of SPD
      break;
    case 2: //1 hit in any layer of SPD
      break;
    case 3: // 4 hits in total in the ITS. At least 1 hit in the first layer of SPD
      break;
    case 4: // 3 hits in total in the ITS. At least 1 hit in any layer of SPD
      break;
    case 5: // 4 hits in total in the ITS. At least 1 hit in any layer of SPD
      break;
    case 6: // 5 hits in total in the ITS. At least 1 hit in any layer of SPD
      break;
    case 7: // no hits on the SPD
      break;
    case 8: // hit in first layer of SDD. No hits on the SPD
      break;
    case 9: // no ITS required
      break;
    default:
      AliError(Form("ITS cluster cut code %d not supported", ITSclscode));
      return kFALSE;
  }
  return kTRUE;
}

/// Configures the actual ITS clusters cut once configured for the data period
/// and according to the corresponding cut parameter
void AliCSTrackCuts::SetActualITSClustersCut(){

  if(!fESDTrackCuts ) {
    AliFatal("AliESDtrackCut is not initialized. ABORTING!!!");
  }

  switch(fParameters[kITSClsCutParam]){
    case 0: // as standard for selected tracks in data period
      /* do nothing */
      break;
    case 1: //1 hit first layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
      break;
    case 2: //1 hit in any layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      break;
    case 3: // 4 hits in total in the ITS. At least 1 hit in the first layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
      fESDTrackCuts->SetMinNClustersITS(4);
      break;
    case 4: // 3 hits in total in the ITS. At least 1 hit in any layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      fESDTrackCuts->SetMinNClustersITS(3);
      break;
    case 5: // 4 hits in total in the ITS. At least 1 hit in any layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      fESDTrackCuts->SetMinNClustersITS(4);
      break;
    case 6: // 5 hits in total in the ITS. At least 1 hit in any layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      fESDTrackCuts->SetMinNClustersITS(5);
      break;
    case 7: // no hits on the SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      break;
    case 8: // hit in first layer of SDD. No hits on any layer of SPD
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
      break;
    case 9:
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
      break;
    default:
      AliError(Form("ITS cluster cut code %d not supported", fParameters[kITSClsCutParam]));
  }
}

/// Prints the ITS clusters cut
void AliCSTrackCuts::PrintITSClustersCut() const {
  switch(fParameters[kITSClsCutParam]){
    case 0:
      printf("as standard for selected tracks in data period\n");
      break;
    case 1: //1 hit first layer of SPD
      printf("first SPD layer required\n");
      break;
    case 2: //1 hit in any layer of SPD
      printf("any of the SPD layers required\n");
      break;
    case 3: // 4 hits in total in the ITS. At least 1 hit in the first layer of SPD
      printf("first SPD layer required; at least four hits in the whole ITS\n");
      break;
    case 4: // 3 hits in total in the ITS. At least 1 hit in any layer of SPD
      printf("any of the SPD layers required; at least three hits in the whole ITS\n");
      break;
    case 5: // 4 hits in total in the ITS. At least 1 hit in any layer of SPD
      printf("any of the SPD layers required; at least four hits in the whole ITS\n");
      break;
    case 6: // 5 hits in total in the ITS. At least 1 hit in any layer of SPD
      printf("any of the SPD layers required; at least five hits in the whole ITS\n");
      break;
    case 7: // no hits on the SPD
      printf("no hit on the SPD\n");
      break;
    case 8: // hit in first layer of SDD. No hists on the SPD
      printf("no hit on the SPD; hit in the first layer of the SDD\n");
      break;
    case 9:
      printf("no SPD hits required\n");
      break;
    default:
      AliError(Form("Wrong ITS cluster cut code %d stored", fParameters[kITSClsCutParam]));
  }
}

/// Configures the TPC clusters cut
/// \param TPCclscode the TPC clusters cut code
/// | code | TPC clusters cut |
/// |:--:|:--|
/// | 0 | as standard for selected tracks in data period |
/// | 1 | minimum 70 clusters |
/// | 2 | minimum 80 clusters |
/// | 3 | minimum 100 clusters |
/// | 4 | minimum 50 crossed rows, minimum crossed rows / findable clusters ratio 60% |
/// | 5 | minimum 70 crossed rows, minimum crossed rows / findable clusters ratio 70% |
/// | 6 | minimum 70 crossed rows, minimum crossed rows / findable clusters ratio 80% |
/// | 7 | minimum 70 crossed rows, minimum crossed rows / findable clusters ratio 90% |
/// | | |
/// | 9 | none |
/// \return kTRUE if proper and supported TPC clusters cut code
///
Bool_t AliCSTrackCuts::SetTPCClustersCut(Int_t TPCclscode){

  switch(TPCclscode){
    case 0: // as standard for selected tracks in data period
      break;
    case 1:  // min 70 clusters
      break;
    case 2:  // min 80 clusters
      break;
    case 3:  // min 100
      break;
    case 4:  // min 50 crossed rows, min 60% crossed rows over findable clusters
      break;
    case 5:  // min 70 crossed rows, min 70% crossed rows over findable clusters
      break;
    case 6:  // min 70 crossed rows, min 80% crossed rows over findable clusters
      break;
    case 7:  // min 70 crossed rows, min 90% crossed rows over findable clusters
      break;
    case 9: // disable it
      break;
    default:
      AliError(Form("TPC cluster cut code %d not supported", TPCclscode));
      return kFALSE;
  }
  return kTRUE;
}

/// Configures the actual TPC clusters cut once configured for the data period
/// and according to the corresponding cut parameter
void AliCSTrackCuts::SetActualTPCClustersCut(){

  if(!fESDTrackCuts ) {
    AliFatal("AliESDtrackCut is not initialized. ABORTING!!!");
  }

  switch(fParameters[kTPCClsCutParam]){
    case 0: // as standard for selected tracks in data period
      /* do nothing */
      break;
    case 1:  // min 70 clusters
      fESDTrackCuts->SetMinNClustersTPC(70);
      break;
    case 2:  // min 80 clusters
      fESDTrackCuts->SetMinNClustersTPC(80);
      break;
    case 3:  // min 100
      fESDTrackCuts->SetMinNClustersTPC(100);
      break;
    case 4:  // min 50 crossed rows, min 60% crossed rows over findable clusters
      fESDTrackCuts->SetMinNCrossedRowsTPC(50);
      fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
      break;
    case 5:  // min 70 crossed rows, min 70% crossed rows over findable clusters
      fESDTrackCuts->SetMinNCrossedRowsTPC(70);
      fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
      break;
    case 6:  // min 70 crossed rows, min 80% crossed rows over findable clusters
      fESDTrackCuts->SetMinNCrossedRowsTPC(70);
      fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      break;
    case 7:  // min 70 crossed rows, min 90% crossed rows over findable clusters
      fESDTrackCuts->SetMinNCrossedRowsTPC(70);
      fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
      break;
    case 9: // disable it
      fESDTrackCuts->SetMinNClustersTPC();
      break;
    default:
      AliError(Form("TPC cluster cut code %d not supported", fParameters[kTPCClsCutParam]));
  }
}

/// Prints the TPC clusters cut
void AliCSTrackCuts::PrintTPCClustersCut() const {

  switch(fParameters[kTPCClsCutParam]){
    case 0:
      printf("TPC clusters cut as standard for selected tracks in data period\n");
      break;
    case 1:  // min 70 clusters
      printf("minimum 70 clusters required\n");
      break;
    case 2:  // min 80 clusters
      printf("minimum 80 clusters required\n");
      break;
    case 3:  // min 100
      printf("minimum 100 clusters required\n");
      break;
    case 4:  // min 50 crossed rows, min 60% crossed rows over findable clusters
      printf("minimum 50 crossed rows required, minimum 60%% crossed rows to findable clusters ratio\n");
      break;
    case 5:  // min 70 crossed rows, min 70% crossed rows over findable clusters
      printf("minimum 70 crossed rows required, minimum 70%% crossed rows to findable clusters ratio\n");
      break;
    case 6:  // min 70 crossed rows, min 80% crossed rows over findable clusters
      printf("minimum 70 crossed rows required, minimum 80%% crossed rows to findable clusters ratio\n");
      break;
    case 7:  // min 70 crossed rows, min 90% crossed rows over findable clusters
      printf("minimum 70 crossed rows required, minimum 90%% crossed rows to findable clusters ratio\n");
      break;
    case 9:
      printf("no TPC clusters cut required\n");
      break;
    default:
      AliError(Form("Wrong TPC cluster cut code %d stored", fParameters[kTPCClsCutParam]));
  }
}

/// Configures the track DCA cut
/// \param dcacode the DCA cut code
/// | code | max \f$ \mbox{DCA}_{XY} \f$ to vertex | max \f$ \mbox{DCA}_{Z} \f$ to vertex |
/// |:--:|:--:|:--:|
/// | 0 | as standard for selected tracks in data period | as standard for selected tracks in data period |
/// | 1 | tight DCA according to the data period | tight DCA according to the data period |
/// | 2 | loose DCA according to the data period | loose DCA according to the data period |
/// | 3 | tight DCA according to the data period | 1 |
/// | 4 | tight DCA according to the data period | 5 |
/// | 5 | 1 cm | 2 cm |
/// | 6 | \f$ 0.0525+\frac{0.175}{p_t^{1.1}} \f$ | 2 |
/// | 9 | not active cut | not active cut |
/// \return kTRUE if proper and supported DCA cut code
///
Bool_t AliCSTrackCuts::SetDCACut(Int_t dcacode)
{
  switch(dcacode){
    case 0:
      break;
    case 1:
      break;
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;
    case 9:
      break;
    default:
      AliError(Form("DCA cut code %d not supported", dcacode));
      return kFALSE;
  }
  return kTRUE;
}

/// Configures the actual DCA cut once configured for the data period
/// and according to the corresponding cut parameter
void AliCSTrackCuts::SetActualDCACut()
{
  if(!fESDTrackCuts ) {
    AliFatal("AliESDtrackCut is not initialized. ABORTING!!!");
  }
  switch(fParameters[kDCACutParam]){
    case 0: // as standard for selected tracks in data period
      /* do nothing */
      break;
    case 1: /* tight, tight */
      switch (fBaseSystem) {
        case k2010based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
          fESDTrackCuts->SetMaxDCAToVertexZ(2);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        case k2011based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
          fESDTrackCuts->SetMaxDCAToVertexZ(2);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        default:
        /* do nothing */
          break;
      }
      break;
    case 2: /* loose, loose */
      switch (fBaseSystem) {
        case k2010based:
        case k2011based:
          fESDTrackCuts->SetMaxDCAToVertexXY(2.4);
          fESDTrackCuts->SetMaxDCAToVertexZ(3.2);
          fESDTrackCuts->SetDCAToVertex2D(kTRUE);
          break;
        default:
        /* do nothing */
          break;
      }
      break;
    case 3: /* tight, 1 cm */
      switch (fBaseSystem) {
        case k2010based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
          fESDTrackCuts->SetMaxDCAToVertexZ(1);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        case k2011based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
          fESDTrackCuts->SetMaxDCAToVertexZ(1);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        default:
        /* do nothing */
          break;
      }
      break;
    case 4: /* tight, 5 cm */
      switch (fBaseSystem) {
        case k2010based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
          fESDTrackCuts->SetMaxDCAToVertexZ(5);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        case k2011based:
          fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
          fESDTrackCuts->SetMaxDCAToVertexZ(5);
          fESDTrackCuts->SetDCAToVertex2D(kFALSE);
          break;
        default:
        /* do nothing */
          break;
      }
      break;
    case 5: /* 1 cm , 2, cm */
      fESDTrackCuts->SetMaxDCAToVertexXY(1);
      fESDTrackCuts->SetMaxDCAToVertexZ(2);
      fESDTrackCuts->SetDCAToVertex2D(kTRUE);
      break;
    case 6: /* specific */
      fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0525+0.175/pt^1.1");
      fESDTrackCuts->SetMaxDCAToVertexZ(2);
      fESDTrackCuts->SetDCAToVertex2D(kFALSE);
      break;
    case 9:
      fESDTrackCuts->SetMaxDCAToVertexZ(1000);
      fESDTrackCuts->SetMaxDCAToVertexXY(1000);
      fESDTrackCuts->SetDCAToVertex2D(kFALSE);
      break;
    default:
      AliError(Form("DCA cut code %d not supported", fParameters[kDCACutParam]));
  }
}

/// Prints the DCA cut information
void AliCSTrackCuts::PrintDCACut() const
{
  switch(fParameters[kDCACutParam]){
    case 0:
      printf("as standard for selected tracks in data period\n");
      break;
    case 1:
      printf("Z: tight according to the data period, XY: tight according to the data period\n");
      break;
    case 2:
      printf("Z: loose according to the data period, XY: loose according to the data period\n");
      break;
    case 3:
      printf("Z: 1 cm, XY: tight according to the data period\n");
      break;
    case 4:
      printf("Z: 5 cm, XY: loose according to the data period\n");
      break;
    case 5:
      printf("Z: 2 cm, XY: 1 cm\n");
      break;
    case 6:
      printf("Z: 2 cm, XY: 0.0525+0.175/pt^1.1\n");
      break;
    case 9:
      printf("none\n");
      break;
    default:
      AliError(Form("Wrong DCA cut code %d stored", fParameters[kDCACutParam]));
  }
}

/// Configures the \f$ p_{T} \f$ cut
/// \param ptcode the \f$ p_{T} \f$ cut code
/// | code | minimum \f$ p_{T} \f$ (GeV/c) | maximum \f$ p_{T} \f$ (GeV/c) |
/// |:--:|:--:|:--:|
/// | 0 | not active cut | not active cut |
/// | 1 | 0.2 | 2.0 |
/// | 2 | 0.2 | 3.0 |
/// | 3 | 0.2 | 5.0 |
/// | 4 | 0.2 | 1.4 |
/// | 5 | 0.2 | 1.6 |
/// | 6 | 0.2 | 1.8 |
/// \return kTRUE if proper and supported \f$ p_{T} \f$ cut code
///
Bool_t AliCSTrackCuts::SetPtCut(Int_t ptcode)
{
  switch(ptcode){
  case 0:
    fCutsEnabledMask.ResetBitNumber(kPtCut);
    fMinPt = 0.0;
    fMaxPt = 9999.0;
    break;
  case 1:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 2.0;
    break;
  case 2:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 3.0;
    break;
  case 3:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 5.0;
    break;
  case 4:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 1.4;
    break;
  case 5:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 1.6;
    break;
  case 6:
    fCutsEnabledMask.SetBitNumber(kPtCut);
    fMinPt = 0.2;
    fMaxPt = 1.8;
    break;
  default:
    AliError(Form("Pt cut code %d not supported", ptcode));
    return kFALSE;
  }
  return kTRUE;
}

/// Initializes the cuts
///
/// Initializes the needed data and allocates the needed histograms list if needed
/// \param name an additional name to precede the cuts string
void AliCSTrackCuts::InitCuts(const char *name){

  if (name == NULL) name = GetName();

  if (fQALevel > kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    if(fHistogramsList != NULL)
      delete fHistogramsList;

    fHistogramsList =new TList();
    fHistogramsList->SetOwner(kTRUE);
    fHistogramsList->SetName(name);

    TH1::AddDirectory(oldstatus);
  }
}

/// Allocates the different histograms if needed
///
/// It is supposed that the current cuts string is the running one
void AliCSTrackCuts::DefineHistograms(){

  if (fQALevel > kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    /* the original name is used as title for the statistics histogram so, preserve it */
    TString originalTempName = fHistogramsList->GetName();
    fHistogramsList->SetName(Form("%s_%s",fHistogramsList->GetName(),GetCutsString()));

    fhCutsStatistics = new TH1F(Form("CutsStatistics_%s",GetCutsString()),Form("%s tracks cuts statistics",originalTempName.Data()),kNCuts+4,-0.5,kNCuts+3.5);
    fhCutsStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
    fhCutsStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");
    for (Int_t i = 0; i < kNCuts; i++)
      fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
    fHistogramsList->Add(fhCutsStatistics);

    if(fQALevel == kQALevelHeavy){
      fhCutsCorrelation = new TH2F(Form("CutCorrelation_%s",GetCutsString()),"Cuts correlation",kNCuts+2,-0.5,kNCuts+1.5,kNCuts+2,-0.5,kNCuts+1.5);
      for (Int_t i=0; i<kNCuts; i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
        fhCutsCorrelation->GetYaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
      }
      fHistogramsList->Add(fhCutsCorrelation);
    }

    fhPtVsDCAxy[0]  = new TH2F(Form("PtVsDCAxyPtB_%s",GetCutsString()),"p_{T} vs DCA_{XY} before;DCA_{XY} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fhPtVsDCAxy[1]  = new TH2F(Form("PtVsDCAxyPtA_%s",GetCutsString()),"p_{T} vs DCA_{XY};DCA_{XY} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsDCAxy[0]);
    fHistogramsList->Add(fhPtVsDCAxy[1]);

    fhPtVsDCAz[0]  = new TH2F(Form("PtVsDCAzPtB_%s",GetCutsString()),"p_{T} vs DCA_{Z} before;DCA_{Z} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fhPtVsDCAz[1]  = new TH2F(Form("PtVsDCAzPtA_%s",GetCutsString()),"p_{T} vs DCA_{Z};DCA_{Z} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsDCAz[0]);
    fHistogramsList->Add(fhPtVsDCAz[1]);

    fhPtVsTPCCls[0]  = new TH2F(Form("PtVsTPCClustersB_%s",GetCutsString()),"p_{T} vs no. of TPC clusters before;no of clusters;p_{T} (GeV/c)",170,0,170,400,0.,10.);
    fhPtVsTPCCls[1]  = new TH2F(Form("PtVsTPCClustersA_%s",GetCutsString()),"p_{T} vs no. of TPC clusters;no of clusters;p_{T} (GeV/c)",170,0,170,400,0.,10.);
    fHistogramsList->Add(fhPtVsTPCCls[0]);
    fHistogramsList->Add(fhPtVsTPCCls[1]);

    fhPtVsTPCRowOverFindCls[0] = new TH2F(Form("PtVsTPCCROFCB_%s",GetCutsString()),
        "p_{T} vs TPC crossed rows findable clusters ratio before;crossed rows / findable clusters;p_{T} (GeV/c)",100,0,1,400,0.,10.);
    fhPtVsTPCRowOverFindCls[1] = new TH2F(Form("PtVsTPCCROFCA_%s",GetCutsString()),
        "p_{T} vs TPC crossed rows findable clusters ratio;crossed rows / findable clusters;p_{T} (GeV/c)",100,0,1,400,0.,10.);
    fHistogramsList->Add(fhPtVsTPCRowOverFindCls[0]);
    fHistogramsList->Add(fhPtVsTPCRowOverFindCls[1]);
    fhPtVsEta[0] = new TH2F(Form("PtVsEtaB_%s",GetCutsString()),"p_{T} vs #eta before;#eta;p_{T} (GeV/c)",100,-2.0,2.0,400,0.,10.);
    fhPtVsEta[1] = new TH2F(Form("PtVsEtaA_%s",GetCutsString()),"p_{T} vs #eta;#eta;p_{T} (GeV/c)",100,-2.0,2.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsEta[0]);
    fHistogramsList->Add(fhPtVsEta[1]);

    if(fQALevel == kQALevelHeavy){
      fhEtaVsPhi[0] = new TH2F(Form("EtaVsPhiB_%s",GetCutsString()),"#eta vs #phi before;#phi;#eta", 360, 0.0, 360.0, 100, -2.0, 2.0);
      fhEtaVsPhi[1] = new TH2F(Form("EtaVsPhiA_%s",GetCutsString()),"#eta vs #phi;#phi;#eta", 360, 0.0, 360.0, 100, -2.0, 2.0);
      fHistogramsList->Add(fhEtaVsPhi[0]);
      fHistogramsList->Add(fhEtaVsPhi[1]);

      fESDTrackCuts->SetHistogramsOn(kTRUE);
      fESDTrackCuts->DefineHistograms();
      fHistogramsList->Add(fESDTrackCuts);
    }

    TH1::AddDirectory(oldstatus);
  }
}

/// \cond CLASSIMP
ClassImp(AliCSTrackCuts);
/// \endcond
