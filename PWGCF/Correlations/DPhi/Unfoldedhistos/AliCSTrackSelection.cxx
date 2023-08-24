/**************************************************************************
* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TObjString.h>
#include <TBits.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include "AliVTrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliCSTrackMaps.h"
#include "AliCSTrackSelection.h"
#include "AliCSTrackCuts.h"
#include "AliPIDResponse.h"
#include "AliCSPIDCuts.h"
#include "AliLog.h"
#include "TH3.h"

/// \file AliCSTrackSelection.cxx
/// \brief Implementation of track selection class within the correlation studies analysis

/// Default constructor for serialization
AliCSTrackSelection::AliCSTrackSelection()
  : TNamed(),
    fPoiNames{},
    fTrackId(kWrongPOIid),
    fParticleId(kWrongPOIid),
    fQALevel(AliCSAnalysisCutsBase::kQALevelNone),
    fInclusiveTrackCuts(),
    fInclusivePIDCuts(),
    fInclusiveCutsStrings(),
    fInclusivePidCutsStrings(),
    fExclusiveCuts(),
    fExclusiveCutsStrings(),
    fExclusivePidCutsStrings(),
    fCutsActivatedMask(NULL),
    fCutsString(""),
    fDataPeriod(AliCSAnalysisCutsBase::kNoPeriod),
    fHighestPrimaryLabel(0),
    fPIDResponse(NULL),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhPtVsDCAxy{nullptr},
    fhPtVsDCAz{nullptr},
    fhPtVsTPCCls{nullptr},
    fhPtVsTPCRows{nullptr},
    fhPtVsTPCRowOverFindCls{nullptr},
    fhEtaVsPhi{nullptr},
    fhPtVsEta{nullptr},
    fhITSdEdxSignalVsP{nullptr},
    fhTPCdEdxSignalVsP{nullptr},
    fhTOFSignalVsP{nullptr},
    fhTOFSelSignalVsP{nullptr},
    fhTPCTOFSigmaVsP{{nullptr}},
    fhTPCdEdxSignalDiffVsP{{nullptr}},
    fhTOFSignalDiffVsP{{nullptr}},
    fhTPCdEdxSelSignalDiffVsP{{nullptr}},
    fhTOFSelSignalDiffVsP{{nullptr}},
    fhPvsTOFMassSq{nullptr},
    fhSelPvsTOFMassSq{nullptr},
    fHistogramsList(NULL)
{
  /* we own the cuts so, they will be destroyed when we are */
  fInclusiveTrackCuts.SetOwner(kTRUE);
  fInclusivePIDCuts.SetOwner(kTRUE);
  fInclusiveCutsStrings.SetOwner(kTRUE);
  fInclusivePidCutsStrings.SetOwner(kTRUE);
  fExclusiveCuts.SetOwner(kTRUE);
  fExclusiveCutsStrings.SetOwner(kTRUE);
  fExclusivePidCutsStrings.SetOwner(kTRUE);
}

/// Constructor
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSTrackSelection::AliCSTrackSelection(const char* name, const char* title)
  : TNamed(name, title),
    fTrackId(kWrongPOIid),
    fParticleId(kWrongPOIid),
    fQALevel(AliCSAnalysisCutsBase::kQALevelNone),
    fInclusiveTrackCuts(),
    fInclusivePIDCuts(),
    fInclusiveCutsStrings(),
    fInclusivePidCutsStrings(),
    fExclusiveCuts(),
    fExclusiveCutsStrings(),
    fExclusivePidCutsStrings(),
    fCutsActivatedMask(NULL),
    fCutsString(""),
    fDataPeriod(AliCSAnalysisCutsBase::kNoPeriod),
    fHighestPrimaryLabel(0),
    fPIDResponse(NULL),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhPtVsDCAxy{nullptr},
    fhPtVsDCAz{nullptr},
    fhPtVsTPCCls{nullptr},
    fhPtVsTPCRows{nullptr},
    fhPtVsTPCRowOverFindCls{nullptr},
    fhEtaVsPhi{nullptr},
    fhPtVsEta{nullptr},
    fhITSdEdxSignalVsP{nullptr},
    fhTPCdEdxSignalVsP{nullptr},
    fhTOFSignalVsP{nullptr},
    fhTOFSelSignalVsP{nullptr},
    fhTPCTOFSigmaVsP{{nullptr}},
    fhTPCdEdxSignalDiffVsP{{nullptr}},
    fhTOFSignalDiffVsP{{nullptr}},
    fhTPCdEdxSelSignalDiffVsP{{nullptr}},
    fhTOFSelSignalDiffVsP{{nullptr}},
    fhPvsTOFMassSq{nullptr},
    fhSelPvsTOFMassSq{nullptr},
    fHistogramsList(NULL)
{
  /* we own the cuts so, they will be destroyed when we are */
  fInclusiveTrackCuts.SetOwner(kTRUE);
  fInclusivePIDCuts.SetOwner(kTRUE);
  fInclusiveCutsStrings.SetOwner(kTRUE);
  fInclusivePidCutsStrings.SetOwner(kTRUE);
  fExclusiveCuts.SetOwner(kTRUE);
  fExclusiveCutsStrings.SetOwner(kTRUE);
  fExclusivePidCutsStrings.SetOwner(kTRUE);
}

/// Destructor
/// We don't own anything, everything we allocate is owned
/// by the output list
AliCSTrackSelection::~AliCSTrackSelection()
{

}

/// \brief Sets the track selection configuration by means of a configuration string
/// \param confstring the configuration string
/// \return kTRUE if the instance was correctly configured kFALSE otherwise
/// The configuration string will have the format
/// Tracks:trackstring;PID:pidstring;
/// In their turn, trackstring and pidstring might have the format
/// cutstring,cutstring+cutstring and cutstring could start with a minus
/// to configure rejection cuts instead of acceptance cuts. For PID
/// strings, the first character/s (second after a potential '-') indicate
/// the family towards the PID cut is applied: 'e' for electrons, 'pi' for
/// pions, 'mu' for muons, 'k' for kaons and 'p' for protons.
///
Bool_t AliCSTrackSelection::InitializeFromString(const char *confstring)
{
  Bool_t accepted = kTRUE;
  fCutsString = confstring;
  TString sztmp = confstring;

  /* the track cuts strings                           */
  if (sztmp.BeginsWith("Tracks:")) {
    sztmp.Remove(0,strlen("Tracks:"));
    TString sztmptrk = sztmp(0,sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);

    TObjArray *tokens = sztmptrk.Tokenize("+");
    for (Int_t icut = 0; icut < tokens->GetEntries(); icut++) {
      if (((TObjString*) tokens->At(icut))->String().BeginsWith("-")) {
        TObjString *posz = new TObjString(((TObjString*) tokens->At(icut))->String().Remove(0,1).Data());
        fExclusiveCutsStrings.Add(posz);
      }
      else {
        TObjString *posz = new TObjString(((TObjString*) tokens->At(icut))->String().Data());
        fInclusiveCutsStrings.Add(posz);
      }
    }
    delete tokens;
  }
  else {
    AliFatal("No tracks selection cuts string. ABORTING!!!");
    return kFALSE;
  }

  /* the PID cuts strings                                                 */
  if (sztmp.BeginsWith("PID:")) {
    sztmp.Remove(0,strlen("PID:"));
    TString sztmpid = sztmp(0, sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);

    /* let's build the inclusive and exclusive PID cuts lists */
    auto populate = [](TObjArray& lst, TString str) {
      TObjArray* tok = str.Tokenize("&");
      for (int icf = 0; icf < tok->GetEntries(); ++icf) {
        TObjString* cfstr = new TObjString(((TObjString*)tok->At(icf))->String());
        if (cfstr->String().BeginsWith("e") || cfstr->String().BeginsWith("mu") || cfstr->String().BeginsWith("pi") || cfstr->String().BeginsWith("p") || cfstr->String().BeginsWith("k")) {
          lst.Add(cfstr);
        } else {
          ::Fatal("AliCSTrackSelection::InitializeFromString()", "Wrong family in track selection PID cuts string. ABORTING!!!");
          return false;
        }
      }
      delete tok;
      return true;
    };
    TObjArray *tokens = sztmpid.Tokenize("+");
    for (Int_t icut = 0; icut < tokens->GetEntries(); icut++) {
      if (((TObjString*) tokens->At(icut))->String().BeginsWith("-")) {
        accepted = accepted && populate(fExclusivePidCutsStrings, ((TObjString*)tokens->At(icut))->String().Remove(0, 1));
      }
      else {
        accepted = accepted && populate(fInclusivePidCutsStrings, ((TObjString*)tokens->At(icut))->String());
      }
    }
    delete tokens;
  }
  else {
    AliFatal("No tracks selection PID cuts string. ABORTING!!!");
    return kFALSE;
  }

  /* let's get the PID response instance if needed */
  if (fInclusivePidCutsStrings.GetEntries() != 0 || fExclusivePidCutsStrings.GetEntries() != 0) {
    AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
    if(manager != NULL) {
      if (!AliCSAnalysisCutsBase::IsOnTheFlyMC()) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*) (manager->GetInputEventHandler());
        fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
        /* if we need PID response instance and it is not there we cannot continue */
        if (fPIDResponse == NULL)
          AliFatal("No PID response instance. ABORTING!!!");
      }
    }
    else {
      AliFatal("No analysis manager instance. ABORTING!!!");
    }
  }

  AliInfo("Stored track selection configuration string:");
  AliInfo(Form("\t%s", fCutsString.Data()));
  return accepted;
}


/// Processes a potential change in the run number
///
/// Checks if the current period under analysis has changes and if so
/// updates the needed members
void AliCSTrackSelection::NotifyRun() {

  /* checks the change in the analysis period */
  if (AliCSAnalysisCutsBase::GetGlobalPeriod() != fDataPeriod) {

    fDataPeriod = AliCSAnalysisCutsBase::GetGlobalPeriod();

    fCutsActivatedMask = new TBits(fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()+fExclusiveCuts.GetEntriesFast());

    /* Inform the different set of cuts */
    for (Int_t ix = 0; ix < fInclusiveTrackCuts.GetEntriesFast(); ix++)
      ((AliCSTrackCutsBase *) fInclusiveTrackCuts[ix])->NotifyRun();
    for (Int_t ix = 0; ix < fInclusivePIDCuts.GetEntriesFast(); ix++)
      ((AliCSTrackCutsBase *) fInclusivePIDCuts[ix])->NotifyRun();
    for (Int_t ix = 0; ix < fExclusiveCuts.GetEntriesFast(); ix++)
      ((AliCSTrackCutsBase *) fExclusiveCuts[ix])->NotifyRun();

    /* and now we ask for histogram allocation */
    DefineHistograms();
  }
}

/// Processes the start of a new event
///
/// Checks if the running analysis is over MC data and stored the
/// needed data for track selection
void AliCSTrackSelection::NotifyEvent() {

  if (AliCSAnalysisCutsBase::IsMC()) {
    AliMCEventHandler *eventHandler = AliCSAnalysisCutsBase::GetMCEventHandler();

    AliGenEventHeader *mainHeader = NULL;
    Int_t nNoOfHeaders = 0;
    fHighestPrimaryLabel = 0;

    if (eventHandler != NULL) {
      /* MC ESD format */
      AliMCEvent *mcEvent = eventHandler->MCEvent();
      AliHeader *header = mcEvent->Header();
      if (header != NULL) {
        AliGenCocktailEventHeader *cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
        if (cocktailHeader != NULL) {
          nNoOfHeaders = cocktailHeader->GetHeaders()->GetEntries();
          mainHeader = dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());
        }
        else {
          if (dynamic_cast<AliGenEventHeader*>(header->GenEventHeader()) != NULL) {
            AliInfo("From header->GenEventHeader");
            mainHeader = dynamic_cast<AliGenEventHeader*>(header->GenEventHeader());
            nNoOfHeaders = 1;
          }
          else {
            AliError("MC analysis but no MC event cocktail header");
          }
        }
      }
      else {
        AliError("MC analysis but no MC event header");
      }
    }
    else {
      /* MC AOD format */
      AliInputEventHandler *handler = AliCSAnalysisCutsBase::GetInputEventHandler();
      AliAODEvent *event = (AliAODEvent *)handler->GetEvent();
      AliAODMCHeader* header = (AliAODMCHeader*) event->GetList()->FindObject(AliAODMCHeader::StdBranchName());

      if (header != NULL) {
        nNoOfHeaders = header->GetNCocktailHeaders();
        mainHeader = header->GetCocktailHeader(0);
      }
      else {
        AliError("MC analysis but no MC event cocktail header");
      }
    }

    if (mainHeader != NULL) {
      if (nNoOfHeaders > 1) {
        AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks with labels above %d.",
            nNoOfHeaders, mainHeader->ClassName(), mainHeader->NProduced()-1));
        fHighestPrimaryLabel = mainHeader->NProduced();
      }
    }
    else {
      AliError("MC anlysis but no proper header");
    }
  }
}

/// get the paritcle of interest internal id for the PID species
/// \param sp the PID species
/// \return the internal POI id
AliCSTrackSelection::poiIds AliCSTrackSelection::poiid(AliPID::EParticleType sp)
{
  switch (sp) {
    case AliPID::kPion:
      return kPOIpi;
    case AliPID::kKaon:
      return kPOIka;
    case AliPID::kProton:
      return kPOIpr;
    default:
      return kWrongPOIid;
  }
};

/// Check whether the passed track is accepted by the different cuts
/// \param trk the track to analyze whether it is accepted or not
/// \return kTRUE if the track  is accepted, kFALSE otherwise
///
/// The track shall be accepted by any of the inclusive set of track cuts
/// and rejected by all the exclusive ones
Bool_t AliCSTrackSelection::IsTrackAccepted(AliVTrack *trk) {

  /* just to be sure */
  if (trk == NULL) return kFALSE;

  /* if MC analysis is ongoing get rid of "ghost" tracks */
  if(AliCSAnalysisCutsBase::IsMC()) if (trk->GetLabel() < 0) return kFALSE;

  /* if MC analysis is ongoing get rid of injected signals */
  if (AliCSAnalysisCutsBase::IsMC()) if (IsFromMCInjectedSignal(trk->GetLabel())) return kFALSE;

  /* for the time being */
  Bool_t accepted = kTRUE;
  Bool_t inclusivetrack = (fInclusiveTrackCuts.GetEntries() == 0);
  Bool_t inclusivepid = true;
  Bool_t exclusive = kFALSE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask->ResetAllBits();

  /* let's get the DCA at this level because it will be latter used */
  bool gooddca = true;
  float dca[2]{0.0f, 0.0f}, bCov[3]{0.0f, 0.0f, 0.0f};
  if (trk->IsA() == AliESDtrack::Class()) {
    trk->GetImpactParameters(dca, bCov);
  } else {
    AliAODTrack* aodt = dynamic_cast<AliAODTrack*>(trk);
    float pos[3] = {0.0f, 0.0f, 0.0f};
    if (aodt->GetPosition(pos)) {
      /* we got a DCA */
      dca[0] = pos[0];
      dca[1] = pos[1];
    } else {
      /* the track was not constrained to the PV */
      double dcap[2]{0.0f, 0.0f}, bCovp[3]{0.0f, 0.0f, 0.0f};
      AliExternalTrackParam etp;
      etp.CopyFromVTrack(aodt);
      if (!etp.PropagateToDCA(AliCSAnalysisCutsBase::GetVertex(), AliCSAnalysisCutsBase::GetMagField(), 3.0, dcap, bCovp)) {
        /* wrong propagation to DCA */
        gooddca = false;
      } else if (bCovp[0] <= 0 || bCovp[2] <= 0) {
        /* wrong DCA resolution */
        gooddca = false;
      } else {
        dca[0] = dcap[0];
        dca[1] = dcap[1];
      }
    }
  }

  /* Check first the inclusive set of track cuts */
  for (Int_t ix = 0; ix < fInclusiveTrackCuts.GetEntriesFast(); ix++) {

    Bool_t cutaccepted = ((AliCSTrackCutsBase*)fInclusiveTrackCuts[ix])->IsTrackAccepted(trk, dca);
    /* if track is not accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix,!cutaccepted);
    inclusivetrack = inclusivetrack || cutaccepted;
  }

  /* now the inclusive set of pid cuts checking their consistency */
  fTrackId = kWrongPOIid;
  if (fInclusivePIDCuts.GetEntriesFast() > 0) {
    for (Int_t ix = 0; ix < fInclusivePIDCuts.GetEntriesFast(); ix++) {
      Bool_t cutaccepted = ((AliCSTrackCutsBase*)fInclusivePIDCuts[ix])->IsTrackAccepted(trk, dca);
      /* if track is not accepted the cut is activated */
      fCutsActivatedMask->SetBitNumber(ix + fInclusiveTrackCuts.GetEntriesFast(), !cutaccepted);
      /* check the consistency */
      if (inclusivepid && cutaccepted) {
        /* still consistent */
        if (fTrackId == kWrongPOIid) {
          fTrackId = poiid(((AliCSPIDCuts*)fInclusivePIDCuts[ix])->GetTargetSpecies());
        } else {
          if (fTrackId != poiid(((AliCSPIDCuts*)fInclusivePIDCuts[ix])->GetTargetSpecies())) {
            /* more than one id associated to the track, track is rejected but we keep on the loop for getting all statistics */
            fTrackId = kWrongPOIid;
            inclusivepid = false;
          }
        }
      }
    }
    /* just in case none were accepted */
    inclusivepid = fTrackId != kWrongPOIid;
  } else {
    /* if the track is reconstructed and accepted then by default is a hadron */
    fTrackId = kPOIh;
  }

  /* and now the exclusive ones */
  for (Int_t ix = 0; ix < fExclusiveCuts.GetEntriesFast(); ix++) {
    Bool_t cutaccepted = ((AliCSTrackCutsBase*)fExclusiveCuts[ix])->IsTrackAccepted(trk, dca);
    /* if track is accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast(),cutaccepted);
    exclusive = exclusive || cutaccepted;
  }

  /* now decide if accepted or not */
  accepted = gooddca && inclusivetrack && inclusivepid && !exclusive;

  if (fQALevel > AliCSAnalysisCutsBase::kQALevelNone) {
    /* let's fill the histograms */
    fhCutsStatistics->Fill("n tracks", 1);
    if (!accepted)
      fhCutsStatistics->Fill("n cut tracks", 1);

    for (Int_t ix=0; ix<fInclusiveTrackCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindFixBin(fInclusiveTrackCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fInclusiveTrackCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix))
        fhCutsStatistics->Fill(fInclusiveTrackCuts[ix]->GetName(), 1);

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx = ix; jx < fInclusiveTrackCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fInclusiveTrackCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusiveTrackCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix) && fCutsActivatedMask->TestBitNumber(jx)) {
            fhCutsCorrelation->Fill(fInclusiveTrackCuts[ix]->GetName(), fInclusiveTrackCuts[jx]->GetName(), 1);
          }
        }
        for (Int_t jx = 0; jx < fInclusivePIDCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fInclusivePIDCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusivePIDCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix) && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast())) {
            fhCutsCorrelation->Fill(fInclusiveTrackCuts[ix]->GetName(), fInclusivePIDCuts[jx]->GetName(), 1);
          }
        }
        for (Int_t jx=0; jx<fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix)
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            fhCutsCorrelation->Fill(fInclusiveTrackCuts[ix]->GetName(), fExclusiveCuts[jx]->GetName(), 1);
          }
        }
      }
    }
    for (Int_t ix=0; ix<fInclusivePIDCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindFixBin(fInclusivePIDCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fInclusivePIDCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()))
        fhCutsStatistics->Fill(fInclusivePIDCuts[ix]->GetName(), 1);

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx = ix; jx < fInclusivePIDCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fInclusivePIDCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusivePIDCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast())) {
            fhCutsCorrelation->Fill(fInclusivePIDCuts[ix]->GetName(), fInclusivePIDCuts[jx]->GetName(), 1);
          }
        }
        for (Int_t jx=0; jx<fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            fhCutsCorrelation->Fill(fInclusivePIDCuts[ix]->GetName(), fExclusiveCuts[jx]->GetName(), 1);
          }
        }
      }
    }
    for (Int_t ix=0; ix<fExclusiveCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindFixBin(fExclusiveCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fExclusiveCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()))
        fhCutsStatistics->Fill(fExclusiveCuts[ix]->GetName(), 1);

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx = ix; jx < fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindFixBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            fhCutsCorrelation->Fill(fExclusiveCuts[ix]->GetName(), fExclusiveCuts[jx]->GetName(), 1);
          }
        }
      }
    }

    /* we now need to consider the potential constrained track */
    AliVTrack *ttrk = trk;
    if (trk->GetID() < 0)
      ttrk = AliCSTrackMaps::GetOriginalTrack(dynamic_cast<AliAODTrack*>(trk));
    Float_t nCrossedRowsTPC = ttrk->GetTPCCrossedRows();
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (ttrk->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / ttrk->GetTPCNclsF();
    }

    for (Int_t i = 0; i < 2; i++) {

      fhPtVsDCAxy[i]->Fill(dca[0],trk->Pt());
      fhPtVsDCAz[i]->Fill(dca[1],trk->Pt());
      fhPtVsTPCCls[i]->Fill(ttrk->GetTPCNcls(),trk->Pt());
      fhPtVsTPCRows[i]->Fill(ttrk->GetTPCCrossedRows(), trk->Pt());
      fhPtVsTPCRowOverFindCls[i]->Fill(ratioCrossedRowsOverFindableClustersTPC,trk->Pt());
      fhPtVsEta[i]->Fill(trk->Eta(),trk->Pt());

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        fhEtaVsPhi[i]->Fill(trk->Phi()*180.0/TMath::Pi(),trk->Eta());
      }

      if (fInclusivePidCutsStrings.GetEntries() != 0 || fExclusivePidCutsStrings.GetEntries() != 0) {
        /* don't fill before if not requested */
        if (i == 1 || (fQALevel > AliCSAnalysisCutsBase::kQALevelLight)) {
          auto poiidx = [&](poiIds id) {
            switch (id) {
              case kPOIpi:
                return 0;
              case kPOIka:
                return 1;
              case kPOIpr:
                return 2;
              default:
                AliFatal(Form("Wrong POI id %d for filling selected tracks histograms", id));
                return -1;
            }
          };
          fhITSdEdxSignalVsP[i]->Fill(trk->P(),ttrk->GetITSsignal());
          fhTPCdEdxSignalVsP[i]->Fill(trk->P(),TMath::Abs(ttrk->GetTPCsignal()));
          if (i == 1 && fInclusivePidCutsStrings.GetEntries() != 0) {
            /* PID selection */
            fhTPCdEdxSelSignalVsP[poiidx(fTrackId)]->Fill(trk->P(), TMath::Abs(ttrk->GetTPCsignal()));
          }
          auto gettofbeta = [&](auto tofsignal) {
            static const Double_t c_cm_ps = TMath::C() * 1.0e2 * 1.0e-12;
            double tracklen_cm = trk->GetIntegratedLength();
            double toftime_ps = tofsignal - fPIDResponse->GetTOFResponse().GetStartTime(trk->P());
            return tracklen_cm / toftime_ps / c_cm_ps;
          };
          auto gettofmasssq = [](auto const& trk, auto beta) {
            return trk->P() * trk->P() * (1.0 / beta / beta - 1.0);
          };
          if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
            double beta = gettofbeta(ttrk->GetTOFsignal());
            double tofmasssq = gettofmasssq(trk, beta);
            fhTOFSignalVsP[i]->Fill(trk->P(), beta);
            fhPvsTOFMassSq[i]->Fill(tofmasssq, trk->P());
            if (i == 1 && fInclusivePidCutsStrings.GetEntries() != 0) {
              /* PID selection */
              fhTOFSelSignalVsP[poiidx(fTrackId)]->Fill(trk->P(), beta);
              fhSelPvsTOFMassSq[poiidx(fTrackId)]->Fill(tofmasssq, trk->P());
            }
          }
          auto spec = [](int ix) {
            switch (ix) {
              case 0:
                return AliPID::kPion;
                break;
              case 1:
                return AliPID::kKaon;
                break;
              case 2:
                return AliPID::kProton;
                break;
            }
            return AliPID::kUnknown;
          };
          for (int j = 0; j < 3; ++j) {
            fhTPCdEdxSignalDiffVsP[j][i]->Fill(ttrk->P(), ttrk->GetTPCsignal() - fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, ttrk, spec(j)));
            if (i == 1 && fInclusivePidCutsStrings.GetEntries() != 0) {
              /* PID selection */
              fhTPCdEdxSelSignalDiffVsP[j][poiidx(fTrackId)]->Fill(ttrk->P(), ttrk->GetTPCsignal() - fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, ttrk, spec(j)));
              if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
                fhTOFSelSignalDiffVsP[j][poiidx(fTrackId)]->Fill(ttrk->P(), gettofbeta(ttrk->GetTOFsignal()) - gettofbeta(fPIDResponse->GetTOFResponse().GetExpectedSignal(ttrk, spec(j))));
              }
            }
            if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
              fhTPCTOFSigmaVsP[j][i]->Fill(fPIDResponse->NumberOfSigmasTPC(ttrk, spec(j)), fPIDResponse->NumberOfSigmasTOF(ttrk, spec(j)), ttrk->P());
              fhTOFSignalDiffVsP[j][i]->Fill(ttrk->P(), gettofbeta(ttrk->GetTOFsignal()) - gettofbeta(fPIDResponse->GetTOFResponse().GetExpectedSignal(ttrk, spec(j))));
            }
          }
        }
      }

      /* don't fill after if event not accepted */
      if (!accepted) break;
    }
  }
  return accepted;
}


/// Check whether the true track associated to the passed track is accepted by the different cuts
/// \param trk the track to analyze whether its associated true track is accepted or not
/// \return kTRUE if the associated true track  is accepted, kFALSE otherwise
///
Bool_t AliCSTrackSelection::IsTrueTrackAccepted(AliVTrack *trk) {

  /* reject ghost tracks */
  if (trk->GetLabel() < 0) return kFALSE;

  return IsTrueTrackAccepted(trk->GetLabel());
}


/// Check whether the passed true track is accepted by the inclusive cuts and rejected by the exclusive cuts
/// \param itrk the index of true track to analyze whether it is accepted or not
/// \return kTRUE if the track  is accepted, kFALSE otherwise
///
Bool_t AliCSTrackSelection::IsTrueTrackAccepted(Int_t itrk) {

  /* discard injected signals */
  if (IsFromMCInjectedSignal(itrk)) return kFALSE;

  /* for the time being */
  Bool_t accepted = kTRUE;
  Bool_t inclusivetrack = (fInclusiveTrackCuts.GetEntriesFast() == 0);
  Bool_t inclusivepid = true;
  Bool_t exclusive = kFALSE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask->ResetAllBits();

  /* Check first the inclusive set of track cuts */
  for (Int_t ix = 0; ix < fInclusiveTrackCuts.GetEntriesFast(); ix++) {

    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fInclusiveTrackCuts[ix])->IsTrueTrackAccepted(itrk);
    /* if track is not accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix,!cutaccepted);
    inclusivetrack = inclusivetrack || cutaccepted;
  }

  /* now the inclusive set of pid cuts checking their consistency */
  fParticleId = kWrongPOIid;
  if (fInclusivePIDCuts.GetEntriesFast() > 0) {
    for (Int_t ix = 0; ix < fInclusivePIDCuts.GetEntriesFast(); ix++) {

      Bool_t cutaccepted = ((AliCSTrackCutsBase*)fInclusivePIDCuts[ix])->IsTrueTrackAccepted(itrk);
      /* if track is not accepted the cut is activated */
      fCutsActivatedMask->SetBitNumber(ix + fInclusiveTrackCuts.GetEntriesFast(), !cutaccepted);
      /* check the consistency */
      if (inclusivepid && cutaccepted) {
        /* still consistent */
        if (fParticleId == kWrongPOIid) {
          fParticleId = poiid(((AliCSPIDCuts*)fInclusivePIDCuts[ix])->GetTargetSpecies());
        } else {
          if (fParticleId != poiid(((AliCSPIDCuts*)fInclusivePIDCuts[ix])->GetTargetSpecies())) {
            /* more than one id associated to the track, track is rejected but we keep on the loop for getting all statistics */
            fParticleId = kWrongPOIid;
            inclusivepid = false;
          }
        }
      }
    }
    inclusivepid = inclusivepid && (fParticleId != kWrongPOIid);
  } else {
    /* if the track is accepted then by default is a hadron */
    fParticleId = kPOIh;
  }

  /* and now the exclusive ones */
  for (Int_t ix = 0; ix < fExclusiveCuts.GetEntriesFast(); ix++) {
    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fExclusiveCuts[ix])->IsTrueTrackAccepted(itrk);
    /* if track is accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast(),cutaccepted);
    exclusive = exclusive || cutaccepted;
  }

  /* now decide if accepted or not */
  accepted = inclusivetrack && inclusivepid && !exclusive;
  return accepted;
}

/// Check whether the true particle associated to a reconstructed track is primary
/// \param trk the reconstructed track
/// \return kTRUE if the associated particle is primary kFALSE otherwise
Bool_t AliCSTrackSelection::IsTruePrimary(AliVTrack *trk) {

  return AliCSTrackCuts::IsTruePrimary(trk);
}


/// Check whether the passed track label is associated to a MC injected signal
/// \param itrk label of the track to analyze whether it is associated to an injected signal
/// \return kTRUE if the track is associated to an injected signal kFALSE otherwise
Bool_t AliCSTrackSelection::IsFromMCInjectedSignal(Int_t itrk) {

  Bool_t isassociated = kFALSE;

  if (fHighestPrimaryLabel > 0) {
    AliMCEventHandler *eventHandler = AliCSAnalysisCutsBase::GetMCEventHandler();

    if (eventHandler != NULL) {
      /* MC ESD data */
      Int_t label = itrk;
      AliMCEvent* mcevent = eventHandler->MCEvent();
      AliVParticle *mother = mcevent->GetTrack(label);

      /* we have to find the primary one */
      while (!mcevent->IsPhysicalPrimary(label)) {
        label = mother->GetMother();
        if (label < 0) break;
        mother = mcevent->GetTrack(label);
        if (mother == NULL) break;
      }

      if (!(label < 0) && (mother != NULL)) {
        if (!(label < fHighestPrimaryLabel)) {
          isassociated = kTRUE;
        }
      }
    }
    else {
      /* MC ESD data */
      Int_t label = itrk;
      TClonesArray *arrayMC = AliCSAnalysisCutsBase::GetMCTrueArray();
      AliVParticle *mother = (AliVParticle *) arrayMC->At(label);

      /* we have to find the primary one */
      while ((mother != NULL) && !(mother->IsPhysicalPrimary())) {
        label = mother->GetMother();
        if (label < 0) break;
        mother = (AliVParticle *) arrayMC->At(label);
      }

      if (!(label < 0) && (mother != NULL)) {
        if (!(label < fHighestPrimaryLabel)) {
          isassociated = kTRUE;
        }
      }
    }
  }
  return isassociated;
}


/// Initializes the cuts
///
/// Inclusive and exclusive set of cuts are created and the needed histograms list
/// is allocated if needed the order is forwarded to the inclusive and exclusive set of cuts.
/// \param name an additional name to precede the cuts string
void AliCSTrackSelection::InitCuts(const char *name){

  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TList* tmplist = new TList();
  fPoiNames.clear();
  fPoiNames.push_back("Ha");

  auto addcut = [&](auto& lst, auto& cut, auto cfstr, std::string kind) {
    cut->InitializeCutsFromCutString(cfstr);
    lst.Add(cut);
    cut->SetQALevelOutput(fQALevel);
    cut->InitCuts(Form("%s cut", kind.c_str()));
    if (cut->GetHistogramsList() != NULL)
      tmplist->Add(cut->GetHistogramsList());
  };

  auto populatetrks = [&](auto& lst, std::string cutrole, auto& cutscf) {
    for (int icut = 0; icut < cutscf.GetEntries(); ++icut) {
      AliCSTrackCuts* cut = new AliCSTrackCuts(Form("CS_%sTrackCut%d", cutrole.c_str(), icut), Form("CS %s Track Cut %d", cutrole.c_str(), icut));
      addcut(lst, cut, ((TObjString*)cutscf.At(icut))->String().Data(), Form("%s track", cutrole.c_str()));
    }
  };

  auto populatepids = [&](auto& lst, std::string cutrole, auto& cutscf, bool updatepoilst = false) {
    int ncuts[AliPID::kSPECIESC] = {0};
    for (int icut = 0; icut < cutscf.GetEntries(); ++icut) {
      TString szcut = ((TObjString*)cutscf.At(icut))->String().Data();
      AliPID::EParticleType target = AliPID::kUnknown;
      Int_t delchars = 0;
      if (szcut.BeginsWith("e")) {
        delchars = 1;
        target = AliPID::kElectron;
      } else if (szcut.BeginsWith("mu")) {
        delchars = 2;
        target = AliPID::kMuon;
      } else if (szcut.BeginsWith("pi")) {
        delchars = 2;
        target = AliPID::kPion;
      } else if (szcut.BeginsWith("p")) {
        delchars = 1;
        target = AliPID::kProton;
      } else if (szcut.BeginsWith("k")) {
        delchars = 1;
        target = AliPID::kKaon;
      } else {
        AliFatal("Inconsistent family in track selection PID cuts string. ABORTING!!!");
      }
      szcut.Remove(0, delchars);
      ncuts[target]++;
      AliCSPIDCuts* cut = new AliCSPIDCuts(Form("CS_%sPidCut%d", cutrole.c_str(), icut), Form("CS %s PID Cut %d", cutrole.c_str(), icut), target, ncuts[target]);
      addcut(lst, cut, szcut.Data(), Form("%s PID", cutrole.c_str()));
    }
    static std::vector<std::string> speciesNames{"El", "Mu", "Pi", "Ka", "Pr"};
    if (updatepoilst) {
      for (uint isp = 0; isp < AliPID::kSPECIES; ++isp) {
        if (ncuts[isp] > 0) {
          fPoiNames.push_back(speciesNames[isp]);
        }
      }
    }
  };

  populatetrks(fInclusiveTrackCuts, "Inclusive", fInclusiveCutsStrings);
  populatetrks(fExclusiveCuts, "Exclusive", fExclusiveCutsStrings);
  populatepids(fInclusivePIDCuts, "Inclusive", fInclusivePidCutsStrings, true);
  populatepids(fExclusiveCuts, "Exclusive", fExclusivePidCutsStrings);

  if (name == NULL) name = GetName();

  if (fQALevel > AliCSAnalysisCutsBase::kQALevelNone || tmplist->GetEntries() != 0) {
    if(fHistogramsList != NULL)
      delete fHistogramsList;

    fHistogramsList = tmplist;
    fHistogramsList->SetOwner(kTRUE);
    fHistogramsList->SetName(name);
  }
  else {
    delete tmplist;
  }
  TH1::AddDirectory(oldstatus);
}

/// Allocates the different histograms if needed
///
/// It is supposed that the current cuts string is the running one
void AliCSTrackSelection::DefineHistograms(){

  if (fQALevel > AliCSAnalysisCutsBase::kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    /* the original name is used as title for the statistics histogram so, preserve it */
    TString originalTempName = fHistogramsList->GetName();
    fHistogramsList->SetName(Form("%s_%s",fHistogramsList->GetName(),fCutsString.Data()));

    Int_t nNoOfCuts = fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()+fExclusiveCuts.GetEntriesFast();
    fhCutsStatistics = new TH1F(Form("CutsStatistics_%s",fCutsString.Data()),Form("%s tracks cuts statistics",originalTempName.Data()),nNoOfCuts+4,-0.5,nNoOfCuts+3.5);
    fhCutsStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
    fhCutsStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");
    for (Int_t i = 0; i < fInclusiveTrackCuts.GetEntriesFast(); i++)
      fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, fInclusiveTrackCuts[i]->GetName());
    for (Int_t i = 0; i < fInclusivePIDCuts.GetEntriesFast(); i++)
      fhCutsStatistics->GetXaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+i+4, fInclusivePIDCuts[i]->GetName());
    for (Int_t i = 0; i < fExclusiveCuts.GetEntriesFast(); i++)
      fhCutsStatistics->GetXaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()+i+4, fExclusiveCuts[i]->GetName());
    fHistogramsList->Add(fhCutsStatistics);

    if(fQALevel == AliCSAnalysisCutsBase::kQALevelHeavy){
      fhCutsCorrelation = new TH2F(Form("CutCorrelation_%s",fCutsString.Data()),"Cuts correlation",nNoOfCuts+2,-0.5,nNoOfCuts+1.5,nNoOfCuts+2,-0.5,nNoOfCuts+1.5);
      for (Int_t i = 0; i < fInclusiveTrackCuts.GetEntriesFast(); i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(i+2, fInclusiveTrackCuts[i]->GetName());
        fhCutsCorrelation->GetYaxis()->SetBinLabel(i+2, fInclusiveTrackCuts[i]->GetName());
      }
      for (Int_t i = 0; i < fInclusivePIDCuts.GetEntriesFast(); i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+i+2, fInclusivePIDCuts[i]->GetName());
        fhCutsCorrelation->GetYaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+i+2, fInclusivePIDCuts[i]->GetName());
      }
      for (Int_t i = 0; i < fExclusiveCuts.GetEntriesFast(); i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()+i+2, fExclusiveCuts[i]->GetName());
        fhCutsCorrelation->GetYaxis()->SetBinLabel(fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()+i+2, fExclusiveCuts[i]->GetName());
      }
      fHistogramsList->Add(fhCutsCorrelation);
    }

    fhPtVsDCAxy[0]  = new TH2F(Form("PtVsDCAxyPtB_%s",fCutsString.Data()),"p_{T} vs DCA_{XY} before;DCA_{XY} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fhPtVsDCAxy[1]  = new TH2F(Form("PtVsDCAxyPtA_%s",fCutsString.Data()),"p_{T} vs DCA_{XY};DCA_{XY} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsDCAxy[0]);
    fHistogramsList->Add(fhPtVsDCAxy[1]);

    fhPtVsDCAz[0]  = new TH2F(Form("PtVsDCAzPtB_%s",fCutsString.Data()),"p_{T} vs DCA_{Z} before;DCA_{Z} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fhPtVsDCAz[1]  = new TH2F(Form("PtVsDCAzPtA_%s",fCutsString.Data()),"p_{T} vs DCA_{Z};DCA_{Z} (cm);p_{T} (GeV/c)",800,-4.0,4.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsDCAz[0]);
    fHistogramsList->Add(fhPtVsDCAz[1]);

    fhPtVsTPCCls[0]  = new TH2F(Form("PtVsTPCClustersB_%s",fCutsString.Data()),"p_{T} vs no. of TPC clusters before;no of clusters;p_{T} (GeV/c)",170,0,170,400,0.,10.);
    fhPtVsTPCCls[1]  = new TH2F(Form("PtVsTPCClustersA_%s",fCutsString.Data()),"p_{T} vs no. of TPC clusters;no of clusters;p_{T} (GeV/c)",170,0,170,400,0.,10.);
    fHistogramsList->Add(fhPtVsTPCCls[0]);
    fHistogramsList->Add(fhPtVsTPCCls[1]);

    fhPtVsTPCRows[0] = new TH2F(Form("PtVsTPCRowsB_%s", fCutsString.Data()), "p_{T} vs no. of TPC crossed rows before;no of crossed rows;p_{T} (GeV/c)", 170, 0, 170, 400, 0., 10.);
    fhPtVsTPCRows[1] = new TH2F(Form("PtVsTPCRowsA_%s", fCutsString.Data()), "p_{T} vs no. of TPC crossed rows;no of crossed rows;p_{T} (GeV/c)", 170, 0, 170, 400, 0., 10.);
    fHistogramsList->Add(fhPtVsTPCRows[0]);
    fHistogramsList->Add(fhPtVsTPCRows[1]);

    fhPtVsTPCRowOverFindCls[0] = new TH2F(Form("PtVsTPCCROFCB_%s", fCutsString.Data()),
                                          "p_{T} vs TPC crossed rows findable clusters ratio before;crossed rows / findable clusters;p_{T} (GeV/c)", 100, 0, 1, 400, 0., 10.);
    fhPtVsTPCRowOverFindCls[1] = new TH2F(Form("PtVsTPCCROFCA_%s",fCutsString.Data()),
        "p_{T} vs TPC crossed rows findable clusters ratio;crossed rows / findable clusters;p_{T} (GeV/c)",100,0,1,400,0.,10.);
    fHistogramsList->Add(fhPtVsTPCRowOverFindCls[0]);
    fHistogramsList->Add(fhPtVsTPCRowOverFindCls[1]);
    fhPtVsEta[0] = new TH2F(Form("PtVsEtaB_%s",fCutsString.Data()),"p_{T} vs #eta before;#eta;p_{T} (GeV/c)",100,-2.0,2.0,400,0.,10.);
    fhPtVsEta[1] = new TH2F(Form("PtVsEtaA_%s",fCutsString.Data()),"p_{T} vs #eta;#eta;p_{T} (GeV/c)",100,-2.0,2.0,400,0.,10.);
    fHistogramsList->Add(fhPtVsEta[0]);
    fHistogramsList->Add(fhPtVsEta[1]);

    if(fQALevel == AliCSAnalysisCutsBase::kQALevelHeavy){
      fhEtaVsPhi[0] = new TH2F(Form("EtaVsPhiB_%s",fCutsString.Data()),"#eta vs #phi before;#phi;#eta", 360, 0.0, 360.0, 100, -2.0, 2.0);
      fhEtaVsPhi[1] = new TH2F(Form("EtaVsPhiA_%s",fCutsString.Data()),"#eta vs #phi;#phi;#eta", 360, 0.0, 360.0, 100, -2.0, 2.0);
      fHistogramsList->Add(fhEtaVsPhi[0]);
      fHistogramsList->Add(fhEtaVsPhi[1]);
    }

    if (fInclusivePidCutsStrings.GetEntries() != 0 || fExclusivePidCutsStrings.GetEntries() != 0) {
      auto name = [](int idx) {
        switch (idx) {
          case 0:
            return "pi";
            break;
          case 1:
            return "k";
            break;
          case 2:
            return "p";
            break;
        }
        return "WRONG";
      };
      auto title = [](int idx) {
        switch (idx) {
          case 0:
            return "#pi";
            break;
          case 1:
            return "k";
            break;
          case 2:
            return "p";
            break;
        }
        return "WRONG";
      };

      /* build the P bins */
      const Int_t nPbins = 150;
      Double_t minP = 0.05;
      Double_t maxP = 20.0;
      Double_t *edges = new Double_t [nPbins + 1];
      Double_t factor = TMath::Power(maxP / minP, 1. / nPbins);
      edges[0] = minP; for (Int_t bin = 0; bin < nPbins; bin++) edges[bin+1] = factor * edges[bin];

      fhITSdEdxSignalVsP[0] = new TH2F(Form("ITSdEdxSignalB_%s",fCutsString.Data()),"ITS dE/dx signal before;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200);
      fhITSdEdxSignalVsP[1] = new TH2F(Form("ITSdEdxSignalA_%s",fCutsString.Data()),"ITS dE/dx signal;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200);
      fHistogramsList->Add(fhITSdEdxSignalVsP[0]);
      fHistogramsList->Add(fhITSdEdxSignalVsP[1]);

      fhTPCdEdxSignalVsP[0] = new TH2F(Form("TPCdEdxSignalB_%s",fCutsString.Data()),"TPC dE/dx signal before;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200.0);
      fhTPCdEdxSignalVsP[1] = new TH2F(Form("TPCdEdxSignalA_%s",fCutsString.Data()),"TPC dE/dx signal;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200.0);
      fHistogramsList->Add(fhTPCdEdxSignalVsP[0]);
      fHistogramsList->Add(fhTPCdEdxSignalVsP[1]);

      if (fInclusivePidCutsStrings.GetEntries() != 0) {
        for (int i = 0; i < 3; ++i) {
          fhTPCdEdxSelSignalVsP[i] = new TH2F(Form("TPCdEdxSel%sSignal_%s", name(i), fCutsString.Data()), Form("TPC dE/dx signal for selected %s;P (GeV/c);#frac{dE}{dx}_{%s} (au)", title(i), title(i)), nPbins, edges, 800, 0.0, 200.0);
          fHistogramsList->Add(fhTPCdEdxSelSignalVsP[i]);
        }
      }

      fhTOFSignalVsP[0] = new TH2F(Form("TOFSignalB_%s",fCutsString.Data()),"TOF signal before;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
      fhTOFSignalVsP[1] = new TH2F(Form("TOFSignalA_%s",fCutsString.Data()),"TOF signal;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
      fHistogramsList->Add(fhTOFSignalVsP[0]);
      fHistogramsList->Add(fhTOFSignalVsP[1]);

      fhPvsTOFMassSq[0] = new TH2F(Form("PvsMsqB_%s", fCutsString.Data()), "Momentum versus #it{m}^{2} before;#it{m}^{2} ((GeV/c^{2})^{2});P (GeV/c)", 140, 0.0, 1.4, nPbins, edges);
      fhPvsTOFMassSq[1] = new TH2F(Form("PvsMsqA_%s", fCutsString.Data()), "Momentum versus #it{m}^{2};#it{m}^{2} ((GeV/c^{2})^{2});P (GeV/c)", 140, 0.0, 1.4, nPbins, edges);
      fHistogramsList->Add(fhPvsTOFMassSq[0]);
      fHistogramsList->Add(fhPvsTOFMassSq[1]);

      if (fInclusivePidCutsStrings.GetEntries() != 0) {
        for (int i = 0; i < 3; ++i) {
          fhTOFSelSignalVsP[i] = new TH2F(Form("TOFSel%sSignal_%s", name(i), fCutsString.Data()), Form("TOF signal for selected %s;P (GeV/c);#beta_{%s}", title(i), title(i)), nPbins, edges, 400, 0.0, 1.1);
          fhSelPvsTOFMassSq[i] = new TH2F(Form("Sel%sPvsMsq_%s", name(i), fCutsString.Data()), Form("Momentum versus #it{m}^{2} for selected %s;#it{m}^{2}_{%s} ((GeV/c^{2})^{2});P (GeV/c)", title(i), title(i)), 140, 0.0, 1.4, nPbins, edges);
          fHistogramsList->Add(fhTOFSelSignalVsP[i]);
          fHistogramsList->Add(fhSelPvsTOFMassSq[i]);
        }
      }

      if (fQALevel == AliCSAnalysisCutsBase::kQALevelHeavy) {
        for (int i = 0; i < 3; ++i) {
          fhTPCTOFSigmaVsP[i][0] = new TH3F(Form("TPCTOFSigma%sVsPB_%s", name(i), fCutsString.Data()), Form("n#sigma to the %s line before;n#sigma_{TPC}^{%s};n#sigma_{TOF}^{%s}", title(i), title(i), title(i)),
                                            120, -6.0, 6.0, 120, -6.0, 6.0, nPbins, minP, maxP);
          fhTPCTOFSigmaVsP[i][1] = new TH3F(Form("TPCTOFSigma%sVsPA_%s", name(i), fCutsString.Data()), Form("n#sigma to the %s line;n#sigma_{TPC}^{%s};n#sigma_{TOF}^{%s}", title(i), title(i), title(i)),
                                            120, -6.0, 6.0, 120, -6.0, 6.0, nPbins, minP, maxP);
          fhTPCTOFSigmaVsP[i][0]->GetZaxis()->Set(nPbins, edges);
          fhTPCTOFSigmaVsP[i][1]->GetZaxis()->Set(nPbins, edges);
          fHistogramsList->Add(fhTPCTOFSigmaVsP[i][0]);
          fHistogramsList->Add(fhTPCTOFSigmaVsP[i][1]);

          fhTPCdEdxSignalDiffVsP[i][0] = new TH2F(Form("TPCdEdxSignal%sDiffVsPB_%s", name(i), fCutsString.Data()), Form("TPC dE/dx to the %s line before;P (GeV/c);dE/dx - <dE/dx>_{%s}", title(i), title(i)),
                                                  nPbins, edges, 800, -200.0, 200.0);
          fhTPCdEdxSignalDiffVsP[i][1] = new TH2F(Form("TPCdEdxSignal%sDiffVsPA_%s", name(i), fCutsString.Data()), Form("TPC dE/dx to the %s line;P (GeV/c);dE/dx - <dE/dx>_{%s}", title(i), title(i)),
                                                  nPbins, edges, 800, -200.0, 200.0);
          fHistogramsList->Add(fhTPCdEdxSignalDiffVsP[i][0]);
          fHistogramsList->Add(fhTPCdEdxSignalDiffVsP[i][1]);

          fhTOFSignalDiffVsP[i][0] = new TH2F(Form("TOFSignal%sDiffVsPB_%s", name(i), fCutsString.Data()), Form("TOF #beta to the %s line before;P (GeV/c);#beta - <#beta>_{%s}", title(i), title(i)),
                                              nPbins, edges, 400, -1.1, 1.1);
          fhTOFSignalDiffVsP[i][1] = new TH2F(Form("TOFSignal%sDiffVsPA_%s", name(i), fCutsString.Data()), Form("TOF #beta to the %s line;P (GeV/c);#beta - <#beta>_{%s}", title(i), title(i)),
                                              nPbins, edges, 400, -1.1, 1.1);
          fHistogramsList->Add(fhTOFSignalDiffVsP[i][0]);
          fHistogramsList->Add(fhTOFSignalDiffVsP[i][1]);

          if (fInclusivePidCutsStrings.GetEntries() != 0) {
            for (int j = 0; j < 3; ++j) {
              fhTPCdEdxSelSignalDiffVsP[i][j] = new TH2F(Form("TPCdEdxSel%sSignal%sDiffVsP_%s", name(j), name(i), fCutsString.Data()), Form("TPC dE/dx of selected %s to the %s line;P (GeV/c);dE/dx_{%s} - <dE/dx>_{%s}", title(j), title(i), title(j), title(i)),
                                                         nPbins, edges, 800, -200.0, 200.0);
              fHistogramsList->Add(fhTPCdEdxSelSignalDiffVsP[i][j]);
              fhTOFSelSignalDiffVsP[i][j] = new TH2F(Form("TOFSel%sSignal%sDiffVsP_%s", name(j), name(i), fCutsString.Data()), Form("TOF #beta of selected %s to the %s line;P (GeV/c);#beta_{%s} - <#beta>_{%s}", title(j), title(i), title(j), title(i)),
                                                     nPbins, edges, 400, -1.1, 1.1);
              fHistogramsList->Add(fhTOFSelSignalDiffVsP[i][j]);
            }
          }
        }
      }
    }

    TH1::AddDirectory(oldstatus);
  }
}

/// \cond CLASSIMP
ClassImp(AliCSTrackSelection);
/// \endcond
