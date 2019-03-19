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

/// \file AliCSTrackSelection.cxx
/// \brief Implementation of track selection class within the correlation studies analysis

/// Default constructor for serialization
AliCSTrackSelection::AliCSTrackSelection() :
  TNamed(),
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

  for (Int_t i = 0; i < 2; i++) {
    fhPtVsDCAxy[i] = NULL;
    fhPtVsDCAz[i] = NULL;
    fhPtVsTPCCls[i] = NULL;
    fhPtVsTPCRowOverFindCls[i] = NULL;
    fhEtaVsPhi[i] = NULL;
    fhPtVsEta[i] = NULL;
    fhITSdEdxSignalVsP[i] = NULL;
    fhTPCdEdxSignalVsP[i] = NULL;
    fhTOFSignalVsP[i] = NULL;
  }
}

/// Constructor
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSTrackSelection::AliCSTrackSelection(const char *name, const char *title) :
  TNamed(name,title),
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

  for (Int_t i = 0; i < 2; i++) {
    fhPtVsDCAxy[i] = NULL;
    fhPtVsDCAz[i] = NULL;
    fhPtVsTPCCls[i] = NULL;
    fhPtVsTPCRowOverFindCls[i] = NULL;
    fhEtaVsPhi[i] = NULL;
    fhPtVsEta[i] = NULL;
    fhITSdEdxSignalVsP[i] = NULL;
    fhTPCdEdxSignalVsP[i] = NULL;
    fhTOFSignalVsP[i] = NULL;
  }
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

    /* so far only reject e- after track selection */
    TObjArray *tokens = sztmpid.Tokenize("+");
    for (Int_t icut = 0; icut < tokens->GetEntries(); icut++) {
      if (((TObjString*) tokens->At(icut))->String().BeginsWith("-")) {
        TObjString *posz = new TObjString(((TObjString*) tokens->At(icut))->String().Remove(0,1).Data());
        if (posz->String().BeginsWith("e") || posz->String().BeginsWith("mu") || posz->String().BeginsWith("pi")
            || posz->String().BeginsWith("p") || posz->String().BeginsWith("k")) {
          fExclusivePidCutsStrings.Add(posz);
        }
        else {
          AliFatal("Wrong family in track selection PID cuts string. ABORTING!!!");
          return kFALSE;
        }
      }
      else {
        TObjString *posz = new TObjString(((TObjString*) tokens->At(icut))->String().Data());
        if (posz->String().BeginsWith("e") || posz->String().BeginsWith("mu") || posz->String().BeginsWith("pi")
            || posz->String().BeginsWith("p") || posz->String().BeginsWith("k")) {
          fInclusivePidCutsStrings.Add(posz);
        }
        else {
          AliFatal("Wrong family in track selection PID cuts string. ABORTING!!!");
          return kFALSE;
        }
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
      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (manager->GetInputEventHandler());
      fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
      /* if we need PID response instance and it is not there we cannot continue */
      if (fPIDResponse == NULL)
        AliFatal("No PID response instance. ABORTING!!!");
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
  Bool_t inclusivepid = (fInclusivePIDCuts.GetEntries() == 0);
  Bool_t exclusive = kFALSE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask->ResetAllBits();

  /* Check first the inclusive set of track cuts */
  for (Int_t ix = 0; ix < fInclusiveTrackCuts.GetEntriesFast(); ix++) {

    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fInclusiveTrackCuts[ix])->IsTrackAccepted(trk);
    /* if track is not accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix,!cutaccepted);
    inclusivetrack = inclusivetrack || cutaccepted;
  }

  /* now the inclusive set of pid cuts */
  for (Int_t ix = 0; ix < fInclusivePIDCuts.GetEntriesFast(); ix++) {

    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fInclusivePIDCuts[ix])->IsTrackAccepted(trk);
    /* if track is not accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast(),!cutaccepted);
    inclusivepid = inclusivepid || cutaccepted;
  }

  /* and now the exclusive ones */
  for (Int_t ix = 0; ix < fExclusiveCuts.GetEntriesFast(); ix++) {
    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fExclusiveCuts[ix])->IsTrackAccepted(trk);
    /* if track is accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast(),cutaccepted);
    exclusive = exclusive || cutaccepted;
  }

  /* now decide if accepted or not */
  accepted = inclusivetrack && inclusivepid && !exclusive;

  if (fQALevel > AliCSAnalysisCutsBase::kQALevelNone) {
    /* let's fill the histograms */
    fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n tracks")));
    if (!accepted)
      fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n cut tracks")));

    for (Int_t ix=0; ix<fInclusiveTrackCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindBin(fInclusiveTrackCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fInclusiveTrackCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix))
        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fInclusiveTrackCuts[ix]->GetName())));

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx=ix; jx<fInclusiveTrackCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fInclusiveTrackCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusiveTrackCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix) && fCutsActivatedMask->TestBitNumber(jx)) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fInclusiveTrackCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fInclusiveTrackCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
        for (Int_t jx=ix; jx<fInclusivePIDCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fInclusivePIDCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusivePIDCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix) && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast())) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fInclusiveTrackCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fInclusivePIDCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
        for (Int_t jx=0; jx<fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix)
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fInclusiveTrackCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fExclusiveCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }
    for (Int_t ix=0; ix<fInclusivePIDCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindBin(fInclusivePIDCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fInclusivePIDCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()))
        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fInclusivePIDCuts[ix]->GetName())));

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx=ix; jx<fInclusivePIDCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fInclusivePIDCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fInclusivePIDCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast())) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fInclusivePIDCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fInclusivePIDCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
        for (Int_t jx=0; jx<fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fInclusivePIDCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fExclusiveCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }
    for (Int_t ix=0; ix<fExclusiveCuts.GetEntriesFast(); ix++) {
      if (fhCutsStatistics->GetXaxis()->FindBin(fExclusiveCuts[ix]->GetName()) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", ix, fExclusiveCuts[ix]->GetName()));

      if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast()))
        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fExclusiveCuts[ix]->GetName())));

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        for (Int_t jx=0; jx<fExclusiveCuts.GetEntriesFast(); jx++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fExclusiveCuts[jx]->GetName()) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", jx, fExclusiveCuts[jx]->GetName()));

          if (fCutsActivatedMask->TestBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())
              && fCutsActivatedMask->TestBitNumber(jx+fInclusiveTrackCuts.GetEntriesFast()+fInclusivePIDCuts.GetEntriesFast())) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fExclusiveCuts[ix]->GetName()));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fExclusiveCuts[jx]->GetName()));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }

    /* get some values needed for histograms filling */
    Float_t dca[2], bCov[3];
    if (trk->IsA() == AliESDtrack::Class()) {
      trk->GetImpactParameters(dca,bCov);
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
      fhPtVsTPCRowOverFindCls[i]->Fill(ratioCrossedRowsOverFindableClustersTPC,trk->Pt());
      fhPtVsEta[i]->Fill(trk->Eta(),trk->Pt());

      if (fQALevel > AliCSAnalysisCutsBase::kQALevelLight) {
        fhEtaVsPhi[i]->Fill(trk->Phi()*180.0/TMath::Pi(),trk->Eta());
      }

      if (fInclusivePidCutsStrings.GetEntries() != 0 || fExclusivePidCutsStrings.GetEntries() != 0) {
        /* don't fill before if not requested */
        if (i == 1 || (fQALevel > AliCSAnalysisCutsBase::kQALevelLight)) {
          fhITSdEdxSignalVsP[i]->Fill(trk->P(),ttrk->GetITSsignal());
          fhTPCdEdxSignalVsP[i]->Fill(trk->P(),TMath::Abs(ttrk->GetTPCsignal()));
          if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
            static const Double_t c_cm_ps = TMath::C() * 1.0e2 * 1.0e-12;
            Double_t tracklen_cm = trk->GetIntegratedLength();
            Double_t toftime_ps = ttrk->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(trk->P());
            Double_t beta = tracklen_cm / toftime_ps / c_cm_ps;
            fhTOFSignalVsP[i]->Fill(trk->P(),beta);
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
  Bool_t inclusivepid = (fInclusivePIDCuts.GetEntriesFast() == 0);
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

  /* now the inclusive set of pid cuts */
  for (Int_t ix = 0; ix < fInclusivePIDCuts.GetEntriesFast(); ix++) {

    Bool_t cutaccepted = ((AliCSTrackCutsBase *) fInclusivePIDCuts[ix])->IsTrueTrackAccepted(itrk);
    /* if track is not accepted the cut is activated */
    fCutsActivatedMask->SetBitNumber(ix+fInclusiveTrackCuts.GetEntriesFast(),!cutaccepted);
    inclusivepid = inclusivepid || cutaccepted;
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

  TList *tmplist = new TList();

  for (Int_t icut = 0; icut < fInclusiveCutsStrings.GetEntries(); icut++) {
    AliCSTrackCuts *cut = new AliCSTrackCuts(Form("CS_InclusiveTrackCut%d",icut),Form("CS Inclusive Track Cut %d",icut));
    cut->InitializeCutsFromCutString(((TObjString*)fInclusiveCutsStrings.At(icut))->String().Data());
    fInclusiveTrackCuts.Add(cut);
    cut->SetQALevelOutput(fQALevel);
    cut->InitCuts("Inclusive track cut");
    if (cut->GetHistogramsList() != NULL)
      tmplist->Add(cut->GetHistogramsList());
  }

  for (Int_t icut = 0; icut < fInclusivePidCutsStrings.GetEntries(); icut++) {
    TString szcut = ((TObjString*)fInclusivePidCutsStrings.At(icut))->String().Data();
    AliPID::EParticleType target = AliPID::kUnknown;
    Int_t delchars = 0;
    if (szcut.BeginsWith("e")) { delchars = 1; target = AliPID::kElectron; }
    else if (szcut.BeginsWith("mu")) { delchars = 2; target = AliPID::kMuon; }
    else if (szcut.BeginsWith("pi")) { delchars = 2; target = AliPID::kPion; }
    else if (szcut.BeginsWith("p")) { delchars = 1; target = AliPID::kProton; }
    else if (szcut.BeginsWith("k")) { delchars = 1; target = AliPID::kKaon; }
    else {
      AliFatal("Inconsistent family in track selection PID cuts string. ABORTING!!!");
    }
    szcut.Remove(0,delchars);
    AliCSPIDCuts *cut = new AliCSPIDCuts(Form("CS_InclusivePidCut%d",icut),Form("CS Inclusive PID Cut %d",icut),target);
    cut->InitializeCutsFromCutString(szcut.Data());
    fInclusivePIDCuts.Add(cut);
    cut->SetQALevelOutput(fQALevel);
    cut->InitCuts("Inclusive pid cut");
    if (cut->GetHistogramsList() != NULL)
      tmplist->Add(cut->GetHistogramsList());
  }

  for (Int_t icut = 0; icut < fExclusiveCutsStrings.GetEntries(); icut++) {
    AliCSTrackCuts *cut = new AliCSTrackCuts(Form("CS_ExclusiveTrackCut%d",icut),Form("CS Exclusive Track Cut %d",icut));
    cut->InitializeCutsFromCutString(((TObjString*)fExclusiveCutsStrings.At(icut))->String().Data());
    fExclusiveCuts.Add(cut);
    cut->SetQALevelOutput(fQALevel);
    cut->InitCuts("Exclusive track cut");
    if (cut->GetHistogramsList() != NULL)
      tmplist->Add(cut->GetHistogramsList());
  }

  for (Int_t icut = 0; icut < fExclusivePidCutsStrings.GetEntries(); icut++) {
    TString szcut = ((TObjString*)fExclusivePidCutsStrings.At(icut))->String().Data();
    AliPID::EParticleType target  = AliPID::kUnknown;
    Int_t delchars = 0;
    if (szcut.BeginsWith("e")) { delchars = 1; target = AliPID::kElectron; }
    else if (szcut.BeginsWith("mu")) { delchars = 2; target = AliPID::kMuon; }
    else if (szcut.BeginsWith("pi")) { delchars = 2; target = AliPID::kPion; }
    else if (szcut.BeginsWith("p")) { delchars = 1; target = AliPID::kProton; }
    else if (szcut.BeginsWith("k")) { delchars = 1; target = AliPID::kKaon; }
    else {
      AliFatal("Inconsistent family in track selection PID cuts string. ABORTING!!!");
    }
    szcut.Remove(0,delchars);
    AliCSPIDCuts *cut = new AliCSPIDCuts(Form("CS_ExclusivePidCut%d",icut),Form("CS Exclusive PID Cut %d",icut),target);
    cut->InitializeCutsFromCutString(szcut.Data());
    fExclusiveCuts.Add(cut);
    cut->SetQALevelOutput(fQALevel);
    cut->InitCuts("Exclusive pid cut");
    if (cut->GetHistogramsList() != NULL)
      tmplist->Add(cut->GetHistogramsList());
  }

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

    fhPtVsTPCRowOverFindCls[0] = new TH2F(Form("PtVsTPCCROFCB_%s",fCutsString.Data()),
        "p_{T} vs TPC crossed rows findable clusters ratio before;crossed rows / findable clusters;p_{T} (GeV/c)",100,0,1,400,0.,10.);
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

      fhTOFSignalVsP[0] = new TH2F(Form("TOFSignalB_%s",fCutsString.Data()),"TOF signal before;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
      fhTOFSignalVsP[1] = new TH2F(Form("TOFSignalA_%s",fCutsString.Data()),"TOF signal;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
      fHistogramsList->Add(fhTOFSignalVsP[0]);
      fHistogramsList->Add(fhTOFSignalVsP[1]);
    }

    TH1::AddDirectory(oldstatus);
  }
}

/// \cond CLASSIMP
ClassImp(AliCSTrackSelection);
/// \endcond
