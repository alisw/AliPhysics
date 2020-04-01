/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskCorrelationsStudies.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#include "AliAnalysisTaskCorrelationsStudies.h"

#include <TBits.h>
#include "Riostream.h"
#include "TGrid.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3I.h"
#include "TProfile3D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TParameter.h"

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliCSTrackMaps.h"
#include "AliCSEventCuts.h"
#include "AliCSPIDCuts.h"
#include "AliDptDptCorrelations.h"
#include "AliCSPairAnalysis.h"
#include "AliVEvent.h"
#include "AliAODHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"

ClassImp(AliAnalysisTaskCorrelationsStudies)

const Int_t AliAnalysisTaskCorrelationsStudies::kgTHnDimension = 4;

AliCSTrackMaps aodTrackMaps("AODtrackMaps", "The AOD tracks id maps"); ///< the id track maps for AOD tracks

/// Default constructor
AliAnalysisTaskCorrelationsStudies::AliAnalysisTaskCorrelationsStudies() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutput(NULL),
    fTaskConfigurationString(""),
    fTaskActivitiesString(""),
    fEventCutsString(""),
    fEfficiencyProfileToEnforce(""),
    fEfficiencyProfileOnTrue(""),
    fAdditionalMCRecOption(""),
    fTrackSelectionString(""),
    fEventCuts(NULL),
    fTrackSelectionCuts(NULL),
    fDataPeriod(AliCSAnalysisCutsBase::kNoPeriod),
    fEnforceEfficiencyProfile(kFALSE),
    fOnTrueEfficiencyProfile(kNoEfficiencyOnTrue),
    fCorrectOnTrueEfficiency(kFALSE),
    fDoProcessCorrelations(kFALSE),
    fDoProcessPairAnalysis(kFALSE),
    fMCRecOption(kNone),
    fProcessCorrelations(),
    fProcessMCRecCorrelationsWithOptions(),
    fProcessTrueCorrelations(),
    fProcessPairAnalysis(),
    fProcessTruePairAnalysis(),
    ffEfficiencyProfile(NULL),
    fRandomGenerator(NULL),
    fhWeightsTrack_1(NULL),
    fhWeightsTrack_2(NULL),
    fhEffCorrTrack_1(NULL),
    fhEffCorrTrack_2(NULL),
    fhPairEfficiency_PP(NULL),
    fhPairEfficiency_PM(NULL),
    fhPairEfficiency_MM(NULL),
    fhPairEfficiency_MP(NULL),
    fPositiveTrackPdf(NULL),
    fNegativeTrackPdf(NULL),
    fTrueToRec(NULL),
    fTrueToRecWrong(NULL),
    fMCRecFlags(NULL),
    fMCTruePrimaryFlags(NULL),
    fMCFlagsStorageSize(30*1024),
    fhOnTrueEfficiencyProfile_1(NULL),
    fhOnTrueEfficiencyProfile_2(NULL),
    fhPt(NULL),
    fhPtPos(NULL),
    fhPtNeg(NULL),
    fhP(NULL),
    fhPPos(NULL),
    fhPNeg(NULL),
    fhTruePt(NULL),
    fhTruePtPos(NULL),
    fhTruePtNeg(NULL),
    fhTrueP(NULL),
    fhTruePPos(NULL),
    fhTruePNeg(NULL),
    fhPurePt(NULL),
    fhPurePtPos(NULL),
    fhPurePtNeg(NULL),
    fhPureP(NULL),
    fhPurePPos(NULL),
    fhPurePNeg(NULL),
    fhUnConstrainedPt(NULL),
    fhPtDifference(NULL),
    fhEtaB(NULL),
    fhEtaA(NULL),
    fhTrueEta(NULL),
    fhPhiB(NULL),
    fhPhiA(NULL),
    fhTruePhi(NULL),
    fhEtaVsPhiB(NULL),
    fhEtaVsPhiA(NULL),
    fhTrueEtaVsPhi(NULL),
    fhPtVsEtaB(NULL),
    fhPtVsEtaA(NULL),
    fhTruePtVsEta(NULL),
    fhPt3DB(NULL),
    fhPt3DA(NULL),
    fhTruePt3D(NULL),
    fh3Dn1B(NULL),
    fh3Dn1A(NULL),
    fhLambdaVsMultiplicity(NULL),
    fhAcceptedVsMultiplicity(NULL),
    fhTrueLambdaVsPrimaries(NULL),
    fhTrueAcceptedVsPrimaries(NULL),
    fhTrueRecNotAccepted(NULL),
    fhTrueNotReconstructed(NULL),
    fhTrueMultiRec(NULL),
    fhPtTrueMultiRec(NULL)
{
    // Dummy constructor ALWAYS needed for I/O.
}

/// Analysis task constructor
/// \param name the name to assign to the task
AliAnalysisTaskCorrelationsStudies::AliAnalysisTaskCorrelationsStudies(const char *name) // All data members should be initialized here
   :AliAnalysisTaskSE(name),
    fOutput(NULL),
    fTaskConfigurationString(""),
    fTaskActivitiesString(""),
    fEventCutsString(""),
    fEfficiencyProfileToEnforce(""),
    fEfficiencyProfileOnTrue(""),
    fAdditionalMCRecOption(""),
    fTrackSelectionString(""),
    fEventCuts(NULL),
    fTrackSelectionCuts(NULL),
    fDataPeriod(AliCSAnalysisCutsBase::kNoPeriod),
    fEnforceEfficiencyProfile(kFALSE),
    fOnTrueEfficiencyProfile(kNoEfficiencyOnTrue),
    fCorrectOnTrueEfficiency(kFALSE),
    fDoProcessCorrelations(kFALSE),
    fDoProcessPairAnalysis(kFALSE),
    fMCRecOption(kNone),
    fProcessCorrelations("DptDptCorrelations"),
    fProcessMCRecCorrelationsWithOptions("DptDptCorrelationsMCRecOptions"),
    fProcessTrueCorrelations("DptDptTrueCorrelations"),
    fProcessPairAnalysis("PairAnalysis"),
    fProcessTruePairAnalysis("TruePairAnalysis"),
    ffEfficiencyProfile(NULL),
    fRandomGenerator(NULL),
    fhWeightsTrack_1(NULL),
    fhWeightsTrack_2(NULL),
    fhEffCorrTrack_1(NULL),
    fhEffCorrTrack_2(NULL),
    fhPairEfficiency_PP(NULL),
    fhPairEfficiency_PM(NULL),
    fhPairEfficiency_MM(NULL),
    fhPairEfficiency_MP(NULL),
    fPositiveTrackPdf(NULL),
    fNegativeTrackPdf(NULL),
    fTrueToRec(NULL),
    fTrueToRecWrong(NULL),
    fMCRecFlags(NULL),
    fMCTruePrimaryFlags(NULL),
    fMCFlagsStorageSize(30*1024),
    fhOnTrueEfficiencyProfile_1(NULL),
    fhOnTrueEfficiencyProfile_2(NULL),
    fhPt(NULL),
    fhPtPos(NULL),
    fhPtNeg(NULL),
    fhP(NULL),
    fhPPos(NULL),
    fhPNeg(NULL),
    fhTruePt(NULL),
    fhTruePtPos(NULL),
    fhTruePtNeg(NULL),
    fhTrueP(NULL),
    fhTruePPos(NULL),
    fhTruePNeg(NULL),
    fhPurePt(NULL),
    fhPurePtPos(NULL),
    fhPurePtNeg(NULL),
    fhPureP(NULL),
    fhPurePPos(NULL),
    fhPurePNeg(NULL),
    fhUnConstrainedPt(NULL),
    fhPtDifference(NULL),
    fhEtaB(NULL),
    fhEtaA(NULL),
    fhTrueEta(NULL),
    fhPhiB(NULL),
    fhPhiA(NULL),
    fhTruePhi(NULL),
    fhEtaVsPhiB(NULL),
    fhEtaVsPhiA(NULL),
    fhTrueEtaVsPhi(NULL),
    fhPtVsEtaB(NULL),
    fhPtVsEtaA(NULL),
    fhTruePtVsEta(NULL),
    fhPt3DB(NULL),
    fhPt3DA(NULL),
    fhTruePt3D(NULL),
    fh3Dn1B(NULL),
    fh3Dn1A(NULL),
    fhLambdaVsMultiplicity(NULL),
    fhAcceptedVsMultiplicity(NULL),
    fhTrueLambdaVsPrimaries(NULL),
    fhTrueAcceptedVsPrimaries(NULL),
    fhTrueRecNotAccepted(NULL),
    fhTrueNotReconstructed(NULL),
    fhTrueMultiRec(NULL),
    fhPtTrueMultiRec(NULL)
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    DefineOutput(1, TList::Class());                                            // for output list
    DefineOutput(2, TList::Class());                                            // for output list
    DefineOutput(3, TList::Class());                                            // for output list
    DefineOutput(4, TList::Class());                                            // for output list
}

/// Default destructor
AliAnalysisTaskCorrelationsStudies::~AliAnalysisTaskCorrelationsStudies()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
    if (fhWeightsTrack_1 != NULL) delete fhWeightsTrack_1;
    if (fhWeightsTrack_2 != NULL) delete fhWeightsTrack_2;
    if (fhEffCorrTrack_1 != NULL) delete fhEffCorrTrack_1;
    if (fhEffCorrTrack_2 != NULL) delete fhEffCorrTrack_2;
    if (fPositiveTrackPdf != NULL) delete fPositiveTrackPdf;
    if (fNegativeTrackPdf != NULL) delete fNegativeTrackPdf;
    if (fTrueToRec != NULL) delete fTrueToRec;
    if (fTrueToRecWrong != NULL) delete fTrueToRecWrong;
    if (fMCRecFlags != NULL) delete [] fMCRecFlags;
    if (fMCTruePrimaryFlags != NULL) delete [] fMCTruePrimaryFlags;
    if (ffEfficiencyProfile != NULL) delete ffEfficiencyProfile;
    if (fhOnTrueEfficiencyProfile_1 != NULL) delete fhOnTrueEfficiencyProfile_1;
    if (fhOnTrueEfficiencyProfile_2 != NULL) delete fhOnTrueEfficiencyProfile_2;
    if (fRandomGenerator != NULL) delete fRandomGenerator;
    delete fEventCuts;
    delete fTrackSelectionCuts;
}

/// \brief Builds the efficiency profiles for both rec and true if needed
/// Information for the profiles is taken from the stored configuration
/// strings for both rec and true.
/// The configuration string could be a single expression depending on x (pT),
/// or a concatenation |bin,bin,expr| terms each of them defining an initial
/// bin, a final bin and an expression to assign the efficiency on that bin
/// interval. In this case the expression will depend on pT bin number. Again
/// the expression  will not depend on the pT value of the concrete bin being
/// handled. If pT value were needed use the corresponding factor within the
/// expression.
Bool_t AliAnalysisTaskCorrelationsStudies::BuildEfficiencyProfiles() {

  Bool_t done = kTRUE;

  /* let's get the limits of the profiles */
  Float_t ptlow = fhPt->GetXaxis()->GetBinLowEdge(1);
  Float_t ptup = fhPt->GetXaxis()->GetBinUpEdge(fhPt->GetNbinsX());

  /* let's build the efficiency profile formula if required */
  if (fEfficiencyProfileToEnforce.Length() != 0) {
    fEnforceEfficiencyProfile = kTRUE;
    if (ffEfficiencyProfile != NULL) delete ffEfficiencyProfile;
    if (fRandomGenerator != NULL) delete fRandomGenerator;

    ffEfficiencyProfile = new TF1("EfficiencyProfile",fEfficiencyProfileToEnforce.Data(),ptlow,ptup);
    fRandomGenerator = new TRandom3();
    AliInfo(Form("Configured efficiency profile enforcement with formula: %s", ffEfficiencyProfile->GetExpFormula().ReplaceAll("x","pT").Data()));
  }

  /* let's build the efficiency profile for true data if required */
  if (fEfficiencyProfileOnTrue.Length() != 0) {
    fOnTrueEfficiencyProfile = kEfficiencyProfile;
    if (fhOnTrueEfficiencyProfile_1 != NULL) delete fhOnTrueEfficiencyProfile_1;
    if (fhOnTrueEfficiencyProfile_2 != NULL) delete fhOnTrueEfficiencyProfile_2;
    if (fRandomGenerator != NULL) delete fRandomGenerator;

    /* let's check for a single expression for the whole pT range */
    TObjArray *intervals = fEfficiencyProfileOnTrue.Tokenize(TString("|"));
    if (intervals->GetEntries() > 1) {
      fhOnTrueEfficiencyProfile_1 = (TH2F*) fhPtVsEtaB->Clone("OnTrueEfficiencyProfileTrackOne");
      fhOnTrueEfficiencyProfile_2 = (TH2F*) fhPtVsEtaB->Clone("OnTrueEfficiencyProfileTrackTwo");
      fhOnTrueEfficiencyProfile_1->Reset();
      fhOnTrueEfficiencyProfile_2->Reset();

      for (Int_t entry = 0; entry < intervals->GetEntries(); entry++) {
        Int_t start, end;
        char expr[128];
        sscanf(((TObjString *)intervals->At(entry))->String().Data(), "%d,%d,%s", &start, &end, expr);
        AliInfo(Form("Assigning efficiency profile from bin %d to bin %d, with values: %s", start+1, end, expr));
        TF1 *formula = new TF1(Form("OnTrue%dFormula", entry),expr,start,end);
        for (Int_t ptbin = start; ptbin < end; ptbin++) {
          for (Int_t etabin = 0; etabin < fhOnTrueEfficiencyProfile_1->GetNbinsX(); etabin++) {
            fhOnTrueEfficiencyProfile_1->SetBinContent(etabin+1, ptbin+1,formula->Eval(ptbin));
            fhOnTrueEfficiencyProfile_2->SetBinContent(etabin+1, ptbin+1,formula->Eval(ptbin));
          }
        }
        delete formula;
      }
    }
    else {
      fhOnTrueEfficiencyProfile_1 = (TH2F*) fhPtVsEtaB->Clone("OnTrueEfficiencyProfileTrackOne");
      fhOnTrueEfficiencyProfile_2 = (TH2F*) fhPtVsEtaB->Clone("OnTrueEfficiencyProfileTrackTwo");
      fhOnTrueEfficiencyProfile_1->Reset();
      fhOnTrueEfficiencyProfile_2->Reset();

      /* we could have a file with two-dimensional histograms for efficiency */
      if (fEfficiencyProfileOnTrue.Contains("/")) {
        Bool_t oldstatus = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kFALSE);
        TFile *effprofile = TFile::Open(fEfficiencyProfileOnTrue.Data());
        if (effprofile != NULL && effprofile->IsOpen()) {
          fhOnTrueEfficiencyProfile_1 = (TH2F*)((TH2F*) effprofile->Get("OnTrueEfficiencyProfile_1"))->Clone("StoredOnTrueEfficiencyProfile_1");
          fhOnTrueEfficiencyProfile_2 = (TH2F*)((TH2F*) effprofile->Get("OnTrueEfficiencyProfile_2"))->Clone("StoredOnTrueEfficiencyProfile_2");
          effprofile->Close();
          delete effprofile;
          if (!(fhOnTrueEfficiencyProfile_1 != NULL) || !(fhOnTrueEfficiencyProfile_1 != NULL)) {
            AliFatal(Form("On true efficiency profiles not found in file %s. ABORTING!!!", fEfficiencyProfileOnTrue.Data()));
          }
        }
        else {
          AliFatal(Form("On true efficiency profile file %s not found. ABORTING!!!", fEfficiencyProfileOnTrue.Data()));
        }
        TH1::AddDirectory(oldstatus);
      }
      else {
        /* we could still have a single term, check it and process it */
        if (((TObjString *)intervals->At(0))->String().Contains(",")) {
          Int_t start, end;
          char expr[128];
          sscanf(((TObjString *)intervals->At(0))->String().Data(), "%d,%d,%s", &start, &end, expr);
          AliInfo(Form("Assigning efficiency profile from bin %d to bin %d, with values: %s", start+1, end, expr));
          TF1 *formula = new TF1(Form("OnTrue%dFormula", 0),expr,start,end);
          for (Int_t ptbin = start; ptbin < end; ptbin++) {
            for (Int_t etabin = 0; etabin < fhOnTrueEfficiencyProfile_1->GetNbinsX(); etabin++) {
              fhOnTrueEfficiencyProfile_1->SetBinContent(etabin+1, ptbin+1,formula->Eval(ptbin));
              fhOnTrueEfficiencyProfile_2->SetBinContent(etabin+1, ptbin+1,formula->Eval(ptbin));
            }
          }
          delete formula;
        }
        else {
          TF1 *formulaEfficiencyProfile = new TF1("OnTrueEfficiencyProfileFormula",fEfficiencyProfileOnTrue.Data(),ptlow,ptup);
          for (Int_t ptbin = 0; ptbin < fhOnTrueEfficiencyProfile_1->GetNbinsY(); ptbin++) {
            for (Int_t etabin = 0; etabin < fhOnTrueEfficiencyProfile_1->GetNbinsX(); etabin++) {
              fhOnTrueEfficiencyProfile_1->SetBinContent(etabin+1, ptbin+1,formulaEfficiencyProfile->Eval(ptbin));
              fhOnTrueEfficiencyProfile_2->SetBinContent(etabin+1, ptbin+1,formulaEfficiencyProfile->Eval(ptbin));
            }
          }
          delete formulaEfficiencyProfile;
        }
      }
    }
    delete intervals;
    AliInfo("Configured two-dimensional efficiency profile on true data for track one");
    for (Int_t ptbin = 0; ptbin < fhOnTrueEfficiencyProfile_1->GetNbinsY(); ptbin++) {
      printf("%f, ", fhOnTrueEfficiencyProfile_1->GetBinContent(4,ptbin+1));
      if ((ptbin+1)%8 == 0) printf("\n");
    }
    printf("\n");
    for (Int_t ptbin = 0; ptbin < fhOnTrueEfficiencyProfile_1->GetNbinsY(); ptbin++) {
      printf("%f, ", fhOnTrueEfficiencyProfile_1->GetBinContent(12,ptbin+1));
      if ((ptbin+1)%8 == 0) printf("\n");
    }
    printf("\n");
    AliInfo("Configured two-dimensional efficiency profile on true data for track two");
    for (Int_t ptbin = 0; ptbin < fhOnTrueEfficiencyProfile_2->GetNbinsY(); ptbin++) {
      printf("%f, ", fhOnTrueEfficiencyProfile_2->GetBinContent(4,ptbin+1));
      if ((ptbin+1)%8 == 0) printf("\n");
    }
    printf("\n");
    for (Int_t ptbin = 0; ptbin < fhOnTrueEfficiencyProfile_2->GetNbinsY(); ptbin++) {
      printf("%f, ", fhOnTrueEfficiencyProfile_2->GetBinContent(12,ptbin+1));
      if ((ptbin+1)%8 == 0) printf("\n");
    }
    printf("\n");
    fRandomGenerator = new TRandom3();

    if (fDoProcessCorrelations) {
      /* should we correct true for this efficiency profile */
      if (fCorrectOnTrueEfficiency) {
        /* yes, so apply it as efficiency correction */
//        fProcessTrueCorrelations.SetEfficiencyCorrection(fhOnTrueEfficiencyProfile_1,fhOnTrueEfficiencyProfile_2);
        AliFatal("Applying two-dimensional profile as efficiency correction on true data not yet supported");
        AliInfo("Applying the same profile as efficiency correction on true data");
      }
      else {
        /* true will not apply any efficiency correction */
        fProcessTrueCorrelations.SetEfficiencyCorrection(NULL,NULL);
      }
    }
  }
  else {
    if (fDoProcessCorrelations) {
      /* true will not apply any efficiency correction */
      fProcessTrueCorrelations.SetEfficiencyCorrection(NULL,NULL);
    }
  }

  if (fDoProcessCorrelations) {
    /* for the time being, rec with true will never apply efficiency correction */
    fProcessMCRecCorrelationsWithOptions.SetEfficiencyCorrection(NULL,NULL);
  }

  return done;
}

/// \brief Return the number of MC reconstructed tracks for the current event
/// Basically this is needed due to the fact that MC information is
/// differently stored in ESD than in AOD
///
/// \return the number of MC reconstructed tracks for the current event
Int_t AliAnalysisTaskCorrelationsStudies::GetNoOfMCRecTracks() {
  if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
    return fInputEvent->GetNumberOfTracks();
  }
  else {
    AliAODHeader *header = (AliAODHeader *) fInputEvent->GetHeader();
    return header->GetNumberOfESDTracks();
  }
}

/// \brief Return the number of true particles for the current event
/// Basically this is needed due to the fact that MC information is
/// differently stored in ESD than in AOD
///
/// \return the number of true particles for the current event
Int_t AliAnalysisTaskCorrelationsStudies::GetNoOfTrueParticles() {
  if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
    return fMCEvent->GetNumberOfTracks();
  }
  else {
    return AliCSAnalysisCutsBase::GetMCTrueArray()->GetEntries();
  }
}

/// \brief Return the number of true physical primary particles for
/// the current event
/// Basically this is needed due to the fact that MC information is
/// differently stored in ESD than in AOD
///
/// \return the number of true physical primary particles for the
/// current event
Int_t AliAnalysisTaskCorrelationsStudies::GetNoOfTruePrimaries() {
  if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
    return fMCEvent->GetNumberOfPrimaries();
  }
  else {
    TClonesArray *arrayMC = AliCSAnalysisCutsBase::GetMCTrueArray();
    Int_t nNoOfPrimaries = 0;

    for (Int_t itrk = 0; itrk < arrayMC->GetEntriesFast(); itrk++) {
      AliAODMCParticle *particle = (AliAODMCParticle *) arrayMC->At(itrk);
      if (particle != NULL && particle->IsPhysicalPrimary()) {
        nNoOfPrimaries++;
      }
    }
    return nNoOfPrimaries;
  }
}

/// \brief Builds the true to reconstructed relation to allow quick access
/// to both data sets
/// The relation is materialized as an array addressed by true label
/// or true track index value in the MC event. The object stored
/// corresponds to the reconstructed track which has as label the
/// index in the array under which it is stored.
/// There are two arrays one with the reconstructed stored when their
/// labels are positive (or zero) and the other with the reconstructed
/// stored when their labels are negative. The index corresponds to the
/// label absolute value
/// The additional array with flags are addressed by the track id. The
/// addressing takes always two steps due to the fact of a potential
/// presence of MC AOD format
void AliAnalysisTaskCorrelationsStudies::BuildTrueRecRelation() {
  /* allocate the structures if they are not there */
  if (fTrueToRec == NULL) {
    fTrueToRec = new TObjArray(30 * 1024); /* initial value for minimize the growing of the array */
    fTrueToRec->SetOwner(kFALSE); /* we will not own the rec tracks */
  }
  if (fTrueToRecWrong == NULL) {
    fTrueToRecWrong = new TObjArray(30 * 1024); /* initial value for minimize the growing of the array */
    fTrueToRecWrong->SetOwner(kFALSE); /* we will not own the rec tracks */
  }
  if (fMCRecFlags == NULL) {
    fMCRecFlags = new UInt_t[fMCFlagsStorageSize];
  }
  if (fMCTruePrimaryFlags == NULL) {
    fMCTruePrimaryFlags = new UInt_t[fMCFlagsStorageSize];
  }

  /* expand the storage if needed */
  Int_t nNoOfTrueTracks = GetNoOfTrueParticles();
  Int_t nNoOfTruePrimaries = GetNoOfTruePrimaries();
  Int_t nNoOfMCRecTracks = GetNoOfMCRecTracks();

  if (fTrueToRec->GetSize() < nNoOfTrueTracks)
    fTrueToRec->Expand((Int_t(nNoOfTrueTracks/1024)+5)*1024);
  if (fTrueToRecWrong->GetSize() < nNoOfTrueTracks)
    fTrueToRecWrong->Expand((Int_t(nNoOfTrueTracks/1024)+5)*1024);
  if (fMCFlagsStorageSize < TMath::Max(nNoOfTruePrimaries,nNoOfMCRecTracks)) {
    fMCFlagsStorageSize = (Int_t(TMath::Max(nNoOfTruePrimaries,nNoOfMCRecTracks)/1024)+5)*1024;
    delete [] fMCRecFlags;
    delete [] fMCTruePrimaryFlags;
    fMCRecFlags = new UInt_t[fMCFlagsStorageSize];
    fMCTruePrimaryFlags = new UInt_t[fMCFlagsStorageSize];
  }

  /* clear previous events information */
  fTrueToRec->Clear();
  fTrueToRecWrong->Clear();
  memset(fMCRecFlags,0,fMCFlagsStorageSize*sizeof(UInt_t));
  memset(fMCTruePrimaryFlags,kNone,fMCFlagsStorageSize*sizeof(UInt_t));

  AliInfo(Form("Rec flags size: %d, True to rec size: %d, True to rec wrong size: %d",
      fMCFlagsStorageSize, fTrueToRec->GetSize(), fTrueToRecWrong->GetSize()));

  AliVEvent *event = InputEvent();

  for (Int_t iaod = 0; iaod < nNoOfMCRecTracks; iaod++) {
    /* just to be sure we do not surpass the limits */
    if (!(iaod < event->GetNumberOfTracks())) break;

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(event->GetTrack(iaod)); // pointer to reconstructed track
    if(!vtrack) {
      fMCRecFlags[iaod] |= kNoTrack;
      continue; /* error is reported in normal rec track processing */
    }
    /* so far let's handle the not constrained tracks */
    if (vtrack->GetID() < 0) continue;

    Int_t iesd = vtrack->GetID();

    if (vtrack->GetLabel() < 0) {

      fMCRecFlags[iesd] |= kNegativeLabel;

      if (fTrackSelectionCuts->IsFromMCInjectedSignal(TMath::Abs(vtrack->GetLabel()))) {
        fMCRecFlags[iesd] |= kFromInjectedSignal;
        continue;
      }

      if (fTrueToRecWrong->At(TMath::Abs(vtrack->GetLabel())) != NULL) {
        AliVTrack *prevvtrack = (AliVTrack *) fTrueToRecWrong->At(TMath::Abs(vtrack->GetLabel()));

        fMCRecFlags[iesd] |= kSameNegativeLabel;
        fMCRecFlags[prevvtrack->GetID()] |= kSameNegativeLabel;
      }
      else
        fTrueToRecWrong->AddAt(vtrack, TMath::Abs(vtrack->GetLabel()));
    }
    else {
      if (fTrackSelectionCuts->IsFromMCInjectedSignal(vtrack->GetLabel())) {
        fMCRecFlags[iesd] |= kFromInjectedSignal;
        continue;
      }

      if (fTrueToRec->At(vtrack->GetLabel()) != NULL) {
        /* we don't build multi reconstruction information yet */
      }
      else
        fTrueToRec->AddAt(vtrack, vtrack->GetLabel());
    }
  }
}

/// \brief Re-builds the true to reconstructed relation based on the
/// accepted reconstructed tracks
/// The most relevant information incorporated is the proper mapping
/// of the true to rec when few rec point to the same true but they
/// are not all accepted
void AliAnalysisTaskCorrelationsStudies::BuildTrueRecAccRelation() {

  AliVEvent *event = InputEvent();
  for (Int_t iaod = 0; iaod < GetNoOfMCRecTracks(); iaod++) {
    /* just to be sure we do not surpass the limits */
    if (!(iaod < event->GetNumberOfTracks())) break;

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(event->GetTrack(iaod)); // pointer to reconstructed track
    if(!vtrack) {
      continue; /* error is reported in normal rec track processing */
    }

    if (vtrack->GetLabel() < 0) {
      Int_t iesd = vtrack->GetID();
      if (iesd < 0) iesd = -1 - iesd;

      /* track should not have been accepted */
      if ((fMCRecFlags[iesd] & kAccepted) == kAccepted) {
        AliInfo(Form("ERROR? - track %d with negative label %d has been accepted. AOD Constrained?",iesd,vtrack->GetLabel()));
      }
    }
    else {
      Int_t iesd = vtrack->GetID();
      if (iesd < 0) iesd = -1 - iesd;

      if (fTrueToRec->At(vtrack->GetLabel()) != NULL) {
        AliVTrack *prevvtrack = (AliVTrack *) fTrueToRec->At(vtrack->GetLabel());
        if ((fMCRecFlags[iesd] & kAccepted) == kAccepted) {
          /* track has been accepted, store multi reconstruction information */
          if (prevvtrack->GetID() != iesd) {
            /* different tracks so update the information accordingly */
            if ((fMCRecFlags[prevvtrack->GetID()] & kAccepted) == kAccepted) {
              /* both tracks are accepted with the same label */
              fMCRecFlags[iesd] |= kSamePositiveLabel;
              fMCRecFlags[prevvtrack->GetID()] |= kSamePositiveLabel;
              AliError(Form("ERROR - track %d and %d accepted with the same label %d",
                  iesd, prevvtrack->GetID(), vtrack->GetLabel()));
            }
            else {
              /* the previous track was not accepted, set this one instead */
              fTrueToRec->AddAt(vtrack, vtrack->GetLabel());
            }
          }
          else {
            /* the same track so, no multi reconstruction so far */
            fMCRecFlags[iesd] &= ~kSamePositiveLabel;
          }
        }
        else {
          /* track has not been accepted, nothing to do */
        }
      }
      else {
        /* not mapped, incorporate it to the map */
        fTrueToRec->AddAt(vtrack, vtrack->GetLabel());
      }
    }
  }
}

/// \brief Sets the task configuration by means of a configuration string
/// \param confstring the configuration string
/// \return kTRUE if the task was correctly configured kFALSE otherwise
/// The configuration string will have the format
/// Task:activities;Events:eventstring;Tracks:trackstring;PID:pidstring;
/// In their turn, eventstring, trackstring and pidstring might have the format
/// cutstring or cutstring+cutstring, and cutstring could start with a '-'
/// to signal an exlusive cut.
/// The activities could be none (empty), any of 'corr' or 'pairs' or both
/// as 'corr+pairs' where 'corr' stands for correlation analysis and 'pairs'
/// for track pairs analysis
///
/// Additionally, between the event and tracks configuration strings there are
/// a set of options which are only effective for MC data analysis
///
/// - Effprof:profile;  the efficiency profile to enforce on MC rec data
/// - Trueeffprof:profile;  the efficiency profile to enforce on MC true date
/// - Recoption:option;   options for the additional MC rec results
///
/// The profile string could be a single expression depending on x (pT),
/// or a concatenation |bin,bin,expr| terms each of them defining an initial
/// bin, a final bin and an expression to assign the efficiency on that bin
/// interval. In this case the expression will depend on pT bin number. Again
/// the expression  will not depend on the pT value of the concrete bin being
/// handled. If pT value were needed use the corresponding factor within the
/// expression.
///
/// The option string could be one of the following
/// - trueval   the reconstructed results will use reconstructed tracks but with true data
/// - primaries   the reconstructed results will use true primary reconstructed tracks
/// It should only be used at initial task configuration
Bool_t AliAnalysisTaskCorrelationsStudies::Configure(const char *confstring)
{
  Bool_t accepted = kTRUE;
  fTaskConfigurationString = confstring;
  TString sztmp = confstring;

  /* the task activities */
  if (sztmp.BeginsWith("Task:")) {
    sztmp.Remove(0,strlen("Task:"));
    fTaskActivitiesString = sztmp(0,sztmp.Index(";"));
    sztmp.Remove(0,sztmp.Index(";")+1);

    TObjArray *tokens = fTaskActivitiesString.Tokenize("+");
    if (tokens->GetEntries() == 2) {
      if (((TObjString*)tokens->At(0))->String().EqualTo("corr")) {
        if (((TObjString*)tokens->At(1))->String().EqualTo("pairs")) {
          fDoProcessCorrelations = kTRUE;
          fDoProcessPairAnalysis = kTRUE;
        }
        else {
          AliFatal("Wrong task activities. ABORTING!!!");
          return kFALSE;
        }
      }
      else if (((TObjString*)tokens->At(1))->String().EqualTo("corr")) {
        if (((TObjString*)tokens->At(0))->String().EqualTo("pairs")) {
          fDoProcessCorrelations = kTRUE;
          fDoProcessPairAnalysis = kTRUE;
        }
        else {
          AliFatal("Wrong task activities. ABORTING!!!");
          return kFALSE;
        }
      }
      else {
        AliFatal("Wrong task activities. ABORTING!!!");
        return kFALSE;
      }
      delete tokens;
    }
    else if (tokens->GetEntries() == 1) {
      if (((TObjString*)tokens->At(0))->String().EqualTo("corr")) {
        fDoProcessCorrelations = kTRUE;
      }
      else if (((TObjString*)tokens->At(0))->String().EqualTo("pairs")) {
          fDoProcessPairAnalysis = kTRUE;
      }
      else {
        AliFatal("Wrong task activities. ABORTING!!!");
        return kFALSE;
      }
    }
    else if (tokens->GetEntries() != 0) {
      AliFatal("Wrong task activities. ABORTING!!!");
      return kFALSE;
    }
  }

  /* the event cuts */
  if (sztmp.BeginsWith("Events:")) {
    sztmp.Remove(0,strlen("Events:"));
    fEventCutsString = sztmp(0,sztmp.Index(";"));
    sztmp.Remove(0,sztmp.Index(";")+1);
  }
  else {
    AliFatal("No event cuts string. ABORTING!!!");
    return kFALSE;
  }

  /* the efficiency profile to enforce */
  if (sztmp.BeginsWith("Effprof:")) {
    sztmp.Remove(0,strlen("Effprof:"));
    fEfficiencyProfileToEnforce = sztmp(0, sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);
  }

  /* the efficiency profile on true data */
  if (sztmp.BeginsWith("Trueeffprof:")) {
    sztmp.Remove(0,strlen("Trueeffprof:"));
    fEfficiencyProfileOnTrue = sztmp(0, sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);
    if (fEfficiencyProfileOnTrue.EqualTo("onlyrec")) {
      fEfficiencyProfileOnTrue.Clear();
      fOnTrueEfficiencyProfile = kOnlyReconstructed;
    }
  }

  /* options for the additional MC rec results */
  if (sztmp.BeginsWith("Recoption:")) {
    sztmp.Remove(0,strlen("Recoption:"));
    fAdditionalMCRecOption = sztmp(0, sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);

    if (fAdditionalMCRecOption.EqualTo("trueval")) {
      fMCRecOption = kRecWithTrue;
    }
    else if(fAdditionalMCRecOption.EqualTo("primaries")) {
      fMCRecOption = kRecTruePrimaries;
    }
    else if(fAdditionalMCRecOption.EqualTo("notacc")) {
      fMCRecOption = kRecWithNotAccepted;
    }
    else {
      AliFatal("Unknown option for additional MC rec results. ABORTING!!!");
      return kFALSE;
    }
  }

  /* the track selection strings                           */
  fTrackSelectionString = sztmp;
  /* let's make a quick survey to detect inconsistent format */
  /* the tracks section */
  if (sztmp.BeginsWith("Tracks:")) {
    sztmp.Remove(0,strlen("Tracks:"));
    TString sztmptrk = sztmp(0,sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);
  }
  else {
    AliFatal("No tracks selection cuts string. ABORTING!!!");
    return kFALSE;
  }
  /* the PID section */
  if (sztmp.BeginsWith("PID:")) {
    sztmp.Remove(0,strlen("PID:"));
    TString sztmpid = sztmp(0, sztmp.Index(";"));
    sztmp.Remove(0, sztmp.Index(";")+1);
  }
  else {
    AliFatal("No tracks selection PID cuts string. ABORTING!!!");
    return kFALSE;
  }

  if (!fTaskConfigurationString.EqualTo(ProduceConfigurationString())) {
    AliError("Produced configuration ");
    AliError(Form("\t%s",ProduceConfigurationString()));
    AliError("does not match the passed one");
    AliError(Form("\t%s",fTaskConfigurationString.Data()));
    AliFatal("ABORTING!!!");
    return kFALSE;
  }

  AliInfo("Task configured with configuration string:");
  AliInfo(Form("\t%s", fTaskConfigurationString.Data()));
  return accepted;
}

/// \brief Builds a configuration string out of the task variables
/// \return the built configuration string
const char *AliAnalysisTaskCorrelationsStudies::ProduceConfigurationString() {
  static char buffer[2048];

  sprintf(buffer, "Task:%s;Events:%s;%s%s%s%s",
      fTaskActivitiesString.Data(),
      fEventCutsString.Data(),
      ((fEfficiencyProfileToEnforce.Length() != 0) ? Form("Effprof:%s;",fEfficiencyProfileToEnforce.Data()) : ""),
      ((fEfficiencyProfileOnTrue.Length() != 0) ? Form("Trueeffprof:%s;",fEfficiencyProfileOnTrue.Data()) :
          ((fOnTrueEfficiencyProfile == kOnlyReconstructed) ? "Trueeffprof:onlyrec;" : "")),
      ((fAdditionalMCRecOption.Length() != 0) ? Form("Recoption:%s;",fAdditionalMCRecOption.Data()) : ""),
      fTrackSelectionString.Data());

  return buffer;
}

///\brief Establishes the processing correlations instance configuration
/// \param confstring the configuration string
/// \param pattern string pattern for the weigths or particle density files
/// \param szContainerPrefix for placing the prefix of the particle and action combination
/// \return kTRUE if everything went OK kFALSE otherwise
Bool_t AliAnalysisTaskCorrelationsStudies::ConfigureCorrelations(const char *confstring, const char *pattern, TString &szContainerPrefix)
{
  TString sztmp = confstring;
  TString szTrack1;
  TString szTrack2;

  /* the correlation configuration */
  if (sztmp.BeginsWith("Corr:")) {
    sztmp.Remove(0,strlen("Corr:"));
  }
  else {
    AliFatal("No correlations configuration string. ABORTING!!!");
    return kFALSE;
  }

  TObjArray *tokens = sztmp.Tokenize(",");
  if ((tokens->GetEntries() == 4) || (tokens->GetEntries() == 5)) {
    /* track polarities */
    if (((TObjString*) tokens->At(0))->String().EqualTo("--")) {
      fProcessCorrelations.SetSameSign(kTRUE);
      fProcessCorrelations.SetRequestedCharge_1(-1);
      fProcessCorrelations.SetRequestedCharge_2(-1);
      fProcessMCRecCorrelationsWithOptions.SetSameSign(kTRUE);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_1(-1);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_2(-1);
      fProcessTrueCorrelations.SetSameSign(kTRUE);
      fProcessTrueCorrelations.SetRequestedCharge_1(-1);
      fProcessTrueCorrelations.SetRequestedCharge_2(-1);
      szTrack1 = "m1";
      szTrack2 = "m2";
      szContainerPrefix = "MM";
    }
    else if (((TObjString*) tokens->At(0))->String().EqualTo("++")) {
      fProcessCorrelations.SetSameSign(kTRUE);
      fProcessCorrelations.SetRequestedCharge_1(1);
      fProcessCorrelations.SetRequestedCharge_2(1);
      fProcessMCRecCorrelationsWithOptions.SetSameSign(kTRUE);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_1(1);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_2(1);
      fProcessTrueCorrelations.SetSameSign(kTRUE);
      fProcessTrueCorrelations.SetRequestedCharge_1(1);
      fProcessTrueCorrelations.SetRequestedCharge_2(1);
      szTrack1 = "p1";
      szTrack2 = "p2";
      szContainerPrefix = "PP";
    }
    else if (((TObjString*) tokens->At(0))->String().EqualTo("+-")) {
      fProcessCorrelations.SetSameSign(kFALSE);
      fProcessCorrelations.SetRequestedCharge_1(1);
      fProcessCorrelations.SetRequestedCharge_2(-1);
      fProcessMCRecCorrelationsWithOptions.SetSameSign(kFALSE);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_1(1);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_2(-1);
      fProcessTrueCorrelations.SetSameSign(kFALSE);
      fProcessTrueCorrelations.SetRequestedCharge_1(1);
      fProcessTrueCorrelations.SetRequestedCharge_2(-1);
      szTrack1 = "p1";
      szTrack2 = "m2";
      szContainerPrefix = "PM";
    }
    else if (((TObjString*) tokens->At(0))->String().EqualTo("-+")) {
      fProcessCorrelations.SetSameSign(kFALSE);
      fProcessCorrelations.SetRequestedCharge_1(-1);
      fProcessCorrelations.SetRequestedCharge_2(1);
      fProcessMCRecCorrelationsWithOptions.SetSameSign(kFALSE);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_1(-1);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_2(1);
      fProcessTrueCorrelations.SetSameSign(kFALSE);
      fProcessTrueCorrelations.SetRequestedCharge_1(-1);
      fProcessTrueCorrelations.SetRequestedCharge_2(1);
      szTrack1 = "m1";
      szTrack2 = "p2";
      szContainerPrefix = "MP";
    }
    else if (((TObjString*) tokens->At(0))->String().EqualTo("**")) {
      fProcessCorrelations.SetSameSign(kFALSE);
      fProcessCorrelations.SetAllCombinations(kTRUE);
      fProcessCorrelations.SetRequestedCharge_1(1);
      fProcessCorrelations.SetRequestedCharge_2(-1);
      fProcessMCRecCorrelationsWithOptions.SetSameSign(kFALSE);
      fProcessMCRecCorrelationsWithOptions.SetAllCombinations(kTRUE);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_1(1);
      fProcessMCRecCorrelationsWithOptions.SetRequestedCharge_2(-1);
      fProcessTrueCorrelations.SetSameSign(kFALSE);
      fProcessTrueCorrelations.SetAllCombinations(kTRUE);
      fProcessTrueCorrelations.SetRequestedCharge_1(1);
      fProcessTrueCorrelations.SetRequestedCharge_2(-1);
      szTrack1 = "p1";
      szTrack2 = "m2";
      szContainerPrefix = "AA";
    }
    else {
      AliFatal("Requested track polarities string not properly configured.ABORTING!!!");
      return kFALSE;
    }

    /* singles or pairs */
    if (((TObjString*) tokens->At(1))->String().EqualTo("singles")) {
      fProcessCorrelations.SetSinglesOnly(kTRUE);
      fProcessMCRecCorrelationsWithOptions.SetSinglesOnly(kTRUE);
      fProcessTrueCorrelations.SetSinglesOnly(kTRUE);
      szContainerPrefix += "_S";
    }
    else if (((TObjString*) tokens->At(1))->String().EqualTo("pairs")) {
      fProcessCorrelations.SetSinglesOnly(kFALSE);
      fProcessMCRecCorrelationsWithOptions.SetSinglesOnly(kFALSE);
      fProcessTrueCorrelations.SetSinglesOnly(kFALSE);
      szContainerPrefix += "_P";
    }
    else {
      AliFatal("Singles / pairs string not properly configured.ABORTING!!!");
      return kFALSE;
    }

    /* resonace rejection configuration */
    if (((TObjString*) tokens->At(2))->String().Contains("resonances:")) {
      fProcessCorrelations.ConfigureResonances(((TObjString*) tokens->At(2))->String().Data());
    }
    else {
      AliFatal("Resonance string not properly configured.ABORTING!!!");
      return kFALSE;
    }

    /* use weights or efficiency corrections or both */
    if (((TObjString*) tokens->At(3))->String().EqualTo("weights") ||
        ((TObjString*) tokens->At(3))->String().EqualTo("effcorr") ||
        ((TObjString*) tokens->At(3))->String().EqualTo("weightseffcorr") ||
        ((TObjString*) tokens->At(3))->String().EqualTo("weightspairseff") ) {
      /* we will use weights or efficiency corrections so the weights filename has to be there */
      if (tokens->GetEntries() == 5) {
        /* get the weights histos */
        TFile *weights;
        Char_t localbuffer[2048];

        TGrid::Connect("alien:");
        weights = TFile::Open(((TObjString*) tokens->At(4))->String(),"OLD");
        if (weights != NULL && weights->IsOpen()) {
          /* the weights correction */
          if (((TObjString*) tokens->At(3))->String().EqualTo("weights") ||
              ((TObjString*) tokens->At(3))->String().EqualTo("weightseffcorr") ||
              ((TObjString*) tokens->At(3))->String().EqualTo("weightspairseff")) {
            sprintf(localbuffer,"%s%s", pattern,szTrack1.Data());
            fhWeightsTrack_1 = (TH3F*) weights->Get(localbuffer);
            sprintf(localbuffer,"%s%s", pattern,szTrack2.Data());
            fhWeightsTrack_2 = (TH3F*) weights->Get(localbuffer);
            if ((fhWeightsTrack_1 != NULL) && (fhWeightsTrack_2 != NULL)) {
              fProcessCorrelations.SetUseWeights(kTRUE);
              AliInfo("===========STORED CORRECTION WEIGHTS====================");
              AliInfo(Form("Track 1: %c charge, weights histogram: %s", ((TObjString*) tokens->At(0))->String()[0],fhWeightsTrack_1->GetName()));
              AliInfo(Form("Track 2: %c charge, weights histogram: %s", ((TObjString*) tokens->At(0))->String()[1],fhWeightsTrack_2->GetName()));
              AliInfo("===========END STORED CORRECTION WEIGHTS================");
              szContainerPrefix += "W";
            }
            else {
              AliFatal(Form("Not able to find weights histograms %s in file %s. ABORTING!!!",localbuffer, ((TObjString*) tokens->At(3))->String().Data()));
              return kFALSE;
            }
          }
          if (((TObjString*) tokens->At(3))->String().EqualTo("effcorr") ||
              ((TObjString*) tokens->At(3))->String().EqualTo("weightseffcorr") ) {
            /* the efficiency correction */
            sprintf(localbuffer,"%seff_%s", pattern,szTrack1.Data());
            fhEffCorrTrack_1 = (TH1F*) weights->Get(localbuffer);
            sprintf(localbuffer,"%seff_%s", pattern,szTrack2.Data());
            fhEffCorrTrack_2 = (TH1F*) weights->Get(localbuffer);
            if ((fhEffCorrTrack_1 != NULL) && (fhEffCorrTrack_2 != NULL)) {
              AliInfo("===========STORED EFFICIENCY CORRECTION====================");
              AliInfo(Form("Track 1: %c charge, efficiency correction histogram: %s", ((TObjString*) tokens->At(0))->String()[0],fhEffCorrTrack_1->GetName()));
              AliInfo(Form("Track 2: %c charge, efficiency correction histogram: %s", ((TObjString*) tokens->At(0))->String()[1],fhEffCorrTrack_2->GetName()));
              AliInfo("===========END STORED EFFICIENCY CORRECTION================");
              szContainerPrefix += "E";
            }
            else {
              AliFatal(Form("Not able to find efficiency correction histograms %s in file %s. ABORTING!!!", localbuffer, ((TObjString*) tokens->At(3))->String().Data()));
              return kFALSE;
            }
          }
          if (((TObjString*) tokens->At(3))->String().EqualTo("pairseff") ||
              ((TObjString*) tokens->At(3))->String().EqualTo("weightspairseff") ) {
            /* the pairs efficiencies */
            sprintf(localbuffer,"%spairseff_PP", pattern);
            fhPairEfficiency_PP = (THn*) weights->Get(localbuffer);
            sprintf(localbuffer,"%spairseff_PM", pattern);
            fhPairEfficiency_PM = (THn*) weights->Get(localbuffer);
            sprintf(localbuffer,"%spairseff_MM", pattern);
            fhPairEfficiency_MM = (THn*) weights->Get(localbuffer);
            sprintf(localbuffer,"%spairseff_MP", pattern);
            fhPairEfficiency_MP = (THn*) weights->Get(localbuffer);
            if ((fhPairEfficiency_PP != NULL) && (fhPairEfficiency_PM != NULL) && (fhPairEfficiency_MM != NULL) && (fhPairEfficiency_MP != NULL)) {
              AliInfo("===========STORED PAIRS EFFICIENCY=========================");
              AliInfo(Form("Pair PP: pairs efficiency histogram: %s", fhPairEfficiency_PP->GetName()));
              AliInfo(Form("Pair PM: pairs efficiency histogram: %s", fhPairEfficiency_PM->GetName()));
              AliInfo(Form("Pair MM: pairs efficiency histogram: %s", fhPairEfficiency_MM->GetName()));
              AliInfo(Form("Pair MP: pairs efficiency histogram: %s", fhPairEfficiency_MP->GetName()));
              AliInfo("===========END STORED PAIRS EFFICIENCY=====================");
              szContainerPrefix += "PE";
            }
            else {
              AliFatal(Form("Not able to find pairs efficiency histograms %s in file %s. ABORTING!!!", localbuffer, ((TObjString*) tokens->At(3))->String().Data()));
              return kFALSE;
            }
          }
        }
        else {
          AliFatal(Form("Weights file %s not available.ABORTING!!!",((TObjString*) tokens->At(3))->String().Data()));
          return kFALSE;
        }
      }
      else {
        AliFatal("Use weights required but the weights filename is not present.ABORTING!!!");
        return kFALSE;
      }
    }
    else if (((TObjString*) tokens->At(3))->String().BeginsWith("simulate")) {
      /* get the number of events to generate per real event */
      Int_t nSimEventsPerEvent;
      if (((TObjString*) tokens->At(3))->String().EqualTo("simulate"))
        nSimEventsPerEvent = 1;
      else
        sscanf(((TObjString*) tokens->At(3))->String().Data(), "simulate-%d", &nSimEventsPerEvent);

      fProcessCorrelations.SetUseSimulation(kTRUE);
      fProcessCorrelations.SetSimEventsPerEvent(nSimEventsPerEvent);
      /* we will use particle profiles so the particle profiles filename has to be there */
      if (tokens->GetEntries() == 5) {
        /* get the weights histos */
        TFile *trkprofiles;
        Char_t localbuffer[2048];

        TGrid::Connect("alien:");
        trkprofiles = TFile::Open(((TObjString*) tokens->At(4))->String(),"OLD");
        if (trkprofiles != NULL && trkprofiles->IsOpen()) {
          /* get the positive tracks density function */
          fPositiveTrackPdf = new TObjArray(64); fPositiveTrackPdf->SetOwner(kTRUE);
          for (Int_t ix = 0; ;ix++) {
            sprintf(localbuffer, "%sp_yield_%d", pattern, ix);
            TH3F *hslot = (TH3F*) trkprofiles->Get(localbuffer);
            if (hslot != NULL)
              fPositiveTrackPdf->Add(hslot);
            else
              break;
          }
          printf("Took %d positive tracks slots with %s pattern\n", fPositiveTrackPdf->GetEntriesFast(), Form("%sp_yield", pattern));
          /* and now the negative one */
          fNegativeTrackPdf = new TObjArray(64); fNegativeTrackPdf->SetOwner(kTRUE);
          for (Int_t ix = 0; ;ix++) {
            sprintf(localbuffer, "%sm_yield_%d", pattern, ix);
            TH3F *hslot = (TH3F*) trkprofiles->Get(localbuffer);
            if (hslot != NULL)
              fNegativeTrackPdf->Add(hslot);
            else
              break;
          }
          printf("Took %d negative tracks slots with %s pattern\n", fNegativeTrackPdf->GetEntriesFast(), Form("%sm_yield", pattern));
          if ((fPositiveTrackPdf->GetEntriesFast() != 0) &&
              (fNegativeTrackPdf->GetEntriesFast() != 0) &&
              (fPositiveTrackPdf->GetEntriesFast() == fNegativeTrackPdf->GetEntriesFast())) {
            AliInfo("=========== STORED SIMULATION TRACKS DENSITY PROFILES =================");
            AliInfo(Form("Positive tracks: %d zvertex slots", fPositiveTrackPdf->GetEntriesFast()));
            AliInfo(Form("Negative tracks: %d zvertex slots", fNegativeTrackPdf->GetEntriesFast()));
            AliInfo("============= END SIMULATION TRACKS DENSITY PROFILES ==================");
            szContainerPrefix += "SIM";
          }
          else {
            AliFatal(Form("Not able to find consistent track density histograms in file %s. ABORTING!!!",((TObjString*) tokens->At(3))->String().Data()));
            return kFALSE;
          }
        }
        else {
          AliFatal(Form("Particle profiles file %s not available.ABORTING!!!",((TObjString*) tokens->At(3))->String().Data()));
          return kFALSE;
        }
      }
      else {
        AliFatal("Use simulation required but the particle profiles filename is not present.ABORTING!!!");
        return kFALSE;
      }
    }
    else if (((TObjString*) tokens->At(3))->String().EqualTo("effcorrtrue")) {
      fCorrectOnTrueEfficiency = kTRUE;
      fProcessCorrelations.SetUseWeights(kFALSE);
      szContainerPrefix += "NW";
    }
    else if (((TObjString*) tokens->At(3))->String().EqualTo("noweights")) {
      fProcessCorrelations.SetUseWeights(kFALSE);
      szContainerPrefix += "NW";
    }
    else {
      AliFatal("(no)weights string not properly configured.ABORTING!!!");
      return kFALSE;
    }

  }
  else {
    AliFatal("Correlations configuration string not properly configured.ABORTING!!!");
    return kFALSE;
  }
  delete tokens;

  /* true will never use weights nor simulation */
  fProcessMCRecCorrelationsWithOptions.SetUseWeights(kFALSE);
  fProcessMCRecCorrelationsWithOptions.SetUseSimulation(kFALSE);
  fProcessTrueCorrelations.SetUseWeights(kFALSE);
  fProcessTrueCorrelations.SetUseSimulation(kFALSE);


  return kTRUE;
}

/// \brief Establishes the bining for the correlation object instance
/// \param confstring the string containing the bins configuration
/// \return kTRUE if everything went OK kFALSE otherwise
Bool_t AliAnalysisTaskCorrelationsStudies::ConfigureCorrelationsBinning(const char *confstring) {

  TString sztmp = confstring;
  TString szTrack1;
  TString szTrack2;

  /* the correlation configuration */
  if (sztmp.BeginsWith("Binning:")) {
    sztmp.Remove(0,strlen("Binning:"));
  }
  else {
    AliFatal("No correlations binning configuration string. ABORTING!!!");
    return kFALSE;
  }

  return (fProcessCorrelations.ConfigureBinning(sztmp.Data()) &&
          fProcessMCRecCorrelationsWithOptions.ConfigureBinning(sztmp.Data()) &&
          fProcessTrueCorrelations.ConfigureBinning(sztmp.Data()) &&
          fProcessPairAnalysis.ConfigureBinning(sztmp.Data()) &&
          fProcessTruePairAnalysis.ConfigureBinning(sztmp.Data()));
}


/// Creates the object to be preserved on the results files
void AliAnalysisTaskCorrelationsStudies::UserCreateOutputObjects()
{
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // Create histograms
  // Called once (on the worker node)
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!

  fEventCuts = new AliCSEventCuts("CS_Event_Cuts", "CS Event Cuts");
  if (!fEventCuts->InitializeCutsFromCutString(fEventCutsString)) {
    AliFatal("Not able to initialize event cuts. ABORTING!!!");
  }

  fEventCuts->SetQALevelOutput(AliCSAnalysisCutsBase::kQALevelHeavy);
  fEventCuts->InitCuts();
  fOutput->Add(fEventCuts->GetHistogramsList());

  fTrackSelectionCuts = new AliCSTrackSelection("CS_Track_Selection", "CS Track Selection");
  fTrackSelectionCuts->InitializeFromString(fTrackSelectionString);
  fTrackSelectionCuts->SetQALevelOutput(AliCSAnalysisCutsBase::kQALevelHeavy);
  fTrackSelectionCuts->InitCuts();
  fOutput->Add(fTrackSelectionCuts->GetHistogramsList());

  /* now initialize the processing correlations instance */
  if (fDoProcessCorrelations) {
    fProcessCorrelations.Initialize();
    fProcessMCRecCorrelationsWithOptions.Initialize();
    fProcessTrueCorrelations.Initialize();
    fProcessCorrelations.SetWeigths(fhWeightsTrack_1, fhWeightsTrack_2);
    fProcessCorrelations.SetEfficiencyCorrection(fhEffCorrTrack_1, fhEffCorrTrack_2);
    fProcessCorrelations.SetPairEfficiencyCorrection(fhPairEfficiency_PP, fhPairEfficiency_PM, fhPairEfficiency_MM, fhPairEfficiency_MP);
    fProcessCorrelations.SetSimultationPdfs(fPositiveTrackPdf, fNegativeTrackPdf);
  }

  /* now initialize the pair analysis instance */
  if (fDoProcessPairAnalysis) {
    fProcessPairAnalysis.Initialize();
    fProcessTruePairAnalysis.Initialize();
    /* we use the weights, they should contain NUA x NUE in four dimensions vtxz, eta, phi and pT */
    fProcessPairAnalysis.SetSinglesEfficiency(fhWeightsTrack_1, fhWeightsTrack_2);
    fProcessPairAnalysis.SetPairEfficiency(fhPairEfficiency_PP, fhPairEfficiency_PM, fhPairEfficiency_MM, fhPairEfficiency_MP);
    /* and incorporate its histograms to the output list */
    /* now initialize the pair analysis instance */
    fOutput->Add(fProcessPairAnalysis.GetHistogramsList());
    fOutput->Add(fProcessTruePairAnalysis.GetHistogramsList());
  }

  // Create histograms
  Int_t ptbins = 30;
  Float_t ptlow = 0.1, ptup = 3.1;
  fhPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
  fhPt->SetMarkerStyle(kFullCircle);
  fhPtPos = new TH1F("fHistPtPos", "P_{T} distribution for reconstructed (#{+});P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
  fhPtNeg = new TH1F("fHistPtNeg", "P_{T} distribution for reconstructed (#{-});P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
  fhP = new TH1F("fHistP", "P distribution for reconstructed;P (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
  fhPPos = new TH1F("fHistPPos", "P distribution for reconstructed (#{+});P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);
  fhPNeg = new TH1F("fHistPNeg", "P distribution for reconstructed (#{-});P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);
  fhUnConstrainedPt = new TH1F("fHistUnConstrainedPt","Unconstrained P_{T} distribution for reconstructed;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
  fhPtDifference = new TH1F("fHistDifferencePt","Difference P_{T} distribution for reconstructed;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, -0.05, 0.05);

  Int_t etabins = 40;
  Float_t etalow = -2.0, etaup = 2.0;
  Int_t zvtxbins = 40;
  Float_t zvtxlow = -10.0, zvtxup = 10.0;
  Int_t phibins = 72;

  fhEtaB = new TH1F("fHistEtaB","#eta distribution for reconstructed before;#eta;counts",etabins, etalow, etaup);
  fhEtaA = new TH1F("fHistEtaA","#eta distribution for reconstructed;#eta;counts",etabins, etalow, etaup);

  fhPhiB = new TH1F("fHistPhiB","#phi distribution for reconstructed before;#phi;counts",360, 0.0, 360.0);
  fhPhiA = new TH1F("fHistPhiA","#phi distribution for reconstructed;#phi;counts",360, 0.0, 360.0);

  fhEtaVsPhiB = new TH2F(Form("CSTaskEtaVsPhiB_%s",fTaskConfigurationString.Data()),"#eta vs #phi before;#phi;#eta", 360, 0.0, 360.0, 100, etalow, etaup);
  fhEtaVsPhiA = new TH2F(Form("CSTaskEtaVsPhiA_%s",fTaskConfigurationString.Data()),"#eta vs #phi;#phi;#eta", 360, 0.0, 360.0, 100, -2.0, 2.0);

  fhPtVsEtaB = new TH2F(Form("fhPtVsEtaB_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta before;#eta;p_{T} (GeV/c)",etabins,etalow,etaup,200,0.0,10.0);
  fhPtVsEtaA = new TH2F(Form("fhPtVsEtaA_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta;#eta;p_{T} (GeV/c)",etabins,etalow,etaup,200,0.0,10.0);

  fhPt3DB = new TProfile3D(Form("fhPt3DB_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta, #phi and vtx_{z} before;#eta;#phi;vtx_{z}",etabins,etalow,etaup,phibins,0.0,2*TMath::Pi(),zvtxbins,zvtxlow,zvtxup);
  fhPt3DA = new TProfile3D(Form("fhPt3DA_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta, #phi and vtx_{z};#eta;#phi;vtx_{z}",etabins,etalow,etaup,phibins,0.0,2*TMath::Pi(),zvtxbins,zvtxlow,zvtxup);

  fh3Dn1B = new TH3I(Form("fh3Dn1B_%s",fTaskConfigurationString.Data()),"n_{1} vs #eta, #phi and p_{T} before;#eta;#varphi;p_{T}",etabins,etalow,etaup,phibins,0.0,2*TMath::Pi(),ptbins,ptlow,ptup);
  fh3Dn1A = new TH3I(Form("fh3Dn1A_%s",fTaskConfigurationString.Data()),"n_{1} vs #eta, #phi and p_{T};#eta;#varphi;p_{T}",etabins,etalow,etaup,phibins,0.0,2*TMath::Pi(),ptbins,ptlow,ptup);

  /* the tracking of the event multiplicity distribution for evaluating the lambda angle */
  Int_t mBins[2] = {4000,4000};
  Double_t mLow[2] = {0.0,0.0};
  Double_t mUp[2] = {4000.0,4000.0};
  fhAcceptedVsMultiplicity = new THnSparseI("AcceptedVsMultiplicity",
      "Accepted vs Reconstructed tracks;reconstructed tracks,accepted tracks,counts",
      2,mBins,mLow,mUp);
  fhLambdaVsMultiplicity = new TH2I("fhLambdaVsMultipliciy","#lambda;reconstructed tracks;#lambda;counts",
      4000,0,4000,72,-TMath::PiOver2(),TMath::PiOver2());


  fOutput->Add(fhPt);
  fOutput->Add(fhPtPos);
  fOutput->Add(fhPtNeg);
  fOutput->Add(fhP);
  fOutput->Add(fhPPos);
  fOutput->Add(fhPNeg);
  fOutput->Add(fhUnConstrainedPt);
  fOutput->Add(fhPtDifference);
  fOutput->Add(fhEtaB);
  fOutput->Add(fhEtaA);
  fOutput->Add(fhPhiB);
  fOutput->Add(fhPhiA);
  fOutput->Add(fhEtaVsPhiB);
  fOutput->Add(fhEtaVsPhiA);
  fOutput->Add(fhPtVsEtaB);
  fOutput->Add(fhPtVsEtaA);
  fOutput->Add(fhPt3DB);
  fOutput->Add(fhPt3DA);
  fOutput->Add(fh3Dn1B);
  fOutput->Add(fh3Dn1A);
  fOutput->Add(fhAcceptedVsMultiplicity);
  fOutput->Add(fhLambdaVsMultiplicity);

  /* incorporate configuration parameters if needed */
  if (fDoProcessCorrelations) {
    fProcessMCRecCorrelationsWithOptions.GetHistogramsList()->AddFirst(new TParameter<Int_t>("AdditionalMCRecOption",fMCRecOption,'f'));
  }

  TH1::AddDirectory(oldstatus);

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
  if (fDoProcessCorrelations) {
    PostData(2, fProcessCorrelations.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
    PostData(3, fProcessMCRecCorrelationsWithOptions.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
    PostData(4, fProcessTrueCorrelations.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
  }
  else {
    /* we post dummy list for avoiding fatal error */
    TList *dummy1 = new TList(); dummy1->SetOwner(kTRUE);
    TList *dummy2 = new TList(); dummy2->SetOwner(kTRUE);
    TList *dummy3 = new TList(); dummy3->SetOwner(kTRUE);
    PostData(2, dummy1); // Post data for ALL output slots >0 here, to get at least an empty histogram
    PostData(3, dummy2); // Post data for ALL output slots >0 here, to get at least an empty histogram
    PostData(4, dummy2); // Post data for ALL output slots >0 here, to get at least an empty histogram
  }
}

/// A new run analysis is starting
/// Check for a potential change in the period (will always happen for the first
/// run in the analysis) and produce all the output depending on the concrete analysis
void AliAnalysisTaskCorrelationsStudies::NotifyRun() {
  AliInfo("");

  AliCSAnalysisCutsBase::NotifyRunGlobal();
  fEventCuts->NotifyRun();
  fTrackSelectionCuts->NotifyRun();

  /* checks the change in the analysis period */
  if (AliCSAnalysisCutsBase::GetGlobalPeriod() == fDataPeriod) return;
  fDataPeriod = AliCSAnalysisCutsBase::GetGlobalPeriod();

  /* now we create additional MC histograms if applicable */
  if (AliCSAnalysisCutsBase::IsMC()) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    Int_t ptbins = fhPt->GetNbinsX();
    Double_t ptlow = fhPt->GetXaxis()->GetBinLowEdge(1), ptup = fhPt->GetXaxis()->GetBinUpEdge(ptbins);

    fhTruePt = new TH1F("fHistTruePt", "P_{T} distribution (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhTruePtPos = new TH1F("fHistTruePtPos", "P_{T} distribution (+) (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhTruePtNeg = new TH1F("fHistTruePtNeg", "P_{T} distribution (-) (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhTrueP = new TH1F("fHistTrueP", "P distribution (truth);P (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhTruePPos = new TH1F("fHistTruePPos", "P distribution (+) (truth);P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);
    fhTruePNeg = new TH1F("fHistTruePNeg", "P distribution (-) (truth);P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);
    fhPurePt = new TH1F("fHistPurePt", "P_{T} distribution for reconstructed (pure);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhPurePtPos = new TH1F("fHistPurePtPos", "P_{T} distribution for reconstructed (+) (pure);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhPurePtNeg = new TH1F("fHistPurePtNeg", "P_{T} distribution for reconstructed (-) (pure);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhPureP = new TH1F("fHistPureP", "P distribution for reconstructed (pure);P (GeV/c);dN/dP_{T} (c/GeV)", ptbins, ptlow, ptup);
    fhPurePPos = new TH1F("fHistPurePPos", "P distribution for reconstructed (+) (pure);P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);
    fhPurePNeg = new TH1F("fHistPurePNeg", "P distribution for reconstructed (-) (pure);P (GeV/c);dN/dP (c/GeV)", ptbins, ptlow, ptup);

    Int_t etabins = fhEtaB->GetNbinsX();
    Double_t etalow = fhEtaB->GetXaxis()->GetBinLowEdge(1), etaup = fhEtaB->GetXaxis()->GetBinUpEdge(etabins);
    Int_t zvtxbins = fhPt3DB->GetNbinsZ();
    Double_t zvtxlow = fhPt3DB->GetZaxis()->GetBinLowEdge(1), zvtxup = fhPt3DB->GetZaxis()->GetBinUpEdge(zvtxbins);
    Int_t phibins = fhPhiB->GetNbinsX();
    Int_t ptvsetabins = fhPtVsEtaB->GetNbinsY();
    Double_t ptvsetalow = fhPtVsEtaB->GetYaxis()->GetBinLowEdge(1), ptvsetaup = fhPtVsEtaB->GetYaxis()->GetBinUpEdge(ptvsetabins);

    fhTrueEta = new TH1F("fHistTrueEta","#eta distribution (truth);#eta;counts",etabins, etalow, etaup);
    fhTruePhi = new TH1F("fHistTruePhi","#phi distribution (truth);#phi;counts",360, 0.0, 360.0);
    fhTrueEtaVsPhi = new TH2F(Form("CSTaskTrueEtaVsPhi_%s",fTaskConfigurationString.Data()),"#eta vs #phi (truth);#phi;#eta", 360, 0.0, 360.0, 100, etalow, etaup);
    fhTruePtVsEta = new TH2F(Form("fhTruePtVsEta_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta (truth);#eta;p_{T} (GeV/c)",etabins,etalow,etaup,ptvsetabins,ptvsetalow,ptvsetaup);

    fhTruePt3D = new TProfile3D(Form("fhTruePt3D_%s",fTaskConfigurationString.Data()),"p_{T} vs #eta, #phi and vtx_{z} (truth);#eta;#phi;vtx_{z}",etabins,etalow,etaup,phibins,0.0,2*TMath::Pi(),zvtxbins,zvtxlow,zvtxup);

    /* the tracking of the event multiplicity distribution for evaluating the lambda angle */
    Int_t mBins[2] = {4000,4000};
    Double_t mLow[2] = {0.0,0.0};
    Double_t mUp[2] = {4000.0,4000.0};
    fhTrueAcceptedVsPrimaries= new THnSparseI("TrueAcceptedVsPrimaries",
        "Accepted vs Total tracks (truth);total tracks,primary tracks,counts",
        2,mBins,mLow,mUp);
    fhTrueLambdaVsPrimaries = new TH2I("fhTrueLambdaVsPrimaries","#lambda (truth);primary tracks;#lambda;counts",
        4000,0,4000,72,-TMath::PiOver2(),TMath::PiOver2());


    fOutput->Add(fhTruePt);
    fOutput->Add(fhTruePtPos);
    fOutput->Add(fhTruePtNeg);
    fOutput->Add(fhTrueP);
    fOutput->Add(fhTruePPos);
    fOutput->Add(fhTruePNeg);
    fOutput->Add(fhPurePt);
    fOutput->Add(fhPurePtPos);
    fOutput->Add(fhPurePtNeg);
    fOutput->Add(fhPureP);
    fOutput->Add(fhPurePPos);
    fOutput->Add(fhPurePNeg);
    fOutput->Add(fhTrueEta);
    fOutput->Add(fhTruePhi);
    fOutput->Add(fhTrueEtaVsPhi);
    fOutput->Add(fhTruePtVsEta);
    fOutput->Add(fhTruePt3D);
    fOutput->Add(fhTrueAcceptedVsPrimaries);
    fOutput->Add(fhTrueLambdaVsPrimaries);

    /* additional histograms if not only truth analysis */
    if (!AliCSAnalysisCutsBase::IsMConlyTruth()) {
      /* the tracking of the true reconstructed but not accepted and of the true not reconstructed */
      Int_t bins[kgTHnDimension] = {ptbins,etabins,phibins,AliPID::kSPECIES+1};
      Double_t low[kgTHnDimension] = {ptlow,etalow,0.0,AliPID::kElectron};
      Double_t up[kgTHnDimension] = {ptup,etaup,360.0,AliPID::kProton+2};
      fhTrueRecNotAccepted = new THnSparseI(Form("fhTrueRecNotAccepted_%s",fTaskConfigurationString.Data()),
          "True reconstructed but not accepted;p_{T} (GeV/c);#eta;#phi;species",
          kgTHnDimension,bins,low,up);
      fhTrueNotReconstructed= new THnSparseI(Form("fhTrueNotReconstructed_%s",fTaskConfigurationString.Data()),
          "True not reconstructed;p_{T} (GeV/c);#eta;#phi;species",
          kgTHnDimension,bins,low,up);

      fhTrueMultiRec = new TH1I(Form("fhTrueMultiRec_%s",fTaskConfigurationString.Data()),
          "Multi reconstructed tracks;no. of multi reconstructed tracks", 20, 0, 20);
      fhPtTrueMultiRec = new TH1I(Form("fhPtTrueMultiRec_%s",fTaskConfigurationString.Data()),
          "p_{T} of multi reconstructed tracks;p_{T} (GeV/c)", ptbins, ptlow, ptup);

      fOutput->Add(fhTrueRecNotAccepted);
      fOutput->Add(fhTrueNotReconstructed);
      fOutput->Add(fhTrueMultiRec);
      fOutput->Add(fhPtTrueMultiRec);
    }

    /* let's build the efficiency profile formula if required */
    if (!BuildEfficiencyProfiles()) {
      AliFatal("Something went wrong when building the efficiency profiles. ABORTING!!!");
    }

    /* if only truth analysis remove the content for the reconstructed output lists */
    if (AliCSAnalysisCutsBase::IsMConlyTruth()) {
      if (fDoProcessCorrelations) {
        fProcessCorrelations.GetHistogramsList()->Clear();
        /* if not additional MC rec results we remove the content from the output list */
        if (fMCRecOption == kNone) {
          fProcessMCRecCorrelationsWithOptions.GetHistogramsList()->Clear();
        }
      }
      if (fDoProcessPairAnalysis) {
        fProcessPairAnalysis.GetHistogramsList()->Clear();
      }
    }
    else {
      /* if not additional MC rec results we remove the content from the output list */
      if (fDoProcessCorrelations) {
        if (fMCRecOption == kNone) {
          fProcessMCRecCorrelationsWithOptions.GetHistogramsList()->Clear();
        }
      }
    }

    TH1::AddDirectory(oldstatus);
  }
  else {
    /* if not MC we remove the content of the output list */
    if (fDoProcessPairAnalysis) {
      fProcessTruePairAnalysis.GetHistogramsList()->Clear();
    }
    if (fDoProcessCorrelations) {
      fProcessMCRecCorrelationsWithOptions.GetHistogramsList()->Clear();
      fProcessTrueCorrelations.GetHistogramsList()->Clear();
    }
  }
}

/// A new event is provided for its analysis
void AliAnalysisTaskCorrelationsStudies::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    AliInfo("Got a new event!");

    // Create pointer to reconstructed event
    AliVEvent *event = InputEvent();
    if (!event) { AliError("ERROR: Could not retrieve event"); return; }

    /* notify the new coming event to our cuts */
    fEventCuts->NotifyEvent();

    /* check if the event is accepted */
    if (!fEventCuts->IsEventAccepted(fInputEvent)) return;

    /* notify the new coming event to our track mapping function */
    aodTrackMaps.NotifyEvent();

    /* notify the new coming event to our cuts */
    fTrackSelectionCuts->NotifyEvent();

    /* build true to rec tracks relation if MC event */
    if (fEventCuts->IsMC()) {
      /* provide some feedback about MC */
      AliInfo(Form("========= MC event is %s", ((fMCEvent != NULL) ? "NOT null" : "NULL")));
      if (fMCEvent != NULL) {
        AliInfo(Form("========= event is 0x%lX while MC event is 0x%lX", (unsigned long) fInputEvent, (unsigned long) fMCEvent));
      }
      if (!fEventCuts->IsMConlyTruth()) {
        /* only if not fast MC */
        BuildTrueRecRelation();
      }
    }

    /* process reconstructed tracks if not a fast MC analysis */
    if (!fEventCuts->IsMConlyTruth()) {
      if (fDoProcessCorrelations && fProcessCorrelations.GetUseSimulation()) {
        /* process rec tracks with potential simulation */
        for (Int_t nsim = 0; nsim < fProcessCorrelations.GetSimEventsPerEvent(); nsim++) {
          /* skip sensible process if additional simulation event is ongoing */
          ProcessTracks(nsim != 0);
        }
      }
      else {
        ProcessTracks(kFALSE);
      }
    }

    /* complete true to rec tracks relation if MC event and not fast MC */
    if (fEventCuts->IsMC() && !fEventCuts->IsMConlyTruth()) {
      BuildTrueRecAccRelation();
    }

    /* process true tracks if MC event */
    if (fEventCuts->IsMC()) {
      ProcessTrueTracks();
    }

    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    AliInfo("Processed event!");
    PostData(1, fOutput);
    if (fDoProcessCorrelations) {
      PostData(2, fProcessCorrelations.GetHistogramsList());
      PostData(3, fProcessMCRecCorrelationsWithOptions.GetHistogramsList());
      PostData(4, fProcessTrueCorrelations.GetHistogramsList());
    }
}

/// \brief Processes the tracks for the current event
/// The processing could be invoked several times for the same event
/// in order to simulate the tracks in the process correlations engine.
/// In that case only the first call should be invoked with kFALSE to fill
/// histograms only once with track real values.
/// \param simulated if kTRUE not fill histograms, tracks will be simulated
void AliAnalysisTaskCorrelationsStudies::ProcessTracks(Bool_t simulated) {
  AliInfo("Start processing event tracks");

  AliVEvent *event = InputEvent();

  Double_t vertexz = fEventCuts->GetVertexZ();
  Double_t centrality = fEventCuts->GetCentrality();

  /* only process pair analysis if this is not an additional simulated event */
  if (fDoProcessPairAnalysis && !simulated) {
    /* initialize pair analysis function activities */
    if (!fProcessPairAnalysis.StartEvent(vertexz)) return;
  }

  /* initialize correlation function activities */
  /* reject event if not correctly initialized */
  if (fDoProcessCorrelations) {
    if (!fProcessCorrelations.StartEvent(centrality,vertexz)) return;
    if (fEventCuts->IsMC()) {
      if (fMCRecOption != kNone) {
        if (!fProcessMCRecCorrelationsWithOptions.StartEvent(centrality,vertexz)) return;
      }
    }
  }

  /* prepare for a potential AOD MC data format */
  TClonesArray *arrayMC = NULL;
  if (fEventCuts->IsMC() && AliCSAnalysisCutsBase::GetMCEventHandler() == NULL) arrayMC = AliCSAnalysisCutsBase::GetMCTrueArray();

  // Track loop for reconstructed event
  Int_t ntracks = event->GetNumberOfTracks();
  Int_t nNoOfAccepted = 0;
  for(Int_t i = 0; i < ntracks; i++) {
    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(event->GetTrack(i)); // pointer to reconstructed to track
    if(!vtrack) {
        AliError(Form("ERROR: Could not retrieve vtrack %d",i));
        continue;
    }

    /* enforce the efficiency profile if required */
    if (fEnforceEfficiencyProfile) {
      if (fRandomGenerator->Uniform(1.0) > ffEfficiencyProfile->Eval(vtrack->Pt()))
        continue;
    }

    Bool_t bTrackAccepted = fTrackSelectionCuts->IsTrackAccepted(vtrack);

    /* only fill histograms if this is not an additional simulated event */
    if (!simulated) {
      fhEtaB->Fill(vtrack->Eta());
      fhPhiB->Fill(vtrack->Phi()*180.0/TMath::Pi());
      fhEtaVsPhiB->Fill(vtrack->Phi()*180.0/TMath::Pi(),vtrack->Eta());
      fhPtVsEtaB->Fill(vtrack->Eta(),vtrack->Pt());
      fhPt3DB->Fill(vtrack->Eta(),vtrack->Phi(),vertexz,vtrack->Pt());
      fh3Dn1B->Fill(vtrack->Eta(),vtrack->Phi(),vtrack->Pt());
    }

    if (bTrackAccepted) {
      /* get the original track in case of constrained track */
      Int_t esdid = vtrack->GetID();
      AliVTrack *origtrack = vtrack;
      if (esdid < 0) origtrack = aodTrackMaps.GetOriginalTrack((AliAODTrack*) vtrack);

      /* additional tracks quality tests for MC reconstructed */
      if (fEventCuts->IsMC()) {
        fMCRecFlags[origtrack->GetID()] |= kAccepted;
        if (fMCRecFlags[origtrack->GetID()] != kAccepted && (origtrack->GetLabel() > 0)) {
          AliError(Form("ERROR: Accepted reconstructed track with esdid %d would be rejected due to quality tests %08X",esdid,fMCRecFlags[origtrack->GetID()]));
//          continue;
        }
      }
      /* only fill histograms if this is not an additional simulated event */
      if (!simulated) {
        nNoOfAccepted++;
        fhPt->Fill(vtrack->Pt());
        fhP->Fill(vtrack->P());
        if (vtrack->Charge() > 0) {
          fhPtPos->Fill(vtrack->Pt());
          fhPPos->Fill(vtrack->P());
        }
        if (vtrack->Charge() < 0) {
          fhPtNeg->Fill(vtrack->Pt());
          fhPNeg->Fill(vtrack->P());
        }
        fhEtaA->Fill(vtrack->Eta());
        fhPhiA->Fill(vtrack->Phi()*180.0/TMath::Pi());
        fhEtaVsPhiA->Fill(vtrack->Phi()*180.0/TMath::Pi(),vtrack->Eta());
        fhPtVsEtaA->Fill(vtrack->Eta(),vtrack->Pt());
        fhPt3DA->Fill(vtrack->Eta(),vtrack->Phi(),vertexz,vtrack->Pt());
        fh3Dn1A->Fill(vtrack->Eta(),vtrack->Phi(),vtrack->Pt());
        fhLambdaVsMultiplicity->Fill(ntracks,TMath::PiOver2() - vtrack->Theta());

        /* fill the constrained vs not constrained histograms */
        if (esdid < 0) {
          fhUnConstrainedPt->Fill(origtrack->Pt());
          fhPtDifference->Fill(vtrack->Pt() - origtrack->Pt());               ///< Constrained - unconstrained \f$ p_{T} \f$ spectrum
        }
      }

      /* only process pair analysis if this is not an additional simulated event */
      if (fDoProcessPairAnalysis && !simulated) {
        /* track accepted, process pair analysis function activities */
        fProcessPairAnalysis.ProcessTrack(i, vtrack);
      }

      /* track accepted, process correlation function activities */
      if (fDoProcessCorrelations) {
        fProcessCorrelations.ProcessTrack(i, vtrack);
        /* only process rec with true values if MC and this is not an additional simulated event */
        if (fEventCuts->IsMC() && !simulated) {
          switch (fMCRecOption) {
          case kRecWithTrue:
            /* additional results with reconstructed with true values */
            if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
              fProcessMCRecCorrelationsWithOptions.ProcessTrack(i, fMCEvent->GetTrack(vtrack->GetLabel()));
            }
            else {
              fProcessMCRecCorrelationsWithOptions.ProcessTrack(i, (AliVParticle *) arrayMC->At(vtrack->GetLabel()));
            }
            break;
          case kRecTruePrimaries:
            /* additional results with true primary reconstructed tracks */
            if (fTrackSelectionCuts->IsTruePrimary(vtrack)) {
              fProcessMCRecCorrelationsWithOptions.ProcessTrack(i, vtrack);
            }
            break;
          case kRecWithNotAccepted:
            /* additional results with accepted and not accepted reconstructed tracks */
            fProcessMCRecCorrelationsWithOptions.ProcessTrack(i, vtrack);
            break;
          default:
            break;
          }
        }
      }

      /* only process purity analysis if this is not an additional simulated event */
      if (!simulated) {
        /* store purity information if applicable */
        if (fEventCuts->IsMC()) {
          if (fTrackSelectionCuts->IsTrueTrackAccepted(vtrack)) {
            fhPurePt->Fill(vtrack->Pt());
            fhPureP->Fill(vtrack->P());
            if (vtrack->Charge() > 0) {
              fhPurePtPos->Fill(vtrack->Pt());
              fhPurePPos->Fill(vtrack->P());
            }
            else {
              fhPurePtNeg->Fill(vtrack->Pt());
              fhPurePNeg->Fill(vtrack->P());
            }
          }
        }
      }
    }
    else {
      /* track not accepted, process required correlation function activities */
      if (fEventCuts->IsMC() && !simulated) {
        if (fDoProcessCorrelations) {
          switch (fMCRecOption) {
          case kRecWithNotAccepted:
            /* additional results with accepted and not accepted reconstructed tracks */
            /* but only with its corresponding true accepted, for the time being */
            if (fTrackSelectionCuts->IsTrueTrackAccepted(vtrack)) {
              fProcessMCRecCorrelationsWithOptions.ProcessTrack(i, vtrack);
            }
            break;
          default:
            break;
          }
        }
      }
    }
  }
  /* only process pair analysis if this is not an additional simulated event */
  if (fDoProcessPairAnalysis && !simulated) {
    /* and finalize pair analysis function activities */
    fProcessPairAnalysis.ProcessEventData();
  }

  /* and finalize correlation function activities */
  if (fDoProcessCorrelations) {
    fProcessCorrelations.ProcessEventData();
    /* only process rec with true if MC and this is not an additional simulated event */
    if (fEventCuts->IsMC() && !simulated) {
      if (fMCRecOption != kNone) {
        fProcessMCRecCorrelationsWithOptions.ProcessEventData();
      }
    }
  }
  /* track the event multiplicity */
  Double_t fillmdata[2] = {Double_t(ntracks+0.5),Double_t(nNoOfAccepted+0.5)};
  fhAcceptedVsMultiplicity->Fill(fillmdata);

  AliInfo(Form("TOTAL REC ACCEPTED: %d", nNoOfAccepted));
}

/// \brief Processes the true tracks for the current event
void AliAnalysisTaskCorrelationsStudies::ProcessTrueTracks() {

  AliInfo("");

  Double_t vertexz = 0.0;
  if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
    vertexz = fMCEvent->GetPrimaryVertex()->GetZ();
  }
  else {
    AliVEvent *event = InputEvent();
    AliAODMCHeader *header = (AliAODMCHeader *) event->FindListObject(AliAODMCHeader::StdBranchName());
    vertexz = header->GetVtxZ();
  }

  Double_t centrality = fEventCuts->GetCentrality();

  if (fDoProcessPairAnalysis) {
    /* initialize pair analysis function activities */
    if (!fProcessTruePairAnalysis.StartEvent(vertexz)) return;
  }

 /* initialize correlation function activities */
  /* reject event if not correctly initialized */
  if (fDoProcessCorrelations) {
    if (!fProcessTrueCorrelations.StartEvent(centrality,vertexz)) return;
  }

  Int_t nNoOfTrueReconstructedAndAccepted = 0;
  Int_t nNoOfTrueAccepted = 0;
  Int_t nNoOfTrueMultiReconstructedAndAccepted = 0;
  Int_t nNoOfReconstructed = InputEvent()->GetNumberOfTracks();
  Int_t nNoOfTrue = GetNoOfTrueParticles();
  Int_t nNoOfPrimaries = GetNoOfTruePrimaries();


  Int_t ntracks = GetNoOfTrueParticles();
  TClonesArray *arrayMC = NULL;
  if (AliCSAnalysisCutsBase::GetMCEventHandler() == NULL) arrayMC = AliCSAnalysisCutsBase::GetMCTrueArray();
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++ ){

    if (fTrackSelectionCuts->IsTrueTrackAccepted(iTrack)) {
      AliVParticle *part = NULL;
      Double_t p = 0.0;
      Double_t pt = 0.0;
      Double_t pz = 0.0;
      Double_t eta = 0.0;
      Double_t phi = 0.0;
      Double_t charge = 0.0;
      AliPID::EParticleType species = AliPID::kUnknown;

      /* extract the true particle according to the ESD or AOD format */
      if (AliCSAnalysisCutsBase::GetMCEventHandler() != NULL) {
        part = fMCEvent->GetTrack(iTrack);
      }
      else {
        part = (AliAODMCParticle *) arrayMC->At(iTrack);
      }

      p = part->P();
      pt = part->Pt();
      pz = part->Pz();
      eta = part->Eta();
      phi = part->Phi();
      charge = part->Charge();
      species = AliCSPIDCuts::GetTrueSpecies(part);

      Double_t filldata[kgTHnDimension] = {pt,eta,phi*180.0/TMath::Pi(),
          (AliPID::kProton < species) ?  AliPID::kSPECIES + 0.5 : species + 0.5};

      /* the analysis just for the charged particles */
      if (charge != 0.0) {
        /* enforce efficiency on true if required */
        switch (fOnTrueEfficiencyProfile) {
        case kEfficiencyProfile:
          /* use the stored efficiency profile */
          if (charge > 0) {
            if (fRandomGenerator->Uniform(1.0) > fhOnTrueEfficiencyProfile_1->GetBinContent(
                fhOnTrueEfficiencyProfile_1->FindBin(eta,pt))) {
              continue;
            }
          }
          else {
            if (fRandomGenerator->Uniform(1.0) > fhOnTrueEfficiencyProfile_2->GetBinContent(
                fhOnTrueEfficiencyProfile_2->FindBin(eta,pt))) {
              continue;
            }
          }
          break;
        case kOnlyReconstructed:
          /* incorporate only tracks which were reconstructed */
          if (fTrueToRec->At(iTrack) == NULL) {
            /* not reconstructed so, skip it */
            continue;
          }
          break;
        default:
          break;
        }

        nNoOfTrueAccepted++;
        if (!fEventCuts->IsMConlyTruth()) {
          /* only if not fast MC */
          if (fTrueToRec->At(iTrack) != NULL) {
            /* true was reconstructed, let's check whether it was accepted */
            if ((fMCRecFlags[((AliVTrack *) fTrueToRec->At(iTrack))->GetID()] & kAccepted) == kAccepted) {
              /* the reconstructed was accepted so, what? */
              nNoOfTrueReconstructedAndAccepted++;
              if ((fMCRecFlags[((AliVTrack *) fTrueToRec->At(iTrack))->GetID()] & kSamePositiveLabel) == kSamePositiveLabel) {
                nNoOfTrueMultiReconstructedAndAccepted++;
                fhPtTrueMultiRec->Fill(pt);
              }
            }
            else {
              /* the reconstructed was not accepted */
              fhTrueRecNotAccepted->Fill(filldata);
            }
          }
          else {
            /* the track was not reconstructed */
            fhTrueNotReconstructed->Fill(filldata);
          }
        }

        fhTruePt->Fill(pt);
        fhTrueP->Fill(p);
        fhTrueEta->Fill(eta);
        fhTruePhi->Fill(phi*180.0/TMath::Pi());
        fhTrueEtaVsPhi->Fill(phi*180.0/TMath::Pi(), eta);
        fhTruePtVsEta->Fill(eta,pt);
        fhTruePt3D->Fill(eta,phi,vertexz,pt);
        if (charge > 0) {
          fhTruePtPos->Fill(pt);
          fhTruePPos->Fill(p);
          Double_t theta = ((pz == 0.0) ? TMath::PiOver2() : TMath::ACos(pz/p));
          fhTrueLambdaVsPrimaries->Fill(nNoOfPrimaries,TMath::PiOver2() - theta);
        }
        else {
          fhTruePtNeg->Fill(pt);
          fhTruePNeg->Fill(p);
          Double_t theta = ((pz == 0) ? TMath::PiOver2() : TMath::ACos(pz/p));
          fhTrueLambdaVsPrimaries->Fill(nNoOfPrimaries,TMath::PiOver2() - theta);
        }
        if (fDoProcessPairAnalysis) {
          /* track accepted, process pair analysis function activities */
           fProcessTruePairAnalysis.ProcessTrack(iTrack, part);
        }
        /* track accepted, process correlation function activities */
        if (fDoProcessCorrelations) {
           fProcessTrueCorrelations.ProcessTrack(iTrack, part);
        }
      }
    }
  }
  if (fDoProcessPairAnalysis) {
    /* and finalize pair analysis function activities */
    fProcessTruePairAnalysis.ProcessEventData();
  }
  /* and finalize correlation function activities */
  if (fDoProcessCorrelations) {
    fProcessTrueCorrelations.ProcessEventData();
  }

  if (!fEventCuts->IsMConlyTruth()) {
    /* only if not fast MC */
    if (nNoOfTrueMultiReconstructedAndAccepted != 0)
      fhTrueMultiRec->Fill(nNoOfTrueMultiReconstructedAndAccepted+0.5);
  }

  /* track the event multiplicity */
  Double_t fillmdata[2] = {Double_t(nNoOfPrimaries+0.5),Double_t(nNoOfTrueAccepted+0.5)};
  fhTrueAcceptedVsPrimaries->Fill(fillmdata);

  AliInfo(Form("TOTAL TRUE: %d; TOTAL TRUE PRIMARIES: %d; TOTAL TRUE ACCEPTED: %d",
      nNoOfTrue, nNoOfPrimaries, nNoOfTrueAccepted));
  AliInfo(Form("TOTAL REC: %d; TOTAL TRUE REC AND ACC: %d; TOTAL TRUE ACC MULTI REC: %d",
      nNoOfReconstructed, nNoOfTrueReconstructedAndAccepted, nNoOfTrueMultiReconstructedAndAccepted));
}

/// The analysis task is going to finish. Produce the final results to be preserved.
void AliAnalysisTaskCorrelationsStudies::FinishTaskOutput() {

  AliInfo("Starting!");

  if (fDoProcessPairAnalysis) {
    if (!fEventCuts->IsMConlyTruth()) {
      fProcessPairAnalysis.FinalizeProcess();
    }
    if (fEventCuts->IsMC()) {
      fProcessTruePairAnalysis.FinalizeProcess();
    }
  }

  if (fDoProcessCorrelations) {
    if (!fEventCuts->IsMConlyTruth()) {
      fProcessCorrelations.FinalizeProcess();
    }
    if (fEventCuts->IsMC()) {
      if (fMCRecOption != kNone) {
        fProcessMCRecCorrelationsWithOptions.FinalizeProcess();
      }
      fProcessTrueCorrelations.FinalizeProcess();
    }
  }

  if (fDoProcessCorrelations) {
    if (!fEventCuts->IsMConlyTruth()) {
      PostData(2, fProcessCorrelations.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
    }
    if (fEventCuts->IsMC()) {
      if (fMCRecOption != kNone) {
        PostData(3, fProcessMCRecCorrelationsWithOptions.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
      }
      PostData(4, fProcessTrueCorrelations.GetHistogramsList()); // Post data for ALL output slots >0 here, to get at least an empty histogram
    }
  }

  AliInfo("Done!");
}

/// The local task is terminated. We do nothing for the time being
void AliAnalysisTaskCorrelationsStudies::Terminate(Option_t *)
{
  /* take it away for the time being */
  if (kFALSE) {
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutput) { AliError("ERROR: could not retrieve TList fOutput"); return; }

    fhPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistPt"));
    if (!fhPt) { AliError("ERROR: could not retrieve fHistPt"); return;}
    fhEtaA = dynamic_cast<TH1F*> (fOutput->FindObject("fHistEtaA"));
    if (!fhEtaA) { AliError("ERROR: could not retrieve fHistEtaA"); return;}

    // Get the physics selection histograms with the selection statistics
    //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    //TH2F *histStat = (TH2F*)inputH->GetStatistics();


    // NEW HISTO should be retrieved from the TList container in the above way,
    // so it is available to draw on a canvas such as below

    TCanvas *c = new TCanvas("AliAnalysisTaskCorrelationsStudies","P_{T} & #eta",10,10,1020,510);
    c->Divide(2,1);
    c->cd(1)->SetLogy();
    fhPt->DrawCopy("E");
    c->cd(2);
    fhEtaA->DrawCopy("E");
  }
}
