/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//                  Base class for DStar Analysis
//
//
//  The D* spectra study is done in pt bins:
//  [0,0.5] [0.5,1] [1,2] [2,3] [3,4] [4,5] [5,6] [6,7] [7,8],
//  [8,10],[10,12], [12,16], [16,20] and [20,24]
//
//  Cuts arew centralized in AliRDHFCutsDStartoKpipi
//  Side Band and like sign background are implemented in the macro
//
//-----------------------------------------------------------------------
//
//                         Author A.Grelli
//              ERC-QGP Utrecht University - a.grelli@uu.nl,
//                         Author Y.Wang
//        University of Heidelberg - yifei@physi.uni-heidelberg.de
//                         Author C.Ivan
//             ERC-QGP Utrecht University - c.ivan@uu.nl,
//
//                  modified for EMCAL production check
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEDStarEMCALProductionCheck.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliInputEventHandler.h"
#include "AliTrackerBase.h"
#include "AliGenPythiaEventHeader.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDStarEMCALProductionCheck);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEDStarEMCALProductionCheck::AliAnalysisTaskSEDStarEMCALProductionCheck():
  AliAnalysisTaskSE(),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE),
  fDoSearch(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPID(0),
  fOutputProductionCheck(0),
  fNSigma(3),
  fCuts(0),
  fCEvents(0),
  fTrueDiff2(0),
  fDeltaMassD1(0),
  fCounter(0),
  fAODProtection(1),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fNPtBins(0),
  fAllhist(0x0),
  fPIDhist(0x0),
  fDoDStarVsY(kFALSE),
  fUseEMCalTrigger(kFALSE),
  fTriggerSelectionString(0),
  fCheckEMCALAcceptance(kFALSE),
  fCheckEMCALAcceptanceNumber(0),
  fApplyEMCALClusterEventCut(kFALSE)
{
  //
  /// Default ctor
  //
  for (Int_t i = 0; i < 5; i++) fHistMassPtImpParTCDs[i] = 0;
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarEMCALProductionCheck::AliAnalysisTaskSEDStarEMCALProductionCheck(const Char_t* name, AliRDHFCutsDStartoKpipi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE),
  fDoSearch(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPID(0),
  fOutputProductionCheck(0),
  fNSigma(3),
  fCuts(0),
  fCEvents(0),
  fTrueDiff2(0),
  fDeltaMassD1(0),
  fCounter(0),
  fAODProtection(1),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fNPtBins(0),
  fAllhist(0x0),
  fPIDhist(0x0),
  fDoDStarVsY(kFALSE),
  fUseEMCalTrigger(kFALSE),
  fTriggerSelectionString(0),
  fCheckEMCALAcceptance(kFALSE),
  fCheckEMCALAcceptanceNumber(0),
  fApplyEMCALClusterEventCut(kFALSE)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarEMCALProductionCheck", "Calling Constructor");

  fCuts = cuts;
  for (Int_t i = 0; i < 5; i++) fHistMassPtImpParTCDs[i] = 0;


  DefineOutput(1, TList::Class()); //counters
  DefineOutput(2, TList::Class()); //All Entries output
  DefineOutput(3, TList::Class()); //3sigma PID output
  DefineOutput(4, AliRDHFCutsDStartoKpipi::Class());  //My private output
  DefineOutput(5, AliNormalizationCounter::Class());  // normalization
  DefineOutput(6, TList::Class()); //production check
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarEMCALProductionCheck::~AliAnalysisTaskSEDStarEMCALProductionCheck() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEDStarEMCALProductionCheck", "Calling Destructor");

  delete fOutput;
  delete fOutputAll;
  delete fOutputPID;
  delete fOutputProductionCheck;
  delete fCuts;
  delete fCEvents;
  delete fDeltaMassD1;
  for (Int_t i = 0; i < 5; i++) {
    delete fHistMassPtImpParTCDs[i];
  }
  for (Int_t i = 0; i < ((fNPtBins + 2) * 18); i++) {
    delete fAllhist[i];
    delete fPIDhist[i];
  }
  delete [] fAllhist;
  delete [] fPIDhist;

}
//_________________________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::Init() {
  //
  /// Initialization
  //

  if (fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");
  AliRDHFCutsDStartoKpipi* copyfCuts = new AliRDHFCutsDStartoKpipi(*fCuts);
  fNPtBins = fCuts->GetNPtBins();
  // Post the data
  PostData(4, copyfCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::UserExec(Option_t *)
{
  /// user exec
  if (!fInputEvent) {
    Error("UserExec", "NO EVENT FOUND!");
    return;
  }

  fCEvents->Fill(0);//all events
  if (fAODProtection >= 0) {
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fCEvents->Fill(8);
      return;
    }
    fCEvents->Fill(1);
  }

  fEvents++;

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi = 0;
  TClonesArray *arrayD0toKpi = 0;

  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
                                ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
      arrayD0toKpi = (TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else {
    arrayDStartoD0pi = (TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
    arrayD0toKpi = (TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
  }


  //objects for production check
  AliAODMCHeader *mcHeader = nullptr;
  AliGenPythiaEventHeader * pythiaHeader = nullptr;
  TClonesArray *mcTrackArray = nullptr;
  Double_t crossSection = 0.0;
  Double_t ptHard = 0.0;
  Int_t nTrials = 0;

  // Int_t nPtBins = fCuts->GetNPtBins();
  // const Int_t nPtBinLimits = nPtBins + 1;
  // Int_t PtBinLimits[nPtBinLimits] = fCuts->GetPtBinLimits();

  if (fUseMCInfo) {
    // load MC header
    mcHeader =  (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      printf("AliAnalysisTaskSEDStarEMCALProductionCheck::UserExec: MC header branch not found!\n");
      return;
    }

    AliGenPythiaEventHeader * pythiaHeader = (AliGenPythiaEventHeader*)mcHeader->GetCocktailHeader(0);
    if (!pythiaHeader) {
      printf("AliAnalysisTaskSEDStarEMCALProductionCheck::UserExec: AliGenPythiaEventHeader not found!\n");
      return;
    }
    crossSection = pythiaHeader->GetXsection();
    ptHard = pythiaHeader->GetPtHard();
    nTrials = pythiaHeader->Trials();

    mcTrackArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcTrackArray) {std::cout << "no track array" << std::endl; return;};
  }

  // check before event cut
  if (fUseMCInfo) {
    for (Int_t j = 0; j < mcTrackArray->GetEntriesFast(); j++) {

      AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(j));
      if (!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
      Int_t pdgCodeMC = TMath::Abs(mcTrackParticle->GetPdgCode());

      if (pdgCodeMC == 413)
      { //if the track is a DStar we check if it comes from charm

        Double_t ptMC = mcTrackParticle->Pt();

        Bool_t fromCharm = kFALSE;
        Int_t mother = mcTrackParticle->GetMother();
        Int_t istep = 0;
        while (mother >= 0 ) {
          istep++;
          AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(mcTrackArray->At(mother));
          if (mcGranma) {
            Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
            if ((abspdgGranma == 4) || (abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)) fromCharm = kTRUE;
            mother = mcGranma->GetMother();
          } else {
            printf("AliVertexingHFUtils::IsTrackFromCharm: Failed casting the mother particle!");
            break;
          }
        }

        if (fromCharm)
        {

          Bool_t mcPionDStarPresent = kFALSE;
          Bool_t mcPionD0Present = kFALSE;
          Bool_t mcKaonPresent = kFALSE;
          Int_t nDaughterDStar = mcTrackParticle->GetNDaughters();

          if (nDaughterDStar == 2) {
            for (Int_t iDaughterDStar = 0; iDaughterDStar < 2; iDaughterDStar++) {

              AliAODMCParticle* daughterDStar = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(iDaughterDStar));
              if (!daughterDStar) break;
              Int_t pdgCodeDaughterDStar = TMath::Abs(daughterDStar->GetPdgCode());

              if (pdgCodeDaughterDStar == 211) { //if the track is a pion we save its monte carlo label
                mcPionDStarPresent = kTRUE;

              } else if (pdgCodeDaughterDStar == 421) { //if the track is a D0 we look at its daughters
                Int_t mcLabelD0 = mcTrackParticle->GetDaughterLabel(iDaughterDStar);
                Int_t nDaughterD0 = daughterDStar->GetNDaughters();

                if (nDaughterD0 == 2) {
                  for (Int_t iDaughterD0 = 0; iDaughterD0 < 2; iDaughterD0++) {

                    AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterDStar->GetDaughterLabel(iDaughterD0));
                    if (!daughterD0) break;
                    Int_t pdgCodeDaughterD0 = TMath::Abs(daughterD0->GetPdgCode());

                    if (pdgCodeDaughterD0 == 211) {
                      mcPionD0Present = kTRUE;
                    } else if (pdgCodeDaughterD0 == 321) {
                      mcKaonPresent = kTRUE;
                    } else break;
                  }
                }
              } else break;
            }
          }

          if (mcPionDStarPresent && mcPionD0Present && mcKaonPresent)
          {
            TString fillthis = "";
            fillthis = "DStarPtTruePreEventSelection";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
            fillthis = "DStarPtTruePreEventSelectionWeighted";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);
            // fillthis = "PtHardPreEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard);
            // fillthis = "PtHardWeightedPreEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard, crossSection);
            // fillthis = "WeightsPreEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(crossSection);
            // fillthis = "TrialsPreEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->AddBinContent(1,nTrials);
            fillthis = "DStar_per_bin_true_PreEventSelection";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
            fillthis = "DStar_per_bin_true_PreEventSelection_weighted";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);          
          }
        }
      }
    }
  }

  TString fillthishist = "";
  fillthishist = "PtHardPreEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard);
  fillthishist = "PtHardWeightedPreEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard, crossSection);
  fillthishist = "WeightsPreEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(crossSection);
  fillthishist = "TrialsPreEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->AddBinContent(1,nTrials);



  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField()) < 0.001) return;
  fCEvents->Fill(2);

  fCounter->StoreEvent(aodEvent, fCuts, fUseMCInfo);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass = aodEvent->GetFiredTriggerClasses();
  if (trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  if (!fCuts->IsEventSelected(aodEvent)) {
    if (fCuts->GetWhyRejection() == 6) // rejected for Z vertex
      fCEvents->Fill(6);
    return;
  }

  // Use simulated EMCal trigger for MC. AliEmcalTriggerMakerTask needs to be run first.
  if (fUseEMCalTrigger)
  {
    auto triggercont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer*>(fInputEvent->FindListObject("EmcalTriggerDecision"));
    if (!triggercont)
    {
      AliErrorStream() <<  "Trigger decision container not found in event - not possible to select EMCAL triggers" << std::endl;
      return;
    }
    if (fTriggerSelectionString == "EG1DG1")
    {
      if (!triggercont->IsEventSelected("EG1") && !triggercont->IsEventSelected("DG1")) return;
    } else if (!triggercont->IsEventSelected(fTriggerSelectionString)) return;
  }

  // Get field for EMCAL acceptance and cut events
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  inputHandler->SetNeedField();

  if (fApplyEMCALClusterEventCut)
  {
    Int_t numberOfCaloClustersEvent = aodEvent->GetNumberOfCaloClusters();

    if (numberOfCaloClustersEvent >= 0)
    {
      Bool_t passClusterCuts = kFALSE;
      for (Int_t iCluster = 0; iCluster < numberOfCaloClustersEvent; ++iCluster)
      {
        AliAODCaloCluster * trackEMCALCluster = (AliAODCaloCluster*)aodEvent->GetCaloCluster(iCluster);
        if (trackEMCALCluster->GetNonLinCorrEnergy() < 9.0) continue;
        if (trackEMCALCluster->GetTOF() > 15e-9) continue;
        if (trackEMCALCluster->GetTOF() < -20e-9) continue;
        if (trackEMCALCluster->GetIsExotic()) continue;
        passClusterCuts = kTRUE;
      }
      if (!passClusterCuts) return;
    } else return;
  }

  Bool_t isEvSel = fCuts->IsEventSelected(aodEvent);
  fCEvents->Fill(3);
  if (!isEvSel) return;

  // Load the event
  //  AliInfo(Form("Event %d",fEvents));
  //if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));

  // counters for efficiencies
  Int_t icountReco = 0;

  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2] = {421, 211};
  Int_t pdgDgD0toKpi[2] = {321, 211};

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!vtx1) return;
  if (vtx1->GetNContributors() < 1) return;
  fCEvents->Fill(4);


  //save cluster information for EMCal trigger selection check

  Int_t numberOfCaloClustersEvent = aodEvent->GetNumberOfCaloClusters();

  if (numberOfCaloClustersEvent >= 0)
  {
    for (Int_t iCluster = 0; iCluster < numberOfCaloClustersEvent; ++iCluster)
    {
      AliAODCaloCluster * trackEMCALCluster = (AliAODCaloCluster*)aodEvent->GetCaloCluster(iCluster);

      //save cluster information
      Float_t pos[3]={0};
      trackEMCALCluster->GetPosition(pos);
      TString fillthis = "";
      fillthis = "fHistClusPosition";
      ((TH3F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(pos[0], pos[1], pos[2]);
    }
  } 


  if (!arrayDStartoD0pi || !arrayD0toKpi) {
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  } else AliDebug(2, Form("Found %d vertices", arrayDStartoD0pi->GetEntriesFast()));

  Int_t nSelectedAna = 0;
  Int_t nSelectedProd = 0;

  // check after event cut
  if (fUseMCInfo) {
    for (Int_t j = 0; j < mcTrackArray->GetEntriesFast(); j++) {

      AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(j));
      if (!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
      Int_t pdgCodeMC = TMath::Abs(mcTrackParticle->GetPdgCode());

      if (pdgCodeMC == 413)
      { //if the track is a DStar we check if it comes from charm

        Double_t ptMC = mcTrackParticle->Pt();

        Bool_t fromCharm = kFALSE;
        Int_t mother = mcTrackParticle->GetMother();
        Int_t istep = 0;
        while (mother >= 0 ) {
          istep++;
          AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(mcTrackArray->At(mother));
          if (mcGranma) {
            Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
            if ((abspdgGranma == 4) || (abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)) fromCharm = kTRUE;
            mother = mcGranma->GetMother();
          } else {
            printf("AliVertexingHFUtils::IsTrackFromCharm: Failed casting the mother particle!");
            break;
          }
        }

        if (fromCharm)
        {

          Bool_t mcPionDStarPresent = kFALSE;
          Bool_t mcPionD0Present = kFALSE;
          Bool_t mcKaonPresent = kFALSE;
          Int_t nDaughterDStar = mcTrackParticle->GetNDaughters();

          if (nDaughterDStar == 2) {
            for (Int_t iDaughterDStar = 0; iDaughterDStar < 2; iDaughterDStar++) {

              AliAODMCParticle* daughterDStar = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(iDaughterDStar));
              if (!daughterDStar) break;
              Int_t pdgCodeDaughterDStar = TMath::Abs(daughterDStar->GetPdgCode());

              if (pdgCodeDaughterDStar == 211) { //if the track is a pion we save its monte carlo label
                mcPionDStarPresent = kTRUE;

              } else if (pdgCodeDaughterDStar == 421) { //if the track is a D0 we look at its daughters
                Int_t mcLabelD0 = mcTrackParticle->GetDaughterLabel(iDaughterDStar);
                Int_t nDaughterD0 = daughterDStar->GetNDaughters();

                if (nDaughterD0 == 2) {
                  for (Int_t iDaughterD0 = 0; iDaughterD0 < 2; iDaughterD0++) {

                    AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterDStar->GetDaughterLabel(iDaughterD0));
                    if (!daughterD0) break;
                    Int_t pdgCodeDaughterD0 = TMath::Abs(daughterD0->GetPdgCode());

                    if (pdgCodeDaughterD0 == 211) {
                      mcPionD0Present = kTRUE;
                    } else if (pdgCodeDaughterD0 == 321) {
                      mcKaonPresent = kTRUE;
                    } else break;
                  }
                }
              } else break;
            }
          }

          if (mcPionDStarPresent && mcPionD0Present && mcKaonPresent)
          {
            TString fillthis = "";
            fillthis = "DStarPtTruePostEventSelection";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
            fillthis = "DStarPtTruePostEventSelectionWeighted";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);
            // fillthis = "PtHardPostEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard);
            // fillthis = "PtHardWeightedPostEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard, crossSection);
            // fillthis = "WeightsPostEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(crossSection);
            // fillthis = "TrialsPostEventSelection";
            // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->AddBinContent(1,nTrials);
            fillthis = "DStar_per_bin_true_PostEventSelection";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
            fillthis = "DStar_per_bin_true_PostEventSelection_weighted";
            ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);        
          }
        }
      }
    }
  }

  fillthishist = "";
  fillthishist = "PtHardPostEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard);
  fillthishist = "PtHardWeightedPostEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard, crossSection);
  fillthishist = "WeightsPostEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(crossSection);
  fillthishist = "TrialsPostEventSelection";
  ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->AddBinContent(1,nTrials);

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  // loop over the tracks to search for candidates soft pion
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi < arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {

    // D* candidates and D0 from D*
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    AliAODRecoDecayHF2Prong *trackD0;
    if (dstarD0pi->GetIsFilled() < 1) {
      trackD0 = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->At(dstarD0pi->GetProngID(1));
    } else {
      trackD0 = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    }

    fCEvents->Fill(10);
    TObjArray arrTracks(3);
    for (Int_t ipr = 0; ipr < 3; ipr++) {
      AliAODTrack *tr;
      if (ipr == 0) tr = vHF->GetProng(aodEvent, dstarD0pi, ipr); //soft pion
      else         tr = vHF->GetProng(aodEvent, trackD0, ipr - 1); //D0 daughters
      arrTracks.AddAt(tr, ipr);
    }
    if (!fCuts->PreSelect(arrTracks)) {
      fCEvents->Fill(13);
      continue;
    }

    Bool_t isDStarCand = kTRUE;
    if (!(vHF->FillRecoCasc(aodEvent, dstarD0pi, isDStarCand))) { //Fill the data members of the candidate only if they are empty.
      fCEvents->Fill(12); //monitor how often this fails
      continue;
    }
    if (!dstarD0pi->GetSecondaryVtx()) continue;
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    if (!theD0particle) continue;

    Int_t isDStar = 0;
    TClonesArray *mcArray = 0;
    // AliAODMCHeader *mcHeader = 0;

    Bool_t isPrimary = kTRUE;
    Float_t pdgCode = -2;
    Float_t trueImpParXY = 0.;

    // mc analysis
    if (fUseMCInfo) {
      //MC array need for maching
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
        AliError("Could not find Monte-Carlo in AOD");
        return;
      }
      // load MC header
      // mcHeader =  (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      // if (!mcHeader) {
      //   printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
      //   return;
      // }
      // find associated MC particle for D* ->D0toKpi
      Int_t mcLabel = dstarD0pi->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, mcArray);
      if (mcLabel >= 0) {

        AliAODMCParticle *partDSt = (AliAODMCParticle*)mcArray->At(mcLabel);
        Int_t checkOrigin = CheckOrigin(mcArray, partDSt);
        if (checkOrigin == 5) isPrimary = kFALSE;
        AliAODMCParticle *dg0 = (AliAODMCParticle*)mcArray->At(partDSt->GetDaughterLabel(0));
        //  AliAODMCParticle *dg01 = (AliAODMCParticle*)mcArray->At(dg0->GetDaughterLabel(0));
        pdgCode = TMath::Abs(partDSt->GetPdgCode());
        if (!isPrimary) {
          trueImpParXY = GetTrueImpactParameterD0(mcHeader, mcArray, dg0) * 1000.;
        }
        isDStar = 1;
      } else {
        pdgCode = -1;
      }
    }

    if (pdgCode == -1) AliDebug(2, "No particle assigned! check\n");

    Double_t Dstarpt = dstarD0pi->Pt();

    // quality selction on tracks and region of interest
    Int_t isTkSelected = fCuts->IsSelected(dstarD0pi, AliRDHFCuts::kTracks); // quality cuts on tracks
    if (!isTkSelected) continue;

    if (!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(), dstarD0pi->YDstar())) continue;


    // EMCAL acceptance check
    if (fCheckEMCALAcceptance)
    {
      Int_t numberInAcc = 0;

      AliAODTrack *track[3];
      for (Int_t iDaught = 0; iDaught < 3; iDaught++) {
        track[iDaught] = (AliAODTrack*)arrTracks.At(iDaught);

        Int_t numberOfCaloClusters = aodEvent->GetNumberOfCaloClusters();

        if (numberOfCaloClusters >= 0)
        {
          Int_t trackEMCALClusterNumber = track[iDaught]->GetEMCALcluster();
          if (!(trackEMCALClusterNumber < 0))
          {
            AliAODCaloCluster * trackEMCALCluster = (AliAODCaloCluster*)aodEvent->GetCaloCluster(trackEMCALClusterNumber);
            if (!trackEMCALCluster) continue;

            if (trackEMCALCluster->GetNonLinCorrEnergy() < 9.0) continue;
            if (trackEMCALCluster->GetTOF() > 15e-9) continue;
            if (trackEMCALCluster->GetTOF() < -20e-9) continue;
            if (trackEMCALCluster->GetIsExotic()) continue;
            numberInAcc++;
          }
        }
      }
      // Cut on number of events in EMCAL acceptance
      if (numberInAcc < fCheckEMCALAcceptanceNumber) continue;
    }

    //histos for impact par studies - D0!!!
    Double_t ptCand = dstarD0pi->Get2Prong()->Pt();
    Double_t invMass = dstarD0pi->InvMassD0();
    Double_t impparXY = dstarD0pi->Get2Prong()->ImpParXY() * 10000.;

    Double_t arrayForSparse[3] = {invMass, ptCand, impparXY};
    Double_t arrayForSparseTrue[3] = {invMass, ptCand, trueImpParXY};

    // set the D0 and D* search window  bin by bin - D* window useful to speed up the reconstruction and D0 window used *ONLY* to calculate side band bkg for the background subtraction methods, for the standard analysis the value in the cut file is considered

    if (0 <= Dstarpt && Dstarpt < 0.5) {
      if (fAnalysis == 1) {
        fD0Window = 0.035;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.020;
        fPeakWindow = 0.0018;
      }
    }
    if (0.5 <= Dstarpt && Dstarpt < 1.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.035;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.020;
        fPeakWindow = 0.0018;
      }
    }
    if (1.0 <= Dstarpt && Dstarpt < 2.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.035;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.020;
        fPeakWindow = 0.0018;
      }
    }
    if (2.0 <= Dstarpt && Dstarpt < 3.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.035;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.022;
        fPeakWindow = 0.0016;
      }
    }
    if (3.0 <= Dstarpt && Dstarpt < 4.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.035;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.026;
        fPeakWindow = 0.0014;
      }
    }
    if (4.0 <= Dstarpt && Dstarpt < 5.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.045;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.026;
        fPeakWindow = 0.0014;
      }
    }
    if (5.0 <= Dstarpt && Dstarpt < 6.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.045;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.026;
        fPeakWindow = 0.006;
      }
    }
    if (6.0 <= Dstarpt && Dstarpt < 7.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.055;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.026;
        fPeakWindow = 0.006;
      }
    }
    if (Dstarpt >= 7.0) {
      if (fAnalysis == 1) {
        fD0Window = 0.074;
        fPeakWindow = 0.03;
      } else {
        fD0Window = 0.026;
        fPeakWindow = 0.006;
      }
    }

    nSelectedProd++;
    nSelectedAna++;

    // check that we are close to signal in the DeltaM - here to save time for PbPb
    Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t invmassDelta = dstarD0pi->DeltaInvMass();

    if (TMath::Abs(invmassDelta - (mPDGDstar - mPDGD0)) > fPeakWindow) continue;
    Int_t isSelected = fCuts->IsSelected(dstarD0pi, AliRDHFCuts::kCandidate, aodEvent); //selected
    if (isSelected > 0) fCEvents->Fill(11);

    // after cuts
    if (fDoImpParDstar && isSelected) {
      fHistMassPtImpParTCDs[0]->Fill(arrayForSparse);
      if (isPrimary) fHistMassPtImpParTCDs[1]->Fill(arrayForSparse);
      else {
        fHistMassPtImpParTCDs[2]->Fill(arrayForSparse);
        fHistMassPtImpParTCDs[3]->Fill(arrayForSparseTrue);
      }
    }

    if (fDoDStarVsY && isSelected) {
      ((TH3F*) (fOutputPID->FindObject("deltamassVsyVsPt")))->Fill(dstarD0pi->DeltaInvMass(), dstarD0pi->YDstar(), dstarD0pi->Pt() );
    }

    // check after cuts
    if (fUseMCInfo) {
      Int_t mcLabel = dstarD0pi->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, mcTrackArray);
      if (mcLabel >= 0) {
        AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(mcLabel));
        if (!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
        Int_t pdgCodeMC = TMath::Abs(mcTrackParticle->GetPdgCode());

        if (pdgCodeMC == 413)
        { //if the track is a DStar we check if it comes from charm

          Double_t ptMC = mcTrackParticle->Pt();

          Bool_t fromCharm = kFALSE;
          Int_t mother = mcTrackParticle->GetMother();
          Int_t istep = 0;
          while (mother >= 0 ) {
            istep++;
            AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(mcTrackArray->At(mother));
            if (mcGranma) {
              Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
              if ((abspdgGranma == 4) || (abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)) fromCharm = kTRUE;
              mother = mcGranma->GetMother();
            } else {
              printf("AliVertexingHFUtils::IsTrackFromCharm: Failed casting the mother particle!");
              break;
            }
          }

          if (fromCharm)
          {

            Bool_t mcPionDStarPresent = kFALSE;
            Bool_t mcPionD0Present = kFALSE;
            Bool_t mcKaonPresent = kFALSE;
            Int_t nDaughterDStar = mcTrackParticle->GetNDaughters();

            if (nDaughterDStar == 2) {
              for (Int_t iDaughterDStar = 0; iDaughterDStar < 2; iDaughterDStar++) {

                AliAODMCParticle* daughterDStar = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(iDaughterDStar));
                if (!daughterDStar) break;
                Int_t pdgCodeDaughterDStar = TMath::Abs(daughterDStar->GetPdgCode());

                if (pdgCodeDaughterDStar == 211) { //if the track is a pion we save its monte carlo label
                  mcPionDStarPresent = kTRUE;

                } else if (pdgCodeDaughterDStar == 421) { //if the track is a D0 we look at its daughters
                  Int_t mcLabelD0 = mcTrackParticle->GetDaughterLabel(iDaughterDStar);
                  Int_t nDaughterD0 = daughterDStar->GetNDaughters();

                  if (nDaughterD0 == 2) {
                    for (Int_t iDaughterD0 = 0; iDaughterD0 < 2; iDaughterD0++) {

                      AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterDStar->GetDaughterLabel(iDaughterD0));
                      if (!daughterD0) break;
                      Int_t pdgCodeDaughterD0 = TMath::Abs(daughterD0->GetPdgCode());

                      if (pdgCodeDaughterD0 == 211) {
                        mcPionD0Present = kTRUE;
                      } else if (pdgCodeDaughterD0 == 321) {
                        mcKaonPresent = kTRUE;
                      } else break;
                    }
                  }
                } else break;
              }
            }

            if (mcPionDStarPresent && mcPionD0Present && mcKaonPresent)
            {
              TString fillthis = "";
              fillthis = "DStarPtTruePostCuts";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
              fillthis = "DStarPtTruePostCutsWeighted";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);
              // fillthis = "PtHardPostCuts";
              // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard);
              // fillthis = "PtHardWeightedPostCuts";
              // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptHard, crossSection);
              // fillthis = "WeightsPostCuts";
              // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(crossSection);
              // fillthis = "TrialsPostCuts";
              // ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->AddBinContent(1,nTrials);
              fillthis = "DStarPtPostCuts";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(Dstarpt);
              fillthis = "DStarPtPostCutsWeighted";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(Dstarpt, crossSection);

              fillthis = "DStar_per_bin_true_PostCuts";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC);
              fillthis = "DStar_per_bin_true_PostCuts_weighted";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(ptMC, crossSection);  

              fillthis = "DStar_per_bin_PostCuts";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(Dstarpt);
              fillthis = "DStar_per_bin_PostCuts_weighted";
              ((TH1F*)(fOutputProductionCheck->FindObject(fillthis)))->Fill(Dstarpt, crossSection);                
            }
          }
        }
      }
    }

    fillthishist = "";
    fillthishist = "PtHardPostCuts";
    ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard);
    fillthishist = "PtHardWeightedPostCuts";
    ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(ptHard, crossSection);
    fillthishist = "WeightsPostCuts";
    ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->Fill(crossSection);
    fillthishist = "TrialsPostCuts";
    ((TH1F*)(fOutputProductionCheck->FindObject(fillthishist)))->AddBinContent(1,nTrials);


    // fill PID
    FillSpectrum(dstarD0pi, isDStar, fCuts, isSelected, fOutputPID, fPIDhist);
    SideBandBackground(dstarD0pi, fCuts, isSelected, fOutputPID, fPIDhist);
    //WrongSignForDStar(dstarD0pi,fCuts,fOutputPID);

    //swich off the PID selection
    fCuts->SetUsePID(kFALSE);
    Int_t isSelectedNoPID = fCuts->IsSelected(dstarD0pi, AliRDHFCuts::kCandidate, aodEvent); //selected
    fCuts->SetUsePID(kTRUE);

    FillSpectrum(dstarD0pi, isDStar, fCuts, isSelectedNoPID, fOutputAll, fAllhist);
    //    SideBandBackground(dstarD0pi,fCuts,isSelectedNoPID, fOutputAll);

    // rare D search ------
    if (fDoSearch) {
      TLorentzVector lorentzTrack1(0, 0, 0, 0); // lorentz 4 vector
      TLorentzVector lorentzTrack2(0, 0, 0, 0); // lorentz 4 vector

      for (Int_t i = 0; i < aodEvent->GetNumberOfTracks(); i++) {

        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
        if (!aodTrack) AliFatal("Not a standard AOD");

        if (dstarD0pi->Charge() == aodTrack->Charge()) continue;
        if ((!(aodTrack->GetStatus()&AliESDtrack::kITSrefit) || (!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        if (TMath::Abs(invmassDelta - (mPDGDstar - mPDGD0)) > 0.02) continue;

        //build the D1 mass
        Double_t mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();

        lorentzTrack1.SetPxPyPzE( dstarD0pi->Px(), dstarD0pi->Py(), dstarD0pi->Pz(), dstarD0pi->E(413) );
        lorentzTrack2.SetPxPyPzE( aodTrack->Px(), aodTrack->Py(), aodTrack->Pz(), aodTrack->E(mass) );

        //D1 mass
        Double_t d1mass = ((lorentzTrack1 + lorentzTrack2).M());
        //mass difference - at 0.4117 and 0.4566
        fDeltaMassD1->Fill(d1mass - dstarD0pi->InvMassDstarKpipi());
      }
    }

    if (isDStar == 1) {
      fTrueDiff2->Fill(dstarD0pi->Pt(), dstarD0pi->DeltaInvMass());
    }

  }

  fCounter->StoreCandidates(aodEvent, nSelectedProd, kTRUE);
  fCounter->StoreCandidates(aodEvent, nSelectedAna, kFALSE);

  delete vHF;

  AliDebug(2, Form("Found %i Reco particles that are D*!!", icountReco));

  PostData(1, fOutput);
  PostData(2, fOutputAll);
  PostData(3, fOutputPID);
  PostData(5, fCounter);
  PostData(6, fOutputProductionCheck);

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::Terminate(Option_t*)
{
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.

  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fDeltaMassD1     = dynamic_cast<TH1F*>(fOutput->FindObject("fDeltaMassD1"));
  fTrueDiff2      = dynamic_cast<TH2F*>(fOutput->FindObject("fTrueDiff2"));

  fOutputAll = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputAll) {
    printf("ERROR: fOutputAll not available\n");
    return;
  }
  fOutputPID = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputPID) {
    printf("ERROR: fOutputPID not available\n");
    return;
  }
  fOutputProductionCheck = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputProductionCheck) {
    printf("ERROR: fOutputProductionCheck not available\n");
    return;
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::UserCreateOutputObjects() {
/// output
  Info("UserCreateOutputObjects", "CreateOutputObjects of task %s\n", GetName());

  //slot #1
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("listAll");

  fOutputPID = new TList();
  fOutputPID->SetOwner();
  fOutputPID->SetName("listPID");

  fOutputProductionCheck = new TList();
  fOutputProductionCheck->SetOwner();
  fOutputProductionCheck->SetName("listPID");


  // define histograms
  DefineHistograms();

  //Counter for Normalization
  fCounter = new AliNormalizationCounter(Form("%s", GetOutputSlot(5)->GetContainer()->GetName()));
  fCounter->Init();

  if (fDoImpParDstar) CreateImpactParameterHistos();

  PostData(1, fOutput);
  PostData(2, fOutputAll);
  PostData(3, fOutputPID);
  PostData(5, fCounter);
  PostData(6, fOutputProductionCheck);

  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSEDStarEMCALProductionCheck::DefineHistograms() {
  /// Create histograms

  fCEvents = new TH1F("fCEvents", "counter", 14, 0, 14);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");
  fCEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
  fCEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
  fCEvents->GetXaxis()->SetBinLabel(3, "good prim vtx and B field");
  fCEvents->GetXaxis()->SetBinLabel(4, "no event selected");
  fCEvents->GetXaxis()->SetBinLabel(5, "no vtx contributors");
  fCEvents->GetXaxis()->SetBinLabel(6, "trigger for PbPb");
  fCEvents->GetXaxis()->SetBinLabel(7, "no z vtx");
  fCEvents->GetXaxis()->SetBinLabel(9, "nEvents Mismatched dAOD");
  fCEvents->GetXaxis()->SetBinLabel(11, "no. of cascade candidates");
  fCEvents->GetXaxis()->SetBinLabel(12, "no. of Dstar after selection cuts");
  fCEvents->GetXaxis()->SetBinLabel(13, "no. of not on-the-fly rec Dstar");
  fCEvents->GetXaxis()->SetBinLabel(14, "no. of Dstar rejected by preselect"); //toadd

  fOutput->Add(fCEvents);

  fTrueDiff2 = new TH2F("DiffDstar_pt", "True Reco diff vs pt", 200, 0, 15, 900, 0, 0.3);
  fOutput->Add(fTrueDiff2);

  fDeltaMassD1 = new TH1F("DeltaMassD1", "delta mass d1", 600, 0, 0.8);
  fOutput->Add(fDeltaMassD1);
  //temp a

  fAllhist = new TH1F*[(fNPtBins + 2) * 18];
  fPIDhist = new TH1F*[(fNPtBins + 2) * 18];

  TString nameMass = " ", nameSgn = " ", nameBkg = " ";

  for (Int_t i = -2; i < fNPtBins; i++) {
    nameMass = "histDeltaMass_";
    nameMass += i + 1;
    nameSgn = "histDeltaSgn_";
    nameSgn += i + 1;
    nameBkg = "histDeltaBkg_";
    nameBkg += i + 1;

    if (i == -2) {
      nameMass = "histDeltaMass";
      nameSgn = "histDeltaSgn";
      nameBkg = "histDeltaBkg";
    }

    TH1F* spectrumMass = new TH1F(nameMass.Data(), "D^{*}-D^{0} invariant mass; #DeltaM [GeV/c^{2}]; Entries", 700, 0.13, 0.2);
    TH1F* spectrumSgn = new TH1F(nameSgn.Data(), "D^{*}-D^{0} Signal invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries", 700, 0.13, 0.2);
    TH1F* spectrumBkg = new TH1F(nameBkg.Data(), "D^{*}-D^{0} Background invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries", 700, 0.13, 0.2);

    nameMass = "histD0Mass_";
    nameMass += i + 1;
    nameSgn = "histD0Sgn_";
    nameSgn += i + 1;
    nameBkg = "histD0Bkg_";
    nameBkg += i + 1;

    if (i == -2) {
      nameMass = "histD0Mass";
      nameSgn = "histD0Sgn";
      nameBkg = "histD0Bkg";
    }

    TH1F* spectrumD0Mass = new TH1F(nameMass.Data(), "D^{0} invariant mass; M(D^{0}) [GeV/c^{2}]; Entries", 200, 1.75, 1.95);
    TH1F* spectrumD0Sgn = new TH1F(nameSgn.Data(), "D^{0} Signal invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries", 200, 1.75, 1.95);
    TH1F* spectrumD0Bkg = new TH1F(nameBkg.Data(), "D^{0} Background invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries", 200, 1.75, 1.95);

    nameMass = "histDstarMass_";
    nameMass += i + 1;
    nameSgn = "histDstarSgn_";
    nameSgn += i + 1;
    nameBkg = "histDstarBkg_";
    nameBkg += i + 1;

    if (i == -2) {
      nameMass = "histDstarMass";
      nameSgn = "histDstarSgn";
      nameBkg = "histDstarBkg";
    }

    TH1F* spectrumDstarMass = new TH1F(nameMass.Data(), "D^{*} invariant mass; M(D^{*}) [GeV/c^{2}]; Entries", 200, 1.9, 2.1);
    TH1F* spectrumDstarSgn = new TH1F(nameSgn.Data(), "D^{*} Signal invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries", 200, 1.9, 2.1);
    TH1F* spectrumDstarBkg = new TH1F(nameBkg.Data(), "D^{*} Background invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries", 200, 1.9, 2.1);

    nameMass = "histSideBandMass_";
    nameMass += i + 1;
    if (i == -2) {
      nameMass = "histSideBandMass";
    }

    TH1F* spectrumSideBandMass = new TH1F(nameMass.Data(), "D^{*}-D^{0} sideband mass; M(D^{*}) [GeV/c^{2}]; Entries", 200, 0.1, 0.2);

    nameMass = "histWrongSignMass_";
    nameMass += i + 1;
    if (i == -2) {
      nameMass = "histWrongSignMass";
    }

    TH1F* spectrumWrongSignMass = new TH1F(nameMass.Data(), "D^{*}-D^{0} wrongsign mass; M(D^{*}) [GeV/c^{2}]; Entries", 200, 0.1, 0.2);


    spectrumMass->Sumw2();
    spectrumSgn->Sumw2();
    spectrumBkg->Sumw2();

    spectrumMass->SetLineColor(6);
    spectrumSgn->SetLineColor(2);
    spectrumBkg->SetLineColor(4);

    spectrumMass->SetMarkerStyle(20);
    spectrumSgn->SetMarkerStyle(20);
    spectrumBkg->SetMarkerStyle(20);
    spectrumMass->SetMarkerSize(0.6);
    spectrumSgn->SetMarkerSize(0.6);
    spectrumBkg->SetMarkerSize(0.6);
    spectrumMass->SetMarkerColor(6);
    spectrumSgn->SetMarkerColor(2);
    spectrumBkg->SetMarkerColor(4);

    spectrumD0Mass->Sumw2();
    spectrumD0Sgn->Sumw2();
    spectrumD0Bkg->Sumw2();

    spectrumD0Mass->SetLineColor(6);
    spectrumD0Sgn->SetLineColor(2);
    spectrumD0Bkg->SetLineColor(4);

    spectrumD0Mass->SetMarkerStyle(20);
    spectrumD0Sgn->SetMarkerStyle(20);
    spectrumD0Bkg->SetMarkerStyle(20);
    spectrumD0Mass->SetMarkerSize(0.6);
    spectrumD0Sgn->SetMarkerSize(0.6);
    spectrumD0Bkg->SetMarkerSize(0.6);
    spectrumD0Mass->SetMarkerColor(6);
    spectrumD0Sgn->SetMarkerColor(2);
    spectrumD0Bkg->SetMarkerColor(4);

    spectrumDstarMass->Sumw2();
    spectrumDstarSgn->Sumw2();
    spectrumDstarBkg->Sumw2();

    spectrumDstarMass->SetLineColor(6);
    spectrumDstarSgn->SetLineColor(2);
    spectrumDstarBkg->SetLineColor(4);

    spectrumDstarMass->SetMarkerStyle(20);
    spectrumDstarSgn->SetMarkerStyle(20);
    spectrumDstarBkg->SetMarkerStyle(20);
    spectrumDstarMass->SetMarkerSize(0.6);
    spectrumDstarSgn->SetMarkerSize(0.6);
    spectrumDstarBkg->SetMarkerSize(0.6);
    spectrumDstarMass->SetMarkerColor(6);
    spectrumDstarSgn->SetMarkerColor(2);
    spectrumDstarBkg->SetMarkerColor(4);

    spectrumSideBandMass->Sumw2();
    spectrumSideBandMass->SetLineColor(4);
    spectrumSideBandMass->SetMarkerStyle(20);
    spectrumSideBandMass->SetMarkerSize(0.6);
    spectrumSideBandMass->SetMarkerColor(4);

    spectrumWrongSignMass->Sumw2();
    spectrumWrongSignMass->SetLineColor(4);
    spectrumWrongSignMass->SetMarkerStyle(20);
    spectrumWrongSignMass->SetMarkerSize(0.6);
    spectrumWrongSignMass->SetMarkerColor(4);

    TH1F* allMass = (TH1F*)spectrumMass->Clone();
    TH1F* allSgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* allBkg  = (TH1F*)spectrumBkg->Clone();

    TH1F* pidMass = (TH1F*)spectrumMass->Clone();
    TH1F* pidSgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* pidBkg  = (TH1F*)spectrumBkg->Clone();

    fOutputAll->Add(allMass);
    fOutputAll->Add(allSgn);
    fOutputAll->Add(allBkg);
    fAllhist[i + 2 + ((fNPtBins + 2)*kDeltaMass)] = allMass;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDeltaSgn)] = allSgn;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDeltaBkg)] = allBkg;

    fOutputPID->Add(pidMass);
    fOutputPID->Add(pidSgn);
    fOutputPID->Add(pidBkg);
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDeltaMass)] = pidMass;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDeltaSgn)] = pidSgn;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDeltaBkg)] = pidBkg;

    TH1F* allD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* allD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* allD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pidD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pidD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pidD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    fOutputAll->Add(allD0Mass);
    fOutputAll->Add(allD0Sgn);
    fOutputAll->Add(allD0Bkg);
    fAllhist[i + 2 + ((fNPtBins + 2)*kDzMass)] = allD0Mass;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDzSgn)] = allD0Sgn;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDzBkg)] = allD0Bkg;

    fOutputPID->Add(pidD0Mass);
    fOutputPID->Add(pidD0Sgn);
    fOutputPID->Add(pidD0Bkg);
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDzMass)] = pidD0Mass;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDzSgn)] = pidD0Sgn;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDzBkg)] = pidD0Bkg;

    TH1F* allDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* allDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* allDstarBkg = (TH1F*)spectrumDstarBkg->Clone();

    TH1F* pidDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pidDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pidDstarBkg = (TH1F*)spectrumDstarBkg->Clone();

    fOutputAll->Add(allDstarMass);
    fOutputAll->Add(allDstarSgn);
    fOutputAll->Add(allDstarBkg);
    fAllhist[i + 2 + ((fNPtBins + 2)*kDstarMass)] = allDstarMass;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDstarSgn)] = allDstarSgn;
    fAllhist[i + 2 + ((fNPtBins + 2)*kDstarBkg)] = allDstarBkg;

    fOutputPID->Add(pidDstarMass);
    fOutputPID->Add(pidDstarSgn);
    fOutputPID->Add(pidDstarBkg);
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDstarMass)] = pidDstarMass;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDstarSgn)] = pidDstarSgn;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kDstarBkg)] = pidDstarBkg;

    TH1F* allSideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pidSideBandMass = (TH1F*)spectrumSideBandMass->Clone();

    fOutputAll->Add(allSideBandMass);
    fOutputPID->Add(pidSideBandMass);
    fAllhist[i + 2 + ((fNPtBins + 2)*kSideBandMass)] = allSideBandMass;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kSideBandMass)] = pidSideBandMass;

    TH1F* allWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pidWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();

    fOutputAll->Add(allWrongSignMass);
    fOutputPID->Add(pidWrongSignMass);
    fAllhist[i + 2 + ((fNPtBins + 2)*kWrongSignMass)] = allWrongSignMass;
    fPIDhist[i + 2 + ((fNPtBins + 2)*kWrongSignMass)] = pidWrongSignMass;

  }

  // pt spectra
  nameMass = "ptMass";
  nameSgn = "ptSgn";
  nameBkg = "ptBkg";

  TH1F* ptspectrumMass = new TH1F(nameMass.Data(), "D^{*} p_{T}; p_{T} [GeV]; Entries", 400, 0, 50);
  TH1F* ptspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal p_{T} - MC; p_{T} [GeV]; Entries", 400, 0, 50);
  TH1F* ptspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background p_{T} - MC; p_{T} [GeV]; Entries", 400, 0, 50);

  ptspectrumMass->Sumw2();
  ptspectrumSgn->Sumw2();
  ptspectrumBkg->Sumw2();

  ptspectrumMass->SetLineColor(6);
  ptspectrumSgn->SetLineColor(2);
  ptspectrumBkg->SetLineColor(4);

  ptspectrumMass->SetMarkerStyle(20);
  ptspectrumSgn->SetMarkerStyle(20);
  ptspectrumBkg->SetMarkerStyle(20);
  ptspectrumMass->SetMarkerSize(0.6);
  ptspectrumSgn->SetMarkerSize(0.6);
  ptspectrumBkg->SetMarkerSize(0.6);
  ptspectrumMass->SetMarkerColor(6);
  ptspectrumSgn->SetMarkerColor(2);
  ptspectrumBkg->SetMarkerColor(4);

  TH1F* ptallMass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptallSgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptallBkg = (TH1F*)ptspectrumBkg->Clone();

  TH1F* ptpidMass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptpidSgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptpidBkg = (TH1F*)ptspectrumBkg->Clone();

  fOutputAll->Add(ptallMass);
  fOutputAll->Add(ptallSgn);
  fOutputAll->Add(ptallBkg);
  fAllhist[((fNPtBins + 2)*kptMass)] = ptallMass;
  fAllhist[((fNPtBins + 2)*kptSgn)] = ptallSgn;
  fAllhist[((fNPtBins + 2)*kptBkg)] = ptallBkg;

  fOutputPID->Add(ptpidMass);
  fOutputPID->Add(ptpidSgn);
  fOutputPID->Add(ptpidBkg);
  fPIDhist[(fNPtBins + 2)*kptMass] = ptpidMass;
  fPIDhist[(fNPtBins + 2)*kptSgn] = ptpidSgn;
  fPIDhist[(fNPtBins + 2)*kptBkg] = ptpidBkg;
  // eta spectra
  nameMass = "etaMass";
  nameSgn = "etaSgn";
  nameBkg = "etaBkg";

  TH1F* etaspectrumMass = new TH1F(nameMass.Data(), "D^{*} #eta; #eta; Entries", 200, -1, 1);
  TH1F* etaspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal #eta - MC; #eta; Entries", 200, -1, 1);
  TH1F* etaspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background #eta - MC; #eta; Entries", 200, -1, 1);

  etaspectrumMass->Sumw2();
  etaspectrumSgn->Sumw2();
  etaspectrumBkg->Sumw2();

  etaspectrumMass->SetLineColor(6);
  etaspectrumSgn->SetLineColor(2);
  etaspectrumBkg->SetLineColor(4);

  etaspectrumMass->SetMarkerStyle(20);
  etaspectrumSgn->SetMarkerStyle(20);
  etaspectrumBkg->SetMarkerStyle(20);
  etaspectrumMass->SetMarkerSize(0.6);
  etaspectrumSgn->SetMarkerSize(0.6);
  etaspectrumBkg->SetMarkerSize(0.6);
  etaspectrumMass->SetMarkerColor(6);
  etaspectrumSgn->SetMarkerColor(2);
  etaspectrumBkg->SetMarkerColor(4);

  TH1F* etaallMass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etaallSgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etaallBkg = (TH1F*)etaspectrumBkg->Clone();

  TH1F* etapidMass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etapidSgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etapidBkg = (TH1F*)etaspectrumBkg->Clone();

  fOutputAll->Add(etaallMass);
  fOutputAll->Add(etaallSgn);
  fOutputAll->Add(etaallBkg);
  fAllhist[(fNPtBins + 2)*ketaMass] = etaallMass;
  fAllhist[(fNPtBins + 2)*ketaSgn] = etaallSgn;
  fAllhist[(fNPtBins + 2)*ketaBkg] = etaallBkg;

  fOutputPID->Add(etapidMass);
  fOutputPID->Add(etapidSgn);
  fOutputPID->Add(etapidBkg);
  fPIDhist[(fNPtBins + 2)*ketaMass] = etapidMass;
  fPIDhist[(fNPtBins + 2)*ketaSgn] = etapidSgn;
  fPIDhist[(fNPtBins + 2)*ketaBkg] = etapidBkg;

  if (fDoDStarVsY) {
    TH3F* deltamassVsyVsPtPID = new TH3F("deltamassVsyVsPt", "delta mass Vs y Vs pT;  #DeltaM [GeV/c^{2}]; y; p_{T} [GeV/c]", 700, 0.13, 0.2, 40, -1, 1, 36, 0., 36.);
    fOutputPID->Add(deltamassVsyVsPtPID);
  }

  TString name_DStarPtTruePreEventSelection = "DStarPtTruePreEventSelection";
  TH1F* hist_DStarPtTruePreEventSelection = new TH1F(name_DStarPtTruePreEventSelection.Data(), "DStarPtTruePreEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePreEventSelection->Sumw2();
  hist_DStarPtTruePreEventSelection->SetLineColor(6);
  hist_DStarPtTruePreEventSelection->SetMarkerStyle(20);
  hist_DStarPtTruePreEventSelection->SetMarkerSize(0.6);
  hist_DStarPtTruePreEventSelection->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePreEventSelection = (TH1F*)hist_DStarPtTruePreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePreEventSelection);

  TString name_DStarPtTruePreEventSelectionWeighted = "DStarPtTruePreEventSelectionWeighted";
  TH1F* hist_DStarPtTruePreEventSelectionWeighted = new TH1F(name_DStarPtTruePreEventSelectionWeighted.Data(), "DStarPtTruePreEventSelectionWeighted; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePreEventSelectionWeighted->Sumw2();
  hist_DStarPtTruePreEventSelectionWeighted->SetLineColor(6);
  hist_DStarPtTruePreEventSelectionWeighted->SetMarkerStyle(20);
  hist_DStarPtTruePreEventSelectionWeighted->SetMarkerSize(0.6);
  hist_DStarPtTruePreEventSelectionWeighted->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePreEventSelectionWeighted = (TH1F*)hist_DStarPtTruePreEventSelectionWeighted->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePreEventSelectionWeighted);

  TString name_DStarPtTruePostEventSelection = "DStarPtTruePostEventSelection";
  TH1F* hist_DStarPtTruePostEventSelection = new TH1F(name_DStarPtTruePostEventSelection.Data(), "DStarPtTruePostEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePostEventSelection->Sumw2();
  hist_DStarPtTruePostEventSelection->SetLineColor(6);
  hist_DStarPtTruePostEventSelection->SetMarkerStyle(20);
  hist_DStarPtTruePostEventSelection->SetMarkerSize(0.6);
  hist_DStarPtTruePostEventSelection->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePostEventSelection = (TH1F*)hist_DStarPtTruePostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePostEventSelection);

  TString name_DStarPtTruePostEventSelectionWeighted = "DStarPtTruePostEventSelectionWeighted";
  TH1F* hist_DStarPtTruePostEventSelectionWeighted = new TH1F(name_DStarPtTruePostEventSelectionWeighted.Data(), "DStarPtTruePostEventSelectionWeighted; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePostEventSelectionWeighted->Sumw2();
  hist_DStarPtTruePostEventSelectionWeighted->SetLineColor(6);
  hist_DStarPtTruePostEventSelectionWeighted->SetMarkerStyle(20);
  hist_DStarPtTruePostEventSelectionWeighted->SetMarkerSize(0.6);
  hist_DStarPtTruePostEventSelectionWeighted->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePostEventSelectionWeighted = (TH1F*)hist_DStarPtTruePostEventSelectionWeighted->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePostEventSelectionWeighted);

  TString name_DStarPtTruePostCuts = "DStarPtTruePostCuts";
  TH1F* hist_DStarPtTruePostCuts = new TH1F(name_DStarPtTruePostCuts.Data(), "DStarPtTruePostCuts; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePostCuts->Sumw2();
  hist_DStarPtTruePostCuts->SetLineColor(6);
  hist_DStarPtTruePostCuts->SetMarkerStyle(20);
  hist_DStarPtTruePostCuts->SetMarkerSize(0.6);
  hist_DStarPtTruePostCuts->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePostCuts = (TH1F*)hist_DStarPtTruePostCuts->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePostCuts);

  TString name_DStarPtTruePostCutsWeighted = "DStarPtTruePostCutsWeighted";
  TH1F* hist_DStarPtTruePostCutsWeighted = new TH1F(name_DStarPtTruePostCutsWeighted.Data(), "DStarPtTruePostCutsWeighted; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtTruePostCutsWeighted->Sumw2();
  hist_DStarPtTruePostCutsWeighted->SetLineColor(6);
  hist_DStarPtTruePostCutsWeighted->SetMarkerStyle(20);
  hist_DStarPtTruePostCutsWeighted->SetMarkerSize(0.6);
  hist_DStarPtTruePostCutsWeighted->SetMarkerColor(6);
  TH1F* histogram_DStarPtTruePostCutsWeighted = (TH1F*)hist_DStarPtTruePostCutsWeighted->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtTruePostCutsWeighted);

  TString name_DStarPtPostCuts = "DStarPtPostCuts";
  TH1F* hist_DStarPtPostCuts = new TH1F(name_DStarPtPostCuts.Data(), "DStarPtPostCuts; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtPostCuts->Sumw2();
  hist_DStarPtPostCuts->SetLineColor(6);
  hist_DStarPtPostCuts->SetMarkerStyle(20);
  hist_DStarPtPostCuts->SetMarkerSize(0.6);
  hist_DStarPtPostCuts->SetMarkerColor(6);
  TH1F* histogram_DStarPtPostCuts = (TH1F*)hist_DStarPtPostCuts->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtPostCuts);

  TString name_DStarPtPostCutsWeighted = "DStarPtPostCutsWeighted";
  TH1F* hist_DStarPtPostCutsWeighted = new TH1F(name_DStarPtPostCutsWeighted.Data(), "DStarPtPostCutsWeighted; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_DStarPtPostCutsWeighted->Sumw2();
  hist_DStarPtPostCutsWeighted->SetLineColor(6);
  hist_DStarPtPostCutsWeighted->SetMarkerStyle(20);
  hist_DStarPtPostCutsWeighted->SetMarkerSize(0.6);
  hist_DStarPtPostCutsWeighted->SetMarkerColor(6);
  TH1F* histogram_DStarPtPostCutsWeighted = (TH1F*)hist_DStarPtPostCutsWeighted->Clone();
  fOutputProductionCheck->Add(histogram_DStarPtPostCutsWeighted);

  TString name_PtHardPreEventSelection = "PtHardPreEventSelection";
  TH1F* hist_PtHardPreEventSelection = new TH1F(name_PtHardPreEventSelection.Data(), "PtHardPreEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardPreEventSelection->Sumw2();
  hist_PtHardPreEventSelection->SetLineColor(6);
  hist_PtHardPreEventSelection->SetMarkerStyle(20);
  hist_PtHardPreEventSelection->SetMarkerSize(0.6);
  hist_PtHardPreEventSelection->SetMarkerColor(6);
  TH1F* histogram_PtHardPreEventSelection = (TH1F*)hist_PtHardPreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_PtHardPreEventSelection);

  TString name_PtHardPostEventSelection = "PtHardPostEventSelection";
  TH1F* hist_PtHardPostEventSelection = new TH1F(name_PtHardPostEventSelection.Data(), "PtHardPostEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardPostEventSelection->Sumw2();
  hist_PtHardPostEventSelection->SetLineColor(6);
  hist_PtHardPostEventSelection->SetMarkerStyle(20);
  hist_PtHardPostEventSelection->SetMarkerSize(0.6);
  hist_PtHardPostEventSelection->SetMarkerColor(6);
  TH1F* histogram_PtHardPostEventSelection = (TH1F*)hist_PtHardPostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_PtHardPostEventSelection);

  TString name_PtHardPostCuts = "PtHardPostCuts";
  TH1F* hist_PtHardPostCuts = new TH1F(name_PtHardPostCuts.Data(), "PtHardPostCuts; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardPostCuts->Sumw2();
  hist_PtHardPostCuts->SetLineColor(6);
  hist_PtHardPostCuts->SetMarkerStyle(20);
  hist_PtHardPostCuts->SetMarkerSize(0.6);
  hist_PtHardPostCuts->SetMarkerColor(6);
  TH1F* histogram_PtHardPostCuts = (TH1F*)hist_PtHardPostCuts->Clone();
  fOutputProductionCheck->Add(histogram_PtHardPostCuts);

  TString name_PtHardWeightedPreEventSelection = "PtHardWeightedPreEventSelection";
  TH1F* hist_PtHardWeightedPreEventSelection = new TH1F(name_PtHardWeightedPreEventSelection.Data(), "PtHardWeightedPreEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardWeightedPreEventSelection->Sumw2();
  hist_PtHardWeightedPreEventSelection->SetLineColor(6);
  hist_PtHardWeightedPreEventSelection->SetMarkerStyle(20);
  hist_PtHardWeightedPreEventSelection->SetMarkerSize(0.6);
  hist_PtHardWeightedPreEventSelection->SetMarkerColor(6);
  TH1F* histogram_PtHardWeightedPreEventSelection = (TH1F*)hist_PtHardWeightedPreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_PtHardWeightedPreEventSelection);

  TString name_PtHardWeightedPostEventSelection = "PtHardWeightedPostEventSelection";
  TH1F* hist_PtHardWeightedPostEventSelection = new TH1F(name_PtHardWeightedPostEventSelection.Data(), "PtHardWeightedPostEventSelection; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardWeightedPostEventSelection->Sumw2();
  hist_PtHardWeightedPostEventSelection->SetLineColor(6);
  hist_PtHardWeightedPostEventSelection->SetMarkerStyle(20);
  hist_PtHardWeightedPostEventSelection->SetMarkerSize(0.6);
  hist_PtHardWeightedPostEventSelection->SetMarkerColor(6);
  TH1F* histogram_PtHardWeightedPostEventSelection = (TH1F*)hist_PtHardWeightedPostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_PtHardWeightedPostEventSelection);

  TString name_PtHardWeightedPostCuts = "PtHardWeightedPostCuts";
  TH1F* hist_PtHardWeightedPostCuts = new TH1F(name_PtHardWeightedPostCuts.Data(), "PtHardWeightedPostCuts; p_{T} [GeV/c]; Entries", 5000, 0, 1000);
  hist_PtHardWeightedPostCuts->Sumw2();
  hist_PtHardWeightedPostCuts->SetLineColor(6);
  hist_PtHardWeightedPostCuts->SetMarkerStyle(20);
  hist_PtHardWeightedPostCuts->SetMarkerSize(0.6);
  hist_PtHardWeightedPostCuts->SetMarkerColor(6);
  TH1F* histogram_PtHardWeightedPostCuts = (TH1F*)hist_PtHardWeightedPostCuts->Clone();
  fOutputProductionCheck->Add(histogram_PtHardWeightedPostCuts);

  TString name_WeightsPreEventSelection = "WeightsPreEventSelection";
  TH1F* hist_WeightsPreEventSelection = new TH1F(name_WeightsPreEventSelection.Data(), "WeightsPreEventSelection; p_{T} [GeV/c]; Entries", 50000, 0, 1000);
  hist_WeightsPreEventSelection->Sumw2();
  hist_WeightsPreEventSelection->SetLineColor(6);
  hist_WeightsPreEventSelection->SetMarkerStyle(20);
  hist_WeightsPreEventSelection->SetMarkerSize(0.6);
  hist_WeightsPreEventSelection->SetMarkerColor(6);
  TH1F* histogram_WeightsPreEventSelection = (TH1F*)hist_WeightsPreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_WeightsPreEventSelection);

  TString name_WeightsPostEventSelection = "WeightsPostEventSelection";
  TH1F* hist_WeightsPostEventSelection = new TH1F(name_WeightsPostEventSelection.Data(), "WeightsPostEventSelection; p_{T} [GeV/c]; Entries", 50000, 0, 1000);
  hist_WeightsPostEventSelection->Sumw2();
  hist_WeightsPostEventSelection->SetLineColor(6);
  hist_WeightsPostEventSelection->SetMarkerStyle(20);
  hist_WeightsPostEventSelection->SetMarkerSize(0.6);
  hist_WeightsPostEventSelection->SetMarkerColor(6);
  TH1F* histogram_WeightsPostEventSelection = (TH1F*)hist_WeightsPostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_WeightsPostEventSelection);

  TString name_WeightsPostCuts = "WeightsPostCuts";
  TH1F* hist_WeightsPostCuts = new TH1F(name_WeightsPostCuts.Data(), "WeightsPostCuts; p_{T} [GeV/c]; Entries", 50000, 0, 1000);
  hist_WeightsPostCuts->Sumw2();
  hist_WeightsPostCuts->SetLineColor(6);
  hist_WeightsPostCuts->SetMarkerStyle(20);
  hist_WeightsPostCuts->SetMarkerSize(0.6);
  hist_WeightsPostCuts->SetMarkerColor(6);
  TH1F* histogram_WeightsPostCuts = (TH1F*)hist_WeightsPostCuts->Clone();
  fOutputProductionCheck->Add(histogram_WeightsPostCuts);

  TString name_TrialsPreEventSelection = "TrialsPreEventSelection";
  TH1F* hist_TrialsPreEventSelection = new TH1F(name_TrialsPreEventSelection.Data(), "TrialsPreEventSelection; p_{T} [GeV/c]; Entries", 1, 0, 1);
  hist_TrialsPreEventSelection->Sumw2();
  hist_TrialsPreEventSelection->SetLineColor(6);
  hist_TrialsPreEventSelection->SetMarkerStyle(20);
  hist_TrialsPreEventSelection->SetMarkerSize(0.6);
  hist_TrialsPreEventSelection->SetMarkerColor(6);
  TH1F* histogram_TrialsPreEventSelection = (TH1F*)hist_TrialsPreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_TrialsPreEventSelection);

  TString name_TrialsPostEventSelection = "TrialsPostEventSelection";
  TH1F* hist_TrialsPostEventSelection = new TH1F(name_TrialsPostEventSelection.Data(), "TrialsPostEventSelection; p_{T} [GeV/c]; Entries", 1, 0, 1);
  hist_TrialsPostEventSelection->Sumw2();
  hist_TrialsPostEventSelection->SetLineColor(6);
  hist_TrialsPostEventSelection->SetMarkerStyle(20);
  hist_TrialsPostEventSelection->SetMarkerSize(0.6);
  hist_TrialsPostEventSelection->SetMarkerColor(6);
  TH1F* histogram_TrialsPostEventSelection = (TH1F*)hist_TrialsPostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_TrialsPostEventSelection);

  TString name_TrialsPostCuts = "TrialsPostCuts";
  TH1F* hist_TrialsPostCuts = new TH1F(name_TrialsPostCuts.Data(), "TrialsPostCuts; p_{T} [GeV/c]; Entries", 1, 0, 1);
  hist_TrialsPostCuts->Sumw2();
  hist_TrialsPostCuts->SetLineColor(6);
  hist_TrialsPostCuts->SetMarkerStyle(20);
  hist_TrialsPostCuts->SetMarkerSize(0.6);
  hist_TrialsPostCuts->SetMarkerColor(6);
  TH1F* histogram_TrialsPostCuts = (TH1F*)hist_TrialsPostCuts->Clone();
  fOutputProductionCheck->Add(histogram_TrialsPostCuts);

  Int_t nPtBins = fCuts->GetNPtBins();
  // const Int_t nPtBinLimits = nPtBins + 1;
  Float_t * PtBinLimits = fCuts->GetPtBinLimits();

  TString name_DStar_per_bin_true_PreEventSelection ="DStar_per_bin_true_PreEventSelection";
  TH1F* hist_DStar_per_bin_true_PreEventSelection = new TH1F(name_DStar_per_bin_true_PreEventSelection.Data(),"DStar_per_bin_true_PreEventSelection; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PreEventSelection = (TH1F*)hist_DStar_per_bin_true_PreEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PreEventSelection);

  TString name_DStar_per_bin_true_PreEventSelection_weighted ="DStar_per_bin_true_PreEventSelection_weighted";
  TH1F* hist_DStar_per_bin_true_PreEventSelection_weighted = new TH1F(name_DStar_per_bin_true_PreEventSelection_weighted.Data(),"DStar_per_bin_true_PreEventSelection_weighted; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PreEventSelection_weighted = (TH1F*)hist_DStar_per_bin_true_PreEventSelection_weighted->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PreEventSelection_weighted);

  TString name_DStar_per_bin_true_PostEventSelection ="DStar_per_bin_true_PostEventSelection";
  TH1F* hist_DStar_per_bin_true_PostEventSelection = new TH1F(name_DStar_per_bin_true_PostEventSelection.Data(),"DStar_per_bin_true_PostEventSelection; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PostEventSelection = (TH1F*)hist_DStar_per_bin_true_PostEventSelection->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PostEventSelection);

  TString name_DStar_per_bin_true_PostEventSelection_weighted ="DStar_per_bin_true_PostEventSelection_weighted";
  TH1F* hist_DStar_per_bin_true_PostEventSelection_weighted = new TH1F(name_DStar_per_bin_true_PostEventSelection_weighted.Data(),"DStar_per_bin_true_PostEventSelection_weighted; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PostEventSelection_weighted = (TH1F*)hist_DStar_per_bin_true_PostEventSelection_weighted->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PostEventSelection_weighted);

  TString name_DStar_per_bin_true_PostCuts ="DStar_per_bin_true_PostCuts";
  TH1F* hist_DStar_per_bin_true_PostCuts = new TH1F(name_DStar_per_bin_true_PostCuts.Data(),"DStar_per_bin_true_PostCuts; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PostCuts = (TH1F*)hist_DStar_per_bin_true_PostCuts->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PostCuts);

  TString name_DStar_per_bin_true_PostCuts_weighted ="DStar_per_bin_true_PostCuts_weighted";
  TH1F* hist_DStar_per_bin_true_PostCuts_weighted = new TH1F(name_DStar_per_bin_true_PostCuts_weighted.Data(),"DStar_per_bin_true_PostCuts_weighted; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_true_PostCuts_weighted = (TH1F*)hist_DStar_per_bin_true_PostCuts_weighted->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_true_PostCuts_weighted);

  TString name_DStar_per_bin_PostCuts ="DStar_per_bin_PostCuts";
  TH1F* hist_DStar_per_bin_PostCuts = new TH1F(name_DStar_per_bin_PostCuts.Data(),"DStar_per_bin_PostCuts; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_PostCuts = (TH1F*)hist_DStar_per_bin_PostCuts->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_PostCuts);

  TString name_DStar_per_bin_PostCuts_weighted ="DStar_per_bin_PostCuts_weighted";
  TH1F* hist_DStar_per_bin_PostCuts_weighted = new TH1F(name_DStar_per_bin_PostCuts_weighted.Data(),"DStar_per_bin_PostCuts_weighted; Entries",nPtBins,PtBinLimits); 
  TH1F* histogram_DStar_per_bin_PostCuts_weighted = (TH1F*)hist_DStar_per_bin_PostCuts_weighted->Clone();
  fOutputProductionCheck->Add(histogram_DStar_per_bin_PostCuts_weighted);

  TString name_fHistClusPosition ="fHistClusPosition";
  TH3F* hist_fHistClusPosition = new TH3F(name_fHistClusPosition.Data(),";#it{x} (cm);#it{y} (cm);#it{z} (cm)", 50, -500, 500, 50, -500, 500, 50, -500, 500);
  TH3F* histogram_fHistClusPosition = (TH3F*)hist_fHistClusPosition->Clone();
  fOutputProductionCheck->Add(histogram_fHistClusPosition);


  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist) {
  //
  /// Fill histos for D* spectrum
  //

  if (!isSel) return;

  // D0 window
  Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();


  Int_t ptbin = cuts->PtBin(part->Pt());
  Double_t pt = part->Pt();
  Double_t eta = part->Eta();

  Double_t invmassDelta = part->DeltaInvMass();
  Double_t invmassDstar = part->InvMassDstarKpipi();

  TString fillthis = "";
  Bool_t massInRange = kFALSE;

  Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();

  // delta M(Kpipi)-M(Kpi)
  if (TMath::Abs(invmassDelta - (mPDGDstar - mPDGD0)) < fPeakWindow) massInRange = kTRUE;

  if (fUseMCInfo) {
    if (isDStar == 1) {
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDzSgn)]->Fill(invmassD0);
      histlist[(fNPtBins + 2)*kDzSgn]->Fill(invmassD0);
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDstarSgn)]->Fill(invmassDstar);
      histlist[(fNPtBins + 2)*kDstarSgn]->Fill(invmassDstar);
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDeltaSgn)]->Fill(invmassDelta);
      histlist[(fNPtBins + 2)*kDeltaSgn]->Fill(invmassDelta);
      if (massInRange) {
        histlist[(fNPtBins + 2)*kptSgn]->Fill(pt);
        histlist[(fNPtBins + 2)*ketaSgn]->Fill(eta);
      }
    }
    else {//background
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDzBkg)]->Fill(invmassD0);
      histlist[(fNPtBins + 2)*kDzBkg]->Fill(invmassD0);
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDstarBkg)]->Fill(invmassDstar);
      histlist[(fNPtBins + 2)*kDstarBkg]->Fill(invmassDstar);
      histlist[ptbin + 1 + ((fNPtBins + 2)*kDeltaBkg)]->Fill(invmassDelta);
      histlist[(fNPtBins + 2)*kDeltaBkg]->Fill(invmassDelta);
      if (massInRange) {
        histlist[(fNPtBins + 2)*kptBkg]->Fill(pt);
        histlist[(fNPtBins + 2)*ketaBkg]->Fill(eta);
      }
    }
  }
  //no MC info, just cut selection
  histlist[ptbin + 1 + ((fNPtBins + 2)*kDzMass)]->Fill(invmassD0);
  histlist[(fNPtBins + 2)*kDzMass]->Fill(invmassD0);
  histlist[ptbin + 1 + ((fNPtBins + 2)*kDstarMass)]->Fill(invmassDstar);
  histlist[(fNPtBins + 2)*kDstarMass]->Fill(invmassDstar);
  histlist[ptbin + 1 + ((fNPtBins + 2)*kDeltaMass)]->Fill(invmassDelta);
  histlist[(fNPtBins + 2)*kDeltaMass]->Fill(invmassDelta);

  if (massInRange) {
    histlist[(fNPtBins + 2)*kptMass]->Fill(pt);
    histlist[(fNPtBins + 2)*ketaMass]->Fill(eta);
  }

  return;
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::SideBandBackground(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist) {

  ///  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas
  /// (expected detector resolution) on the left and right frm the D0 mass. Each band
  ///  has a width of ~5 sigmas. Two band needed  for opening angle considerations

  if (!isSel) return;

  Int_t ptbin = cuts->PtBin(part->Pt());

  // select the side bands intervall
  Double_t invmassD0    = part->InvMassD0();
  if (TMath::Abs(invmassD0 - 1.865) > 4 * fD0Window && TMath::Abs(invmassD0 - 1.865) < 8 * fD0Window) {

    // for pt and eta
    Double_t invmassDelta = part->DeltaInvMass();

    histlist[ptbin + 1 + ((fNPtBins + 2)*kSideBandMass)]->Fill(invmassDelta);
    histlist[(fNPtBins + 2)*kSideBandMass]->Fill(invmassDelta);


  }
}
//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarEMCALProductionCheck::WrongSignForDStar(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, TList *listout) {
  //
  /// assign the wrong charge to the soft pion to create background
  //
  Int_t ptbin = cuts->PtBin(part->Pt());

  Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();
  if (TMath::Abs(invmassD0 - mPDGD0) > fD0Window) return;

  AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)part->Get2Prong();

  Int_t okDzWrongSign;
  Double_t wrongMassD0 = 0.;

  Int_t isSelected = cuts->IsSelected(part, AliRDHFCuts::kCandidate); //selected
  if (!isSelected) {
    return;
  }

  okDzWrongSign =  1;

  //if is D*+ than assume D0bar
  if (part->Charge() > 0 && (isSelected == 1)) {
    okDzWrongSign = 0;
  }

  // assign the wrong mass in case the cuts return both D0 and D0bar
  if (part->Charge() > 0 && (isSelected == 3)) {
    okDzWrongSign = 0;
  }

  //wrong D0 inv mass
  if (okDzWrongSign != 0) {
    wrongMassD0 = theD0particle->InvMassD0();
  } else if (okDzWrongSign == 0) {
    wrongMassD0 = theD0particle->InvMassD0bar();
  }

  if (TMath::Abs(wrongMassD0 - 1.865) < fD0Window) {

    // wrong D* inv mass
    Double_t e[3];
    if (part->Charge() > 0) {
      e[0] = theD0particle->EProng(0, 321);
      e[1] = theD0particle->EProng(1, 211);
    } else {
      e[0] = theD0particle->EProng(0, 211);
      e[1] = theD0particle->EProng(1, 321);
    }
    e[2] = part->EProng(0, 211);

    Double_t esum = e[0] + e[1] + e[2];
    Double_t pds = part->P();

    Double_t   wrongMassDstar = TMath::Sqrt(esum * esum - pds * pds);

    TString fillthis = "";
    fillthis = "histWrongSignMass_";
    fillthis += ptbin;
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar - wrongMassD0);
    fillthis = "histWrongSignMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar - wrongMassD0);

  }
}

//-------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEDStarEMCALProductionCheck::CheckOrigin(TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const {
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma = 0;
  Bool_t isFromB = kFALSE;
  while (mother > 0 ) {
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma) {
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
        isFromB = kTRUE;
      }
      mother = mcGranma->GetMother();
    } else {
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  if (isFromB) return 5;
  else return 4;
}
//-------------------------------------------------------------------------------------
Float_t AliAnalysisTaskSEDStarEMCALProductionCheck::GetTrueImpactParameterD0(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const {
  /// true impact parameter calculation

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partDp->XvYvZv(origD);
  Short_t charge = partDp->Charge();
  Double_t pXdauTrue[3], pYdauTrue[3], pZdauTrue[3];
  Int_t labelFirstDau = partDp->GetDaughterLabel(0);

  Int_t nDau = partDp->GetNDaughters();

  Int_t theDau = 0;
  if (nDau == 2) {
    for (Int_t iDau = 0; iDau < 2; iDau++) {
      Int_t ind = labelFirstDau + iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if (!part) {
        AliError("Daughter particle not found in MC array");
        return 99999.;
      }
      Int_t pdgCode = TMath::Abs(part->GetPdgCode());
      if (pdgCode == 211 || pdgCode == 321) {
        pXdauTrue[theDau] = part->Px();
        pYdauTrue[theDau] = part->Py();
        pZdauTrue[theDau] = part->Pz();
        ++theDau;
      }
    }
  }
  if (theDau != 2) {
    AliError("Wrong number of decay prongs");
    return 99999.;
  }

  Double_t d0dummy[3] = {0., 0., 0.};
  AliAODRecoDecayHF aodD0MC(vtxTrue, origD, 3, charge, pXdauTrue, pYdauTrue, pZdauTrue, d0dummy);
  return aodD0MC.ImpParXY();

}
//______________________________________________________-
void AliAnalysisTaskSEDStarEMCALProductionCheck::CreateImpactParameterHistos() {
  /// Histos for impact paramter study

  Int_t nbins[3] = {400, 200, fNImpParBins};
  Double_t xmin[3] = {1.75, 0., fLowerImpPar};
  Double_t xmax[3] = {1.98, 20., fHigherImpPar};

  fHistMassPtImpParTCDs[0] = new THnSparseF("hMassPtImpParAll",
      "Mass vs. pt vs.imppar - All",
      3, nbins, xmin, xmax);
  fHistMassPtImpParTCDs[1] = new THnSparseF("hMassPtImpParPrompt",
      "Mass vs. pt vs.imppar - promptD",
      3, nbins, xmin, xmax);
  fHistMassPtImpParTCDs[2] = new THnSparseF("hMassPtImpParBfeed",
      "Mass vs. pt vs.imppar - DfromB",
      3, nbins, xmin, xmax);
  fHistMassPtImpParTCDs[3] = new THnSparseF("hMassPtImpParTrueBfeed",
      "Mass vs. pt vs.true imppar -DfromB",
      3, nbins, xmin, xmax);
  fHistMassPtImpParTCDs[4] = new THnSparseF("hMassPtImpParBkg",
      "Mass vs. pt vs.imppar - backgr.",
      3, nbins, xmin, xmax);

  for (Int_t i = 0; i < 5; i++) {
    fOutput->Add(fHistMassPtImpParTCDs[i]);
  }
}

