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

//-------------------------------------------------------------------------
// 1) Analysis task to identify the jet whose cluster fired the trigger patch
// 2) perform some QA on the patch / cluster / jet
// 3) and pass the saved out collection of jet(s) to other tasks
//
// currently set up for GA trigger
//
// Author: Joel Mazer (joel.mazer@cern.ch)
//-------------------------------------------------------------------------

// task head include
//#include "AliAnalysisTaskEmcalTriggerPatchJetMatch.h"
#include <iostream>

// ROOT includes
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

// event handler (and pico's) includes
#include "AliAnalysisManager.h"
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliESDInputHandler.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVVZERO.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalTriggerPatchInfoAP.h"
#include "AliAODHeader.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEmcalTriggerPatchJetMatch.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalTriggerPatchJetMatch)

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerPatchJetMatch::AliAnalysisTaskEmcalTriggerPatchJetMatch() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalTriggerPatchJetMatch", kTRUE),
  fDebug(kFALSE), fAttachToEvent(kTRUE),
  fTriggerClass(""),
  fJetContainer(0),
  fMaxPatchEnergy(0), fMaxPatchADCEnergy(0),
  fTriggerType(1), //-1),
  fNFastOR(16),
  fMainTrigCat(kTriggerLevel1Jet), // options: kTriggerLevel0, kTriggerLevel1Jet, kTriggerLevel1Gamma, kTriggerRecalcJet, kTriggerRecalGamma
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),// Recalculated max trigger patch; does not need to be above trigger threshold
  fTriggerCategory(kTriggerRecalcGamma),
  fMainTrigSimple(kFALSE), // (kTRUE)
  fJetPtCut(10.0), // jet pt cut at 10 by default
  fPatchECut(10.0), // patch E cut at 10 by default
  fTrkBias(0), fClusBias(0),
  doComments(0), fUseALLrecalcPatches(0), // defaults to max patch only
  fJetTriggeredEventname("JETthatTriggeredEvent"), fCaloClustersName(""),
  fMaxPatch(0), // max patch object
  fMainPatchType(kManual), // (kEmcalJet)
  fhNEvents(0),
  fhTriggerbit(0), 
  fHistRhovsCent(0),
  fh3PtEtaPhiTracks(0), fh3PtEtaPhiTracksOnEmcal(0), fh3PtEtaPhiTracksToProp(0), fh3PtEtaPhiTracksProp(0), fh3PtEtaPhiTracksNoProp(0),
  fh2CentPtJet(0),
  fh3PtEtaPhiJet(0),
  fh2NJetsPt(0),
  fh3PtEtaAreaJet(0),
  fh2PtNConstituentsCharged(0), fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0), fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0), fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0), fh3PtLeadJet2VsPatchEnergy(0),
  fh3PtLeadJet1PatchEnergyVZEROAmp(0), fh3PtLeadJet1RawPatchEnergyVZEROAmp(0),
  fh3PatchEnergyEtaPhiCenterJ1(0), fh3PatchEnergyEtaPhiCenterJ2(0), fh3PatchEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1(0), fh3PatchADCEnergyEtaPhiCenterJ2(0), fh3PatchADCEnergyEtaPhiCenterJ1J2(0),
  fh3PatchEnergyEtaPhiCenterG1(0), fh3PatchEnergyEtaPhiCenterG2(0), fh3PatchEnergyEtaPhiCenterG1G2(0),
  fh3PatchADCEnergyEtaPhiCenterG1(0), fh3PatchADCEnergyEtaPhiCenterG2(0), fh3PatchADCEnergyEtaPhiCenterG1G2(0),
  fh3PatchADCEnergyEtaPhiCenterAll(0),
  fh3EEtaPhiCell(0), fh2ECellVsCent(0), fh2CellEnergyVsTime(0), fh3EClusELeadingCellVsTime(0),
  fh3JetReacCent(0),
  fHistClusEnergy(0), fHistClusofJetEnergy(0),
  fHistEventSelectionQA(0), fhQAinfoAllPatchesCounter(0), fhQAinfoCounter(0), fhQAmaxinfoCounter(0),
  fhRecalcGammaPatchEnergy(0), fhGammaLowPatchEnergy(0), fhGammaLowSimplePatchEnergy(0),
  fhRecalcJetPatchEnergy(0), fhJetLowPatchEnergy(0), fhJetLowSimplePatchEnergy(0),
  fhMainTriggerPatchEnergy(0),
  fJetTriggeredEvent(0), fRecalcTriggerPatches(0),
  fhnPatchMaxClus(0x0), fhnPatchMatch(0x0), fhnPatchMatch2(0x0), fhnPatchMatchJetLeadClus(0x0)
{
  // Default constructor.
  for(Int_t j=0; j<16; j++) {
    fHistdPhidEtaPatchJetCluster[j] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerPatchJetMatch::AliAnalysisTaskEmcalTriggerPatchJetMatch(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(kFALSE), fAttachToEvent(kTRUE),
  fTriggerClass(""),
  fJetContainer(0),
  fMaxPatchEnergy(0), fMaxPatchADCEnergy(0),
  fTriggerType(1), //-1),
  fNFastOR(16),
  fMainTrigCat(kTriggerLevel1Jet), // options: kTriggerLevel0, kTriggerLevel1Jet, kTriggerLevel1Gamma, kTriggerRecalcJet, kTriggerRecalGamma
  //fMainTrigCat(kTriggerRecalcJet),  // Recalculated max trigger patch; does not need to be above trigger threshold
  //fMainTrigCat(kTriggerRecalcGamma),// Recalculated max trigger patch; does not need to be above trigger threshold
  fTriggerCategory(kTriggerRecalcGamma),
  fMainTrigSimple(kFALSE), // (kTRUE)
  fJetPtCut(10.0), // jet pt cut at 10 by default
  fPatchECut(10.0), // patch E cut at 10 by default
  fTrkBias(0), fClusBias(0),
  doComments(0), fUseALLrecalcPatches(0), // defaults to max patch only
  fJetTriggeredEventname("JETthatTriggeredEvent"), fCaloClustersName(""),
  fMaxPatch(0),
  fMainPatchType(kManual), //kEmcalJet
  fhNEvents(0),
  fhTriggerbit(0),
  fHistRhovsCent(0),
  fh3PtEtaPhiTracks(0), fh3PtEtaPhiTracksOnEmcal(0), fh3PtEtaPhiTracksToProp(0), fh3PtEtaPhiTracksProp(0), fh3PtEtaPhiTracksNoProp(0),
  fh2CentPtJet(0),
  fh3PtEtaPhiJet(0),
  fh2NJetsPt(0),
  fh3PtEtaAreaJet(0),
  fh2PtNConstituentsCharged(0), fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0), fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0), fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0), fh3PtLeadJet2VsPatchEnergy(0),
  fh3PtLeadJet1PatchEnergyVZEROAmp(0), fh3PtLeadJet1RawPatchEnergyVZEROAmp(0),
  fh3PatchEnergyEtaPhiCenterJ1(0), fh3PatchEnergyEtaPhiCenterJ2(0), fh3PatchEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1(0), fh3PatchADCEnergyEtaPhiCenterJ2(0), fh3PatchADCEnergyEtaPhiCenterJ1J2(0),
  fh3PatchEnergyEtaPhiCenterG1(0), fh3PatchEnergyEtaPhiCenterG2(0), fh3PatchEnergyEtaPhiCenterG1G2(0),
  fh3PatchADCEnergyEtaPhiCenterG1(0), fh3PatchADCEnergyEtaPhiCenterG2(0), fh3PatchADCEnergyEtaPhiCenterG1G2(0),
  fh3PatchADCEnergyEtaPhiCenterAll(0),
  fh3EEtaPhiCell(0), fh2ECellVsCent(0), fh2CellEnergyVsTime(0), fh3EClusELeadingCellVsTime(0),
  fh3JetReacCent(0),
  fHistClusEnergy(0), fHistClusofJetEnergy(0),
  fHistEventSelectionQA(0), fhQAinfoAllPatchesCounter(0), fhQAinfoCounter(0), fhQAmaxinfoCounter(0),
  fhRecalcGammaPatchEnergy(0), fhGammaLowPatchEnergy(0), fhGammaLowSimplePatchEnergy(0),
  fhRecalcJetPatchEnergy(0), fhJetLowPatchEnergy(0), fhJetLowSimplePatchEnergy(0),
  fhMainTriggerPatchEnergy(0),
  fJetTriggeredEvent(0), fRecalcTriggerPatches(0),
  fhnPatchMaxClus(0x0), fhnPatchMatch(0x0), fhnPatchMatch2(0x0), fhnPatchMatchJetLeadClus(0x0)
{
  // Standard constructor.
  for(Int_t j=0; j<16; j++) {
    fHistdPhidEtaPatchJetCluster[j] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerPatchJetMatch::~AliAnalysisTaskEmcalTriggerPatchJetMatch()
{
  // Destructor.
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTriggerPatchJetMatch::ExecOnce(){
  AliAnalysisTaskEmcalJet::ExecOnce();

  // Init the analysis
  if(fJetTriggeredEventname=="") fJetTriggeredEventname = Form("fJetTriggeredEventname_%s", GetName());
  fJetTriggeredEvent = new TClonesArray("AliEmcalJet");
  fJetTriggeredEvent->SetName(fJetTriggeredEventname);    
  fJetTriggeredEvent->SetOwner(kTRUE);

  fRecalcTriggerPatches = new TClonesArray("AliEmcalTriggerPatchInfo");
  fRecalcTriggerPatches->SetName("RecalcTriggerPatches");
  fRecalcTriggerPatches->SetOwner(kTRUE);

  // add jet object (of jets triggering event) to event if not yet there
  if(fAttachToEvent) {
    if (InputEvent()->FindListObject(fJetTriggeredEventname)) {
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fJetTriggeredEventname.Data()));
    } else {
      InputEvent()->AddObject(fJetTriggeredEvent);
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTriggerPatchJetMatch::SelectEvent() {
  // Decide if event should be selected for analysis
  fhNEvents->Fill(3.5);

  // this section isn't needed for LHC11h
  if(!fTriggerClass.IsNull()) {
    //Check if requested trigger was fired
    TString trigType1 = "J1";
    TString trigType2 = "J2";
    if(fTriggerClass.Contains("G")) {
      trigType1 = "G1";
      trigType2 = "G2";
    }

    // get fired trigger classes
    TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();

    if(fTriggerClass.Contains(trigType1.Data()) && fTriggerClass.Contains(trigType2.Data())) { //if events with J1&&J2 are requested
      if(!firedTrigClass.Contains(trigType1.Data()) || !firedTrigClass.Contains(trigType2.Data()) ) //check if both are fired
        return kFALSE;
    } else {
      if(!firedTrigClass.Contains(fTriggerClass)) return kFALSE;
      else if(fTriggerClass.Contains(trigType1.Data()) && firedTrigClass.Contains(trigType2.Data())) //if J2 is requested also add triggers which have J1&&J2. Reject if J1 is requested and J2 is fired
	    return kFALSE;
    }
  }

  fhNEvents->Fill(1.5);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerPatchJetMatch::FillTriggerPatchHistos() {
  // Fill trigger patch histos for main trigger
  if(!fTriggerPatchInfo) {
    AliFatal(Form("%s: TriggerPatchInfo object %s does not exist. Aborting", GetName(), fJetTriggeredEventname.Data()));
    return;
  } 

  // see if event was selected
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  ExtractMainPatch();
  //if(fMainPatchType == kManual) ExtractMainPatch();
  // not set up yet for jet patch
  //else if(fMainPatchType == kEmcalJet) cout<<"kEmcalJet"<<endl; //fMaxPatch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);
  //AliEmcalTriggerPatchInfo *patch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);

  fMaxPatchEnergy = 0;
  fMaxPatchADCEnergy = 0;

  if(fMaxPatch) {
    fMaxPatchEnergy = fMaxPatch->GetPatchE();
    fMaxPatchADCEnergy = fMaxPatch->GetADCAmpGeVRough();
    fh3PatchADCEnergyEtaPhiCenterAll->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());

    // the following don't exist for LHC11h data
    //Check if requested trigger was fired, relevant for data sets with 2 EMCal trigger thresholds
    // JET trigger
    if(fMaxPatch->IsJetLow() && !fMaxPatch->IsJetHigh()) { //main patch only fired low threshold trigger
      fh3PatchEnergyEtaPhiCenterJ2->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ2->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    }
    else if(fMaxPatch->IsJetHigh() && !fMaxPatch->IsJetLow()) { //main patch only fired high threshold trigger - should never happen
      fh3PatchEnergyEtaPhiCenterJ1->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ1->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    }
    else if(fMaxPatch->IsJetHigh() && fMaxPatch->IsJetLow()) { //main patch fired both triggers
      fh3PatchEnergyEtaPhiCenterJ1J2->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ1J2->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    } // JE

    // GAMMA trigger
    if(fMaxPatch->IsGammaLow() && !fMaxPatch->IsGammaHigh()) { //main patch only fired low threshold trigger
      fh3PatchEnergyEtaPhiCenterG2->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterG2->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    }
    else if(fMaxPatch->IsGammaHigh() && !fMaxPatch->IsGammaLow()) { //main patch only fired high threshold trigger - should never happen
      fh3PatchEnergyEtaPhiCenterG1->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterG1->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    }
    else if(fMaxPatch->IsGammaHigh() && fMaxPatch->IsGammaLow()) { //main patch fired both triggers
      fh3PatchEnergyEtaPhiCenterG1G2->Fill(fMaxPatch->GetPatchE(),fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterG1G2->Fill(fMaxPatchADCEnergy,fMaxPatch->GetEtaGeo(),fMaxPatch->GetPhiGeo());
    } // GA

    // fill max patch counters 
    FillTriggerPatchQA(fhQAmaxinfoCounter, trig, fMaxPatch);

  } // have patch

}

//_________________________________________________________________________________________
void AliAnalysisTaskEmcalTriggerPatchJetMatch::ExtractMainPatch() {
  //Find main trigger
  if(!fTriggerPatchInfo) return;

  // reset array of patches to save
  fRecalcTriggerPatches->Clear();

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntriesFast();

  //loop over patches to define trigger type of event
  Int_t nG1 = 0, nG2 = 0, nJ1 = 0, nJ2 = 0, nL0 = 0;

  // see if event was selected
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  // check the Fired Trigger Classes
  TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();

  //extract main trigger patch
  AliEmcalTriggerPatchInfo *patch;
  Double_t emax = -1.;
  for (Int_t iPatch = 0, patchacc = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEmcalTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (!patch) continue;

    // count trigger types
    if (patch->IsGammaHigh()) nG1++;
    if (patch->IsGammaLow())  nG2++;
    if (patch->IsJetHigh()) nJ1++;
    if (patch->IsJetLow())  nJ2++;
    if (patch->IsLevel0())  nL0++;

    // fill Energy spectra of recalculated Jet and GA patches
    if(patch->IsRecalcGamma()) { fhRecalcGammaPatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsGammaLow()) { fhGammaLowPatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsGammaLowSimple()) { fhGammaLowSimplePatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsRecalcJet()) { fhRecalcJetPatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsJetLow()) { fhJetLowPatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsJetLowSimple()) { fhJetLowSimplePatchEnergy->Fill(patch->GetPatchE()); }
    if(patch->IsMainTrigger()) { fhMainTriggerPatchEnergy->Fill(patch->GetPatchE()); }

    // fill QA counter histo for all 'trig class' patches
    FillTriggerPatchQA(fhQAinfoAllPatchesCounter, trig, patch);

    // make sure trigger "fTriggerClass" actually was fired in event
    if(!firedTrigClass.Contains(fTriggerClass)) continue;

    // check that we have a recalculated (OFFLINE) trigger patch of type fTriggerCategory
    if(fTriggerCategory == kTriggerRecalcGamma) {
      if(!patch->IsRecalcGamma()) continue;
    } else if (fTriggerCategory == kTriggerRecalcJet) {
      if(!patch->IsRecalcGamma()) continue;
    }

    // method for filling collection output array of 'saved' recalculated patches
    (*fRecalcTriggerPatches)[patchacc] = patch;
    ++patchacc;

    // fill QA counter histo for Recalculated 'trig class' patches
    FillTriggerPatchQA(fhQAinfoCounter, trig, patch);

    // find max patch energy
    if(patch->GetPatchE()>emax) {
      fMaxPatch = patch;
      emax = patch->GetPatchE();
    } // max patch
  } // loop over patches

  // check on patch Energy
  if(firedTrigClass.Contains(fTriggerClass) && fMaxPatch && fMaxPatch->GetPatchE() > fPatchECut) {
    // get cluster container and find leading cluster
    AliClusterContainer  *clusContp = GetClusterContainer(0);
    if(!clusContp) {
      AliError(Form("ERROR: Cluster container doesn't exist\n"));
      return;
    }
    AliVCluster* leadclus = clusContp->GetLeadingCluster();
    if(!leadclus) return;

    // initialize variables and get leading cluster parameters
    double leadclusEta = 0, leadclusPhi = 0, leadclusE = 0;
    TLorentzVector clusvect;
    leadclus->GetMomentum(clusvect, const_cast<Double_t*>(fVertex));
    leadclusEta = clusvect.Eta();
    leadclusPhi = clusvect.Phi();
    leadclusE = leadclus->E();

    // get patch variables
    double fMaxPatchPhiCM = fMaxPatch->GetPhiCM();
    double fMaxPatchPhiGeo = fMaxPatch->GetPhiGeo();
    double fMaxPatchEtaCM = fMaxPatch->GetEtaCM();
    double fMaxPatchEtaGeo = fMaxPatch->GetEtaGeo();
    double fMaxPatchPhiMin = fMaxPatch->GetPhiMin();
    double fMaxPatchPhiMax = fMaxPatch->GetPhiMax();
    double fMaxPatchEtaMin = fMaxPatch->GetEtaMin();
    double fMaxPatchEtaMax = fMaxPatch->GetEtaMax();
    Int_t nTrigBit = fMaxPatch->GetTriggerBits();

    // positional variables
    double dEtaGeo = 1.0*TMath::Abs(leadclusEta - fMaxPatchEtaGeo);
    double dEtaCM = 1.0*TMath::Abs(leadclusEta - fMaxPatchEtaCM);
    double dPhiGeo = 1.0*TMath::Abs(leadclusPhi - fMaxPatchPhiGeo);
    double dPhiCM = 1.0*TMath::Abs(leadclusPhi - fMaxPatchPhiCM);
    double maxPatchE = fMaxPatch->GetPatchE();
    double maxPatchADC = fMaxPatch->GetADCAmpGeVRough();

    // patch summary (meeting energy cut requirement)
    if(doComments) {
      cout<<endl<<"Patch summary: "<<endl;
      cout<<Form("Number of patches = %d, Max Patch Energy = %f GeV, MAXClusterE = %f, Phi = %f, Eta = %f", nPatch, fMaxPatch->GetPatchE(), leadclus->E(), leadclusPhi, leadclusEta)<<endl;
      cout<<Form("CM in Phi = %f, in Eta = %f, Geo Center in Phi = %f, in Eta = %f, TriggerBits = %d", fMaxPatchPhiCM, fMaxPatchEtaCM, fMaxPatchPhiGeo, fMaxPatchEtaGeo, nTrigBit)<<endl;
      cout<<Form("phi min = %f, phi max = %f, eta min = %f, eta max = %f", fMaxPatchPhiMin, fMaxPatchPhiMax, fMaxPatchEtaMin, fMaxPatchEtaMax)<<endl;
    }

    Double_t fill[7] = {dEtaGeo, dEtaCM, dPhiGeo, dPhiCM, maxPatchADC, maxPatchE, leadclusE};
    fhnPatchMaxClus->Fill(fill);

  } // patch energy cut

//  cout<<Form("Jet:   low: %d, high: %d," ,nJ2, nJ1)<<"   ";//<<endl;
//  cout<<Form("Gamma: low: %d, high: %d," ,nG2, nG1)<<"   ";//<<endl;
//  cout<<Form("L0: %d", nL0)<<endl;
//  cout<<Form("Max Patch Energy: %f GeV", fMaxPatch->GetPatchE())<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerPatchJetMatch::UserCreateOutputObjects()
{
  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhNEvents = new TH1F("fhNEvents","fhNEvents;selection;N_{evt}",5,0,5);
  fOutput->Add(fhNEvents);

  fhTriggerbit = new TProfile("fhTriggerbit","fhTriggerbit;;TriggerBit",1,0,1);
  fOutput->Add(fhTriggerbit);

  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, 300, 0., 300.);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("s#rho_{ch} (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCent);
    
  Int_t fgkNCentBins = 21;
  Float_t kMinCent   = 0.;
  Float_t kMaxCent   = 105.;
  Double_t *binsCent = new Double_t[fgkNCentBins+1];
  for(Int_t i=0; i<=fgkNCentBins; i++) binsCent[i]=(Double_t)kMinCent + (kMaxCent-kMinCent)/fgkNCentBins*(Double_t)i ;
  binsCent[fgkNCentBins-1] = 100.5;
  binsCent[fgkNCentBins] = 101.5;
    
  Int_t fgkNdEPBins = 18*8;
  Float_t kMindEP   = 0.;
  Float_t kMaxdEP   = 1.*TMath::Pi()/2.;
  Double_t *binsdEP = new Double_t[fgkNdEPBins+1];
  for(Int_t i=0; i<=fgkNdEPBins; i++) binsdEP[i]=(Double_t)kMindEP + (kMaxdEP-kMindEP)/fgkNdEPBins*(Double_t)i ;

  Int_t fgkNPtBins = 200;
  Float_t kMinPt   = -50.;
  Float_t kMaxPt   = 150.;
  Double_t *binsPt = new Double_t[fgkNPtBins+1];
  for(Int_t i=0; i<=fgkNPtBins; i++) binsPt[i]=(Double_t)kMinPt + (kMaxPt-kMinPt)/fgkNPtBins*(Double_t)i ;

  Int_t fgkNPhiBins = 18*8;
  Float_t kMinPhi   = 0.;
  Float_t kMaxPhi   = 2.*TMath::Pi();
  Double_t *binsPhi = new Double_t[fgkNPhiBins+1];
  for(Int_t i=0; i<=fgkNPhiBins; i++) binsPhi[i]=(Double_t)kMinPhi + (kMaxPhi-kMinPhi)/fgkNPhiBins*(Double_t)i ;

  Int_t fgkNEtaBins = 100;
  Float_t fgkEtaMin = -1.;
  Float_t fgkEtaMax =  1.;
  Double_t *binsEta=new Double_t[fgkNEtaBins+1];
  for(Int_t i=0; i<=fgkNEtaBins; i++) binsEta[i]=(Double_t)fgkEtaMin + (fgkEtaMax-fgkEtaMin)/fgkNEtaBins*(Double_t)i ;

  Int_t fgkNAreaBins = 100;
  Float_t kMinArea   = 0.;
  Float_t kMaxArea   = 1.;
  Double_t *binsArea = new Double_t[fgkNAreaBins+1];
  for(Int_t i=0; i<=fgkNAreaBins; i++) binsArea[i]=(Double_t)kMinArea + (kMaxArea-kMinArea)/fgkNAreaBins*(Double_t)i ;

  Int_t fgkNConstBins = 100;
  Float_t kMinConst   = 0.;
  Float_t kMaxConst   = 100.;
  Double_t *binsConst = new Double_t[fgkNConstBins+1];
  for(Int_t i=0; i<=fgkNConstBins; i++) binsConst[i]=(Double_t)kMinConst + (kMaxConst-kMinConst)/fgkNConstBins*(Double_t)i ;

  Int_t fgkNMeanPtBins = 100;
  Float_t kMinMeanPt   = 0.;
  Float_t kMaxMeanPt   = 20.;
  Double_t *binsMeanPt = new Double_t[fgkNMeanPtBins+1];
  for(Int_t i=0; i<=fgkNMeanPtBins; i++) binsMeanPt[i]=(Double_t)kMinMeanPt + (kMaxMeanPt-kMinMeanPt)/fgkNMeanPtBins*(Double_t)i ;

  Int_t fgkNNEFBins = 101;
  Float_t kMinNEF   = 0.;
  Float_t kMaxNEF   = 1.01;
  Double_t *binsNEF = new Double_t[fgkNNEFBins+1];
  for(Int_t i=0; i<=fgkNNEFBins; i++) binsNEF[i]=(Double_t)kMinNEF + (kMaxNEF-kMinNEF)/fgkNNEFBins*(Double_t)i ;

  Int_t fgkNzBins = 101;
  Float_t kMinz   = 0.;
  Float_t kMaxz   = 1.01;
  Double_t *binsz = new Double_t[fgkNzBins+1];
  for(Int_t i=0; i<=fgkNzBins; i++) binsz[i]=(Double_t)kMinz + (kMaxz-kMinz)/fgkNzBins*(Double_t)i ;

  Int_t fgkNJetTypeBins = 2;
  Float_t kMinJetType   = -0.5;
  Float_t kMaxJetType   = 1.5;
  Double_t *binsJetType = new Double_t[fgkNJetTypeBins+1];
  for(Int_t i=0; i<=fgkNJetTypeBins; i++) binsJetType[i]=(Double_t)kMinJetType + (kMaxJetType-kMinJetType)/fgkNJetTypeBins*(Double_t)i ;

  Int_t fgkNTimeBins = 100;
  Float_t kMinTime   = -200.;
  Float_t kMaxTime   = 200;
  Double_t *binsTime = new Double_t[fgkNTimeBins+1];
  for(Int_t i=0; i<=fgkNTimeBins; i++) binsTime[i]=(Double_t)kMinTime + (kMaxTime-kMinTime)/fgkNTimeBins*(Double_t)i ;

  Int_t fgkNVZEROBins = 100;
  Float_t kMinVZERO   = 0.;
  Float_t kMaxVZERO   = 25000;
  Double_t *binsVZERO = new Double_t[fgkNVZEROBins+1];
  for(Int_t i=0; i<=fgkNVZEROBins; i++) binsVZERO[i]=(Double_t)kMinVZERO + (kMaxVZERO-kMinVZERO)/fgkNVZEROBins*(Double_t)i ;

  Double_t enBinEdges[3][2];
  enBinEdges[0][0] = 1.; //10 bins
  enBinEdges[0][1] = 0.1;
  enBinEdges[1][0] = 5.; //8 bins
  enBinEdges[1][1] = 0.5;
  enBinEdges[2][0] = 100.;//95 bins
  enBinEdges[2][1] = 1.;

  const Float_t enmin1 =  0;
  const Float_t enmax1 =  enBinEdges[0][0];
  const Float_t enmin2 =  enmax1 ;
  const Float_t enmax2 =  enBinEdges[1][0];
  const Float_t enmin3 =  enmax2 ;
  const Float_t enmax3 =  enBinEdges[2][0];//fgkEnMax;
  const Int_t nbin11 = (int)((enmax1-enmin1)/enBinEdges[0][1]);
  const Int_t nbin12 = (int)((enmax2-enmin2)/enBinEdges[1][1])+nbin11;
  const Int_t nbin13 = (int)((enmax3-enmin3)/enBinEdges[2][1])+nbin12;

  Int_t fgkNEnBins=nbin13;
  Double_t *binsEn=new Double_t[fgkNEnBins+1];
  for(Int_t i=0; i<=fgkNEnBins; i++) {
    if(i<=nbin11) binsEn[i]=(Double_t)enmin1 + (enmax1-enmin1)/nbin11*(Double_t)i ;
    if(i<=nbin12 && i>nbin11) binsEn[i]=(Double_t)enmin2 + (enmax2-enmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;
    if(i<=nbin13 && i>nbin12) binsEn[i]=(Double_t)enmin3 + (enmax3-enmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;
  }

  fh3PtEtaPhiTracks = new TH3F("fh3PtEtaPhiTracks","fh3PtEtaPhiTracks;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracks);

  fh3PtEtaPhiTracksOnEmcal = new TH3F("fh3PtEtaPhiTracksOnEmcal","fh3PtEtaPhiTracksOnEmcal;#it{p}_{T}^{track}_{emc};#eta_{emc};#varphi_{emc}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksOnEmcal);

  fh3PtEtaPhiTracksToProp = new TH3F("fh3PtEtaPhiTracksToProp","fh3PtEtaPhiTracksToProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksToProp);

  fh3PtEtaPhiTracksProp = new TH3F("fh3PtEtaPhiTracksProp","fh3PtEtaPhiTracksProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksProp);

  fh3PtEtaPhiTracksNoProp = new TH3F("fh3PtEtaPhiTracksNoProp","fh3PtEtaPhiTracksNoProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksNoProp);

  fh2CentPtJet = new TH2F("fh2CentPtJet","fh2CentPtJet;cent;#it{p}_{T}^{jet}",fgkNCentBins,binsCent,fgkNPtBins,binsPt);
  fOutput->Add(fh2CentPtJet);

  fh3PtEtaPhiJet = new TH3F("fh3PtEtaPhiJet","fh3PtEtaPhiJet;#it{p}_{T}^{jet};#eta;#varphi",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiJet);

  fh2NJetsPt = new TH2F("fh2NJetsPt","fh2NJetsPt;N_{jets};#it{p}_{T}^{jet}",20,-0.5,19.5,fgkNPtBins,binsPt);
  fOutput->Add(fh2NJetsPt);

  fh3PtEtaAreaJet = new TH3F("fh3PtEtaAreaJet","fh3PtEtaAreaJet;#it{p}_{T}^{jet};#eta;A",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNAreaBins,binsArea);
  fOutput->Add(fh3PtEtaAreaJet);

  fh2PtNConstituentsCharged = new TH2F("fh2PtNConstituentsCharged","fh2PtNConstituentsCharged;#it{p}_{T}^{jet};N_{charged constituents}",fgkNPtBins,binsPt,fgkNConstBins,binsConst);
  fOutput->Add(fh2PtNConstituentsCharged);

  fh2PtNConstituents = new TH2F("fh2PtNConstituents","fh2PtNConstituents;#it{p}_{T}^{jet};N_{constituents}",fgkNPtBins,binsPt,fgkNConstBins,binsConst);
  fOutput->Add(fh2PtNConstituents);

  fh2PtMeanPtConstituentsCharged = new TH2F("fh2PtMeanPtConstituentsCharged","fh2PtMeanPtConstituentsCharged;#it{p}_{T}^{jet};charged #langle #it{p}_{T} #rangle",fgkNPtBins,binsPt,fgkNMeanPtBins,binsMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsCharged);

  fh2PtMeanPtConstituentsNeutral = new TH2F("fh2PtMeanPtConstituentsNeutral","fh2PtMeanPtConstituentsNeutral;#it{p}_{T}^{jet};neutral langle #it{p}_{T} #rangle",fgkNPtBins,binsPt,fgkNMeanPtBins,binsMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsNeutral);

  fh2PtNEF = new TH2F("fh2PtNEF","fh2PtNEF;#it{p}_{T}^{jet};NEF",fgkNPtBins,binsPt,fgkNNEFBins,binsNEF);
  fOutput->Add(fh2PtNEF);

  fh3NEFEtaPhi = new TH3F("fh3NEFEtaPhi","fh3NEFEtaPhi;NEF;#eta;#varphi",fgkNNEFBins,binsNEF,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3NEFEtaPhi);

  fh2NEFNConstituentsCharged = new TH2F("fh2NEFNConstituentsCharged","fh2NEFNConstituentsCharged;NEF;N_{charged constituents}",fgkNNEFBins,binsNEF,fgkNConstBins,binsConst);
  fOutput->Add(fh2NEFNConstituentsCharged);

  fh2NEFNConstituentsNeutral = new TH2F("fh2NEFNConstituentsNeutral","fh2NEFNConstituentsNeutral;NEF;N_{clusters}",fgkNNEFBins,binsNEF,fgkNConstBins,binsConst);
  fOutput->Add(fh2NEFNConstituentsNeutral);

  fh2Ptz = new TH2F("fh2Ptz","fh2Ptz;#it{p}_{T}^{jet};z=p_{t,trk}^{proj}/p_{jet}",fgkNPtBins,binsPt,fgkNzBins,binsz);
  fOutput->Add(fh2Ptz);

  fh2PtLeadJet1VsLeadJet2 = new TH2F("fh2PtLeadJet1VsLeadJet2","fh2PtLeadJet1VsLeadJet2;#it{p}_{T}^{jet 1};#it{p}_{T}^{jet 2}",fgkNPtBins,binsPt,fgkNPtBins,binsPt);
  fOutput->Add(fh2PtLeadJet1VsLeadJet2);

  fh3EEtaPhiCluster = new TH3F("fh3EEtaPhiCluster","fh3EEtaPhiCluster;E_{clus};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3EEtaPhiCluster);

  fh3PtLeadJet1VsPatchEnergy = new TH3F("fh3PtLeadJet1VsPatchEnergy","fh3PtLeadJet1VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNJetTypeBins,binsJetType);
  fOutput->Add(fh3PtLeadJet1VsPatchEnergy);
  fh3PtLeadJet2VsPatchEnergy = new TH3F("fh3PtLeadJet2VsPatchEnergy","fh3PtLeadJet2VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNJetTypeBins,binsJetType);
  fOutput->Add(fh3PtLeadJet2VsPatchEnergy);

  fh3PtLeadJet1PatchEnergyVZEROAmp = new TH3F("fh3PtLeadJet1PatchEnergyVZEROAmp","fh3PtLeadJet1VsPatchEnergyVZEROAmp;#it{p}_{T}^{jet 1};Amplitude_{patch};VZERO amp",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNVZEROBins,binsVZERO);
  fOutput->Add(fh3PtLeadJet1PatchEnergyVZEROAmp);
  fh3PtLeadJet1RawPatchEnergyVZEROAmp = new TH3F("fh3PtLeadJet1RawPatchEnergyVZEROAmp","fh3PtLeadJet1RawPatchEnergyVZEROAmp;#it{p}_{T}^{jet 1};ADC Amplitude_{patch} (GeV);VZERO amp",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNVZEROBins,binsVZERO);
  fOutput->Add(fh3PtLeadJet1RawPatchEnergyVZEROAmp);

  fh3PatchEnergyEtaPhiCenterJ1 = new TH3F("fh3PatchEnergyEtaPhiCenterJ1","fh3PatchEnergyEtaPhiCenterJ1;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ1);

  fh3PatchEnergyEtaPhiCenterJ2 = new TH3F("fh3PatchEnergyEtaPhiCenterJ2","fh3PatchEnergyEtaPhiCenterJ2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ2);

  fh3PatchEnergyEtaPhiCenterJ1J2 = new TH3F("fh3PatchEnergyEtaPhiCenterJ1J2","fh3PatchEnergyEtaPhiCenterJ1J2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ1J2);

  fh3PatchADCEnergyEtaPhiCenterJ1 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ1","fh3PatchADCEnergyEtaPhiCenterJ1;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ1);

  fh3PatchADCEnergyEtaPhiCenterJ2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ2","fh3PatchADCEnergyEtaPhiCenterJ2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ2);

  fh3PatchADCEnergyEtaPhiCenterJ1J2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ1J2","fh3PatchADCEnergyEtaPhiCenterJ1J2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ1J2);

  fh3PatchEnergyEtaPhiCenterG1 = new TH3F("fh3PatchEnergyEtaPhiCenterG1","fh3PatchEnergyEtaPhiCenterG1;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterG1);

  fh3PatchEnergyEtaPhiCenterG2 = new TH3F("fh3PatchEnergyEtaPhiCenterG2","fh3PatchEnergyEtaPhiCenterG2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterG2);

  fh3PatchEnergyEtaPhiCenterG1G2 = new TH3F("fh3PatchEnergyEtaPhiCenterG1G2","fh3PatchEnergyEtaPhiCenterG1G2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterG1G2);

  fh3PatchADCEnergyEtaPhiCenterG1 = new TH3F("fh3PatchADCEnergyEtaPhiCenterG1","fh3PatchADCEnergyEtaPhiCenterG1;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterG1);

  fh3PatchADCEnergyEtaPhiCenterG2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterG2","fh3PatchADCEnergyEtaPhiCenterG2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterG2);

  fh3PatchADCEnergyEtaPhiCenterG1G2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterG1G2","fh3PatchADCEnergyEtaPhiCenterG1G2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterG1G2);

  fh3PatchADCEnergyEtaPhiCenterAll = new TH3F("fh3PatchADCEnergyEtaPhiCenterAll","fh3PatchADCEnergyEtaPhiCenterAll;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterAll);

  fh3EEtaPhiCell = new TH3F("fh3EEtaPhiCell","fh3EEtaPhiCell;E_{cell};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3EEtaPhiCell);

  fh2ECellVsCent = new TH2F("fh2ECellVsCent","fh2ECellVsCent;centrality;E_{cell}",101,-1,100,500,0.,5.);
  fOutput->Add(fh2ECellVsCent);

  fh2CellEnergyVsTime = new TH2F("fh2CellEnergyVsTime","fh2CellEnergyVsTime;E_{cell};time",fgkNEnBins,binsEn,fgkNTimeBins,binsTime);
  fOutput->Add(fh2CellEnergyVsTime);

  fh3EClusELeadingCellVsTime = new TH3F("fh3EClusELeadingCellVsTime","fh3EClusELeadingCellVsTime;E_{cluster};E_{leading cell};time_{leading cell}",fgkNEnBins,binsEn,fgkNEnBins,binsEn,fgkNTimeBins,binsTime);
  fOutput->Add(fh3EClusELeadingCellVsTime);
    
  fh3JetReacCent = new TH3F("fh3JetReacCent","fh3JetReacCent;E_{Jet};Centrality;dEP",fgkNEnBins,binsEn,fgkNCentBins,binsCent,fgkNdEPBins,binsdEP);
  fOutput->Add(fh3JetReacCent);

  fhQAinfoAllPatchesCounter = new TH1F("fhQAinfoAllPatchesCounter", "QA trigger info counters for all patches", 20, 0.5, 20.5);
  fOutput->Add(fhQAinfoAllPatchesCounter);

  fhQAinfoCounter = new TH1F("fhQAinfoCounter", "QA trigger info counters", 20, 0.5, 20.5);
  fOutput->Add(fhQAinfoCounter);
    
  fhQAmaxinfoCounter = new TH1F("fhQAmaxinfoCounter", "QA Max patch trigger info counters", 20, 0.5, 20.5);
  fOutput->Add(fhQAmaxinfoCounter);

  fhRecalcGammaPatchEnergy = new TH1F("fhRecalcGammaPatchEnergy", "Recalculated Gamma Patch Energy", 200, 0, 50); //100
  fOutput->Add(fhRecalcGammaPatchEnergy);

  fhGammaLowPatchEnergy = new TH1F("fhGammaLowPatchEnergy", "Gamma Low Patch Energy", 200, 0, 50); //100
  fOutput->Add(fhGammaLowPatchEnergy);

  fhGammaLowSimplePatchEnergy = new TH1F("fhGammaLowSimplePatchEnergy", "Gamma Low Simple Patch Energy", 200, 0, 50); //100
  fOutput->Add(fhGammaLowSimplePatchEnergy);

  fhRecalcJetPatchEnergy = new TH1F("fhRecalcJetPatchEnergy", "Recalculated Jet Patch Energy", 200, 0, 100);
  fOutput->Add(fhRecalcJetPatchEnergy);

  fhJetLowPatchEnergy = new TH1F("fhJetLowPatchEnergy", "Jet Low Patch Energy", 200, 0, 100);
  fOutput->Add(fhJetLowPatchEnergy);

  fhJetLowSimplePatchEnergy = new TH1F("fhJetLowSimplePatchEnergy", "Jet Low Simple Patch Energy", 200, 0, 100);
  fOutput->Add(fhJetLowSimplePatchEnergy);

  fhMainTriggerPatchEnergy = new TH1F("fhMainTriggerPatchEnergy", "Main Trigger Patch Energy", 200, 0, 100);
  fOutput->Add(fhMainTriggerPatchEnergy);

  fHistClusEnergy = new TH1F("fHistClusEnergy", "Cluster Energy distribution", 200, 0, 50);
  fOutput->Add(fHistClusEnergy);


  fHistClusofJetEnergy = new TH1F("fHistClusofJetEnergy", "Cluster of Jet Energy distribution", 200, 0, 20);
  fOutput->Add(fHistClusofJetEnergy);

  for(Int_t j=0; j<16; j++) {
    fHistdPhidEtaPatchJetCluster[j] = new TH2F(Form("fHistdPhidEtaPatchJetCluster_%d",j), "dPhi-dEta distribution between max recalculated Gamma patch and most energetic jet cluster", 144, 0, 2.016, 144, 0, 2.016);
    fOutput->Add(fHistdPhidEtaPatchJetCluster[j]);
  }

  Int_t nDim=7;
  Int_t *nbins = new Int_t[nDim];
  Double_t *xmin = new Double_t[nDim];
  Double_t *xmax = new Double_t[nDim]; 
  for (Int_t i=0; i<nDim; i++){
    nbins[i]=144;
    xmin[i]=0;
    xmax[i]=2.016;
  }

  nbins[4]=100; xmax[4]=100;
  nbins[5]=100; xmax[5]=100.;
  nbins[6]=100; xmax[6]=100.;

  fhnPatchMaxClus = new THnSparseF("fhnPatchMaxClus","fhn Patch Max Cluster Distributions", nDim,nbins,xmin,xmax);
  fhnPatchMaxClus->GetAxis(0)->SetTitle("#Delta#etaGeo");          // 0
  fhnPatchMaxClus->GetAxis(1)->SetTitle("#Delta#etaCM");           // 1
  fhnPatchMaxClus->GetAxis(2)->SetTitle("#Delta#phiGeo");          // 2
  fhnPatchMaxClus->GetAxis(3)->SetTitle("#Delta#phiCM");           // 3
  fhnPatchMaxClus->GetAxis(4)->SetTitle("Max Patch ADC");          // 4
  fhnPatchMaxClus->GetAxis(5)->SetTitle("Max Patch Energy");       // 5
  fhnPatchMaxClus->GetAxis(6)->SetTitle("Leading Cluster Energy"); // 6
  fOutput->Add(fhnPatchMaxClus);

  // QA before matching
  Int_t nDim1=8;
  Int_t *nbins1 = new Int_t[nDim1];
  Double_t *xmin1 = new Double_t[nDim1];
  Double_t *xmax1 = new Double_t[nDim1]; 
  nbins1[0]=10; xmin1[0]=0.; xmax1[0]=100.;
  nbins1[1]=250; xmin1[1]=0.; xmax1[1]=250.;
  nbins1[2]=50; xmin1[2]=0.; xmax1[2]=50;
  nbins1[3]=144; xmin1[3]=0.; xmax1[3]=2.016;
  nbins1[4]=144; xmin1[4]=0.; xmax1[4]=2.016;
  nbins1[5]=300; xmin1[5]=0.; xmax1[5]=300;
  nbins1[6]=500; xmin1[6]=0.; xmax1[6]=500;
  nbins1[7]=3; xmin1[7]=0.0; xmax1[7]=1.0*TMath::Pi()/2.0;

  // before cuts to perform match
  fhnPatchMatch = new THnSparseF("fhnPatchMatch","fhn Patch Match before cuts", nDim1,nbins1,xmin1,xmax1);
  fhnPatchMatch->GetAxis(0)->SetTitle("Centrality %");              // 0
  fhnPatchMatch->GetAxis(1)->SetTitle("Jet p_{T}");                 // 1
  fhnPatchMatch->GetAxis(2)->SetTitle("Max Cluster Energy");        // 2
  fhnPatchMatch->GetAxis(3)->SetTitle("#Delta#phi Geo");            // 3
  fhnPatchMatch->GetAxis(4)->SetTitle("#Delta#eta Geo");            // 4
  fhnPatchMatch->GetAxis(5)->SetTitle("Max Patch Energy");          // 5
  fhnPatchMatch->GetAxis(6)->SetTitle("Max Patch ADC");             // 6
  fhnPatchMatch->GetAxis(7)->SetTitle("Event Plane - Jet Angle");   // 7
  fOutput->Add(fhnPatchMatch);

  // after cuts to perform match
  fhnPatchMatch2 = new THnSparseF("fhnPatchMatch2","fhn Patch Match after cuts", nDim1,nbins1,xmin1,xmax1);
  fhnPatchMatch2->GetAxis(0)->SetTitle("Centrality %");              // 0
  fhnPatchMatch2->GetAxis(1)->SetTitle("Jet p_{T}");                 // 1
  fhnPatchMatch2->GetAxis(2)->SetTitle("Max Cluster Energy");        // 2
  fhnPatchMatch2->GetAxis(3)->SetTitle("#Delta#phi Geo");            // 3
  fhnPatchMatch2->GetAxis(4)->SetTitle("#Delta#eta Geo");            // 4
  fhnPatchMatch2->GetAxis(5)->SetTitle("Max Patch Energy");          // 5
  fhnPatchMatch2->GetAxis(6)->SetTitle("Max Patch ADC");             // 6
  fhnPatchMatch2->GetAxis(7)->SetTitle("Event Plane - Jet Angle");   // 7
  fOutput->Add(fhnPatchMatch2);

  // for jet cluster matched to patch
  Int_t nDim2=18;
  Int_t *nbins2 = new Int_t[nDim2];
  Double_t *xmin2 = new Double_t[nDim2];
  Double_t *xmax2 = new Double_t[nDim2]; 
  for (Int_t i=0; i<nDim2; i++){
    nbins2[i]=144;
    xmin2[i]=0;
  }

  nbins2[0]=10; xmax2[0]=100;
  nbins2[1]=3; xmax2[1]=1.0*TMath::Pi()/2.;
  nbins2[2]=250; xmax2[2]=250; 
  nbins2[3]=72; xmin2[3]=1.2; xmax2[3]=3.4; //*TMath::Pi();
  nbins2[4]=56; xmin2[4]=-0.7; xmax2[4]=0.7;
  nbins2[5]=40; xmax2[5]=2.;
  nbins2[6]=50; xmax2[6]=50;
  nbins2[7]=100; xmax2[7]=100.;
  nbins2[8]=50; xmax2[8]=50;
  nbins2[9]=100; xmax2[9]=100.;
  nbins2[10]=72; xmin2[10]=1.2; xmax2[10]=3.4; //2.0*TMath::Pi();
  nbins2[11]=56; xmin2[11]=-0.7; xmax2[11]=0.7;
  nbins2[12]=500; xmax2[12]=500.;
  nbins2[13]=500; xmax2[13]=500.;
  nbins2[14]=300; xmax2[14]=300.;
  nbins2[15]=300; xmax2[15]=300.;
  nbins2[16]=72; xmin2[16]=1.4; xmax2[16]=3.2;
  nbins2[17]=56; xmin2[17]=-0.7; xmax2[17]=0.7;

  //Double_t fill[18] = {fCent, dEPJet, jet->Pt(), jet->Phi(), jet->Eta(), jet->Area(), jet->GetNumberOfTracks(), jet->MaxTrackPt(), jet->GetNumberOfClusters(), maxClusterE, maxClusterPhi, maxClusterEta, kAmplitudeOnline, kAmplitudeOnline, kEnergyOnline, kEnergyOffline, fMaxPatchPhiGeo, fMaxPatchEtaGeo};                

  fhnPatchMatchJetLeadClus = new THnSparseF("fhnPatchMatchJetLeadClus","fhn Patch Match to Jet Leading Cluster", nDim2, nbins2,xmin2,xmax2);
  fhnPatchMatchJetLeadClus->GetAxis(0)->SetTitle("Centrality %");                          // 0
  fhnPatchMatchJetLeadClus->GetAxis(1)->SetTitle("Event Plane Orientation");               // 1
  fhnPatchMatchJetLeadClus->GetAxis(2)->SetTitle("Jet p_{T}");                             // 2
  fhnPatchMatchJetLeadClus->GetAxis(3)->SetTitle("Jet #phi");                              // 3
  fhnPatchMatchJetLeadClus->GetAxis(4)->SetTitle("Jet #eta");                              // 4
  fhnPatchMatchJetLeadClus->GetAxis(5)->SetTitle("Jet Area");                              // 5
  fhnPatchMatchJetLeadClus->GetAxis(6)->SetTitle("Number of Jet Tracks");                  // 6
  fhnPatchMatchJetLeadClus->GetAxis(7)->SetTitle("Max Track p_{T} constituent of Jet");    // 7
  fhnPatchMatchJetLeadClus->GetAxis(8)->SetTitle("Number of Jet Clusters");                // 8
  fhnPatchMatchJetLeadClus->GetAxis(9)->SetTitle("Max Cluster Energy constituent of Jet"); // 9
  fhnPatchMatchJetLeadClus->GetAxis(10)->SetTitle("Cluster #phi");                         // 10
  fhnPatchMatchJetLeadClus->GetAxis(11)->SetTitle("Cluster #eta");                         // 11
  fhnPatchMatchJetLeadClus->GetAxis(12)->SetTitle("Max Patch Amplitude Online");           // 12
  fhnPatchMatchJetLeadClus->GetAxis(13)->SetTitle("Max Patch Amplitude Offline");          // 13
  fhnPatchMatchJetLeadClus->GetAxis(14)->SetTitle("Max Patch Energy Online");              // 14
  fhnPatchMatchJetLeadClus->GetAxis(15)->SetTitle("Max Patch Energy Offline");             // 15
  fhnPatchMatchJetLeadClus->GetAxis(16)->SetTitle("Max Patch Geometric Center in #phi");   // 16
  fhnPatchMatchJetLeadClus->GetAxis(17)->SetTitle("Max Patch Geometric Center in #eta");   // 17
  fOutput->Add(fhnPatchMatchJetLeadClus);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fOutput->Add(fHistEventSelectionQA);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
    if (h2){
      h2->Sumw2();
      continue;
    }
    TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
    if (h3){
      h3->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  if(binsCent)              delete [] binsCent;
  if(binsdEP)               delete [] binsdEP;
  if(binsEn)                delete [] binsEn;
  if(binsPt)                delete [] binsPt;
  if(binsPhi)               delete [] binsPhi;
  if(binsEta)               delete [] binsEta;
  if(binsArea)              delete [] binsArea;
  if(binsConst)             delete [] binsConst; 
  if(binsMeanPt)            delete [] binsMeanPt;  
  if(binsNEF)               delete [] binsNEF;
  if(binsz)                 delete [] binsz;
  if(binsJetType)           delete [] binsJetType;
  if(binsTime)              delete [] binsTime;
  if(binsVZERO)             delete [] binsVZERO;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTriggerPatchJetMatch::FillHistograms() {
  // Fill histograms.

  // check and fill a Event Selection QA histogram for different trigger selections
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  // get fired trigger classes
  TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();

  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA, trig);

  // create pointer to list of input event
  TList *list = InputEvent()->GetList();
  if(!list) {
    AliError(Form("ERROR: list not attached\n"));
    return kTRUE;
  }

  // reset collect at start of event
  fJetTriggeredEvent->Clear();

  //Tracks
  AliParticleContainer *partCont = GetParticleContainer(0);
  if (partCont) {
    //AliPicoTrack *track = dynamic_cast<AliPicoTrack*>(partCont->GetNextAcceptParticle(0));
    AliVTrack *track = dynamic_cast<AliVTrack*>(partCont->GetNextAcceptParticle(0));
    //if(!track) cout<<"No tracks... "<<endl;
    while(track) {
      Double_t trkphi = track->Phi()*TMath::RadToDeg();
      fh3PtEtaPhiTracks->Fill(track->Pt(),track->Eta(),track->Phi());

      //Select tracks which should be propagated
      if(track->Pt()>=0.350) {
        if (TMath::Abs(track->Eta())<=0.9 && trkphi > 10 && trkphi < 250) {
          fh3PtEtaPhiTracksOnEmcal->Fill(track->GetTrackPtOnEMCal(),track->GetTrackEtaOnEMCal(),track->GetTrackPhiOnEMCal());
          fh3PtEtaPhiTracksToProp->Fill(track->Pt(),track->Eta(),track->Phi());
          if(track->GetTrackPtOnEMCal()>=0) fh3PtEtaPhiTracksProp->Fill(track->Pt(),track->Eta(),track->Phi());
        else
          fh3PtEtaPhiTracksNoProp->Fill(track->Pt(),track->Eta(),track->Phi());
        } // acceptance cut
      } // track pt cut
      track = dynamic_cast<AliPicoTrack*>(partCont->GetNextAcceptParticle());
    } // track exist
  } // particle container

  //Clusters
  AliClusterContainer  *clusCont = GetClusterContainer(0);
  if (clusCont) {
    Int_t nclusters = clusCont->GetNClusters();
    for (Int_t ic = 0; ic < nclusters; ic++) {
      AliVCluster *cluster = static_cast<AliVCluster*>(clusCont->GetCluster(ic));
      if (!cluster) {
	    AliDebug(2,Form("Could not receive cluster %d", ic));
	    continue;
      }

      if (!cluster->IsEMCAL()) {
        AliDebug(2,Form("%s: Cluster is not emcal",GetName()));
	    continue;
      }

      TLorentzVector lp;
      cluster->GetMomentum(lp, const_cast<Double_t*>(fVertex));
      fh3EEtaPhiCluster->Fill(lp.E(),lp.Eta(),lp.Phi());
      if(fCaloCells) {
        Double_t leadCellE = GetEnergyLeadingCell(cluster);
        Double_t leadCellT = cluster->GetTOF();
        fh3EClusELeadingCellVsTime->Fill(lp.E(),leadCellE,leadCellT*1e9);
      } // if calo cells exist
    } // cluster loop
  } // if cluster container exist

  //Clusters - get cluster collection
  TClonesArray *clusters = 0x0;
  clusters = dynamic_cast<TClonesArray*>(list->FindObject(fCaloClustersName));  // CHANGED
  if (!clusters) {
    AliError(Form("Pointer to clusters %s == 0", fCaloClustersName.Data()));  // CHANGED
    return kTRUE;
  } // verify existence of clusters

  // loop over clusters
  const Int_t Nclusters = clusters->GetEntries();
  for (Int_t iclus = 0; iclus < Nclusters; iclus++){
    AliVCluster* clus = static_cast<AliVCluster*>(clusters->At(iclus));
    if(!clus){
      AliError(Form("Couldn't get AliVCluster %d\n", iclus));
      continue;
    }

    // is cluster in EMCal?
    if (!clus->IsEMCAL()) {
      AliDebug(2,Form("%s: Cluster is not emcal",GetName()));
      continue;
    }

    // get some info on cluster and fill histo
    TLorentzVector lp;
    clus->GetMomentum(lp, const_cast<Double_t*>(fVertex));
    fHistClusEnergy->Fill(lp.E());
  }

  //Jets
  Double_t ptLeadJet1 = 0.;
  Double_t ptLeadJet2 = 0.;

  // background rho
  fHistRhovsCent->Fill(fCent,GetRhoVal(fJetContainer));

  TArrayI *nJetsArr = new TArrayI(fh2NJetsPt->GetNbinsY()+1);
  nJetsArr->Reset(0);
  nJetsArr->Set(fh2NJetsPt->GetNbinsY()+1);

  // get jet container
  if (GetJetContainer(fJetContainer)) {
    const Int_t njets = GetNJets(fJetContainer);

    // loop over jets
    for (Int_t ij = 0, jacc = 0; ij < njets; ij++) {
      AliEmcalJet* jet = GetAcceptJetFromArray(ij,fJetContainer);
      if (!jet) continue; //jet not selected

      // get jetpt after background subtraction
      Double_t jetPt = jet->Pt() - GetRhoVal(fJetContainer)*jet->Area();
      if(jetPt>ptLeadJet1) ptLeadJet1=jetPt;

      // fill histos
      Double_t dEPJet = RelativeEP(jet->Phi(), fEPV0);
      fh3JetReacCent->Fill(jet->E(),fCent,dEPJet);
      fh2CentPtJet->Fill(fCent,jetPt);
      fh3PtEtaPhiJet->Fill(jetPt,jet->Eta(),jet->Phi());
      fh3PtEtaAreaJet->Fill(jetPt,jet->Eta(),jet->Area());

      //count jets above certain pT threshold
      Int_t ptbin = fh2NJetsPt->GetYaxis()->FindBin(jetPt);
      for(Int_t iptbin = ptbin; iptbin<=fh2NJetsPt->GetNbinsY(); iptbin++)
        nJetsArr->AddAt(nJetsArr->At(iptbin)+1,iptbin);

      // fill histos
      fh2PtNConstituentsCharged->Fill(jetPt,jet->GetNumberOfTracks());
      fh2PtNConstituents->Fill(jetPt,jet->GetNumberOfConstituents());
      fh2PtNEF->Fill(jetPt,jet->NEF());
      fh3NEFEtaPhi->Fill(jet->NEF(),jet->Eta(),jet->Phi());
      fh2NEFNConstituentsCharged->Fill(jet->NEF(),jet->GetNumberOfTracks());
      fh2NEFNConstituentsNeutral->Fill(jet->NEF(),jet->GetNumberOfClusters());

      // sum up charged pt of jet
      AliVParticle *vp;
      Double_t sumPtCh = 0.;
      for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
        vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
        if(!vp) continue;
        fh2Ptz->Fill(jetPt,GetZ(vp,jet));
        sumPtCh+=vp->Pt();
      }
      
      if(jet->GetNumberOfTracks()>0) fh2PtMeanPtConstituentsCharged->Fill(jetPt,sumPtCh/(double)(jet->GetNumberOfTracks()) );

      AliVCluster *vc = 0x0;
      Double_t sumPtNe = 0.;
      Double_t maxClusterEta = 0., maxClusterPhi = 0., maxClusterE = 0.;

      // apply jet pt cut
      if(jet->Pt() > fJetPtCut) {
        // check for cluster containers
        if (clusCont && clusCont->GetArray()) {
          // get leading cluster of jet
          AliVCluster *clustermax = jet->GetLeadingCluster(clusCont->GetArray());

          // check that we have a leading cluster and it is above a set bias - for MAX patch
          if(clustermax && clustermax->E() > fClusBias  && fMaxPatch && firedTrigClass.Contains(fTriggerClass) && (!fUseALLrecalcPatches)) {
            TLorentzVector nPart;
            clustermax->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
            maxClusterEta = nPart.Eta();
            maxClusterPhi = nPart.Phi();
            maxClusterE = clustermax->E();

            // get max patch location
            double fMaxPatchPhiGeo = fMaxPatch->GetPhiGeo();
            double fMaxPatchEtaGeo = fMaxPatch->GetEtaGeo();
            double dPhiPatchLeadCl = 1.0*TMath::Abs(fMaxPatchPhiGeo - nPart.Phi());
            double dEtaPatchLeadCl = 1.0*TMath::Abs(fMaxPatchEtaGeo - nPart.Eta());

            // get some info to fill in sparse
            double dEPJet = RelativeEP(jet->Phi(), fEPV0);
            double kAmplitudeOnline = fMaxPatch->GetADCAmp();
            double kAmplitudeOffline = fMaxPatch->GetADCOfflineAmp();
            double kEnergyOnline = fMaxPatch->GetADCAmpGeVRough();
            double kEnergyOffline = fMaxPatch->GetPatchE();

            // get patch variables
            double etamin = TMath::Min(fMaxPatch->GetEtaMin(), fMaxPatch->GetEtaMax());
            double etamax = TMath::Max(fMaxPatch->GetEtaMin(), fMaxPatch->GetEtaMax());
            double phimin = TMath::Min(fMaxPatch->GetPhiMin(), fMaxPatch->GetPhiMax());
            double phimax = TMath::Max(fMaxPatch->GetPhiMin(), fMaxPatch->GetPhiMax());

            for(int maxbinE = 0; maxbinE<16; maxbinE++) {
               if(maxClusterE > maxbinE) fHistdPhidEtaPatchJetCluster[maxbinE]->Fill(dPhiPatchLeadCl, dEtaPatchLeadCl);
            }

            // fill sparse array before and after match
            Double_t fillarr[8] = {fCent, jet->Pt(), maxClusterE, dPhiPatchLeadCl, dEtaPatchLeadCl, kEnergyOffline, kAmplitudeOnline, dEPJet};
            fhnPatchMatch->Fill(fillarr);
            if(maxClusterEta > etamin && maxClusterEta < etamax && maxClusterPhi > phimin && maxClusterPhi < phimax) fhnPatchMatch2->Fill(fillarr);

            // patch meeting offline energy cut
            if(fMaxPatch->GetPatchE() > fPatchECut) {

              // look to geometrically match patch to leading cluster of jet
              if(maxClusterEta > etamin && maxClusterEta < etamax && maxClusterPhi > phimin && maxClusterPhi < phimax){
                if(doComments) {
                  cout<<"*********************************************"<<endl;
                  cout<<"Proper match (cluster constituent of jet to fired max patch: ";
                  cout<<Form("Jet # = %i, Jet Pt = %f, Jet Phi =  %f, Jet Eta = %f", ij, jet->Pt(), jet->Phi(), jet->Eta())<<endl;
                  cout<<Form("NClus in jet = %i, MaxClusterE = %f, Phi = %f, Eta = %f", jet->GetNumberOfClusters(), maxClusterE, maxClusterPhi, maxClusterEta)<<endl;
                  cout<<Form("LeadClusE = %f, Phi = %f, Eta = %f", maxClusterE, maxClusterPhi, maxClusterEta)<<endl;
                  cout<<"*********************************************"<<endl;
                } // do comments

                // fill sparse for match
                Double_t fill[18] = {fCent, dEPJet, jet->Pt(), jet->Phi(), jet->Eta(), jet->Area(), (Double_t)jet->GetNumberOfTracks(), jet->MaxTrackPt(), (Double_t)jet->GetNumberOfClusters(), maxClusterE, maxClusterPhi, maxClusterEta, kAmplitudeOnline, kAmplitudeOffline, kEnergyOnline, kEnergyOffline, fMaxPatchPhiGeo, fMaxPatchEtaGeo};
                fhnPatchMatchJetLeadClus->Fill(fill);

                // method for filling collection output array of 'saved' jets
                (*fJetTriggeredEvent)[jacc] = jet;
                ++jacc;

              } // MATCH!
            } // patch > Ecut GeV
          } // have max cluster


// ====================================================
          // check that we have a leading cluster and it is above a set bias - for ALL recalculated (OFFLINE) patches
          if(clustermax && clustermax->E() > fClusBias && firedTrigClass.Contains(fTriggerClass) && (fUseALLrecalcPatches)) {
            TLorentzVector nPart;
            clustermax->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
            maxClusterEta = nPart.Eta();
            maxClusterPhi = nPart.Phi();
            maxClusterE = clustermax->E();

            // number of recalculated patches in event
            Int_t nPatch = fRecalcTriggerPatches->GetEntriesFast();

            AliEmcalTriggerPatchInfo *patch;
            for(Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
              patch = (AliEmcalTriggerPatchInfo*)fRecalcTriggerPatches->At(iPatch);
              if(!patch) continue;

              // get max patch location
              double fPatchPhiGeo = patch->GetPhiGeo();
              double fPatchEtaGeo = patch->GetEtaGeo();
              double dPhiPatchLeadCl = 1.0*TMath::Abs(fPatchPhiGeo - maxClusterPhi);
              double dEtaPatchLeadCl = 1.0*TMath::Abs(fPatchEtaGeo - maxClusterEta);

              // get some info to fill in sparse
              double dEPJet = RelativeEP(jet->Phi(), fEPV0);
              double kAmplitudeOnline = patch->GetADCAmp();
              double kAmplitudeOffline = patch->GetADCOfflineAmp();
              double kEnergyOnline = patch->GetADCAmpGeVRough();
              double kEnergyOffline = patch->GetPatchE();

              // get patch variables
              double etamin = TMath::Min(patch->GetEtaMin(), patch->GetEtaMax());
              double etamax = TMath::Max(patch->GetEtaMin(), patch->GetEtaMax());
              double phimin = TMath::Min(patch->GetPhiMin(), patch->GetPhiMax());
              double phimax = TMath::Max(patch->GetPhiMin(), patch->GetPhiMax());

              for(int maxbinE = 0; maxbinE<16; maxbinE++) {
                 if(maxClusterE > maxbinE) fHistdPhidEtaPatchJetCluster[maxbinE]->Fill(dPhiPatchLeadCl, dEtaPatchLeadCl);
              }

              // fill sparse array before and after match
              Double_t fillarr[8] = {fCent, jet->Pt(), maxClusterE, dPhiPatchLeadCl, dEtaPatchLeadCl, kEnergyOffline, kAmplitudeOnline, dEPJet};
              fhnPatchMatch->Fill(fillarr);
              if(maxClusterEta > etamin && maxClusterEta < etamax && maxClusterPhi > phimin && maxClusterPhi < phimax) fhnPatchMatch2->Fill(fillarr);

              // patch meeting offline energy cut
              if(patch->GetPatchE() > fPatchECut) {
                // look to geometrically match patch to leading cluster of jet
                if(maxClusterEta > etamin && maxClusterEta < etamax && maxClusterPhi > phimin && maxClusterPhi < phimax){
                  if(doComments) {
                    cout<<"*********************************************"<<endl;
                    cout<<"Proper match (cluster constituent of jet to fired ALL recalculated (OFFLINE) patch: ";
                    cout<<Form("Jet # = %i, Jet Pt = %f, Jet Phi =  %f, Jet Eta = %f", ij, jet->Pt(), jet->Phi(), jet->Eta())<<endl;
                    cout<<Form("NClus in jet = %i, MaxClusterE = %f, Phi = %f, Eta = %f", jet->GetNumberOfClusters(), maxClusterE, maxClusterPhi, maxClusterEta)<<endl;
                    cout<<Form("LeadClusE = %f, Phi = %f, Eta = %f", maxClusterE, maxClusterPhi, maxClusterEta)<<endl;
                    cout<<"*********************************************"<<endl;
                  } // do comments

                  // fill sparse for match
                  Double_t fill[18] = {fCent, dEPJet, jet->Pt(), jet->Phi(), jet->Eta(), jet->Area(), (Double_t)jet->GetNumberOfTracks(), jet->MaxTrackPt(), (Double_t)jet->GetNumberOfClusters(), maxClusterE, maxClusterPhi, maxClusterEta, kAmplitudeOnline, kAmplitudeOffline, kEnergyOnline, kEnergyOffline, fPatchPhiGeo, fPatchEtaGeo};
                  fhnPatchMatchJetLeadClus->Fill(fill);

                  // method for filling collection output array of 'saved' jets
                  (*fJetTriggeredEvent)[jacc] = jet;
                  ++jacc;

                } // MATCH!
              } // patch > Ecut GeV
            }// patch loop
          } // have max cluster
//=======================================================
        } // have cluster container
      } // jet->pt() > cut

      // sum up neutral energy of jet
      if (clusCont) {
        for(Int_t icc=0; icc<jet->GetNumberOfClusters(); icc++) {
          vc = static_cast<AliVCluster*>(clusCont->GetCluster(icc));
          if(!vc) continue;
          TLorentzVector lp;
          vc->GetMomentum(lp, const_cast<Double_t*>(fVertex));
          sumPtNe+=lp.Pt();
        } // cluster loop

        if(jet->GetNumberOfClusters()>0) fh2PtMeanPtConstituentsNeutral->Fill(jetPt,sumPtNe/(double)(jet->GetNumberOfClusters()) );
      } // cluster container exists
    } // jet loop

    for(Int_t i=1; i<=fh2NJetsPt->GetNbinsY(); i++) {
      Int_t nJetsInEvent = nJetsArr->At(i);
      fh2NJetsPt->Fill(nJetsInEvent,fh2NJetsPt->GetYaxis()->GetBinCenter(i));
    }
  }

  // leading jet comparison with each other and max patches
  if(GetJetContainer(fJetContainer)) fh2PtLeadJet1VsLeadJet2->Fill(ptLeadJet1,ptLeadJet2);
  fh3PtLeadJet1VsPatchEnergy->Fill(ptLeadJet1,fMaxPatchEnergy,fTriggerType);
  fh3PtLeadJet2VsPatchEnergy->Fill(ptLeadJet2,fMaxPatchEnergy,fTriggerType);

  //Get VZERO amplitude
  Float_t VZEROAmp = InputEvent()->GetVZEROData()->GetTriggerChargeA() + InputEvent()->GetVZEROData()->GetTriggerChargeC();
  fh3PtLeadJet1PatchEnergyVZEROAmp->Fill(ptLeadJet1,fMaxPatchEnergy,VZEROAmp);
  fh3PtLeadJet1RawPatchEnergyVZEROAmp->Fill(ptLeadJet1,fMaxPatchADCEnergy,VZEROAmp);

  delete nJetsArr;

  //Cells - some QA plots
  if(fCaloCells) {
    const Short_t nCells   = fCaloCells->GetNumberOfCells();

    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = fCaloCells->GetCellNumber(iCell);
      Double_t cellE = fCaloCells->GetCellAmplitude(cellId);
      Double_t cellT = fCaloCells->GetCellTime(cellId);
      TVector3 pos;
      fGeom->GetGlobal(cellId, pos);
      TLorentzVector lv(pos,cellE);
      Double_t cellEta = lv.Eta();
      Double_t cellPhi = lv.Phi();
      if(cellPhi<0.) cellPhi+=TMath::TwoPi();
      if(cellPhi>TMath::TwoPi()) cellPhi-=TMath::TwoPi();

      AliDebug(2,Form("cell energy = %f  time = %f",cellE,cellT*1e9));
      fh2CellEnergyVsTime->Fill(cellE,cellT*1e9);
      fh3EEtaPhiCell->Fill(cellE,cellEta,cellPhi);
      fh2ECellVsCent->Fill(fCent,cellE);
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTriggerPatchJetMatch::Run() {
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  fhTriggerbit->Fill(0.5,GetCollisionCandidates());

  // just in case..
  fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

  //Check if event is selected (vertex & pile-up)
  //if(!SelectEvent()) return kFALSE;

  // when we have the patch object to match, peform analysis
  if(fTriggerPatchInfo) FillTriggerPatchHistos();

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalTriggerPatchJetMatch::Terminate(Option_t *) {
  // Called once at the end of the analysis.
}
//________________________________________________________________________
Double_t AliAnalysisTaskEmcalTriggerPatchJetMatch::GetZ(const AliVParticle *trk, const AliEmcalJet *jet) const {  
  // Get Z of constituent trk
  return GetZ(trk->Px(),trk->Py(),trk->Pz(),jet->Px(),jet->Py(),jet->Pz());
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalTriggerPatchJetMatch::GetZ(Double_t trkPx, Double_t trkPy, Double_t trkPz, Double_t jetPx, Double_t jetPy, Double_t jetPz) const {
  // Get the z of a constituent inside of a jet
  Double_t pJetSq = jetPx*jetPx+jetPy*jetPy+jetPz*jetPz;
  if(pJetSq>0.)
    return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/pJetSq;
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f",GetName(), pJetSq)); 
    return 0;
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalTriggerPatchJetMatch::GetLeadingCellId(const AliVCluster *clus) const {
  //Get energy of leading cell in cluster
  if(!fCaloCells) return -1;

  Double_t emax = -1.;
  Int_t iCellAbsIdMax = -1;
  Int_t nCells = clus->GetNCells();
  for(Int_t i = 0; i<nCells; i++) {
    Int_t absId = clus->GetCellAbsId(i);
    Double_t cellE = fCaloCells->GetCellAmplitude(absId);
    if(cellE>emax) {
      emax = cellE;
      iCellAbsIdMax = absId;
    }
  }
  return iCellAbsIdMax;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalTriggerPatchJetMatch::GetEnergyLeadingCell(const AliVCluster *clus) const {
  //Get energy of leading cell in cluster
  if(!fCaloCells) return -1.;

  Int_t absID = GetLeadingCellId(clus);
  if(absID>-1) return fCaloCells->GetCellAmplitude(absID);
  else return -1.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalTriggerPatchJetMatch::GetECross(Int_t absID) const {
  //Get Ecross = sum of energy of neighbouring cells (using uncalibrated energy)
  if(!fCaloCells) return -1.;

  Double_t ecross = -1.;

  Int_t absID1 = -1;
  Int_t absID2 = -1;
  Int_t absID3 = -1;
  Int_t absID4 = -1;

  Int_t imod = -1, iphi =-1, ieta=-1, iTower = -1, iIphi = -1, iIeta = -1;
  fGeom->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

  if( iphi < AliEMCALGeoParams::fgkEMCALRows-1)
    absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if( iphi > 0 )
    absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  if( ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) ) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  }
  else if( ieta == 0 && imod%2 ) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  }
  else  {
    if( ieta < AliEMCALGeoParams::fgkEMCALCols-1 )
      absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if( ieta > 0 )
      absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
  }

  Double_t ecell1 = fCaloCells->GetCellAmplitude(absID1);
  Double_t ecell2 = fCaloCells->GetCellAmplitude(absID2);
  Double_t ecell3 = fCaloCells->GetCellAmplitude(absID3);
  Double_t ecell4 = fCaloCells->GetCellAmplitude(absID4);

  ecross = ecell1+ecell2+ecell3+ecell4;

  return ecross;
}

//_________________________________________________________________________
Float_t AliAnalysisTaskEmcalTriggerPatchJetMatch::RelativeEP(Double_t objAng, Double_t EPAng) const {
  // function to calculate angle between object and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi = EPAng - objAng;

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() )
    dphi = dphi + 1*TMath::Pi();
  if( dphi>1*TMath::Pi())
    dphi = dphi - 1*TMath::Pi();

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  }

  return dphi;   // dphi in [0, Pi/2]
}

TH1* AliAnalysisTaskEmcalTriggerPatchJetMatch::FillTriggerPatchQA(TH1* h, UInt_t trig, AliEmcalTriggerPatchInfo *fPatch) {
  // check and fill a QA histogram
  if(fPatch->IsLevel0()) h->Fill(1);
  if(fPatch->IsJetLow()) h->Fill(2);
  if(fPatch->IsJetHigh()) h->Fill(3);
  if(fPatch->IsGammaLow()) h->Fill(4);
  if(fPatch->IsGammaHigh()) h->Fill(5);
  if(fPatch->IsMainTrigger()) h->Fill(6);
  if(fPatch->IsJetLowSimple()) h->Fill(7);
  if(fPatch->IsJetHighSimple()) h->Fill(8);
  if(fPatch->IsGammaLowSimple()) h->Fill(9);
  if(fPatch->IsGammaHighSimple()) h->Fill(10);
  if(fPatch->IsMainTriggerSimple()) h->Fill(11);
  if(fPatch->IsOfflineSimple()) h->Fill(12);
  if(fPatch->IsRecalcJet()) h->Fill(13);
  if(fPatch->IsRecalcGamma()) h->Fill(14);
  if(trig & AliVEvent::kEMCEJE) h->Fill(19);
  if(trig & AliVEvent::kEMCEGA) h->Fill(20);

  h->GetXaxis()->SetBinLabel(1, "Level0");   
  h->GetXaxis()->SetBinLabel(2, "JetLow");
  h->GetXaxis()->SetBinLabel(3, "JetHigh");
  h->GetXaxis()->SetBinLabel(4, "GammaLow");
  h->GetXaxis()->SetBinLabel(5, "GammaHigh");
  h->GetXaxis()->SetBinLabel(6, "MainTrigger");
  h->GetXaxis()->SetBinLabel(7, "JetLowSimple");
  h->GetXaxis()->SetBinLabel(8, "JetHighSimple");
  h->GetXaxis()->SetBinLabel(9, "GammaLowSimple");
  h->GetXaxis()->SetBinLabel(10, "GammaHighSimple");
  h->GetXaxis()->SetBinLabel(11, "MainTriggerSimple");
  h->GetXaxis()->SetBinLabel(12, "OfflineSimple");
  h->GetXaxis()->SetBinLabel(13, "RecalcJet");
  h->GetXaxis()->SetBinLabel(14, "RecalcGamma");
  h->GetXaxis()->SetBinLabel(15, "");
  h->GetXaxis()->SetBinLabel(16, "");
  h->GetXaxis()->SetBinLabel(17, "");
  h->GetXaxis()->SetBinLabel(18, "");
  h->GetXaxis()->SetBinLabel(19, "kEMCEJE");
  h->GetXaxis()->SetBinLabel(20, "kEMCEGA");

  // set x-axis labels vertically
  //h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}

Bool_t AliAnalysisTaskEmcalTriggerPatchJetMatch::CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const {
  Bool_t hasfound = kFALSE;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *mypatch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    Double_t etamin = TMath::Min(mypatch->GetEtaMin(), mypatch->GetEtaMax()),
    etamax = TMath::Max(mypatch->GetEtaMin(), mypatch->GetEtaMax()),
    phimin = TMath::Min(mypatch->GetPhiMin(), mypatch->GetPhiMax()),
    phimax = TMath::Max(mypatch->GetPhiMin(), mypatch->GetPhiMax());
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      hasfound = kTRUE;
      break;
    }
  }
  return hasfound;
}

TH1* AliAnalysisTaskEmcalTriggerPatchJetMatch::FillEventTriggerQA(TH1* h, UInt_t trig) {
  // check and fill a Event Selection QA histogram for different trigger selections after cuts
  if(trig == 0) h->Fill(1);
  if(trig & AliVEvent::kAny) h->Fill(2);
  if(trig & AliVEvent::kAnyINT) h->Fill(3);
  if(trig & AliVEvent::kMB) h->Fill(4);
  if(trig & AliVEvent::kINT7) h->Fill(5);
  if(trig & AliVEvent::kEMC1) h->Fill(6);
  if(trig & AliVEvent::kEMC7) h->Fill(7);
  if(trig & AliVEvent::kEMC8) h->Fill(8);
  if(trig & AliVEvent::kEMCEJE) h->Fill(9);
  if(trig & AliVEvent::kEMCEGA) h->Fill(10);
  if(trig & AliVEvent::kCentral) h->Fill(11);
  if(trig & AliVEvent::kSemiCentral) h->Fill(12);
  if(trig & AliVEvent::kINT8) h->Fill(13);

  if(trig & (AliVEvent::kEMCEJE | AliVEvent::kMB)) h->Fill(14);
  if(trig & (AliVEvent::kEMCEGA | AliVEvent::kMB)) h->Fill(15);
  if(trig & (AliVEvent::kAnyINT | AliVEvent::kMB)) h->Fill(16);

  if(trig & (AliVEvent::kEMCEJE & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(17);
  if(trig & (AliVEvent::kEMCEGA & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(18);
  if(trig & (AliVEvent::kAnyINT & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(19);

  // label bins of the analysis trigger selection summary
  h->GetXaxis()->SetBinLabel(1, "no trigger");
  h->GetXaxis()->SetBinLabel(2, "kAny");
  h->GetXaxis()->SetBinLabel(3, "kAnyINT");
  h->GetXaxis()->SetBinLabel(4, "kMB");
  h->GetXaxis()->SetBinLabel(5, "kINT7");
  h->GetXaxis()->SetBinLabel(6, "kEMC1");
  h->GetXaxis()->SetBinLabel(7, "kEMC7");
  h->GetXaxis()->SetBinLabel(8, "kEMC8");
  h->GetXaxis()->SetBinLabel(9, "kEMCEJE");
  h->GetXaxis()->SetBinLabel(10, "kEMCEGA");
  h->GetXaxis()->SetBinLabel(11, "kCentral");
  h->GetXaxis()->SetBinLabel(12, "kSemiCentral");
  h->GetXaxis()->SetBinLabel(13, "kINT8");
  h->GetXaxis()->SetBinLabel(14, "kEMCEJE or kMB");
  h->GetXaxis()->SetBinLabel(15, "kEMCEGA or kMB");
  h->GetXaxis()->SetBinLabel(16, "kAnyINT or kMB");
  h->GetXaxis()->SetBinLabel(17, "kEMCEJE & (kMB or kCentral or kSemiCentral)");
  h->GetXaxis()->SetBinLabel(18, "kEMCEGA & (kMB or kCentral or kSemiCentral)");
  h->GetXaxis()->SetBinLabel(19, "kAnyINT & (kMB or kCentral or kSemiCentral)");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
