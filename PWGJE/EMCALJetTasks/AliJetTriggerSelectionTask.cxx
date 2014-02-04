// $Id$
//
// Jet trigger selection task.
//
// Author: S.Aiola

#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>

#include "AliVEvent.h"
#include "AliVCluster.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"

#include "AliJetTriggerSelectionTask.h"

ClassImp(AliJetTriggerSelectionTask)

//________________________________________________________________________
AliJetTriggerSelectionTask::AliJetTriggerSelectionTask() : 
  AliAnalysisTaskEmcalJet("AliJetTriggerSelectionTask", kFALSE),
  fEnergyThreshold(0),
  fMaxDistance2(0.0225),
  fTriggerBits(AliVEvent::kEMCEJE),
  fTaskSettingsOk(kFALSE),
  fNTriggers(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0)
{
  // Default constructor.
  
  for (Int_t i = 0; i < 999; i++) {
    fTrigPos[i][0] = -999;
    fTrigPos[i][1] = -999;
  }
}

//________________________________________________________________________
AliJetTriggerSelectionTask::AliJetTriggerSelectionTask(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kFALSE),
  fEnergyThreshold(0),
  fMaxDistance2(0.0225),
  fTriggerBits(AliVEvent::kEMCEJE),
  fTaskSettingsOk(kFALSE),
  fNTriggers(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 999; i++) {
    fTrigPos[i][0] = -999;
    fTrigPos[i][1] = -999;
  }
}

//________________________________________________________________________
void AliJetTriggerSelectionTask::ExecOnce()
{
  // Initialize the task.

  AliAnalysisTaskEmcalJet::ExecOnce();

  fTaskSettingsOk = kTRUE;

  fVZERO = InputEvent()->GetVZEROData();
  if (!fVZERO) {
    AliError(Form("%s: AliVVZERO not available, task will not be executed!",GetName()));
    fTaskSettingsOk = kFALSE;
  }
  
  if (GetClusterArray() == 0) {
    AliError(Form("%s: No cluster collection provided, task will not be executed!",GetName()));
    fTaskSettingsOk = kFALSE;
  }

  if (GetJetArray() == 0) {
    AliError(Form("%s: No jet collection provided, task will not be executed!",GetName()));
    fTaskSettingsOk = kFALSE;
  }
 
  if (!fEnergyThreshold) {
    AliError(Form("%s: No threshold function provided, task will not be executed!",GetName()));
    fTaskSettingsOk = kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliJetTriggerSelectionTask::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  if (fVZERO) {
    fV0ATotMult = fVZERO->GetMTotV0A();
    fV0CTotMult = fVZERO->GetMTotV0C();
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetTriggerSelectionTask::Run()
{
  // Run the analysis.

  if (!fTaskSettingsOk) return kFALSE;

  FindTriggers();
  SelectJets();
  
  return kTRUE;
}

//________________________________________________________________________
void AliJetTriggerSelectionTask::FindTriggers()
{
  fNTriggers = 0;

  Int_t nclusters = GetNClusters();

  Double_t th = fEnergyThreshold->Eval(fV0ATotMult+fV0CTotMult);

  for (Int_t i = 0; i < nclusters; i++) {
    if (fNTriggers >= 999) {
      AliError("More than 999 triggers found!");
      break;
    }

    AliVCluster *cluster = GetAcceptClusterFromArray(i);
    if (!cluster)
      continue;

    //Printf("Cluster energy=%.3f, th=%.3f",cluster->E(), th);

    if (cluster->E() > th) {
      TLorentzVector vect;
      cluster->GetMomentum(vect,fVertex);
      fTrigPos[fNTriggers][0] = vect.Eta();
      fTrigPos[fNTriggers][1] = vect.Phi();
      fNTriggers++;
    }
  }

  AliDebug(2,Form("%s: %d triggers found among %d candidates (cent=%.1f, mult=%.1f, th=%.2f)!",GetName(),fNTriggers,nclusters,fCent,fV0ATotMult+fV0CTotMult,th));
}

//________________________________________________________________________
void AliJetTriggerSelectionTask::SelectJets()
{
  for (Int_t c = 0; c < fJetCollArray.GetEntriesFast(); c++) {

    Int_t njets = GetNJets(c);

    for (Int_t i = 0; i < njets; i++) {
      AliEmcalJet *jet = GetAcceptJetFromArray(i,c);
      if (IsTriggerJet(jet)) jet->AddTrigger(fTriggerBits);
    }
  }
}

//________________________________________________________________________
Bool_t AliJetTriggerSelectionTask::IsTriggerJet(AliEmcalJet *jet)
{
  if (!jet) return kFALSE;

  for (Int_t i = 0; i < fNTriggers; i++) {
    Double_t deta = jet->Eta() - fTrigPos[i][0];
    Double_t dphi = jet->Phi() - fTrigPos[i][1];
    
    Double_t d2 = deta * deta + dphi * dphi;

    if (d2 < fMaxDistance2)
      return kTRUE;
  }
  
  return kFALSE;
}
