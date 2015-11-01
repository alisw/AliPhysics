/**
 * @file   AliCopyHeaderTask.cxx
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:59:32 2011
 * 
 * @brief  Task to copy ESD header to AOD 
 * 
 * @ingroup pwglf_forward_tasks 
 */

#include "AliCopyHeaderTask.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "TFile.h"
#include "AliEventplane.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"

ClassImp(AliCopyHeaderTask)
#if 0 
; // for emacs - do not remove 
#endif


//____________________________________________________________________
void
AliCopyHeaderTask::UserExec(Option_t*)
{
  // 
  // Called at every event 
  //
  // Copies information from ESD header to AOD header
  //

  // --- Get the I/O objects -----------------------------------------
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(AODEvent());
  if (!esd) { 
    AliWarning("Missing ESD event");
    return;
  }
  if (!aod) { 
    AliWarning("Missing AOD event");
    return;
  }

  // --- Load the data -----------------------------------------------
  LoadBranches();

  // --- Get or create the header ------------------------------------
  AliAODHeader* aodHeader = dynamic_cast<AliAODHeader*>(aod->GetHeader());
  if(!aodHeader) AliFatal("Not a standard AOD");
  if (!aodHeader) { 
    AliWarning("Missing AOD header");
    aodHeader = new AliAODHeader(esd->GetRunNumber(),
				 esd->GetBunchCrossNumber(),
				 esd->GetOrbitNumber(),
				 esd->GetPeriodNumber());
    aod->AddHeader(aodHeader);
  }

  // --- Various event/run stuff -------------------------------------
  aodHeader->SetRunNumber(esd->GetRunNumber());
  aodHeader->SetOfflineTrigger(fInputHandler->IsEventSelected());
  aodHeader->SetBunchCrossNumber(esd->GetBunchCrossNumber());
  aodHeader->SetOrbitNumber(esd->GetOrbitNumber());
  aodHeader->SetPeriodNumber(esd->GetPeriodNumber());
  aodHeader->SetEventType(esd->GetEventType());
  aodHeader->SetEventNumberESDFile(esd->GetHeader()->GetEventNumberInFile());

  // --- Centrality --------------------------------------------------
  if(esd->GetCentrality())
    aodHeader->SetCentrality(esd->GetCentrality());
  else
    aodHeader->SetCentrality(0);

  // --- Trigger -----------------------------------------------------
  aodHeader->SetFiredTriggerClasses(esd->GetFiredTriggerClasses());
  aodHeader->SetTriggerMask(esd->GetTriggerMask()); 
  aodHeader->SetTriggerCluster(esd->GetTriggerCluster());
  aodHeader->SetL0TriggerInputs(esd->GetHeader()->GetL0TriggerInputs());    
  aodHeader->SetL1TriggerInputs(esd->GetHeader()->GetL1TriggerInputs());    
  aodHeader->SetL2TriggerInputs(esd->GetHeader()->GetL2TriggerInputs());    

  // --- Magnetic field, ZDC signal ----------------------------------
  aodHeader->SetMagneticField(esd->GetMagneticField());
  aodHeader->SetMuonMagFieldScale(esd->GetCurrentDip()/6000.);
  aodHeader->SetZDCN1Energy(esd->GetZDCN1Energy());
  aodHeader->SetZDCP1Energy(esd->GetZDCP1Energy());
  aodHeader->SetZDCN2Energy(esd->GetZDCN2Energy());
  aodHeader->SetZDCP2Energy(esd->GetZDCP2Energy());
  aodHeader->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));

  // --- Interacting beams -------------------------------------------
  AliESDHeader* esdHeader = esd->GetHeader();
  if (esdHeader) { 
    aodHeader->SetIRInt2InteractionMap(esdHeader->GetIRInt2InteractionMap());
    aodHeader->SetIRInt1InteractionMap(esdHeader->GetIRInt1InteractionMap());
  }

  // --- ESD file name -----------------------------------------------
  TTree* tree = fInputHandler->GetTree();
  if (tree) {
    TFile* file = tree->GetCurrentFile();
    if (file) aodHeader->SetESDFileName(file->GetName());
  }

  // --- TPC event plane ---------------------------------------------
  AliEventplane* ep = esd->GetEventplane();
  if (ep) aodHeader->SetEventplane(ep);

  // --- Copy primary vertices ---------------------------------------
  CopyVertex(*aod, esd->GetPrimaryVertex(),    AliAODVertex::kPrimary);
  CopyVertex(*aod, esd->GetPrimaryVertexSPD(), AliAODVertex::kMainSPD);
  CopyVertex(*aod, esd->GetPrimaryVertexTPC(), AliAODVertex::kMainTPC);
  
  // --- Loop over pile-ups vertices ---------------------------------
  for (Int_t i = 0; i < esd->GetNumberOfPileupVerticesSPD(); i++) 
    CopyVertex(*aod, esd->GetPileupVertexSPD(i), AliAODVertex::kPileupSPD);
  for (Int_t i = 0; i < esd->GetNumberOfPileupVerticesTracks(); i++) 
    CopyVertex(*aod, esd->GetPileupVertexTracks(i),AliAODVertex::kPileupTracks);


  // --- Reference multiplicities ------------------------------------
  if (fCalculateRefMult) {
    AliESDtrackCuts::MultEstTrackType estType =
      esd->GetPrimaryVertexTracks()->GetStatus()
      ? AliESDtrackCuts::kTrackletsITSTPC
      : AliESDtrackCuts::kTracklets;
    Int_t mult05 = AliESDtrackCuts::GetReferenceMultiplicity(esd,estType,0.5);
    Int_t mult08 = AliESDtrackCuts::GetReferenceMultiplicity(esd,estType,0.8);
    aodHeader->SetRefMultiplicityComb05(mult05);
    aodHeader->SetRefMultiplicityComb08(mult08);
  }
  // --- Copy multiplicity object ------------------------------------
  // AliAnalysisTaskESDfilter::ConvertTracklets(const AliESDEvent& esd)
  AliMultiplicity* mul       = esd->GetMultiplicity();
  AliAODTracklets* tracklets = (aod->GetTracklets());
  if (mul && mul->GetNumberOfTracklets() > 0 && tracklets) {
    Int_t nTracklets =  mul->GetNumberOfTracklets();
    tracklets->CreateContainer(nTracklets);
    for (Int_t i = 0; i < nTracklets; i++) {
      tracklets->SetTracklet(i,
			     mul->GetTheta(i),
			     mul->GetPhi(i),
			     mul->GetDeltaPhi(i),
			     mul->GetLabel(i, 0),
			     mul->GetLabel(i, 1));
    }
    tracklets->SetFiredChipMap      (mul->GetFiredChipMap());
    tracklets->SetFastOrFiredChipMap(mul->GetFastOrFiredChipMap());
    tracklets->SetFiredChips        (0, mul->GetNumberOfFiredChips(0));
    tracklets->SetFiredChips        (1, mul->GetNumberOfFiredChips(1));
  }
}

//____________________________________________________________________
void
AliCopyHeaderTask::CopyVertex(AliAODEvent&        aod, 
			      const AliESDVertex* vtx, 
			      Int_t               type) 
{
  if (!vtx) return;

  // --- Get array of V0's from AOD ----------------------------------
  TClonesArray* arr = aod.GetVertices();
  if (!arr) return;

  // --- Now get some stuff vertex to copy to AOD --------------------
  Int_t    n     = arr->GetEntriesFast(); 
  Double_t pos[] = { 0., 0., 0. };
  Double_t cov[] = { 0., 0., 0., 0., 0., 0. }; 
  Double_t chi2  = vtx->GetChi2toNDF();
  vtx->GetXYZ(pos);
  vtx->GetCovMatrix(cov);

  // --- Create new AOD vertex and set stuff -------------------------
  AliAODVertex* out = new((*arr)[n]) AliAODVertex(pos, cov, chi2, 0, -1, type);
  out->SetName(vtx->GetName());
  out->SetTitle(vtx->GetTitle());
  out->SetBC(vtx->GetBC());

  // --- If from a track vertex, make sure we set the contributors ---
  TString tit(out->GetTitle());
  if (!tit.Contains("VertexerTracks")) 
    out->SetNContributors(vtx->GetNContributors());
}

//____________________________________________________________________
void
AliCopyHeaderTask::Terminate(Option_t*)
{
  // Called at the end of the job  - does nothing 
}

//____________________________________________________________________
Bool_t
AliCopyHeaderTask::Connect()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("Connect", "No analysis manager to connect to.");
    return false;
  }   

  // Add to the manager 
  mgr->AddTask(this);
  
  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());

  return true;
}

//
// EOF
//
