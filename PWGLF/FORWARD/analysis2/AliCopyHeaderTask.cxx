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

ClassImp(AliCopyHeaderTask)
#if 0 
; // for emacs - do not remove 
#endif


void
AliCopyHeaderTask::UserExec(Option_t*)
{
  // 
  // Called at every event 
  //
  // Copies information from ESD header to AOD header
  // 
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

  LoadBranches();

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

  aodHeader->SetRunNumber(esd->GetRunNumber());
  aodHeader->SetOfflineTrigger(fInputHandler->IsEventSelected());
  aodHeader->SetBunchCrossNumber(esd->GetBunchCrossNumber());
  aodHeader->SetOrbitNumber(esd->GetOrbitNumber());
  aodHeader->SetPeriodNumber(esd->GetPeriodNumber());
  aodHeader->SetEventType(esd->GetEventType());
  aodHeader->SetEventNumberESDFile(esd->GetHeader()->GetEventNumberInFile());
  if(esd->GetCentrality())
    aodHeader->SetCentrality(esd->GetCentrality());
  else
    aodHeader->SetCentrality(0);

  aodHeader->SetFiredTriggerClasses(esd->GetFiredTriggerClasses());
  aodHeader->SetTriggerMask(esd->GetTriggerMask()); 
  aodHeader->SetTriggerCluster(esd->GetTriggerCluster());
  aodHeader->SetL0TriggerInputs(esd->GetHeader()->GetL0TriggerInputs());    
  aodHeader->SetL1TriggerInputs(esd->GetHeader()->GetL1TriggerInputs());    
  aodHeader->SetL2TriggerInputs(esd->GetHeader()->GetL2TriggerInputs());    
  
  aodHeader->SetMagneticField(esd->GetMagneticField());
  aodHeader->SetMuonMagFieldScale(esd->GetCurrentDip()/6000.);
  aodHeader->SetZDCN1Energy(esd->GetZDCN1Energy());
  aodHeader->SetZDCP1Energy(esd->GetZDCP1Energy());
  aodHeader->SetZDCN2Energy(esd->GetZDCN2Energy());
  aodHeader->SetZDCP2Energy(esd->GetZDCP2Energy());
  aodHeader->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));

  AliESDHeader* esdHeader = esd->GetHeader();
  if (esdHeader) { 
    aodHeader->SetIRInt2InteractionMap(esdHeader->GetIRInt2InteractionMap());
    aodHeader->SetIRInt1InteractionMap(esdHeader->GetIRInt1InteractionMap());
  }
  
  TTree* tree = fInputHandler->GetTree();
  if (tree) {
    TFile* file = tree->GetCurrentFile();
    if (file) aodHeader->SetESDFileName(file->GetName());
  }

  AliEventplane* ep = esd->GetEventplane();
  if (ep) aodHeader->SetEventplane(ep);

  // Copy primary vertices
  CopyVertex(*aod, esd->GetPrimaryVertex(),    AliAODVertex::kPrimary);
  CopyVertex(*aod, esd->GetPrimaryVertexSPD(), AliAODVertex::kMainSPD);
  CopyVertex(*aod, esd->GetPrimaryVertexTPC(), AliAODVertex::kMainTPC);
  
  // Loop over pile-ups 
  for (Int_t i = 0; i < esd->GetNumberOfPileupVerticesSPD(); i++) 
    CopyVertex(*aod, esd->GetPileupVertexSPD(i), AliAODVertex::kPileupSPD);
  for (Int_t i = 0; i < esd->GetNumberOfPileupVerticesTracks(); i++) 
    CopyVertex(*aod, esd->GetPileupVertexTracks(i),AliAODVertex::kPileupTracks);
    
}

void
AliCopyHeaderTask::CopyVertex(AliAODEvent&        aod, 
			      const AliESDVertex* vtx, 
			      Int_t               type) 
{
  if (!vtx) return;

  TClonesArray* arr = aod.GetVertices();
  if (!arr) return;

  Int_t    n     = arr->GetEntriesFast(); 
  Double_t pos[] = { 0., 0., 0. };
  Double_t cov[] = { 0., 0., 0., 0., 0., 0. }; 
  Double_t chi2  = vtx->GetChi2toNDF();
  vtx->GetXYZ(pos);
  vtx->GetCovMatrix(cov);
  
  AliAODVertex* out = new((*arr)[n]) AliAODVertex(pos, cov, chi2, 0, -1, type);
  out->SetName(vtx->GetName());
  out->SetTitle(vtx->GetTitle());
  out->SetBC(vtx->GetBC());
  TString tit(out->GetTitle());
  if (!tit.Contains("VertexerTracks")) 
    out->SetNContributors(vtx->GetNContributors());
}

void
AliCopyHeaderTask::Terminate(Option_t*)
{
  // Called at the end of the job  - does nothing 
}

//
// EOF
//
