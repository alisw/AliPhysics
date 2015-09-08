#include "AliTestAD.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include <TList.h>
#include <TH2.h>

//====================================================================
namespace {
  AliInputEventHandler* GetInputHandler()
  {
    AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
    if (!am) {
      Warning("GetPhysicsSelection", "No analysis manager");
      return 0;
    }
    
    // Get the input handler - should always be there 
    AliInputEventHandler* ih = 
      dynamic_cast<AliInputEventHandler*>(am->GetInputEventHandler());
    if (!ih) { 
      Warning("GetPhysicsSelection", "No input handler");
      return 0;
    }
    return ih;
  }
  AliPhysicsSelection* GetPhysicsSelection()
  {
    AliInputEventHandler* ih = GetInputHandler();
    if (!ih) return 0;
    
    // Get the physics selection - should always be there 
    AliPhysicsSelection* ps = 
      dynamic_cast<AliPhysicsSelection*>(ih->GetEventSelection());
    if (!ps) {
      Warning("GetPhysicsSelection", "No physics selection");
      return 0;
    }
    return ps;
  }
}

//====================================================================
AliTestAD::AliTestAD(const char*)
  : AliAnalysisTaskSE("ad"),
    fList(0),
    fAD(),
    fIPz0(0),
    fIPz1(0),
    fIPz2(0),
    fIPzSat0(0),
    fIPzSat1(0),
    fIPzSat2(0),
    fIPzSatP(0),
    fIsInit()
{
  DefineOutput(1, TList::Class());
}
//--------------------------------------------------------------------
void AliTestAD::Connect()
{
  // Get the analysis manager and connect this task 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  // Add to the manager 
  mgr->AddTask(this);

  // Connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());

  // Connect output 
  AliAnalysisDataContainer* sumCon = 
    mgr->CreateContainer("ad", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,
			 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(this, 1, sumCon);
}
//--------------------------------------------------------------------
void AliTestAD::UserCreateOutputObjects()
{
  // Define output list 
  fList = new TList;
  fList->SetName(GetName());
  fList->SetOwner();

  Int_t    nSatIPz   = 4000;
  Double_t maxSatIPz = 400;
  // MAke correlation histogram 
  fIPz0 = new TH2D("ipZ0", "AD vs SPD IP_{#it{z}} (>0)",
		   100, -50, 50, nSatIPz, -maxSatIPz, +maxSatIPz);
  fIPz0->SetDirectory(0);
  fIPz0->SetXTitle("SPD IP_{#it{z}} [cm]");
  fIPz0->SetYTitle("AD IP_{#it{z}} [cm]");
  fIPz0->SetZTitle("Events");
  fIPz0->SetDirectory(0);

  fIPz1 = static_cast<TH2*>(fIPz0->Clone("ipZ1"));
  fIPz1->SetTitle("AD vs SPD IP_{#it{z}} (>1)");
  fIPz1->SetDirectory(0);

  fIPz2 = static_cast<TH2*>(fIPz0->Clone("ipZ2"));
  fIPz2->SetTitle("AD vs SPD IP_{#it{z}} (>2)");
  fIPz2->SetDirectory(0);
  
  fIPzSat0 = new TH1D("ipZSat0", "Satellite IP_{#it{z}} (SPD>0 veto)",
		      nSatIPz, -maxSatIPz, +maxSatIPz);
  fIPzSat0->SetXTitle("AD IP_{#it{z}} [cm]");
  fIPzSat0->SetYTitle("Events");
  fIPzSat0->SetFillStyle(3001);
  fIPzSat0->SetFillColor(kRed-3);
  fIPzSat0->SetLineColor(kRed-3);
  fIPzSat0->SetDirectory(0);

  fIPzSat1 = static_cast<TH1*>(fIPzSat0->Clone("ipZSat1"));
  fIPzSat1->SetTitle("Satellite IP_{#it{z}} (SPD>1 veto)");
  fIPzSat1->SetFillColor(kGreen-3);
  fIPzSat1->SetLineColor(kGreen-3);
  fIPzSat1->SetDirectory(0);

  fIPzSat2 = static_cast<TH1*>(fIPzSat0->Clone("ipZSat2"));
  fIPzSat2->SetTitle("Satellite IP_{#it{z}} (SPD>2 veto)");
  fIPzSat2->SetFillColor(kBlue-3);
  fIPzSat2->SetLineColor(kBlue-3);
  fIPzSat2->SetDirectory(0);
  

  fIPzSatP = static_cast<TH1*>(fIPzSat0->Clone("ipZSatP"));
  fIPzSatP->SetTitle("Satellite IP_{#it{z}} (SPD pile-up veto)");
  fIPzSatP->SetFillColor(kMagenta-3);
  fIPzSatP->SetLineColor(kMagenta-3);
  fIPzSatP->SetDirectory(0);
    
  fList->Add(fIPz0);
  fList->Add(fIPz1);
  fList->Add(fIPz2);
  fList->Add(fIPzSat0); 
  fList->Add(fIPzSat1); 
  fList->Add(fIPzSat2); 
  fList->Add(fIPzSatP); 

  PostData(1,fList);
}		 
//--------------------------------------------------------------------
void AliTestAD::UserExec(Option_t*)
{
  // Check for an event 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return;
  }

  // IF we have a physis selection, check that we have a trigger
  UInt_t trg = (GetPhysicsSelection() ?
		GetInputHandler()->IsEventSelected() :
		AliVEvent::kMB);
  if (!trg) return;
  
  // If not initialized yet, do so here. 
  if (!fIsInit) {
    fAD.SetupForData(fList, "", false);
    fIsInit = true;
  }

  // Let sub-procedure test the event 
  fAD.Process(esd);

  // Get and check SPD ip 
  Double_t spdZ0 = AliDisplacedVertexSelectionAD::kInvalidVtxZ;
  Double_t spdZ1 = AliDisplacedVertexSelectionAD::kInvalidVtxZ;
  Double_t spdZ2 = AliDisplacedVertexSelectionAD::kInvalidVtxZ;
  const AliESDVertex* ip = esd->GetPrimaryVertexSPD();
  if (ip && ip->GetZRes() <= 0.2) {
    if (ip->GetNContributors() >= 1) spdZ0 = ip->GetZ();
    if (ip->GetNContributors() >= 2) spdZ1 = ip->GetZ();
    if (ip->GetNContributors() >= 3) spdZ2 = ip->GetZ();
  }

  // Fill Correlation of vertex 
  if (spdZ0 !=  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
    fIPz0->Fill(spdZ0, fAD.GetIPz());
  if (spdZ1 !=  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
    fIPz1->Fill(spdZ1, fAD.GetIPz());
  if (spdZ2 !=  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
    fIPz2->Fill(spdZ2, fAD.GetIPz());

  Bool_t pileup = esd->IsPileupFromSPDInMultBins();
    
  if (fAD.IsSatellite()) {
    if (spdZ0 ==  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
      fIPzSat0->Fill(fAD.GetIPz());
    if (spdZ1 ==  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
      fIPzSat1->Fill(fAD.GetIPz());
    if (spdZ2 ==  AliDisplacedVertexSelectionAD::kInvalidVtxZ)
      fIPzSat2->Fill(fAD.GetIPz());
    if (!pileup)
      fIPzSatP->Fill(fAD.GetIPz());
  }
  
  // Post to output 
  PostData(1,fList);
}
//--------------------------------------------------------------------
//
// EOF
// 
