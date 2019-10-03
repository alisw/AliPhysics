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
#include "AliMultSelection.h"
#include "AliAODHandler.h"

ClassImp(AliCopyHeaderTask)
#if 0 
; // for emacs - do not remove 
#endif


//____________________________________________________________________
AliCopyHeaderTask&
AliCopyHeaderTask::operator=(const AliCopyHeaderTask& other) 
{
  if (this == &other) return *this;
  AliAnalysisTaskSE::operator=(other);
  fCalculateRefMult = other.fCalculateRefMult;
  fCopyCentrality   = other.fCopyCentrality;
  fCopyTracklets    = other.fCopyTracklets;
  fCopyV0           = other.fCopyV0;
  fCopyAD           = other.fCopyAD;
  fCopyZDC          = other.fCopyZDC;
  if (fMultSelection) {
    delete fMultSelection;
    fMultSelection = 0;
  }
  if (other.fMultSelection)
    fMultSelection = new AliMultSelection(*other.fMultSelection);
  return *this;
}

//____________________________________________________________________
AliCopyHeaderTask::AliCopyHeaderTask(const AliCopyHeaderTask& other) 
  : AliAnalysisTaskSE(other),
    fCalculateRefMult(other.fCalculateRefMult),
    fCopyCentrality(other.fCopyCentrality),
    fCopyTracklets(other.fCopyTracklets),
    fCopyV0(other.fCopyV0),
    fCopyAD(other.fCopyAD),
    fCopyZDC(other.fCopyZDC),
    fMultSelection(0)
{
  if (other.fMultSelection) 
    fMultSelection = new AliMultSelection(*other.fMultSelection);
}

//____________________________________________________________________
void
AliCopyHeaderTask::SetCopyOptions(const TString& what)
{
  TObjArray* tokens = what.Tokenize(",");
  TObject*   token  = 0;
  TIter      next(tokens);
  while ((token = next())) {
    TString opt(token->GetName());
    Bool_t  enable = true;
    Int_t   rem    = 0;
    if      (opt.BeginsWith("+")) { rem = 1; enable = true;  }
    else if (opt.BeginsWith("-")) { rem = 1; enable = false; }
    if (rem > 0) opt.Remove(0,1);

    opt.ToLower();

    if      (opt.BeginsWith("cent")) SetCopyCentrality(enable);
    else if (opt.BeginsWith("trac")) SetCopyTracklets(enable);
    else if (opt.BeginsWith("v0") || opt.BeginsWith("vzero"))
      SetCopyV0(enable);
    else if (opt.BeginsWith("ad"))   SetCopyAD(enable);
    else if (opt.BeginsWith("zdc"))  SetCopyZDC(enable);
    else if (opt.BeginsWith("ref")) SetCalculateRefMult(enable);
  }
  tokens->Delete();
}
    

  
//____________________________________________________________________
void
AliCopyHeaderTask::UserCreateOutputObjects()
{
  if (!fCopyCentrality) return;
  
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) {
    AliWarning("No AOD output handler set in analysis manager");
    return;
  }

  fMultSelection = new AliMultSelection("MultSelection");
  ah->AddBranch("AliMultSelection", &fMultSelection);
}
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

  // --- V channel equalization factors ------------------------------
  aodHeader->SetVZEROEqFactors(esd->GetVZEROEqFactors());
  
  
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
  if (fCopyTracklets) {
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

  // --- Copy V0 data ------------------------------------------------
  if (fCopyV0) {
    AliAODVZERO* vzeroData = aod->GetVZEROData();
    *vzeroData = *(esd->GetVZEROData());
  }
  // --- Copy V0 data ------------------------------------------------
  if (fCopyAD) {
    AliAODAD* adData = aod->GetADData();
    *adData = *(esd->GetADData());
  }
  // --- Copy ZDC data ------------------------------------------------
  if (fCopyZDC) {
    AliESDZDC* esdZDC = esd->GetZDCData();    
    AliAODZDC* zdcAOD = aod->GetZDCData();
    zdcAOD->SetZEM1Energy(esdZDC->GetZEM1Energy());
    zdcAOD->SetZEM2Energy(esdZDC->GetZEM2Energy());
    zdcAOD->SetZNCTowers(esdZDC->GetZNCTowerEnergy(),
			 esdZDC->GetZNCTowerEnergyLR());
    zdcAOD->SetZNATowers(esdZDC->GetZNATowerEnergy(),
			 esdZDC->GetZNATowerEnergyLR());
    zdcAOD->SetZPCTowers(esdZDC->GetZPCTowerEnergy(),
			 esdZDC->GetZPCTowerEnergyLR());
    zdcAOD->SetZPATowers(esdZDC->GetZPATowerEnergy(),
			 esdZDC->GetZPATowerEnergyLR());
    zdcAOD->SetZDCParticipants(esdZDC->GetZDCParticipants(),
			       esdZDC->GetZDCPartSideA(),
			       esdZDC->GetZDCPartSideC());
    zdcAOD->SetZDCImpactParameter(esdZDC->GetImpactParameter(),
				  esdZDC->GetImpactParamSideA(),
				  esdZDC->GetImpactParamSideC());
    zdcAOD->SetZDCTDCSum(esdZDC->GetZNTDCSum(0));	
    zdcAOD->SetZDCTDCDiff(esdZDC->GetZNTDCDiff(0));	
    if(esdZDC->IsZNChit()){
      Int_t cable = 10;
      if(esdZDC->IsZDCTDCcablingSet()){ // RUN2
	if(esdZDC->GetZNCTDCChannel()>0)
	  cable = esdZDC->GetZNCTDCChannel();
	else
	  cable = 16;
      }
      zdcAOD->SetZNCTDC(esdZDC->GetZDCTDCCorrected(cable, 0)); //RUN1
    }
    if(esdZDC->IsZNAhit()){
      Int_t cable = 12; 
      if(esdZDC->IsZDCTDCcablingSet()){ // RUN2
	if(esdZDC->GetZNATDCChannel()>0)
	  cable = esdZDC->GetZNATDCChannel();
	else
	  cable = 18;
      }
      zdcAOD->SetZNATDC(esdZDC->GetZDCTDCCorrected(cable, 0));
    }
    if(esdZDC->IsZPChit()){
      Int_t cable = 11;
      if(esdZDC->IsZDCTDCcablingSet()){ // RUN2
	if(esdZDC->GetZPCTDCChannel()>0)
	  cable = esdZDC->GetZPCTDCChannel();
	else
	  cable = 17;
      }
      zdcAOD->SetZPCTDC(esdZDC->GetZDCTDCCorrected(cable, 0));
    }
    if(esdZDC->IsZPAhit()){
      Int_t cable = 13;
      if(esdZDC->IsZDCTDCcablingSet()){ // RUN2
	if(esdZDC->GetZPATDCChannel()>0)
	  cable = esdZDC->GetZPATDCChannel();
	else
	  cable = 19; 
      }
      zdcAOD->SetZPATDC(esdZDC->GetZDCTDCCorrected(cable, 0));
    }
  }
  // --- Copy centrality estimates -----------------------------------
  if (fMultSelection) {
    TObject*          outO  = InputEvent()->FindListObject("MultSelection");
    if (outO) {
      AliMultSelection* outS = static_cast<AliMultSelection*>(outO);
      fMultSelection->Set(outS);
    }
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
