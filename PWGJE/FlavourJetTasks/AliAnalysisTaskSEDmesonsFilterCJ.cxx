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
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
//
// Task for selecting D mesons to be used as an input for D-Jet correlations
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL) XMZhang@lbl.gov
// S. Aiola (Yale University) salvatore.aiola@cern.ch
// S. Antônio (University of São Paulo / Utrecht University) antonio.silva@cern.ch
// B. Trzeciak (Utrecht University) barbara.antonia.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TH3F.h>

#include "AliLog.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliStack.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"

#include "AliAnalysisTaskSEDmesonsFilterCJ.h"

ClassImp(AliAnalysisTaskSEDmesonsFilterCJ)

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ() :
  AliAnalysisTaskEmcal(),
  fUseMCInfo(kFALSE),
  fBuildRMEff(kFALSE),
  fUsePythia(kFALSE),
  fUseReco(kTRUE),
  fCandidateType(0),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fCuts(0),
  fMinMass(0.),
  fMaxMass(0.),
  fInhibitTask(kFALSE),
  fCombineDmesons(kFALSE),
  fMultCand(kFALSE),
  fAnalyseCand(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fMCarray(0),
  fMCHeader(0),
  fCandidateArray(0),
  fSideBandArray(0),
  fCombinedDmesons(0),
  fCombinedDmesonsBkg(0),
  fMCCombinedDmesons(0),
  fNCand(0),
  fNSBCand(0),
  fHistStat(0),
  fHistNSBCandEv(0),
  fHistNCandEv(0),
  fHistImpParS(0),
  fHistImpParB(0),
  fHistPtPion(0),
  fHistInvMassPtD(0),
  fHistInvMassS(0),
  fHistInvMassB(0),
  fHistAlphaDDS(0),
  fHistAlphaDpisS(0),
  fHistAlphaDpiS(0),
  fHistAlphaDKS(0),
  fHistAlphaDDB(0),
  fHistAlphaDpisB(0),
  fHistAlphaDpiB(0),
  fHistAlphaDKB(0),
  fHistDeltaRDDS(0),
  fHistDeltaRDpisS(0),
  fHistDeltaRDpiS(0),
  fHistDeltaRDKS(0),
  fHistDeltaRDDB(0),
  fHistDeltaRDpisB(0),
  fHistDeltaRDpiB(0),
  fHistDeltaRDKB(0),
  fHistAlphaDpiR(0),
  fHistAlphaDKR(0),
  fHistDeltaRDpiR(0),
  fHistDeltaRDKR(0)
{
  //
  // Default constructor
  //

  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
  for (Int_t i=30; i--;) fSigmaD0[i] = 0.;

  fNeedEmcalGeom = kFALSE;
}

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ(const char *name, AliRDHFCuts *cuts, ECandidateType candtype) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fUseMCInfo(kFALSE),
  fBuildRMEff(kFALSE),
  fUsePythia(kFALSE),
  fUseReco(kTRUE),
  fCandidateType(candtype),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fCuts(cuts),
  fMinMass(0.),
  fMaxMass(0.),
  fInhibitTask(kFALSE),
  fCombineDmesons(kFALSE),
  fMultCand(kFALSE),
  fAnalyseCand(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fMCarray(0),
  fMCHeader(0),
  fCandidateArray(0),
  fSideBandArray(0),
  fCombinedDmesons(0),
  fCombinedDmesonsBkg(0),
  fMCCombinedDmesons(0),
  fNCand(0),
  fNSBCand(0),
  fHistStat(0),
  fHistNSBCandEv(0),
  fHistNCandEv(0),
  fHistImpParS(0),
  fHistImpParB(0),
  fHistPtPion(0),
  fHistInvMassPtD(0),
  fHistInvMassS(0),
  fHistInvMassB(0),
  fHistAlphaDDS(0),
  fHistAlphaDpisS(0),
  fHistAlphaDpiS(0),
  fHistAlphaDKS(0),
  fHistAlphaDDB(0),
  fHistAlphaDpisB(0),
  fHistAlphaDpiB(0),
  fHistAlphaDKB(0),
  fHistDeltaRDDS(0),
  fHistDeltaRDpisS(0),
  fHistDeltaRDpiS(0),
  fHistDeltaRDKS(0),
  fHistDeltaRDDB(0),
  fHistDeltaRDpisB(0),
  fHistDeltaRDpiB(0),
  fHistDeltaRDKB(0),
  fHistAlphaDpiR(0),
  fHistAlphaDKR(0),
  fHistDeltaRDpiR(0),
  fHistDeltaRDKR(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //

  Info("AliAnalysisTaskSEDmesonsFilterCJ","Calling Constructor");

  fNeedEmcalGeom = kFALSE;

  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
  for (Int_t i=30; i--;) fSigmaD0[i] = 0.;

  const Int_t nptbins = fCuts->GetNPtBins();
  Float_t defaultSigmaD013[20]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.040};

  switch (fCandidateType) {
    case 0 :
      fCandidateName = "D0";
      fPDGmother = 421;
      fNProngs = 2;
      fPDGdaughters[0] = 211;  // pi
      fPDGdaughters[1] = 321;  // K
      fPDGdaughters[2] = 0;    // empty
      fPDGdaughters[3] = 0;    // empty
      fBranchName = "D0toKpi";
      break;
    case 1 :
      fCandidateName = "DStar";
      fPDGmother = 413;
      fNProngs = 3;
      fPDGdaughters[1] = 211; // pi soft
      fPDGdaughters[0] = 421; // D0
      fPDGdaughters[2] = 211; // pi fromD0
      fPDGdaughters[3] = 321; // K from D0
      fBranchName = "Dstar";

      if (nptbins<20) {
      	 for (Int_t ipt=0; ipt<nptbins;ipt++) fSigmaD0[ipt] = defaultSigmaD013[ipt];
      }
      else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
    default :
      Warning("AliAnalysisTaskSEDmesonsFilterCJ", "%d not accepted!!", fCandidateType);
      break;
  }

  if (fCandidateType == kD0toKpi) SetMassLimits(0.15, fPDGmother);
  if (fCandidateType == kDstartoKpipi) SetMassLimits(0.015, fPDGmother);

  DefineOutput(1, TList::Class());       // histos
  DefineOutput(2, AliRDHFCuts::Class()); // my cuts
  DefineOutput(3, TClonesArray::Class()); //array of candidates
  DefineOutput(4, TClonesArray::Class()); //array of SB candidates
  DefineOutput(5, TClonesArray::Class()); //array of candidates and event tracks
  DefineOutput(6, TClonesArray::Class()); //array of SB candidates and event tracks
  DefineOutput(7, TClonesArray::Class()); //array of MC D and event MC particles
}

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::~AliAnalysisTaskSEDmesonsFilterCJ()
{
  //
  // Destructor
  //

  Info("~AliAnalysisTaskSEDmesonsFilterCJ","Calling Destructor");

  if (fCuts) {
    delete fCuts;
    fCuts   = 0;
  }

  if (fCandidateArray) {
    delete fCandidateArray;
    fCandidateArray = 0;
  }

  if (fSideBandArray) {
    delete fSideBandArray;
    fSideBandArray = 0;
  }

  if (fCombinedDmesons) {
    delete fCombinedDmesons;
    fCombinedDmesons = 0;
    delete fMCCombinedDmesons;
    fMCCombinedDmesons = 0;
  }

  if (fCombinedDmesonsBkg) {
    delete fCombinedDmesonsBkg;
    fCombinedDmesonsBkg = 0;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Init()
{
  //
  // Initialization
  //

  Info("AnalysisTaskSEDmesonsForJetCorrelations::Init()", "Entering method");

  switch (fCandidateType) {
    case 0:
    {
      AliRDHFCutsD0toKpi* copyfCutsDzero = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      copyfCutsDzero->SetName("AnalysisCutsDzero");
      PostData(2, copyfCutsDzero);  // Post the data
    } break;
    case 1:
    {
      AliRDHFCutsDStartoKpipi* copyfCutsDstar = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      copyfCutsDstar->SetName("AnalysisCutsDStar");
      PostData(2, copyfCutsDstar); // Post the cuts
    } break;
    default:
      return;
  }

  return;
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects()
{
  //
  // Create output objects
  //

  Info("UserCreateOutputObjects","CreateOutputObjects of task %s", GetName());

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  DefineHistoForAnalysis(); // define histograms

  if (fUseReco) {
    if (fCandidateType == kD0toKpi){
      fCandidateArray = new TClonesArray("AliAODRecoDecayHF2Prong",10);
      fSideBandArray = new TClonesArray("AliAODRecoDecayHF2Prong",10);
    }
    else if (fCandidateType == kDstartoKpipi) {
      fCandidateArray = new TClonesArray("AliAODRecoCascadeHF",10);
      fSideBandArray = new TClonesArray("AliAODRecoCascadeHF",10);
    }
    else {
      AliWarning(Form("Candidate type %d not recognized!", fCandidateType));
      return;
    }
  }
  else {
    fCandidateArray = new TClonesArray("AliAODMCParticle",10);
    fSideBandArray = new TClonesArray("TObject",0); // not used
  }

  if (fCombineDmesons) {
    fCombinedDmesons = new TClonesArray("AliEmcalParticle",50);
    fCombinedDmesonsBkg = new TClonesArray("AliEmcalParticle",50);
    fMCCombinedDmesons = new TClonesArray("AliAODMCParticle",50);
  }
  else {
    fCombinedDmesons = new TClonesArray("TObject",0); // not used
    fCombinedDmesonsBkg = new TClonesArray("TObject",0); // not used
    fMCCombinedDmesons = new TClonesArray("TObject",0); // not used
  }

  fCandidateArray->SetOwner();
  fCandidateArray->SetName(GetOutputSlot(3)->GetContainer()->GetName());

  //this is used for the DStar side bands and MC!
  fSideBandArray->SetOwner();
  fSideBandArray->SetName(GetOutputSlot(4)->GetContainer()->GetName());

  fCombinedDmesons->SetOwner();
  fCombinedDmesons->SetName(GetOutputSlot(5)->GetContainer()->GetName());

  fCombinedDmesonsBkg->SetOwner();
  fCombinedDmesonsBkg->SetName(GetOutputSlot(6)->GetContainer()->GetName());

  fMCCombinedDmesons->SetOwner();
  fMCCombinedDmesons->SetName(GetOutputSlot(7)->GetContainer()->GetName());

  PostData(1, fOutput);
  PostData(3, fCandidateArray);
  PostData(4, fSideBandArray);
  PostData(5, fCombinedDmesons);
  PostData(6, fCombinedDmesonsBkg);
  PostData(7, fMCCombinedDmesons);

  Info("UserCreateOutputObjects","Data posted for task %s", GetName());
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ExecOnce()
{
  //
  // To be executed only once, for the first event
  //

  AliDebug(2, "Entering ExecOnce()");

  if (fInhibitTask) return;

  // Load the event
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  if (fAodEvent) {
    fArrayDStartoD0pi = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject(fBranchName.Data()));
  }
  else {
    if (AODEvent() && IsStandardAOD()) {

      // In case there is an AOD handler writing a standard AOD, use the AOD
      // event in memory rather than the input (ESD) event.
      fAodEvent = dynamic_cast<AliAODEvent*>(AODEvent());

      // in this case the branches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
      if(aodHandler->GetExtensions()) {
        AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
        AliAODEvent *aodFromExt = ext->GetAOD();
        fArrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      }
      else {
        AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
        fInhibitTask = kTRUE;
        return;
      }
    }
  }

  if (fArrayDStartoD0pi) {
    TString objname(fArrayDStartoD0pi->GetClass()->GetName());
    TClass cls(objname);
    if (!cls.InheritsFrom("AliAODRecoDecayHF2Prong")) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from AliAODRecoDecayHF2Prong! Task will be disabled!",
                    GetName(), cls.GetName(), fArrayDStartoD0pi->GetName()));
      fInhibitTask = kTRUE;
      fArrayDStartoD0pi = 0;
      return;
    }
  }
  else {
    AliError(Form("Could not find array %s, skipping the event. Task '%s' will be disabled!", fBranchName.Data(), GetName()));
    fInhibitTask = kTRUE;
    return;
  }

  if (fUseMCInfo) {
    fMCarray = dynamic_cast<TClonesArray*>(fAodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMCarray) {
      AliError(Form("MC particles not found! Task '%s' will be disabled!", GetName()));
      fInhibitTask = kTRUE;
      return;
    }
    
    fMCHeader = (AliAODMCHeader*)fAodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!fMCHeader) {
      AliError(Form("MC header not found! Task '%s' will be disabled!", GetName()));
      return;
    }
  }

  if (fCombineDmesons) {
    AddObjectToEvent(fCombinedDmesons);
    AddObjectToEvent(fCombinedDmesonsBkg);
    if(fUseMCInfo) AddObjectToEvent(fMCCombinedDmesons);
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::Run()
{
  //
  // Analysis execution
  //

  if (fInhibitTask) return kFALSE;

  Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
  if (matchingAODdeltaAODlevel<=0) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return kFALSE;
  }

  AliDebug(2, "Entering Run()");

  //clear the TClonesArray from the previous event
  fCandidateArray->Clear();
  fSideBandArray->Clear();
  fCombinedDmesons->Clear();
  fCombinedDmesonsBkg->Clear();
  AliDebug(2, "TClonesArray cleared");

  fHistStat->Fill(0);

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return kFALSE;

  //Event selection
  Bool_t iseventselected = fCuts->IsEventSelected(fAodEvent);
  if (!iseventselected) return kFALSE;
  fHistStat->Fill(1);

  AliDebug(2, "Event selected");

  const Int_t nD = fArrayDStartoD0pi->GetEntriesFast();
  AliDebug(2, Form("Found %d vertices", nD));
  if (!fUseMCInfo) fHistStat->Fill(2, nD);
 

  Int_t pdgMeson = 413;
  if (fCandidateType == kD0toKpi) pdgMeson = 421;

  fNCand = 0;
  fNSBCand = 0;

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

      
// for MC response matrix of efficiency studies, fMultCand option only
if(fUseMCInfo && fBuildRMEff){
  
  const Int_t nDMC = fMCarray->GetEntriesFast(); 
  
  for (Int_t iMCcharm = 0; iMCcharm < nDMC; iMCcharm++){ //loop over MC D
    
    AliAODMCParticle* charmPart = static_cast<AliAODMCParticle*>(fMCarray->At(iMCcharm));
    if(TMath::Abs(charmPart->GetPdgCode()) != pdgMeson) continue;
    if(!charmPart) continue;
     fHistStat->Fill(2);
     
    Int_t origin = CheckOrigin(charmPart, fMCarray);
    if (origin < 0) continue;
    if (fRejectQuarkNotFound && origin == kQuarkNotFound) {
      fHistStat->Fill(6);
      continue;
    }
    if (fRejectDfromB && origin == kFromBottom) {
      fHistStat->Fill(7);
      continue;
    }
    if (fKeepOnlyDfromB && origin != kFromBottom) {
      fHistStat->Fill(8);
      continue;
    }

     Int_t decay = CheckDecayChannel(charmPart, fMCarray);
    if( TMath::Abs(charmPart->GetPdgCode()) == 413 && decay != kDecayDStartoKpipi) continue;
    if( TMath::Abs(charmPart->GetPdgCode()) == 421 && decay != kDecayD0toKpi) continue;

    fHistStat->Fill(3);
    
    if (fNCand==fAnalyseCand){  
        
        new ((*fMCCombinedDmesons)[0]) AliAODMCParticle(*charmPart); 
    
        AliAODRecoDecayHF2Prong* charmCand = 0;
        AliAODRecoCascadeHF* dstar = 0;
        
        // loop over reco D candidates to find a match to MC
        Int_t isRecoD = kFALSE;
        for (Int_t icharm = 0; icharm < nD; icharm++) {  
 
          charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fArrayDStartoD0pi->At(icharm)); // D candidates
          if (!charmCand) continue;
          if(!(vHF->FillRecoCand(fAodEvent,charmCand))) continue;
          
          Int_t nprongs = charmCand->GetNProngs();
          AliDebug(2, Form("Candidate is %d, and nprongs = %d", fCandidateType, nprongs));

          if (fCandidateType == kDstartoKpipi) {
              dstar = dynamic_cast<AliAODRecoCascadeHF*>(charmCand);
              if (!dstar) {
                Error("AliAnalysisTaskSEDmesonsFilterCJ::UserExec","Candidate type is D* but object type is wrong (should be AliAODRecoCascadeHF)");
                continue;
              }
          }
    
          Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
          Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi
          
          Int_t mcLabel = NULL;
          if (fCandidateType == kDstartoKpipi) mcLabel = dstar->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, fMCarray);
          else mcLabel = charmCand->MatchToMC(421, fMCarray, fNProngs, fPDGdaughters);

          if(mcLabel == iMCcharm) { isRecoD = kTRUE; break; }
        }
        
          if (!isRecoD) break;
         if (!charmCand) break;
         if (fCandidateType == kDstartoKpipi &&  !dstar) break;
        
         fHistStat->Fill(4);

        // region of interest + cuts
        if (!fCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(pdgMeson))) break;
        
        //candidate selected by cuts and PID
        Int_t isSelected = 0;
        isSelected = fCuts->IsSelected(charmCand, AliRDHFCuts::kAll, fAodEvent); //selected
        if (!isSelected) break;

        fHistStat->Fill(5);

        if (fCandidateType == kDstartoKpipi) {
          AliAODRecoDecayHF2Prong* D0fromDstar = dstar->Get2Prong();
          fHistInvMassS->Fill(dstar->DeltaInvMass());
          fHistImpParS->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
          fHistImpParS->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
          fHistImpParS->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));
        }
        else {
          fHistImpParS->Fill(charmCand->Getd0Prong(0), charmCand->PtProng(0));
          fHistImpParS->Fill(charmCand->Getd0Prong(1), charmCand->PtProng(1));
          fHistInvMassS->Fill(charmCand->InvMassD0());
          fHistInvMassS->Fill(charmCand->InvMassD0bar());
        }
        
        if (fCandidateType == kDstartoKpipi) {
          Int_t isDstar = charmPart == NULL ? 0 : 1;
           //fill histograms of kinematics, using MC truth
          FillDStarMCTruthKinHistos(dstar, isSelected, isDstar);

          new ((*fCandidateArray)[0]) AliAODRecoCascadeHF(*dstar);
          new ((*fCombinedDmesons)[0]) AliEmcalParticle(dstar);
         
        }
        else {
         
          Int_t isD0 = 0;
          Int_t pdgCode = charmPart->GetPdgCode();
          if (pdgCode ==  421)  { isD0 = 1; }
          else if (pdgCode == -421)  { isD0 = -1; }
           //fill histograms of kinematics, using MC truth
          FillD0MCTruthKinHistos(charmCand, isSelected, isD0);
  
          new ((*fCandidateArray)[0]) AliAODRecoDecayHF2Prong(*charmCand);
          new ((*fCombinedDmesons)[0]) AliEmcalParticle(charmCand);
      }
        
        break;
   }
   else { fNCand++; }         
  
} 
}
// for data, or MC without RM or efficiency studies
else {
  for (Int_t icharm = 0; icharm < nD; icharm++) {   //loop over D candidates
    Int_t isSelected = 0;

    AliAODRecoDecayHF2Prong* charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fArrayDStartoD0pi->At(icharm)); // D candidates
    if (!charmCand) continue;
    if(!(vHF->FillRecoCand(fAodEvent,charmCand))) continue;

    Int_t nprongs = charmCand->GetNProngs();
    AliDebug(2, Form("Candidate is %d, and nprongs = %d", fCandidateType, nprongs));

    AliAODRecoCascadeHF* dstar = 0;

    if (fCandidateType == kDstartoKpipi) {
      dstar = dynamic_cast<AliAODRecoCascadeHF*>(charmCand);
      if (!dstar) {
        Error("AliAnalysisTaskSEDmesonsFilterCJ::UserExec","Candidate type is D* but object type is wrong (should be AliAODRecoCascadeHF)");
        continue;
      }
    }

    // region of interest + cuts
    if (!fCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(pdgMeson))) continue;

    if (!fUseMCInfo && fCandidateType == kDstartoKpipi) {
      FillDstarSideBands(dstar);
    }

    
    //candidate selected by cuts and PID
    isSelected = fCuts->IsSelected(charmCand, AliRDHFCuts::kAll, fAodEvent); //selected
    if (!isSelected) continue;

    if (fCandidateType == kDstartoKpipi) {
      ProcessDstar(dstar, isSelected);
    }
    else {
      ProcessD0(charmCand, isSelected);
    }
  } // end of D cand loop
}


  delete vHF;
 

  AliDebug(2, "Loop done");

  if (fCombineDmesons) {
    if (fCombinedDmesons->GetEntriesFast() > 0) {
      AddEventTracks(fCombinedDmesons, GetParticleContainer(0));
    }

    if (fMCCombinedDmesons->GetEntriesFast() > 0) {
        AddMCEventTracks(fMCCombinedDmesons, GetParticleContainer(1));
    }

    if (fCombinedDmesonsBkg->GetEntriesFast() > 0) {
      AddEventTracks(fCombinedDmesonsBkg, GetParticleContainer(0));
    }
  }

  fHistNCandEv->Fill(fCandidateArray->GetEntriesFast());
  if (fCandidateType == kDstartoKpipi || fUseMCInfo) {
    Int_t nsbcand = fSideBandArray->GetEntriesFast();
    fHistStat->Fill(4, nsbcand);
    fHistNSBCandEv->Fill(nsbcand);
  }

  PostData(1, fOutput);
  PostData(3, fCandidateArray);
  PostData(4, fSideBandArray);
  PostData(5, fCombinedDmesons);
  PostData(6, fCombinedDmesonsBkg);
  PostData(7, fMCCombinedDmesons);

  AliDebug(2, "Exiting method");

  return kTRUE;
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ProcessD0(AliAODRecoDecayHF2Prong* charmCand, Int_t isSelected)
{
  // Process the D0 candidate.

  // For MC analysis look for the AOD MC particle associated with the charm candidate
  AliAODMCParticle* charmPart = 0;
  if (fUseMCInfo) {
    Int_t mcLabel = charmCand->MatchToMC(421, fMCarray, fNProngs, fPDGdaughters);

    if (mcLabel >= 0) {
      charmPart = static_cast<AliAODMCParticle*>(fMCarray->At(mcLabel));
    }
  }

  if (charmPart) {

    Int_t origin = CheckOrigin(charmPart, fMCarray);
    if (origin < 0) return;

    if (fRejectQuarkNotFound && origin == kQuarkNotFound) {
      fHistStat->Fill(6);
      return;
    }
    if (fRejectDfromB && origin == kFromBottom) {
      fHistStat->Fill(7);
      return;
    }
    if (fKeepOnlyDfromB && origin != kFromBottom) {
      fHistStat->Fill(8);
      return;
    }

    fHistStat->Fill(2);
  }
  else {
    fHistStat->Fill(5);
  }

  // For MC background fill fSideBandArray
  if (fUseMCInfo && !charmPart) {
    if (fUseReco) {
      if(!fMultCand) new ((*fSideBandArray)[fNSBCand]) AliAODRecoDecayHF2Prong(*charmCand);
      else if(fNSBCand==fAnalyseCand) new ((*fSideBandArray)[0]) AliAODRecoDecayHF2Prong(*charmCand);

      if (fCombineDmesons) {
        if(!fMultCand) new ((*fCombinedDmesonsBkg)[fNSBCand]) AliEmcalParticle(charmCand);
        else if(fNSBCand==fAnalyseCand) new ((*fCombinedDmesonsBkg)[0]) AliEmcalParticle(charmCand);
      }

      fHistImpParB->Fill(charmCand->Getd0Prong(0), charmCand->PtProng(0));
      fHistImpParB->Fill(charmCand->Getd0Prong(1), charmCand->PtProng(1));
      fHistInvMassB->Fill(charmCand->InvMassD0());
      fHistInvMassB->Fill(charmCand->InvMassD0bar());

      fNSBCand++;
    }
  }
  // For data or MC signal fill fCandidateArray
  else {
    // For data or MC with the requirement fUseReco fill with candidates
    if (fUseReco) {
      if(!fMultCand) new ((*fCandidateArray)[fNCand]) AliAODRecoDecayHF2Prong(*charmCand);
      else if(fNCand==fAnalyseCand) new ((*fCandidateArray)[0]) AliAODRecoDecayHF2Prong(*charmCand);

      if (fCombineDmesons) {
        if(!fMultCand) new ((*fCombinedDmesons)[fNCand]) AliEmcalParticle(charmCand);
        else if(fNCand==fAnalyseCand)
        {
            new ((*fCombinedDmesons)[0]) AliEmcalParticle(charmCand);
            if(charmPart) new ((*fMCCombinedDmesons)[0]) AliAODMCParticle(*charmPart);
        }
      }

      fHistImpParS->Fill(charmCand->Getd0Prong(0), charmCand->PtProng(0));
      fHistImpParS->Fill(charmCand->Getd0Prong(1), charmCand->PtProng(1));
      fHistInvMassS->Fill(charmCand->InvMassD0());
      fHistInvMassS->Fill(charmCand->InvMassD0bar());
    }
    // For MC with requirement particle level fill with AliAODMCParticle
    else {
      if(!fMultCand) new ((*fCandidateArray)[fNCand]) AliAODMCParticle(*charmPart);
      else if(fNCand==fAnalyseCand) new ((*fCandidateArray)[0]) AliAODMCParticle(*charmPart);
    }
    fHistStat->Fill(3);
    fNCand++;
  }

  // Now filling some more histograms

  // mass vs pt
  if (isSelected == 1 || isSelected == 3) fHistInvMassPtD->Fill(charmCand->InvMassD0(), charmCand->Pt());
  if (isSelected >= 2) fHistInvMassPtD->Fill(charmCand->InvMassD0bar(), charmCand->Pt());

  if (fUseMCInfo) {  //fill histograms of kinematics, using MC truth
    Int_t isD0 = 0;
    if (charmPart) {
      Int_t pdgCode = charmPart->GetPdgCode();

      if (pdgCode ==  421)  {
        isD0 = 1;
      }
      else if (pdgCode == -421)  {
        isD0 = -1;
      }
      else {
        AliDebug(2, "Not a D0/D0bar meson!");
        return;
      }
    }

    FillD0MCTruthKinHistos(charmCand, isSelected, isD0);
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ProcessDstar(AliAODRecoCascadeHF* dstar, Int_t isSelected)
{
  // Process the D* candidate.

  AliDebug(2, "Entering method");

  // For MC analysis look for the AOD MC particle associated with the charm candidate
  AliAODMCParticle* charmPart = 0;
  if (fUseMCInfo) {
    //D* and D0 prongs needed to MatchToMC method
    Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
    Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi

    Int_t mcLabel = dstar->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, fMCarray);

    if (mcLabel >= 0) {
      charmPart = static_cast<AliAODMCParticle*>(fMCarray->At(mcLabel));
    }
  }

  if (charmPart) {

    Int_t origin = CheckOrigin(charmPart, fMCarray);
    if (origin < 0) return;

    if (fRejectQuarkNotFound && origin == kQuarkNotFound) {
      fHistStat->Fill(6);
      return;
    }
    if (fRejectDfromB && origin == kFromBottom) {
      fHistStat->Fill(7);
      return;
    }
    if (fKeepOnlyDfromB && origin != kFromBottom) {
      fHistStat->Fill(8);
      return;
    }

    fHistStat->Fill(2);
  }
  else {
    fHistStat->Fill(5);
  }

  AliAODRecoDecayHF2Prong* D0fromDstar = dstar->Get2Prong();

  // For MC background fill fSideBandArray
  if (fUseMCInfo && !charmPart) {
    if (fUseReco) {
      if(!fMultCand) new ((*fSideBandArray)[fNSBCand]) AliAODRecoCascadeHF(*dstar);
      else if(fNSBCand==fAnalyseCand) new ((*fSideBandArray)[0]) AliAODRecoCascadeHF(*dstar);

      if (fCombineDmesons) {
        if(!fMultCand) new ((*fCombinedDmesonsBkg)[fNSBCand]) AliEmcalParticle(dstar);
        else if(fNSBCand==fAnalyseCand) new ((*fCombinedDmesonsBkg)[0]) AliEmcalParticle(dstar);
      }

      fHistInvMassB->Fill(dstar->DeltaInvMass());
      fHistImpParB->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
      fHistImpParB->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
      fHistImpParB->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));

      fNSBCand++;
    }
  }
  // For data and MC signal fill fCandidateArray
  else {
    // For data or MC signal with the requirement fUseReco fill with candidates
    if (fUseReco) {
      if(!fMultCand) new ((*fCandidateArray)[fNCand]) AliAODRecoCascadeHF(*dstar);
      else if(fNCand==fAnalyseCand) new ((*fCandidateArray)[0]) AliAODRecoCascadeHF(*dstar);

      if (fCombineDmesons) {
        if(!fMultCand) new ((*fCombinedDmesons)[fNCand]) AliEmcalParticle(dstar);
        else if(fNCand==fAnalyseCand)
        {
            new ((*fCombinedDmesons)[0]) AliEmcalParticle(dstar);
            if(charmPart) new ((*fMCCombinedDmesons)[0]) AliAODMCParticle(*charmPart);
        }
      }
    }
    // For MC signal with requirement particle level fill with AliAODMCParticle
    else {
      if(!fMultCand) new ((*fCandidateArray)[fNCand]) AliAODMCParticle(*charmPart);
      else if(fNCand==fAnalyseCand) new ((*fCandidateArray)[0]) AliAODMCParticle(*charmPart);
    }
    fNCand++;

    fHistInvMassS->Fill(dstar->DeltaInvMass());
    fHistImpParS->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
    fHistImpParS->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
    fHistImpParS->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));
    fHistStat->Fill(3);
  }

  // Now filling some more histograms.

  // select D* in the D0 window.
  // In the cut object window is loose to allow for side bands

  // retrieve the corresponding pt bin for the candidate
  // and set the expected D0 width (x3)

  Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t invMassPDG = mPDGDstar - mPDGD0;

  Double_t ptD = dstar->Pt();
  Int_t ptDbin = fCuts->PtBin(ptD);
  if (ptDbin < 0 || ptDbin >= fCuts->GetNPtBins()) {
    AliError(Form("Pt %.3f out of bounds", ptD));
    return;
  }

  AliDebug(1, Form("Pt bin %d and sigma D0 %.4f", ptDbin, fSigmaD0[ptDbin]));
  //consider the D* candidates only if the mass of the D0 is in 3 sigma wrt the PDG value
  if ((dstar->InvMassD0()>=(mPDGD0-3.*fSigmaD0[ptDbin])) && (dstar->InvMassD0()<=(mPDGD0+3.*fSigmaD0[ptDbin]))) {
    AliDebug(2, "D0 mass within 3 sigma");

    // D* delta mass
    AliDebug(2, "Filling invariant mass vs pt histogram");
    fHistInvMassPtD->Fill(dstar->DeltaInvMass(), ptD); // 2 D slice for pt bins

    // Soft pion pt for good candidates
    Double_t invmassDelta = dstar->DeltaInvMass();
    if (TMath::Abs(invmassDelta - invMassPDG) < 0.0021) {
      AliDebug(2, "Filling pion pt histogram");
      //softpion from D* decay
      AliAODTrack *softPionTrack = static_cast<AliAODTrack*>(dstar->GetBachelor());
      if (softPionTrack) fHistPtPion->Fill(softPionTrack->Pt());
    }
  }

  if (fUseMCInfo) { //fill histograms of kinematics, using MC truth
    Int_t isDstar = charmPart == NULL ? 0 : 1;
    FillDStarMCTruthKinHistos(dstar, isDstar, isSelected);
  }

  AliDebug(2, "Exiting method");
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::FillD0MCTruthKinHistos(AliAODRecoDecayHF2Prong* charmCand, Int_t isSelected, Int_t isD0)
{
  // Fill some histogram with kinematic information of the D0 candidate.

  Double_t ptD = charmCand->Pt();
  Double_t aD = charmCand->Phi();
  Double_t adaugh[2] = {charmCand->PhiProng(0), charmCand->PhiProng(1)};
  AliAODTrack* p0 = static_cast<AliAODTrack*>(charmCand->GetDaughter(0));
  AliAODTrack* p1 = static_cast<AliAODTrack*>(charmCand->GetDaughter(1));
  Float_t dR0 = DeltaR(charmCand, p0);
  Float_t dR1 = DeltaR(charmCand, p1);

  if (isD0 == 0) {  //background
    if (isSelected == 1 || isSelected == 3) { // selected as D0
      fHistAlphaDKB->Fill(aD-adaugh[0],ptD);
      fHistAlphaDpiB->Fill(aD-adaugh[1],ptD);

      fHistDeltaRDKB->Fill(dR0,ptD);
      fHistDeltaRDpiB->Fill(dR1,ptD);
    }
    if (isSelected >= 2) { //selected as D0bar
      fHistAlphaDpiB->Fill(aD-adaugh[0],ptD);
      fHistAlphaDKB->Fill(aD-adaugh[1],ptD);

      fHistDeltaRDpiB->Fill(dR0,ptD);
      fHistDeltaRDKB->Fill(dR1,ptD);
    }
  }
  else if (isD0 == 1) { //D0
    fHistAlphaDKS->Fill(aD-adaugh[0],ptD);
    fHistAlphaDpiS->Fill(aD-adaugh[1],ptD);

    fHistDeltaRDKS->Fill(dR0,ptD);
    fHistDeltaRDpiS->Fill(dR1,ptD);

    if (isSelected == 3) { // selected as both D0bar/D0
      fHistAlphaDpiR->Fill(aD-adaugh[0],ptD);
      fHistAlphaDKR->Fill(aD-adaugh[1],ptD);

      fHistDeltaRDpiR->Fill(dR0,ptD);
      fHistDeltaRDKR->Fill(dR1,ptD);
    }
  }
  else if (isD0 == -1) { //D0bar
    fHistAlphaDKS->Fill(aD-adaugh[1],ptD);
    fHistAlphaDpiS->Fill(aD-adaugh[0],ptD);

    fHistDeltaRDKS->Fill(dR1,ptD);
    fHistDeltaRDpiS->Fill(dR0,ptD);

    if (isSelected == 3) { // selected as both D0bar/D0
      fHistAlphaDpiR->Fill(aD-adaugh[1],ptD);
      fHistAlphaDKR->Fill(aD-adaugh[0],ptD);

      fHistDeltaRDpiR->Fill(dR1, ptD);
      fHistDeltaRDKR->Fill(dR0, ptD);
    }
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::FillDStarMCTruthKinHistos(AliAODRecoCascadeHF* dstar, Int_t isSelected, Int_t isDstar)
{
  // Fill some histogram with kinematic information of the D0 candidate.

  AliAODTrack *softPionTrack = static_cast<AliAODTrack*>(dstar->GetBachelor());
  Double_t ptD  = dstar->Pt();
  Double_t aD  = dstar->Phi();
  Double_t apis= softPionTrack->Phi();

  AliAODRecoDecayHF2Prong* D0fromDstar = dstar->Get2Prong();
  Double_t aD0    = D0fromDstar->Phi();
  Int_t    isD0   = D0fromDstar->Charge()>0 ? kTRUE : kFALSE;
  Double_t aK     = isD0 ? D0fromDstar->PhiProng(0) : D0fromDstar->PhiProng(1);
  Double_t api    = isD0 ? D0fromDstar->PhiProng(1) : D0fromDstar->PhiProng(0);
  Double_t dRDD0  = DeltaR(dstar,D0fromDstar);
  Double_t dRDpis = DeltaR(dstar,softPionTrack);
  Double_t dRDpi  = DeltaR(dstar, isD0 ? static_cast<AliVParticle*>(D0fromDstar->GetDaughter(1)) : static_cast<AliVParticle*>(D0fromDstar->GetDaughter(0)));
  Double_t dRDK   = DeltaR(dstar, isD0 ? static_cast<AliVParticle*>(D0fromDstar->GetDaughter(0)) : static_cast<AliVParticle*>(D0fromDstar->GetDaughter(1)));

  if (isDstar) {
    fHistAlphaDDS  ->Fill(aD-aD0, ptD);
    fHistAlphaDpisS->Fill(aD-apis, ptD);
    fHistAlphaDpiS ->Fill(aD-api, ptD);
    fHistAlphaDKS  ->Fill(aD-aK, ptD);

    fHistDeltaRDDS  ->Fill(dRDD0, ptD);
    fHistDeltaRDpisS->Fill(dRDpis, ptD);
    fHistDeltaRDpiS ->Fill(dRDpi, ptD);
    fHistDeltaRDKS  ->Fill(dRDK, ptD);
  }
  else {
    fHistAlphaDDB  ->Fill(aD-aD0, ptD);
    fHistAlphaDpisB->Fill(aD-apis, ptD);
    fHistAlphaDpiB ->Fill(aD-api, ptD);
    fHistAlphaDKB  ->Fill(aD-aK, ptD);

    fHistDeltaRDDB  ->Fill(dRDD0, ptD);
    fHistDeltaRDpisB->Fill(dRDpis, ptD);
    fHistDeltaRDpiB ->Fill(dRDpi, ptD);
    fHistDeltaRDKB  ->Fill(dRDK, ptD);
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::FillDstarSideBands(AliAODRecoCascadeHF* dstar)
{
  // Fills the array of Dstar side band candidates.

  //select by track cuts the side band candidates (don't want mass cut)
  Int_t isSelected = fCuts->IsSelected(dstar, AliRDHFCuts::kTracks, fAodEvent);
  if (!isSelected) return;

  //add a reasonable cut on the invariant mass (e.g. (+-2\sigma, +-10 \sigma), with \sigma = fSigmaD0[bin])
  Int_t bin = fCuts->PtBin(dstar->Pt());
  if (bin < 0 || bin >= fCuts->GetNPtBins()) {
    AliError(Form("Pt %.3f out of bounds", dstar->Pt()));
    return;
  }

  const Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

  //if data and Dstar from D0 side band
  if (((dstar->InvMassD0() < (mPDGD0-3.*fSigmaD0[bin])) && (dstar->InvMassD0() > (mPDGD0-10.*fSigmaD0[bin]))) /*left side band*/   ||
      ((dstar->InvMassD0() > (mPDGD0+3.*fSigmaD0[bin])) && (dstar->InvMassD0() < (mPDGD0+10.*fSigmaD0[bin]))) /*right side band*/) {

    if(!fMultCand) new ((*fSideBandArray)[fNSBCand]) AliAODRecoCascadeHF(*dstar);
    else if(fNSBCand==fAnalyseCand) new ((*fSideBandArray)[0]) AliAODRecoCascadeHF(*dstar);
    fNSBCand++;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t range, Int_t pdg)
{
  // Set the mass limits using the PDG code and the given range.

  Float_t mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();

  // compute the Delta mass for the D*
  if (fCandidateType==kDstartoKpipi) mass -= TDatabasePDG::Instance()->GetParticle(421)->Mass();

  fMinMass = mass - range;
  fMaxMass = mass + range;

  AliInfo(Form("Setting mass limits to %f, %f", fMinMass, fMaxMass));
  if ((fMinMass<0.) || (fMaxMass<=0.) || (fMaxMass<fMinMass)) {
    AliError(Form("Wrong mass limits! Task '%s' will be disabled!", GetName()));
    fInhibitTask = kTRUE;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t lowlimit, Double_t uplimit)
{
  // Set the mass limits.

  if (uplimit>lowlimit) {
    fMinMass = lowlimit;
    fMaxMass = uplimit;
  }
  else {
    printf("Error! Lower limit larger than upper limit!\n");
    fMinMass = uplimit - uplimit*0.2;
    fMaxMass = uplimit;
  }
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar(Int_t nptbins, Float_t* width)
{
  // Set the D0 width for the D* analysis.

  if (nptbins > 30) {
    AliWarning("Maximum number of bins allowed is 30!");
    return kFALSE;
  }

  if (!width) return kFALSE;
  for (Int_t ipt=0; ipt < nptbins; ipt++) fSigmaD0[ipt] = width[ipt];

  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis()
{
  // Allocate the histograms for the analysis.

  // Statistics
  fHistStat = new TH1I("fHistStat", "Statistics", 9, -0.5, 8.5);
  fHistStat->GetXaxis()->SetBinLabel(1, "N ev anal");
  fHistStat->GetXaxis()->SetBinLabel(2, "N ev sel");
  if(fUseMCInfo && fBuildRMEff){
    fHistStat->GetXaxis()->SetBinLabel(3, "N Gen D");
    fHistStat->GetXaxis()->SetBinLabel(4, "N Gen Sel D");
    fHistStat->GetXaxis()->SetBinLabel(5, "N D");
    fHistStat->GetXaxis()->SetBinLabel(6, "N cand sel cuts");
    fHistStat->GetXaxis()->SetBinLabel(7, "N rej no quark");
    fHistStat->GetXaxis()->SetBinLabel(8, "N rej from B");
    fHistStat->GetXaxis()->SetBinLabel(9, "N rej from D");
  }
  else {
  if (fUseMCInfo) {
    if (fUseReco) {
      fHistStat->GetXaxis()->SetBinLabel(3, "N D");
    }
    else {
      fHistStat->GetXaxis()->SetBinLabel(3, "N Gen D");
    }
  }
  else {
    fHistStat->GetXaxis()->SetBinLabel(3, "N Cand");
  }
  fHistStat->GetXaxis()->SetBinLabel(4, "N cand sel cuts");
  if (fCandidateType == kDstartoKpipi) {
    fHistStat->GetXaxis()->SetBinLabel(5, "N side band cand");
  }
  if (fUseMCInfo) {
    fHistStat->GetXaxis()->SetBinLabel(6, "N Background");
    fHistStat->GetXaxis()->SetBinLabel(7, "N rej no quark");
    fHistStat->GetXaxis()->SetBinLabel(8, "N rej from B");
    fHistStat->GetXaxis()->SetBinLabel(9, "N rej from D");
  }
  }
  
  fHistStat->SetNdivisions(1);
  fOutput->Add(fHistStat);

  fHistNCandEv = new TH1F("fHistNCandEv", "Number of candidates per event (after cuts);# cand/ev", 100, 0., 100.);
  fOutput->Add(fHistNCandEv);

  if ((fCandidateType == kDstartoKpipi) || fUseMCInfo) {
    fHistNSBCandEv = new TH1F("hnSBCandEv", "Number of side bands candidates per event (after cuts);# cand/ev", 100, 0.,100.);
    fOutput->Add(fHistNSBCandEv);
  }

  if (fCandidateType == kDstartoKpipi) {
    fHistPtPion = new TH1F("fHistPtPion", "Primary pions candidates pt", 500, 0., 10.);
    fHistPtPion->SetStats(kTRUE);
    fHistPtPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPtPion->GetYaxis()->SetTitle("entries");
    fOutput->Add(fHistPtPion);
  }

  // Invariant mass related histograms
  const Int_t nbinsmass = 200;
  const Int_t ptbinsD = 100;
  Float_t ptmin =0.;
  Float_t ptmax =50.;
  fHistInvMassPtD = new TH2F("fHistInvMassPtD", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass, ptbinsD, ptmin, ptmax);
  fHistInvMassPtD->SetStats(kTRUE);
  fHistInvMassPtD->GetXaxis()->SetTitle("mass (GeV/c)");
  fHistInvMassPtD->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fOutput->Add(fHistInvMassPtD);

  fHistImpParS = new TH2F("fHistImpParS", "Impact parameter of daughter tracks; Getd0Prong();#it{p}^{daugh}_{T} (GeV/c)",200, -0.1,0.1,ptbinsD, ptmin, ptmax); //same range of pt of the D, but pt daughter used
  fOutput->Add(fHistImpParS);

  fHistInvMassS = new TH1F("fHistInvMassS", "D invariant mass distribution (filled with fCandidateArray)", nbinsmass, fMinMass, fMaxMass);
  fHistInvMassS->SetStats(kTRUE);
  fHistInvMassS->GetXaxis()->SetTitle("mass (GeV/c)");
  fHistInvMassS->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fOutput->Add(fHistInvMassS);

  const Int_t nbinsalpha = 200;
  Float_t minalpha = -TMath::Pi();
  Float_t maxalpha = TMath::Pi();
  const Int_t nbinsdeltaR = 200;
  Float_t mindeltaR = 0.;
  Float_t maxdeltaR = 10.;



  if (fUseMCInfo) {
    fHistImpParB = new TH2F("fHistImpParB", "Impact parameter of daughter tracks (Background); Getd0Prong();#it{p}^{daugh}_{T} (GeV/c)",200, -0.1,0.1,ptbinsD, ptmin, ptmax); //same range of pt of the D, but pt daughter used
    fOutput->Add(fHistImpParB);

    fHistInvMassB = new TH1F("fHistInvMassB", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass);
    fHistInvMassB->SetStats(kTRUE);
    fHistInvMassB->GetXaxis()->SetTitle("mass (GeV/c)");
    fHistInvMassB->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fOutput->Add(fHistInvMassB);

    if (fCandidateType == kDstartoKpipi) {

      fHistAlphaDDS    = new TH2F("fHistAlphaDDS","Angle D^{*}-D^{0} (Signal);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpisS  = new TH2F("fHistAlphaDpisS","Angle D^{*}-#pi_{soft} (Signal);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiS   = new TH2F("fHistAlphaDpiS","Angle D^{*}-#pi (Signal);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKS    = new TH2F("fHistAlphaDKS","Angle D^{*}-K (Signal);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistAlphaDDB    = new TH2F("fHistAlphaDDB","Angle D^{*}-D^{0} (Background);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpisB  = new TH2F("fHistAlphaDpisB","Angle D^{*}-#pi_{soft} (Background);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiB   = new TH2F("fHistAlphaDpiB","Angle D^{*}-#pi (Background);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKB    = new TH2F("fHistAlphaDKB","Angle D^{*}-K (Background);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistDeltaRDDS   = new TH2F("fHistDeltaRDDS","Angle D^{*}-D^{0} (Signal);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpisS = new TH2F("fHistDeltaRDpisS","Angle D^{*}-#pi_{soft} (Signal);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiS  = new TH2F("fHistDeltaRDpiS","Angle D^{*}-#pi (Signal);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKS   = new TH2F("fHistDeltaRDKS","Angle D^{*}-K (Signal);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fHistDeltaRDDB   = new TH2F("fHistDeltaRDDB","Angle D^{*}-D^{0} (Background);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpisB = new TH2F("fHistDeltaRDpisB","Angle D^{*}-#pi_{soft} (Background);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiB  = new TH2F("fHistDeltaRDpiB","Angle D^{*}-#pi (Background);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKB   = new TH2F("fHistDeltaRDKB","Angle D^{*}-K (Background);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fOutput->Add(fHistAlphaDDS);
      fOutput->Add(fHistAlphaDpisS);
      fOutput->Add(fHistAlphaDpiS);
      fOutput->Add(fHistAlphaDKS);
      fOutput->Add(fHistAlphaDDB);
      fOutput->Add(fHistAlphaDpisB);
      fOutput->Add(fHistAlphaDpiB);
      fOutput->Add(fHistAlphaDKB);

      fOutput->Add(fHistDeltaRDDS);
      fOutput->Add(fHistDeltaRDpisS);
      fOutput->Add(fHistDeltaRDpiS);
      fOutput->Add(fHistDeltaRDKS);
      fOutput->Add(fHistDeltaRDDB);
      fOutput->Add(fHistDeltaRDpisB);
      fOutput->Add(fHistDeltaRDpiB);
      fOutput->Add(fHistDeltaRDKB);
    }
    else if (fCandidateType == kD0toKpi) {

      fHistAlphaDpiS  = new TH2F("fHistAlphaDpiS","Angle D^{0}-#pi (Signal);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKS   = new TH2F("fHistAlphaDKS","Angle D^{0}-K (Signal);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiR  = new TH2F("fHistAlphaDpiR","Angle D^{0}-#pi (Reflections);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKR   = new TH2F("fHistAlphaDKR","Angle D^{0}-K (Reflections);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistAlphaDpiB  = new TH2F("fHistAlphaDpiB","Angle D^{0}-#pi (Background);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKB   = new TH2F("fHistAlphaDKB","Angle D^{0}-K (Background);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistDeltaRDpiS = new TH2F("fHistDeltaRDpiS","Angle D^{0}-#pi (Signal);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKS  = new TH2F("fHistDeltaRDKS","Angle D^{0}-K (Signal);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiR = new TH2F("fHistDeltaRDpiR","Angle D^{0}-#pi (Reflections);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKR  = new TH2F("fHistDeltaRDKR","Angle D^{0}-K (Reflections);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fHistDeltaRDpiB = new TH2F("fHistDeltaRDpiB","Angle D^{0}-#pi (Background);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKB  = new TH2F("fHistDeltaRDKB","Angle D^{0}-K (Background);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fOutput->Add(fHistAlphaDpiS);
      fOutput->Add(fHistAlphaDKS);
      fOutput->Add(fHistAlphaDpiR);
      fOutput->Add(fHistAlphaDKR);
      fOutput->Add(fHistAlphaDpiB);
      fOutput->Add(fHistAlphaDKB);

      fOutput->Add(fHistDeltaRDpiS);
      fOutput->Add(fHistDeltaRDKS);
      fOutput->Add(fHistDeltaRDpiR);
      fOutput->Add(fHistDeltaRDKR);
      fOutput->Add(fHistDeltaRDpiB);
      fOutput->Add(fHistDeltaRDKB);
    }
  }

  return kTRUE;
}

//_______________________________________________________________________________
Float_t AliAnalysisTaskSEDmesonsFilterCJ::DeltaR(AliVParticle *p1, AliVParticle *p2) const
{
  // Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)

  if (!p1 || !p2) return -1;

  Double_t phi1 = p1->Phi();
  Double_t eta1 = p1->Eta();
  Double_t phi2 = p2->Phi();
  Double_t eta2 = p2->Eta();

  Double_t dPhi = phi1 - phi2;
  if(dPhi < -(TMath::Pi())/2)    dPhi = dPhi + TMath::TwoPi();
  if(dPhi > (3*(TMath::Pi()))/2) dPhi = dPhi - TMath::TwoPi();

  Double_t dEta = eta1 - eta2;

  Double_t deltaR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

  return deltaR;
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::AddEventTracks(TClonesArray* coll, AliParticleContainer* tracks)
{
  //
  // Add event tracks to a collection that already contains the D candidates, excluding the daughters of the D candidates
  //

  if (!tracks) return;

  TObjArray allDaughters(10);
  allDaughters.SetOwner(kFALSE);

  TIter next(coll);
  AliEmcalParticle* emcpart = 0;
  while ((emcpart = static_cast<AliEmcalParticle*>(next()))) {
    AliAODRecoDecay* reco = dynamic_cast<AliAODRecoDecay*>(emcpart->GetTrack());
    AliDebug(2, Form("Found a D meson candidtate with pT = %.3f, eta = %.3f, phi = %.3f", reco->Pt(), reco->Eta(), reco->Phi()));
    if (reco) AddDaughters(reco, allDaughters);
  }

  tracks->ResetCurrentID();
  AliVTrack* track = 0;
  Int_t n = coll->GetEntriesFast();
 
  while ((track = static_cast<AliVTrack*>(tracks->GetNextAcceptParticle()))) {
    if(fUseMCInfo && fUsePythia){
          AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(track);
          bool isInj = IsTrackInjected(aodtrack, fMCHeader, fMCarray);
          if(!isInj) continue;
    }
   
    if (allDaughters.Remove(track) == 0) {
      new ((*coll)[n]) AliEmcalParticle(track);
      n++;
      AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is included", tracks->GetCurrentID(), track->Pt(), track->Eta(), track->Phi()));
    }
    else {
      AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is excluded", tracks->GetCurrentID(), track->Pt(), track->Eta(), track->Phi()));
    }
  }
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskSEDmesonsFilterCJ::AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters)
{
  // Add all the dauthers of cand in an array. Follows all the decay cascades.

  Int_t n = cand->GetNDaughters();

  //Printf("AddDaughters: the number of dauhters is %d", n);

  Int_t ntot = 0;
  Double_t pt = 0;
  for (Int_t i = 0; i < n; i++) {
    AliVTrack* track = dynamic_cast<AliVTrack*>(cand->GetDaughter(i));
    if (!track) continue;

    AliAODRecoDecay* cand2 = dynamic_cast<AliAODRecoDecay*>(track);

    if (cand2) {
      //Printf("Daughter pT = %.3f --> ", track->Pt());
      pt += AddDaughters(cand2, daughters);
    }
    else {
      if (!track->InheritsFrom("AliAODTrack")) {
        Printf("Warning: One of the daughters is not of type 'AliAODTrack' nor 'AliAODRecoDecay'.");
        continue;
      }
      //Printf("Daughter pT = %.3f", track->Pt());
      daughters.AddLast(track);
      pt += track->Pt();
      ntot++;
    }
  }

  //Printf("Total pt of the daughters = %.3f", pt);

  return pt;
}
//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::AddMCEventTracks(TClonesArray* coll, AliParticleContainer* mctracks)
{
    //
    // Add event tracks to a collection that already contains the D candidates, excluding the daughters of the D candidates
    //

    if (!mctracks) return;

    TObjArray allMCDaughters(10);
    allMCDaughters.SetOwner(kFALSE);

    AliAODMCParticle* mcD = (AliAODMCParticle*)coll->At(0);

    AddMCDaughters(mcD,allMCDaughters,fMCarray);

    mctracks->ResetCurrentID();
    AliAODMCParticle* mcpart = 0;
    Int_t n = coll->GetEntriesFast();
    while ((mcpart = static_cast<AliAODMCParticle*>(mctracks->GetNextAcceptParticle()))) {
        if(TMath::Abs(mcpart->Charge())==0) continue;
        if(fUsePythia){
          bool isInj = IsMCTrackInjected(mcpart, fMCHeader, fMCarray);
          if(!isInj) continue;
        }
        if (allMCDaughters.Remove(mcpart) == 0) {
            new ((*coll)[n]) AliAODMCParticle(*mcpart);
            n++;
            AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is included", mctracks->GetCurrentID(), mcpart->Pt(), mcpart->Eta(), mcpart->Phi()));
        }
        else {
            AliDebug(2, Form("Track %d (pT = %.3f, eta = %.3f, phi = %.3f) is excluded", mctracks->GetCurrentID(), mcpart->Pt(), mcpart->Eta(), mcpart->Phi()));
        }
    }
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskSEDmesonsFilterCJ::AddMCDaughters(AliAODMCParticle* mcDmeson, TObjArray& mcdaughters, TClonesArray* mcArray)
{
    // Add all the dauthers of cand in an array. Follows all the decay cascades.

    Int_t n = mcDmeson->GetNDaughters();

    //Printf("AddDaughters: the number of dauhters is %d", n);
    Double_t pt = 0;

    for (Int_t i = 0; i < n; i++) {
        AliAODMCParticle* DDaughter = static_cast<AliAODMCParticle*>(mcArray->At(mcDmeson->GetDaughter(i)));
        if (!DDaughter) continue;

        if (DDaughter->GetNDaughters()>0) {
            //Printf("Daughter pT = %.3f --> ", track->Pt());
            pt += AddMCDaughters(DDaughter, mcdaughters, mcArray);
            mcdaughters.AddLast(DDaughter);
        }
        else {
            //Printf("Daughter pT = %.3f", track->Pt());
            mcdaughters.AddLast(DDaughter);
            pt += DDaughter->Pt();
        }
    }

    //Printf("Total pt of the daughters = %.3f", pt);

    return pt;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(AliAODRecoDecay* cand, TClonesArray* mcArray)
{
  // Checks whether the mother of the D meson candidate comes from a charm or a bottom quark.

  if (!mcArray) return -1;

  Int_t labDau0 = static_cast<AliVTrack*>(cand->GetDaughter(0))->GetLabel();
  if (labDau0 < 0) return -1;

  AliAODMCParticle* part = static_cast<AliAODMCParticle*>(mcArray->At(labDau0));
  return CheckOrigin(part, mcArray);
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(AliAODRecoDecay* cand, AliStack* stack)
{
  // Checks whether the mother of the D meson candidate comes from a charm or a bottom quark.

  if (!stack) return -1;

  Int_t labDau0 = static_cast<AliVTrack*>(cand->GetDaughter(0))->GetLabel();
  if (labDau0 < 0) return -1;

  return CheckOrigin(labDau0, stack);
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray)
{
  // Checks whether the mother of the particle comes from a charm or a bottom quark.

  if (!part) return -1;
  if (!mcArray) return -1;

  Int_t pdgGranma = 0;
  Int_t mother = part->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma = 0;
  Bool_t isFromB = kFALSE;
  Bool_t isQuarkFound = kFALSE;

  while (mother >= 0) {
    istep++;
    AliAODMCParticle* mcGranma = static_cast<AliAODMCParticle*>(mcArray->At(mother));
    if (mcGranma >= 0) {
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
        isFromB = kTRUE;
      }

      if (abspdgGranma == 4 || abspdgGranma == 5) isQuarkFound = kTRUE;
      mother = mcGranma->GetMother();
    }
    else {
      ::Error("AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin", "Could not retrieve mother particle %d!", mother);
      break;
    }
  }

  if (isQuarkFound) {
    if (isFromB) {
      return kFromBottom;
    }
    else {
      return kFromCharm;
    }
  }
  else {
    return kQuarkNotFound;
  }
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(Int_t ipart, AliStack* stack)
{
  // Checks whether the mother of the particle comes from a charm or a bottom quark.

  if (!stack) return -1;

  TParticle* part = stack->Particle(ipart);
  if (!part) return -1;

  Int_t pdgGranma = 0;
  Int_t mother = part->GetFirstMother();
  Int_t istep = 0;
  Int_t abspdgGranma = 0;
  Bool_t isFromB = kFALSE;
  Bool_t isQuarkFound = kFALSE;

  while (mother >= 0) {
    istep++;
    TParticle* mcGranma = stack->Particle(mother);
    if (mcGranma >= 0) {
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
        isFromB = kTRUE;
      }

      if (abspdgGranma == 4 || abspdgGranma == 5) isQuarkFound = kTRUE;
      mother = mcGranma->GetFirstMother();
    }
    else {
      ::Error("AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin", "Could not retrieve mother particle %d!", mother);
      break;
    }
  }

  if (isQuarkFound) {
    if (isFromB) {
      return kFromBottom;
    }
    else {
      return kFromCharm;
    }
  }
  else {
    return kQuarkNotFound;
  }
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckDecayChannel(AliAODMCParticle* part, TClonesArray* mcArray)
{
  // Determine the decay channel

  if (!part) return -1;
  if (!mcArray) return -1;

  Int_t decay = kDecayOther;

  Int_t absPdgPart = TMath::Abs(part->GetPdgCode());

  if (part->GetNDaughters() == 2) {

    AliAODMCParticle* d1 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughter(0)));
    AliAODMCParticle* d2 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughter(1)));

    if (!d1 || !d2) {
      return decay;
    }

    Int_t absPdg1 = TMath::Abs(d1->GetPdgCode());
    Int_t absPdg2 = TMath::Abs(d2->GetPdgCode());

    if (absPdgPart == 421) { // D0 -> K pi

      if ((absPdg1 == 211 && absPdg2 == 321) || // pi K
          (absPdg1 == 321 && absPdg2 == 211)) { // K pi
        decay = kDecayD0toKpi;
      }
    }

    if (absPdgPart == 413) { // D* -> D0 pi

      if (absPdg1 == 421 && absPdg2 == 211) {  // D0 pi
        Int_t D0decay = CheckDecayChannel(d1, mcArray);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }

      if (absPdg1 == 211 && absPdg2 == 421) {  // pi D0
        Int_t D0decay = CheckDecayChannel(d2, mcArray);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }
    }
  }

  return decay;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDmesonsFilterCJ::CheckDecayChannel(Int_t ipart, AliStack* stack)
{
  // Determine the decay channel

  if (!stack) return -1;

  TParticle* part = stack->Particle(ipart);

  if (!part) return -1;

  Int_t decay = kDecayOther;

  if (part->GetNDaughters() == 2) {

    Int_t id1 = part->GetDaughter(0);
    Int_t id2 = part->GetDaughter(1);

    TParticle* d1 = stack->Particle(id1);
    TParticle* d2 = stack->Particle(id2);

    if (!d1 || !d2) {
      return decay;
    }

    Int_t absPdg1 = TMath::Abs(d1->GetPdgCode());
    Int_t absPdg2 = TMath::Abs(d2->GetPdgCode());


    if (part->GetPdgCode() == 421) { // D0 -> K pi

      if ((absPdg1 == 211 && absPdg2 == 321) || // pi K
          (absPdg1 == 321 && absPdg2 == 211)) { // K pi
        decay = kDecayD0toKpi;
      }
    }

    if (part->GetPdgCode() == 413) { // D* -> D0 pi

      if (absPdg1 == 421 && absPdg2 == 211) {  // D0 pi
        Int_t D0decay = CheckDecayChannel(id1, stack);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }

      if (absPdg1 == 211 && absPdg2 == 421) {  // pi D0
        Int_t D0decay = CheckDecayChannel(id2, stack);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }
    }
  }
  
  return decay;
  
}


//_____________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  /// method to check if a track comes from a given generator

  Int_t lab=TMath::Abs(track->GetLabel());
  nameGen=AliVertexingHFUtils::GetGenerator(lab,header);
  
  //  Int_t countControl=0;
  
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){
      printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=AliVertexingHFUtils::GetGenerator(mother,header);
    // countControl++;
    // if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
    //   printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
    //   break;
    // }
  }
  
  return;
}
//----------------------------------------------------------------------
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::IsTrackInjected(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a track comes from the signal event or from the underlying Hijing event
  TString nameGen;
  GetTrackPrimaryGenerator(track,header,arrayMC,nameGen);
  
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  
  return kTRUE;
}


//_____________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::GetMCTrackPrimaryGenerator(AliAODMCParticle *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  /// method to check if a track comes from a given generator

  Int_t lab=TMath::Abs(track->GetLabel());
  nameGen=AliVertexingHFUtils::GetGenerator(lab,header);
  
  //  Int_t countControl=0;
  
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){
      printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=AliVertexingHFUtils::GetGenerator(mother,header);
    // countControl++;
    // if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
    //   printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
    //   break;
    // }
  }
  
  return;
}
//----------------------------------------------------------------------
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::IsMCTrackInjected(AliAODMCParticle *track,AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a track comes from the signal event or from the underlying Hijing event
  TString nameGen;
  GetMCTrackPrimaryGenerator(track,header,arrayMC,nameGen);
  
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  
  return kTRUE;
}
