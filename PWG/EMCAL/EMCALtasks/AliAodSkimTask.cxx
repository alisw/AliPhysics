#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSystem.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODMCHeader.h>
#include <AliAODMCParticle.h>
#include <AliAODVertex.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include "AliAodSkimTask.h"
#include "TObjectTable.h"
using namespace std;
ClassImp(AliAodSkimTask)

AliAodSkimTask::AliAodSkimTask(const char* name) :
  AliAnalysisTaskSE(name), fClusMinE(-1), fTrackMinPt(-1), fTrackMaxPt(-1), fDoBothMinTrackAndClus(0), fCutMinPt(0), fCutFilterBit(-1), fGammaBr(""),
  fDoCopyHeader(1),  fDoCopyVZERO(1),  fDoCopyTZERO(1),  fDoCopyVertices(1),  fDoCopyTOF(1), fDoCopyTracklets(1), fDoCopyTracks(1), fDoRemoveTracks(0), fDoCleanTracks(0),
  fDoRemCovMat(0), fDoRemPid(0), fDoCopyTrigger(1), fDoCopyPTrigger(0), fDoCopyCells(1), fDoCopyPCells(0), fDoCopyClusters(1), fDoCopyDiMuons(0),  fDoCopyTrdTracks(0),
  fDoCopyV0s(0), fDoCopyCascades(0), fDoCopyZDC(1), fDoCopyConv(0), fDoCopyKinks(0), fDoCopyMC(1), fDoCopyMCHeader(1), fDoVertWoRefs(0), fDoVertMain(0), fDoCleanTracklets(0),
  fDoCopyUserTree(0), fDoPhosFilt(0), fDoRemoveMcParts(0), fCutMcIsPrimary(0), fCutMcIsPhysPrimary(0), fCutMcPt(-1), fCutMcY(1.0), fCutMcPhos(0), fCutMcEmcal(0),
  fTrials(0), fPyxsec(0), fPytrials(0), fPypthardbin(0), fAOD(0), fAODMcHeader(0), fOutputList(0), fHevs(0), fHclus(0), fHtrack(0)
{
  if (name) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
  }
}

AliAodSkimTask::~AliAodSkimTask()
{
  if (fOutputList) {
    delete fOutputList;
  }
  delete fHevs;
  delete fHclus;
  delete fHtrack;
}

void AliAodSkimTask::CleanTrack(AliAODTrack *t)
{
  if (!fDoCleanTracks)
    return;
  //cout << "Clean tracks " << endl;
  t->SetRAtAbsorberEnd(0);
  t->SetChi2MatchTrigger(0);
  t->SetMuonClusterMap(0);
  t->SetITSMuonClusterMap(0);
  t->SetMUONtrigHitsMapTrg(0);
  t->SetMUONtrigHitsMapTrk(0);
  t->SetMFTClusterPattern(0);
  t->SetMatchTrigger(0);
  t->SetIsMuonGlobalTrack(0);
  t->SetXYAtDCA(-999., -999.);
  t->SetPxPyPzAtDCA(-999., -999., -999.);
  if (fDoRemCovMat)
    t->RemoveCovMatrix();
  if (fDoRemPid) {
    AliAODPid *pid = t->GetDetPid();
    delete pid;
    t->SetDetPID(0);
  } else if (t->GetDetPid()) {
    AliAODPid *pid = t->GetDetPid();
    AliAODPid *nid = new AliAODPid;
    nid->SetTPCsignal(pid->GetTPCsignal());
    nid->SetTPCsignalN(pid->GetTPCsignalN());
    nid->SetTPCmomentum(pid->GetTPCmomentum());
    nid->SetTPCTgl(pid->GetTPCTgl());
    /* Not used and getter not implemented
       AliTPCdEdxInfo *dedx = new AliTPCdEdxInfo(pid->GetTPCdEdxInfo());
       nid->SetTPCdEdxInfo(dedx); */
    nid->SetTOFsignal(pid->GetTOFsignal());
    Double_t val[5];
    pid->GetTOFpidResolution(val);
    nid->SetTOFpidResolution(val);
    pid->GetIntegratedTimes(val,5);
    nid->SetIntegratedTimes(val);
    delete pid;
    t->SetDetPID(nid);
  }
}

void AliAodSkimTask::CopyCascades()
{
  if (!fDoCopyCascades)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetCascades();
  TClonesArray *in  = evin->GetCascades();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyCells()
{
  if (!fDoCopyCells)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloCells *out = eout->GetEMCALCells();
  AliAODCaloCells *in  = evin->GetEMCALCells();
  *out = *in;
}

void AliAodSkimTask::CopyCellsP()
{
  if (!fDoCopyPCells)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloCells *out = eout->GetPHOSCells();
  AliAODCaloCells *in  = evin->GetPHOSCells();
  *out = *in;
}

void AliAodSkimTask::CopyClusters()
{
  if (!fDoCopyClusters)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetCaloClusters();
  TClonesArray *in  = evin->GetCaloClusters();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous clusters not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyConv()
{
  if (!fDoCopyConv)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = dynamic_cast<TClonesArray*>(eout->FindListObject(fGammaBr));
  TClonesArray *in  = dynamic_cast<TClonesArray*>(evin->FindListObject(fGammaBr));
  if (!in) {
    evin->GetList()->ls();
    AliFatal(Form("%s: Could not find conversion branch with name %s!",GetName(), fGammaBr.Data()));
  }
  if (in && !out) {
    out = new TClonesArray("AliAODConversionPhoton",2*in->GetEntries());
    out->SetName(fGammaBr);
    oh->AddBranch("TClonesArray", &out);
  }
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    out->Delete();
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyDimuons()
{
  if (!fDoCopyDiMuons)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetDimuons();
  TClonesArray *in  = evin->GetDimuons();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous dimuons not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyHeader()
{
  if (!fDoCopyHeader)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODHeader *out = (AliAODHeader*)eout->GetHeader();
  AliAODHeader *in  = (AliAODHeader*)evin->GetHeader();
  *out = *in;
  out->SetUniqueID(fTrials);
}

void AliAodSkimTask::CopyKinks()
{
  if (!fDoCopyKinks)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetKinks();
  TClonesArray *in  = evin->GetKinks();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous kinks not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyMc()
{
  if (!fDoCopyMC)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
  TClonesArray *in  = static_cast<TClonesArray*>(evin->FindListObject(AliAODMCParticle::StdBranchName()));
  if (in && !out) {
    fgAODMCParticles = new TClonesArray("AliAODMCParticle",2*in->GetEntries());
    fgAODMCParticles->SetName(AliAODMCParticle::StdBranchName());
    oh->AddBranch("TClonesArray", &fgAODMCParticles);
    out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
  }
  if (in && out) {
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous mcparticles not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
    for (Int_t i=0;i<out->GetEntriesFast();++i) {
      AliAODMCParticle *mc = static_cast<AliAODMCParticle*>(out->At(i));
      if (!KeepMcPart(mc))
        new ((*out)[i]) AliAODMCParticle;
    }
  }
}

void AliAodSkimTask::CopyMcHeader()
{
  if (!fDoCopyMCHeader)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODMCHeader *out = static_cast<AliAODMCHeader*>(eout->FindListObject(AliAODMCHeader::StdBranchName()));
  AliAODMCHeader *in  = static_cast<AliAODMCHeader*>(evin->FindListObject(AliAODMCHeader::StdBranchName()));
  if (in && !out) {
    fAODMcHeader = new AliAODMCHeader();
    fAODMcHeader->SetName(AliAODMCHeader::StdBranchName());
    oh->AddBranch("AliAODMCHeader",&fAODMcHeader);
    out = static_cast<AliAODMCHeader*>(eout->FindListObject(AliAODMCHeader::StdBranchName()));
  }
  if (in && out) {
    *out = *in;
    if ((in->GetCrossSection()==0) && (fPyxsec>0)) {
      out->SetCrossSection(fPyxsec);
      out->SetTrials(fPytrials);
      out->SetPtHard(fPypthardbin);
    }
  }
}

void AliAodSkimTask::CopyUserTree()
{
  if (!fDoCopyUserTree)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  TTree *tout = oh->GetTree();
  TList *lout = tout->GetUserInfo();
  if (lout->FindObject("alirootVersion")==0) {
    TList *lin = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetUserInfo();
    TString apver(gSystem->BaseName(gSystem->Getenv("ALICE_PHYSICS")));
    lout->Add(new TObjString(Form("AodSkim: ver %s, tag %s with settings %s",GetVersion(),apver.Data(),Str())));
    AliInfo(Form("%s: Set user info %s", GetName(), lout->At(0)->GetName()));
    for (Int_t jj=0;jj<lin->GetEntries()-1;++jj) {
      lout->Add(lin->At(jj)->Clone(lin->At(jj)->GetName()));
    }
  }
}

void AliAodSkimTask::CopyTof()
{
  if (!fDoCopyTOF)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliTOFHeader *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader());
  const AliTOFHeader *in = evin->GetTOFHeader();
  *out = *in;
}

void AliAodSkimTask::CopyTracklets()
{
  if (!fDoCopyTracklets)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODTracklets *out = eout->GetTracklets();
  AliAODTracklets *in  = evin->GetTracklets();
  *out = *in;
  if (fDoCleanTracklets) {
    Int_t n=in->GetNumberOfTracklets();
    out->SetTitle(Form("Ntracklets=%d",n));
    out->DeleteContainer();
  }
}

void AliAodSkimTask::CopyTracks()
{
  if (!fDoCopyTracks)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetTracks();
  TClonesArray *in  = evin->GetTracks();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous tracks not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
  for (Int_t i=0;i<out->GetEntries();++i) {
    AliAODTrack *t = static_cast<AliAODTrack*>(out->At(i));
    if (KeepTrack(t)) {
      CleanTrack(t);
      if (fDoVertMain)
        t->SetProdVertex(0);
    } else {
      t->~AliAODTrack();
      new ((*out)[i]) AliAODTrack;
    }
  }
}

void AliAodSkimTask::CopyTrdTracks()
{
  if (!fDoCopyTrdTracks)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = static_cast<TClonesArray*>(eout->FindListObject("trdTracks"));
  TClonesArray *in  = static_cast<TClonesArray*>(evin->FindListObject("trdTracks"));
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous trdtracks not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyTrigger()
{
  if (!fDoCopyTrigger)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloTrigger *out = eout->GetCaloTrigger("EMCAL");
  AliAODCaloTrigger *in  = evin->GetCaloTrigger("EMCAL");
  *out = *in;
}

void AliAodSkimTask::CopyTriggerP()
{
  if (!fDoCopyPTrigger)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloTrigger *out = eout->GetCaloTrigger("PHOS");
  AliAODCaloTrigger *in  = evin->GetCaloTrigger("PHOS");
  *out = *in;
}

void AliAodSkimTask::CopyTZero()
{
  if (!fDoCopyTZERO)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODTZERO *out = eout->GetTZEROData();
  AliAODTZERO *in  = evin->GetTZEROData();
  *out = *in;
}

void AliAodSkimTask::CopyV0s()
{
  if (!fDoCopyV0s)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetV0s();
  TClonesArray *in  = evin->GetV0s();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous v0s not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
}

void AliAodSkimTask::CopyVertices()
{
  if (!fDoCopyVertices)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TClonesArray *out = eout->GetVertices();
  TClonesArray *in  = evin->GetVertices();
  if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
    AliFatal(Form("%s: Previous vertices not deleted. This should not happen!",GetName()));
  }
  out->AbsorbObjects(in);
  Int_t marked=-1;
  for (Int_t i=0; i<out->GetEntries(); ++i) {
    AliAODVertex *v = static_cast<AliAODVertex*>(out->At(i));
    Int_t nc = v->CountRealContributors();
    Int_t nd = v->GetNDaughters();
    if (fDoVertWoRefs) {
      if (nc>0)
        v->SetNContributors(nc);
      else
        v->SetNContributors(nd);
    }
    if (fDoVertMain) {
      TString tmp(v->GetName());
      if (tmp.Length()==0)
        continue;
      if (!tmp.Contains("PrimaryVertex")&&!tmp.Contains("SPDVertex")&&!tmp.Contains("TPCVertex"))
        continue;
      marked=i+1;
    }
  }
  if (marked>0) {
    out->RemoveRange(marked,out->GetEntries()-1);
    out->Compress();
  }
}

void AliAodSkimTask::CopyVZero()
{
  if (!fDoCopyVZERO)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODVZERO *out = eout->GetVZEROData();
  AliAODVZERO *in  = evin->GetVZEROData();
  *out = *in;
}

void AliAodSkimTask::CopyZdc()
{
  if (!fDoCopyZDC)
    return;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODZDC *out = eout->GetZDCData();
  AliAODZDC *in  = evin->GetZDCData();
  *out = *in;
}

Bool_t AliAodSkimTask::IsDcalAcc(Double_t phi, Double_t eta)
{
  const Double_t etamin=0.22*0.9;
  const Double_t etamax=0.7*1.1;
  const Double_t phimin=260./360*TMath::TwoPi()*0.9;
  const Double_t phimax=327./360*TMath::TwoPi()*1.1;
  const Double_t phiin=TVector2::Phi_0_2pi(phi);
  if ((eta>etamin)&&(eta<etamax) && (phiin>phimin) && (phiin<phimax))
    return kTRUE;
  const Double_t phimin2=320./360*TMath::TwoPi()*0.9;
  if ((TMath::Abs(eta)<etamax) && (phiin>phimin2) && (phiin<phimax))
    return kTRUE;
  return kFALSE;
}

Bool_t AliAodSkimTask::IsEmcalAcc(Double_t phi, Double_t eta)
{
  const Double_t etaabs=0.7*1.1;
  const Double_t phimin=80./360*TMath::TwoPi()*0.9;
  const Double_t phimax=187./360*TMath::TwoPi()*1.1;
  const Double_t phiin=TVector2::Phi_0_2pi(phi);
  if ((TMath::Abs(eta)<etaabs) && (phiin>phimin) && (phiin<phimax))
    return kTRUE;
  return kFALSE;
}

Bool_t AliAodSkimTask::IsPhosAcc(Double_t phi, Double_t eta)
{
  const Double_t etaabs=0.12*1.1;
  const Double_t phimin=250./360*TMath::TwoPi()*0.9;
  const Double_t phimax=320./360*TMath::TwoPi()*1.1;
  const Double_t phiin=TVector2::Phi_0_2pi(phi);
  if ((TMath::Abs(eta)<etaabs) && (phiin>phimin) && (phiin<phimax))
    return kTRUE;
  return kFALSE;
}

Bool_t AliAodSkimTask::KeepMcPart(AliAODMCParticle *p)
{
  if (!fDoRemoveMcParts)
    return kTRUE;
  if (fCutMcIsPrimary && !p->IsPrimary())
    return kFALSE;
  if (fCutMcIsPhysPrimary && !p->IsPhysicalPrimary())
    return kFALSE;
  Double_t pt  = p->Pt();
  if (pt<fCutMcPt)
    return kFALSE;
  Double_t y   = p->Y();
  if (TMath::Abs(y)>fCutMcY)
    return kFALSE;
  Double_t phi = p->Phi();
  Double_t eta = p->Eta();
  if (fCutMcPhos) {
    if (!IsPhosAcc(phi,eta))
      return kFALSE;
  }
  if (fCutMcEmcal) {
    if (!IsEmcalAcc(phi,eta))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t AliAodSkimTask::KeepTrack(AliAODTrack *t)
{
  if (!fDoRemoveTracks)
    return kTRUE;
  //cout << "Keep tracks " << endl;
  if (t->IsMuonGlobalTrack())
    return kFALSE;
  if (t->Pt()<fCutMinPt)
    return kFALSE;
  if (fCutFilterBit!=(UInt_t)-1) {
    if (t->TestFilterBit(fCutFilterBit)==0)
      return kFALSE;
  }
  return kTRUE;
}

Bool_t AliAodSkimTask::PythiaInfoFromFile(const char* currFile, Float_t &xsec, Float_t &trials, Int_t &pthard)
{
  TString file(currFile);
  xsec = 0;
  trials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1,Form("File name: %s",file.Data()));

  // Get the pt hard bin
  TString strPthard(file);

  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec()) {
    pthard = strPthard.Atoi();
  }
  else {
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));
    pthard = -1;
  }

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root"));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = static_cast<TKey*>(fxsec->GetListOfKeys()->At(0));
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      xsec = static_cast<TProfile*>(list->FindObject("h1Xsec"))->GetBinContent(1);
      trials = static_cast<TH1F*>(list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = static_cast<TTree*>(fxsec->Get("Xsection"));
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

Bool_t AliAodSkimTask::SelectEvent()
{
  // Accept event only if an EMCal-cluster with a minimum energy of fClusMinE has been found in the event
  Bool_t storeE = kFALSE;
  if (fClusMinE>0) {
    TClonesArray *cls  = fAOD->GetCaloClusters();
    for (Int_t i=0; i<cls->GetEntriesFast(); ++i) {
      AliAODCaloCluster *clus = static_cast<AliAODCaloCluster*>(cls->At(i));
      if (fDoPhosFilt && clus->IsEMCAL())
        continue;
      if (!fDoPhosFilt && !clus->IsEMCAL())
        continue;
      if (clus->IsEMCAL() && fCutMcEmcal) {
        Double_t vec[3]={0,0,0};
        TLorentzVector p;
        clus->GetMomentum(p,vec);
        if (!IsEmcalAcc(p.Phi(),p.Eta()))
          continue;
      }
      Double_t e = clus->E();
      fHclus->Fill(e);
      if (e>fClusMinE) {
        storeE = kTRUE;
      }
    }
  } else {
    storeE = kTRUE;
  }

  // Accept event only if an track with a minimum pT of fTrackMinPt has been found in the event
  Bool_t storePt = kFALSE;
  if (fTrackMinPt > 0 ){
    TClonesArray *tracks = fAOD->GetTracks();
    for (Int_t i=0;i<tracks->GetEntries();++i) {
      AliAODTrack *t = static_cast<AliAODTrack*>(tracks->At(i));
      Double_t pt = t->Pt();
      fHtrack->Fill(pt);
      if (pt>fTrackMinPt) {
        storePt = kTRUE;
      }
      if(fTrackMaxPt > 0  &&  fTrackMaxPt < pt){
        storePt = kFALSE;
      } 
    }
  } else {
    storePt = kTRUE;
  }

  Bool_t store = kFALSE;
  if (fDoBothMinTrackAndClus && fClusMinE>0 && fTrackMinPt > 0){
    // request that both conditions are full-filled for propagating the event
    store     = (storeE && storePt);
  } else if (!fDoBothMinTrackAndClus && fClusMinE>0 && fTrackMinPt > 0){
    // request that at least one of the conditions is fullfilled
    store     = (storeE || storePt);
  } else if ( fClusMinE>0 ){
    store     = storeE;
  } else if ( fTrackMinPt>0 ){
    store     = storePt;
  } else {
    store     = kTRUE;
  }

  return store;
}

const char *AliAodSkimTask::Str() const
{
  return Form("%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
              fDoCopyHeader,
              fDoCopyVZERO,
              fDoCopyTZERO,
              fDoCopyVertices,
              fDoCopyTOF,
              fDoCopyTracklets,
              fDoCopyTracks,
              fDoCopyTrigger,
              fDoCopyPTrigger,
              fDoCopyCells,
              fDoCopyPCells,
              fDoCopyClusters,
              fDoCopyDiMuons,
              fDoCopyTrdTracks,
              fDoCopyV0s,
              fDoCopyCascades,
              fDoCopyZDC,
              fDoCopyConv,
              fDoCopyKinks,
              fDoCopyMC,
              fDoCopyMCHeader);
}

void AliAodSkimTask::Terminate(Option_t *)
{
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man->GetAnalysisType()!=0)
     return;
   if (fHevs==0)
     return;
   Int_t norm = fHevs->GetEntries();
   if (norm<1)
     norm=1;
   cout << "AliAodSkimTask " << GetName() << " terminated with accepted fraction of events: " << fHevs->GetBinContent(2)/norm
	<< " (" << fHevs->GetBinContent(2) << "/" << fHevs->GetEntries() << ")" << endl;
}

void AliAodSkimTask::UserCreateOutputObjects()
{
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  if (oh) {
    TFile *fout = oh->GetTree()->GetCurrentFile();
    fout->SetCompressionLevel(2);
  }

  cout << "AliAodSkimTask " << GetName() << " version " << GetVersion() << " running with " << Str() << endl;
      
  fOutputList = new TList;
  fOutputList->SetOwner();
  fHevs = new TH1F("hEvs","",2,-0.5,1.5);
  fOutputList->Add(fHevs);
  fHclus = new TH1F("hClus",";E (GeV)",200,0,100);
  fOutputList->Add(fHclus);
  fHtrack = new TH1F("hTrack",";p_{T} (GeV/c)",200,0,100);
  fOutputList->Add(fHtrack);
  PostData(1, fOutputList);
}

void AliAodSkimTask::UserExec(Option_t *)
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD)
    return;

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  if (!oh) {
    AliFatal(Form("%s: No output handler found", GetName()));
    return;
  }
  oh->SetFillAOD(kFALSE);

  Bool_t store = SelectEvent();
  
  if (!store) {
    ++fTrials;
    fHevs->Fill(0);
    return;
  }

  fHevs->Fill(1);
  oh->SetFillAOD(kTRUE);

  CopyUserTree();
  CopyHeader();
  CopyVZero();
  CopyTZero();
  CopyVertices();
  CopyTof();
  CopyTracklets();
  CopyTracks();
  CopyTrigger();
  CopyTriggerP();
  CopyCells();
  CopyCellsP();
  CopyClusters();
  CopyDimuons();
  CopyTrdTracks();
  CopyV0s();
  CopyCascades();
  CopyZdc();
  CopyConv();
  CopyKinks();
  CopyMc();
  CopyMcHeader();
  CopyMore(); // Overload to do more in derived classes

  if (gDebug>10) {
    AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
    Int_t  run = eout->GetRunNumber();
    AliAODVertex *v=(AliAODVertex*)eout->GetVertices()->At(0);
    Int_t   vzn = v->GetNContributors();
    Double_t vz = v->GetZ();
    cout << "AliAodSkimTask " << GetName() << " debug run " << run << " " << vzn << " " << vz << endl;
  }

  fTrials = 0;
  PostData(1, fOutputList);
}

Bool_t AliAodSkimTask::UserNotify()
{
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();
  Int_t nevents = tree->GetEntriesFast();

  Int_t slevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  Bool_t res = PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);
  gErrorIgnoreLevel=slevel;

  if (res) {
    cout << "AliAodSkimTask " << GetName() << " found xsec info: " << xsection << " " << trials << " " << pthardbin << " " << nevents << endl;
    fPyxsec      = xsection;
    fPytrials    = trials;
    fPypthardbin = pthardbin;
  }
  return res;
}
