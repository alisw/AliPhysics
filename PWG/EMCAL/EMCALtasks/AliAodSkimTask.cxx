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
  AliAnalysisTaskSE(name), fClusMinE(-1), fTrackMinPt(-1), fTrackMaxPt(-1), fDoBothMinTrackAndClus(0), fCutMC(1), fYCutMC(0.7), fCutMinPt(0), fCutFilterBit(-1), fGammaBr(""),
  fDoCopyHeader(1),  fDoCopyVZERO(1),  fDoCopyTZERO(1),  fDoCopyVertices(1),  fDoCopyTOF(1), fDoCopyTracklets(1), fDoCopyTracks(1), fDoRemoveTracks(0), fDoCleanTracks(0),
  fDoRemCovMat(0), fDoRemPid(0), fDoCopyTrigger(1), fDoCopyPTrigger(0), fDoCopyCells(1), fDoCopyPCells(0), fDoCopyClusters(1), fDoCopyDiMuons(0),  fDoCopyTrdTracks(0),
  fDoCopyV0s(0), fDoCopyCascades(0), fDoCopyZDC(1), fDoCopyConv(0), fDoCopyMC(1), fDoCopyMCHeader(1), fDoVertWoRefs(0), fDoVertMain(0), fDoCleanTracklets(0), fDoCopyUserTree(0),
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
      if (!clus->IsEMCAL())
        continue;
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
  return Form("mine%.2f_%dycut%.2f_%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
              fClusMinE,
              fCutMC,
              fYCutMC,
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
  AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
  AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
  TTree *tout = oh->GetTree();
  if (fDoCopyUserTree) {
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

  if (fDoCopyHeader) {
    AliAODHeader *out = (AliAODHeader*)eout->GetHeader();
    AliAODHeader *in  = (AliAODHeader*)evin->GetHeader();
    *out = *in;
    out->SetUniqueID(fTrials);
  }

  if (fDoCopyVZERO) {
    AliAODVZERO *out = eout->GetVZEROData();
    AliAODVZERO *in  = evin->GetVZEROData();
    *out = *in;
  }

  if (fDoCopyTZERO) {
    AliAODTZERO *out = eout->GetTZEROData();
    AliAODTZERO *in  = evin->GetTZEROData();
    *out = *in;
  }

  if (fDoCopyVertices) {
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
        if (!tmp.Contains("PrimaryVertex")&&!tmp.Contains("SPDVertex")&&!tmp.Contains("TPCVertex"))
          continue;
        marked=i;
        break;
      }
    }
    if (marked>0) {
      out->RemoveRange(marked,out->GetEntries());
      out->Compress();
    }
  }
  
  if (fDoCopyTOF) {
    AliTOFHeader *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader());
    const AliTOFHeader *in = evin->GetTOFHeader();
    *out = *in;
  }

  if (fDoCopyTracks) {
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

  if (fDoCopyTracklets) {
    AliAODTracklets *out = eout->GetTracklets();
    AliAODTracklets *in  = evin->GetTracklets();
    *out = *in;
    if (fDoCleanTracklets) {
      Int_t n=in->GetNumberOfTracklets();
      out->SetTitle(Form("Ntracklets=%d",n));
      out->DeleteContainer();
    }
  }

  if (fDoCopyTrigger) {
    AliAODCaloTrigger *out = eout->GetCaloTrigger("EMCAL");
    AliAODCaloTrigger *in  = evin->GetCaloTrigger("EMCAL");
    *out = *in;
  }

  if (fDoCopyPTrigger) {
    AliAODCaloTrigger *out = eout->GetCaloTrigger("PHOS");
    AliAODCaloTrigger *in  = evin->GetCaloTrigger("PHOS");
    *out = *in;
  }

  if (fDoCopyCells) {
    AliAODCaloCells *out = eout->GetEMCALCells();
    AliAODCaloCells *in  = evin->GetEMCALCells();
      *out = *in;
  }

  if (fDoCopyPCells) {
    AliAODCaloCells *out = eout->GetPHOSCells();
    AliAODCaloCells *in  = evin->GetPHOSCells();
    *out = *in;
  }

  if (fDoCopyClusters) {
    TClonesArray *out = eout->GetCaloClusters();
    TClonesArray *in  = evin->GetCaloClusters();
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous clusters not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyTrdTracks) {
    TClonesArray *out = static_cast<TClonesArray*>(eout->FindListObject("trdTracks"));
    TClonesArray *in  = static_cast<TClonesArray*>(eout->FindListObject("trdTracks"));
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous trdtracks not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyV0s) {
    TClonesArray *out = eout->GetV0s();
    TClonesArray *in  = evin->GetV0s();
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous v0s not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyCascades) {
    TClonesArray *out = eout->GetCascades();
    TClonesArray *in  = evin->GetCascades();
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyZDC) {
    AliAODZDC *out = eout->GetZDCData();
    AliAODZDC *in  = evin->GetZDCData();
    *out = *in;
  }

  if (fDoCopyDiMuons) {
    TClonesArray *out = eout->GetDimuons();
    TClonesArray *in  = evin->GetDimuons();
    if (out->GetEntries()>0) { // just checking if the deletion of previous event worked
      AliFatal(Form("%s: Previous dimuons not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyConv) {
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

  if (fDoCopyMC) {
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
      if (fCutMC) {
        for (Int_t i=0;i<out->GetEntriesFast();++i) {
          AliAODMCParticle *mc = static_cast<AliAODMCParticle*>(in->At(i));
          if ((mc==0)&&(i==0)) {
            AliError(Form("%s: No MC info, skipping this event!",GetName()));
            oh->SetFillAOD(kFALSE);
            return;
          }
          if ((mc==0)||(TMath::Abs(mc->Y())>fYCutMC))
            new ((*out)[i]) AliAODMCParticle;
        }
      }
    }
  }

  if (fDoCopyMCHeader) {
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

  if (gDebug>10) {
    Int_t  run = eout->GetRunNumber();
    AliAODVertex *v=(AliAODVertex*)eout->GetVertices()->At(0);
    Int_t   vzn = v->GetNContributors();
    Double_t vz = v->GetZ();
    cout << "debug run " << run << " " << vzn << " " << vz << endl;
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

