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
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include "AliAodSkimTask.h"

using namespace std;
ClassImp(AliAodSkimTask)

AliAodSkimTask::AliAodSkimTask() : 
  AliAnalysisTaskSE(), fClusMinE(-1), fCutMC(1), fYCutMC(0.7),
  fDoCopyHeader(1),  fDoCopyVZERO(1),  fDoCopyTZERO(1),  fDoCopyVertices(1),  fDoCopyTOF(1), fDoCopyTracks(1), fDoCopyTrigger(1), fDoCopyPTrigger(0), 
  fDoCopyCells(1), fDoCopyPCells(0), fDoCopyClusters(1), fDoCopyDiMuons(0), fDoCopyZDC(1), fDoCopyMC(1), fDoCopyMCHeader(1), fTrials(0), fPyxsec(0), 
  fPytrials(0), fPypthardbin(0), fAOD(0), fAODMcHeader(0), fOutputList(0), fHevs(0), fHclus(0)
{
}

AliAodSkimTask::AliAodSkimTask(const char* name) : 
  AliAnalysisTaskSE(name), fClusMinE(-1), fCutMC(1), fYCutMC(0.7),
  fDoCopyHeader(1),  fDoCopyVZERO(1),  fDoCopyTZERO(1),  fDoCopyVertices(1),  fDoCopyTOF(1), fDoCopyTracks(1), fDoCopyTrigger(1), fDoCopyPTrigger(0), 
  fDoCopyCells(1), fDoCopyPCells(0), fDoCopyClusters(1), fDoCopyDiMuons(0), fDoCopyZDC(1), fDoCopyMC(1), fDoCopyMCHeader(1), fTrials(0), fPyxsec(0), 
  fPytrials(0), fPypthardbin(0), fAOD(0), fAODMcHeader(0), fOutputList(0), fHevs(0), fHclus(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAodSkimTask::~AliAodSkimTask()
{
  if (fOutputList) {
    delete fOutputList;
  }
  delete fHevs;
  delete fHclus;
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

  Bool_t store = kFALSE;
  if (fClusMinE>0) {
    TClonesArray *cls  = fAOD->GetCaloClusters();  
    for (Int_t i=0; i<cls->GetEntriesFast(); ++i) {
      AliAODCaloCluster *clus = static_cast<AliAODCaloCluster*>(cls->At(i));
      if (!clus->IsEMCAL())
	continue;
      Double_t e = clus->E();
      fHclus->Fill(e);
      if (e>fClusMinE) {
	store = kTRUE;
      }
    }
  } else { 
    store = kTRUE;
  }

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
  if (tout) {
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
      AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
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
      AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
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
      AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
    }
    out->AbsorbObjects(in);
  }

  if (fDoCopyDiMuons) { 
    TClonesArray *out = eout->GetDimuons();
    TClonesArray *in  = evin->GetDimuons();
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
	AliFatal(Form("%s: Previous event not deleted. This should not happen!",GetName()));
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

void AliAodSkimTask::Terminate(Option_t *)
{
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man->GetAnalysisType()!=0)
     return;
   if (fHevs==0)
     return;
   cout << "AliAodSkimTask " << GetName() << " terminated with accepted fraction of events: " << fHevs->GetBinContent(2)/fHevs->GetEntries() 
	<< " (" << fHevs->GetBinContent(2) << "/" << fHevs->GetEntries() << ")" << endl;
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

const char *AliAodSkimTask::Str() const
{
  return Form("mine%.2f_%dycut%.2f_%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
	      fClusMinE,
	      fCutMC,
	      fYCutMC,
	      fDoCopyHeader,
	      fDoCopyVZERO,
	      fDoCopyTZERO,
	      fDoCopyVertices,
	      fDoCopyTOF,
	      fDoCopyTracks,
	      fDoCopyTrigger,
	      fDoCopyPTrigger,
	      fDoCopyCells,
	      fDoCopyPCells,
	      fDoCopyClusters,
	      fDoCopyDiMuons,
	      fDoCopyZDC,
	      fDoCopyMC,
	      fDoCopyMCHeader);
}
