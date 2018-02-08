#include <Riostream.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TFile.h>
#include <TSystem.h>
#include <TKey.h>
#include <AliAODMCHeader.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODMCParticle.h>
#include <AliAnalysisManager.h>
#include "AliAodSkimTask.h"
#include <AliLog.h>

using namespace std;
ClassImp(AliAodSkimTask)

AliAodSkimTask::AliAodSkimTask() : AliAnalysisTaskSE(), fClusMinE(-1), fCutMC(1), fTrials(0), fPyxsec(0), fPytrials(0), 
                                   fPypthardbin(0), fAOD(0), fAODMcHeader(0), fOutputList(0)
{
}

AliAodSkimTask::AliAodSkimTask(const char* name) : AliAnalysisTaskSE(name), fClusMinE(-1), fCutMC(1), fTrials(0), fPyxsec(0), fPytrials(0),
                                                   fPypthardbin(0), fAOD(0), fAODMcHeader(0), fOutputList(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAodSkimTask::~AliAodSkimTask()
{
  if (fOutputList) {
    delete fOutputList;
  }
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

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  if (oh)
    oh->SetFillAOD(kFALSE);

  if (1) {
    oh->SetFillAOD(kTRUE);
    AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
    AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
    TTree *tout = oh->GetTree();
    if (tout) {
      TList *lout = tout->GetUserInfo();
      if (lout->FindObject("alirootVersion")==0) {
	TList *lin = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetUserInfo();
	for (Int_t jj=0;jj<lin->GetEntries()-1;++jj) { 
	  lout->Add(lin->At(jj)->Clone(lin->At(jj)->GetName()));
	}
      }
    }
	
    if (1) {
      AliAODHeader *out = (AliAODHeader*)eout->GetHeader();                      
      AliAODHeader *in  = (AliAODHeader*)evin->GetHeader();  	    
      *out = *in;
      out->SetUniqueID(fTrials);
    }
    if (1) {   
      AliAODVZERO *out = eout->GetVZEROData();	                 
      AliAODVZERO *in  = evin->GetVZEROData();
      *out = *in;                  
    }
    if (1) {   
      AliAODTZERO *out = eout->GetTZEROData();	                 
      AliAODTZERO *in  = evin->GetTZEROData(); 	    
      *out = *in; 
    }
    if (1) {   
      TClonesArray *out = eout->GetVertices(); 
      TClonesArray *in  = evin->GetVertices();      
      for (Int_t i=0;i<in->GetEntriesFast();++i) {
	AliAODVertex *v = static_cast<AliAODVertex*>(in->At(i));
	new ((*out)[i]) AliAODVertex(*v);
      }
    }
    if (1) {   
      AliTOFHeader *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader()); 
      const AliTOFHeader *in = evin->GetTOFHeader();	    
      *out = *in;                  
    }
    if (1) {
      TClonesArray *out = eout->GetTracks();	                 
      TClonesArray *in  = evin->GetTracks();	
      for (Int_t i=0;i<in->GetEntriesFast();++i) {
	AliAODTrack *t = static_cast<AliAODTrack*>(in->At(i));
	new ((*out)[i]) AliAODTrack(*t);
      }
    }
    if (1) { 
      AliAODCaloTrigger *out = eout->GetCaloTrigger("EMCAL");
      AliAODCaloTrigger *in  = evin->GetCaloTrigger("EMCAL");
      *out = *in;
    }
    if (1) { 
      AliAODCaloCells *out = eout->GetEMCALCells();                  
      AliAODCaloCells *in  = evin->GetEMCALCells();    
      *out = *in;
    }
    if (1) { 
      TClonesArray *out = eout->GetCaloClusters();	         
      TClonesArray *in  = evin->GetCaloClusters();  
      for (Int_t i=0;i<in->GetEntriesFast();++i) {
	AliAODCaloCluster *c = static_cast<AliAODCaloCluster*>(in->At(i));
	new ((*out)[i]) AliAODCaloCluster(*c);
      }
    }

    if (1) {
      TClonesArray *out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
      TClonesArray *in  = static_cast<TClonesArray*>(evin->FindListObject(AliAODMCParticle::StdBranchName()));
      if (in && !out) {
	fgAODMCParticles = new TClonesArray("AliAODMCParticle",1000);
	fgAODMCParticles->SetName(AliAODMCParticle::StdBranchName());
	oh->AddBranch("TClonesArray", &fgAODMCParticles);
	out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
      } 
      if (in && out) {
	if (fCutMC) {
	  for (Int_t i=0;i<in->GetEntriesFast();++i) {
	    AliAODMCParticle *mc = static_cast<AliAODMCParticle*>(in->At(i));
	    if (TMath::Abs(mc->Y())>1.2)
	      new ((*out)[i]) AliAODMCParticle;
	    else
	      new ((*out)[i]) AliAODMCParticle(*mc);
	  }
	}
      }      
    }

    if (1) {
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
    fTrials = 0;
    PostData(1, fOutputList);
  }
}

Bool_t AliAodSkimTask::UserNotify()
{
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
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
    cout << "AliAodSkimTask " << GetName() << " found " << xsection << " " << trials << " " << pthardbin << " " << nevents << endl;
    fPyxsec      = xsection;      
    fPytrials    = trials;      
    fPypthardbin = pthardbin;
  }

  return res;
}

void AliAodSkimTask::Terminate(Option_t *)
{
  cout << "AliAodSkimTask " << GetName() << " terminated with accepted fraction of events: " << fHevs->GetBinContent(2)/fHevs->GetEntries() << endl;
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
