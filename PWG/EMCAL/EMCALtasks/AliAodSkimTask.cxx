#include <Riostream.h>
#include <TChain.h>
#include <TGraph.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <AliAODMCHeader.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODMCParticle.h>
#include <AliAnalysisManager.h>
#include "AliAodSkimTask.h"


using namespace std;
ClassImp(AliAodSkimTask)

AliAodSkimTask::AliAodSkimTask() : AliAnalysisTaskSE(), fClusMinE(-1), fCutMC(1), fTrials(0), fAOD(0), fAODMcHeader(0), fOutputList(0)
{
}

AliAodSkimTask::AliAodSkimTask(const char* name) : AliAnalysisTaskSE(name), fClusMinE(-1), fCutMC(1), fTrials(0), fAOD(0), fAODMcHeader(0), fOutputList(0)
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
      if (clus->E()>fClusMinE) {
	store = kTRUE;
	break;
      }
    }
  } else { 
    store = kTRUE;
  }

  if (!store) {
    ++fTrials;
    return;
  }

  //AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
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
      new (out) TClonesArray(*in); 
    }
    if (1) {   
      AliTOFHeader *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader()); 
      const AliTOFHeader *in = evin->GetTOFHeader();	    
      *out = *in;                  
    }
    if (1) {
      TClonesArray *out = eout->GetTracks();	                 
      TClonesArray *in  = evin->GetTracks();	
      new (out) TClonesArray(*in);     
    }
    if (1) { 
      AliAODCaloTrigger *out = eout->GetCaloTrigger("EMCAL");
      AliAODCaloTrigger *in  = evin->GetCaloTrigger("EMCAL");
      *out = *in;
    }
    if (1) { 
      AliAODCaloCells *out = eout->GetEMCALCells();                  
      AliAODCaloCells *in  = evin->GetEMCALCells();    
      new (out) AliAODCaloCells(*in);  
    }
    if (1) { 
      TClonesArray *out = eout->GetCaloClusters();	         
      TClonesArray *in  = evin->GetCaloClusters();  
      new (out) TClonesArray(*in); 
    }

    if (1) {
      TClonesArray *out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
      TClonesArray *in  = static_cast<TClonesArray*>(evin->FindListObject(AliAODMCParticle::StdBranchName()));
      if (in && !out) {
	fgAODMCParticles = new TClonesArray("AliAODMCParticle",1000);
	fgAODMCParticles->SetName(AliAODMCParticle::StdBranchName());
	oh->AddBranch("TClonesArray", &fgAODMCParticles);
	out = static_cast<TClonesArray*>(eout->FindListObject(AliAODMCParticle::StdBranchName()));
      } else if (in && out) {
	new (out) TClonesArray(*in); 
	if (fCutMC) {
	  for (Int_t i=0;i<out->GetEntriesFast();++i) {
	    AliAODMCParticle *mc = static_cast<AliAODMCParticle*>(out->At(i));
	    if (TMath::Abs(mc->Y())>1.2)
	      new ((*out)[i]) AliAODMCParticle;
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
      } else if (in && out) {
	*out = *in;
      }
    }

    fTrials = 0;
    PostData(1, fOutputList);
  }
}

void AliAodSkimTask::Terminate(Option_t *)
{
}
