#include "AliFMDMCHitHandler.h"
#include "AliAnalysisManager.h"
#include <TError.h>
#include <AliLog.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TROOT.h>
#include <AliStack.h>
#include <AliMCEvent.h>

ClassImp(AliFMDMCHitHandler)
#if 0 // For Emacs - do not remove
;
#endif

//____________________________________________________________________
AliFMDMCHitHandler::AliFMDMCHitHandler(const char* name,
				       const char* what,
				       AliMCEventHandler* parent)
  : AliMCEventHandler(name, what),
    fParent(parent), 
    fFile(0),
    fTree(0),
    fDir(0),
    fArray(0),
    fNEvents(0), 
    fNEventsPerFile(0),
    fNEventsInContainer(0),
    fEvent(0), 
    fFileNumber(0),
    fTreeName(""),
    fFileBase("")
{
  // Constructor 
  // 
  // Parameters: 
  //   name The name 
}

//____________________________________________________________________
TString*
AliFMDMCHitHandler::GetParentPath() const
{
  if (!fParent) { 
    AliWarning("No parent");
    return 0;
  }
  return fParent->GetInputPath();
}

//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::Init(Option_t* opt)
{
  // Initialize the input
  // 
  // @param opt Options 
  // 
  // @return true on success 
  AliDebugF(10,"AliFMDMCHitHandler::Init(\"%s\")", opt);

  TString option(opt);
  if (option.EqualTo("proof") || option.EqualTo("local")) return true;

  TString t = "Tree";
  TString b = "";
  TClass* cl = gROOT->GetClass(GetTitle());
  if (cl) { 
    if (cl->InheritsFrom("AliHit")) { 
      t += "H";
      b =  "Hits";
    }
    else if (cl->InheritsFrom("AliSDigit")) {
      t += "S";
      b =  "SDigits";
    }
    else if (cl->InheritsFrom("AliDigit")) {
      t += "D";
      b =  "Digits";
    }
    else 
      t = "";
  }
  if (!t.IsNull()) fTreeName = t;
  if (!b.IsNull()) fFileBase = b;


  fArray = new TClonesArray(GetTitle());
  
  TTree* treeE = fParent->GetTree();
  if (!treeE) { 
    AliError("Parent does not have an events tree");
    return false;
  }

  // Get number of events in this directory 
  fNEventsPerFile = -1;
  fNEvents        = treeE->GetEntries();
  fEvent          = 0;
  fFileNumber     = 0;

  if (!OpenFile(fFileNumber)) return false;

  return true;
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::BeginEvent(Long64_t entry)
{
  // Called at the beginning of an event 
  // 
  // @param entry Entry in tree 
  // 
  // @return true on success
  AliDebugF(10,"AliFMDMCHitHandler::BeginEvent(%lld)", entry);

  if (entry == -1) 
    fEvent++;
  else 
    fEvent = entry;

  if (fEvent >= fNEvents) { 
    AliWarningF("Event number out of range %d/%d", fEvent, fNEvents);
    return false;
  }

  if (fNEventsPerFile < 0) {
    TTree* treeK = fParent->TreeK();
    if (!treeK) { 
      AliError("Parent does not have a kinematics tree");
      return false;
    }
  
    TFile* fileK = treeK->GetCurrentFile();
    if (!fileK) { 
      AliError("Kinematics tree has no associated file");
      return false;
    }
    // Get the number of events per file 
    fNEventsPerFile = fileK->GetNkeys() - fileK->GetNProcessIDs();
  }
  return LoadEvent(fEvent);
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::Notify(const char* path)
{
  // Called when the input file is changed 
  // 
  // @param path New path 
  //
  // @return true on success
  AliDebugF(10,"AliFMDMCHitHandler::Notify(\"%s\")", path);
  return true;
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::FinishEvent()
{
  // Called at the end of an event 
  // 
  // @return true on success
  AliDebug(10,"AliFMDMCHitHandler::FinishEvent()");
  return true;
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::Terminate()
{
  // Called at the end of a job 
  // 
  // @return true on success 
  AliDebug(10,"AliFMDMCHitHandler::Terminate()");
  return true;
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::TerminateIO()
{
  // Called at the end of a sub-job
  // 
  // @return true on success
  AliDebug(10,"AliFMDMCHitHandler::TerminateIO()");
  return true;
}

//____________________________________________________________________
void
AliFMDMCHitHandler::ResetIO()
{
  // Reset the I/O
  // 
  //
  AliDebug(10,"AliFMDMCHitHandler::ResetIO()");

  TString* path = GetParentPath();
  AliDebugF(10,"Got parent path %s", path ? path->Data() : "null");

  if (fFile) { 
    delete fFile;
    fFile = 0;
  }
}
//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::OpenFile(Int_t fileNo)
{
  TString* path = GetParentPath();
  AliDebugF(10,"Got parent path %s", path ? path->Data() : "null");
  if (!path) return false;

  TString ext("");
  if (fileNo > 0) ext = TString::Format("%d", fileNo);

  TString w(GetTitle());
  if (w.EndsWith("s")) w.Chop();
    
  TString fn = TString::Format("%s%s.%s%s.root", 
			       path->Data(), GetName(), 
			       fFileBase.Data(), ext.Data());
  Info("Init", "Opening %s", fn.Data());
  fFile = TFile::Open(fn, "READ");
  if (!fFile) { 
    AliErrorF("Failed to open %s", fn.Data());
    return false;
  }

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDMCHitHandler::LoadEvent(Int_t iev)
{
  // Load an event 
  // 
  // @param iev Event number 
  // 
  // @return true on success 
  AliDebugF(10,"AliFMDMCHitHandler::LoadEvent(%d)", iev);

  Int_t iNew = iev / fNEventsPerFile;
  if (iNew != fFileNumber) { 
    fFileNumber = iNew;
    if (!OpenFile(fFileNumber)) return false;
  }
  if (!fFile) return false;

  TString folder = TString::Format("Event%d", iev);
  fFile->GetObject(folder, fDir);
  if (!fDir) { 
    AliWarningF("Folder %s not found in file", folder.Data());
    return false;
  }

  fDir->GetObject(fTreeName, fTree);
  if (!fTree) { 
    AliWarningF("Folder %s does not contain the %s tree %s", 
		folder.Data(), GetTitle(), fTreeName.Data());
    return false;
  }

  fTree->SetBranchAddress(GetName(), &fArray);
  return true;
}


//____________________________________________________________________
AliFMDMCHitHandler*
AliFMDMCHitHandler::Create(const char* name, const char* what)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AliFMDMCHitHandler::Create", "No analysis manager");
    return 0;
  }
  
  AliVEventHandler* vmc = mgr->GetMCtruthEventHandler();
  if (!vmc) { 
    ::Error("AliFMDMCHitHandler::Create", "No MC truth handler");
    return 0;
  }
  
  AliMCEventHandler* mc = dynamic_cast<AliMCEventHandler*>(vmc);
  if (!mc) { 
    ::Error("AliFMDMCHitHandler::Create", 
	    "MC truth handler not a AliMCEventHandler, but %s", 
	    vmc->ClassName());
    return 0;
  }

  AliFMDMCHitHandler* ret = new AliFMDMCHitHandler(name, what, mc);
  mc->AddSubsidiaryHandler(ret);
  
  return ret;
}

//____________________________________________________________________
TClonesArray*
AliFMDMCHitHandler::GetParticleArray(AliFMDMCHitHandler* handler, 
				     Int_t particle)
{
  if (!handler) { 
    ::Error("AliFMDMCHitHandler::GetArray", "No handler passed");
    return 0;
  }

  AliMCEventHandler* mc = handler->GetParent();
  if (!mc) {
    ::Error("AliFMDMCHitHandler::GetArray", "Handler has no parent");
    return 0;
  }
  
  AliMCEvent* event = mc->MCEvent();
  if (!event) { 
    ::Error("AliFMDMCHitHandler::GetArray", "No MC event");
    return 0;
  }
  
  AliStack* stack = event->Stack();
  if (!event) { 
    ::Error("AliFMDMCHitHandler::GetArray", "Event has no stack");
    return 0;
  }
  
  TTree* tree = handler->GetTree();
  if (!tree) {
    ::Error("AliFMDMCHitHandler::GetArray", "Handler has no tree");
    return 0;
  }
    
  Int_t treeIdx = stack->TreeKEntry(particle);
  if (treeIdx < 0 || treeIdx >= tree->GetEntries()) { 
    ::Error("AliFMDMCHitHandler::GetArray", 
	    "Index %d of %d out of bounds [0,%lld]", treeIdx, particle, 
	    tree->GetEntries()-1);
    return 0;
  }
  
  tree->GetEntry(treeIdx);
  
  return handler->GetArray();
}

  
    
//____________________________________________________________________
//
// EOF
//
