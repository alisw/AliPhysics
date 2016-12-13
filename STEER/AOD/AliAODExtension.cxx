#include "AliAODExtension.h"

//-------------------------------------------------------------------------
//     Support class for AOD extensions. This is created by the user analysis
//     that requires a separate file for some AOD branches. The name of the 
//     AliAODExtension object is the file name where the AOD branches will be
//     stored.
//-------------------------------------------------------------------------

#include "AliAODBranchReplicator.h"
#include "AliAODEvent.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "Riostream.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TList.h"
#include "TMap.h"
#include "TMap.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

using std::endl;
using std::cout;
ClassImp(AliAODExtension)

//______________________________________________________________________________
AliAODExtension::AliAODExtension() : TNamed(), 
fAODEvent(0), fTreeE(0), fFileE(0), fNtotal(0), fNpassed(0), 
fSelected(kFALSE), fTreeBuffSize(30000000), fMemCountAOD(0),
fRepFiMap(0x0), fRepFiList(0x0), fEnableReferences(kTRUE), fObjectList(0)
{
  // default ctor
}

//______________________________________________________________________________
AliAODExtension::AliAODExtension(const char* name, const char* title, Bool_t isfilter)
:TNamed(name,title), 
fAODEvent(0), 
fTreeE(0), 
fFileE(0), 
fNtotal(0), 
fNpassed(0),
fSelected(kFALSE),
fTreeBuffSize(30000000),
fMemCountAOD(0),
fRepFiMap(0x0),
fRepFiList(0x0),
fEnableReferences(kTRUE),
fObjectList(0x0)
{
  // Constructor.
  if (isfilter) {
    TObject::SetBit(kFilteredAOD);
    printf("####### Added AOD filter %s\n", name);
  } else printf("####### Added AOD extension %s\n", name);
  KeepUnspecifiedBranches();
}   

//______________________________________________________________________________
AliAODExtension::~AliAODExtension()
{
  // Destructor.
  if(fFileE){
    // is already handled in TerminateIO
    fFileE->Close();
    delete fFileE;
    fTreeE = 0;
    fAODEvent = 0;
  }
  if (fTreeE) delete fTreeE;
  if (fRepFiMap) fRepFiMap->DeleteAll();
  delete fRepFiMap; // the map is owner
  delete fRepFiList; // the list is not
  delete fObjectList; // not owner
}

//______________________________________________________________________________
void AliAODExtension::AddBranch(const char* cname, void* addobj)
{
  // Add a new branch to the aod 
  
  if (!fAODEvent) {
    char type[20];
    gROOT->ProcessLine(Form("{TString s_tmp; AliAnalysisManager::GetAnalysisManager()->GetAnalysisTypeString(s_tmp); sprintf((char*)%p, \"%%s\", s_tmp.Data());}", type));
    Init(type);
  }
  TDirectory *owd = gDirectory;
  if (fFileE) {
    fFileE->cd();
  }
  char** apointer = (char**) addobj;
  TObject* obj = (TObject*) *apointer;
  
  fAODEvent->AddObject(obj);
  
  TString bname(obj->GetName());
  
  if (!fTreeE->FindBranch(bname.Data())) 
  {
    Bool_t acceptAdd(kTRUE);
    
    if ( TestBit(kDropUnspecifiedBranches) )
    {
      // check that this branch is in our list of specified ones...
      // otherwise do not add it !
      TIter next(fRepFiMap);
      TObjString* p;
      
      acceptAdd=kFALSE;
      
      while ( ( p = static_cast<TObjString*>(next()) ) && !acceptAdd )
      {
        if ( p->String() == bname ) acceptAdd=kTRUE;
      }
    }
    
    if ( acceptAdd ) 
    {
      // Do the same as if we book via 
      // TTree::Branch(TCollection*)
      
      fObjectList->Add(obj);

      const Int_t kSplitlevel = 99; // default value in TTree::Branch()
      const Int_t kBufsize = 32000; // default value in TTree::Branch()
      
      fTreeE->Bronch(bname.Data(), cname, 
                     fAODEvent->GetList()->GetObjectRef(obj),
                     kBufsize, kSplitlevel - 1);
    }
  }
  owd->cd();
}

//______________________________________________________________________________
Bool_t AliAODExtension::FinishEvent()
{
  // Fill current event.
  fNtotal++;
  if (!IsFilteredAOD()) {
    fAODEvent->MakeEntriesReferencable();
    FillTree();
    return kTRUE;
  }  
  // Filtered AOD. Fill only if event is selected.
  if (!fSelected) return kTRUE;
  
  TIter next(fRepFiList);
  
  AliAODBranchReplicator* repfi;
  
  while ( ( repfi = static_cast<AliAODBranchReplicator*>(next()) ) )
  {
    repfi->ReplicateAndFilter(*fAODEvent);
  }
  fNpassed++;
  FillTree();
  fSelected = kFALSE; // so that next event will not be selected unless demanded
  return kTRUE;
}  

//______________________________________________________________________________
void AliAODExtension::FillTree() 
{
  //
  //   Fill AOD extension tree and check AutoFlush settings
  //
  
  Long64_t nbf = fTreeE->Fill();
  
  // Check buffer size and set autoflush if fTreeBuffSize is reached
  if (fTreeBuffSize>0 && fTreeE->GetAutoFlush()<0 && 
      (fMemCountAOD += nbf)>fTreeBuffSize ) { // default limit is still not reached
    nbf = fTreeE->GetZipBytes();
    if (nbf>0) nbf = -nbf;
    else       nbf = fTreeE->GetEntries();
    fTreeE->SetAutoFlush(nbf);
    AliInfo(Form("Calling fTreeE->SetAutoFlush(%lld) | W:%lld T:%lld Z:%lld", 
		 nbf,fMemCountAOD,fTreeE->GetTotBytes(),fTreeE->GetZipBytes()));  
    
  }
}

//______________________________________________________________________________
Bool_t AliAODExtension::Init(Option_t *option)
{
  // Initialize IO.
  
  AliCodeTimerAuto(GetName(),0);
  
  if(!fAODEvent) 
  {
    fAODEvent = new AliAODEvent();    
  }
  
  TDirectory *owd = gDirectory;
  TString opt(option);
  opt.ToLower();
  
  if (opt.Contains("proof")) 
  {
    // proof
    // Merging via files. Need to access analysis manager via interpreter.
    gROOT->ProcessLine(Form("AliAnalysisDataContainer *c_common_out = AliAnalysisManager::GetAnalysisManager()->GetCommonOutputContainer();"));
    gROOT->ProcessLine(Form("AliAnalysisManager::GetAnalysisManager()->OpenProofFile(c_common_out, \"RECREATE\", \"%s\");", fName.Data()));
    fFileE = gFile;
  } 
  else 
  {
    fFileE = new TFile(GetName(), "RECREATE");
  }  
  fTreeE = new TTree("aodTree", "AliAOD tree");
  
  delete fObjectList;
  fObjectList = new TList;
  fObjectList->SetOwner(kFALSE); // be explicit we're not the owner...
  TList* inputList = fAODEvent->GetList();
  TIter next(inputList);
  TObject* o;
  
  while ( ( o = next() ) )
  {
    // Loop on the objects that are within the main AOD, and see what to do with them :
    // - transmit them to our AOD as they are
    // - filter them (by means of an AliAODBranchReplicator)
    // - drop them completely
    
    Bool_t mustKeep(kFALSE);
    
    TString test(o->ClassName());
    test.ToUpper();
    // check if there is a replicator for the header
    Bool_t headerHasReplicator = fRepFiMap && (fRepFiMap->FindObject(o->GetName())!=0x0);
    if (test.BeginsWith("ALIAODHEADER") && !headerHasReplicator)
    {
      // do not allow to drop header branch
      mustKeep=kTRUE;
    }
    
    if ( fRepFiMap && !mustKeep )
    {
      // we have some replicators, so let's see what the relevant one decides about this object
      TObject* specified = fRepFiMap->FindObject(o->GetName()); // FindObject finds key=o->GetName() in the map
      if (specified)
      {
        AliAODBranchReplicator* repfi = dynamic_cast<AliAODBranchReplicator*>(fRepFiMap->GetValue(o->GetName())); // GetValue gets the replicator corresponding to key=o->GetName()
        if ( repfi ) 
        {        
          TList* replicatedList = repfi->GetList();
          if (replicatedList)
          {
            AliAODEvent::AssignIDtoCollection(replicatedList);
            TIter nextRep(replicatedList);
            TObject* objRep;
            while ( ( objRep = nextRep() ) )
            {
              if ( !fObjectList->FindObject(objRep) ) // insure we're not adding several times the same object
              {                
                fObjectList->Add(objRep);                  
              }
            }
          }
          else
          {
            AliError(Form("replicatedList from %s is null !",repfi->GetName()));
          }
        }
      }
      else
      {
        if ( !TestBit(kDropUnspecifiedBranches) ) 
        {
          // object o will be transmitted to the output AOD, unchanged
          fObjectList->Add(o);
        }
      }
    } 
    else
    {
      // no replicator, so decide based on the policy about dropping unspecified branches
      if ( mustKeep || !TestBit(kDropUnspecifiedBranches) )
      {
        // object o will be transmitted to the output AOD, unchanged
        fObjectList->Add(o);
      }
    }
  }
    
  if (fEnableReferences) 
  {
    fTreeE->BranchRef();    
  }
    
  fTreeE->Branch(fObjectList);
  
  owd->cd();
  
  return kTRUE;
}

//______________________________________________________________________________
void AliAODExtension::Print(Option_t* opt) const
{
  // Print info about this extension
  
  cout << opt << Form("%s - %s - %s - aod %p",IsFilteredAOD() ? "FilteredAOD" : "Extension",
                      GetName(),GetTitle(),GetAOD()) << endl;
  if ( !fEnableReferences ) 
  {
    cout << opt << opt << "References are disabled ! Hope you know what you are doing !" << endl;
  }
  if ( TestBit(kDropUnspecifiedBranches) )
  {
    cout << opt << opt << "All branches not explicitely specified will be dropped" << endl;
  }
  
  TIter next(fRepFiMap);
  TObjString* s;
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    AliAODBranchReplicator* br = static_cast<AliAODBranchReplicator*>(fRepFiMap->GetValue(s->String().Data()));
    
    cout << opt << opt << "Branch " << s->String();
    if (br)
    {
      cout << " will be filtered by class " << br->ClassName();
    }
    else
    {
      cout << " will be transmitted as is";
    }
    cout << endl;
  }
}

//______________________________________________________________________________
void AliAODExtension::SetEvent(AliAODEvent* event)
{
  // Connects to an external event
  if (!IsFilteredAOD()) {
    Error("SetEvent", "Not allowed to set external event for non filtered AOD's");   
    return;
  }
  fAODEvent = event;
}

//______________________________________________________________________________
void AliAODExtension::AddAODtoTreeUserInfo()
{
  // Add aod event to tree user info
  
  if (!fTreeE) return;
  
  AliAODEvent* aodEvent(fAODEvent);
  
  if ( IsFilteredAOD() )
  {
    // cannot attach fAODEvent (which is shared with our AliAODHandler mother)
    // so we create a custom (empty) AliAODEvent 
    // Has also the advantage we can specify only the list of objects
    // that are actually there in this filtered aod
    //
    aodEvent = new AliAODEvent;
    TIter nextObj(fObjectList);
    TObject* o;
    while ( ( o = nextObj() ) ) 
    {
      aodEvent->AddObject(o);
    }    
  }

  TList *l = aodEvent->GetList();
  if (l) {
    for(int i = 0;i < l->GetEntries(); ++i){
      TObject *pObject = l->At(i);
      if(pObject->InheritsFrom(TClonesArray::Class())){
       ((TClonesArray*)pObject)->Delete();
      } else if(!pObject->InheritsFrom(TCollection::Class())){
       TClass *pClass = TClass::GetClass(pObject->ClassName());
       if (pClass && pClass->GetListOfMethods()->FindObject("Clear")) {
         AliDebug(1, Form("Clear for object %s class %s", pObject->GetName(), pObject->ClassName()));
         pObject->Clear();
       }
      } else {
         AliWarning(Form("No method to clear for object %s class %s", pObject->GetName(), pObject->ClassName()));
      }
    }
  }
  fTreeE->GetUserInfo()->Add(aodEvent);
}

//______________________________________________________________________________
Bool_t AliAODExtension::TerminateIO()
{
  // Terminate IO
  if (TObject::TestBit(kFilteredAOD))
    printf("AOD Filter %s: events processed: %d   passed: %d\n", GetName(), fNtotal, fNpassed);
  else
    printf("AOD extension %s: events processed: %d\n", GetName(), fNtotal);
  if (fFileE) 
  {
    fFileE->Write();
    fFileE->Close();
    delete fFileE;
    fFileE = 0;
    fTreeE = 0;
    fAODEvent = 0;
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliAODExtension::FilterBranch(const char* branchName, AliAODBranchReplicator* repfi)
{
  // Specify a filter/replicator for a given branch
  //
  // If repfi=0x0, this will disable the branch (in the output) completely.
  //
  // repfi is adopted by this class, i.e. user should not delete it.
  //
  // WARNING : branch name must be exact.
  //
  // See also the documentation for AliAODBranchReplicator class.
  //
  
  if (!fRepFiMap)
  {
    fRepFiMap = new TMap;
    fRepFiMap->SetOwnerKeyValue(kTRUE,kTRUE);
    fRepFiList = new TList;
    fRepFiList->SetOwner(kFALSE);
  }
  
  fRepFiMap->Add(new TObjString(branchName),repfi);
  
  if (repfi && !fRepFiList->FindObject(repfi))
  {
    // insure we get unique and non-null replicators in this list
    fRepFiList->Add(repfi);
  }
}

