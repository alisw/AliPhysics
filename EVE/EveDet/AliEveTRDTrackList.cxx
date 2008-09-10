// Uncomment to display debugging infos
//#define ALIEVETRDTRACKLIST_DEBUG

#include "AliEveTRDTrackList.h"

#include <AliTRDtrackV1.h>
#include <TFile.h>
#include <TFunction.h>
#include <TList.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeStream.h>
#include <EveDet/AliEveTRDData.h>

ClassImp(AliEveTRDTrackList)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackList ////////////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackList::AliEveTRDTrackList(const Text_t* n, const Text_t* t, Bool_t doColor):
  TEveElementList(n, t, doColor),
  fMacroList(0),
  fMacroSelList(0),
  fDataFromMacroList(0),
  fDataTree(0)
{
  // Only accept childs of type AliEveTRDTrack
  SetChildClass(AliEveTRDTrack::Class());

  fMacroList = new TList();
  fMacroSelList = new TList();
  fDataFromMacroList = new TList();

  fDataTree = new TTreeSRedirector("TRD.TrackListMacroData.root");

  AddStandardMacros();
}

//__________________________________________________________
AliEveTRDTrackList::~AliEveTRDTrackList()
{
  if (fMacroList != 0)
  {
    fMacroList->Clear();
    delete fMacroList;
    fMacroList = 0;
  }
  if (fMacroSelList != 0)
  {
    fMacroSelList->Clear();
    delete fMacroSelList;
    fMacroSelList = 0;
  } 
  if (fDataFromMacroList != 0)
  {
    fDataFromMacroList->Clear();
    delete fDataFromMacroList;
    fDataFromMacroList = 0;
  } 
  if (fDataTree != 0)
  {
    delete fDataTree;
    fDataTree = 0;
  } 
}

//__________________________________________________________
Int_t AliEveTRDTrackList::AddMacro(const Char_t* path, const Char_t* nameC)
{
  // First check the type of the macro:
  // If it has the signature of a selection macro:
  // Bool_t MacroName(AliTRDtrackV1*)
  // it is assumed to be a selection macro.
  // If it has the signature of a process macro:
  // void MacroName(AliTRDtrackV1*, Double_t*&, Int_t&)
  // it is assumed to be a process macro.
  // In all other cases: Macro is rejected
  Bool_t isSelectionMacro = kFALSE;
  Bool_t hasCorrectSignature = kFALSE;
  

  Char_t* entryName = MakeMacroEntry(path, nameC);

  Char_t* pathname = new Char_t[300];
  memset(pathname, '\0', sizeof(Char_t) * 300);

  // Expand the path and create the pathname
  Char_t* systemPath = gSystem->ExpandPathName(path);
  sprintf(pathname, "%s/%s", systemPath, nameC);
  delete systemPath;
  systemPath = 0;

  // Delete ".C" from filename
  Char_t* name = new Char_t[strlen(nameC)];
  memset(name, '\0', sizeof(Char_t) * strlen(nameC));
  strncpy(name, nameC, strlen(nameC) - 2);

  // Check, if files exists
  FILE* fp = 0;

  fp = fopen(pathname, "rb");
  if (fp != 0)
  {
    fclose(fp);
    fp = 0;
  }
  else
  {
    if (name != 0)  delete name;
    name = 0;
    if (pathname != 0)  delete pathname;
    pathname = 0;
    if (entryName != 0)  delete entryName;
    entryName = 0;

    return NOT_EXIST_ERROR;
  }

  // Clean up root, load the desired macro and then check the type of the macro
  //gROOT->Reset("a");
  gROOT->Reset();
  gROOT->ProcessLineSync(Form(".L %s", pathname));

  // Selection macro?
  TFunction* f = gROOT->GetGlobalFunctionWithPrototype(name, "AliTRDtrackV1*", kTRUE);
  if (f != 0x0)
  {
    if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
    {
      hasCorrectSignature = kTRUE;
      isSelectionMacro = kTRUE;
    }
  }
  // Process macro?
  else
  {
    f = gROOT->GetGlobalFunctionWithPrototype(name, "AliTRDtrackV1*, Double_t*&, Int_t&", kTRUE);
    if (f != 0x0)
    {
      if (!strcmp(f->GetReturnTypeName(), "void"))
      {
        hasCorrectSignature = kTRUE;
        isSelectionMacro = kFALSE;
      }
    }
  }

  // Clean up again / unload this function
  gROOT->ProcessLineSync(Form(".U %s", pathname));
  //gROOT->Reset("a");
  gROOT->Reset();
  if (name != 0)  delete name;
  name = 0;
  if (pathname != 0)  delete pathname;
  pathname = 0;

  // Has not the correct signature!
  if (!hasCorrectSignature) 
  {
    if (entryName != 0)  delete entryName;
    entryName = 0;
    return SIGNATURE_ERROR;
  }

  Int_t returnValue = WARNING;

  // Only add macro, if it is not already in the list
  if (!isSelectionMacro && fMacroList->FindObject(entryName) == 0)
  {
    fMacroList->Add(new TObjString(entryName));
    fMacroList->Sort();

    returnValue = SUCCESS;
  }
  else if (isSelectionMacro && fMacroSelList->FindObject(entryName) == 0)
  {
    fMacroSelList->Add(new TObjString(entryName));
    fMacroSelList->Sort();
    
    returnValue = SUCCESS;
  }
  else  returnValue = WARNING;

  if (entryName != 0)  delete entryName;
  entryName = 0;

  return returnValue;
}

//__________________________________________________________
void AliEveTRDTrackList::AddMacroFast(const Char_t* path, const Char_t* name, Bool_t toSelectionList)
{
  Char_t* entry = MakeMacroEntry(path, name);
  if (entry != 0) {
    if (toSelectionList)  fMacroSelList->Add(new TObjString(entry));
    else                  fMacroList->Add(new TObjString(entry));
    
    delete entry;
    entry = 0;

#ifdef ALIEVETRDTRACKLIST_DEBUG
    // Successfull add will only be displayed in debug mode
    printf("#AliEveTRDTrackList: Standard macros: Added macro %s/%s to %s list\n", path, name, 
           (toSelectionList ? "selection" : "process"));
#endif
  } else {
    // Error will always be displayed
    printf("#AliEveTRDTrackList: Standard macros: ERROR: Could not add macro %s/%s to %s list\n", path, name, 
           (toSelectionList ? "selection" : "process"));
  } 
}

//__________________________________________________________
void AliEveTRDTrackList::AddStandardMacros()
{
  // Add your standard macros here, e.g.: 
  // To add a macro without any checks (very fast, but unsafe):
  // AddMacroFast("$(ALICE_ROOT)/myFolder", "myMacroName.C", isSelMacro);
  // To add a macro with checks (slower, but safe):
  // AddMacro("$(ALICE_ROOT)/myFolder", "myMacroName.C");
  // -> If the file does not exist, nothing happens. So if you want to handle this,
  // use the return value of AddMacro (NOT_EXIST_ERROR is returned, if file does not exist)
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "clusterSelection.C");
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "chargeDistr.C");
}

//__________________________________________________________
void AliEveTRDTrackList::ApplyProcessMacros(TList* iterator)
{
  Char_t name[100];
  Char_t path[300];
  Char_t pathname[400];
  Char_t cmd[430];

  AliEveTRDTrack* track = 0;
  AliTRDtrackV1 *trackv1 = 0;

  // Clear root
  gROOT->Reset();
  
  // Clear old data and re-allocate
  if (fDataFromMacroList != 0) delete fDataFromMacroList;
  fDataFromMacroList = new TList();

  if (fDataTree == 0) fDataTree = new TTreeSRedirector("TRD.TrackListMacroData.root");
  
  // Walk through the list of tracks
  Int_t trackNum = 0;
  
  //Double_t* results = new Double_t[iterator->GetEntries()];
  //for (Int_t i = 0; i < iterator->GetEntries(); i++)  results[i] = 0;
  //Double_t result = 0;
    
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter, ++trackNum)
  {
    track = dynamic_cast<AliEveTRDTrack*>(*iter);

    if (!track)  continue;
      
    // Skip tracks that have not been selected
    if (!track->GetRnrState())  continue;
      
    trackv1 = (AliTRDtrackV1*)track->GetUserData();

    track->ExportToCINT((Text_t*)"automaticTrack");
    // Cast to AliTRDtrackV1
    gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1 = (AliTRDtrackV1*)automaticTrack->GetUserData();");

    // Collect data for each macro
    for (Int_t i = 0; i < iterator->GetEntries(); i++)
    {
      memset(name, '\0', sizeof(Char_t) * 100);
      memset(path, '\0', sizeof(Char_t) * 300);
      memset(pathname, '\0', sizeof(Char_t) * 400);
      memset(cmd, '\0', sizeof(Char_t) * 430);

#ifdef ALIEVETRDTRACKLIST_DEBUG
      // Display this message only once
      if (iter == this->BeginChildren()) 
        printf("AliEveTRDTrackList: Applying process macro: %s\n", iterator->At(i)->GetTitle());
#endif

      // Extract path and name -> Make pathname
      sscanf(iterator->At(i)->GetTitle(), "%s (Path: %s)", name, path);
      // Delete ")" at the end of path
      path[strlen(path)] = '\0';
      path[strlen(path) - 1] = '\0';
      sprintf(pathname, "%s/%s", path, name);
      sprintf(cmd, ".x %s(automaticTrackV1, results, n);", pathname);

      // Add to "data-from-list", but only once!
      if (iter == this->BeginChildren())  fDataFromMacroList->Add(new TObjString(name));

      //results[i] += (Double_t)gROOT->ProcessLineSync(cmd); 
      //result = (Double_t)gROOT->ProcessLineSync(cmd);

      //(*fDataTree) << Form("TrackData%d", i)
      //  << Form("Macro%d=", i) << result << "\n";   
      
      // Create data pointers in CINT, execute the macro and get the data
      gROOT->ProcessLineSync("Double_t* results = 0;");
      gROOT->ProcessLineSync("Int_t n = 0;");
      gROOT->ProcessLineSync(cmd);
      Double_t* results = (Double_t*)gROOT->ProcessLineSync("results;");
      Int_t nResults = (Int_t)gROOT->ProcessLineSync("n;");
      
      if (results == 0)
      {
        Error("Apply macros", Form("Error reading data from macro \"%s\"", name));
        continue;
      }
      for (Int_t resInd = 0; resInd < nResults; resInd++)
      {
        (*fDataTree) << Form("TrackData%d", i) << Form("Macro%d=", i) << results[resInd] << (Char_t*)"\n";   
      }

      delete results;
      results = 0;
    }
  }    
/*  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    (*fDataTree) << Form("TrackData") << Form("Macro%d=", i) << results[i] << "\n";    
  }

  if (results != 0) delete results;
  results = 0;
*/

  delete fDataTree;
  fDataTree = 0;

  // Clear root
  gROOT->Reset();
  
  // Now the data is stored in "TRD.TrackListMacroData.root"
  // The editor will access this file to display the data
}


//__________________________________________________________
void AliEveTRDTrackList::ApplySelectionMacros(TList* iterator)
{
  Char_t name[100];
  Char_t path[300];
  Char_t pathname[400];
  Char_t cmd[430];

  AliEveTRDTrack* track = 0;
  AliTRDtrackV1 *trackv1 = 0;
  Bool_t selectedByMacro = kFALSE;

  // Clear root
  gROOT->Reset();

  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    memset(name, '\0', sizeof(Char_t) * 100);
    memset(path, '\0', sizeof(Char_t) * 300);
    memset(pathname, '\0', sizeof(Char_t) * 400);
    memset(cmd, '\0', sizeof(Char_t) * 430);

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Applying selection macro: %s\n", iterator->At(i)->GetTitle());
#endif
      
    // Extract path and name -> Make pathname
    sscanf(iterator->At(i)->GetTitle(), "%s (Path: %s)", name, path);
    // Delete ")" at the end of path.
    path[strlen(path)] = '\0';
    path[strlen(path) - 1] = '\0';
    
    sprintf(pathname, "%s/%s", path, name);
    sprintf(cmd, ".x %s(automaticTrackV1);", pathname);

    // Walk through the list of tracks
    for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter)
    {
      track = dynamic_cast<AliEveTRDTrack*>(*iter);

      if (!track) continue;
      
      trackv1 = (AliTRDtrackV1*)track->GetUserData();

      track->ExportToCINT((Text_t*)"automaticTrack");
      // Cast to AliTRDtrackV1
      gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1 = (AliTRDtrackV1*)automaticTrack->GetUserData();");
      selectedByMacro = (Bool_t)gROOT->ProcessLineSync(cmd);
      track->SetRnrState(selectedByMacro);         
    }
  }

  // Clear root
  gROOT->Reset();  
}


//__________________________________________________________
Char_t* AliEveTRDTrackList::MakeMacroEntry(const Char_t* path, const Char_t* name)
{
  Char_t* entry = new Char_t[400];
  memset(entry, '\0', sizeof(Char_t) * 400);

  Char_t* systemPath = gSystem->ExpandPathName(path);
  sprintf(entry, "%s (Path: %s)", name, systemPath);
  delete systemPath;
  systemPath = 0;

  return entry;
}


//__________________________________________________________
void AliEveTRDTrackList::RemoveProcessMacros(TList* iterator) 
{
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    fMacroList->Remove(fMacroList->FindObject(iterator->At(i)->GetTitle()));
  }
}


//__________________________________________________________
void AliEveTRDTrackList::RemoveSelectionMacros(TList* iterator) 
{
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    fMacroSelList->Remove(fMacroSelList->FindObject(iterator->At(i)->GetTitle()));
  }
}
