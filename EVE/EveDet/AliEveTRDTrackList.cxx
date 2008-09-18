// Uncomment to display debugging infos
//#define ALIEVETRDTRACKLIST_DEBUG

#include "AliEveTRDTrackList.h"

#include <AliTRDtrackV1.h>
#include <TFile.h>
#include <TFunction.h>
#include <TH1.h>
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
  fDataTree(0),
  fHistoDataSelected(0),
  fMacroListSelected(0),
  fMacroSelListSelected(0),
  fSelectedTab(1)             // Standard tab: "Apply macros" (index 1)
{
  // Only accept childs of type AliEveTRDTrack
  SetChildClass(AliEveTRDTrack::Class());

  fMacroList = new TList();
  fMacroSelList = new TList();
  fDataFromMacroList = new TList();

  //fDataTree = new TTreeSRedirector("TRD.TrackListMacroData.root");

  AddStandardMacros();
}

//______________________________________________________
AliEveTRDTrackList::~AliEveTRDTrackList()
{
  if (fMacroList != 0)
  {
    fMacroList->Delete();
    delete fMacroList;
    fMacroList = 0;
  }
  if (fMacroSelList != 0)
  {
    fMacroSelList->Delete();
    delete fMacroSelList;
    fMacroSelList = 0;
  } 
  if (fDataFromMacroList != 0)
  {
    fDataFromMacroList->Delete();
    delete fDataFromMacroList;
    fDataFromMacroList = 0;
  } 
  if (fDataTree != 0)
  {
    delete fDataTree;
    fDataTree = 0;
  } 
}

//______________________________________________________
Int_t AliEveTRDTrackList::AddMacro(const Char_t* path, const Char_t* nameC, Bool_t forceReload)
{
  // First check the type of the macro:
  // If it has the signature of a selection macro:
  // Bool_t MacroName(AliTRDtrackV1*)
  // it is assumed to be a selection macro.
  // If it has the signature of a process macro:
  // void MacroName(AliTRDtrackV1*, Double_t*&, Int_t&)
  // it is assumed to be a process macro.
  // In all other cases: Macro is rejected
  Bool_t isHistoMacro = kFALSE;
  Bool_t isSelectionMacro = kFALSE;
  Bool_t hasCorrectSignature = kFALSE;
  

  Char_t* entryName = MakeMacroEntry(path, nameC);

  Char_t pathname[fkMaxMacroPathNameLength];
  memset(pathname, '\0', sizeof(Char_t) * fkMaxMacroPathNameLength);

  // Expand the path and create the pathname
  Char_t* systemPath = gSystem->ExpandPathName(path);
  sprintf(pathname, "%s/%s", systemPath, nameC);
  delete systemPath;
  systemPath = 0;

  // Delete ".C" from filename
  Char_t name[fkMaxMacroNameLength];
  memset(name, '\0', sizeof(Char_t) * fkMaxMacroNameLength);
  
  for (UInt_t ind = 0; ind < fkMaxMacroNameLength && ind < strlen(nameC) - 2; ind++)  name[ind] = nameC[ind];

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
    if (entryName != 0)  delete entryName;
    entryName = 0;

    return NOT_EXIST_ERROR;
  }

  // Clean up root, load the desired macro and then check the type of the macro
  //gROOT->Reset("a");
  gROOT->Reset();
 
  if (forceReload)  gROOT->ProcessLineSync(Form(".L %s++", pathname));
  else              gROOT->ProcessLineSync(Form(".L %s+", pathname));

  // Selection macro or process macro of type 2 (histo)?
  TFunction* f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*", kTRUE);
  if (f != 0x0)
  {
    // Some additional check (is the parameter EXACTLY of the desired type?)
    if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0)
    {
      // Selection macro?
      if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
      {      
        hasCorrectSignature = kTRUE;
        isSelectionMacro = kTRUE;
        isHistoMacro = kFALSE;
      }
      // Process macro of type 2 (histo)?
      else if (!strcmp(f->GetReturnTypeName(), "TH1*"))
      {
        hasCorrectSignature = kTRUE;
        isSelectionMacro = kFALSE;
        isHistoMacro = kTRUE;
      }
    }
  }
  // Process macro of type 1?
  else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*, Double_t*&, Int_t&", kTRUE)) != 0x0)
  {
    if (!strcmp(f->GetReturnTypeName(), "void"))
    {
      // Some additional check (are the parameters EXACTLY of the desired type?)
      if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0 &&
          strstr(f->GetMangledName(), "cODouble_tmUaNsP") != 0x0 &&
          strstr(f->GetMangledName(), "cOInt_taNsP") != 0x0)
      {
        hasCorrectSignature = kTRUE;
        isSelectionMacro = kFALSE;
        isHistoMacro = kFALSE;
      }
    }
  }

  //// Clean up again / unload this function
  //gROOT->ProcessLineSync(Form(".U %s", pathname));
  //gROOT->Reset("a");
  // Clean up again
  gROOT->Reset();
  
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

    // We do not know, where the element has been inserted - deselect this list
    fMacroListSelected = 0;

    returnValue = SUCCESS;
  }
  else if (isSelectionMacro && fMacroSelList->FindObject(entryName) == 0)
  {
    fMacroSelList->Add(new TObjString(entryName));
    fMacroSelList->Sort();
  
    // We do not know, where the element has been inserted - deselect this list
    fMacroSelListSelected = 0;
    
    returnValue = SUCCESS;
  }
  else  returnValue = WARNING;

  if (entryName != 0)  delete entryName;
  entryName = 0;

  return returnValue;
}

//______________________________________________________
void AliEveTRDTrackList::AddMacroFast(const Char_t* entry, Bool_t toSelectionList)
{
  if (toSelectionList)
  {
    fMacroSelList->Add(new TObjString(entry));
    fMacroSelList->Sort();

    // We do not know, where the element has been inserted - deselect this list
    fMacroSelListSelected = 0;
  }
  else 
  {
    fMacroList->Add(new TObjString(entry));
    fMacroList->Sort();

    // We do not know, where the element has been inserted - deselect this list
    fMacroListSelected = 0;
  }
}

//______________________________________________________
void AliEveTRDTrackList::AddMacroFast(const Char_t* path, const Char_t* name, Bool_t toSelectionList)
{
  Char_t* entry = MakeMacroEntry(path, name);
  if (entry != 0)
  {
    AddMacroFast(entry, toSelectionList);

#ifdef ALIEVETRDTRACKLIST_DEBUG
    // Successfull add will only be displayed in debug mode
    printf("#AliEveTRDTrackList::AddMacroFast: Added macro \"%s/%s\" to %s list\n", path, name, 
           (toSelectionList ? "selection" : "process"));
#endif
    
    delete entry;
    entry = 0;
  }
  else
  {
    // Error will always be displayed
    printf("#AliEveTRDTrackList::AddMacroFast: ERROR: Could not add macro \"%s/%s\" to %s list\n", path, name, 
           (toSelectionList ? "selection" : "process"));
  }
}

//______________________________________________________
void AliEveTRDTrackList::AddStandardMacros()
{
  // Add your standard macros here, e.g.: 
  // To add a macro without any checks (very fast, but unsafe):
  // AddMacroFast("$(ALICE_ROOT)/myFolder", "myMacroName.C", isSelMacro);
  // To add a macro with checks (slower, but safe):
  // AddMacro("$(ALICE_ROOT)/myFolder", "myMacroName.C");
  // -> If the file does not exist, nothing happens. So if you want to handle this,
  // use the return value of AddMacro (NOT_EXIST_ERROR is returned, if file does not exist)
  // (-> You can also check for other return values (see AddMacro(...)))
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "clusterSelection.C");
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "chargeDistr.C");
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "clusterResiduals.C");
  AddMacro("$(ALICE_ROOT)/TRD/qaRec/macros", "PH.C");
}

//______________________________________________________
void AliEveTRDTrackList::ApplyProcessMacros(TList* iterator)
{
  if (iterator->GetEntries() <= 0)  return;

  Char_t name[fkMaxMacroNameLength];
  Char_t** cmds = new Char_t*[iterator->GetEntries()];
  Bool_t* isHistoMacro = new Bool_t[iterator->GetEntries()];

  Int_t numHistoMacros = 0;
  TH1** histos = 0;

  AliEveTRDTrack* track = 0;
  AliTRDtrackV1 *trackv1 = 0;
  TH1* returnedHist = 0x0;

  // Clear root
  gROOT->Reset();
  
  // Clear old data and re-allocate
  if (fDataFromMacroList != 0) delete fDataFromMacroList;
  fDataFromMacroList = new TList();

  fHistoDataSelected = 0;

  if (fDataTree == 0) fDataTree = new TTreeSRedirector("TRD.TrackListMacroData.root");

  // Collect the commands for each macro and add them to "data-from-list"
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    memset(name, '\0', sizeof(Char_t) * fkMaxMacroNameLength);
    
    cmds[i] = new Char_t[(fkMaxMacroPathNameLength + fkMaxApplyCommandLength)];
    memset(cmds[i], '\0', sizeof(Char_t) * (fkMaxMacroNameLength + fkMaxApplyCommandLength));

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Applying process macro: %s\n", iterator->At(i)->GetTitle());
#endif
 
    // Extract the name
    sscanf(iterator->At(i)->GetTitle(), "%s (Path: %*s)", name);
   
    // Delete ".C" at the end 
    // -> Note: Physical address pointer, do NOT delete. / Changes "name" as well!
    Char_t* dotC = (Char_t*)strrchr(name, '.');
    if (dotC != 0)
    {
      *dotC = '\0';
      dotC++;
      *dotC = '\0';
    }

    // Add to "data-from-list"
    fDataFromMacroList->Add(new TObjString(name));

    // Find the type of the process macro
    if (!IsHistogramMacro(name))
    {
      // Type 1
      isHistoMacro[i] = kFALSE;
      // Create the command 
      sprintf(cmds[i], "%s(automaticTrackV1, results, n);", name);
    }
    else
    {
      // Type 2 (histo)
      isHistoMacro[i] = kTRUE;
      numHistoMacros++;
      // Create the command 
      sprintf(cmds[i], "%s(automaticTrackV1);", name);
    } 
  }  

  // Allocate memory for the histograms
  if (numHistoMacros > 0)  histos = new TH1*[numHistoMacros];
  for (Int_t i = 0; i < numHistoMacros; i++)  histos[i] = 0;
  
  // Walk through the list of tracks     
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter)
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
    for (Int_t i = 0, histoIndex = 0; i < iterator->GetEntries(); i++)
    {
      // Process for macro type 2 (histo)
      if (isHistoMacro[i])
      {
        returnedHist = (TH1*)gROOT->ProcessLineSync(cmds[i]);
        if (returnedHist != 0x0)
        {
          if (histos[histoIndex] == 0)  histos[histoIndex] = (TH1*)gROOT->ProcessLineSync(cmds[i]);
          else  histos[histoIndex]->Add((const TH1*)gROOT->ProcessLineSync(cmds[i]));
      
          delete returnedHist;
          returnedHist = 0;
        }
        histoIndex++;
      }
      // Process for macro type 1
      else
      {
        // Create data pointers in CINT, execute the macro and get the data
        gROOT->ProcessLineSync("Double_t* results = 0;");
        gROOT->ProcessLineSync("Int_t n = 0;");
        gROOT->ProcessLineSync(cmds[i]);
        Double_t* results = (Double_t*)gROOT->ProcessLineSync("results;");
        Int_t nResults = (Int_t)gROOT->ProcessLineSync("n;");
        
        if (results == 0)
        {
          Error("Apply macros", Form("Error reading data from macro \"%s\"", iterator->At(i)->GetTitle()));
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
  }    

  for (Int_t i = 0, histoIndex = 0; i < iterator->GetEntries() && histoIndex < numHistoMacros; i++)
  {
    if (isHistoMacro[i])
    {
      // Might be empty (e.g. no tracks have been selected)!
      if (histos[histoIndex] != 0)
      {
        (*fDataTree) << Form("TrackData%d", i) << Form("Macro%d=", i) << histos[histoIndex] << (Char_t*)"\n";
      }
      histoIndex++;
    }
  }

  if (fDataTree != 0) delete fDataTree;
  fDataTree = 0;

  if (cmds != 0)  delete [] cmds;
  if (isHistoMacro != 0)  delete isHistoMacro;
  isHistoMacro = 0;

  if (histos != 0)  delete [] histos;
  histos = 0;

  // Clear root
  gROOT->Reset();
  
  // If there is data, select the first data set
  if (iterator->GetEntries() > 0) SETBIT(fHistoDataSelected, 0);

  // Now the data is stored in "TRD.TrackListMacroData.root"
  // The editor will access this file to display the data
}

//______________________________________________________
void AliEveTRDTrackList::ApplySelectionMacros(TList* iterator)
{
  Char_t name[fkMaxMacroNameLength];
  Char_t cmd[(fkMaxMacroNameLength + fkMaxApplyCommandLength)];

  AliEveTRDTrack* track = 0;
  AliTRDtrackV1 *trackv1 = 0;
  Bool_t selectedByMacro = kFALSE;

  // Clear root
  gROOT->Reset();

  // Select all tracks at first. A track is then deselect, if at least one selection macro
  // returns kFALSE for this track
  // Enable all tracks (Note: EnableListElements(..) will call "ElementChanged", which will cause unforseen behaviour!)
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter)
  {
    ((TEveElement*)(*iter))->SetRnrState(kTRUE);
  }
  SetRnrState(kTRUE);
  
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {

    memset(name, '\0', sizeof(Char_t) * fkMaxMacroNameLength);
    memset(cmd, '\0', sizeof(Char_t) * (fkMaxMacroNameLength + fkMaxApplyCommandLength));

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Applying selection macro: %s\n", iterator->At(i)->GetTitle());
#endif
    
    // Extract the name
    sscanf(iterator->At(i)->GetTitle(), "%s (Path: %*s)", name);
    // Delete ".C" at the end 
    // -> Note: Physical address pointer, do NOT delete. / Changes "name" as well!
    Char_t* dotC = (Char_t*)strrchr(name, '.');
    if (dotC != 0)
    {
      *dotC = '\0';
      dotC++;
      *dotC = '\0';
    }

    // Create the command
    sprintf(cmd, "%s(automaticTrackV1);", name);

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
      track->SetRnrState(selectedByMacro && track->GetRnrState());         
    }
  }

  // Clear root
  gROOT->Reset();  
}

//______________________________________________________
Char_t* AliEveTRDTrackList::MakeMacroEntry(const Char_t* path, const Char_t* name)
{
  Char_t* entry = new Char_t[(fkMaxMacroPathNameLength + 30)];
  memset(entry, '\0', sizeof(Char_t) * (fkMaxMacroPathNameLength + 30));

  Char_t* systemPath = gSystem->ExpandPathName(path);
  sprintf(entry, "%s (Path: %s)", name, systemPath);
  delete systemPath;
  systemPath = 0;

  return entry;
}

//______________________________________________________
Bool_t AliEveTRDTrackList::IsHistogramMacro(const Char_t* name) 
{
  TFunction* f = 0x0;
  if ((f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*", kTRUE)) != 0x0)
    if (strcmp(f->GetReturnTypeName(), "TH1*") == 0)  return kTRUE;

  return kFALSE;  
}

//______________________________________________________
void AliEveTRDTrackList::RemoveProcessMacros(TList* iterator) 
{
  TObjString* obj = 0;
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    obj = (TObjString*)fMacroList->Remove(fMacroList->FindObject(iterator->At(i)->GetTitle()));
    
    if (obj != 0) delete obj;
  }
  obj = 0;
}

//______________________________________________________
void AliEveTRDTrackList::RemoveSelectionMacros(TList* iterator) 
{
  TObjString* obj = 0;
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    obj = (TObjString*)fMacroSelList->Remove(fMacroSelList->FindObject(iterator->At(i)->GetTitle()));
    if (obj != 0) delete obj;
  }
  obj = 0;
}

