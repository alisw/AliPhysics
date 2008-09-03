#include "AliEveTRDTrackListEditor.h"

ClassImp(AliEveTRDTrackListEditor)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackListEditor //////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackListEditor::AliEveTRDTrackListEditor(const TGWindow* p, Int_t width, Int_t height,
				                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options, back),
  fM(0),
  fMainFrame(0),
  fMemberFrame(0),
  fBrowseFrame(0),
  bBrowse(0),
  bApplyMacros(0),
  bRemoveMacros(0),
  teField(0),
  tvMemberList(0),
  tlMacroList(0),
  tlMacroSelList(0),
  fileInfo(0),
  fileTypes(0),
  fLabel1(0), fLabel2(0), fLabel3(0),
  fLine1(0), fLine2(0), fLine3(0)
{
  // Constructor.
  fMainFrame = CreateEditorTabSubFrame("Apply macros");
 
  // Functionality for adding macros 
  fLabel1 = new TGLabel(fMainFrame,"Add macro(s):");
  fMainFrame->AddFrame(fLabel1);
  fBrowseFrame = new TGHorizontalFrame(fMainFrame);

  teField = new TGTextEntry(fBrowseFrame);
  teField->Connect("ReturnPressed()","AliEveTRDTrackListEditor", this, "HandleMacroPathSet()"); 
  fBrowseFrame->AddFrame(teField);
  
  bBrowse = new TGTextButton(fBrowseFrame, "Browse");
  bBrowse->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "BrowseMacros()");
  fBrowseFrame->AddFrame(bBrowse);
  fMainFrame->AddFrame(fBrowseFrame);

  fLine1 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel2 = new TGLabel(fMainFrame,"Selection macros:");
  fMainFrame->AddFrame(fLabel2);

  tlMacroSelList = new TGListBox(fMainFrame);
  tlMacroSelList->Resize(194, 94);
  tlMacroSelList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(tlMacroSelList);

  fLine2 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel3 = new TGLabel(fMainFrame,"Process macros:");
  fMainFrame->AddFrame(fLabel3);

  tlMacroList = new TGListBox(fMainFrame);
  tlMacroList->Resize(194, 94);
  tlMacroList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(tlMacroList);

  fLine3 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine3, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));  

  bApplyMacros = new TGTextButton(fMainFrame, "Apply selected macro(s)");
  bApplyMacros->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "ApplyMacros()");
  bApplyMacros->SetRightMargin(12);
  fMainFrame->AddFrame(bApplyMacros);

  bRemoveMacros = new TGTextButton(fMainFrame, "Remove selected macro(s)");
  bRemoveMacros->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "RemoveMacros()");
  fMainFrame->AddFrame(bRemoveMacros);

  // List members
  fMemberFrame = CreateEditorTabSubFrame("Members");
  
  tvMemberList = new TGTextView(fMemberFrame, 220, 120, "");
  fMemberFrame->AddFrame(tvMemberList);

  // Set up file dialog
  fileInfo = new TGFileInfo();
  fileInfo->SetMultipleSelection(kTRUE);

  fileTypes = new Char_t*[6];
  fileTypes[0] = "All files"; fileTypes[1] = "*";
  fileTypes[2] = "ROOT macros"; fileTypes[3] = "*.C";
  fileTypes[4] = 0; fileTypes[5] = 0;
  fileInfo->fFileTypes = (const Char_t**)fileTypes;
  fileInfo->fFileTypeIdx = 2;
  fileInfo->fMultipleSelection = kTRUE;
}


//________________________________________________________
void AliEveTRDTrackListEditor::AddMacro(const Char_t* entryName, const Char_t* nameC, const Char_t* pathname)
{
  // First check the type of the macro:
  // If it has the signature of a selection macro:
  // Bool_t MacroName(AliEveTRDTrack
  // it is assumed to be a selection macro. In all other cases: Process macro
  TEveMacro* macro = NULL;
  Bool_t isSelectionMacro = kFALSE;
  Char_t signature[120];
  memset(signature, '\0', sizeof(Char_t) * 120);

  // Delete ".C" from filename
  Char_t* name = new Char_t[strlen(nameC)];
  memset(name, '\0', sizeof(Char_t) * strlen(nameC));
  strncpy(name, nameC, strlen(nameC) - 2);
  
  // Create signature
  sprintf(signature, "Bool_t %s(AliEveTRDTrack", name);
  
  macro = new TEveMacro(pathname);
  if (!macro) {
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Cannot access file!", kMBIconExclamation, kMBOk);
    return;
  }
  if (macro->TMacro::GetLineWith(signature) != 0)   isSelectionMacro = kTRUE;
  else                                              isSelectionMacro = kFALSE; 
  delete macro;
  macro = 0;
  
  // Only add macro, if it is not already in the list
  if (!isSelectionMacro && fM->macroList->FindObject(entryName) == 0){
    fM->macroList->Add(new TObjString(entryName));
    fM->macroList->Sort();

    UpdateMacroList();
  } else if (isSelectionMacro && fM->macroSelList->FindObject(entryName) == 0) {
    fM->macroSelList->Add(new TObjString(entryName));
    fM->macroSelList->Sort();

    UpdateMacroList();
  } else {
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Warning", "Macro is already in list (won't be added again)!", kMBIconExclamation, kMBOk);
  }
}

//________________________________________________________
void AliEveTRDTrackListEditor::ApplyMacros()
{
  TList* iterator = new TList();
  TEveMacro* macro = NULL;
  Char_t name[100];
  Char_t path[300];
  Char_t pathname[400];
  //Char_t cmd[430];

  // First apply the selection macros
  tlMacroSelList->GetSelectedEntries(iterator);
  AliEveTRDTrack* track = 0;
  //AliTRDtrackV1 *trackv1 = 0;
  Bool_t selectedByMacro = kFALSE;

  for (Int_t i = 0; i < iterator->GetEntries(); i++){
    memset(name, '\0', sizeof(Char_t) * 100);
    memset(path, '\0', sizeof(Char_t) * 300);
    memset(pathname, '\0', sizeof(Char_t) * 400);
    
    // Extract path and name -> Make pathname
    sscanf(iterator->At(i)->GetTitle(), "%s (Path: %s)", name, path);
    // Delete ")" at the end of path
    path[strlen(path) - 1] = '\0';
    sprintf(pathname, "%s/%s", path, name);
    
    TString sPathname(pathname);
    gSystem->ExpandPathName(sPathname);

    // Load and execute macro
    macro = new TEveMacro(sPathname);

    // Walk through the list of tracks
    for (TEveElement::List_i iter = fM->BeginChildren(); iter != fM->EndChildren(); ++iter){
      track = dynamic_cast<AliEveTRDTrack*>(*iter);

      if (!track) continue;
  
      //trackv1 = (AliTRDtrackV1*)track->GetUserData();
      
      // Make this track available
      track->ExportToCINT("automaticTrack");
      selectedByMacro = (Bool_t)macro->Exec("automaticTrack");
      track->SetRnrState(selectedByMacro);         

      if (macro != NULL)    delete macro;
      return;
    }
  }  

  if (iterator != NULL) delete iterator;  

  // Now apply the process macros
  iterator = new TList();
  tlMacroList->GetSelectedEntries(iterator);

  // Make tracklist availabe
  fM->ExportToCINT("trackList");
  for (Int_t i = 0; i < iterator->GetEntries(); i++){
    memset(name, '\0', sizeof(Char_t) * 100);
    memset(path, '\0', sizeof(Char_t) * 300);
    memset(pathname, '\0', sizeof(Char_t) * 400);
    //memset(cmd, '\0', sizeof(Char_t) * 430);

    // Extract path and name -> Make pathname
    sscanf(iterator->At(i)->GetTitle(), "%s (Path: %s)", name, path);
    // Delete ")" at the end of path
    path[strlen(path) - 1] = '\0';
    sprintf(pathname, "%s/%s", path, name);

    TString sPathname(pathname);
    gSystem->ExpandPathName(sPathname);

    //sprintf(cmd, ".x %s(trackList)", pathname);
    //gROOT->ProcessLine(cmd);
    // Load and execute macro
    macro = new TEveMacro(sPathname);

    macro->Exec("trackList");
    if (macro != NULL)    delete macro;
  }

  if (iterator != NULL) delete iterator;    
}


//________________________________________________________
void AliEveTRDTrackListEditor::BrowseMacros()
{
  new TGFileDialog(gClient->GetRoot(), GetMainFrame(), kFDOpen, fileInfo);
  
  if (fileInfo->fIniDir != 0 && fileInfo->fFileNamesList != 0)
  {       
    Char_t entryName[300];
    memset(entryName, '\0', sizeof(Char_t) * 300);
    
    // Extract filenames
    TObject* iter = fileInfo->fFileNamesList->First();
 
    Char_t* name = 0;

    while (iter != 0)
    {
      name = strrchr(iter->GetName(), '/');
      // Delete '"' at the end
      name[strlen(name)] = '\0';
      sprintf(entryName, "%s (Path: %s)", name + 1, fileInfo->fIniDir);
        
      AddMacro(entryName, name + 1, iter->GetName()); 
      iter = (TObjString*)fileInfo->fFileNamesList->After(iter);
    }
  }

  // Some error occurs, when one ends the filedialog with "cancel": fileInfo->fFileNamesList is set to 0x0, but
  // in the next launch no new memory is allocated. So do this manually.
  if (fileInfo->fFileNamesList == 0)  fileInfo->fFileNamesList = new TList();
}


//________________________________________________________
void AliEveTRDTrackListEditor::HandleMacroPathSet()
{
  if (strlen(teField->GetText()) != 0){         			
    // Check if file exists
    FILE* fp = NULL;

    fp = fopen(teField->GetText(), "rb");
    if (fp != NULL){
      fclose(fp);
      Char_t entryName[300];
      memset(entryName, '\0', sizeof(Char_t) * 300);
  
      // Extract filename
      Char_t* name = strrchr(teField->GetText(), '/');

      // Current path
      if (name == NULL){
        sprintf(entryName, "%s (Path: .)", teField->GetText());
        sprintf(name, "%s", teField->GetText());

        // Add path to textfield
        Char_t pathname[100];
        memset(pathname, '\0', sizeof(Char_t) * 100);
        sprintf(pathname, "./%s", teField->GetText());
        teField->SetText(pathname);
      } else {// Different path
        // Extract path
        Char_t* path = new Char_t[240];
        memset(path, '\0', sizeof(Char_t) * 240);
        strncpy(path, teField->GetText(), strlen(teField->GetText()) - strlen(name));
        sprintf(entryName, "%s (Path: %s)", name + 1, path);

        if (path != NULL)  delete path;

        // Ignore the slash "/" in the following
        name++;
      }
      AddMacro(entryName, name, teField->GetText());   
    } else {
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
      "File does not exist or you do not have read permission!", kMBIconExclamation, kMBOk);
    }
  }
}


//________________________________________________________
void AliEveTRDTrackListEditor::RemoveMacros()
{
  TList* iterator = new TList();
  
  tlMacroList->GetSelectedEntries(iterator);

  for (Int_t i = 0; i < iterator->GetEntries(); i++){
      fM->macroList->Remove(fM->macroList->FindObject(iterator->At(i)->GetTitle()));
  }

  tlMacroSelList->GetSelectedEntries(iterator);

  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
      fM->macroSelList->Remove(fM->macroSelList->FindObject(iterator->At(i)->GetTitle()));
  }

  UpdateMacroList();
}

void AliEveTRDTrackListEditor::SetModel(TObject* obj)
{
  // Set model object.
  fM = dynamic_cast<AliEveTRDTrackList*>(obj);

  // Add members to a list
  tvMemberList->Clear();
  // In order to prevent the first line from being empty, do the following:
  TEveElement::List_i iterator = fM->BeginChildren();
  if (iterator != fM->EndChildren())
  {
    tvMemberList->SetText(new TGText(((AliEveTRDTrack*)(*iterator))->GetName()));
    iterator++;
    for ( ; iterator != fM->EndChildren(); ++iterator)
            tvMemberList->AddLineFast(((AliEveTRDTrack*)(*iterator))->GetName());
  }
  
  tvMemberList->ShowTop();
  tvMemberList->Update();

  UpdateMacroList();
}

void AliEveTRDTrackListEditor::UpdateMacroList()
{
  tlMacroList->RemoveAll();
 
  TObjString* iter = (TObjString*)fM->macroList->First();

  while (iter != 0)
  {
    tlMacroList->AddEntry(iter->GetName(), -1);
    iter = (TObjString*)fM->macroList->After(iter);
  }

  tlMacroList->SortByName();


  tlMacroSelList->RemoveAll();
 
  iter = (TObjString*)fM->macroSelList->First();

  while (iter != 0)
  {
    tlMacroSelList->AddEntry(iter->GetName(), -1);
    iter = (TObjString*)fM->macroSelList->After(iter);
  }

  tlMacroSelList->SortByName();  
}
