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
  fileInfo(0),
  fileTypes(0)
{
  // Constructor.
  fMainFrame = CreateEditorTabSubFrame("Apply macros");
  
  // Functionality for adding macros 
  fBrowseFrame = new TGHorizontalFrame(fMainFrame);

  teField = new TGTextEntry(fBrowseFrame);
  teField->Connect("ReturnPressed()","AliEveTRDTrackListEditor", this, "HandleMacroPathSet()"); 
  fBrowseFrame->AddFrame(teField);
  
  bBrowse = new TGTextButton(fBrowseFrame, "Browse");
  bBrowse->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "BrowseMacros()");
  fBrowseFrame->AddFrame(bBrowse);
  fMainFrame->AddFrame(fBrowseFrame);

  tlMacroList = new TGListBox(fMainFrame);
  tlMacroList->Resize(194, 120);
  tlMacroList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(tlMacroList);

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

void AliEveTRDTrackListEditor::AddMacro(const Char_t* pathname)
{
  // Only add macro, if it is not already in the list
  if (fM->macroList->FindObject(pathname) == 0)
  {
      fM->macroList->Add(new TObjString(pathname));
      fM->macroList->Sort();

      UpdateMacroList();
  }
  else
  {
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Warning", "Macro is already in list (won't be added again)!",
                   kMBIconExclamation, kMBOk);
  }
}

void AliEveTRDTrackListEditor::ApplyMacros()
{
  TList* iterator = new TList();
  //TEveMacro* macro = NULL;
  Char_t name[100];
  Char_t path[300];
  Char_t pathname[400];
  Char_t cmd[430];

  memset(name, '\0', sizeof(Char_t) * 100);
  memset(path, '\0', sizeof(Char_t) * 300);
  memset(pathname, '\0', sizeof(Char_t) * 400);
  memset(cmd, '\0', sizeof(Char_t) * 430);

  tlMacroList->GetSelectedEntries(iterator);

  // Make tracklist availabe
  fM->ExportToCINT("trackList");
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
      // Extract path and name -> Make pathname
      sscanf(iterator->At(i)->GetTitle(), "%s (Path: %s)", name, path);
      // Delete ")" at the end of path
      path[strlen(path) - 1] = '\0';
      sprintf(pathname, "%s/%s", path, name);

      TString sPathname(pathname);
      gSystem->ExpandPathName(sPathname);
 
      sprintf(cmd, ".x %s(trackList)", pathname);
      gROOT->ProcessLine(cmd);
      // Load and execute macro
      //macro = new TEveMacro(sPathname);

      //macro->Exec("trackList");
      //if (macro != NULL)    delete macro;
  }

  if (iterator != NULL) delete iterator;    
}
/*
void AliEveTRDTrackListEditor::BrowseMacros()
{
  new TGFileDialog(gClient->GetRoot(), new TGWindow(), kFDOpen, fileInfo);
  
  if (fileInfo->fFilename != 0 && fileInfo->fIniDir != 0)
  {        			
    Char_t pathname[300];
    
    // Extract filename
    Char_t* name = strrchr(fileInfo->fFilename, '/');
    // Delete '"' at the end
    name[strlen(name)] = '\0';
    sprintf(pathname, "%s (Path: %s)", name + 1, fileInfo->fIniDir);
  
    AddMacro(pathname); 
  }
}
*/

void AliEveTRDTrackListEditor::BrowseMacros()
{
  new TGFileDialog(gClient->GetRoot(), GetMainFrame(), kFDOpen, fileInfo);
  
  if (fileInfo->fIniDir != 0)
  {        			
    Char_t pathname[300];
    
    // Extract filenames
    TObject* iter = fileInfo->fFileNamesList->First();
 
    Char_t* name = 0;

    while (iter != 0)
    {
      name = strrchr(iter->GetName(), '/');
      // Delete '"' at the end
      name[strlen(name)] = '\0';
      sprintf(pathname, "%s (Path: %s)", name + 1, fileInfo->fIniDir);
  
      AddMacro(pathname); 
      iter = (TObjString*)fileInfo->fFileNamesList->After(iter);
    }
  }
}

void AliEveTRDTrackListEditor::HandleMacroPathSet()
{
  if (strlen(teField->GetText()) != 0)
  {         			
    // Check if file exists
    FILE* fp = NULL;

    fp = fopen(teField->GetText(), "rb");
    if (fp != NULL)
    {
        fclose(fp);
        Char_t pathname[300];
    
        // Extract filename
        Char_t* name = strrchr(teField->GetText(), '/');

        // Current path
        if (name == NULL)
        {
            sprintf(pathname, "%s (Path: .)", teField->GetText());
        }
        // Different path
        else
        {
            // Extract path
            Char_t* path = new Char_t[240];
            memset(path, '\0', sizeof(Char_t) * 240);
            strncpy(path, teField->GetText(), strlen(teField->GetText()) - strlen(name));
            sprintf(pathname, "%s (Path: %s)", name + 1, path);

            if (path != NULL)  delete path;
        }

        AddMacro(pathname);   
    }
    else
    {
        new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                     "File does not exist or you do not have read permission!", kMBIconExclamation, kMBOk);
    }
  }
}

void AliEveTRDTrackListEditor::RemoveMacros()
{
  TList* iterator = new TList();
  
  tlMacroList->GetSelectedEntries(iterator);

  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
      fM->macroList->Remove(fM->macroList->FindObject(iterator->At(i)->GetTitle()));
  }

  UpdateMacroList();
}

void AliEveTRDTrackListEditor::SetModel(TObject* obj)
{
  // Set model object.
  fM = dynamic_cast<AliEveTRDTrackList*>(obj);

  // Add members to a list
  tvMemberList->Clear();
  for (TEveElement::List_i iterator = fM->BeginChildren(); iterator != fM->EndChildren(); ++iterator)
  {
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
}
