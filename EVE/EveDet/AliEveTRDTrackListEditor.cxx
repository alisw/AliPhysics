#include "TTreeStream.h"
#include <EveDet/AliEveTRDTrackList.h>
#include "AliEveTRDTrackListEditor.h"

#include <TGFileDialog.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGedEditor.h>     ////// MAYBE THIS CAN BE REMOVED
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGListBox.h>
#include <TGMsgBox.h>
#include <TGLabel.h>
#include <TG3DLine.h>
#include <TEveMacro.h>
#include <TEveManager.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <AliTRDtrackV1.h>
#include <EveDet/AliEveTRDData.h>

ClassImp(AliEveTRDTrackListEditor)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackListEditor //////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackListEditor::AliEveTRDTrackListEditor(const TGWindow* p, Int_t width, Int_t height,
				                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options, back),
  fM(0),
  fMainFrame(0),
  fHistoFrame(0),
  fHistoSubFrame(0),
  fBrowseFrame(0),
  fbBrowse(0),
  fbApplyMacros(0),
  fbRemoveMacros(0),
  fbDrawHisto(0),
  fteField(0),
  ftlMacroList(0),
  ftlMacroSelList(0),
  fFileInfo(0),
  fFileTypes(0),
  fLabel1(0), fLabel2(0), fLabel3(0), fLabel4(0),
  fLine1(0), fLine2(0), fLine3(0), fLine4(0),
  fCheckButtons(0)
{  
  fMainFrame = CreateEditorTabSubFrame("Apply macros");
 
  // Functionality for adding macros 
  fLabel1 = new TGLabel(fMainFrame,"Add macro(s):");
  fMainFrame->AddFrame(fLabel1);
  fBrowseFrame = new TGHorizontalFrame(fMainFrame);

  fteField = new TGTextEntry(fBrowseFrame);
  fteField->Connect("ReturnPressed()","AliEveTRDTrackListEditor", this, "HandleMacroPathSet()"); 
  fBrowseFrame->AddFrame(fteField);
  
  fbBrowse = new TGTextButton(fBrowseFrame, "Browse");
  fbBrowse->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "BrowseMacros()");
  fBrowseFrame->AddFrame(fbBrowse);
  fMainFrame->AddFrame(fBrowseFrame);

  fLine1 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine1, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));
  fLabel2 = new TGLabel(fMainFrame,"Selection macros:");
  fMainFrame->AddFrame(fLabel2);

  ftlMacroSelList = new TGListBox(fMainFrame);
  ftlMacroSelList->Resize(194, 94);
  ftlMacroSelList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(ftlMacroSelList);

  fLine2 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine2, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));
  fLabel3 = new TGLabel(fMainFrame,"Process macros:");
  fMainFrame->AddFrame(fLabel3);

  ftlMacroList = new TGListBox(fMainFrame);
  ftlMacroList->Resize(194, 94);
  ftlMacroList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(ftlMacroList);

  fLine3 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine3, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));  

  fbApplyMacros = new TGTextButton(fMainFrame, "Apply selected macro(s)");
  fbApplyMacros->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "ApplyMacros()");
  fbApplyMacros->SetRightMargin(12);
  fMainFrame->AddFrame(fbApplyMacros);

  fbRemoveMacros = new TGTextButton(fMainFrame, "Remove selected macro(s)");
  fbRemoveMacros->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "RemoveMacros()");
  fMainFrame->AddFrame(fbRemoveMacros);

  // Stuff for displaying histograms
  fHistoFrame = CreateEditorTabSubFrame("Histograms");  
  fHistoFrame->SetMapSubwindows(kTRUE);
  fLabel4 = new TGLabel(fHistoFrame,"Data from applied macros:");
  fHistoFrame->AddFrame(fLabel4);

  fHistoSubFrame = new TGVerticalFrame(fHistoFrame);
  fHistoSubFrame->SetMapSubwindows(kTRUE);
  fHistoSubFrame->Resize(194, 200);
  fHistoFrame->AddFrame(fHistoSubFrame);

  fLine4 = new TGHorizontal3DLine(fHistoFrame, 194, 8);
  fHistoFrame->AddFrame(fLine4, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));  

  fbDrawHisto = new TGTextButton(fHistoFrame, "Draw histogram");
  fbDrawHisto->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "DrawHistos()");
  fHistoFrame->AddFrame(fbDrawHisto);

  // Set up file dialog
  fFileInfo = new TGFileInfo();
  fFileInfo->SetMultipleSelection(kTRUE);

  fFileTypes = new Char_t*[6];
  fFileTypes[0] = (Char_t*)"All files"; fFileTypes[1] = (Char_t*)"*";
  fFileTypes[2] = (Char_t*)"ROOT macros"; fFileTypes[3] = (Char_t*)"*.C";
  fFileTypes[4] = 0; fFileTypes[5] = 0;
  fFileInfo->fFileTypes = (const Char_t**)fFileTypes;
  fFileInfo->fFileTypeIdx = 2;
  fFileInfo->fMultipleSelection = kTRUE;

  // Set focus on "Apply macros" tab
  //fMainFrame->TGWindow::RequestFocus();
  //Update();
}

//______________________________________________________
AliEveTRDTrackListEditor::~AliEveTRDTrackListEditor()
{
  if (fFileTypes != 0)
  {
    delete [] fFileTypes;
    fFileTypes = 0;
  }

  if (fFileInfo != 0)
  {
    delete fFileInfo; 
    fFileInfo = 0;
  }
}

//______________________________________________________
void AliEveTRDTrackListEditor::AddMacro(const Char_t* path, const Char_t* name)
{
  Int_t result = fM->AddMacro(path, name);

  switch (result)
  {
  case SUCCESS:
    UpdateMacroList();
    break;
  case WARNING:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Warning", "Macro is already in list (won't be added again)!",
                 kMBIconExclamation, kMBOk);
    break;
  case SIGNATURE_ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Macro has not the signature of...\n...a selection macro: Bool_t YourMacro(AliTRDtrackV1*)\n...a process macro: void YourMacro(AliTRDtrackV1*, Double_t*&, Int_t&)", 
                 kMBIconExclamation, kMBOk);
    break;               
  case NOT_EXIST_ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "File does not exist or you do not have read permission!", kMBIconExclamation, kMBOk);
    break;
  default:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 Form("AliEveTRDTrackList::AddMacro exited with unknown return value: %d", result),
                 kMBIconExclamation, kMBOk);
    break;
  }
}

//______________________________________________________
void AliEveTRDTrackListEditor::ApplyMacros()
{
  // First apply the selection macros
  TList* iterator = new TList();
  ftlMacroSelList->GetSelectedEntries(iterator);
  fM->ApplySelectionMacros(iterator);
  
  // Update view
  gEve->Redraw3D();

  if (iterator != 0) delete iterator;  

  // Now apply the process macros
  iterator = new TList();
  ftlMacroList->GetSelectedEntries(iterator);
  fM->ApplyProcessMacros(iterator);
  
  // Now the data is stored in the following file and can easily be accessed
  TFile* file = new TFile("TRD.TrackListMacroData.root", "READ");
  if (!file)  
  {
    Error("Apply macros", "Cannot open file \"TRD.TrackListMacroData.root\"");
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Apply macros", 
                 "Cannot open file \"TRD.TrackListMacroData.root\"", kMBIconExclamation, kMBOk);
    
    if (iterator != 0)  delete iterator;  
    iterator = 0;

    return;  
  }
    
  TTree* t = 0;
  for (Int_t i = 0; i < iterator->GetEntries(); i++) {
    t = (TTree*)file->Get(Form("TrackData%d", i));
    if (t != 0) {
      gEve->AddCanvasTab(Form("Macro%d", i));
      t->Draw(Form("Macro%d", i), "1");
 
      delete t;
      t = 0;

      // ONLY DISPLAY ONE MACRO (the first one possible) -> Remove the next line to display all
      break;
    } else {
      Error("Apply macros", Form("No data for macro%d found!", i));
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Apply macros", 
                   Form("No data for macro%d found!", i), kMBIconExclamation, kMBOk);  
    }
  }
  file->Close("R");
  delete file;
  file = 0;

  if (iterator != 0)  delete iterator;  
  iterator = 0;  

  // Update histogram tab (data has to be reloaded)
  //fHistoFrame->TGWindow::RequestFocus();
  //fGedEditor->TGCompositeFrame::ShowFrame(fHistoFrame);
  SetModel(fM);
  Update();
}


//______________________________________________________
void AliEveTRDTrackListEditor::BrowseMacros()
{
  new TGFileDialog(gClient->GetRoot(), GetMainFrame(), kFDOpen, fFileInfo);
  
  if (fFileInfo->fIniDir != 0 && fFileInfo->fFileNamesList != 0)
  {       
    // Extract filenames
    TObject* iter = fFileInfo->fFileNamesList->First();
 
    Char_t* name = 0;

    while (iter != 0)
    {
      // NOTE: fileInfo->fFileNamesList will be changed by that, too!
      name = strrchr(iter->GetName(), '/');
      // Delete '"' at the end
      name[strlen(name)] = '\0';
              
      AddMacro(fFileInfo->fIniDir, name + 1); 
      iter = (TObjString*)fFileInfo->fFileNamesList->After(iter);
    }
  }

  // -> The following problem has been fixed (trunk -> Changes according to 03 September 2008):
  // Some error occurs, when one ends the filedialog with "cancel": fileInfo->fFileNamesList is set to 0x0, but
  // in the next launch no new memory is allocated. So do this manually.
  //if (fileInfo->fFileNamesList == 0)  fileInfo->fFileNamesList = new TList();
}

//______________________________________________________
void AliEveTRDTrackListEditor::DrawHistos()
{
  Int_t nHistograms = GetNSelectedHistograms();
  if (nHistograms <= 0)
  {
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "No data selected. Please select the data you want to plot!", kMBIconExclamation, kMBOk);
    return;
  }
  if (nHistograms > 3)
  {
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), 
                 "Error", "Only histograms with up to 3 dimensions supported. Please select 1,2 or 3 data macros!",
                 kMBIconExclamation, kMBOk);
    return;
  }

  TFile* file = new TFile("TRD.TrackListMacroData.root", "READ");
  if (!file)  
  {
    Error("Draw histograms", "Cannot open file \"TRD.TrackListMacroData.root\"");
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                 "Cannot open file \"TRD.TrackListMacroData.root\"", kMBIconExclamation, kMBOk);
    return;
  }
  
  TTree* t = 0;
  TTree* tFriend1 = 0;
  TTree* tFriend2 = 0;

  Int_t indexOfMacro1 = 0;
  Int_t indexOfMacro2 = 0;
  Int_t indexOfMacro3 = 0;

  // Load the trees in succession and remember the entries
  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    if (fCheckButtons[i]->TGButton::GetState() == kButtonDown)
    {
      if (t == 0)
      {
        indexOfMacro1 = i;
        if (!(t = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro%d found!", i));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro%d found!", i), kMBIconExclamation, kMBOk);
          break;   
        }
     
        // 1d histogram   
        if (nHistograms == 1) 
        {
          gEve->AddCanvasTab(Form("Macro%d", indexOfMacro1));
          t->Draw(Form("Macro%d", indexOfMacro1), "1");

          break;     
        }
      }
      else if (tFriend1 == 0)
      {
        indexOfMacro2 = i;
        if (!(tFriend1 = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro%d found!", i));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro%d found!", i), kMBIconExclamation, kMBOk);
          break;   
        }
        
        // 2d histogram
        if (nHistograms == 2) 
        {
          gEve->AddCanvasTab(Form("Macro%d - Macro%d", indexOfMacro1, indexOfMacro2));
          t->AddFriend(tFriend1);
          t->Draw(Form("Macro%d:Macro%d", indexOfMacro1, indexOfMacro2), "1");
 
          break;     
        }
      }    
      // 3d histogram
      else
      {
        indexOfMacro3 = i;
        if (!(tFriend2 = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro%d found!", i));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro%d found!", i), kMBIconExclamation, kMBOk);
          break;   
        }
        
        gEve->AddCanvasTab(Form("Macro%d - Macro%d - Macro%d", indexOfMacro1, indexOfMacro2, indexOfMacro3));
        t->AddFriend(tFriend1);
        t->AddFriend(tFriend2);
        t->Draw(Form("Macro%d:Macro%d:Macro%d", indexOfMacro1, indexOfMacro2, indexOfMacro3), "1");
 
        break;     
      }
    }
  }

  if (t != 0) delete t;
  t = 0;
  if (tFriend1 != 0)  delete tFriend1;
  tFriend1 = 0;
  if (tFriend2 != 0)  delete tFriend2;
  tFriend2 = 0;

  file->Close("R");
  delete file;
  file = 0;
}

//______________________________________________________
Int_t AliEveTRDTrackListEditor::GetNSelectedHistograms()
{
  Int_t count = 0;
  
  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    if (fCheckButtons[i]->TGButton::GetState() == kButtonDown)  count++;
  }

  return count;
}

//______________________________________________________
void AliEveTRDTrackListEditor::HandleMacroPathSet()
{
  if (strlen(fteField->GetText()) != 0)
  {         			
    // Check if file exists
    FILE* fp = NULL;

    fp = fopen(fteField->GetText(), "rb");
    if (fp != NULL)
    {
      fclose(fp);

      // Extract filename
      Char_t* name = strrchr(fteField->GetText(), '/');

      // Current path
      if (name == NULL)
      {
        name = new Char_t[100];
        memset(name, '\0', sizeof(Char_t) * 100);
        sprintf(name, "%s", fteField->GetText());

        // Add path to textfield
        Char_t pathname[100];
        memset(pathname, '\0', sizeof(Char_t) * 100);
        sprintf(pathname, "./%s", fteField->GetText());
        fteField->SetText(pathname);

        AddMacro(".", name);  
        if (name != 0)  delete name;
        name = 0;
      }
      // Different path
      else
      {
        // Extract path
        Char_t* path = new Char_t[240];
        memset(path, '\0', sizeof(Char_t) * 240);
        strncpy(path, fteField->GetText(), strlen(fteField->GetText()) - strlen(name));
        
        // Ignore the slash "/" in name
        AddMacro(path, name + 1);  
  
        if (path != 0)  delete path;
        path = 0;
      }       
    }
    else
    {
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                   "File does not exist or you do not have read permission!", kMBIconExclamation, kMBOk);
    }
  }
}

//______________________________________________________
void AliEveTRDTrackListEditor::RemoveMacros()
{
  TList* iterator = new TList();
  
  ftlMacroList->GetSelectedEntries(iterator);
  fM->RemoveProcessMacros(iterator);

  if (iterator != 0)  delete iterator;

  iterator = new TList();
  ftlMacroSelList->GetSelectedEntries(iterator);
  fM->RemoveSelectionMacros(iterator);

  UpdateMacroList();

  if (iterator != 0)  delete iterator;
  iterator = 0;
}


//______________________________________________________
void AliEveTRDTrackListEditor::SetModel(TObject* obj)
{
  // Set model object
  fM = dynamic_cast<AliEveTRDTrackList*>(obj);

  UpdateMacroList();
  UpdateHistoList();  
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateHistoList()
{
  fHistoSubFrame->TGCompositeFrame::RemoveAll();
  
  // Set buttons for histograms
  if (fCheckButtons != 0) delete fCheckButtons;
  fCheckButtons = new TGCheckButton*[fM->fDataFromMacroList->GetEntries()];
  
  TObjString* iter = (TObjString*)fM->fDataFromMacroList->First();
  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries() && iter != 0; i++)
  {
    fCheckButtons[i] = new TGCheckButton(fHistoSubFrame, iter->GetName());
    fHistoSubFrame->AddFrame(fCheckButtons[i]);
    
    fCheckButtons[i]->SetState(kButtonUp, kFALSE);
    fCheckButtons[i]->MapRaised();
            
    iter = (TObjString*)fM->fDataFromMacroList->After(iter);
  }  
}


//______________________________________________________
void AliEveTRDTrackListEditor::UpdateMacroList()
{
  ftlMacroList->RemoveAll();
 
  TObjString* iter = (TObjString*)fM->fMacroList->First();

  while (iter != 0)
  {
    ftlMacroList->AddEntry(iter->GetName(), -1);
    iter = (TObjString*)fM->fMacroList->After(iter);
  }

  ftlMacroList->SortByName();


  ftlMacroSelList->RemoveAll();
 
  iter = (TObjString*)fM->fMacroSelList->First();

  while (iter != 0)
  {
    ftlMacroSelList->AddEntry(iter->GetName(), -1);
    iter = (TObjString*)fM->fMacroSelList->After(iter);
  }

  ftlMacroSelList->SortByName(); 
}
