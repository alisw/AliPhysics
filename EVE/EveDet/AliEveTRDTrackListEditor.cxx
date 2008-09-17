#include <EveDet/AliEveTRDData.h>
#include <EveDet/AliEveTRDTrackList.h>
#include "AliEveTRDTrackListEditor.h"

#include <EveBase/AliEveEventManager.h>
#include <AliTRDtrackV1.h>
#include <TGButton.h>
#include <TCanvas.h>     
#include <TEveBrowser.h>
#include <TEveGedEditor.h> 
#include <TEveMacro.h>
#include <TEveManager.h>
#include <TFile.h>
#include <TG3DLine.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGMsgBox.h>
#include <TGTab.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TH1.h>
#include <TTreeStream.h>


ClassImp(AliEveTRDTrackListEditor)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackListEditor //////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackListEditor::AliEveTRDTrackListEditor(const TGWindow* p, Int_t width, Int_t height,
				                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options, back),
  fM(0),
  fHistoCanvas(0),
  fHistoCanvasName(0),
  fInheritMacroList(kFALSE),
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

  fHistoCanvasName = new TGString("");

  // Handle the signal "Selected(Int_t ind)"
  ftlMacroList->Connect("Selected(Int_t)", "AliEveTRDTrackListEditor", this, "UpdateMacroListSelection(Int_t)");
  ftlMacroSelList->Connect("Selected(Int_t)", "AliEveTRDTrackListEditor", this, "UpdateMacroSelListSelection(Int_t)");

  // Handle the signal "NewEventLoaded"
  gAliEveEvent->Connect("NewEventLoaded()", "AliEveTRDTrackListEditor", this, "HandleNewEventLoaded()");

  // Handle the signal "Selected" (another tab has been selected)
  GetGedEditor()->GetTab()->Connect("Selected(Int_t)", "AliEveTRDTrackListEditor", 
                                    this, "HandleTabChangedToIndex(Int_t)");
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
  // Close and delete all tabs that have been created by this class
  CloseTabs();

  if (fHistoCanvasName != 0)
  {
    delete fHistoCanvasName;
    fHistoCanvasName = 0;
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
                 "Macro has not the signature of...\n...a selection macro: Bool_t YourMacro(const AliTRDtrackV1*)\n...a process macro (type 1): void YourMacro(const AliTRDtrackV1*, Double_t*&, Int_t&)\n...a process macro (type2): TH1* YourMacro(const AliTRDtrackV1*)", 
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

  // Update histogram tab (data has to be reloaded)
  SetModel(fM);
  Update();

  // AlieveTRDTrackList::ApplyProcessMacros() automatically selects a macro -> Draw the histogram for it,
  // if a process macro has been applied
  if (iterator->GetEntries() > 0) 
  {
    // Set focus on "Histograms" tab
    GetGedEditor()->GetTab()->SetTab("Histograms");

    DrawHistos();
  }

  if (iterator != 0)  delete iterator;  
  iterator = 0;  
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
void AliEveTRDTrackListEditor::CloseTabs()
{
  if (fHistoCanvas != 0)
  {
    // Close the created tab, if it exists
    if (fHistoCanvasName != 0)
    {
      if (gEve->GetBrowser()->GetTab(1)->SetTab(fHistoCanvasName->GetString()))
      {
        // Now the created tab is the current one and can be deleted
        gEve->GetBrowser()->GetTab(1)->RemoveTab();
      }
    }
    // With the tab removal, the canvas will be deleted automatically!
    fHistoCanvas = 0;
  }
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

  // Check, if a histo macro shall be drawn
  Int_t indexOfHistoMacro = -1;
  Int_t selectedChecked = 0;
  for (Int_t j = 0; j < fM->fDataFromMacroList->GetEntries(); j++)
  {
    if (fCheckButtons[j]->TGButton::GetState() == kButtonDown)
    {
      selectedChecked++;

      // Histo macro?
      if (fM->IsHistogramMacro(fM->fDataFromMacroList->At(j)->GetName()))
      {
        // Is also another macro selected?
        if (nHistograms > 1)
        {
          // Histo macros cannot(!) be correlated!
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       "Process macros of type 2 (return value \"TH1*\") cannot be correlated with other macros", 
                       kMBIconExclamation, kMBOk);
          return;        
        }

        // Mark this histo macro for drawing
        indexOfHistoMacro = j;

        // Have all selected macros been checked? -> If yes, we are done with this
        if (selectedChecked == nHistograms)  break;
      }
    }
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

  // Variable for the loop below -> Will be set to aborting value, if a histo macro is drawn
  Int_t i = 0;
  
  // Draw histo macro?
  if (indexOfHistoMacro >= 0)
  {
    if ((t = (TTree*)file->Get(Form("TrackData%d", indexOfHistoMacro))))
    {
      SetDrawingToHistoCanvasTab();
 
      TH1* myHist = 0;
      t->SetBranchAddress(Form("Macro%d", indexOfHistoMacro), &myHist);
      t->GetEntry(0);
      if (myHist != 0)  myHist->Draw();
      else
      {
        Error("Draw histograms", Form("No histogram for histo macro \"%s\" found!", 
                                      fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()));
        new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                     Form("No histogram for histo macro \"%s\" found!", 
                          fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()), kMBIconExclamation, kMBOk);
               
      }

      UpdateHistoCanvasTab();      
    }
    else
    {
      Error("Draw histograms", Form("No data for histo macro \"%s\" found!", 
                                    fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()));
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                   Form("No data for histo macro \"%s\" found!", 
                        fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()), kMBIconExclamation, kMBOk);
    }

    // Skip the loop below
    i = fM->fDataFromMacroList->GetEntries();
  }

  // Load the trees in succession and remember the entries -> Plot the process macros of type 1
  for ( ; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    if (fCheckButtons[i]->TGButton::GetState() == kButtonDown)
    {
      if (t == 0)
      {
        indexOfMacro1 = i;
        if (!(t = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()),
                            kMBIconExclamation, kMBOk);
          break;   
        }

        // 1d histogram - macro of type 1   
        if (nHistograms == 1) 
        {
          SetDrawingToHistoCanvasTab();

          t->Draw(Form("Macro%d", indexOfMacro1), "1");
          ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(Form("%s;%s",
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName()));

          UpdateHistoCanvasTab();

          break;     
        }
      }
      else if (tFriend1 == 0)
      {
        indexOfMacro2 = i;
        if (!(tFriend1 = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()),
                            kMBIconExclamation, kMBOk);
          break;   
        }
        
        // 2d histogram
        if (nHistograms == 2) 
        {
          SetDrawingToHistoCanvasTab();

          t->AddFriend(tFriend1);
          t->Draw(Form("Macro%d:Macro%d", indexOfMacro1, indexOfMacro2), "1");
          ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(Form("%s - %s;%s;%s",
            fM->fDataFromMacroList->At(indexOfMacro2)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro2)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName()));

          UpdateHistoCanvasTab();
 
          break;     
        }
      }    
      // 3d histogram
      else
      {
        indexOfMacro3 = i;
        if (!(tFriend2 = (TTree*)file->Get(Form("TrackData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!", fM->fDataFromMacroList->At(i)->GetName()),
                            kMBIconExclamation, kMBOk);
          break;   
        }

        SetDrawingToHistoCanvasTab();

        t->AddFriend(tFriend1);
        t->AddFriend(tFriend2);
        t->Draw(Form("Macro%d:Macro%d:Macro%d", indexOfMacro1, indexOfMacro2, indexOfMacro3), "1");
        ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(Form("%s - %s - %s;%s;%s;%s",
            fM->fDataFromMacroList->At(indexOfMacro3)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro2)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro3)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro2)->GetName(),
            fM->fDataFromMacroList->At(indexOfMacro1)->GetName()));
        
        UpdateHistoCanvasTab();
 
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
        name = new Char_t[AliEveTRDTrackList::fkMaxMacroNameLength];
        memset(name, '\0', sizeof(Char_t) * AliEveTRDTrackList::fkMaxMacroNameLength);
        sprintf(name, "%s", fteField->GetText());

        // Add path to textfield -> Path is "./" -> Use length for the name + 2
        Char_t pathname[AliEveTRDTrackList::fkMaxMacroNameLength + 2];
        memset(pathname, '\0', sizeof(Char_t) * (AliEveTRDTrackList::fkMaxMacroNameLength + 2));
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
        Char_t* path = new Char_t[AliEveTRDTrackList::fkMaxMacroPathLength];
        memset(path, '\0', sizeof(Char_t) * AliEveTRDTrackList::fkMaxMacroPathLength);
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
void AliEveTRDTrackListEditor::HandleNewEventLoaded()
{
  // Inherit the macro list for the next track list!
  fInheritMacroList = kTRUE;

  // Close the tabs
  CloseTabs();
}

//______________________________________________________
void AliEveTRDTrackListEditor::HandleTabChangedToIndex(Int_t index)
{
  fM->SetSelectedTab(index);
}

//______________________________________________________
void AliEveTRDTrackListEditor::InheritMacroList()
{
  // The old macro lists are stored in the corresponding list boxes -> add them to the track list
    
  // Selection macros
  fM->fMacroSelList->Delete();
  for (Int_t i = 0; i < ftlMacroSelList->GetNumberOfEntries(); i++)
  {
    fM->AddMacroFast(ftlMacroSelList->GetEntry(i)->GetTitle(), kTRUE);
  }

  // Process macros
  fM->fMacroList->Delete();
  for (Int_t i = 0; i < ftlMacroList->GetNumberOfEntries(); i++)
  {
    fM->AddMacroFast(ftlMacroList->GetEntry(i)->GetTitle(), kFALSE);
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

  // Selected macros are deleted from the list -> No selected entries left
  fM->fMacroListSelected = 0;
  fM->fMacroSelListSelected = 0;

  UpdateMacroList();

  if (iterator != 0)  delete iterator;
  iterator = 0;
}

//______________________________________________________
void AliEveTRDTrackListEditor::SetDrawingToHistoCanvasTab()
{
  // If the tab with the canvas has been closed, the canvas will be deleted.
  // So, if there is no tab, set the canvas pointer to zero and recreate it in a new tab.
  if (fHistoCanvas != 0) 
  {
    if (gEve->GetBrowser()->GetTab(1)->SetTab(fHistoCanvasName->GetString()) == 0)
    {
      fHistoCanvas = 0;
    }
  }

  if (!fHistoCanvas)
  {
    fHistoCanvas = gEve->AddCanvasTab(fM->GetName());     
  }
                           
  gPad = fHistoCanvas;
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateHistoCanvasTab()
{
  // Update name of the tab (tab has been set to current tab!)
  fHistoCanvasName->SetString(fM->GetName());  

  // Use a copy of fHistoCanvasName!! -> If the user closes a tab manually, the TGString
  // will be deleted -> Error might occur, when accessing the pointer   
  gEve->GetBrowser()->GetTab(1)->GetCurrentTab()->SetText(new TGString(fHistoCanvasName));

  // Switch tabs to force redrawing
  gEve->GetBrowser()->GetTab(1)->SetTab(0);
  gEve->GetBrowser()->GetTab(1)->SetTab(fHistoCanvasName->GetString());
  fHistoCanvas->Update();
}

//______________________________________________________
void AliEveTRDTrackListEditor::SetModel(TObject* obj)
{  
  // Set model object
  fM = dynamic_cast<AliEveTRDTrackList*>(obj);

  if (fM == 0) 
  {
    Error("SetModel", "Parameter is zero pointer");
    return;
  }

  // If macro list shall be inherited from previously loaded track list, do so
  if (fInheritMacroList)
  {
    InheritMacroList();
    fInheritMacroList = kFALSE;
  }
  
  UpdateMacroList();
  UpdateHistoList(); 

  // View correct tab
  GetGedEditor()->GetTab()->SetTab(fM->fSelectedTab); 
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateDataFromMacroListSelection()
{
  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    fM->SetHistoDataSelection(i, fCheckButtons[i]->IsOn());
  }
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateHistoList()
{
  fHistoSubFrame->TGCompositeFrame::Cleanup();
  
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
    fCheckButtons[i]->SetOn(fM->HistoDataIsSelected(i));
    fCheckButtons[i]->Connect("Clicked()", "AliEveTRDTrackListEditor", this, "UpdateDataFromMacroListSelection()");
            
    iter = (TObjString*)fM->fDataFromMacroList->After(iter);
  }  
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateMacroList()
{
  ftlMacroList->RemoveAll();
 
  TObjString* iter = (TObjString*)fM->fMacroList->First();

  Int_t ind = 0;
  while (iter != 0)
  {
    ftlMacroList->AddEntry(iter->GetName(), ind++);
    iter = (TObjString*)fM->fMacroList->After(iter);
  }

  ftlMacroList->SortByName();

  // Select, what has been selected before
  for (Int_t i = 0; i < fM->fMacroList->GetEntries(); i++)
  {
    ftlMacroList->Select(i, fM->MacroListIsSelected(i));
  }



  ftlMacroSelList->RemoveAll();
 
  iter = (TObjString*)fM->fMacroSelList->First();

  ind = 0;
  while (iter != 0)
  {
    ftlMacroSelList->AddEntry(iter->GetName(), ind++);
    iter = (TObjString*)fM->fMacroSelList->After(iter);
  }

  ftlMacroSelList->SortByName(); 

  // Select, what has been selected before
  for (Int_t i = 0; i < fM->fMacroSelList->GetEntries(); i++)
  {
    ftlMacroSelList->Select(i, fM->MacroSelListIsSelected(i));
  }
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateMacroListSelection(Int_t ind)
{
  // Toggle selected item
  fM->SetMacroListSelection(ind, !fM->MacroListIsSelected(ind));
}

//______________________________________________________
void AliEveTRDTrackListEditor::UpdateMacroSelListSelection(Int_t ind)
{
  // Toggle selected item
  fM->SetMacroSelListSelection(ind, !fM->MacroSelListIsSelected(ind));
}
