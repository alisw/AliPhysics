// Author: Benjamin Hess   06/11/2009

/*************************************************************************
 * Copyright (C) 2009, Alexandru Bercuci, Benjamin Hess.                 *
 * All rights reserved.                                                  *
 *************************************************************************/

//TODO: Documentation
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveListAnalyserEditor                                     //
//                                                                      //
// The AliEveListAnalyserEditor provides the graphical          //
// functionality for the AliEveTRDAnalyseObject. It creates the tabs    //
// and canvases, when they are needed and, as well, frees allocated     //
// memory on destruction (or if new events are loaded and thus some     //
// tabs are closed). The function DrawHistos() accesses the temporary   //
// file created by the AliEveListAnalyser and draws the desired //
// data (the file will be created within the call of ApplyMacros()).    //
// Have a look at this function to learn more about the structure of    //
// the file and how to access the data.                                 //
//////////////////////////////////////////////////////////////////////////

#include <EveDet/AliEveListAnalyser.h>
#include "AliEveListAnalyserEditor.h"

#include <EveBase/AliEveEventManager.h>
#include <TCanvas.h>     
#include <TEveBrowser.h>
#include <TEveGedEditor.h> 
#include <TEveMacro.h>
#include <TEveManager.h>
#include <TFile.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGMsgBox.h>
#include <TGTab.h>
#include <TMap.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TGTextEntry.h>
#include <TGTextEdit.h>
#include <TGComboBox.h>
#include <TGTextView.h>
#include <TH1.h>
#include <TTreeStream.h>


ClassImp(AliEveListAnalyserEditor)

///////////////////////////////////////////////////////////
/////////////   AliEveListAnalyserEditor //////////
///////////////////////////////////////////////////////////
AliEveListAnalyserEditor::AliEveListAnalyserEditor(const TGWindow* p, Int_t width, Int_t height,
				                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options, back),
  fM(0),
  fHistoCanvas(0),
  fHistoCanvasName(0),
  fInheritedMacroList(0),
  fInheritSettings(kFALSE),
  fMainFrame(0),
  fHistoFrame(0),
  fHistoSubFrame(0),
  fBrowseFrame(0),
  fbBrowse(0),
  fbNew(0),
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
  // Creates the AliEveListAnalyserEditor.

  // Functionality for adding macros  
  fMainFrame = CreateEditorTabSubFrame("Process");
   
  fLabel1 = new TGLabel(fMainFrame,"Add macro(s):");
  fMainFrame->AddFrame(fLabel1);
  fBrowseFrame = new TGHorizontalFrame(fMainFrame);

  fteField = new TGTextEntry(fBrowseFrame);
  fteField->SetToolTipText("Enter the pathname of the macro you want to add here and press \"Enter\"");
  fteField->Connect("ReturnPressed()","AliEveListAnalyserEditor", this, "HandleMacroPathSet()"); 
  fBrowseFrame->AddFrame(fteField);
  
  fbBrowse = new TGTextButton(fBrowseFrame, "Browse");
  fbBrowse->SetToolTipText("Browse the macro you want to add");
  fbBrowse->Connect("Clicked()", "AliEveListAnalyserEditor", this, "BrowseMacros()");
  fBrowseFrame->AddFrame(fbBrowse);
  
  fbNew = new TGTextButton(fBrowseFrame, "New");
  fbNew->SetToolTipText("Start macro creation wizard");
  fbNew->Connect("Clicked()", "AliEveListAnalyserEditor", this, "NewMacros()");
  fBrowseFrame->AddFrame(fbNew);
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
  fLabel3 = new TGLabel(fMainFrame,"Process plugins:");
  fMainFrame->AddFrame(fLabel3);

  ftlMacroList = new TGListBox(fMainFrame);
  ftlMacroList->Resize(194, 94);
  ftlMacroList->SetMultipleSelections(kTRUE);
  fMainFrame->AddFrame(ftlMacroList);

  fLine3 = new TGHorizontal3DLine(fMainFrame, 194, 8);
  fMainFrame->AddFrame(fLine3, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));  

  fbApplyMacros = new TGTextButton(fMainFrame, "Apply plugin(s)");
  fbApplyMacros->SetToolTipText("Apply all selected macros/class functins to the list of objects -> A data file will be generated");
  fbApplyMacros->Connect("Clicked()", "AliEveListAnalyserEditor", this, "ApplyMacros()");
  fbApplyMacros->SetRightMargin(12);
  fMainFrame->AddFrame(fbApplyMacros);

  fbRemoveMacros = new TGTextButton(fMainFrame, "Remove plugin(s)");
  fbRemoveMacros->SetToolTipText("Remove the selected macros/class functions from the list(s)");
  fbRemoveMacros->Connect("Clicked()", "AliEveListAnalyserEditor", this, "RemoveMacros()");
  fMainFrame->AddFrame(fbRemoveMacros);

  // Stuff for displaying histograms
  fHistoFrame = CreateEditorTabSubFrame("Results");  
  fHistoFrame->SetMapSubwindows(kTRUE);
  fLabel4 = new TGLabel(fHistoFrame,"Data from plugins:");
  fHistoFrame->AddFrame(fLabel4);

  fHistoSubFrame = new TGVerticalFrame(fHistoFrame);
  fHistoSubFrame->SetMapSubwindows(kTRUE);
  fHistoSubFrame->Resize(194, 200);
  fHistoFrame->AddFrame(fHistoSubFrame);

  fLine4 = new TGHorizontal3DLine(fHistoFrame, 194, 8);
  fHistoFrame->AddFrame(fLine4, new TGLayoutHints(kLHintsLeft  | kLHintsTop, 2, 2, 8, 2));  

  fbDrawHisto = new TGTextButton(fHistoFrame, "Draw projections");
  fbDrawHisto->SetToolTipText("Uses the data file created by the last \"Apply selected macro(s)\".\nClick here to display the data histograms of the selected macros.\nSelect multiple macros to create multi-dimensional plots.\nHisto macros cannot be used for multi-dimensional plots!");
  fbDrawHisto->Connect("Clicked()", "AliEveListAnalyserEditor", this, "DrawHistos()");
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
  ftlMacroList->Connect("Selected(Int_t)", "AliEveListAnalyserEditor", this, "UpdateMacroListSelection(Int_t)");
  ftlMacroSelList->Connect("Selected(Int_t)", "AliEveListAnalyserEditor", this, "UpdateMacroListSelection(Int_t)");

  // Handle the signal "NewEventLoaded"
  AliEveEventManager::GetMaster()->Connect("NewEventLoaded()", "AliEveListAnalyserEditor", this, "HandleNewEventLoaded()");

  // Handle the signal "Selected" (another tab has been selected)
  GetGedEditor()->GetTab()->Connect("Selected(Int_t)", "AliEveListAnalyserEditor", this, "HandleTabChangedToIndex(Int_t)");
}

//______________________________________________________
AliEveListAnalyserEditor::~AliEveListAnalyserEditor()
{
  // Destructor: Closes all tabs created by this object and
  // frees the corresponding memory.

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
  
  if (fInheritedMacroList != 0)
  {
    fInheritedMacroList->Delete();
    delete fInheritedMacroList;
    fInheritedMacroList = 0;
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::AddMacro(const Char_t* name, const Char_t* path)
{
  // Adds the macro path/name to the macro list. A warning is provided, if there is
  // something wrong, e.g. if the macro does not have the correct signature.
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
  case ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", "Fail to load the macro (check messages in the terminal)!",
                 kMBIconExclamation, kMBOk);
    break;
  case SIGNATURE_ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Macro has not the signature of...\n...a single object selection macro: Bool_t YourMacro(const YourObjectType*)\n...a correlated objects selection macro: Bool_t YourMacro(const YourObjectType*, const YourObjectType*)\n...a single object analyse macro: void YourMacro(const YourObjectType*, Double_t*&, Int_t&)\n...a correlated objects analyse macro: void YourMacro(const YourObjectType*, const YourObjectType*, Double_t*&, Int_t&)\n...a single object histo macro: TH1* YourMacro(const YourObjectType*)\n...a correlated objects histo macro: TH1* YourMacro(const YourObjectType*, const YourObjectType*)", 
                 kMBIconExclamation, kMBOk);
    break;               
  case NOT_EXIST_ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "File does not exist or you do not have read permission!", kMBIconExclamation, kMBOk);
    break;
  case UNKNOWN_OBJECT_TYPE_ERROR:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Unknown object type of macro parameter!", kMBIconExclamation, kMBOk);
    break;
  default:
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 Form("AliEveListAnalyser::AddMacro exited with unknown return value: %d", result),
                 kMBIconExclamation, kMBOk);
    break;
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::ApplyMacros()
{
  // Applies the selected macros and updates the view.

  Bool_t success = kFALSE;

  // First apply the single object selection macros
  TList* selIterator = new TList();
  ftlMacroSelList->GetSelectedEntries(selIterator);
  fM->ApplySTSelectionMacros(selIterator);
  
  // Update view
  gEve->Redraw3D();

  // Now apply the process macros
  TList* procIterator = new TList();
  ftlMacroList->GetSelectedEntries(procIterator);
  success = fM->ApplyProcessMacros(selIterator, procIterator);

  // Update histogram tab (data has to be reloaded)
  SetModel(fM);
  Update();

  // AliEveListAnalyser::ApplyProcessMacros() automatically selects a macro -> Draw the histogram for it,
  // if a process macro has been applied
  if (success && procIterator->GetEntries() > 0) 
  {
    // Set focus on "Histograms" tab
    GetGedEditor()->GetTab()->SetTab("Results");

    DrawHistos();
  }

  if (selIterator != 0) delete selIterator;
  selIterator = 0;  
  if (procIterator != 0)  delete procIterator;  
  procIterator = 0;  
  
  if (!success)
  {
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "AliEveListAnalyser::ApplyProcessMacros experienced an error (cf. CINT-output)!", 
                 kMBIconExclamation, kMBOk);  
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::BrowseMacros()
{
  // Creates a file-dialog. The selected files will be added to the macro list
  // via AddMacro(...).

  new TGFileDialog(gClient->GetRoot(), GetMainFrame(), kFDOpen, fFileInfo);
  
  if (fFileInfo->fIniDir != 0 && fFileInfo->fFileNamesList != 0)
  {       
    // Extract filenames
    TObject* iter = fFileInfo->fFileNamesList->First();
 
    Char_t* name = 0;

    while (iter != 0)
    {
      // NOTE: fileInfo->fFileNamesList will be changed by that, too!
      name = (Char_t*)strrchr(iter->GetName(), '/');
      // Delete '"' at the end
      name[strlen(name)] = '\0';
              
      AddMacro(name + 1, fFileInfo->fIniDir); 
      iter = (TObjString*)fFileInfo->fFileNamesList->After(iter);
    }
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::CloseTabs()
{
  // Closes + deletes the tabs created by this object

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
void AliEveListAnalyserEditor::DrawHistos()
{
  // Accesses the temporary data file created by the last call of ApplyMacros() and draws
  // histograms according to the selection in the "Histograms"-tab.
 
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

      // Histo macro? -> To check this, look for the substring "(histo macro)"
      if (strstr(fM->fDataFromMacroList->At(j)->GetName(), "(histo macro)") != 0)
      {
        // Is also another macro selected?
        if (nHistograms > 1)
        {
          // Histo macros cannot(!) be correlated!
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       "Histo macros (return value \"TH1*\") cannot be combined with other macros", 
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

  TFile* file = new TFile(Form("/tmp/ListAnalyserMacroData_%s.root", gSystem->Getenv("USER")), "READ");
  if (!file)  
  {
    Error("Draw histograms", Form("Cannot open file \"/tmp/ListAnalyserMacroData_%s.root\"", 
                                  gSystem->Getenv("USER")));
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                 Form("Cannot open file \"/tmp/ListAnalyserMacroData_%s.root\"", gSystem->Getenv("USER")),
                 kMBIconExclamation, kMBOk);
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
    if ((t = (TTree*)file->Get(Form("ObjectData%d", indexOfHistoMacro))))
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
      Error("Draw histograms", Form("No data for histo macro \"%s\" found!\nMaybe no objects have been selected.", 
                                    fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()));
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                   Form("No data for histo macro \"%s\" found!\nMaybe no objects have been selected.", 
                        fM->fDataFromMacroList->At(indexOfHistoMacro)->GetName()), kMBIconExclamation, kMBOk);
    }

    // Skip the loop below
    i = fM->fDataFromMacroList->GetEntries();
  }

  // Load the trees in succession and remember the entries -> Plot the analyse macros
  for ( ; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    if (fCheckButtons[i]->TGButton::GetState() == kButtonDown)
    {
      if (t == 0)
      {
        indexOfMacro1 = i;
        if (!(t = (TTree*)file->Get(Form("ObjectData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                                        fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                            fM->fDataFromMacroList->At(i)->GetName()), kMBIconExclamation, kMBOk);
          break;   
        }

        // 1d histogram
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
        if (!(tFriend1 = (TTree*)file->Get(Form("ObjectData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                                        fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                            fM->fDataFromMacroList->At(i)->GetName()),
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
        if (!(tFriend2 = (TTree*)file->Get(Form("ObjectData%d", i))))
        { 
          Error("Draw histograms", Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                                        fM->fDataFromMacroList->At(i)->GetName()));
          new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error - Draw histograms", 
                       Form("No data for macro \"%s\" found!\nMaybe no objects have been selected.", 
                            fM->fDataFromMacroList->At(i)->GetName()), kMBIconExclamation, kMBOk);
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
Int_t AliEveListAnalyserEditor::GetNSelectedHistograms() const
{
  // Returns the number of selected macros (or rather: of their selected data) in the "Histograms"-tab

  Int_t count = 0;
  
  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    if (fCheckButtons[i]->TGButton::GetState() == kButtonDown)  count++;
  }

  return count;
}

//______________________________________________________
void AliEveListAnalyserEditor::HandleMacroPathSet()
{
  // Takes the input of the text field (adding a macro), checks if the macro can be
  // accessed (and that it exists) and adds the macro to the macro list via AddMacro(...).
  // You can use environment variables in the text field, e.g. "$ALICE_ROOT/Eve/alice-macro/myMacro.C".

  if (strlen(fteField->GetText()) != 0)
  {  
    // Expand the pathname
    Char_t* systemPath = gSystem->ExpandPathName(fteField->GetText());
    fteField->SetText(systemPath);
    delete systemPath;
    systemPath = 0;
       			
    // Check if file exists
    FILE* fp = NULL;

    fp = fopen(fteField->GetText(), "rb");
    if (fp != NULL)
    {
      fclose(fp);

      // Extract filename
      Char_t* name = (Char_t*)strrchr(fteField->GetText(), '/');

      // Current path
      if (name == NULL)
      {
        name = new Char_t[AliEveListAnalyser::fkMaxMacroNameLength];
        memset(name, '\0', sizeof(Char_t) * AliEveListAnalyser::fkMaxMacroNameLength);
        sprintf(name, "%s", fteField->GetText());

        // Add path to textfield -> Path is "./" -> Use length for the name + 2
        Char_t pathname[AliEveListAnalyser::fkMaxMacroNameLength + 2];
        memset(pathname, '\0', sizeof(Char_t) * (AliEveListAnalyser::fkMaxMacroNameLength + 2));
        sprintf(pathname, "./%s", fteField->GetText());
        fteField->SetText(pathname);

        AddMacro(name);  
        if (name != 0)  delete name;
        name = 0;
      }
      // Different path
      else
      {
        // Extract path
        Char_t* path = new Char_t[AliEveListAnalyser::fkMaxMacroPathLength];
        memset(path, '\0', sizeof(Char_t) * AliEveListAnalyser::fkMaxMacroPathLength);
        strncpy(path, fteField->GetText(), strlen(fteField->GetText()) - strlen(name));
        
        // Ignore the slash "/" in name
        AddMacro(name + 1, path);  
  
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
void AliEveListAnalyserEditor::HandleNewEventLoaded()
{
  // Closes the tabs created by this object and sets a flag that will
  // cause the function SetModel() to inherit the macro lists
  // for the next AliEveListAnalyser from the current one.

  // Inherit the macro list for the next analyse object list!
  fInheritSettings = kTRUE;

  // Close the tabs
  CloseTabs();
}

//______________________________________________________
void AliEveListAnalyserEditor::HandleTabChangedToIndex(Int_t index)
{
  // Saves the current tab in the current AliEveListAnalyser.

  fM->SetSelectedTab(index);
}

//______________________________________________________
void AliEveListAnalyserEditor::InheritMacroList()
{
  // The old macro list is possibly stored in the corresponding interior map. This function will 
  // use this interior map to move the data from the interior map to the newly loaded AliEveListAnalyser. 
  // Then the interior map will be cleaned up. With this, the settings will be inherited from the previously 
  // loaded AliEveListAnalyser.

  if (fInheritedMacroList == 0)  return;

  // Clear list  
  fM->fMacroList->Delete();

  // Store data from interior list in the analyse object list's map
  TMapIter* iter = (TMapIter*)fInheritedMacroList->MakeIterator();
  
  TObject* key = 0;
  TGeneralMacroData* macro = 0;
  
  while ((key = iter->Next()) != 0)
  {
    macro = (TGeneralMacroData*)fInheritedMacroList->GetValue(key);
    if (macro != 0)  fM->fMacroList->Add(new TObjString(key->GetName()), 
                                         new TGeneralMacroData(macro->GetName(), macro->GetPath(), macro->GetType()));
    else
    {
      Error("AliEveListAnalyserEditor::InheritMacroList", Form("Failed to inherit the macro \"%s\"!", key));
    }
  }
  
  fInheritedMacroList->Delete();
  delete fInheritedMacroList;
  fInheritedMacroList = 0;
}

//______________________________________________________
void AliEveListAnalyserEditor::NewMacros()
{
  // Start the macro creation wizard.
  // thanks to Jacek Otwinowski<J.Otwinowski@GSI.DE> for this suggestion

  AliEveGeneralMacroWizard *wizz = new AliEveGeneralMacroWizard();
  wizz->Connect("Create(Char_t*)", "AliEveListAnalyserEditor", this, "AddMacro(Char_t*)");
}

//______________________________________________________
void AliEveListAnalyserEditor::RemoveMacros()
{
  // Removes the selected macros from the corresponding list.

  TList* iterator = new TList();
  
  ftlMacroList->GetSelectedEntries(iterator);
  fM->RemoveSelectedMacros(iterator);

  if (iterator != 0)  delete iterator;

  iterator = new TList();
  ftlMacroSelList->GetSelectedEntries(iterator);
  fM->RemoveSelectedMacros(iterator);

  // Selected macros are deleted from the list -> No selected entries left
  fM->fMacroListSelected = 0;

  UpdateMacroList();

  if (iterator != 0)  delete iterator;
  iterator = 0;
}

//______________________________________________________
void AliEveListAnalyserEditor::SaveMacroList(TMap* list)
{
  // Saves the provided macro list in an interior list. This list will be used by
  // InheritMacroList() to restore the data in "list". With this method one is able
  // to inherit the macro list from analyse object list to analyse object list (i.e. from event to event).

  if (fInheritedMacroList != 0)
  {
    fInheritedMacroList->Delete();
    delete fInheritedMacroList;
  }
  fInheritedMacroList = new TMap();
  fInheritedMacroList->SetOwnerKeyValue(kTRUE, kTRUE);

  TMapIter* iter = (TMapIter*)list->MakeIterator();
  TObject* key = 0;
  TGeneralMacroData* macro = 0;
  
  while ((key = iter->Next()) != 0)
  {
    macro = (TGeneralMacroData*)fM->fMacroList->GetValue(key);
    if (macro != 0) fInheritedMacroList->Add(new TObjString(key->GetName()), 
                                             new TGeneralMacroData(macro->GetName(), macro->GetPath(), macro->GetType()));
    else
    {
      Error("AliEveListAnalyserEditor::SaveMacroList", Form("Failed to inherit the macro \"%s\"!", key));
    }
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::SetDrawingToHistoCanvasTab()
{
  // Sets gPad to the tab with the name of the current AliEveListAnalyser. If this tab does
  // not exist, it will be created. Otherwise, it is re-used.

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
void AliEveListAnalyserEditor::SetModel(TObject* obj)
{  
  // Sets the model object, updates the related data in the GUI and
  // inherits settings (cf. Inherit*(...)), if the flag fInheritSettings is set to kTRUE.

  fM = dynamic_cast<AliEveListAnalyser*>(obj);

  if (fM == 0) 
  {
    Error("SetModel", "Parameter is zero pointer");
    return;
  }

  // Provide a pointer to this editor
  fM->fEditor = this;

  // If macro list + track style shall be inherited from previously loaded track list, do so
  if (fInheritSettings)
  {
    InheritMacroList();

    fInheritSettings = kFALSE;
  }

  UpdateMacroList();
  UpdateHistoList(); 

  // View correct tab
  GetGedEditor()->GetTab()->SetTab(fM->GetSelectedTab()); 
}

//______________________________________________________
void AliEveListAnalyserEditor::UpdateDataFromMacroListSelection()
{
  // Saves the current selection in the "Histograms"-tab to the current
  // AliEveListAnalyser. This means that the selection is updated and won't
  // get lost, if another editor is loaded in Eve.

  for (Int_t i = 0; i < fM->fDataFromMacroList->GetEntries(); i++)
  {
    fM->SetHistoDataSelection(i, fCheckButtons[i]->IsOn());
  }
}

//______________________________________________________
void AliEveListAnalyserEditor::UpdateHistoCanvasTab()
{
   // Updates the histogram and the corresponding tab (including titles).

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
void AliEveListAnalyserEditor::UpdateHistoList()
{
  // Reloads (updates) the buttons in the "Histograms"-tab via
  // the current AliEveListAnalyser (data).

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
    fCheckButtons[i]->Connect("Clicked()", "AliEveListAnalyserEditor", this, "UpdateDataFromMacroListSelection()");
            
    iter = (TObjString*)fM->fDataFromMacroList->After(iter);
  }  
}

//______________________________________________________
void AliEveListAnalyserEditor::UpdateMacroList()
{
  // Reloads (updates) the macro list (selection AND process macros) via
  // the current AliEveListAnalyser (data).

  ftlMacroList->RemoveAll();
  ftlMacroSelList->RemoveAll();
   
  TMapIter* iter = (TMapIter*)fM->fMacroList->MakeIterator();
  TObject* key = 0;
  TGeneralMacroData* macro = 0;

  Int_t ind = 0;
  while ((key = iter->Next()) != 0)
  {
    macro = (TGeneralMacroData*)fM->fMacroList->GetValue(key);
    if (macro != 0)
    {
      if (macro->IsProcessMacro())
      {
        ftlMacroList->AddEntry(macro->GetName(), ind);
        // Select, what has been selected before
        ftlMacroList->Select(ind, fM->MacroListIsSelected(ind));
        ind++;
      }
      else if (macro->IsSelectionMacro())
      {
        ftlMacroSelList->AddEntry(macro->GetName(), ind);
        // Select, what has been selected before
        ftlMacroSelList->Select(ind, fM->MacroListIsSelected(ind));
        ind++;
      }
      else
      {
        Error("AliEveListAnalyserEditor::UpdateMacroList()", 
              Form("Macro \"%s/%s.C\" is neither a selection macro nor a process macro!",
                   macro->GetPath(), macro->GetName()));                                        
      }
    }
    else
    {
      Error("AliEveListAnalyserEditor::UpdateMacroList()", 
              Form("Macro list is corrupted: Macro \"%s\" not found!", key->GetName()));            
    }     
  }

  ftlMacroList->SortByName(); 
  ftlMacroSelList->SortByName(); 
}

//______________________________________________________
void AliEveListAnalyserEditor::UpdateMacroListSelection(Int_t ind)
{
  // Saves the current selection in the macro listS to the current
  // AliEveListAnalyser. This means that the selection is updated and won't
  // get lost, if another editor is loaded in Eve.
  // NOTE: The indices in BOTH lists will be unique!

  // Toggle selected item
  fM->SetMacroListSelection(ind, !fM->MacroListIsSelected(ind));
}


//______________________________________________________
//______________________________________________________
//______________________________________________________


/////////////////////////////////////////////////
ClassImp(AliEveGeneralMacroWizard)
/////////////////////////////////////////////////

//______________________________________________________
AliEveGeneralMacroWizard::AliEveGeneralMacroWizard(const TGWindow* p)
  :TGMainFrame(p ? p : gClient->GetRoot(), 10, 10, kMainFrame | kVerticalFrame)
  ,fTextName(0x0)
  ,fTextObjectType(0x0)
  ,fCombo(0x0)
  ,fTextEdit(0x0)
  ,fbCreate(0x0)
  ,fbCancel(0x0)
{
  const Int_t width = 300;

  // horizontal frame
  TGHorizontalFrame *fFrameName = new TGHorizontalFrame(this, 10, 10, kHorizontalFrame);
  TGLabel *fLabel = new TGLabel(fFrameName, "Name*");
  fLabel->SetTextJustify(36);
  fLabel->SetMargins(0,0,0,0);
  fLabel->SetWrapLength(-1);
  fFrameName->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  fTextName = new TGTextEntry(fFrameName);
  fTextName->SetMaxLength(255);
  fTextName->SetAlignment(kTextLeft);
  fTextName->SetText("");
  fTextName->SetToolTipText("The name of your macro");
  fTextName->Resize(width, fTextName->GetDefaultHeight());
  fFrameName->AddFrame(fTextName, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2));

  // horizontal frame
  TGHorizontalFrame *fFrameObjectType = new TGHorizontalFrame(this, 10, 10, kHorizontalFrame);
  fLabel = new TGLabel(fFrameObjectType, "Object type of macro");
  fLabel->SetTextJustify(36);
  fLabel->SetMargins(0,0,0,0);
  fLabel->SetWrapLength(-1);
  fFrameObjectType->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  fTextObjectType = new TGTextEntry(fFrameObjectType);
  fTextObjectType->SetMaxLength(255);
  fTextObjectType->SetAlignment(kTextLeft);
  fTextObjectType->SetText("");
  // Limit max.length to 80 characters
  fTextObjectType->SetMaxLength(80);
  fTextObjectType->SetToolTipText("The type of objects, your macro will work with");
  fTextObjectType->Resize(width, fTextObjectType->GetDefaultHeight());
  fFrameObjectType->AddFrame(fTextObjectType, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2));

  // horizontal frame
  TGHorizontalFrame *fFrameComment = new TGHorizontalFrame(this,10,10,kHorizontalFrame);
  fLabel = new TGLabel(fFrameComment, "Comment");
  fLabel->SetTextJustify(36);
  fLabel->SetMargins(0,0,0,0);
  fLabel->SetWrapLength(-1);
  fFrameComment->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  fTextEdit = new TGTextEdit(fFrameComment, width, 5*fTextName->GetDefaultHeight());
  fFrameComment->AddFrame(fTextEdit, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2));

  // horizontal frame
  TGHorizontalFrame *fFrameType = new TGHorizontalFrame(this,10,10,kHorizontalFrame);
  fLabel = new TGLabel(fFrameType, "Type*");
  fLabel->SetTextJustify(36);
  fLabel->SetMargins(0,0,0,0);
  fLabel->SetWrapLength(-1);
  fFrameType->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  fCombo = new TGComboBox(fFrameType, -1, kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
  fCombo->AddEntry("Single Object Selection", AliEveListAnalyser::kSingleObjectSelect);
  fCombo->AddEntry("Pair Objects Selection", AliEveListAnalyser::kCorrelObjectSelect);
  fCombo->AddEntry("Single Object Analyse", AliEveListAnalyser::kSingleObjectAnalyse);
  fCombo->AddEntry("Single Object Histo", AliEveListAnalyser::kSingleObjectHisto);
  fCombo->AddEntry("Pair Objects Analyse", AliEveListAnalyser::kCorrelObjectAnalyse);
  fCombo->AddEntry("Pair Objects Histo", AliEveListAnalyser::kCorrelObjectHisto);
  fCombo->Select(-1);
  fCombo->Resize(width, fTextName->GetDefaultHeight());
  fFrameType->AddFrame(fCombo, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2));

  // horizontal frame
  TGHorizontalFrame *fFrameAction = new TGHorizontalFrame(this,10,10,kHorizontalFrame);
  fbCancel = new TGTextButton(fFrameAction, "Cancel");
  fbCancel->SetToolTipText("Exit macro creation wizard");
  fFrameAction->AddFrame(fbCancel, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2)); 
  fbCreate = new TGTextButton(fFrameAction, "Done");
  fbCreate->SetToolTipText("Use settings to create the macro");
  fFrameAction->AddFrame(fbCreate, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,2)); 


  // horizontal frame
  TGHorizontalFrame *fFrameText = new TGHorizontalFrame(this,10,10,kHorizontalFrame);
  fLabel = new TGLabel(fFrameText, "(*) Mandatory fields");
  fLabel->SetTextJustify(36);
  fLabel->SetMargins(0,0,0,0);
  fLabel->SetWrapLength(-1);
  fFrameText->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  // put things together  
  AddFrame(fFrameName, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));
  AddFrame(fFrameObjectType, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));
  AddFrame(fFrameComment, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));
  AddFrame(fFrameType, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));
  AddFrame(fFrameAction, new TGLayoutHints(kLHintsRight | kLHintsTop | kLHintsExpandX,2,2,2,2));

  TGHorizontal3DLine *fLine = new TGHorizontal3DLine(this, 281, 2);
  AddFrame(fLine, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));
  AddFrame(fFrameText, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));


  SetWindowName("Macro Wizard");
  SetMWMHints(kMWMDecorAll,
              kMWMFuncAll,
              kMWMInputModeless);
  MapSubwindows();

  Resize(GetDefaultSize());
  MapWindow();

  // Do the linking
  //fCombo->Connect("Selected(Int_t)", "AliEveGeneralMacroWizard", this, "Create(Int_t)");
  fbCreate->Connect("Clicked()", "AliEveGeneralMacroWizard", this, "HandleCreate()");
  fbCancel->Connect("Clicked()", "AliEveGeneralMacroWizard", this, "CloseWindow()");

  // Standard choice
  fCombo->Select(1, kFALSE);
}  

const Char_t *fGeneralIncludes = 
// TODO: Remove include files corresponding to a track list
"#if !defined(__CINT__) || defined(__MAKECINT__)\n"
"#include <TROOT.h>\n"
"#include <TH1.h>\n"
"#include <TRD/AliTRDgeometry.h>\n"
"#include <TRD/AliTRDcluster.h>\n"
"#include <TRD/AliTRDseedV1.h>\n"
"#include <TRD/AliTRDtrackV1.h>\n"
"#endif\n";

const Char_t *fGeneralMacroTemplate[7] = {
""
,"  if (!object) return kFALSE;\n"

,"  n = 0;\n"
"  r = 0x0;\n"
"  if (!object) return;\n"

,"  if (!object) return 0x0;\n"
"  TH1* h = 0x0;\n\n"
"// Set bins, xmin and xmax here\n"
"  Int_t n = 1;\n"
"  Double_t xmin = 0;\n"
"  Double_t xmax = 100;\n\n" 
"  if(!(h = (TH1*)gROOT->FindObject(\"h\"))){\n"
"    h = new TH1(\"h\", \"Title\", n, xmin, xmax);\n"
"    h->GetXaxis()->SetTitle("");\n"
"    h->GetYaxis()->SetTitle("");\n"
"  } else h->Reset();\n"

,"  if (!object) return kFALSE;\n"
"  if (!object2) return kFALSE;\n"

,"  n = 0;\n"
"  r = 0x0;\n"
"  if (!object) return;\n"
"  if (!object2) return;\n"

,"  if (!object) return 0x0;\n"
"  if (!object2) return 0x0;\n"
"  TH1* h = 0x0;\n\n"
"// Set bins, xmin and xmax here\n"
"  Int_t n = 1;\n"
"  Double_t xmin = 0;\n"
"  Double_t xmax = 100;\n\n" 
"  if(!(h = (TH1*)gROOT->FindObject(\"h\"))){\n"
"    h = new TH1(\"h\", \"Title\", n, xmin, xmax);\n"
"    h->GetXaxis()->SetTitle("");\n"
"    h->GetYaxis()->SetTitle("");\n"
"  } else h->Reset();\n"
};


//TODO: Needed?
/*
const Char_t *fGeneralMacroTemplate_WithType[7] = {
""
,"  if (!object) return kFALSE;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return kFALSE;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n" 
"  if (!myObject) return kFALSE;\n"

,"  n = 0;\n"
"  r = 0x0;\n"
"  if (!object) return;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n" 
"  if (!myObject) return;\n"

,"  if (!object) return 0x0;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return 0x0;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n" 
"  if (!myObject) return 0x0;\n\n"
"  TH1* h = 0x0;\n\n"
"// Set bins, xmin and xmax here\n"
"  Int_t n = 1;\n"
"  Double_t xmin = 0;\n"
"  Double_t xmax = 100;\n\n" 
"  if(!(h = (TH1*)gROOT->FindObject(\"h\"))){\n"
"    h = new TH1(\"h\", \"Title\", n, xmin, xmax);\n"
"    h->GetXaxis()->SetTitle("");\n"
"    h->GetYaxis()->SetTitle("");\n"
"  } else h->Reset();\n"

,"  if (!object) return kFALSE;\n"
"  if (!object2) return kFALSE;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return kFALSE;\n"
"  if (object2->IsA() != OBJECTTYPE::Class()) return kFALSE;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n"
"  const OBJECTTYPE* myObject2 = dynamic_cast<const OBJECTTYPE*>(object2);\n" 
"  if (!myObject) return kFALSE;\n"
"  if (!myObject2) return kFALSE;\n"

,"  n = 0;\n"
"  r = 0x0;\n"
"  if (!object) return;\n"
"  if (!object2) return;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return;\n"
"  if (object2->IsA() != OBJECTTYPE::Class()) return;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n"
"  const OBJECTTYPE* myObject2 = dynamic_cast<const OBJECTTYPE*>(object2);\n" 
"  if (!myObject) return;\n"
"  if (!myObject2) return;\n"

,"  if (!object) return 0x0;\n"
"  if (!object2) return 0x0;\n"
"  if (object->IsA() != OBJECTTYPE::Class()) return 0x0;\n"
"  if (object2->IsA() != OBJECTTYPE::Class()) return 0x0;\n\n"
"  const OBJECTTYPE* myObject = dynamic_cast<const OBJECTTYPE*>(object);\n"
"  const OBJECTTYPE* myObject2 = dynamic_cast<const OBJECTTYPE*>(object2);\n" 
"  if (!myObject) return 0x0;\n"
"  if (!myObject2) return 0x0;\n"
"  TH1* h = 0x0;\n\n"
"// Set bins, xmin and xmax here\n"
"  Int_t n = 1;\n"
"  Double_t xmin = 0;\n"
"  Double_t xmax = 100;\n\n" 
"  if(!(h = (TH1*)gROOT->FindObject(\"h\"))){\n"
"    h = new TH1(\"h\", \"Title\", n, xmin, xmax);\n"
"    h->GetXaxis()->SetTitle("");\n"
"    h->GetYaxis()->SetTitle("");\n"
"  } else h->Reset();\n"
};
*/

//______________________________________________________
void AliEveGeneralMacroWizard::Create(Int_t type)
{
  const Char_t* name = fTextName->GetText();
  if(strcmp(name,"") == 0)
  {
    Error("AliEveGeneralMacroWizard::Create", "Please specify a name for your macro.");
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Please specify a name for your macro.", kMBIconExclamation, kMBOk);
    //fCombo->Select(-1);
    return;
  }

  Bool_t useGivenType = kFALSE;

  // Remove white-spaces
  TString* typeStr = new TString();

  typeStr->Append(fTextObjectType->GetText());
  typeStr->ReplaceAll(" ", "");
  fTextObjectType->SetText(typeStr->Data(), kFALSE);

  // If an object type is provided by the user, use it!
  if (strlen(typeStr->Data()) > 0)
  {
    // Check, if the class really exists
    if (TClass::GetClass(typeStr->Data()) != 0x0)
    {
      useGivenType = kTRUE; 
    }
    else
    {
      Int_t buttonsPressed = 0;
      new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Unknown object type", 
        Form("The class of your object, \"%s\" has not been found. Do you really want to create your macro with this object type?", 
             typeStr->Data()), kMBIconExclamation, kMBYes | kMBNo, &buttonsPressed);

      if (buttonsPressed & kMBYes)  useGivenType = kTRUE;
      else useGivenType = kFALSE;

      // Cancel creation
      if (!useGivenType)
      {
        typeStr->Clear();
        if (typeStr != 0) delete typeStr;
        typeStr = 0;
        return;
      }
    }
  }

  // Note: gSystem->AccessPathName(...) returns kTRUE, if the access FAILED!
  if(!gSystem->AccessPathName(Form("./%s.C", name)))
  {
    // If there is already a file with this name -> Error
    Error("AliEveGeneralMacroWizard::Create", Form("A macro \"%s.C\" already exists in the current directory!\nPlease choose another name!", name));
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 Form("A macro \"%s.C\" already exists in the current directory!\nPlease choose another name!", name), kMBIconExclamation, kMBOk);
    //fCombo->Select(-1);
    return;
  }

  FILE* fp = 0x0;
  if(!(fp = fopen(Form("%s.C", name), "wt"))){
    Error("AliEveGeneralMacroWizard::Create", "Couldn't create macro file.");
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 "Couldn't create macro file.", kMBIconExclamation, kMBOk);
    //fCombo->Select(-1);
    return;
  }

  TGText* comment = fTextEdit->GetText();
  Char_t* line = 0x0; Int_t iline = 0;
  while((line = comment->GetLine(TGLongPosition(0,iline++), 200))) fprintf(fp, "// %s\n", line);

  fprintf(fp, "\n%s\n", fGeneralIncludes);

  TString* tempStr = new TString();

  // Use default type
  if (!useGivenType)
  {
    typeStr->Clear();
    typeStr->Append("TObject");
  }

  switch(type){
  case AliEveListAnalyser::kSingleObjectSelect:
    // Use "Bool_t 'NAME'(const 'OBJECTTYPE' *object)\n"
    tempStr->Append("Bool_t ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object)\n");
    fprintf(fp, tempStr->Data());
    break;
  case AliEveListAnalyser::kCorrelObjectSelect:
    // Use "Bool_t 'NAME'(const 'OBJECTTYPE' *object, const 'OBJECTTYPE' *object2)\n"
    tempStr->Append("Bool_t ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object, const ").Append(typeStr->Data()).Append(" *object2)\n");
    fprintf(fp, tempStr->Data());
    break;
  case AliEveListAnalyser::kSingleObjectAnalyse:    
    // Use "void 'NAME'(const 'OBJECTTYPE' *object, Double_t*& r, Int_t& n)\n"
    tempStr->Append("void ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object, Double_t*& r, Int_t& n)\n");
    fprintf(fp, tempStr->Data());
    break;
  case AliEveListAnalyser::kSingleObjectHisto:
    // Use "TH1* 'NAME'(const 'OBJECTTYPE' *object)\n"
    tempStr->Append("TH1* ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object)\n");
    fprintf(fp, tempStr->Data());
    break;
  case AliEveListAnalyser::kCorrelObjectAnalyse:
    // Use "void 'NAME'(const TObject *object, const TObject *object2, Double_t*& r, Int_t& n)\n"
    tempStr->Append("void ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object, const ").Append(typeStr->Data()).Append(" *object2, Double_t*& r, Int_t& n)\n");
    fprintf(fp, tempStr->Data());
    break;
  case AliEveListAnalyser::kCorrelObjectHisto:
    // Use "TH1* 'NAME'(const 'OBJECTTYPE' *object, const 'OBJECTTYPE' *object2)\n"
    tempStr->Append("TH1* ").Append(name).Append("(const ").Append(typeStr->Data()).Append(" *object, const ").Append(typeStr->Data()).Append(" *object2)\n");
    fprintf(fp, tempStr->Data());
    break;
  default:
    Error("AliEveGeneralMacroWizard::Create", Form("Unknown type[%d]", type));
    new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error", 
                 Form("Unknown type[%d]", type), kMBIconExclamation, kMBOk);
    fclose(fp);
    gSystem->Exec(Form("rm -f %s.C", name));
    //fCombo->Select(-1);

    tempStr->Clear();
    if (tempStr != 0) delete tempStr;
    tempStr = 0;

    typeStr->Clear();
    if (typeStr != 0) delete typeStr;
    typeStr = 0;

    return;
  }

  tempStr->Clear();
  if (tempStr != 0) delete tempStr;
  tempStr = 0;

  typeStr->Clear();
  if (typeStr != 0) delete typeStr;
  typeStr = 0;

//TODO: Version below?!
  fprintf(fp, "{\n%s\n", fGeneralMacroTemplate[type]);
/*      
  if (useGivenType)
  {
    // Replace "OBJECTTYPE" with the class name
    TString* tempStr = new TString();
    tempStr->Append(fGeneralMacroTemplate_WithType[type]);

    tempStr->ReplaceAll("OBJECTTYPE", fTextObjectType->GetText());

    fprintf(fp, "{\n%s\n", tempStr->Data());

    if (tempStr != 0)
    {
      tempStr->Clear();
      delete tempStr;
      tempStr = 0;
    }
  }
  else
  {
    fprintf(fp, "{\n%s\n", fGeneralMacroTemplate[type]);
  }
*/
  fprintf(fp, "// add your own code here\n\n\n}\n");
  fclose(fp);

  Emit("Create(Int_t)", type);
  Create((Char_t*)name);
  CloseWindow();
}

//______________________________________________________
void AliEveGeneralMacroWizard::Create(Char_t *name)
{
  Emit("Create(Char_t*)", Form("%s.C", name));
}

//______________________________________________________
void AliEveGeneralMacroWizard::HandleCreate()
{
  Create(fCombo->GetSelected());
}
