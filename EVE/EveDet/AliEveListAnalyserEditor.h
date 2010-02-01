// Author: Benjamin Hess   29/01/2010

/*************************************************************************
 * Copyright (C) 2009-2010, Alexandru Bercuci, Benjamin Hess.            *
 * All rights reserved.                                                  *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveListAnalyserEditor                                             //
//                                                                      //
// The AliEveListAnalyserEditor provides the graphical func-            //
// tionality for the AliEveListAnalyser. It creates the tabs            //
// and canvases, when they are needed and, as well, frees allocated     //
// memory on destruction (or if new events are loaded and thus some     //
// tabs are closed).                                                    //
// The function DrawHistos() accesses the temporary file created by the //
// AliEveListAnalyser and draws the desired data (the file will         //
// be created within the call of ApplyMacros()). Have a look at this    //
// function to learn more about the structure of the file and how to    //
// access the data.                                                     //
// You can add objects to the list (of analysis objects) "by clicking"! //
// To do this, click the "start" button in the "list" tab. Pressing it, //
// connects the class to signals of objects in the viewer.              //
// You have to kinds of selection:                                      //
//                                                                      //
// Secondary selection:                                                 //
// You can hold "CTRL"+"ALT" (depending on your system, "ALT" alone can //
// also be fine) and click an single object (e.g. a single cluster of a //
// TEvePointSet) in the viewer to add it to the list. If the object is  //
// already in the list, it will be removed from it!                     //
//                                                                      //
// Primary selection:                                                   //
// Just click the object you want to add in the viewer (or as well in   //
// the browser (left panel)). If the object is already in the list, it  //
// will be removed from it!                                             //
//                                                                      //
// For both cases: Note:                                                //
// If you have added all the desired objects, please press the "stop"   //
// button in the "list" tab to disconnect the class from the signals.   //
// If you want to remove an object, you HAVE to use the same procedure  //
// that you have used for adding it. e.g. you cannot(!) remove an       //
// object added by the secondary selection method by using the primary  //
// selection method!                                                    //
//////////////////////////////////////////////////////////////////////////

#ifndef AliEveListAnalyserEditor_H
#define AliEveListAnalyserEditor_H

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

#ifndef ROOT_TGFrame
#include <TGFrame.h>
#endif

class AliEveListAnalyser;
class AliTRDReconstructor;
class TCanvas;     
class TEveBrowser;           
class TEveGedEditor;
class TEveManager;
class TFile;
class TGButtonGroup;
class TGCheckButton;
class TGFileInfo;
class TGGroupFrame;
class TGHorizontal3DLine;
class TGHorizontalFrame;
class TGLabel;
class TGListBox;
class TGRadioButton;
class TGString;
class TGTab;
class TGTextButton;
class TGTextEntry;
class TGVerticalFrame;
class TH1;
class TMacroData;
class TMap;
class TMapIter;
class TTree;

class AliEveListAnalyserEditor: public TGedFrame
{
public:
  AliEveListAnalyserEditor(const TGWindow* p = 0, Int_t width = 170, Int_t height = 30,
		                       UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  virtual ~AliEveListAnalyserEditor();
  virtual void SetModel(TObject* obj);

  void    AddMacro(const Char_t* name, const Char_t* path = ".");  
  void    ApplyMacros();
  void    BrowseMacros();
  void    CloseTabs();
  //void    DoAddPrimSelectedObjects();
  //void    DoRemovePrimSelectedObjects();
  void    DoResetObjectList();
  void    DoStartAddingObjects();
  void    DoStopAddingObjects();
  void    DrawHistos();
  Int_t   GetNSelectedHistograms() const;
  void    HandleMacroPathSet();
  void    HandleNewEventLoaded();
  void    HandleTabChangedToIndex(Int_t);
  void    NewMacros();
  void    RemoveMacros();
  void    SaveMacroList(TMap* list);
  void    UpdateDataFromMacroListSelection();
  void    UpdateHistoList();
  void    UpdateMacroList();
  void    UpdateMacroListSelection(Int_t ind);
  
protected:
  AliEveListAnalyser* fM;                                               // Model object

  void InheritMacroList();                                                                

private:
  AliEveListAnalyserEditor(const AliEveListAnalyserEditor&);            // Not implemented
  AliEveListAnalyserEditor& operator=(const AliEveListAnalyserEditor&); // Not implemented 

  // Help functions
  void SetDrawingToHistoCanvasTab();        
  void UpdateHistoCanvasTab();             

  TCanvas*          fHistoCanvas;            // Canvas for the histograms
  TGString*         fHistoCanvasName;        // Name of the histogram canvas

  TMap*             fInheritedMacroList;     // Stores the from the analyse object list inherited macro list

  Bool_t            fInheritSettings;        // Flag indicating, whether the macro list will be inherited from
                                             // the previously loaded analyse object list within the next call of SetModel

  TGHorizontalFrame* fBrowseFrame;           // Frame for features corresponding to searching macros
  TGVerticalFrame*   fHistoFrame;            // Top frame for the histogram stuff
  TGVerticalFrame*   fHistoSubFrame;         // Frame for the histogram buttons themselves
  TGVerticalFrame*   fMainFrame;             // Top frame for macro functionality.
  TGVerticalFrame*   fObjectFrame;           // Frame for features corresponding to adding objects to the list
  
  //TGTextButton*   fbAddPrimObjects;          // "Add selected object(s)" button
  TGTextButton*   fbApplyMacros;             // "Apply macros" button
  TGTextButton*   fbBrowse;                  // "Browse" button
  TGTextButton*   fbDrawHisto;               // "Draw histogram" button
  TGTextButton*   fbNew;                     // "New" button  
  //TGTextButton*   fbRemovePrimObjects;       // "Remove selected object(s)" button
  TGTextButton*   fbRemoveMacros;            // "Remove macros" button
  TGTextButton*   fbReset;                   // "Reset" (list of added objects) button
  TGTextButton*   fbStart;                   // "Start" (adding objects to list) button
  TGTextButton*   fbStop;                    // "Stop" (adding objects to list) button
  TGTextEntry*    fteField;                  // Text field to insert macro path manually
  TGListBox*      ftlMacroList;              // To display the list of (process) macros
  TGListBox*      ftlMacroSelList;           // To display the list of (selection) macros

  TGFileInfo*     fFileInfo;                 // Holds data about opening macros
  Char_t**        fFileTypes;                // File types (for macros)

  // Some labels
  TGLabel* fLabel1;
  TGLabel* fLabel2;
  TGLabel* fLabel3;
  TGLabel* fLabel4;
     
  // Some lines
  TGHorizontal3DLine *fLine1;
  TGHorizontal3DLine *fLine2;
  TGHorizontal3DLine *fLine3;
  TGHorizontal3DLine *fLine4; 

  TGCheckButton** fCheckButtons;            // Check buttons for histograms 

  ClassDef(AliEveListAnalyserEditor, 0);    // Editor for AliEveListAnalyser.
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveGeneralMacroWizard                                             //
//                                                                      //
// Wizard for creating new macros.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TGTextEdit;
class TGComboBox;
class AliEveGeneralMacroWizard : public TGMainFrame
{
public:
  AliEveGeneralMacroWizard(const TGWindow* p = 0);
  void Create(Int_t type); //*SIGNAL*
  void Create(Char_t *pname); //*SIGNAL*
  void HandleCreate();
  void HandleSelectionChanged(Int_t sel);

private:
  AliEveGeneralMacroWizard(const AliEveGeneralMacroWizard&);
  AliEveGeneralMacroWizard& operator=(const AliEveGeneralMacroWizard&);

  TGTextButton *fbCancel;                  // "Cancel" button
  TGTextButton *fbCreate;                  // "Done" button
  TGComboBox   *fCombo;                    // "Type"
  TGTextEdit   *fTextEdit;                 // "Comments"
  TGTextEntry  *fTextIncludes;             // "Includes"
  TGTextEntry  *fTextName;                 // "Name"
  TGTextEntry  *fTextObjectType;           // "1st object type"  
  TGTextEntry  *fTextObjectType2;          // "2nd object type"
  
  ClassDef(AliEveGeneralMacroWizard, 0);      // Helper class to create macro templates 
};

#endif
