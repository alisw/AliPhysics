// Author: Benjamin Hess   25/09/2008

/*************************************************************************
 * Copyright (C) 2008, Alexandru Bercuci, Benjamin Hess.                 *
 * All rights reserved.                                                  *
 *************************************************************************/

#ifndef AliEveTRDTrackListEditor_H
#define AliEveTRDTrackListEditor_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveTRDTrackListEditor                                             //
//                                                                      //
// The AliEveTRDTrackListEditor provides the graphical functionality    //
// for the AliEveTRDTrackList. It creates the tabs and canvases, when   //
// they are needed and, as well, frees allocated memory on destruction  //
// (or if new events are loaded and thus some tabs are closed).         //
// The function DrawHistos() accesses the temporary file created by the //
// AliEveTRDTrackList and draws the desired data (the file will be      //
// created within the call of ApplyMacros()). Have a look at this       //
// function to learn more about the structure of the file and how to    //
// access the data.                                                     //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

#ifndef ROOT_TGFrame
#include <TGFrame.h>
#endif

class AliEveTRDTrack;
class AliEveTRDTrackList;
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

class AliEveTRDTrackListEditor: public TGedFrame
{
public:
  AliEveTRDTrackListEditor(const TGWindow* p = 0, Int_t width = 170, Int_t height = 30,
		                       UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  virtual ~AliEveTRDTrackListEditor();
  virtual void SetModel(TObject* obj);

  void    AddMacro(const Char_t* name, const Char_t* path = ".");  
  void    ApplyMacros();
  void    BrowseMacros();
  void    CloseTabs();
  void    DrawHistos();
  Int_t   GetNSelectedHistograms() const;
  void    HandleMacroPathSet();
  void    HandleNewEventLoaded();
  void    HandleTabChangedToIndex(Int_t);
  void    NewMacros();
  void    RemoveMacros();
  void    SaveMacroList(TMap* list);
  void    SetTrackColor(Int_t ind);
  void    SetTrackModel(Int_t ind);
  void    UpdateDataFromMacroListSelection();
  void    UpdateHistoList();
  void    UpdateMacroList();
  void    UpdateMacroListSelection(Int_t ind);
  
protected:
  AliEveTRDTrackList* fM;                                               // Model object

  void InheritMacroList();                               
  void InheritStyle();                                    

private:
  AliEveTRDTrackListEditor(const AliEveTRDTrackListEditor&);            // Not implemented
  AliEveTRDTrackListEditor& operator=(const AliEveTRDTrackListEditor&); // Not implemented 

  TCanvas*          fHistoCanvas;            // Canvas for the histograms
  TGString*         fHistoCanvasName;        // Name of the histogram canvas

  TMap*             fInheritedMacroList;     // Stores the from the track list inherited macro list

  Bool_t            fInheritSettings;        // Flag indicating, whether the macro list and the style settings will be 
                                             // inherited from the previously loaded track list within the next call 
                                             // of SetModel

  TGHorizontalFrame* fStyleFrame;            // Frame for the style stuff
  TGVerticalFrame*   fMainFrame;             // Top frame for macro functionality.
  TGVerticalFrame*   fHistoFrame;            // Top frame for the histogram stuff
  TGVerticalFrame*   fHistoSubFrame;         // Frame for the histogram buttons themselves
  TGHorizontalFrame* fBrowseFrame;           // Frame for features corresponding to searching macros
  TGButtonGroup*     fbgStyleColor;          // Button group for the color model
  TGButtonGroup*     fbgStyleTrack;          // Button group for the track model
  
  TGRadioButton**    frbColor;               // Radio buttons for the color model
  TGRadioButton**    frbTrack;               // Radio buttons for the track model

  TGTextButton*   fbBrowse;                  // "Browse" button
  TGTextButton*   fbNew;                     // "New" button
  TGTextButton*   fbApplyMacros;             // "Apply macros" button
  TGTextButton*   fbRemoveMacros;            // "Remove macros" button
  TGTextButton*   fbDrawHisto;               // "Draw histogram" button
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
  TGHorizontal3DLine *fLine5;  

  TGCheckButton** fCheckButtons;            // Check buttons for histograms

  // Help functions
  void SetDrawingToHistoCanvasTab();        
  void UpdateHistoCanvasTab();              

  ClassDef(AliEveTRDTrackListEditor, 0);    // Editor for AliEveTRDTrackList.
};

class TGTextEdit;
class TGComboBox;
class AliEveTRDMacroWizzard : public TGMainFrame
{
public:
  AliEveTRDMacroWizzard(const TGWindow* p = 0);
  void Create(Int_t typ);
  void Done(Char_t *macro)
    { Emit("Done(Char_t*)", macro); } //*SIGNAL*
private:
  AliEveTRDMacroWizzard(const AliEveTRDMacroWizzard&);
  AliEveTRDMacroWizzard& operator=(const AliEveTRDMacroWizzard&);

  TGTextEntry *fText;
  TGComboBox  *fCombo;
  TGTextEdit  *fTextEdit;
  
  ClassDef(AliEveTRDMacroWizzard, 0);    // Helper class to create macro templates 
};

#endif
