#ifndef AliEveTRDTrackListEditor_H
#define AliEveTRDTrackListEditor_H

#include <TGedFrame.h>

class AliEveTRDTrack;
class AliEveTRDTrackList;
class TFile;
class TGCheckButton;
class TGFileInfo;
class TGHorizontal3DLine;
class TGHorizontalFrame;
class TGLabel;
class TGListBox;
class TGTextButton;
class TGTextEntry;
class TGVerticalFrame;
class TTree;

class AliEveTRDTrackListEditor: public TGedFrame
{
public:
  AliEveTRDTrackListEditor(const TGWindow* p = 0, Int_t width = 170, Int_t height = 30,
		                       UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  virtual ~AliEveTRDTrackListEditor();

  virtual void SetModel(TObject* obj);

  void ApplyMacros();               // Apply macros
  void BrowseMacros();              // Browse macros
  void DrawHistos();                // Draw histograms
  void HandleMacroPathSet();        // Handles the "macro path set"-event 
  void RemoveMacros();              // Removes the selected macros from the lists

protected:
  AliEveTRDTrackList* fM;           // Model object.

private:
  AliEveTRDTrackListEditor(const AliEveTRDTrackListEditor&);            // Not implemented
  AliEveTRDTrackListEditor& operator=(const AliEveTRDTrackListEditor&); // Not implemented 

  void AddMacro(const Char_t* path, const Char_t* name);  // Adds macro to the macro list
  Int_t GetNSelectedHistograms();                         // Get the number of selected histograms for drawing
  void UpdateHistoList();                                 // Updates the histogram list
  void UpdateMacroList();                                 // Updates the macro list


  TGVerticalFrame*  fMainFrame;              // Top frame for macro functionality.
  TGVerticalFrame*  fHistoFrame;             // Top frame for the histogram stuff
  TGVerticalFrame*  fHistoSubFrame;          // Frame for the histogram buttons themselves
  TGHorizontalFrame* fBrowseFrame;           // For searching macros

  TGTextButton*   fbBrowse;                  // "Browse" button
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

  // Check buttons for histograms
  TGCheckButton** fCheckButtons;

  ClassDef(AliEveTRDTrackListEditor, 0);    // Editor for AliEveTRDTrackList.
};

#endif
