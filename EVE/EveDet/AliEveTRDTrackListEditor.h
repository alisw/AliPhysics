#ifndef AliEveTRDTrackListEditor_H
#define AliEveTRDTrackListEditor_H

#include <TGedFrame.h>

class AliEveTRDTrack;
class AliEveTRDTrackList;
class TCanvas;     
class TEveBrowser;           
class TEveGedEditor;
class TEveManager;
class TFile;
class TGCheckButton;
class TGFileInfo;
class TGHorizontal3DLine;
class TGHorizontalFrame;
class TGLabel;
class TGListBox;
class TGString;
class TGTab;
class TGTextButton;
class TGTextEntry;
class TGVerticalFrame;
class TH1;
class TTree;

class AliEveTRDTrackListEditor: public TGedFrame
{
public:
  AliEveTRDTrackListEditor(const TGWindow* p = 0, Int_t width = 170, Int_t height = 30,
		                       UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  virtual ~AliEveTRDTrackListEditor();

  virtual void SetModel(TObject* obj);

  void AddMacro(const Char_t* path, const Char_t* name);  // Adds macro to the macro list
  void ApplyMacros();                                     // Apply macros
  void BrowseMacros();                                    // Browse macros
  void CloseTabs();                                       // Closes + deletes all the tabs opened created by this class
  void DrawHistos();                                      // Draw histograms
  Int_t GetNSelectedHistograms();                         // Get the number of selected histograms for drawing
  void HandleMacroPathSet();                              // Handles the "macro path set"-signal
  void HandleNewEventLoaded();                            // Handles the "NewEventLoaded()"-signal
  void HandleTabChangedToIndex(Int_t);                    // Handles the "Selected(Int_t id)"-signal (tab changed)
  void RemoveMacros();                                    // Removes the selected macros from the lists
  void UpdateDataFromMacroListSelection();                // Updates the selection in the "data from macro"-list
  void UpdateHistoList();                                 // Updates the histogram list
  void UpdateMacroList();                                 // Updates the macro list
  void UpdateMacroListSelection(Int_t ind);               // Updates the selection of the process macro list
  void UpdateMacroSelListSelection(Int_t ind);            // Updates the selection of the selection macro list
  
protected:
  AliEveTRDTrackList* fM;                                 // Model object

  void InheritMacroList();                                // Inherits macro list from the previously loaded track list

private:
  AliEveTRDTrackListEditor(const AliEveTRDTrackListEditor&);            // Not implemented
  AliEveTRDTrackListEditor& operator=(const AliEveTRDTrackListEditor&); // Not implemented 

  TCanvas*          fHistoCanvas;            // Canvas for the histograms
  TGString*         fHistoCanvasName;        // Name of the histogram canvas

  Bool_t            fInheritMacroList;       // Flag indicating, whether the macro list will be inherited from the
                                             // previously loaded track list within the next call of SetModel

  TGVerticalFrame*   fMainFrame;             // Top frame for macro functionality.
  TGVerticalFrame*   fHistoFrame;            // Top frame for the histogram stuff
  TGVerticalFrame*   fHistoSubFrame;         // Frame for the histogram buttons themselves
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

  // Help functions
  void SetDrawingToHistoCanvasTab();        // Sets focus on the tab for histograms and makes fHistoCanvas be the
                                            // current tab
  void UpdateHistoCanvasTab();              // Updates the histogram and the corresponding tab (including titles)

  ClassDef(AliEveTRDTrackListEditor, 0);    // Editor for AliEveTRDTrackList.
};

#endif
