#ifndef AliEveTRDTrackListEditor_H
#define AliEveTRDTrackListEditor_H

#include <TGedFrame.h>
#include <TGFileDialog.h>
#include <TGButton.h>
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
#include <EveDet/AliEveTRDTrackList.h>

class AliEveTRDTrackListEditor: public TGedFrame
{
public:
  AliEveTRDTrackListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		           UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDTrackListEditor() {};

  virtual void SetModel(TObject* obj);

  void ApplyMacros();               // Apply macros
  void BrowseMacros();              // Browse macros
  void HandleMacroPathSet();        // Handle "macro path set"-event 
  void RemoveMacros();              // Remove macros

protected:
  AliEveTRDTrackList* fM;           // Model object.

private:
  AliEveTRDTrackListEditor(const AliEveTRDTrackListEditor&);            // Not implemented
  AliEveTRDTrackListEditor& operator=(const AliEveTRDTrackListEditor&); // Not implemented 

  void AddMacro(const Char_t* Entryame, const Char_t* name,     // Add macro to the macro list
                const Char_t* pathname);
  void UpdateMacroList();                                       // Updates the macro list


  TGVerticalFrame*  fMainFrame;             // Top frame for macro functionality.
  TGVerticalFrame*  fMemberFrame;           // Top frame for member list
  TGHorizontalFrame* fBrowseFrame;          // For searching macros

  TGTextButton*   bBrowse;                  // Browse button
  TGTextButton*   bApplyMacros;             // Apply macros button
  TGTextButton*   bRemoveMacros;            // Remove macros button
  TGTextEntry*    teField;                  // Text field to insert macro path manually
  TGTextView*     tvMemberList;             // To display the list of members
  TGListBox*      tlMacroList;              // To display the list of (process) macros
  TGListBox*      tlMacroSelList;           // To display the list of (selection) macros

  TGFileInfo*     fileInfo;                 // Holds data about opening macros
  Char_t**    fileTypes;                    // File types (for macros)

  // Some labels
  TGLabel* fLabel1;
  TGLabel* fLabel2;
  TGLabel* fLabel3;
     
  // Some lines
  TGHorizontal3DLine *fLine1;
  TGHorizontal3DLine *fLine2;
  TGHorizontal3DLine *fLine3;

  ClassDef(AliEveTRDTrackListEditor, 0);    // Editor for AliEveTRDTrackList.
};

#endif
