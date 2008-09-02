#ifndef AliEveTRDTrackListEditor_H
#define AliEveTRDTrackListEditor_H

#include <TGedFrame.h>
#include <TGFileDialog.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGListBox.h>
#include <TGMsgBox.h>
//#include <TEveMacro.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <EveDet/AliEveTRDTrackList.h>

class AliEveTRDTrackListEditor: public TGedFrame
{
public:
  AliEveTRDTrackListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		           UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDTrackListEditor() {};

  virtual void SetModel(TObject* obj);

  void ApplyMacros();         // Apply macros
  void BrowseMacros();        // Browse macros
  void HandleMacroPathSet();  // Handle "macro path set"-event 
  void RemoveMacros();        // Remove macros

protected:
  AliEveTRDTrackList* fM;     // Model object.

private:
  AliEveTRDTrackListEditor(const AliEveTRDTrackListEditor&);            // Not implemented
  AliEveTRDTrackListEditor& operator=(const AliEveTRDTrackListEditor&); // Not implemented 

  void AddMacro(const Char_t* pathname);   // Add macro to the macro list
  void UpdateMacroList();                  // Updates the macro list


  TGVerticalFrame*  fMainFrame; // Top frame for macro functionality.
  TGVerticalFrame*  fMemberFrame; // Top frame for member list
  TGHorizontalFrame* fBrowseFrame; // For searching macros

  TGTextButton*   bBrowse; // Browse button
  TGTextButton*   bApplyMacros; // Apply macros button
  TGTextButton*   bRemoveMacros; // Remove macros button
  TGTextEntry*    teField; // Text field to insert macro path manually
  TGTextView*     tvMemberList; // To display the list of members
  TGListBox*      tlMacroList; // To display the list of macros

  TGFileInfo*     fileInfo; // Holds data about opening macros
  Char_t**    fileTypes; // File types (for macros)

  ClassDef(AliEveTRDTrackListEditor, 0); // Editor for AliEveTRDTrackList.
};

#endif
