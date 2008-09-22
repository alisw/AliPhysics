#ifndef AliEveTRDTrackList_H
#define AliEveTRDTrackList_H

#include <TEveElement.h>
#include <EveDet/AliEveTRDData.h>

#define SIGNATURE_ERROR   -1
#define NOT_EXIST_ERROR   -2
#define WARNING           0
#define SUCCESS           1

#define MAX_MACRO_NAME_LENGTH     100
#define MAX_MACRO_PATH_LENGTH     300
#define MAX_APPLY_COMMAND_LENGTH   50

class AliEveTRDTrack;
class AliTRDtrackV1;
class TFile;
class TFunction;
class TH1;
class TObjString;
class TList;
class TTreeSRedirector;

class AliEveTRDTrackList: public TEveElementList
{
  friend class AliEveTRDTrackListEditor;

public:
  enum
  {
    // Maximum length (number of characters) for a macro name:
    fkMaxMacroNameLength = MAX_MACRO_NAME_LENGTH,    
    // Maximum length (number of characters) for a macro path:
    fkMaxMacroPathLength = MAX_MACRO_PATH_LENGTH,    
    // Maximum length (number of characters) for a macro pathname:
    fkMaxMacroPathNameLength = MAX_MACRO_NAME_LENGTH + MAX_MACRO_PATH_LENGTH,  
    // Maximum length (number of characters) for "apply macro" commands in addition to the length of the name, path... 
    fkMaxApplyCommandLength = MAX_APPLY_COMMAND_LENGTH  
  };

  AliEveTRDTrackList(const Text_t* n = "AliEveTRDTrackList", const Text_t* t = "", Bool_t doColor = kFALSE);
  virtual ~AliEveTRDTrackList();

  Int_t AddMacro(const Char_t* path, const Char_t* name,        // Adds a macro (path/name) to the corresponding list
                 Bool_t forceReload = kFALSE);                  // (automatic recognition / check) -> library will be
                                                                // built, if it does not exist, or updated, if the 
                                                                // macro code has been changed. If forceReload is
                                                                // kTRUE, the library will always be (re-)built!
  void AddMacroFast(const Char_t* entry,                        // Adds an entry to the corresponding list (cf. below)
                    Bool_t toSelectionList);     
  void AddMacroFast(const Char_t* path, const Char_t* name,     // Adds a macro (path/name) to the selection (process)
                    Bool_t toSelectionList);                    // macro list, if second parameter is kTRUE (kFALSE).
                                                                // No checks are performed (fast) and no libraries are
                                                                // loaded. Do use only, if library already exists!
  virtual void AddStandardMacros();                             // Adds standard macros to the lists
  Bool_t ApplyProcessMacros(TList* iterator);                   // Uses the iterator (for the selected process 
                                                                // macros) to apply the selected macros to the data.
                                                                // Return kTRUE on success, otherwise kFALSE. If there
                                                                // no process macros selected, kTRUE is returned (no 
                                                                // error!).
  void ApplySelectionMacros(TList* iterator);                   // Uses the iterator (for the selected selection
                                                                // macros) to apply the selected macros to the data
  Char_t* MakeMacroEntry(const Char_t* path, const Char_t* name);  // Constructs an entry for the macro
                                                                   // lists with path and name   
  void RemoveProcessMacros(TList* iterator);                    // Uses the iterator (for the selected process
                                                                // macros) to remove the process macros from 
                                                                // the corresponding list.    
  void RemoveSelectionMacros(TList* iterator);                  // Uses the iterator (for the selected selection
                                                                // macros) to remove the selection macros from 
                                                                // the corresponding list.  
  //void SetTrackStyleState(UChar_t c);   // Sets the track model and the color model for each track in the list

protected:
  TList*  fMacroList;                 // List of (process) macros
  TList*  fMacroSelList;              // List of (selection) macros
  TList*  fDataFromMacroList;         // List of macros that currently have data for histograms

  TTreeSRedirector *fDataTree;       // Tree containing data for histograms

  Int_t   fHistoDataSelected;          // Stores the selection for the data of the histograms
  Int_t   fMacroListSelected;          // Stores the selection of the process macro list
  Int_t   fMacroSelListSelected;       // Stores the selection of the selection macro list

  Char_t  fSelectedTab;                                            // Holds the index of the selected tab
  UChar_t fSelectedStyle;        // Holds the selected track style

  Char_t GetSelectedTab()                                         // Gets the selected tab
    { return fSelectedTab;  }

  UChar_t GetSelectedTrackStyle()     // Gets the selected track style
    { return fSelectedStyle;  }
  Bool_t HistoDataIsSelected(Int_t index)               // Is entry in list selected?
    { return TESTBIT(fHistoDataSelected, index);  }  
   
  Bool_t MacroListIsSelected(Int_t index)               // Is entry in list selected?
    { return TESTBIT(fMacroListSelected, index);  }     

  Bool_t MacroSelListIsSelected(Int_t index)            // Is entry in list selected?
    { return TESTBIT(fMacroSelListSelected, index);  }  

  Bool_t IsHistogramMacro(const Char_t* name);          // Returns kTRUE, if a macro with name "name" has been
                                                        // loaded into a shared library and has the signature
                                                        // of a process macro of type 2 (histogram) -> NO
                                                        // additional check with mangled name!!

  void SetHistoDataSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fHistoDataSelected, index); else CLRBIT(fHistoDataSelected, index);  }  

  void SetMacroListSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fMacroListSelected, index); else CLRBIT(fMacroListSelected, index);  }  

  void SetMacroSelListSelection(Int_t index, Bool_t set)    // Set selection of entry in list
    { if (set) SETBIT(fMacroSelListSelected, index); else CLRBIT(fMacroSelListSelected, index);  }   
    
  void SetSelectedTab(Int_t index)                                          // Sets the selected tab
    { fSelectedTab = (Char_t)index; }  

  void SetSelectedTrackStyle(UChar_t index)     // Sets the selected track style
    { fSelectedStyle = index;  }

  void UpdateTrackStyle(AliEveTRDTrack::AliEveTRDTrackState s, UChar_t ss=0);

private:
  AliEveTRDTrackList(const AliEveTRDTrackList&);            // Not implemented
  AliEveTRDTrackList& operator=(const AliEveTRDTrackList&); // Not implemented             

  ClassDef(AliEveTRDTrackList, 0);  // Class containing a list of tracks
};

#endif
