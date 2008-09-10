#ifndef AliEveTRDTrackList_H
#define AliEveTRDTrackList_H

#include <TEveElement.h>

#define SIGNATURE_ERROR   -1
#define NOT_EXIST_ERROR   -2
#define WARNING           0
#define SUCCESS           1

class AliEveTRDTrack;
class AliTRDtrackV1;
class TFile;
class TFunction;
class TList;
class TTreeSRedirector;

class AliEveTRDTrackList: public TEveElementList
{
  friend class AliEveTRDTrackListEditor;

public:
  AliEveTRDTrackList(const Text_t* n = "AliEveTRDTrackList", const Text_t* t = "", Bool_t doColor = kFALSE);
  virtual ~AliEveTRDTrackList();

protected:
  TList* fMacroList;                 // List of (process) macros
  TList* fMacroSelList;              // List of (selection) macros
  TList* fDataFromMacroList;         // List of macros that currently have data for histograms

  TTreeSRedirector *fDataTree;       // Tree containing data for histograms

  Int_t AddMacro(const Char_t* path, const Char_t* name);       // Adds a macro (path/name) to the corresponding list
                                                                // (automatic recognition / check)
  void AddMacroFast(const Char_t* path, const Char_t* name,     // Adds a macro (path/name) to the selection (process)
                    Bool_t toSelectionList);                    // macro list, if second parameter is kTRUE (kFALSE).
                                                                // No checks are performed (fast).
  virtual void AddStandardMacros();                             // Adds standard macros to the lists
  void ApplyProcessMacros(TList* iterator);                     // Uses the iterator (for the selected process 
                                                                // macros) to apply the selected macros to the data
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
private:
  AliEveTRDTrackList(const AliEveTRDTrackList&);            // Not implemented
  AliEveTRDTrackList& operator=(const AliEveTRDTrackList&); // Not implemented             

  ClassDef(AliEveTRDTrackList, 0);  // Class containing a list of tracks
};

#endif
