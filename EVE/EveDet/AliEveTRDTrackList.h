// Author: Benjamin Hess   23/09/2008

/*************************************************************************
 * Copyright (C) 2008, Alexandru Bercuci, Benjamin Hess.                 *
 * All rights reserved.                                                  *
 *************************************************************************/

#ifndef AliEveTRDTrackList_H
#define AliEveTRDTrackList_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveTRDTrackList                                                   //
//                                                                      //
// An AliEveTRDTrackList is, in principal, a TEveElementList with some  //
// sophisticated features. You can add macros to this list, which then  //
// can be applied to the list of tracks (these tracks can be added to   //
// the list in the same way as for the TEveElementList). In general,    //
// please use AddMacro(...) for this purpose.                           //
// Macros that are no longer needed can be removed from the list via    //
// RemoveSelectionMacros(...) or RemoveProcessMacros(...) respectively. //
// This function takes an iterator of the list of entries that are to   //
// be removed. An entry looks like:                                     //
// "MacroName.C (Path: MacroPath)". This is the way, the information    //
// about a macro is stored in the AliEveTRDTrackList. If you have path  //
// and name of a macro, use MakeMacroEntry(...) to get the corresponding//
// entry. The type of the macros is stored in a map. You can get the    //
// macro type via GetMacroType(...).                                    //
// With ApplySTSelectionMacros(...) or ApplyProcessMacros(...)          //
// respectively you can apply the macros to the track list via          //
// iterators (same style like for RemoveProcessMacros(...) (cf. above)).//
// Selection macros (de-)select macros according to a selection rule    //
// by setting the rnr-state of the tracks.                              //
// If multiple selection macros are applied, a track is selected, if    //
// all selection macros select the track.                               //
// Process macros create data or histograms, which will be stored in    //
// a temporary file. The editor of this class will access this file     //
// and draw all the stuff within it's DrawHistos() function. The file   //
// will be deleted by the destructor.                                   //
//                                                                      //
// Currently, the following macro types are supported:                  //
// Selection macros:                                                    //
// Bool_t YourMacro(const AliTRDtrackV1*);                              //
// Bool_t YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*);        //
//                                                                      //
// Process macros:                                                      //
// void YourMacro(const AliTRDtrackV1*, Double_t*&, Int_t&);            //
// void YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*,           //
//                Double_t*&, Int_t&);                                  //
// TH1* YourMacro(const AliTRDtrackV1*);                                //
// TH1* YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*);          //
//                                                                      //
// The macros which take 2 tracks are applied to all pairs              //
// fullfilling the selection criteria.                                  //
//////////////////////////////////////////////////////////////////////////


#include <TEveElement.h>
#include <EveDet/AliEveTRDData.h>

#define SIGNATURE_ERROR   -1
#define NOT_EXIST_ERROR   -2
#define WARNING           0
#define SUCCESS           1

#define MAX_MACRO_NAME_LENGTH     100
#define MAX_MACRO_PATH_LENGTH     300
#define MAX_APPLY_COMMAND_LENGTH  120

#define UNSETBIT(n,i)  ((n) &= ~BIT(i))

class AliEveTRDTrack;
class AliTRDReconstructor;
class TFile;
class TFunction;
class TH1;
class TObjString;
class TList;
class TMap;
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

  // Macro types
  enum AliEveTRDTrackListMacroType
  {
    kUnknown = 0,
    kSingleTrackSelect = 1,
    kSingleTrackAnalyse = 2,
    kSingleTrackHisto = 3,
    kCorrelTrackSelect = 4,
    kCorrelTrackAnalyse = 5,
    kCorrelTrackHisto = 6
  };
  

  AliEveTRDTrackList(const Text_t* n = "AliEveTRDTrackList", const Text_t* t = "", Bool_t doColor = kFALSE);
  virtual ~AliEveTRDTrackList();

  Int_t AddMacro(const Char_t* path, const Char_t* name, Bool_t forceReload = kFALSE);                  
  void AddMacroFast(const Char_t* entry, AliEveTRDTrackListMacroType type);          
  void AddMacroFast(const Char_t* path, const Char_t* name, AliEveTRDTrackListMacroType type);        
  virtual void AddStandardMacros();                           
  Bool_t ApplyProcessMacros(const TList* selIterator, const TList* procIterator);               
  void ApplySTSelectionMacros(const TList* iterator);

  // Returns the type of the macro of the corresponding entry (i.e. "macro.C (Path: path)"). 
  // If you have only the name and the path, you can simply use MakeMacroEntry.
  // If "UseList" is kTRUE, the type will be looked up in the internal list (very fast). But if this list
  // does not exist, you have to use kFALSE for this parameter. Then the type will be determined by the
  // prototype! NOTE: It is assumed that the macro has been compiled! If not, the return value is not
  // predictable, but normally will be kUnknown.
  // Note: AddMacro(Fast) will update the internal list and RemoveProcess(/Selection)Macros respectively.
  AliEveTRDTrackListMacroType GetMacroType(const Char_t* entry, Bool_t UseList = kTRUE) const; 
  Char_t* MakeMacroEntry(const Char_t* path, const Char_t* name) const;  
  void RemoveProcessMacros(const TList* iterator);                   
  void RemoveSelectionMacros(const TList* iterator);                  

protected:
  TList* fMacroList;                 // List of (process) macros
  TList* fMacroSelList;              // List of (selection) macros
  TList* fDataFromMacroList;         // List of macros that currently have data for histograms

  TMap*  fMacroTypes;                // Contains the type of each macro

  TTreeSRedirector *fDataTree;       // Tree containing data for histograms

  Int_t fHistoDataSelected;          // Stores the selection for the data of the histograms
  Int_t fMacroListSelected;          // Stores the selection of the process macro list
  Int_t fMacroSelListSelected;       // Stores the selection of the selection macro list

  Char_t fSelectedTab;               // Holds the index of the selected tab
  UChar_t fSelectedStyle;            // Holds the selected track style

  Char_t GetSelectedTab() const                          // Gets the selected tab
    { return fSelectedTab;  }

  UChar_t GetSelectedTrackStyle() const                  // Gets the selected track style
    { return fSelectedStyle;  }

  Bool_t HistoDataIsSelected(Int_t index) const          // Is entry in list selected?
    { return TESTBIT(fHistoDataSelected, index);  }  
   
  Bool_t MacroListIsSelected(Int_t index) const          // Is entry in list selected?
    { return TESTBIT(fMacroListSelected, index);  }     

  Bool_t MacroSelListIsSelected(Int_t index) const       // Is entry in list selected?
    { return TESTBIT(fMacroSelListSelected, index);  }  

  void SetHistoDataSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fHistoDataSelected, index); else UNSETBIT(fHistoDataSelected, index);  }  

  void SetMacroListSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fMacroListSelected, index); else UNSETBIT(fMacroListSelected, index);  }  

  void SetMacroSelListSelection(Int_t index, Bool_t set)    // Set selection of entry in list
    { if (set) SETBIT(fMacroSelListSelected, index); else UNSETBIT(fMacroSelListSelected, index);  }   
    
  void SetSelectedTab(Int_t index)                          // Sets the selected tab
    { fSelectedTab = (Char_t)index; }  

  void SetSelectedTrackStyle(UChar_t index)                 // Sets the selected track style
    { fSelectedStyle = index;  }

  void UpdateTrackStyle(AliEveTRDTrack::AliEveTRDTrackState s, UChar_t ss = 0); 

private:
  AliEveTRDTrackList(const AliEveTRDTrackList&);            // Not implemented
  AliEveTRDTrackList& operator=(const AliEveTRDTrackList&); // Not implemented             

  ClassDef(AliEveTRDTrackList, 0);  // Class containing a list of tracks
};

#endif
