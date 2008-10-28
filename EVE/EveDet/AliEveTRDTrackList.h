// Author: Benjamin Hess   25/09/2008

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
// RemoveSelectedMacros(...).This function takes an iterator of the     //
// list of macros that are to be removed.                               //
// be removed. An entry looks like:                                     //
// The data for each macro consists of path, name, type and the command //
// that will be used to apply the macro. This stuff is stored in a map  //
// which takes the macro name for the key and the above mentioned data  //
// in a TMacroData-object for the value.                                //
// You can get the macro type via GetMacroType(...).                    //
// With ApplySTSelectionMacros(...) or ApplyProcessMacros(...)          //
// respectively you can apply the macros to the track list via          //
// iterators (same style like for RemoveSelectedMacros(...)(cf. above)).//
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
// The macros which take 2 tracks are applied to all track pairs        //
// (whereby BOTH tracks of the pair have to be selected by the single   //
// track selection macros and have to be unequal, otherwise they will   //
// be skipped) that have been selected by ALL correlated tracks         //
// selection macros. The selection macros with 2 tracks do NOT affect   //
// process macros that process only a single track!                     //
//////////////////////////////////////////////////////////////////////////


#include <TEveElement.h>
#include <EveDet/AliEveTRDData.h>

#define SIGNATURE_ERROR   -1
#define NOT_EXIST_ERROR   -2
#define ERROR            -3
#define WARNING           0
#define SUCCESS           1

#define MAX_MACRO_NAME_LENGTH     100
#define MAX_MACRO_PATH_LENGTH     300
#define MAX_APPLY_COMMAND_LENGTH  120

#define UNSETBIT(n,i)  ((n) &= ~BIT(i))

class AliEveTRDTrack;
class AliEveTRDTrackListEditor;
class AliTRDReconstructor;
class TFile;
class TFunction;
class TH1;
class TObjString;
class TList;
class TMap;
class TPair;
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
  Bool_t AddMacroFast(const Char_t* path, const Char_t* name, AliEveTRDTrackListMacroType type);        
  virtual void AddStandardContent();                           
  Bool_t ApplyProcessMacros(const TList* selIterator, const TList* procIterator);               
  void ApplySTSelectionMacros(const TList* iterator);

  // Returns the type of the macro of the corresponding entry (i.e. "macro.C (Path: path)"). 
  // If you have only the name and the path, you can simply use MakeMacroEntry.
  // If "UseList" is kTRUE, the type will be looked up in the internal list (very fast). But if this list
  // does not exist, you have to use kFALSE for this parameter. Then the type will be determined by the
  // prototype! NOTE: It is assumed that the macro has been compiled! If not, the return value is not
  // predictable, but normally will be kUnknown.
  // Note: AddMacro(Fast) will update the internal list and RemoveProcess(/Selection)Macros respectively.
  AliEveTRDTrackListMacroType GetMacroType(const Char_t* name, Bool_t UseList = kTRUE) const; 
  void RemoveSelectedMacros(const TList* iterator);                                    

protected:
  AliEveTRDTrackListEditor* fEditor; // Pointer to the editor of this list             

  TList* fDataFromMacroList;         // List of macros that currently have data for histograms

  TMap*  fMacroList;                 // Stores the names, paths, types and commands of all macros added to this list

  TTreeSRedirector *fDataTree;       // Tree containing data for histograms

  Int_t fHistoDataSelected;          // Stores the selection for the data of the histograms
  Int_t fMacroListSelected;          // Stores the selection of the macro list

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

  void SetHistoDataSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fHistoDataSelected, index); else UNSETBIT(fHistoDataSelected, index);  }  

  void SetMacroListSelection(Int_t index, Bool_t set)       // Set selection of entry in list
    { if (set) SETBIT(fMacroListSelected, index); else UNSETBIT(fMacroListSelected, index);  }  
    
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMacroData                                                           //
//                                                                      //
// Stores macro data which will be used by AliEveTRDTrackList.          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TMacroData: public TObject
{
public:
  TMacroData(const Char_t* name, const Char_t* path,                                                  // Constructor
             AliEveTRDTrackList::AliEveTRDTrackListMacroType type = AliEveTRDTrackList::kUnknown):
    TObject(),
    fIsSelected(kFALSE),
    fType(AliEveTRDTrackList::kUnknown)
  {
    // The command is automatically set via type.

    SetName(name);
    SetPath(path);
    SetType(type);

    // Register the commands for each type here
    switch (type)
    {
      case AliEveTRDTrackList::kSingleTrackSelect:
      case AliEveTRDTrackList::kSingleTrackHisto:
        SetCmd(Form("%s(automaticTrackV1_1);", name));
        break;

      case AliEveTRDTrackList::kSingleTrackAnalyse:
        SetCmd(Form("%s(automaticTrackV1_1, results, n);", name));
        break;

      case AliEveTRDTrackList::kCorrelTrackSelect:
      case AliEveTRDTrackList::kCorrelTrackHisto:
        SetCmd(Form("%s(automaticTrackV1_1, automaticTrackV1_2);", name));
        break;

      case AliEveTRDTrackList::kCorrelTrackAnalyse:
        SetCmd(Form("%s(automaticTrackV1_1, automaticTrackV1_2, results, n);", name));
        break;

      default:
        SetCmd("");
        break;
    }
  }

  const Char_t* GetCmd() const                // Returns the command that will be used to call this macro
    { return fCmd;  }
  const Char_t* GetName() const               // Returns the macro name (without ".C")
    { return fName;  }
  const Char_t* GetPath() const               // Returns the path of the macro
    { return fPath;  }
  AliEveTRDTrackList::AliEveTRDTrackListMacroType GetType() const   // Returns the type of the macro
    { return fType;  }
  Bool_t IsProcessMacro() const             // Returns whether the macro is a process type macro or not
  {
    switch (fType)
    {
      case AliEveTRDTrackList::kSingleTrackAnalyse:
      case AliEveTRDTrackList::kSingleTrackHisto:
      case AliEveTRDTrackList::kCorrelTrackAnalyse:
      case AliEveTRDTrackList::kCorrelTrackHisto:
        return kTRUE;
        break;
      default:
        break;
    }
    
    return kFALSE;
  }

  Bool_t IsSelected() const                   // Returns whether the macro is selected or not
    { return fIsSelected; }
  Bool_t IsSelectionMacro() const             // Returns whether the macro is a selection type macro or not
  {
    switch (fType)
    {
      case AliEveTRDTrackList::kSingleTrackSelect:
      case AliEveTRDTrackList::kCorrelTrackSelect:
        return kTRUE;
        break;
      default:
        break;
    }
    
    return kFALSE;
  }

  void SetCmd(const char* newCmd)             // Sets the command that will be used to call this macro
  { 
    memset(fCmd, '\0', sizeof(Char_t) * MAX_APPLY_COMMAND_LENGTH);
    sprintf(fCmd, "%s", newCmd);
  }
  void SetName(const char* newName)           // Sets the macro name (please use without ".C")
  { 
    memset(fName, '\0', sizeof(Char_t) * MAX_MACRO_NAME_LENGTH);
    sprintf(fName, "%s", newName);
  }
  void SetPath(const char* newPath)           // Sets the path of the macro
  { 
    memset(fPath, '\0', sizeof(Char_t) * MAX_MACRO_PATH_LENGTH);
    sprintf(fPath, "%s", newPath);
  }
  void SetSelected(Bool_t selection)          // Sets whether the macro is selected or not
  {
    fIsSelected = selection;
  }
  void SetType(AliEveTRDTrackList::AliEveTRDTrackListMacroType newType)  // Sets the type of the macro
  {
    fType = newType;
  }

private:
  Char_t fCmd[MAX_APPLY_COMMAND_LENGTH];                  // Command that will be used to call this macro
  Char_t fName[MAX_MACRO_NAME_LENGTH];                    // Macro name (without ".C"!)
  Char_t fPath[MAX_MACRO_PATH_LENGTH];                    // Path of the macro
  Bool_t fIsSelected;                                     // Is macro selected (e.g. in the editor's list)?
  AliEveTRDTrackList::AliEveTRDTrackListMacroType fType;  // Type of the macro  
};

#endif
