// Author: Benjamin Hess   29/01/2010

/*************************************************************************
 * Copyright (C) 2009-2010, Alexandru Bercuci, Benjamin Hess.            *
 * All rights reserved.                                                  *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEveListAnalyser                                                   //
//                                                                      //
// An AliEveListAnalyser is, in principal, a TEveElementList with some  //
// sophisticated features. You can add macros to this list, which then  //
// can be applied to the list of analysis objects (these objects can be //
// added to the list in the same way as for the TEveElementList, but    //
// also "by clicking" (cf. AliEveListAnaLyserEditor)).                  //
// In general, please use AddMacro(...) for this purpose.               //
// Macros that are no longer needed can be removed from the list via    //
// RemoveSelectedMacros(...). This function takes an iterator of the    //
// list of macros that are to be removed.                               //
// An entry looks like:                                                 //
// The data for each macro consists of path, name, type and the command //
// that will be used to apply the macro. This stuff is stored in a map  //
// which takes the macro name for the key and the above mentioned data  //
// in a TGeneralMacroData-object for the value.                         //
// You can get the macro type via GetMacroType(...).                    //
// To find the type of objects the macro will deal with (corresponds to //
// "YourObjectType" in the examples below) please use                   //
// GetMacroObjectType(...).                                             //
// With ApplySOSelectionMacros(...) or ApplyProcessMacros(...)          //
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
// Bool_t YourMacro(const YourObjectType*);                             //
// Bool_t YourMacro(const YourObjectType*, const YourObjectType2*);     //
//                                                                      //
// Process macros:                                                      //
// void YourMacro(const YourObjectType*, Double_t*&, Int_t&);           //
// void YourMacro(const YourObjectType*, const YourObjectType2*,        //
//                Double_t*&, Int_t&);                                  //
// TH1* YourMacro(const YourObjectType*);                               //
// TH1* YourMacro(const YourObjectType*, const YourObjectType2*);       //
//                                                                      //
// The macros which take 2 tracks are applied to all track pairs        //
// (whereby BOTH tracks of the pair have to be selected by the single   //
// track selection macros and have to be unequal, otherwise they will   //
// be skipped) that have been selected by ALL correlated tracks         //
// selection macros. The selection macros with 2 tracks do NOT affect   //
// process macros that process only a single track!                     //
//////////////////////////////////////////////////////////////////////////

#ifndef AliEveListAnalyser_H
#define AliEveListAnalyser_H

#include <TClass.h>
#include <TEveElement.h>
#include <AliEveTRDData.h>

#define SIGNATURE_ERROR           -1
#define NOT_EXIST_ERROR           -2
#define ERROR                     -3
#define UNKNOWN_OBJECT_TYPE_ERROR -4
#define WARNING                   0
#define SUCCESS                   1

#define MAX_MACRO_NAME_LENGTH     100
#define MAX_MACRO_PATH_LENGTH     300
#define MAX_APPLY_COMMAND_LENGTH  120

#define UNSETBIT(n,i)  ((n) &= ~BIT(i))


class AliEveListAnalyserEditor;
class AliTRDReconstructor;
class TClass;
class TEveManager;
class TEveSelection;
class TFile;
class TFunction;
class TH1;
class TList;
class TMap;
class TObjString;
class TPair;
class TString;
class TTreeSRedirector;

class AliEveListAnalyser: public TEveElementList
{
  friend class AliEveListAnalyserEditor;

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
  enum AliEveListAnalyserMacroType
  {
    kUnknown = 0,
    kSingleObjectSelect = 1,
    kSingleObjectAnalyse = 2,
    kSingleObjectHisto = 3,
    kCorrelObjectSelect = 4,
    kCorrelObjectAnalyse = 5,
    kCorrelObjectHisto = 6
  };
  

  AliEveListAnalyser(const Text_t* n = "AliEveListAnalyser", const Text_t* t = "", Bool_t doColor = kFALSE);
  virtual ~AliEveListAnalyser();

  Int_t AddMacro(const Char_t* path, const Char_t* name, Bool_t forceReload = kFALSE);                      
  Bool_t AddMacroFast(const Char_t* path, const Char_t* name, AliEveListAnalyserMacroType type, TClass* objectType, TClass* objectType2);     
  Int_t AddPrimSelectedObject(TEveElement* el);
  //void AddPrimSelectedObjects();
  void AddSecSelectedSingleObjectToList(Int_t pointId);  
  virtual void AddStandardContent();                             
  Bool_t ApplyProcessMacros(const TList* selIterator, const TList* procIterator);               
  void ApplySOSelectionMacros(const TList* iterator);
  Bool_t GetConnected()             // Returns whether "adding objects by clicking" is enabled or not.
    { return fConnected;  };
  TClass* GetMacroObjectType(const Char_t* name, Int_t argNum = 1) const;

  // Returns the type of the macro of the corresponding entry (i.e. "macro.C (Path: path)"). 
  // If you have only the name and the path, you can simply use MakeMacroEntry.
  // If "UseList" is kTRUE, the type will be looked up in the internal list (very fast). But if this list
  // does not exist, you have to use kFALSE for this parameter. Then the type will be determined by the
  // prototype! NOTE: It is assumed that the macro has been compiled! If not, the return value is not
  // predictable, but normally will be kUnknown.
  // "objectType" gives the class/"type of object" of the pointer(s), the macro accepts as a parametre.
  // Note: AddMacro(Fast) will update the internal list and RemoveProcess(/Selection)Macros respectively.
  AliEveListAnalyserMacroType GetMacroType(const Char_t* name, const Char_t* objectType = "TObject", 
                                           const Char_t* objectType2 = "TObject", Bool_t UseList = kTRUE) const; 
  
  //void RemovePrimSelectedObjects();
  void RemoveSelectedMacros(const TList* iterator);
  void ResetObjectList();
  Bool_t StartAddingObjects();
  Bool_t StopAddingObjects();

protected:     
  Bool_t fConnected;                 // Connection to the TEvePointSet signal

  TList* fDataFromMacroList;         // List of macros that currently have data for histograms

  AliEveListAnalyserEditor* fEditor; // Pointer to the editor of this list

  TMap*  fMacroList;                 // Stores the names, paths, types and commands of all macros added to this list

  TTreeSRedirector *fDataTree;       // Tree containing data for histograms

  Int_t fHistoDataSelected;          // Stores the selection for the data of the histograms
  Int_t fMacroListSelected;          // Stores the selection of the macro list

  Char_t fSelectedTab;               // Holds the index of the selected tab

  Char_t GetSelectedTab() const                          // Gets the selected tab
    { return fSelectedTab;  }

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

private:
  AliEveListAnalyser(const AliEveListAnalyser&);            // Not implemented
  AliEveListAnalyser& operator=(const AliEveListAnalyser&); // Not implemented             

  ClassDef(AliEveListAnalyser, 0);  // Class containing a list of analyse objects
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGeneralMacroData                                                    //
//                                                                      //
// Stores macro data which will be used by AliEveListAnalyser.          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TGeneralMacroData: public TObject
{
public:
  TGeneralMacroData(const Char_t* name, const Char_t* path,                                                        // Constructor
             AliEveListAnalyser::AliEveListAnalyserMacroType type = AliEveListAnalyser::kUnknown,
             TClass* objectType = TObject::Class(), TClass* objectType2 = TObject::Class()):
    TObject(),
    fObjectType(0x0),
    fObjectType2(0x0),
    fIsSelected(kFALSE),
    fType(AliEveListAnalyser::kUnknown)
    
  {
    // The command is automatically set via type.

    SetName(name);
    SetObjectType(objectType);
    SetObjectType2(objectType2);
    SetPath(path);
    SetType(type);

    // Register the commands for each type here
    // Note: The framework will cast all data objects to "TObject*"
    // -> We need to cast to the correct type for the macro!
    switch (type)
    {
      case AliEveListAnalyser::kSingleObjectSelect:
      case AliEveListAnalyser::kSingleObjectHisto:
        SetCmd(Form("%s((%s*)automaticObject_1);", name, fObjectType->GetName()));
        break;

      case AliEveListAnalyser::kSingleObjectAnalyse:
        SetCmd(Form("%s((%s*)automaticObject_1, results, n);", name, fObjectType->GetName()));
        break;

      case AliEveListAnalyser::kCorrelObjectSelect:
      case AliEveListAnalyser::kCorrelObjectHisto:
        SetCmd(Form("%s((%s*)automaticObject_1, (%s*)automaticObject_2);", name, fObjectType->GetName(), fObjectType2->GetName()));
        break;

      case AliEveListAnalyser::kCorrelObjectAnalyse:
        SetCmd(Form("%s((%s*)automaticObject_1, (%s*)automaticObject_2, results, n);", name, fObjectType->GetName(), fObjectType2->GetName()));
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
  TClass* GetObjectType() const               // Returns the object type of the macro parameter (1st pointer)
    { return fObjectType;  }
  TClass* GetObjectType2() const              // Returns the object type of the macro parameter (2nd pointer)
    { return fObjectType2;  }
  const Char_t* GetPath() const               // Returns the path of the macro
    { return fPath;  }
  AliEveListAnalyser::AliEveListAnalyserMacroType GetType() const   // Returns the type of the macro
    { return fType;  }
  Bool_t IsProcessMacro() const               // Returns whether the macro is a process type macro or not
  {
    switch (fType)
    {
      case AliEveListAnalyser::kSingleObjectAnalyse:
      case AliEveListAnalyser::kSingleObjectHisto:
      case AliEveListAnalyser::kCorrelObjectAnalyse:
      case AliEveListAnalyser::kCorrelObjectHisto:
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
      case AliEveListAnalyser::kSingleObjectSelect:
      case AliEveListAnalyser::kCorrelObjectSelect:
        return kTRUE;
        break;
      default:
        break;
    }
    
    return kFALSE;
  }

  void SetCmd(const char* newCmd)               // Sets the command that will be used to call this macro
  { 
    memset(fCmd, '\0', sizeof(Char_t) * MAX_APPLY_COMMAND_LENGTH);
    snprintf(fCmd, MAX_APPLY_COMMAND_LENGTH, "%s", newCmd);
  }
  void SetName(const char* newName)             // Sets the macro name (please use without ".C")
  { 
    memset(fName, '\0', sizeof(Char_t) * MAX_MACRO_NAME_LENGTH);
    snprintf(fName, MAX_MACRO_NAME_LENGTH, "%s", newName);
  }
  void SetObjectType(TClass* newObjectType)     // Sets the object type of the macro parameter (1st pointer)
  { 
    fObjectType = newObjectType;
  }
  void SetObjectType2(TClass* newObjectType2)   // Sets the object type of the macro parameter (2nd pointer)
  { 
    fObjectType2 = newObjectType2;
  }
  void SetPath(const char* newPath)             // Sets the path of the macro
  { 
    memset(fPath, '\0', sizeof(Char_t) * MAX_MACRO_PATH_LENGTH);
    snprintf(fPath, MAX_MACRO_PATH_LENGTH, "%s", newPath);
  }
  void SetSelected(Bool_t selection)            // Sets whether the macro is selected or not
  {
    fIsSelected = selection;
  }
  void SetType(AliEveListAnalyser::AliEveListAnalyserMacroType newType)  // Sets the type of the macro
  {
    fType = newType;
  }

private:
  TGeneralMacroData(const TGeneralMacroData&);            // Not implemented
  TGeneralMacroData& operator=(const TGeneralMacroData&); // Not implemented 

  Char_t fCmd[MAX_APPLY_COMMAND_LENGTH];                  // Command that will be used to call this macro
  Char_t fName[MAX_MACRO_NAME_LENGTH];                    // Macro name (without ".C"!)
  TClass* fObjectType;                                    // Type of object of the macro parameter (1st pointer)
  TClass* fObjectType2;                                   // Type of object of the macro parameter (2nd pointer)
  Char_t fPath[MAX_MACRO_PATH_LENGTH];                    // Path of the macro
  Bool_t fIsSelected;                                     // Is macro selected (e.g. in the editor's list)?
  AliEveListAnalyser::AliEveListAnalyserMacroType fType;  // Type of the macro  
};

#endif
