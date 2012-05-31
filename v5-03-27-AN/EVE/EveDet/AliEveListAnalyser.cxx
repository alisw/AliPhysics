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


// Uncomment to display debugging infos
//#define AliEveListAnalyser_DEBUG

#include <TEveManager.h>
#include <TEveSelection.h>
#include <TFile.h>
#include <TFunction.h>
#include <TH1.h>
#include <TList.h>
#include <TMap.h>
#include <TMethodArg.h>
#include <TMethodCall.h>
#include <TObjString.h>
#include <TQObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeStream.h>

#include <AliTRDReconstructor.h>

#include <EveDet/AliEveListAnalyser.h>
#include <EveDet/AliEveListAnalyserEditor.h>

ClassImp(AliEveListAnalyser)

///////////////////////////////////////////////////////////
/////////////   AliEveListAnalyser ////////////////////////
///////////////////////////////////////////////////////////
AliEveListAnalyser::AliEveListAnalyser(const Text_t* n, const Text_t* t, Bool_t doColor):
  TEveElementList(n, t, doColor),
  fConnected(kFALSE),
  fDataFromMacroList(0x0),
  fEditor(0x0),
  fMacroList(0x0),
  fDataTree(0x0),
  fHistoDataSelected(0),
  fMacroListSelected(0),
  fSelectedTab(2)                               // Standard tab: "Apply macros" (index 2)
{
  // Creates the AliEveListAnalyser.

  // Only accept childs of type TEveElement
  SetChildClass(TEveElement::Class());

  // Allocate memory for the lists and declare them as owners of their contents
  fDataFromMacroList = new TList();
  fDataFromMacroList->TCollection::SetOwner(kTRUE);

  fMacroList = new TMap();
  // Set map to owner of it's objects to delete them, if they are removed from the map
  fMacroList->SetOwnerKeyValue(kTRUE, kTRUE);

  // Set the build directory for AClic
  if(gSystem->AccessPathName(Form("%s/.QArec" , gSystem->Getenv("HOME")))) gSystem->Exec("mkdir $HOME/.QArec");
  gSystem->SetBuildDir(Form("%s/.QArec", gSystem->Getenv("HOME")));

  AddStandardContent();
}

//______________________________________________________
AliEveListAnalyser::~AliEveListAnalyser()
{
  // Frees allocated memory (lists etc.).

  // Stop adding objects
  StopAddingObjects();

  // Let the editor know that the list will be destroyed -> The editor will save the data
  if (fEditor != 0)
  {
    fEditor->SaveMacroList(fMacroList);
    fEditor = 0;
  }

  if (fDataFromMacroList != 0)
  {
    fDataFromMacroList->Delete();
    delete fDataFromMacroList;
    fDataFromMacroList = 0;
  } 
  if (fDataTree != 0)
  {
    delete fDataTree;
    fDataTree = 0;
  } 
  if (fMacroList != 0)
  {
    fMacroList->DeleteAll();
    delete fMacroList;
    fMacroList = 0;
  }
  // Note: gSystem->AccessPathName(...) returns kTRUE, if the access FAILED!
  if(!gSystem->AccessPathName(Form("/tmp/ListAnalyserMacroData_%s.root", gSystem->Getenv("USER")))) 
    gSystem->Exec(Form("rm /tmp/ListAnalyserMacroData_%s.root", gSystem->Getenv("USER")));
}

//______________________________________________________
Int_t AliEveListAnalyser::AddMacro(const Char_t* path, const Char_t* nameC, Bool_t forceReload)
{
  // Checks, if the file exists and if the signature is correct.
  // If these criteria are fullfilled, the library for this macro is built
  // and the macro is added to the corresponding list.
  // Supported macro types:
  // Selection macros:                                                    
  // Bool_t YourMacro(const YourObjectType*)
  // Bool_t YourMacro(const YourObjectType*, const YourObjectType2*)
  //
  // Process macros:                                                      
  // void YourMacro(const YourObjectType*, Double_t*&, Int_t&)             
  // void YourMacro(const YourObjectType*, const YourObjectType2*, Double_t*&, Int_t&)                                   
  // TH1* YourMacro(const YourObjectType*)                                 
  // TH1* YourMacro(const YourObjectType*, const YourObjectType2*)                              

  Char_t pathname[fkMaxMacroPathNameLength];
  memset(pathname, '\0', sizeof(Char_t) * fkMaxMacroPathNameLength);

  // Expand the path and create the pathname
  Char_t* systemPath = gSystem->ExpandPathName(path);
  snprintf(pathname, fkMaxMacroPathNameLength, "%s/%s", systemPath, nameC);
  delete systemPath;
  systemPath = 0;

  // Delete ".C" from filename
  Char_t name[fkMaxMacroNameLength];
  memset(name, '\0', sizeof(Char_t) * fkMaxMacroNameLength);
  
  for (UInt_t ind = 0; ind < fkMaxMacroNameLength && ind < strlen(nameC) - 2; ind++)  name[ind] = nameC[ind];

  // Check, if files exists
  FILE* fp = 0x0;
  if((fp = fopen(pathname, "rb"))){
    fclose(fp);
    fp = 0x0;
  } else  return NOT_EXIST_ERROR;
  
  // Clean up root, load the desired macro and then check the type of the macro
  // A.B. gROOT->Reset();
 
  gROOT->ProcessLineSync(Form(".L %s+%c", pathname, forceReload ? '+' : ' '));

  TClass* objectType;
  TClass* objectType2;

  objectType = GetMacroObjectType(name, 1);
  objectType2 = GetMacroObjectType(name, 2);

  // We need this line... otherwise, in some cases, there will be problems concerning ACLIC
  gROOT->ProcessLineSync(Form(".L %s", pathname));

  if (!objectType)  return UNKNOWN_OBJECT_TYPE_ERROR;

  // This might be a macro dealing with only 1 object... test this afterwards!
  Bool_t testSecondObj = kFALSE;
  if (!objectType2) 
  {
    objectType2 = TObject::Class();
    testSecondObj = kTRUE;
  }
  AliEveListAnalyserMacroType type = GetMacroType(name, objectType->GetName(), objectType2->GetName(), kFALSE);

  if (testSecondObj)
  {
    switch (type)
    {
    case AliEveListAnalyser::kCorrelObjectSelect:
    case AliEveListAnalyser::kCorrelObjectAnalyse:
    case AliEveListAnalyser::kCorrelObjectHisto:
      // There must be a second type -> Error!
      return UNKNOWN_OBJECT_TYPE_ERROR;
      break;
    default:
      // Ok, single object macro!
      break;
    }
  }


  // Clean up again
  // A.B. gROOT->Reset();
 
  // Has not the correct signature!
  if (type == kUnknown) return SIGNATURE_ERROR;

  // Only add macro, if it is not already in the list
  Int_t returnValue = WARNING;
  if(fMacroList->GetValue(name) == 0) 
  {
    returnValue = AddMacroFast(path, name, type, objectType, objectType2) ? SUCCESS : ERROR;
  }

  return returnValue;
}

//______________________________________________________
Bool_t AliEveListAnalyser::AddMacroFast(const Char_t* path, const Char_t* name, AliEveListAnalyserMacroType type, 
                                        TClass* objectType, TClass* objectType2)
{
  // Adds a macro (path/name) to the corresponding list. No checks are performed (file exists, 
  // macro already in list/map, signature correct),  no libraries are created!
  // You can use this function only, if the macro has been added successfully before 
  // (and then maybe was removed). The function is very fast. On success kTRUE is returned, otherwise: kFALSE;
  // Note: If your macro takes only 1 pointer as a parameter, just use "0x0" for objectType2!

  Bool_t success = kFALSE;

  switch (type)
  {
    case kSingleObjectSelect:
    case kCorrelObjectSelect:
    case kSingleObjectAnalyse:
    case kSingleObjectHisto:
    case kCorrelObjectAnalyse:
    case kCorrelObjectHisto:
      fMacroList->Add(new TObjString(name), new TGeneralMacroData(name, path, type, objectType, objectType2));

      // We do not know, where the element has been inserted - deselect this list
      fMacroListSelected = 0;
    
      success = kTRUE;

#ifdef AliEveListAnalyser_DEBUG
      // Successfull add will only be displayed in debug mode
      printf("AliEveListAnalyser::AddMacroFast: Added macro \"%s/%s\" with object types \"%s\" and \"%s\" to the corresponding list\n", 
             path, name, objectType->GetName(), objectType2->GetName());
#endif

      break;

    default:
      // Error will always be displayed
      printf("AliEveListAnalyser::AddMacroFast: ERROR: Could not add macro \"%s/%s\" with object types \"%s\" and \"%s\" to the corresponding list\n", path, name, objectType->GetName(), objectType2->GetName());

      success = kFALSE;

      break;
  }

  return success;
}


//______________________________________________________
Int_t AliEveListAnalyser::AddPrimSelectedObject(TEveElement* el)
{
  // Adds the TEveElement el to the list. If it is already in the list, it is removed.
  // If the list is the only parent of the clicked object, the object is moved outside the list in the browser (not deleted!).
  // If you want to delete the object, just select it there and choose "Destroy" in the menu.
  // This function is designed to be used together with a signal:
  // It adds the (primarily) selected objects in the viewer to the list (objects that are already in the list are removed!).
  // Returns "ERROR" (cf. defines) on error, "WARNING" if the element does not contain any user data and else "SUCCESS" (i.e.
  // the element has been added successfully or the element is the list itself and therefore ignored, or the element is ignored 
  // because it has been added via secondary selection).
 
  if (!el)
  {
    Error("AliEveListAnalyser::AddPrimSelectedObject", "Zero pointer!");

    return ERROR;
  }

  // If the clicked element is the list itself, just do nothing.
  if (el == this)  return SUCCESS;

  if (!this->HasChild(el))
  {

    // Ignore objects that do not have any user data, since these cannot be used in the analysis!
    if (el->GetUserData() == 0x0)
    {
      Warning("AddPrimSelectedObject", "Selected object does not contain any \"user data\" and therefore is ignored!");

      return WARNING;
    }

    // Element clicked that is not in the list (and is not the list itself!) -> Add this element to the list
    this->AddElement(el);
    this->SetTitle(Form("Objects %d", this->NumChildren()));
    gEve->Redraw3D();
  }
  else
  {
    // Element clicked that is already in the list. Remove it. But: Only take care of objects that have been added
    // via primary selection (name does not start with "[sec")
    if (TString(el->GetElementName()).BeginsWith("[sec:"))  return SUCCESS;


    // Element is a child of this list. So, if there are only 2 parents, we know them: list + eve selection. In this case,
    // the element needs to be destroyed. If there are more parents, just remove the element from the list.
    // Since the elements editor will be opened, the element is not deleted, but "moved" outside the list (in the browser).
    if (el->NumParents() > 2) 
    {
      this->RemoveElement(el);
    }
    else
    {
      // There must be at least 2 parents!
      if (el->NumParents() <= 1)  return ERROR;

      TEveElement* listObj = 0x0;
      listObj = this->FindChild(el->GetElementName());
      if (!listObj)  return ERROR;

      gEve->AddElement(listObj, 0);
      // Alternatively: Switch on that the element is NOT destroyed, instead of adding it outside the list. Warning: Memory leaks possible.
      //listObj->SetDestroyOnZeroRefCnt(kFALSE);
      this->RemoveElement(listObj);
      //gEve->RemoveElement(listObj, 0);
    }
  } 

  this->SetTitle(Form("Objects %d", this->NumChildren()));
  gEve->Redraw3D();

  return SUCCESS;
}


/*
//______________________________________________________
void AliEveListAnalyser::AddPrimSelectedObjects()
{
  // Adds the (primarily) selected objects in the viewer to the list (objects that are already in the list are ignored).
  // Hold the CTRL-key for multiple selection.

  TEveSelection* eveSel = gEve->GetSelection();
  if (!eveSel)
  {
    Error("AliEveListAnalyser::AddPrimSelectedObjects", "Failed to get the selection!\n");
    return;
  }
  
  TEveElement* elem = 0x0;
  Bool_t changedSomething = kFALSE;

  for (TEveElement::List_i iter = eveSel->BeginChildren(); iter != eveSel->EndChildren(); ++iter)
  {
    if(!(elem = dynamic_cast<TEveElement*>(*iter))) continue;

    if (!this->HasChild(elem) && elem != this)
    {
      // Element clicked that is not in the list (and is not the list itself!) -> Add this element to list
      this->AddElement(elem);
      this->SetTitle(Form("Objects %d", this->NumChildren()));
      changedSomething = kTRUE;
    }
  }

  if (changedSomething) gEve->Redraw3D();
}
*/

//______________________________________________________
void AliEveListAnalyser::AddSecSelectedSingleObjectToList(Int_t pointId)
{
  // This function adds the selected object (secondary selection in the viewer) to the list
  // of analysis objects. If the object is already in the list, it will be removed from it.
  // This function is used to add single objects of a TEvePointset, e.g. single clusters.

  TEvePointSet* ps = dynamic_cast<TEvePointSet*>((TQObject*) gTQSender);
  if (!ps)
  {
    Error("AliEveListAnalyser::AddSecSelectedSingleObjectToList", "Zero pointer!");
    return;
  }

  // Check, if object is already there. If so, remove it!
  
  // 1st possibility: Object of the list clicked. But: Only take care of objects that have been added
  // via secondary selection (name starts with "[sec"). Note: HasChild will also return kTRUE, if e.g.
  // the whole TEvePointSet of clusters is in the last (but maybe another point of it has been clicked!)

  if (this->HasChild(ps))
  {
    if (TString(ps->GetName()).BeginsWith("[sec:"))
    {
      // I don't know why, but you get a crash, if you try this->RemoveElement(ps) (in some cases).
      // So, use this way instead.
      TEveElement* listObj = this->FindChild(ps->GetName());
      if (listObj)
      {
        listObj->SetUserData(0x0);
        this->RemoveElement(listObj);
        this->SetTitle(Form("Objects %d", this->NumChildren()));
      }

      return;
    }
  }

  TObject* obj = ps->GetPointId(pointId);
  if (obj)
  {
    // 2nd possibility: Same object clicked again
    TEveElement* listObj = 0x0;
    listObj = this->FindChild(Form("[sec:%d] %s%d", obj->GetUniqueID(), obj->GetName(), pointId));
    if (listObj)
    {
      listObj->SetUserData(0x0);
      this->RemoveElement(listObj); 
      this->SetTitle(Form("Objects %d", this->NumChildren()));
      return;
    }

    // Object clicked that is not in the list -> Add this object to list
    TEvePointSet* newPS = new TEvePointSet(Form("[sec:%d] %s%d", obj->GetUniqueID(), obj->GetName(), pointId));
    Double_t x = 0, y = 0, z = 0;
    ps->GetPoint(pointId, x, y, z);
    newPS->SetPoint(0, x, y, z);
    newPS->SetUserData(obj);
    // Choose yellow for the added points and inherit style and size for the marker
    newPS->SetMarkerColor(5);
    newPS->SetMarkerStyle(ps->GetMarkerStyle());
    newPS->SetMarkerSize(ps->GetMarkerSize());
    // Own points -> Will be cleared, if this object is removed
    newPS->SetOwnIds(kTRUE);

    this->AddElement(newPS);
    this->SetTitle(Form("Objects %d", this->NumChildren()));
    gEve->Redraw3D();
  }
  else
  {
    Error("AliEveListAnalyser::AddSecSelectedSingleObjectToList", "Selected object is NULL and therefore ignored!");
  }
}

//______________________________________________________
void AliEveListAnalyser::AddStandardContent()
{
  // Adds standard macros to the macro list.

  // Add your standard macros here, e.g.: 
  // To add a macro use:
  // AddMacro("$(ALICE_ROOT)/myFolder", "myMacroName.C");
  // -> If the file does not exist, nothing happens. So if you want to handle this,
  // use the return value of AddMacro (NOT_EXIST_ERROR is returned, if file does not exist)
  // (-> You can also check for other return values (see AddMacro(...)))

}

//______________________________________________________
Bool_t AliEveListAnalyser::ApplyProcessMacros(const TList* selIterator, const TList* procIterator)
{
  // Uses the procIterator (for the selected process macros) to apply the selected macros to the data.
  // Returns kTRUE on success, otherwise kFALSE. If there no process macros selected, kTRUE is returned 
  // (this is no error!).
  // The single object process macros are applied to all selected objects.
  // The selIterator (for the selected selection macros) will be used to apply the correlated objects selection
  // macros to all object pairs (whereby BOTH objects have to be selected, otherwise they will be skipped).
  // All object pairs that have been selected by ALL correlated objects selection macros will be processed by
  // the correlated objects process macros.

  // No process macros need to be processed
  if (procIterator->GetEntries() <= 0)  return kTRUE;

  // Clear root
  // A.B. gROOT->Reset();
  
  // Clear old data and re-allocate
  if (fDataTree == NULL){
    TDirectory *cwd = gDirectory;
    fDataTree = new TTreeSRedirector(Form("/tmp/ListAnalyserMacroData_%s.root", gSystem->Getenv("USER")));
    cwd->cd();
  }
  if (!fDataTree){
    Error("Apply process macros", "File \"/tmp/ListAnalyserMacroData_%s.root\" could not be accessed properly!",
          gSystem->Getenv("USER"));
    return kFALSE;
  }
  
  if (fDataFromMacroList != 0) {
    fDataFromMacroList->Delete();
    delete fDataFromMacroList;
  }
  fDataFromMacroList = new TList();
  fDataFromMacroList->TCollection::SetOwner(kTRUE);

  fHistoDataSelected = 0;
  TGeneralMacroData* macro(NULL);

  TString* procCmds                      = new TString[procIterator->GetEntries()];
  AliEveListAnalyserMacroType* mProcType = new AliEveListAnalyserMacroType[procIterator->GetEntries()];
  TClass** mProcObjectType               = new TClass*[procIterator->GetEntries()];
  TClass** mProcObjectType2              = new TClass*[procIterator->GetEntries()];

  TString* selCmds(NULL);
  AliEveListAnalyserMacroType* mSelType(NULL);
  TClass** mSelObjectType(NULL);
  TClass** mSelObjectType2(NULL);
  
  Bool_t selectedByCorrSelMacro = kFALSE;

  AliEveListAnalyserMacroType macroType = kUnknown;
  Int_t numHistoMacros = 0;
  TH1** histos(NULL);

  TEveElement* object1(NULL);
  TEveElement* object2(NULL);
  TH1* returnedHist(NULL);

  // Collect the commands for each process macro and add them to "data-from-list"
  for (Int_t i = 0; i < procIterator->GetEntries(); i++){
    macro = (TGeneralMacroData*)fMacroList->GetValue(procIterator->At(i)->GetTitle());

    if (!macro){
      Error("Apply process macros", 
        "Macro list is corrupted: Macro \"%s\" is not registered!",
        procIterator->At(i)->GetTitle());
      continue;
    }

#ifdef AliEveListAnalyser_DEBUG
    printf("AliEveListAnalyser: Checking process macro: %s\n", macro->GetName());
#endif 
           
    // Find the object types of the macro
    mProcObjectType[i] = macro->GetObjectType();
    mProcObjectType2[i] = macro->GetObjectType2();

    // Find the type of the process macro
    macroType = macro->GetType();
    if (macroType == kSingleObjectHisto || macroType == kCorrelObjectHisto){
      mProcType[i] = macroType;
      numHistoMacros++;
      // Create the command 
      procCmds[i] = macro->GetCmd();

      // Add to "data-from-list" -> Mark as a histo macro with the substring "(histo macro)"
      fDataFromMacroList->Add(new TObjString(Form("%s (histo macro)", macro->GetName())));
    } else if (macroType == kSingleObjectAnalyse || macroType == kCorrelObjectAnalyse) {
      mProcType[i] = macroType;
      // Create the command 
      procCmds[i] = macro->GetCmd();

      // Add to "data-from-list"
      fDataFromMacroList->Add(new TObjString(macro->GetName()));
    } else {
      Error("Apply process macros", 
        "Macro list corrupted: Macro \"%s/%s.C\" is not registered as a process macro!",
        macro->GetPath(), macro->GetName());
      mProcType[i] = kUnknown;
    } 
  }  


  Int_t selEntries = selIterator->GetEntries();
  // Collect the commands for each selection macro and add them to "data-from-list"
  if (selEntries > 0) {
    selCmds         = new TString[selEntries];
    mSelType        = new AliEveListAnalyserMacroType[selEntries];
    mSelObjectType  = new TClass*[selEntries];
    mSelObjectType2 = new TClass*[selEntries];
    for (Int_t i = 0; i < selEntries; i++){
      macro = (TGeneralMacroData*)fMacroList->GetValue(selIterator->At(i)->GetTitle());

      if (!macro){
        Error("Apply process macros",
          "Macro list is corrupted: Macro \"%s\" is not registered!",
          selIterator->At(i)->GetTitle());
        continue;
      }

#ifdef AliEveListAnalyser_DEBUG
      printf("AliEveListAnalyser: Checking selection macro: %s\n", macro->GetName());
#endif

      // Find the object types of the macro
      mSelObjectType[i] = macro->GetObjectType();
      mSelObjectType2[i] = macro->GetObjectType2();

      // Find the type of the process macro
      macroType = macro->GetType();

      // Single Object select macro
      if (macroType == kSingleObjectSelect) {
        // Has already been processed by ApplySOSelectionMacros(...)
        mSelType[i] = macroType;
      }
      // Correlated Objects select macro
      else if (macroType == kCorrelObjectSelect) {
        mSelType[i] = macroType;

        // Create the command
        selCmds[i] = macro->GetCmd();
      } else {
        Error("Apply process macros",
          "Macro list corrupted: Macro \"%s/%s.C\" is not registered as a selection macro!",
          macro->GetPath(), macro->GetName());
        mSelType[i] = kUnknown;
      }
    }
  }

  // Allocate memory for the histograms
  if (numHistoMacros > 0){
    histos = new TH1*[numHistoMacros];
    memset(histos, 0, numHistoMacros*sizeof(TH1*));
  }
  Bool_t secondBeforeFirstObject = kTRUE;
  

  //////////////////////////////////////
  // WALK THROUGH THE LIST OF OBJECTS //
  //////////////////////////////////////     
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter){
    if(!(object1 = dynamic_cast<TEveElement*>(*iter))) continue;

    // Skip objects that have not been selected
    if (!object1->GetRnrState())  continue;
    
    // Cast to the "real" object behind
    gROOT->ProcessLineSync(Form("TEveElement *automaticEveElement = (TEveElement*)%p;", (void*)object1));
    gROOT->ProcessLineSync("TObject* automaticObject_1 = (TObject*)automaticEveElement->GetUserData();");

    // Collect data for each macro
    for (Int_t i = 0, histoIndex = 0; i < procIterator->GetEntries(); i++){
      // Find the type of the object and relate it to the macro object type
      // Only apply macro to this object, if...
      // ... the macro takes objects of exactly this type.
      // ... the macro object type is a child of this object's type.
      // Otherwise: Continue

      // Finally, via procCmds[i], the automatic objects are casted to the correct type and analysed by each macro!
      if (((TObject*)object1->GetUserData())->IsA() != mProcObjectType[i] && 
          !((TObject*)object1->GetUserData())->InheritsFrom(mProcObjectType[i]))  continue;

        
      // Single object histo
      if (mProcType[i] == kSingleObjectHisto){
        returnedHist = (TH1*)gROOT->ProcessLineSync(procCmds[i]);
        if (histos && returnedHist)
        {
          if (!histos[histoIndex])  histos[histoIndex] = returnedHist;
          else  
          {
            histos[histoIndex]->Add((const TH1*)returnedHist);
            delete returnedHist;
            returnedHist = 0;
          }
        }
        histoIndex++;
       // Correlated Objects histo
      } else if (mProcType[i] == kCorrelObjectHisto) {
        // To get all pairs, do the second loop over all objects.
        // But: If a macro takes 2 pointers of the same type, we must take care that one gets the same pair, when we exchange the objects
        // (this is not true, if we have different types - even if they inherit from the same classes!).
        // Thus: If the latter case occurs, we ignore an object pair, if the second object is BEFORE the first object in the list.
        // Since then the pair has already been taken into account.
        // Furthermore, we have a pair of objects, if and only if both objects of the pair are selected (Rnr-state)
        // and are not equal.
        // The correlated objects process macro will be applied to all pairs that will be additionally selected by
        // all correlated objects selection macros.

        secondBeforeFirstObject = kTRUE;
        for (TEveElement::List_i iter2 = this->BeginChildren(); iter2 != this->EndChildren(); ++iter2)
        {
          // If the objects are the same, it is not a pair -> continue. From now on: 2nd object BEHIND the 1st object in the list!
          if (iter == iter2)
          {
            secondBeforeFirstObject = kFALSE;
            continue;
          }
          if(!(object2 = dynamic_cast<TEveElement*>(*iter2))) continue;

          // Skip objects that have not been selected
          if (!object2->GetRnrState())  continue;

          // Same check of the macro object type as before
          if (((TObject*)object2->GetUserData())->IsA() != mProcObjectType2[i] && 
              !((TObject*)object2->GetUserData())->InheritsFrom(mProcObjectType2[i]))  continue;
          // Do not process object pairs twice
          if (secondBeforeFirstObject)
          {
            if (mProcObjectType[i] == mProcObjectType2[i]) continue;
          }
      
          // Cast to the "real" object behind
          gROOT->ProcessLineSync(Form("TEveElement *automaticEveElement = (TEveElement*)%p;", (void*)object2));
          gROOT->ProcessLineSync("TObject* automaticObject_2 = (TObject*)automaticEveElement->GetUserData();");

          // Select object by default (so it will be processed, if there are no correlated objects selection macros!)
          selectedByCorrSelMacro = kTRUE;
          for (Int_t j = 0; j < selEntries; j++){
            if (mSelType[j] == kCorrelObjectSelect){
          // Check, whether the macro can deal with both objects. If not, skip it.
          // Note: Again, via selCmds[i], the automatic objects are casted to the correct type!
          if (((TObject*)object1->GetUserData())->IsA() != mSelObjectType[j] && 
              !((TObject*)object1->GetUserData())->InheritsFrom(mSelObjectType[j]))  continue;
          if (((TObject*)object2->GetUserData())->IsA() != mSelObjectType2[j] && 
              !((TObject*)object2->GetUserData())->InheritsFrom(mSelObjectType2[j]))  continue;

              selectedByCorrSelMacro = (Bool_t)gROOT->ProcessLineSync(selCmds[j]);
              if (!selectedByCorrSelMacro)  break;
            }
          }       

          // If the pair has not been selected by the correlated objects selection macros, skip it!
          if (!selectedByCorrSelMacro) continue;

          returnedHist = (TH1*)gROOT->ProcessLineSync(procCmds[i]);
          if (returnedHist && numHistoMacros)
          {
            if (!histos[histoIndex])  histos[histoIndex] = returnedHist;
            else  
            {
              histos[histoIndex]->Add((const TH1*)returnedHist);

              delete returnedHist;
              returnedHist = 0;
            }
          }
        }
        histoIndex++;
      }
      // Single object analyse
      else if (mProcType[i] == kSingleObjectAnalyse) {
        // Create data pointers in CINT, execute the macro and get the data
        gROOT->ProcessLineSync("Double_t* results = 0;");
        gROOT->ProcessLineSync("Int_t n = 0;");
        gROOT->ProcessLineSync(procCmds[i]);
        Double_t* results = (Double_t*)gROOT->ProcessLineSync("results;");
        Int_t nResults = (Int_t)gROOT->ProcessLineSync("n;");
        
        if (results == 0) {
          Error("Apply macros", "Error reading data from macro \"%s\"", procIterator->At(i)->GetTitle());
          continue;
        }
        for (Int_t resInd = 0; resInd < nResults; resInd++){
          (*fDataTree) << Form("ObjectData%d", i) << Form("Macro%d=", i) << results[resInd] << (Char_t*)"\n";   
        }

        delete results;
        results = 0;
      }
      // Correlated objects analyse
      else if (mProcType[i] == kCorrelObjectAnalyse){
        // To get all pairs, do the second loop over all objects.
        // But: If a macro takes 2 pointers of the same type, we must take care that one gets the same pair, when we exchange the objects
        // (this is not true, if we have different types - even if they inherit from the same classes!).
        // Thus: If the latter case occurs, we ignore an object pair, if the second object is BEFORE the first object in the list.
        // Since then the pair has already been taken into account.
        // Furthermore, we have a pair of objects, if and only if both objects of the pair are selected (Rnr-state)
        // and are not equal.
        // The correlated objects process macro will be applied to all pairs that will be additionally selected by
        // all correlated objects selection macros.

        secondBeforeFirstObject = kTRUE;
        for (TEveElement::List_i iter2 = this->BeginChildren(); iter2 != this->EndChildren(); ++iter2)
        {
          // If the objects are the same, it is not a pair -> continue. From now on: 2nd object BEHIND the 1st object in the list!
          if (iter == iter2)
          {
            secondBeforeFirstObject = kFALSE;
            continue;
          }
          if(!(object2 = dynamic_cast<TEveElement*>(*iter2))) continue;
 
          // Skip objects that have not been selected
          if (!object2->GetRnrState())  continue;

          // Same check of the macro object type as before
          if (((TObject*)object2->GetUserData())->IsA() != mProcObjectType2[i] && 
              !((TObject*)object2->GetUserData())->InheritsFrom(mProcObjectType2[i]))  continue;
          // Do not process object pairs twice
          if (secondBeforeFirstObject)
          {
            if (mProcObjectType[i] == mProcObjectType2[i]) continue;
          }
    
          // Cast to the "real" object behind
          gROOT->ProcessLineSync(Form("TEveElement *automaticEveElement = (TEveElement*)%p;", (void*)object2));
          gROOT->ProcessLineSync("TObject* automaticObject_2 = (TObject*)automaticEveElement->GetUserData();");

          // Select object by default (so it will be processed, if there are no correlated objects selection macros!)
          selectedByCorrSelMacro = kTRUE;
          for (Int_t j = 0; j < selEntries; j++) {
            if (mSelType[j] == kCorrelObjectSelect) {
              // Check, whether the macro can deal with both objects. If not, skip it.
              // Note: Again, via selCmds[i], the automatic objects are casted to the correct type! 
              if (((TObject*)object1->GetUserData())->IsA() != mSelObjectType[j] && 
                  !((TObject*)object1->GetUserData())->InheritsFrom(mSelObjectType[j]))  continue;
              if (((TObject*)object2->GetUserData())->IsA() != mSelObjectType2[j] && 
                  !((TObject*)object2->GetUserData())->InheritsFrom(mSelObjectType2[j]))  continue;

              selectedByCorrSelMacro = (Bool_t)gROOT->ProcessLineSync(selCmds[j]);
              if (!selectedByCorrSelMacro)  break;
            }
          }       

          // If the pair has not been selected by the correlated objects selection macros, skip it!
          if (!selectedByCorrSelMacro) continue;
          
          // Create data pointers in CINT, execute the macro and get the data
          gROOT->ProcessLineSync("Double_t* results = 0;");
          gROOT->ProcessLineSync("Int_t n = 0;");
          gROOT->ProcessLineSync(procCmds[i]);
          Double_t* results = (Double_t*)gROOT->ProcessLineSync("results;");
          Int_t nResults = (Int_t)gROOT->ProcessLineSync("n;");
     
          if (results == 0) {
            Error("Apply macros", "Error reading data from macro \"%s\"", procIterator->At(i)->GetTitle());
            continue;
          }
          for (Int_t resInd = 0; resInd < nResults; resInd++) {
            (*fDataTree) << Form("ObjectData%d", i) << Form("Macro%d=", i) << results[resInd] << (Char_t*)"\n";   
          }

          delete results;
          results = 0;
        }
      }
    }
  }    

  for (Int_t i = 0, histoIndex = 0; i < procIterator->GetEntries() && histoIndex < numHistoMacros; i++) {
    if (mProcType[i] == kSingleObjectHisto || mProcType[i] == kCorrelObjectHisto) {
      // Might be empty (e.g. no objects have been selected)!
      if (histos[histoIndex]) {
        (*fDataTree) << Form("ObjectData%d", i) << Form("Macro%d=", i) << histos[histoIndex] << (Char_t*)"\n";
      }
      histoIndex++;
    }
  }

  if (fDataTree != 0) delete fDataTree;
  fDataTree = 0;

  if (procCmds != 0)  delete [] procCmds;
  procCmds = 0;
  if (mProcObjectType != 0) delete [] mProcObjectType;
  mProcObjectType = 0;
  if (mProcObjectType2 != 0) delete [] mProcObjectType2;
  mProcObjectType2 = 0;
  if (mProcType != 0)  delete [] mProcType;
  mProcType = 0;

  if (selCmds != 0)  delete [] selCmds;
  selCmds = 0;
  if (mSelObjectType != 0)  delete [] mSelObjectType;
  mSelObjectType = 0;
  if (mSelObjectType2 != 0)  delete [] mSelObjectType2;
  mSelObjectType2 = 0;
  if (mSelType != 0)  delete [] mSelType;
  mSelType = 0;

  if (histos != 0)  delete [] histos;
  histos = 0;

  // Clear root
  // A.B. gROOT->Reset();
  
  // If there is data, select the first data set
  if (procIterator->GetEntries() > 0) SETBIT(fHistoDataSelected, 0);

  // Now the data is stored in "/tmp/ListAnalyserMacroData_$USER.root"
  // The editor will access this file to display the data
  return kTRUE;
}

//______________________________________________________
void AliEveListAnalyser::ApplySOSelectionMacros(const TList* iterator)
{
  // Uses the iterator (for the selected selection macros) to apply the selected macros to the data.
  // The rnr-states of the objects are set according to the result of the macro calls (kTRUE, if all
  // macros return kTRUE for this object, otherwise: kFALSE).
  // "SO" stands for "single object". This means that only single object selection macros are applied.
  // Correlated objects selection macros will be used inside the call of ApplyProcessMacros(...)!

  TGeneralMacroData* macro = 0;
  AliEveListAnalyserMacroType macroType = kUnknown;
  TEveElement* object1 = 0;
  Bool_t selectedByMacro = kFALSE;

  // Clear root
  // A.B. gROOT->Reset();

  // Select all objecs at first. A object is then deselected, if at least one selection macro
  // returns kFALSE for this object.
  // Enable all objects (Note: EnableListElements(..) will call "ElementChanged", which will cause unforeseen behaviour!)
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter) ((TEveElement*)(*iter))->SetRnrState(kTRUE);
  SetRnrState(kTRUE);
  
  for (Int_t i = 0; i < iterator->GetEntries(); i++){
    macro = (TGeneralMacroData*)fMacroList->GetValue(iterator->At(i)->GetTitle());

    if (!macro){
      Error("Apply selection macros", 
            "Macro list is corrupted: Macro \"%s\" is not registered!", iterator->At(i)->GetTitle());
      continue;
    }

#ifdef AliEveListAnalyser_DEBUG
    printf("AliEveListAnalyser: Applying selection macro: %s\n", macro->GetName());
#endif
    
    // Determine macro type
    macroType = macro->GetType();

    // Single object select macro
    if (macroType == kSingleObjectSelect){
      // Walk through the list of objects
      for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter)
      {
        object1 = dynamic_cast<TEveElement*>(*iter);

        if (!object1) continue;

        // If the object has already been deselected, nothing is to do here
        if (!object1->GetRnrState()) continue;

        // Find the type of the object and relate it to the macro object type
        // Only apply macro to this object, if...
        // ... the macro takes objects of exactly this type.
        // ... the macro object type is a child of this object's type.
        // Otherwise: Continue
        if (((TObject*)object1->GetUserData())->IsA() != macro->GetObjectType() && 
            !((TObject*)object1->GetUserData())->InheritsFrom(macro->GetObjectType()))  continue;

        // Cast to the "real" object behind
        gROOT->ProcessLineSync(Form("TEveElement *automaticEveElement = (TEveElement*)%p;", (void*)object1));
        gROOT->ProcessLineSync("TObject* automaticObject_1 = (TObject*)automaticEveElement->GetUserData();");

        // GetCmd() will cast the automatic objects to the correct type for each macro!
        selectedByMacro = (Bool_t)gROOT->ProcessLineSync(macro->GetCmd());
        object1->SetRnrState(selectedByMacro && object1->GetRnrState());               
      }
    }
    // Correlated objects select macro
    else if (macroType == kCorrelObjectSelect){
      // Will be processed in ApplyProcessMacros(...)
      continue;
    } else {
      Error("Apply selection macros", 
        "Macro list corrupted: Macro \"%s/%s.C\" is not registered as a selection macro!",
        macro->GetPath(), macro->GetName());
    } 
  }

  // Clear root
  // A.B. gROOT->Reset();  
}

//______________________________________________________
TClass* AliEveListAnalyser::GetMacroObjectType(const Char_t* name, Int_t argNum) const
{
  // Returns the type of object (of argument argNum) the macro with name "name" is dealing with; 
  // e.g. if you have the signature:
  // void MyMacro(const AliTRDtrackV1* track, Double_t* &results, Int_t& nResults)
  // the call 'GetMacroObjectType("MyMacro")' yields the AliTRDtrackV1-class.
  // If the macro is not found (or there is an error), 0x0 is returned.

  if (argNum - 1 < 0) return 0x0;

  TFunction* f = gROOT->GetGlobalFunction(name, 0 , kTRUE);
  TMethodArg* m = 0;
  TList* list = 0;

  if (f)
  {
    list = f->GetListOfMethodArgs();
    
    if (!list->IsEmpty())
    {
      m = (TMethodArg*)list->At(argNum - 1);

      if (m)  return TClass::GetClass(m->GetTypeName());
    }
  }  

  // Error
  return 0x0;
}

//______________________________________________________
AliEveListAnalyser::AliEveListAnalyserMacroType AliEveListAnalyser::GetMacroType(const Char_t* name, const Char_t* objectType, 
                                                                                 const Char_t* objectType2, Bool_t UseList) const
{
  // Returns the type of the corresponding macro, that accepts pointers of the classes "objectType" (first pointer) and
  // objectType2" (second pointer) as parametres. 
  // If "UseList" is kTRUE, the type will be looked up in the internal list (very fast). But if this list
  // does not exist, you have to use kFALSE for this parameter. Then the type will be determined by the
  // prototype! NOTE: It is assumed that the macro has been compiled! If not, the return value is not
  // predictable, but normally will be kUnknown.
  // Note: AddMacro(Fast) will update the internal list and RemoveMacros respectively.

  AliEveListAnalyserMacroType type = kUnknown;

  TString* typeStr = 0;
  TString* typeStr2 = 0;
  
  if (objectType != 0) 
  {
    typeStr = new TString(objectType);
    // Remove white-spaces
    typeStr->ReplaceAll(" ", "");
  }
  else
  {
    typeStr = new TString("TObject");
  }
  if (objectType2 != 0) 
  {
    typeStr2 = new TString(objectType2);
    // Remove white-spaces
    typeStr2->ReplaceAll(" ", "");
  }
  else
  {
    typeStr2 = new TString("TObject");
  }

  TString* mangled1Str = new TString();
  TString* mangled2Str = new TString();
  TString* mangled3Str = new TString();
  TString* mangled4Str = new TString();
  TString* mangledArg1Str = new TString();
  TString* mangledArg2Str = new TString();

  // We want "const 'OBJECTTYPE'*"
  mangled1Str->Form("const %s*", typeStr->Data());

  // We want "const 'OBJECTTYPE'*, Double_t*&, Int_t&"
  mangled2Str->Form("const %s*, Double_t*&, Int_t&", typeStr->Data());

  // We want "const 'OBJECTTYPE'*, const 'OBJECTTYPE2'*"
  mangled3Str->Form("const %s*, const %s*", typeStr->Data(), typeStr2->Data());

  // We want "const 'OBJECTTYPE'*, const 'OBJECTTYPE2'*, Double_t*&, Int_t&"
  mangled4Str->Form("const %s*, const %s*, Double_t*&, Int_t&", typeStr->Data(), typeStr2->Data());

  // We want "oPconstsP'OBJECTTYPE'mUsP"
  mangledArg1Str->Form("oPconstsP%smUsP", typeStr->Data());

  // We want "cOconstsP'OBJECTTYPE2'mUsP"
  mangledArg2Str->Form("cOconstsP%smUsP", typeStr2->Data());
  
  // Re-do the check of the macro type
  if (!UseList){
    // Single object select macro or single object histo macro?
    TFunction* f = gROOT->GetGlobalFunctionWithPrototype(name, mangled1Str->Data(), kTRUE);

    if (f != 0x0)
    {
      // Some additional check (is the parameter EXACTLY of the desired type?)
      if (strstr(f->GetMangledName(), mangledArg1Str->Data()) != 0x0)
      {
        // Single object select macro?
        if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
        { 
          type = kSingleObjectSelect;     
        }
        // single object histo macro?
        else if (!strcmp(f->GetReturnTypeName(), "TH1*"))
        {
          type = kSingleObjectHisto;
        }
      }
    }
    // Single object analyse macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, mangled2Str->Data(), kTRUE)) 
             != 0x0)
    {
      if (!strcmp(f->GetReturnTypeName(), "void"))
      {
        // Some additional check (are the parameters EXACTLY of the desired type?)
        if (strstr(f->GetMangledName(), mangledArg1Str->Data()) != 0x0 &&
            strstr(f->GetMangledName(), "cODouble_tmUaNsP") != 0x0 &&
            strstr(f->GetMangledName(), "cOInt_taNsP") != 0x0)
        {
          type = kSingleObjectAnalyse;
        }
      }
    }    
    // Correlated objects select macro or correlated objects histo macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, mangled3Str->Data(), kTRUE)) 
             != 0x0)
    {
      // Some additional check (is the parameter EXACTLY of the desired type?)
      if (strstr(f->GetMangledName(), mangledArg1Str->Data()) != 0x0 &&
          strstr(f->GetMangledName(), mangledArg2Str->Data()) != 0x0)
      {
        // Correlated objects select macro?
        if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
        { 
          type = kCorrelObjectSelect;     
        }
        // Correlated objects histo macro?
        else if (!strcmp(f->GetReturnTypeName(), "TH1*"))
        {
          type = kCorrelObjectHisto;
        }
      }
    }    
    // Correlated objects analyse macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, mangled4Str->Data(), kTRUE)) 
             != 0x0)
    {
      if (!strcmp(f->GetReturnTypeName(), "void"))
      {
        // Some additional check (is the parameter EXACTLY of the desired type?)
        if (strstr(f->GetMangledName(), mangledArg1Str->Data()) != 0x0 &&
            strstr(f->GetMangledName(), mangledArg2Str->Data()) != 0x0 &&
            strstr(f->GetMangledName(), "cODouble_tmUaNsP") != 0x0 &&
            strstr(f->GetMangledName(), "cOInt_taNsP") != 0x0)
        {
          type = kCorrelObjectAnalyse;
        }
      }
    }    
  }
  // Use list to look up the macro type
  else
  {
    TGeneralMacroData* macro = 0;
    macro = (TGeneralMacroData*)fMacroList->GetValue(name);
    if (macro == 0)  return kUnknown; 
    
    type = macro->GetType();
    switch (type)
    {
      case kSingleObjectSelect:
      case kSingleObjectAnalyse:
      case kSingleObjectHisto:
      case kCorrelObjectSelect:
      case kCorrelObjectAnalyse:
      case kCorrelObjectHisto:      
        break;
    default:
      type = kUnknown;
      break;
    }
  }

  // Clean up
  if (mangled1Str != 0)
  {
    mangled1Str->Clear();
    delete mangled1Str;
    mangled1Str = 0;
  }
  if (mangled2Str != 0)
  {
    mangled2Str->Clear();
    delete mangled2Str;
    mangled2Str = 0;
  }
  if (mangled3Str != 0)
  {
    mangled3Str->Clear();
    delete mangled3Str;
    mangled3Str = 0;
  }
  if (mangled4Str != 0)
  {
    mangled4Str->Clear();
    delete mangled4Str;
    mangled4Str = 0;
  }
  if (mangledArg1Str != 0)
  {
    mangledArg1Str->Clear();
    delete mangledArg1Str;
    mangledArg1Str = 0;
  }
  if (mangledArg2Str != 0)
  {
    mangledArg2Str->Clear();
    delete mangledArg2Str;
    mangledArg2Str = 0;
  }

  typeStr->Clear();
  delete typeStr;
  typeStr = 0;

  typeStr2->Clear();
  delete typeStr2;
  typeStr2 = 0;


  return type;
}


/*
//______________________________________________________
void AliEveListAnalyser::RemovePrimSelectedObjects()
{
  // Removes the (primarily) selected objects in the viewer from the list (objects that are already in the list are ignored).
  // Hold the CTRL-key for multiple selection.

  TEveSelection* eveSel = gEve->GetSelection();

  if (!eveSel)
  {
    Error("AliEveListAnalyser::RemovePrimSelectedObjects", "Failed to get the selection!\n");
    return;
  }
  
  TEveElement* elem = 0x0;
  Bool_t changedSomething = kFALSE;

  for (TEveElement::List_i iter = eveSel->BeginChildren(); iter != eveSel->EndChildren(); ++iter)
  {
    if(!(elem = dynamic_cast<TEveElement*>(*iter))) continue;

    // Check, if element is already there. If so, remove it!
    if (this->HasChild(elem) && elem != this)
    {
      this->RemoveElement(elem);
      this->SetTitle(Form("Objects %d", this->NumChildren()));
      changedSomething = kTRUE;
    }
  }

  if (changedSomething) gEve->Redraw3D();
}
*/

//______________________________________________________
void AliEveListAnalyser::RemoveSelectedMacros(const TList* iterator) 
{
  // Uses the iterator (for the selected macros) to remove the selected macros from 
  // the corresponding list.
   
  TObject* key = 0;
  TPair*   entry = 0;
  for (Int_t i = 0; i < iterator->GetEntries(); i++)
  {
    entry = (TPair*)fMacroList->FindObject(iterator->At(i)->GetTitle());

    if (entry == 0)
    {
      Error("AliEveListAnalyser::RemoveSelectedMacros", "Macro \"%s\" not found in list!",
                                                                     iterator->At(i)->GetTitle());
      continue;
    }
    key = entry->Key();

    if (key == 0)   
    {
      Error("AliEveListAnalyser::RemoveSelectedMacros", "Key for macro \"%s\" not found in list!",
                                                                     iterator->At(i)->GetTitle());
      continue;
    }

    // Key and value will be deleted, too, since fMacroList is the owner of them
    Bool_t rem = fMacroList->DeleteEntry(key);

    if (rem)
    {
#ifdef AliEveListAnalyser_DEBUG
    printf("AliEveListAnalyser::RemoveSelectedMacros(): Removed macro: %s\n", iterator->At(i)->GetTitle());
#endif
    }
    else
    {
      Error("AliEveListAnalyser::RemoveSelectedMacros", "Macro \"%s\" could not be removed from the list!",
                                                                     iterator->At(i)->GetTitle());
    }
  }
}

//______________________________________________________
void AliEveListAnalyser::ResetObjectList()
{
  // Removes all objects from the list.

  RemoveElements();
  this->SetTitle(Form("Objects %d", this->NumChildren()));
}

//______________________________________________________
Bool_t AliEveListAnalyser::StartAddingObjects()
{ 
  // Starts adding objects for the analysis. Returns kTRUE on success.

  if (fConnected == kFALSE)
  {
    fConnected = TQObject::Connect("TEvePointSet", "PointSelected(Int_t)", "AliEveListAnalyser", this, "AddSecSelectedSingleObjectToList(Int_t)");
    if (fConnected)  fConnected = TQObject::Connect(gEve->GetSelection(), "SelectionAdded(TEveElement*)", "AliEveListAnalyser", this, "AddPrimSelectedObject(TEveElement*)");

    if (fConnected) return kTRUE;
    
    Error("AliEveListAnalyser::StartAddingObjects", "Connection failed!");
    
    // Connection of 2nd signal failed, but first connection succeeded -> Disconnect 1st signal.
    TQObject::Disconnect("TEvePointSet", "PointSelected(Int_t)", this, "AddObjectToList(Int_t)");
  }

  return kFALSE;
}

//______________________________________________________
Bool_t AliEveListAnalyser::StopAddingObjects()
{
  // Stops adding objects for the analysis. Returns kTRUE on success.

  if (fConnected)
  {
    Bool_t dis1 = kFALSE, dis2 = kFALSE;
    dis1 = TQObject::Disconnect("TEvePointSet", "PointSelected(Int_t)", this, "AddSecSelectedSingleObjectToList(Int_t)");
    dis2 = TQObject::Disconnect(gEve->GetSelection(), "SelectionAdded(TEveElement*)", this, "AddPrimSelectedObject(TEveElement*)");

    if (dis1 || dis2) fConnected = kFALSE;
    if (dis1 && dis2) return kTRUE;
    else
    {
      Error("AliEveListAnalyser::StopAddingObjects", "Disconnection failed!");

      return kFALSE;
    }
  }

  return kTRUE;
}
