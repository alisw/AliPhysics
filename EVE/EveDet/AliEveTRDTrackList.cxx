// Author: Benjamin Hess   29/01/2010

/*************************************************************************
 * Copyright (C) 2009-2010, Alexandru Bercuci, Benjamin Hess.            *
 * All rights reserved.                                                  *
 *************************************************************************/


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


// Uncomment to display debugging infos
//#define ALIEVETRDTRACKLIST_DEBUG

#include <TFile.h>
#include <TFunction.h>
#include <TH1.h>
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeStream.h>
#include <TMethodCall.h>

#include <AliTRDReconstructor.h>

#include <EveDet/AliEveTRDTrackList.h>
#include <EveDet/AliEveTRDTrackListEditor.h>

#include <../PWG1/TRD/AliTRDrecoTask.h>
#include <../PWG1/TRD/AliTRDpwg1Helper.h>

ClassImp(AliEveTRDTrackList)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrackList ////////////////////////
///////////////////////////////////////////////////////////
AliEveTRDTrackList::AliEveTRDTrackList(const Text_t* n, const Text_t* t, Bool_t doColor):
  TEveElementList(n, t, doColor),
  fEditor(0x0),
  fDataFromMacroList(0x0),
  fMacroList(0x0),
  fDataTree(0x0),
  fHistoDataSelected(0),
  fMacroListSelected(0),
  fSelectedTab(1),                              // Standard tab: "Apply macros" (index 1)
  fSelectedStyle(0)
{
  // Creates the AliEveTRDTrackList.

  // Only accept childs of type AliEveTRDTrack
  SetChildClass(AliEveTRDTrack::Class());

  // Allocate memory for the lists and declare them as owners of their contents
  fDataFromMacroList = new TList();
  fDataFromMacroList->TCollection::SetOwner(kTRUE);

  fMacroList = new TMap();
  // Set map to owner of it's objects to delete them, if they are removed from the map
  fMacroList->SetOwnerKeyValue(kTRUE, kTRUE);

  // Set the build directory for AClic
  if(gSystem->AccessPathName(Form("%s/.trdQArec" , gSystem->Getenv("HOME")))) gSystem->Exec("mkdir $HOME/.trdQArec");
  gSystem->SetBuildDir(Form("%s/.trdQArec", gSystem->Getenv("HOME")));

  AddStandardContent();
}

//______________________________________________________
AliEveTRDTrackList::~AliEveTRDTrackList()
{
  // Frees allocated memory (lists etc.).

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
  if(!gSystem->AccessPathName(Form("/tmp/TRD.TrackListMacroData_%s.root", gSystem->Getenv("USER")))) 
    gSystem->Exec(Form("rm /tmp/TRD.TrackListMacroData_%s.root", gSystem->Getenv("USER")));
}

//______________________________________________________
Int_t AliEveTRDTrackList::AddMacro(const Char_t* path, const Char_t* nameC, Bool_t forceReload)
{
  // Checks, if the file exists and if the signature is correct.
  // If these criteria are fullfilled, the library for this macro is built
  // and the macro is added to the corresponding list.
  // Supported macro types:
  // Selection macros:                                                    
  // Bool_t YourMacro(const AliTRDtrackV1*)
  // Bool_t YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*)
  //
  // Process macros:                                                      
  // void YourMacro(const AliTRDtrackV1*, Double_t*&, Int_t&)             
  // void YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*,           
  //                Double_t*&, Int_t&)                                   
  // TH1* YourMacro(const AliTRDtrackV1*)                                 
  // TH1* YourMacro(const AliTRDtrackV1*, const AliTRDtrackV1*)                              

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

  // We need this line... otherwise, in some cases, there will be problems concerning ACLIC
  gROOT->ProcessLineSync(Form(".L %s", pathname));

  AliEveTRDTrackListMacroType type = GetMacroType(name, kFALSE);

  // Clean up again
  // A.B. gROOT->Reset();
  
  // Has not the correct signature!
  if (type == kUnknown)  return SIGNATURE_ERROR;

  // Only add macro, if it is not already in the list
  Int_t returnValue = WARNING;
  if(fMacroList->GetValue(name) == 0) {
    returnValue = AddMacroFast(path, name, type) ? SUCCESS : ERROR;
  }
  return returnValue;
}

//______________________________________________________
Bool_t AliEveTRDTrackList::AddMacroFast(const Char_t* path, const Char_t* name, AliEveTRDTrackListMacroType type)
{
  // Adds a macro (path/name) to the corresponding list. No checks are performed (file exist, 
  // macro already in list/map, signature correct),  no libraries are created!
  // You can use this function only, if the macro has been added successfully before 
  // (and then maybe was removed). The function is very fast. On success kTRUE is returned, otherwise: kFALSE;

  Bool_t success = kFALSE;

  switch (type)
  {
    case kSingleTrackSelect:
    case kCorrelTrackSelect:
    case kSingleTrackAnalyse:
    case kSingleTrackHisto:
    case kCorrelTrackAnalyse:
    case kCorrelTrackHisto:
      fMacroList->Add(new TObjString(name), new TMacroData(name, path, type));

      // We do not know, where the element has been inserted - deselect this list
      fMacroListSelected = 0;
    
      success = kTRUE;

#ifdef ALIEVETRDTRACKLIST_DEBUG
      // Successfull add will only be displayed in debug mode
      printf("AliEveTRDTrackList::AddMacroFast: Added macro \"%s/%s\" to the corresponding list\n", path, name);
#endif

      break;

    default:
      // Error will always be displayed
      printf("AliEveTRDTrackList::AddMacroFast: ERROR: Could not add macro \"%s/%s\" to the corresponding list\n", 
             path, name);

      success = kFALSE;

      break;
  }

  return success;
}

//______________________________________________________
void AliEveTRDTrackList::AddStandardContent()
{
  // Adds standard macros to the macro list.

  // Add your standard macros here, e.g.:
  // To add a macro use:
  // AddMacro("$(ALICE_ROOT)/myFolder", "myMacroName.C");
  // -> If the file does not exist, nothing happens. So if you want to handle this,
  // use the return value of AddMacro (NOT_EXIST_ERROR is returned, if file does not exist)
  // (-> You can also check for other return values (see AddMacro(...)))

  const Char_t *libs[] = {"libANALYSIS.so", "libANALYSISalice.so", "libTENDER.so", "libPWG1.so"};
  Int_t nlibs = static_cast<Int_t>(sizeof(libs)/sizeof(Char_t *));
  for(Int_t ilib=0; ilib<nlibs; ilib++){
    if(gSystem->Load(libs[ilib]) >= 0) continue;
    AliError(Form("Fail loading %s.", libs[ilib]));
    return;
  }

  const Char_t *fgkTRDPWG1taskClassName[AliTRDpwg1Helper::kNTRDQATASKS] = {
    "AliTRDcheckESD"
    ,"AliTRDinfoGen"
    ,"AliTRDcheckDET"
    ,"AliTRDefficiency"
    ,"AliTRDresolution"
    ,"AliTRDcheckPID"
    ,"AliTRDv0Monitor"
  };
  AliTRDrecoTask *task(NULL);
  TList *fPlots(NULL);
  for(Int_t it=2; it<AliTRDpwg1Helper::kNTRDQATASKS; it++){
    TClass c(fgkTRDPWG1taskClassName[it]);
    task = (AliTRDrecoTask*)c.New();
    task->SetMCdata(kFALSE);
    if(!(fPlots = task->GetPlotFunctors())){
      //AliWarning(Form("No Plot functors defined for task \"%s\"", fgkTRDtaskClassName[it]));
      delete task;
      continue;
    }
    if(!(task->Histos())){
      //AliWarning(Form("No Ref Histograms defined for task \"%s\"", fgkTRDtaskClassName[it]));
      delete task;
      continue;
    }

    // export task to CINT and add functions
    gROOT->ProcessLine(Form("%s* %s = (%s*)%p;", fgkTRDPWG1taskClassName[it], task->GetName(), fgkTRDPWG1taskClassName[it], (void*)task));
    TIter iter(fPlots); TMethodCall *m = 0x0;
    while((m = dynamic_cast<TMethodCall*>(iter()))){
      AddMacroFast("", Form("%s->%s", task->GetName(), m->GetMethodName()), kSingleTrackHisto);
    }
  }
}


//______________________________________________________
Bool_t AliEveTRDTrackList::ApplyProcessMacros(const TList* selIterator, const TList* procIterator)
{
  // Uses the procIterator (for the selected process macros) to apply the selected macros to the data.
  // Returns kTRUE on success, otherwise kFALSE. If there no process macros selected, kTRUE is returned 
  // (this is no error!).
  // The single track process macros are applied to all selected tracks.
  // The selIterator (for the selected selection macros) will be used to apply the correlated tracks selection
  // macros to all track pairs (whereby BOTH tracks have to be selected, otherwise they will be skipped).
  // All track pairs that have been selected by ALL correlated tracks selection macros will be processed by
  // the correlated tracks process macros.

  // No process macros need to be processed
  if (procIterator->GetEntries() <= 0)  return kTRUE;

  // Clear root
  // A.B. gROOT->Reset();
  
  // Clear old data and re-allocate
  if (!fDataTree){
    TDirectory *cwd = gDirectory;
    fDataTree = new TTreeSRedirector(Form("/tmp/TRD.TrackListMacroData_%s.root", gSystem->Getenv("USER")));
    cwd->cd();
  }
  if (!fDataTree){
    Error("Apply process macros", "File \"/tmp/TRD.TrackListMacroData_%s.root\" could not be accessed properly!", gSystem->Getenv("USER"));
    return kFALSE;
  }
  
  if (fDataFromMacroList != 0) {
    fDataFromMacroList->Delete();
    delete fDataFromMacroList;
  }
  fDataFromMacroList = new TList();
  fDataFromMacroList->TCollection::SetOwner(kTRUE);

  fHistoDataSelected = 0;


  TMacroData* macro = 0;

  TString* procCmds = new TString[procIterator->GetEntries()];
  AliEveTRDTrackListMacroType* mProcType = new AliEveTRDTrackListMacroType[procIterator->GetEntries()];

  TString* selCmds(NULL);
  AliEveTRDTrackListMacroType* mSelType(NULL);
  if (selIterator->GetEntries() > 0) {
    selCmds = new TString[selIterator->GetEntries()];
    mSelType = new AliEveTRDTrackListMacroType[selIterator->GetEntries()];
  }
  
  Bool_t selectedByCorrSelMacro = kFALSE;

  AliEveTRDTrackListMacroType macroType = kUnknown;
  Int_t numHistoMacros = 0;
  TH1** histos(NULL);

  AliEveTRDTrack* track1(NULL);
  AliEveTRDTrack* track2(NULL);

  // Collect the commands for each process macro and add them to "data-from-list"
  for (Int_t i = 0; i < procIterator->GetEntries(); i++){
    macro = (TMacroData*)fMacroList->GetValue(procIterator->At(i)->GetTitle());

    if (!macro){
      Error("Apply process macros", 
        "Macro list is corrupted: Macro \"%s\" is not registered!",
        procIterator->At(i)->GetTitle());
      continue;
    }

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Checking process macro: %s\n", macro->GetName());
#endif 
           
    // Find the type of the process macro
    macroType = macro->GetType();
    if (macroType == kSingleTrackHisto || macroType == kCorrelTrackHisto){
      mProcType[i] = macroType;
      numHistoMacros++;
      // Create the command 
      procCmds[i] = macro->GetCmd();

      // Add to "data-from-list" -> Mark as a histo macro with the substring "(histo macro)"
      fDataFromMacroList->Add(new TObjString(Form("%s (histo macro)", macro->GetName())));
    } else if (macroType == kSingleTrackAnalyse || macroType == kCorrelTrackAnalyse) {
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

  // Collect the commands for each selection macro and add them to "data-from-list"
  for (Int_t i = 0; i < selIterator->GetEntries(); i++){
    macro = (TMacroData*)fMacroList->GetValue(selIterator->At(i)->GetTitle());

    if (!macro){
      Error("Apply process macros", 
        "Macro list is corrupted: Macro \"%s\" is not registered!",
        selIterator->At(i)->GetTitle());
      continue;
    }

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Checking selection macro: %s\n", macro->GetName());
#endif
       
    // Find the type of the process macro
    macroType = macro->GetType();

    // Single track select macro
    if (macroType == kSingleTrackSelect) {
      // Has already been processed by ApplySTSelectionMacros(...)
      if(mSelType) mSelType[i] = macroType;
    }
    // Correlated tracks select macro
    else if (macroType == kCorrelTrackSelect) {
      if(mSelType) mSelType[i] = macroType;  
 
      // Create the command
      if(selCmds) selCmds[i] = macro->GetCmd();
    } else {
      Error("Apply process macros", 
        "Macro list corrupted: Macro \"%s/%s.C\" is not registered as a selection macro!",
        macro->GetPath(), macro->GetName());
      if(mSelType) mSelType[i] = kUnknown;
    } 
  }  

  // Allocate memory for the histograms
  if (numHistoMacros > 0){
    histos = new TH1*[numHistoMacros];
    memset(histos, 0, numHistoMacros*sizeof(TH1*));
  }

  //////////////////////////////////
  // WALK THROUGH THE LIST OF TRACKS
  //////////////////////////////////     
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter){
    if(!(track1 = dynamic_cast<AliEveTRDTrack*>(*iter))) continue;

    // Skip tracks that have not been selected
    if (!track1->GetRnrState())  continue;
    
    // Cast to AliTRDtrackV1
    gROOT->ProcessLineSync(Form("AliEveTRDTrack *automaticTrack = (AliEveTRDTrack*)%p;", (void*)track1));
    gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1_1 = (AliTRDtrackV1*)automaticTrack->GetUserData();");

    // Collect data for each macro
    for (Int_t i = 0, histoIndex = 0; i < procIterator->GetEntries(); i++){
      // Single track histo
      if (mProcType[i] == kSingleTrackHisto){
        if(histos) histos[histoIndex++] = (TH1*)gROOT->ProcessLineSync(procCmds[i]);
       // Correlated tracks histo
      } else if (mProcType[i] == kCorrelTrackHisto) {
        // Loop over all pairs behind the current one - together with the other loop this will be a loop
        // over all pairs. We have a pair of tracks, if and only if both tracks of the pair are selected (Rnr-state)
        // and are not equal.
        // The correlated tracks process macro will be applied to all pairs that will be additionally selected by
        // all correlated tracks selection macros.
        TEveElement::List_i iter2 = iter;
        iter2++;
        for ( ; iter2 != this->EndChildren(); ++iter2)
        {
          if(!(track2 = dynamic_cast<AliEveTRDTrack*>(*iter2))) continue;

          // Skip tracks that have not been selected
          if (!track2->GetRnrState())  continue;
      
          // Cast to AliTRDtrackV1
          gROOT->ProcessLineSync(Form("AliEveTRDTrack *automaticTrack = (AliEveTRDTrack*)%p;", (void*)track2));
          gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1_2 = (AliTRDtrackV1*)automaticTrack->GetUserData();");

          // Select track by default (so it will be processed, if there are no correlated tracks selection macros!)
          selectedByCorrSelMacro = kTRUE;
          for (Int_t j = 0; j < selIterator->GetEntries(); j++){
            if (mSelType && mSelType[j] == kCorrelTrackSelect){
              selectedByCorrSelMacro = (Bool_t)gROOT->ProcessLineSync(selCmds[j]);
              if (!selectedByCorrSelMacro)  break;
            }
          }       

          // If the pair has not been selected by the correlated tracks selection macros, skip it!
          if (!selectedByCorrSelMacro) continue;
          
          if(histos) histos[histoIndex] = (TH1*)gROOT->ProcessLineSync(procCmds[i]);
        }
        histoIndex++;
      }
      // Single track analyse
      else if (mProcType[i] == kSingleTrackAnalyse) {
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
          (*fDataTree) << Form("TrackData%d", i) << Form("Macro%d=", i) << results[resInd] << (Char_t*)"\n";   
        }

        delete results;
        results = 0;
      }
      // Correlated tracks analyse
      else if (mProcType[i] == kCorrelTrackAnalyse){
        // Loop over all pairs behind the current one - together with the other loop this will be a loop
        // over all pairs. We have a pair of tracks, if and only if both tracks of the pair are selected (Rnr-state)
        // and are not equal.
        // The correlated tracks process macro will be applied to all pairs that will be additionally selected by
        // all correlated tracks selection macros.
        TEveElement::List_i iter2 = iter;
        iter2++;

        for ( ; iter2 != this->EndChildren(); ++iter2) {
          if(!(track2 = dynamic_cast<AliEveTRDTrack*>(*iter2))) continue;
 
          // Skip tracks that have not been selected
          if (!track2->GetRnrState())  continue;
    
          // Cast to AliTRDtrackV1
          gROOT->ProcessLineSync(Form("AliEveTRDTrack *automaticTrack = (AliEveTRDTrack*)%p;", (void*)track2));
          gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1_2 = (AliTRDtrackV1*)automaticTrack->GetUserData();");

          // Select track by default (so it will be processed, if there are no correlated tracks selection macros!)
          selectedByCorrSelMacro = kTRUE;
          for (Int_t j = 0; j < selIterator->GetEntries(); j++) {
            if (mSelType && mSelType[j] == kCorrelTrackSelect) {
              selectedByCorrSelMacro = (Bool_t)gROOT->ProcessLineSync(selCmds[j]);
              if (!selectedByCorrSelMacro)  break;
            }
          }       

          // If the pair has not been selected by the correlated tracks selection macros, skip it!
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
            (*fDataTree) << Form("TrackData%d", i) << Form("Macro%d=", i) << results[resInd] << (Char_t*)"\n";   
          }

          delete results;
          results = 0;
        }
      }
    }
  }    

  for (Int_t i = 0, histoIndex = 0; i < procIterator->GetEntries() && histoIndex < numHistoMacros; i++) {
    if (mProcType[i] == kSingleTrackHisto || mProcType[i] == kCorrelTrackHisto) {
      // Might be empty (e.g. no tracks have been selected)!
      if (histos[histoIndex]) {
        (*fDataTree) << Form("TrackData%d", i) << Form("Macro%d=", i) << histos[histoIndex] << (Char_t*)"\n";
      }
      histoIndex++;
    }
  }

  if (fDataTree) delete fDataTree;
  fDataTree = NULL;

  if (procCmds)  delete [] procCmds;
  procCmds = NULL;
  if (mProcType)  delete [] mProcType;
  mProcType = NULL;

  if (selCmds)  delete [] selCmds;
  selCmds = NULL;
  if (mSelType)  delete [] mSelType;
  mSelType = NULL;

  if (histos)  delete [] histos;
  histos = NULL;

  // Clear root
  // A.B. gROOT->Reset();
  
  // If there is data, select the first data set
  if (procIterator->GetEntries() > 0) SETBIT(fHistoDataSelected, 0);

  // Now the data is stored in "/tmp/TRD.TrackListMacroData_$USER.root"
  // The editor will access this file to display the data
  return kTRUE;
}

//______________________________________________________
void AliEveTRDTrackList::ApplySTSelectionMacros(const TList* iterator)
{
  // Uses the iterator (for the selected selection macros) to apply the selected macros to the data.
  // The rnr-states of the tracks are set according to the result of the macro calls (kTRUE, if all
  // macros return kTRUE for this track, otherwise: kFALSE).
  // "ST" stands for "single track". This means that only single track selection macros are applied.
  // Correlated tracks selection macros will be used inside the call of ApplyProcessMacros(...)!

  TMacroData* macro = 0;
  AliEveTRDTrackListMacroType macroType = kUnknown;
  AliEveTRDTrack* track1 = 0;
  Bool_t selectedByMacro = kFALSE;

  // Clear root
  // A.B. gROOT->Reset();

  // Select all tracks at first. A track is then deselected, if at least one selection macro
  // returns kFALSE for this track.
  // Enable all tracks (Note: EnableListElements(..) will call "ElementChanged", which will cause unforeseen behaviour!)
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter) ((TEveElement*)(*iter))->SetRnrState(kTRUE);
  SetRnrState(kTRUE);
  
  for (Int_t i = 0; i < iterator->GetEntries(); i++){
    macro = (TMacroData*)fMacroList->GetValue(iterator->At(i)->GetTitle());

    if (!macro){
      Error("Apply selection macros", 
            "Macro list is corrupted: Macro \"%s\" is not registered!", iterator->At(i)->GetTitle());
      continue;
    }

#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList: Applying selection macro: %s\n", macro->GetName());
#endif
    
    // Determine macro type
    macroType = macro->GetType();

    // Single track select macro
    if (macroType == kSingleTrackSelect){
      // Walk through the list of tracks
      for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter)
      {
        track1 = dynamic_cast<AliEveTRDTrack*>(*iter);

        if (!track1) continue;

        // If the track has already been deselected, nothing is to do here
        if (!track1->GetRnrState()) continue;

        // Cast to AliTRDtrackV1
        gROOT->ProcessLineSync(Form("AliEveTRDTrack *automaticTrack = (AliEveTRDTrack*)%p;", (void*)track1));
        gROOT->ProcessLineSync("AliTRDtrackV1* automaticTrackV1_1 = (AliTRDtrackV1*)automaticTrack->GetUserData();");
        selectedByMacro = (Bool_t)gROOT->ProcessLineSync(macro->GetCmd());
        track1->SetRnrState(selectedByMacro && track1->GetRnrState());               
      }
    }
    // Correlated tracks select macro
    else if (macroType == kCorrelTrackSelect){
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
AliEveTRDTrackList::AliEveTRDTrackListMacroType AliEveTRDTrackList::GetMacroType(const Char_t* name, Bool_t UseList) const
{
  // Returns the type of the corresponding macro. 
  // If "UseList" is kTRUE, the type will be looked up in the internal list (very fast). But if this list
  // does not exist, you have to use kFALSE for this parameter. Then the type will be determined by the
  // prototype! NOTE: It is assumed that the macro has been compiled! If not, the return value is not
  // predictable, but normally will be kUnknown.
  // Note: AddMacro(Fast) will update the internal list and RemoveMacros respectively.

  AliEveTRDTrackListMacroType type = kUnknown;

  // Re-do the check of the macro type
  if (!UseList){
    // Single track select macro or single track histo macro?
    TFunction* f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*", kTRUE);
    if (f != 0x0)
    {
      // Some additional check (is the parameter EXACTLY of the desired type?)
      if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0)
      {
        // Single track select macro?
        if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
        { 
          type = kSingleTrackSelect;     
        }
        // single track histo macro?
        else if (!strcmp(f->GetReturnTypeName(), "TH1*"))
        {
          type = kSingleTrackHisto;
        }
      }
    }
    // Single track analyse macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*, Double_t*&, Int_t&", kTRUE)) 
             != 0x0)
    {
      if (!strcmp(f->GetReturnTypeName(), "void"))
      {
        // Some additional check (are the parameters EXACTLY of the desired type?)
        if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0 &&
            strstr(f->GetMangledName(), "cODouble_tmUaNsP") != 0x0 &&
            strstr(f->GetMangledName(), "cOInt_taNsP") != 0x0)
        {
          type = kSingleTrackAnalyse;
        }
      }
    }    
    // Correlated tracks select macro or correlated tracks histo macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, "const AliTRDtrackV1*, const AliTRDtrackV1*", kTRUE)) 
             != 0x0)
    {
      // Some additional check (is the parameter EXACTLY of the desired type?)
      if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0 &&
          strstr(f->GetMangledName(), "cOconstsPAliTRDtrackV1mUsP") != 0x0)
      {
        // Correlated track select macro?
        if (!strcmp(f->GetReturnTypeName(), "Bool_t")) 
        { 
          type = kCorrelTrackSelect;     
        }
        // Correlated track histo macro?
        else if (!strcmp(f->GetReturnTypeName(), "TH1*"))
        {
          type = kCorrelTrackHisto;
        }
      }
    }    
    // Correlated tracks analyse macro?
    else if ((f = gROOT->GetGlobalFunctionWithPrototype(name, 
                              "const AliTRDtrackV1*, const AliTRDtrackV1*, Double_t*&, Int_t&", kTRUE)) 
             != 0x0)
    {
      if (!strcmp(f->GetReturnTypeName(), "void"))
      {
        // Some additional check (is the parameter EXACTLY of the desired type?)
        if (strstr(f->GetMangledName(), "oPconstsPAliTRDtrackV1mUsP") != 0x0 &&
            strstr(f->GetMangledName(), "cOconstsPAliTRDtrackV1mUsP") != 0x0 &&
            strstr(f->GetMangledName(), "cODouble_tmUaNsP") != 0x0 &&
            strstr(f->GetMangledName(), "cOInt_taNsP") != 0x0)
        {
          type = kCorrelTrackAnalyse;
        }
      }
    }    
  }
  // Use list to look up the macro type
  else
  {
    TMacroData* macro = 0;
    macro = (TMacroData*)fMacroList->GetValue(name);
    if (macro == 0)  return kUnknown; 
    
    type = macro->GetType();
    switch (type)
    {
      case kSingleTrackSelect:
      case kSingleTrackAnalyse:
      case kSingleTrackHisto:
      case kCorrelTrackSelect:
      case kCorrelTrackAnalyse:
      case kCorrelTrackHisto:      
        break;
    default:
      type = kUnknown;
      break;
    }
  }

  return type;
}

//______________________________________________________
void AliEveTRDTrackList::RemoveSelectedMacros(const TList* iterator) 
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
      Error("AliEveTRDTrackList::RemoveSelectedMacros", "Macro \"%s\" not found in list!",
                                                             iterator->At(i)->GetTitle());
      continue;
    }
    key = entry->Key();

    if (key == 0)   
    {
      Error("AliEveTRDTrackList::RemoveSelectedMacros", "Key for macro \"%s\" not found in list!",
                                                             iterator->At(i)->GetTitle());
      continue;
    }

    // Key and value will be deleted, too, since fMacroList is the owner of them
    Bool_t rem = fMacroList->DeleteEntry(key);

    if (rem)
    {
#ifdef ALIEVETRDTRACKLIST_DEBUG
    printf("AliEveTRDTrackList::RemoveSelectedMacros(): Removed macro: %s\n", iterator->At(i)->GetTitle());
#endif
    }
    else
    {
      Error("AliEveTRDTrackList::RemoveSelectedMacros", "Macro \"%s\" could not be removed from the list!",
                                                             iterator->At(i)->GetTitle());
    }
  }
}

//______________________________________________________
void AliEveTRDTrackList::UpdateTrackStyle(AliEveTRDTrack::AliEveTRDTrackState s, UChar_t ss)
{
  // Updates the track style and sets this style for each track.

  switch(s)
  {
    case AliEveTRDTrack::kSource:
      SETBIT(fSelectedStyle, AliEveTRDTrack::kSource);
      break;  
    case AliEveTRDTrack::kPID:
      CLRBIT(fSelectedStyle, AliEveTRDTrack::kSource);
      switch(ss)
      {
      case AliTRDpidUtil::kLQ:
        CLRBIT(fSelectedStyle, AliEveTRDTrack::kPID);
        break;
      case AliTRDpidUtil::kNN:
        SETBIT(fSelectedStyle, AliEveTRDTrack::kPID);
        break;
      }
      break;  
    case AliEveTRDTrack::kTrackCosmics:
      SETBIT(fSelectedStyle, AliEveTRDTrack::kTrackCosmics);
      break;  
    case AliEveTRDTrack::kTrackModel:
      CLRBIT(fSelectedStyle, AliEveTRDTrack::kTrackCosmics);
      switch(ss)
      {
      case AliEveTRDTrack::kRieman:
        CLRBIT(fSelectedStyle, AliEveTRDTrack::kTrackModel);
        break;
      case AliEveTRDTrack::kKalman:
        //AliWarning("Kalman fit under testing for the moment.");
        SETBIT(fSelectedStyle, AliEveTRDTrack::kTrackModel);
        break;
      }
      break;  
  }


  // Walk through the list of tracks     
  AliEveTRDTrack* track = 0x0;
  for (TEveElement::List_i iter = this->BeginChildren(); iter != this->EndChildren(); ++iter) 
  {
    if (!(track = dynamic_cast<AliEveTRDTrack*>(*iter)))  continue;

    track->SetStatus(fSelectedStyle);
  }
}
