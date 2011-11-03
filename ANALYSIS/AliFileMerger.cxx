/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//  marian.ivanov@cern.ch
//  Utilities for file merging.
//  Additional functionality on top of the standard TFileMerger:
//
//  1. Possibility to Set the reject/accept list.
//     1.a)  Only entries selected in accept list are merged. By default all entries are selected
//           use AddAccept 0 to specify your desired entry
//     1.b)  Entries selected in reject list are not merged. By default the reject list is empty.
//
//  2. syswatch.log is created diring mergin procedure. 
//     Memeory consumption - for reading and for merging can be monitored

//  RS: Changed merger to respect the structure of files being merged (directories, collections...)
//      Additional option: SetNoTrees (default false) to not merge any tree
//      The code mostly taken from root's hadd.cxx
/*
  Usage:
  // Libraries for all classes to be merged should be loaded before using the class
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libTPCcalib"); 
  TH1::AddDirectory(0);

  //Example usage starting from the input data list in text file:
  //
  AliFileMerger merger;
  merger.AddReject("esdFriend");
  merger.IterTXT("calib.list","CalibObjects.root",kFALSE);
  //

*/
//////////////////////////////////////////////////////////////////////////
 

#include <fstream>
#include <THashList.h>
#include <TChain.h>
#include <TKey.h>
#include <TH1.h>
#include <THStack.h>
#include "TSystem.h"
#include "TFile.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMethodCall.h"
#include "Riostream.h"
#include "AliSysInfo.h"
#include "AliFileMerger.h"

ClassImp(AliFileMerger)

////////////////////////////////////////////////////////////////////////

AliFileMerger::AliFileMerger():
  TNamed(),
  fRejectMask(0),
  fAcceptMask(0),
  fNoTrees(kFALSE)
{
  //
  // Default constructor
  //
}

//______________________________________________________________________

AliFileMerger::AliFileMerger(const char* name):
  TNamed(name,name),
  fRejectMask(0),
  fAcceptMask(0),
  fNoTrees(kFALSE)
{
  //
  // 
  //
}


void AliFileMerger::IterAlien(const char* outputDir, const char* outputFileName, const char* pattern, Bool_t dontOverwrite){

  //
  // Merge the files coming out of the calibration job
  // 
  TString outputFile(outputFileName);
  gSystem->ExpandPathName(outputFile);
  TString command;
  // looking for files to be merged in the output directory
  command = Form("find %s/ *%s", outputDir, pattern);
  printf("command: %s\n", command.Data());
  TGrid::Connect("alien://");
  TGridResult *res = gGrid->Command(command);
  if (!res) return;
  TIter nextmap(res);
  TMap *map = 0;
  // loop over the results
  TList sourcelist;  
  sourcelist.SetOwner(kTRUE);
  //
  while((map=(TMap*)nextmap())) {
    // getting the turl
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) {
      // Nothing found - skip this output
      delete res;
      break;
    } 
    printf("looking for file %s\n",(objs->GetString()).Data());
    AddFile(&sourcelist, (objs->GetString()).Data());;
  }
  printf("List of files to be merged:\n");
  sourcelist.Print();
  //  
  TFile* target = TFile::Open( outputFile.Data(), (dontOverwrite ? "CREATE":"RECREATE") );
  if (!target || target->IsZombie()) {
    cerr << "Error opening target file (does " << outputFileName << " exist?)." << endl;
    cerr << "Use force = kTRUE to re-creation of output file." << endl;
    return;
  }  
  MergeRootfile( target, &sourcelist);
  delete res;
  delete target;
}

void AliFileMerger::IterTXT( const char * fileList,  const char* outputFileName, Bool_t dontOverwrite){
  
  // Merge the files indicated in the list - fileList
  // ASCII file option example: 
  // find `pwd`/ | grep AliESDfriends_v1.root > calib.list
  
  // Open the input stream
  ifstream in;
  in.open(fileList);
  // Read the input list of files 
  TString objfile;
  Int_t counter=0;
  TList sourcelist;  
  sourcelist.SetOwner(kTRUE);
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains(".root")) continue; // protection
    gSystem->ExpandPathName(objfile);
    printf("Add file:Counter\t%d\tMerging file %s\n",counter++,objfile.Data());
    AddFile(&sourcelist, objfile.Data());
  }
  //
  printf("List of files to be merged:\n");
  sourcelist.Print();
  //
  TString outputFile(outputFileName);
  gSystem->ExpandPathName(outputFile);
  TFile* target = TFile::Open( outputFile.Data(), (dontOverwrite ? "CREATE":"RECREATE") );
  if (!target || target->IsZombie()) {
    cerr << "Error opening target file (does " << outputFileName << " exist?)." << endl;
    cerr << "Use force = kTRUE to re-creation of output file." << endl;
    return;
  }

  MergeRootfile( target, &sourcelist);
  delete target;
}

void AliFileMerger::StoreResults(TObjArray * array, const char* outputFileName){
  //
  // Storing the results in one single file
  //
  TFile *f = new TFile(outputFileName,"recreate");
  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    object0->Write();
  }
  f->Close();
  delete f;
}


void AliFileMerger::StoreSeparateResults(TObjArray * array, const char* outputFileName){
  //
  // Store the results in separate files (one per object)
  //
  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    TFile *f = new TFile(Form("%s_%s.root",outputFileName,object0->GetName()),"recreate");
    object0->Write();
    f->Close();
    delete f;
  }
}

void AliFileMerger::Merge(TFile* fileIn, TObjArray * array){
  //
  // Merging procedure
  //
  if (!array) return;
  static Int_t counter=-1;
  counter++;
  TObjArray *carray = new TObjArray;   //array of the objects inside current file
  carray->SetOwner(kTRUE);
  
  // load all objects to  memory
  
  TList *farr = fileIn->GetListOfKeys();
  if (!farr) { 
    delete carray;
    return;
  }
  for (Int_t ical=0; ical<farr->GetEntries(); ical++){
    if (!farr->At(ical)) continue;
    TString name(farr->At(ical)->GetName());
    if (!IsAccepted(name)) continue;                        // skip not accepted entries
    TObject *obj = fileIn->Get(name.Data());
    if (obj) carray->AddLast(obj);
    AliSysInfo::AddStamp(name.Data(),1,ical,counter);  
  }
  
  if (carray->GetEntries()==0)  { 
    delete carray;
    return;
  }
  TMethodCall callEnv;
  Int_t entries =carray->GetEntriesFast();
  for (Int_t i=0; i<entries; i++){
    
    TObjArray *templist = new TObjArray(1);
    templist->SetOwner(kFALSE);
    TObject *currentObject = carray->At(i);
    if (!currentObject) { 
      delete templist;
      continue;
    }
    printf("%s\n",currentObject->GetName());
    callEnv.InitWithPrototype(currentObject->IsA(), "Merge", "TCollection*");
    if (!callEnv.IsValid()) {
      delete templist; 
      continue;
    }
    TString oname=currentObject->GetName();
    TObject *mergedObject = array->FindObject(currentObject->GetName());
    if (!mergedObject) {
      array->AddLast(currentObject);
      carray->RemoveAt(i);
      delete templist; 
      continue;
    }
    templist->AddLast(currentObject);
    callEnv.SetParam((Long_t) templist);
    callEnv.Execute(mergedObject);
    AliSysInfo::AddStamp(currentObject->GetName(),2,i,counter);  
    delete templist;
  }
  carray->Delete();
  delete carray;
}

Bool_t AliFileMerger::IsAccepted(TString name){
  //
  // Accept/reject logic
  // name - name of the entry
  //
  //  if fAcceptMask specified   - entry has to be in list of selected
  //  if fRejectMask speciefied  - entry with name speciief in the list are rejected 
  //
  Bool_t accept=kTRUE;
  if (fAcceptMask){
    //
    accept=kFALSE;
    for (Int_t iaccept=0; iaccept<fAcceptMask->GetEntries(); iaccept++){
      if (name.Contains(fAcceptMask->At(iaccept)->GetName())) accept=kTRUE;   // entry was selected
    }
  }
  if (!accept) return kFALSE;

  if (fRejectMask){
    //
    for (Int_t ireject=0; ireject<fRejectMask->GetEntries(); ireject++){
      if (name.Contains(fRejectMask->At(ireject)->GetName())) accept=kFALSE;   // entry was rejected
    }
  }
  return accept;
}




void AliFileMerger::AddReject(const char *reject){
  //
  // add reject string to the list of entries to be rejected for merging
  //
  if (!fRejectMask) fRejectMask = new TObjArray;
  fRejectMask->AddLast(new TObjString(reject));
}
void AliFileMerger::AddAccept(const char *accept){
  //
  // add reject string to the list of entries to be rejected for merging
  //
  if (!fAcceptMask) fAcceptMask = new TObjArray;
  fAcceptMask->AddLast(new TObjString(accept));


}

//___________________________________________________________________________
int AliFileMerger::MergeRootfile( TDirectory *target, TList *sourcelist)
{
  // Merge all objects in a directory
  // modified version of root's hadd.cxx
  Int_t counterF = -1;
  int status = 0;
  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  //
  // find 1st valid file
  TDirectory *first_source = (TDirectory*)sourcelist->First();
  //
  Int_t nguess = sourcelist->GetSize()+1000;
  THashList allNames(nguess);
  ((THashList*)target->GetList())->Rehash(nguess);
  ((THashList*)target->GetListOfKeys())->Rehash(nguess);
  TList listH;
  TString listHargs;
  listHargs.Form("((TCollection*)0x%lx)", (ULong_t)&listH);
  //
  while(first_source) {
    counterF++;
    TDirectory *current_sourcedir = first_source->GetDirectory(path);
    if (!current_sourcedir) {
      first_source = (TDirectory*)sourcelist->After(first_source);
      continue;
    }
    // loop over all keys in this directory
    TChain *globChain = 0;
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    //gain time, do not add the objects in the list in memory
    TH1::AddDirectory(kFALSE);
    //
    int counterK = 0;
    //
    while ( (key = (TKey*)nextkey())) {
      if (current_sourcedir == target) break;
      //
      // check if we don't reject this name
      TString nameK(key->GetName());
      if (!IsAccepted(nameK)) {
	if (!counterF) printf("Object %s is in rejection list, skipping...\n",nameK.Data());
	continue;
      }
      //
      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      if (!strcmp(key->GetClassName(),"TProcessID")) {key->ReadObj(); continue;}
      if (allNames.FindObject(key->GetName())) continue;
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (!cl || !cl->InheritsFrom(TObject::Class())) {
	cout << "Cannot merge object type, name: "
	     << key->GetName() << " title: " << key->GetTitle() << endl;
	continue;
      }
      allNames.Add(new TObjString(key->GetName()));
      AliSysInfo::AddStamp(nameK.Data(),1,counterK++,counterF-1); 
      // read object from first source file
      //current_sourcedir->cd();
      TObject *obj = key->ReadObj();
      //printf("keyname=%s, obj=%x\n",key->GetName(),obj);
      
      if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {
	
	// loop over all source files create a chain of Trees "globChain"
	if (!fNoTrees) { // 
	  TString obj_name;
	  if (path.Length()) {
	    obj_name = path + "/" + obj->GetName();
	  } else {
	    obj_name = obj->GetName();
	  }
	  globChain = new TChain(obj_name);
	  globChain->Add(first_source->GetName());
	  TFile *nextsource = (TFile*)sourcelist->After( first_source );
	  while ( nextsource ) {
	    //do not add to the list a file that does not contain this Tree
	    TFile *curf = TFile::Open(nextsource->GetName());
	    if (curf) {
	      Bool_t mustAdd = kFALSE;
	      if (curf->FindKey(obj_name)) {
		mustAdd = kTRUE;
	      } else {
		//we could be more clever here. No need to import the object
		//we are missing a function in TDirectory
		TObject *aobj = curf->Get(obj_name);
		if (aobj) { mustAdd = kTRUE; delete aobj;}
	      }
	      if (mustAdd) {
		globChain->Add(nextsource->GetName());
	      }
	    }
	    delete curf;
	    nextsource = (TFile*)sourcelist->After( nextsource );
	  }
	}
      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
	// it's a subdirectory
	
	cout << "Found subdirectory " << obj->GetName() << endl;
	// create a new subdir of same name and title in the target file
	target->cd();
	TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
	
	// newdir is now the starting point of another round of merging
	// newdir still knows its depth within the target file via
	// GetPath(), so we can still figure out where we are in the recursion
	status = MergeRootfile( newdir, sourcelist);
	if (status) return status;
	
      } else if ( obj->InheritsFrom(TObject::Class())
		  && obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*") ) {
	// object implements Merge(TCollection*)
	
	// loop over all source files and merge same-name object
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while ( nextsource ) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	    if (key2) {
	      TObject *hobj = key2->ReadObj();
	      hobj->ResetBit(kMustCleanup);
	      listH.Add(hobj);
	      Int_t error = 0;
	      obj->Execute("Merge", listHargs.Data(), &error);
	      if (error) {
		cerr << "Error calling Merge() on " << obj->GetName()
		     << " with the corresponding object in " << nextsource->GetName() << endl;
	      }
	      listH.Delete();
	    }
	  }
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      } else if ( obj->IsA()->InheritsFrom( THStack::Class() ) ) {
	THStack *hstack1 = (THStack*) obj;
	TList* l = new TList();
	
	// loop over all source files and merge the histos of the
	// corresponding THStacks with the one pointed to by "hstack1"
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while ( nextsource ) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(hstack1->GetName());
	    if (key2) {
	      THStack *hstack2 = (THStack*) key2->ReadObj();
	      l->Add(hstack2->GetHists()->Clone());
	      delete hstack2;
	    }
	  }
	  
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
	hstack1->GetHists()->Merge(l);
	l->Delete();
      } else {
	// object is of no type that we can merge
	cout << "Cannot merge object type, name: "
	     << obj->GetName() << " title: " << obj->GetTitle() << endl;
	
	// loop over all source files and write similar objects directly to the output file
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while ( nextsource ) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	    if (key2) {
	      TObject *nobj = key2->ReadObj();
	      nobj->ResetBit(kMustCleanup);
	      int nbytes1 = target->WriteTObject(nobj, key2->GetName(), "SingleKey" );
	      if (nbytes1 <= 0) status = -1;
	      delete nobj;
	    }
	  }
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      }
      
      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      target->cd();
      
      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( TDirectory::Class() )) {
	//printf("cas d'une directory\n");
      } else if(obj->IsA()->InheritsFrom( TTree::Class() )) {
	if (!fNoTrees) {
	  globChain->ls();
	  globChain->Merge(target->GetFile(),0,"keep fast");
	  delete globChain;
	}
      } else {
	int nbytes2 = obj->Write( key->GetName(), TObject::kSingleKey );
	if (nbytes2 <= 0) status = -1;
      }
      oldkey = key;
      delete obj;
    } // while ( ( TKey *key = (TKey*)nextkey() ) )
    first_source = (TDirectory*)sourcelist->After(first_source);
  }
  // save modifications to target file
  target->SaveSelf(kTRUE);
  //
  return status;
}

//___________________________________________________________________________
int AliFileMerger::AddFile(TList* sourcelist, std::string entry)
{
  // add a new file to the list of files
  //  static int count(0);
  if( entry.empty() ) return 0;
  size_t j =entry.find_first_not_of(' ');
  if( j==std::string::npos ) return 0;
  entry = entry.substr(j);
  if( entry.substr(0,1)=="@") {
    std::ifstream indirect_file(entry.substr(1).c_str() );
    if( ! indirect_file.is_open() ) {
      std::cerr<< "Could not open indirect file " << entry.substr(1) << std::endl;
      return 1;
    }
    while( indirect_file ){
      std::string line;
      std::getline(indirect_file, line);
      if( AddFile(sourcelist, line)!=0 )return 1;;
    }
    return 0;
  }
  //  cout << "Source file " << (++count) << ": " << entry << endl;
  
  TFile* source = TFile::Open( entry.c_str());
  if( source==0 ) {
    cout << "Failed to open " << entry << " will skip" << endl;
    return 0;
  }
  sourcelist->Add(source);
  return 0;
}
