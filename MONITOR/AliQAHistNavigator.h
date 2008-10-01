/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//     (see AliQAHistNavigator.cxx for details)
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIQAHISTNAVIGATOR_H
#define ALIQAHISTNAVIGATOR_H

#include "TSystem.h"
#include "Riostream.h"
#include "TH1D.h"
#include "TF1.h"
#include "TList.h"
#include "TObjString.h"
#include "TString.h"
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TKey.h"
#include "TText.h"
#include <list>
#include <string>

class AliQADirList;
class AliQADirListItem;
class AliQAHistNavigator;

class AliQADirList : public TNamed{
public:
    AliQADirList();
    virtual ~AliQADirList();
    TList* GetItems() {return fPItems;}
    TList* GetDirs() {return fPDirs;}
    AliQADirList* GetParent() {return fPParent;}
    void SetParent(AliQADirList* l) {fPParent = l;}
    
private:
    AliQADirList* fPParent;          //pointer to parent folder
    TList* fPItems;     //List of items contained in the list
    TList* fPDirs;      //List of dirs
    AliQADirList(const AliQADirList&);            // Not implemented
    AliQADirList& operator=(const AliQADirList&); // Not implemented

    ClassDef(AliQADirList,999)  //AliQADirListDir
};

class AliQADirListItem : public TObjString {
public:
    AliQADirListItem(const char* s="");
    virtual ~AliQADirListItem();
    AliQADirList* GetParent() {return fPParent;}
    void SetParent(AliQADirList* parent) {fPParent=parent;}

private:
    AliQADirList* fPParent;
    AliQADirListItem(const AliQADirListItem&);            // Not implemented
    AliQADirListItem& operator=(const AliQADirListItem&); // Not implemented

    ClassDef(AliQADirListItem,999)
};

class AliQAHistNavigator {

public:
    AliQAHistNavigator( Int_t run=0 );
    virtual ~AliQAHistNavigator();

    Bool_t GetHistogram(TH1*& hist);
    Bool_t Next();
    Bool_t Prev();

    Bool_t SetFile (TString file ); 
    Bool_t SetFile (Int_t file ); 
    Bool_t SetDetector( TString detector );
    Bool_t SetDetector( Int_t detector );
    Bool_t SetLevel( TString type );
    Bool_t SetLevel( Int_t type );
    Bool_t SetItem( TString histo );
    Bool_t SetItem( Int_t histo );

    void SetLoopAllFiles( const Bool_t s=kTRUE ) {fLoopAllFiles=s;}
    void SetLoopAllDetectors( const Bool_t s=kTRUE ) {fLoopAllDetectors=s;}
    void SetLoopAllLevels( const Bool_t s=kTRUE ) {fLoopAllLevels=s;}
    
    TString GetDetectorName();
    TString GetLevelName();
    TString GetItemName();
    TString GetFileName();
    TString GetDirName();
    TString GetPath(AliQADirListItem* item);
    
    AliQADirList* GetFileList() {return fPListOfFiles;}
    AliQADirList* GetDetectorList() {return fPCurrFile;}
    AliQADirList* GetLevelList() {return fPCurrDetector;}
    TList*    GetItemList(); 
    AliQADirList* GetCurrListOfFiles() {return fPListOfFiles;}
    AliQADirList* GetCurrFile() {return fPCurrFile;}
    AliQADirList* GetCurrDetector() {return fPCurrDetector;}
    AliQADirList* GetCurrLevel() {return fPCurrLevel;}
    AliQADirListItem* GetCurrItem() {return fPCurrItem;}
    
    Bool_t InitOK() {return fInitOK;}
    Bool_t ReReadFiles();
    void SetExpertMode(Bool_t mode);

    Bool_t CloneDirStructure();

private:

    Bool_t OpenCurrentFile();
    Bool_t OpenCurrentDirectory();
    Bool_t GetListOfFiles();
    Bool_t Crawl(AliQADirList* parent);
    
    TFile* fPFile;  //pointer to current open file
    TFile* fPCORRFile; //pointer to file with ntuple
    TFile* fPQAResultFile; //pointer to file with AliQA object
    Int_t fRun;     //runnumber

    //The state of the navigator, these help navigate the "tree"
    AliQADirList* fPCurrFile;  //current list holding detectors
    AliQADirList* fPCurrDetector; //current list holding levels
    AliQADirList* fPCurrLevel;  //current list holding histograms
    AliQADirListItem* fPCurrItem;  //current histogram name
    
    AliQADirList* fPListOfFiles; //Tree-like structure of lists within lists mirroring the layout of histogtams in files

    Bool_t fLoopAllFiles;  //whether to loop over all files
    Bool_t fLoopAllDetectors;  //whether to loop over all detectors
    Bool_t fLoopAllLevels;   //whether to loop over all levels
    
    Bool_t fInitOK;  //whether there is data to navigate
    Bool_t fExpertMode; //expert histogram mode
    TString fExpertDirName; //expert dir name
    TList* fPEmptyList;

    AliQAHistNavigator(const AliQAHistNavigator&);            // Not implemented
    AliQAHistNavigator& operator=(const AliQAHistNavigator&); // Not implemented

    ClassDef(AliQAHistNavigator,999)     //AliQAHistNavigator class
};

#endif

