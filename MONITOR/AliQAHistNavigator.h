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

class AliQAHistNavigator {

public:
    AliQAHistNavigator( Int_t run=0, Int_t rev=0 );

    Bool_t GetHistogram(TH1*& hist);
    Bool_t GetNextHistogram(TH1*& hist);
    Bool_t GetPrevHistogram(TH1*& hist);
    Bool_t Next();
    Bool_t Prev();
    void PrintDebugInfo();
    Bool_t DumpList(TString file="AliQAHistNavigator.conf");
    Bool_t ReadList(TString file="AliQAHistNavigator.conf");

    Bool_t SetFile (TString file ); 
    Bool_t SetFile (Int_t file ); 
    Bool_t SetDetector( TString detector );
    Bool_t SetDetector( Int_t detector );
    Bool_t SetLevel( TString type );
    Bool_t SetLevel( Int_t type );
    Bool_t SetHist( TString histo );
    Bool_t SetHist( Int_t histo );
    void SetLoopAllFiles( const Bool_t s=kTRUE ) {fLoopAllFiles=s;}
    void SetLoopAllDetectors( const Bool_t s=kTRUE ) {fLoopAllDetectors=s;}
    void SetLoopAllLevels( const Bool_t s=kTRUE ) {fLoopAllLevels=s;}
    TString GetDetectorName();
    TString GetLevelName();
    TString GetHistName();
    TString GetFileName();
    TString GetDirName();
    TList* GetFileList() {return fPListOfFiles;}
    TList* GetDetectorList() {return fPCurrFile;}
    TList* GetLevelList() {return fPCurrDetector;}
    TList* GetHistList() {return fPCurrLevel;}
    TList* GetCurrListOfFiles() {return fPListOfFiles;}
    TList* GetCurrFile() {return fPCurrFile;}
    TList* GetCurrDetector() {return fPCurrDetector;}
    TList* GetCurrLevel() {return fPCurrLevel;}
    TObjString* GetCurrHistName() {return fPCurrHistName;}

    Bool_t CloneDirStructure();

private:

    Bool_t OpenCurrentFile();
    Bool_t OpenCurrentDirectory();
    Bool_t GetListOfFiles();
    Bool_t Crawl(TList* parent);
    
    TFile* fPFile;  //pointer to current open file
    Int_t fRun;     //runnumber
    Int_t fCyc;     //Cycle number

    //The state of the navigator, these help navigate the "tree"
    TList* fPCurrFile;  //current list holding detectors
    TList* fPCurrDetector; //current list holding levels
    TList* fPCurrLevel;  //current list holding histograms
    TObjString* fPCurrHistName;  //current histogram name
    
    TList* fPListOfFiles; //Tree-like structure of lists within lists mirroring the layout of histogtams in files

    Bool_t fLoopAllFiles;  //whether to loop over all files
    Bool_t fLoopAllDetectors;  //whether to loop over all detectors
    Bool_t fLoopAllLevels;   //whether to loop over all levels

    ClassDef(AliQAHistNavigator,999)     //AliQAHistNavigator class
};

#endif

