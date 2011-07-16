#ifndef ALIESDINPUTHANDLERRP_H
#define ALIESDINPUTHANDLERRP_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliESDInputHandler.h 24521 2008-03-14 16:43:54Z morsch $ */

//-------------------------------------------------------------------------
//     ESD Input Handler realisation of the AliVEventHandler interface
//     Automatic loading of RecPoint Trees
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliESDInputHandler.h"
class TList;
class TTree;
class TDirectoryFile;
class TString;


class AliESDInputHandlerRP : public AliESDInputHandler {

 public:
    AliESDInputHandlerRP();
    AliESDInputHandlerRP(const char* name, const char* title);
    virtual ~AliESDInputHandlerRP();
    virtual Bool_t       Init(Option_t* opt);
    virtual Bool_t       Init(TTree* tree, Option_t* opt) {return AliESDInputHandler::Init(tree, opt);}
    virtual Bool_t       InitIO(Option_t* opt) {return Init(opt);};
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       FinishEvent();
    virtual Bool_t       LoadEvent(Int_t iev);
    virtual Bool_t       Notify() { return AliESDInputHandler::Notify();}
    virtual Bool_t       Notify(const char* path);
    virtual void         ResetIO();
    //
    virtual TTree*       GetTreeR(const char* det);
 private:
    Bool_t      OpenFile(Int_t i);
    AliESDInputHandlerRP(const AliESDInputHandlerRP& handler);             
    AliESDInputHandlerRP& operator=(const AliESDInputHandlerRP& handler);  
 private:
    TObjArray*      fRTrees;           // List of RecPoint Trees
    TObjArray*      fRDirs;            // List of RecPoint directories
    TList*          fRFiles;           // List of RecPoint Files
    TList*          fDetectors;        // List of detector names
    TDirectoryFile *fDirR;             //! Directory for RP Tree
    Int_t           fEventNumber;      //! Current event number
    Int_t           fFileNumber;       //! Input file number
    Int_t           fEventsPerFile;    //! Number of events per file
    const Char_t   *fExtension;        //! File name extension
    TString        *fPathName;         //! Input file path
    Bool_t          fIsArchive;        //! True if directory is an archive
    ClassDef(AliESDInputHandlerRP, 1);
};

#endif
