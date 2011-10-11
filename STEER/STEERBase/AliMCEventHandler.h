  // -*- mode: C++ -*-
#ifndef ALIMCEVENTHANDLER_H
#define ALIMCEVENTHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliMCEvent
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is contained in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//-------------------------------------------------------------------------
#include "AliInputEventHandler.h"
#include "AliHeader.h"
#include <TExMap.h>

class TFile;
class TTree;
class TList;

class TParticle;
class TString;
class TClonesArray;
class TDirectoryFile;

class AliMCEvent;



class AliMCEventHandler : public AliInputEventHandler
{
public:

    enum PreReadMode_t {kNoPreRead = 0, kLmPreRead = 1, kHmPreRead = 2};

    AliMCEventHandler();
    AliMCEventHandler(const char* name, const char* title);
    virtual ~AliMCEventHandler();
    virtual void         SetOutputFileName(const char* /* fname */) {;}
    virtual const char*  GetOutputFileName() {return 0;}
    virtual void         SetInputPath(const char* fname); 
    virtual void         SetInputTree(TTree* /*tree*/) {;}
    virtual TString*     GetInputPath() const {return fPathName;}
    virtual Bool_t       Init(Option_t* opt);
    virtual Bool_t       GetEntry() {return kTRUE;}
    virtual Bool_t       InitIO(Option_t* opt) {return Init(opt);};
    virtual Bool_t       Init(TTree* /*tree*/, Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    virtual void         ResetIO();
    virtual Bool_t       LoadEvent(Int_t iev);
    virtual void         SetReadTR(Bool_t flag) { fReadTR = flag; }
    virtual void         AddSubsidiaryHandler(AliMCEventHandler* handler);
    virtual void         SetNumberOfEventsInContainer(Int_t nev) {fEventsInContainer = nev;}
    virtual void         SetPreReadMode(PreReadMode_t mode) {fPreReadMode = mode;}
    //
    AliMCEvent* MCEvent() const {return fMCEvent;}
    TTree*      TreeTR()  const {return fTreeTR;}
    TTree*      TreeK()   const {return fTreeK;}
    virtual TTree*      GetTree() const {return fTreeE;}
    Int_t       GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs);
    void        DrawCheck(Int_t i, Int_t search=0);
    Bool_t      InitOk() const {return fInitOk;}
    // Label manipulation
    void   SelectParticle(Int_t i);
    Bool_t IsParticleSelected(Int_t i);
    void   CreateLabelMap();
    Int_t  GetNewLabel(Int_t i);

private:
    Bool_t      OpenFile(Int_t i);
    void  VerifySelectedParticles();
    AliMCEventHandler(const AliMCEventHandler& handler);             
    AliMCEventHandler& operator=(const AliMCEventHandler& handler);  
private:
    AliMCEvent            *fMCEvent;            //! MC Event
    TFile                 *fFileE;              //! File with TreeE
    TFile                 *fFileK;              //! File with TreeK
    TFile                 *fFileTR;             //! File with TreeTR
    TTree                 *fTreeE;              //! TreeE  (Event Headers)
    TTree                 *fTreeK;              //! TreeK  (kinematics tree)
    TTree                 *fTreeTR;             //! TreeTR (track references tree)
    TDirectoryFile        *fDirK;               //! Directory for Kine Tree
    TDirectoryFile        *fDirTR;              //! Directory for TR Tree
    TExMap                 fParticleSelected;   //! List of selected MC particles for t
    TExMap                 fLabelMap;           //! Stores the Map of MC (ESDLabel,AODlabel)  
    Int_t                  fNEvent;             //! Number of events
    Int_t                  fEvent;              //! Current event
    TString               *fPathName;           //! Input file path 
    const Char_t          *fExtension;          //! File name extension 
    Int_t                  fFileNumber;         //! Input file number
    Int_t                  fEventsPerFile;      //! Number of events per file
    Bool_t                 fReadTR;             // determines if TR shall be read
    Bool_t                 fInitOk;             // Initialization ok
    TList                 *fSubsidiaryHandlers; //! List of subsidiary MC handlers (for example for Background)
    Int_t                  fEventsInContainer;  //! Number of events in container class
    PreReadMode_t          fPreReadMode;        //! Pre reading mode
    
    ClassDef(AliMCEventHandler,1)  //MC Truth EventHandler class
};
#endif 

