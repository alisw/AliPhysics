// -*- mode: C++ -*- 
#ifndef ALIMCEVENTHANDLER_H
#define ALIMCEVENTHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliMCEvent
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is containe in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//-------------------------------------------------------------------------

#include "AliVirtualEventHandler.h"
class TFile;
class TTree;
class TParticle;
class TClonesArray;
class AliHeader;
class AliStack;


class AliMCEventHandler : public AliVirtualEventHandler 
{
public:
    AliMCEventHandler();
    AliMCEventHandler(const char* name, const char* title);
    virtual ~AliMCEventHandler();
    virtual void         SetOutputFileName(char* /* fname */) {;}
    virtual char*        GetOutputFileName() {return 0;}
    virtual void         SetInputPath(char* fname) {fPathName = fname;}
    virtual char*        GetInputPath() {return fPathName;}
    virtual Bool_t       InitIO(Option_t* opt);
    virtual Bool_t       BeginEvent();
    virtual Bool_t       Notify();
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    virtual void         ResetIO();
    //
    AliStack* Stack()  {return fStack;}
    TTree*    TreeTR() {return fTreeTR;}
    Int_t     GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs);
    void      DrawCheck(Int_t i);
private:
    void      ReorderAndExpandTreeTR();
    
private:
    TFile            *fFileE;            //! File with TreeE
    TFile            *fFileK;            //! File with TreeK
    TFile            *fFileTR;           //! File with TreeTR
    TFile            *fTmpFileTR;        //! Temporary file with TreeTR to read old format
    TTree            *fTreeE;            //! TreeE  (Event Headers)
    TTree            *fTreeK;            //! TreeK  (kinematics tree)
    TTree            *fTreeTR;           //! TreeTR (track references tree)
    TTree            *fTmpTreeTR;        //! Temporary tree TR to read old format
    AliStack         *fStack;            //! Current pointer to stack
    AliHeader        *fHeader;           //! Current pointer to header
    TClonesArray     *fTrackReferences;  //! Current list of track references
    Int_t             fNEvent;           //! Number of events
    Int_t             fEvent;            //! Current event
    Int_t             fNprimaries;       //! Number of primaries
    Int_t             fNparticles;       //! Number of particles
    char             *fPathName;         //! Input file path 
    ClassDef(AliMCEventHandler,1)  //MC Truth EventHandler class 
};
#endif 

