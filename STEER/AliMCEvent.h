// -*- mode: C++ -*- 
#ifndef ALIMCEVENT_H
#define ALIMCEVENT_H
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


class AliMCEvent : public AliVirtualEventHandler 
{
public:
    AliMCEvent();
    AliMCEvent(const char* name, const char* title);
    virtual ~AliMCEvent();
    virtual void         SetOutputFileName(char* /* fname */) {;}
    virtual char*        GetOutputFileName() {return 0;}
    virtual Bool_t       InitIO(Option_t* opt);
    virtual Bool_t       BeginEvent();
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    //
    AliStack* Stack()  {return fStack;}
    TTree*    TreeTR() {return fTreeTR;}
    Int_t     GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs);
    void      DrawCheck(Int_t i);
	    
private:
    TFile            *fFileE;            //! File with TreeE
    TFile            *fFileK;            //! File with TreeK
    TFile            *fFileTR;           //! File with TreeTR
    TTree            *fTreeE;            //! TreeE  (Event Headers)
    TTree            *fTreeK;            //! TreeK  (kinematics tree)
    TTree            *fTreeTR;           //! TreeTR (track references tree)
    AliStack         *fStack;            //! Current pointer to stack
    AliHeader        *fHeader;           //! Current pointer to header
    TClonesArray     *fTrackReferences;  //! Current list of tarck references
    Int_t             fNEvent;           //! Number of events
    Int_t             fEvent;            //! Current event
    Int_t             fNprimaries;       //! Number of primaries
    Int_t             fNparticles;       //! Number of particles 
    ClassDef(AliMCEvent,1)  //MCEvent class 
};
#endif 

