// -*- mode: C++ -*- 
#ifndef ALIMCEVENT_H
#define ALIMCEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliMCEvent
//      
// Origin: Andreas.Morsch, CERN, andreas.morsch@cern.ch 
//-------------------------------------------------------------------------


#include <TTree.h>
#include <TRefArray.h>

#include <AliVEvent.h>
#include "AliVHeader.h"
#include "AliMCParticle.h"

class AliStack;
class AliHeader;
class AliGenEventHeader;

class TClonesArray;

class AliMCEvent : public AliVEvent {

public:

    AliMCEvent();
    virtual ~AliMCEvent() {;} 
    AliMCEvent(const AliMCEvent& mcEvnt); 
    AliMCEvent& operator=(const AliMCEvent& mcEvnt);
    //
    // Methods implementing the interface
    //
    // Services
    virtual void AddObject(TObject* /*obj*/)               {;}
    virtual TObject* FindListObject(const char */*name*/)  {return 0;}
    virtual TList* GetList() const                         {return 0;}
    virtual void CreateStdContent()                        {;} 
    virtual void GetStdContent()                           {;}
    virtual void ReadFromTree(TTree * /*tree*/, Option_t* /*opt*/) {;}
    virtual void WriteToTree(TTree* /*tree*/)  const {;}

    virtual void SetStdNames()                             {;}
    virtual void Print(Option_t */*option=""*/)  const     {;}
    

    // Header
    virtual AliVHeader* GetHeader()          const         {return 0;}

    // Delegated methods for fESDRun or AODHeader
  
    virtual void     SetRunNumber(Int_t /*n*/)             {;}
    virtual void     SetPeriodNumber(UInt_t /*n*/)         {;}
    virtual void     SetMagneticField(Double_t /*mf*/)     {;}
    
  
    virtual Int_t    GetRunNumber()          const         {return 0;}
    virtual UInt_t   GetPeriodNumber()       const {return 0;}
    virtual Double_t GetMagneticField()      const         {return 0.;}

    // Setters not needed
    virtual void      SetOrbitNumber(UInt_t /*n*/)         {;}
    virtual void      SetBunchCrossNumber(UShort_t /*n*/)  {;}
    virtual void      SetEventType(UInt_t /*eventType*/)   {;}
    virtual void      SetTriggerMask(ULong64_t /*n*/)      {;}
    virtual void      SetTriggerCluster(UChar_t /*n*/)     {;} 

    virtual UInt_t    GetOrbitNumber()       const {return 0;}
    virtual UShort_t  GetBunchCrossNumber()  const {return 0;}
    
    virtual UInt_t    GetEventType()         const {return 0;}

    virtual ULong64_t GetTriggerMask()        const {return 0;}
    virtual UChar_t   GetTriggerCluster()     const {return 0;}
    virtual Double_t  GetZDCN1Energy()        const {return 0.;}
    virtual Double_t  GetZDCP1Energy()        const {return 0.;}
    virtual Double_t  GetZDCN2Energy()        const {return 0.;}
    virtual Double_t  GetZDCP2Energy()        const {return 0.;}
    virtual Double_t  GetZDCEMEnergy(Int_t /*i*/) 
                                              const {return 0.;}
    // Tracks
    virtual AliMCParticle *GetTrack(Int_t i) const;
    virtual Int_t     GetNumberOfTracks()    const {return fNparticles;}
    virtual Int_t     GetNumberOfV0s()       const {return -1;}
    virtual Int_t     GetNumberOfCascades()  const {return -1;}

    //
    // MC Specific methods
    //
    // Getters
    AliStack*    Stack()   {return fStack;}
    AliHeader*   Header()  {return fHeader;}
    AliGenEventHeader* GenEventHeader();
    // Services
    virtual void      ConnectTreeE (TTree* tree);
    virtual void      ConnectTreeK (TTree* tree);
    virtual void      ConnectTreeTR(TTree* tree);
    virtual void      Clean();
    virtual void      FinishEvent();
    virtual Int_t     GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs);
    virtual void      DrawCheck(Int_t i, Int_t search);
private:
    virtual void      ReorderAndExpandTreeTR();
    
private:
    AliStack         *fStack;           // Current pointer to stack
    TClonesArray     *fMCParticles;     // Pointer to list of particles
    TRefArray        *fMCParticleMap;   // Map of MC Particles
    AliHeader        *fHeader;          // Current pointer to header
    TClonesArray     *fTRBuffer;        // Track reference buffer    
    TClonesArray     *fTrackReferences; // Array of track references
    TTree            *fTreeTR;          // Pointer to Track Reference Tree
    TTree            *fTmpTreeTR;       // Temporary tree TR to read old format
    TFile            *fTmpFileTR;       // Temporary file with TreeTR to read old format
    Int_t             fNprimaries;      // Number of primaries
    Int_t             fNparticles;      // Number of particles
    ClassDef(AliMCEvent, 1)  // AliVEvent realisation for MC data
};


#endif 

