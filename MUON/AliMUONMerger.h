#ifndef ALIMUONMERGER_H
#define ALIMUONMERGER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// #include "AliMerger.h"
// #include "AliMergable.h"

class AliMUONPadHit;
class AliHitMap;

typedef enum {kDigitize=0, kMerge = 1} MergeMode_t;

class AliMUONMerger {
 public:
    
    AliMUONMerger();
    virtual ~AliMUONMerger();
    
    // Compare pad hits
    virtual Bool_t Exists(const AliMUONPadHit * sdigit);
    // Update a pad hit
    virtual  void Update(AliMUONPadHit *sdigit);
    // Create a new hit
    virtual  void CreateNew(AliMUONPadHit *sdigit);

    // Initialize merging and digitisation
    virtual void Init();

    // Do the main work
    void Digitise() ;
    
    // Setters -> Later Communication with gAlice 
    void SetSignalEventNumber(Int_t i)     {fEvNrSig = i;}
    void SetBackgroundEventNumber(Int_t i) {fEvNrBgr = i;}    
    void SetBackgroundFileName(char* file) {fFnBgr = file;}        
    void SetMode(MergeMode_t mode) {fMerge = mode;}
	
    enum {kBgTag = -1};
    
 private:    
    // Open the bgr file
    TFile *InitBgr();
    void SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr);
    
 private:
    TTree *fTrH1;                   // ! Hits Tree for background event
    TClonesArray *fHitsBgr;         // ! List of hits for one track only
    TClonesArray *fPadHitsBgr;      // ! List of clusters for one track only
    AliHitMap **fHitMap;            // ! pointer to array of pointers to hitmaps
    Int_t fNch;                     // ! chamber nr (loop variable)
    Int_t fTrack;                   // ! track nr (loop variable)
    TObjArray *fList;               // ! list of AliMUONTransientDigit
    TObjArray *fTrList;             // ! list of tracks
    TClonesArray *fAddress;         // ! pointer to TClonesArray of TVectors with trackinfo
    Int_t fCounter;                 // ! nr. of AliMUONTransientDigit
    Int_t fCountadr;                // ! counter for trinfo
    Int_t fDigits[6];               // ! array with digits
    Int_t fEvNrSig;                 // signal     event number
    Int_t fEvNrBgr;                 // background event number    
    MergeMode_t fMerge;             // merging type kDigitize, kMerge
    char  *fFnBgr;                  // background file name
    TFile *fBgrFile;                // Pointer to background file
    
    ClassDef(AliMUONMerger,0)
};    
#endif

