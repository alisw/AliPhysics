#ifndef ALIMUONDIGITIZER_H
#define ALIMUONDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "AliDigitizer.h"

class AliRunDigitizer;
class AliMUONPadHit;
class AliMUONHitMapA1;

class AliMUONDigitizer : public AliDigitizer {
 public:
    
    AliMUONDigitizer();
    AliMUONDigitizer(AliRunDigitizer * manager);
    virtual ~AliMUONDigitizer();

    // Compare pad hits
    virtual Bool_t Exists(const AliMUONPadHit * sdigit) const;
    // Update a pad hit
    virtual  void Update(AliMUONPadHit *sdigit);
    // Create a new hit
    virtual  void CreateNew(AliMUONPadHit *sdigit);

    // Initialize merging and digitization
    virtual Bool_t Init();

    // Do the main work
    virtual void Exec(Option_t* option=0);
    
    Int_t GetDebug() const {return fDebug;}       // get debug level
    void SetDebug(Int_t level){fDebug = level;}   // set debug level    
    enum {kBgTag = -1};
    
 private:    
    void SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr);
    
 private:
    TClonesArray* fHits;
    TClonesArray* fPadHits;
    AliMUONHitMapA1 **fHitMap;      //! pointer to array of pointers to hitmaps
    Int_t fNch;                     //! chamber nr (loop variable)
    Int_t fTrack;                   //! track nr (loop variable)
    TObjArray *fTDList;             //! list of AliMUONTransientDigits
    Int_t fCounter;                 //! nr. of AliMUONTransientDigit
    Bool_t fSignal;                 //! kTRUE if signal file is processed
    Int_t fMask;                    //! mask dependent on input file
    Int_t fDigits[6];               //! array with digits
    Int_t fDebug;                   //! debug level

    ClassDef(AliMUONDigitizer,1)  // MUON merging/digitization
};    
#endif

