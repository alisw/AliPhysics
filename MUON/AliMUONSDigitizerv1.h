#ifndef ALIMUONSDIGITIZERV1_H
#define ALIMUONSDIGITIZERV1_H

// The AliMUONSDigitizer produces
// SDigits from Hits 
// J.P Cussonneau Subatech Feb 2004

class AliMUONPadHit;
class AliMUONHitMapA1;

class AliMUONSDigitizerv1 {

 public:    
    AliMUONSDigitizerv1();
    virtual ~AliMUONSDigitizerv1();

    // Create a new TransientDigit
    virtual void   AddTransientDigit(AliMUONTransientDigit * mTD);
    // Do the main work
    virtual void   Exec(Option_t* option=0);
    // Verifying a TransientDigit
    virtual Bool_t ExistTransientDigit(AliMUONTransientDigit * mTD); 
    // Getting debug level 
    Int_t          GetDebug() const {return fDebug;}       // get debug level
    // Initialize merging and digitization
    virtual Bool_t Init();
    // Generation of a TransientDigit : Response function of the chamber
    virtual void   MakeTransientDigit(Int_t itrack, Int_t ihit, AliMUONHit * mHit);
    // Setting debug level
    void           SetDebug(Int_t level){fDebug = level;}   // set debug level    
    enum {kBgTag = -1};
    // Updating a TransientDigit
    virtual void   UpdateTransientDigit(Int_t itrack, AliMUONTransientDigit * mTD);
    
    private:    
    void           SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr);
    
    private:
    AliMUONHitMapA1 **fHitMap;      //! pointer to array of pointers to hitmaps
    TObjArray *fTDList;             //! list of AliMUONTransientDigits
    Int_t fTDCounter;                 //! nr. of AliMUONTransientDigits
    Int_t fDebug;                   //! debug level
    Int_t fMask;                    //! mask dependent on input file
    Bool_t fSignal;                 //! kTRUE if signal file is processed


    ClassDef(AliMUONSDigitizerv1,0) 
};    
#endif

