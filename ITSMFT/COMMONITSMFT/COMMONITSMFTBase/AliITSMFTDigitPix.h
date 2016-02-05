#ifndef ALIITSDIGITSPD_H
#define ALIITSDIGITSPD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////
// Digit class for SPD                                     //
/////////////////////////////////////////////////////////////
#include <TObject.h>

class TArrayI;

//______________________________________________________________________
class AliITSMFTDigitPix: public TObject {

 public:
    AliITSMFTDigitPix(); //default creator
    AliITSMFTDigitPix(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSMFTDigitPix(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    //
    AliITSMFTDigitPix(const AliITSMFTDigitPix &h);
    AliITSMFTDigitPix& operator=(const AliITSMFTDigitPix &h);
    //
    virtual ~AliITSMFTDigitPix();

    // returns the array size used to store Tracks and Hits
    static Int_t GetNTracks() {return fgkSize;}
    //returns pointer to array of tracks numbers
    Int_t *GetTracks()  {return &fTracks[0];}
    // returns track number kept in the array element i of fTracks 
    Int_t GetTrack(Int_t i) const {return fTracks[i];}

    Int_t GetHit(Int_t i) const {return fHits[i];}
    Int_t GetCoord1() const {return fCoord1;} // returns fCoord1
    Int_t GetCoord2() const {return fCoord2;} // returns fCoord2
    Int_t GetSignal() const {return fSignal;} // returns fSignal
    // returns the signal in electrons
    Int_t GetSignalPix() const {return fSignalPix;}
    Int_t GetROCycle()   const {return fROCycle;}
    Int_t GetListOfTracks(TArrayI &t);

    void SetCoord1(Int_t i){fCoord1 = i;} // Sets fCoord1 value
    void SetCoord2(Int_t i){fCoord2 = i;} // Sets fCoord12value
    void SetSignal(Int_t i){fSignal = i;} // Sets fSignal value
    void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}

    // set signal in electrons
    void SetSignalPix(Int_t sig) {fSignalPix = sig;}
    void SetROCycle(Int_t cycle) {fROCycle = cycle;}
    void Print(std::ostream *os) const; // Class ascii print function
    void Print(Option_t *option="") const {TObject::Print(option);}
    void Read(std::istream *os);  // Class ascii read function
    Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    static const Int_t fgkSize = 10;//array size
    Int_t   fTracks[fgkSize];   //[fgkSize] tracks making this digit 
    Int_t   fHits[fgkSize];     //[fgkSize] hits associated to the tracks

    Int_t fCoord1; // Cell number on Z axis (SPD+SDD), flag for side type (SSD)
    Int_t fCoord2; // Cell number on X axis (SPD+SDD), strip number (SSD)
    Int_t fSignal; // Signal in ADC counts

    Int_t fSignalPix;   // Signal in electrons
    Int_t fROCycle;     // readout cycle
    ClassDef(AliITSMFTDigitPix,1)   // Simulated digit object for Pixels

};
#endif
