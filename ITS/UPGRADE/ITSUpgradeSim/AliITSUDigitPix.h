#ifndef ALIITSDIGITSPD_H
#define ALIITSDIGITSPD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////
// Digit class for SPD                                     //
/////////////////////////////////////////////////////////////
#include <AliITSdigit.h>

//______________________________________________________________________
class AliITSUDigitPix: public AliITSdigit {

 public:
    AliITSUDigitPix(); //default creator
    AliITSUDigitPix(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSUDigitPix(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSUDigitPix(){/*destructor*/}
    // returns the signal in electrons
    Int_t GetSignalPix() const {return fSignalPix;}
    Int_t GetROCycle()   const {return fROCycle;}
    virtual Int_t GetListOfTracks(TArrayI &t);
    // set signal in electrons
    void SetSignalPix(Int_t sig) {fSignalPix = sig;}
    void SetROCycle(Int_t cycle) {fROCycle = cycle;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual void Read(istream *os);  // Class ascii read function
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    
    Int_t fSignalPix;   // Signal in electrons
    Int_t fROCycle;     // readout cycle
    ClassDef(AliITSUDigitPix,1)   // Simulated digit object for Pixels

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSUDigitPix &source);
istream &operator>>(istream &os,AliITSUDigitPix &source);

#endif
