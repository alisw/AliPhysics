#ifndef ALIITSDIGITSPD_H
#define ALIITSDIGITSPD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////
// Digit class for SPD                                     //
/////////////////////////////////////////////////////////////
#include <AliITSdigit.h>

//______________________________________________________________________
class AliITSdigitPixUpg: public AliITSdigit {

 public:
    AliITSdigitPixUpg(); //default creator
    AliITSdigitPixUpg(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSdigitPixUpg(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSdigitPixUpg(){/*destructor*/}
    // returns the signal in electrons
    Int_t GetSignalPix() const {return fSignalPix;}
    virtual Int_t GetListOfTracks(TArrayI &t);
    // set signal in electrons
    void SetSignalPix(Int_t sig) {fSignalPix = sig;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual void Read(istream *os);  // Class ascii read function
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    
    Int_t fSignalPix;   // Signal in electrons

    ClassDef(AliITSdigitPixUpg,1)   // Simulated digit object for Pixels

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitPixUpg &source);
istream &operator>>(istream &os,AliITSdigitPixUpg &source);

#endif
