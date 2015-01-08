#ifndef ALIITSDIGITSPD_H
#define ALIITSDIGITSPD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/////////////////////////////////////////////////////////////
// Digit class for SPD                                     //
/////////////////////////////////////////////////////////////
#include <AliITSdigit.h>

//______________________________________________________________________
class AliITSdigitSPD: public AliITSdigit {

 public:
    AliITSdigitSPD(); //default creator
    AliITSdigitSPD(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSdigitSPD(){/*destructor*/}
    // returns the signal in electrons
    Int_t GetSignalSPD() const {return fSignalSPD;}
    virtual Int_t GetListOfTracks(TArrayI &t);
    // set signal in electrons
    void SetSignalSPD(Int_t sig) {fSignalSPD = sig;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual void Read(istream *os);  // Class ascii read function
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    
    Int_t fSignalSPD;   // Signal in electrons

    ClassDef(AliITSdigitSPD,3)   // Simulated digit object for SPD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSPD &source);
istream &operator>>(istream &os,AliITSdigitSPD &source);

#endif
