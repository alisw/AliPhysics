#ifndef ALIITSDIGITSSD_H
#define ALIITSDIGITSSD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
////////////////////////////////////////////////////////
// Digit class for SSD                                //
////////////////////////////////////////////////////////

#include <AliITSdigit.h>

//______________________________________________________________________
class AliITSdigitSSD: public AliITSdigit {

 public:
    AliITSdigitSSD(); //default constructor
    //Standard constructor with digits
    AliITSdigitSSD(const Int_t *digits);
    //Standard constructor with digits, tracks, and hits
    AliITSdigitSSD(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSdigitSSD(){/* destructor */}
    // returns the array size used to store Tracks and Hits
    // static Int_t GetNTracks() {return fgkSssd;}
    Int_t  GetSignal() const {/* returns signal*/return fSignal;}
    Int_t  GetStripNumber() const {/* returns strip number*/return fCoord2;}
    //returns 1  when side P and 0 when side N
    Int_t  IsSideP() const {if(fCoord1==0) return 1; else return 0; }
    // returns the pointer to the array of hits which made this digit
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t);
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}


 protected:
    
    
    ClassDef(AliITSdigitSSD,3)   // Simulated digit object for SSD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSSD &source);
istream &operator>>(istream &os,AliITSdigitSSD &source);

#endif
