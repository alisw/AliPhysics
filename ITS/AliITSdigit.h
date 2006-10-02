#ifndef ALIITSDIGIT_H
#define ALIITSDIGIT_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

////////////////////////////////////////////////
//  Digits classes for all ITS detectors      //
////////////////////////////////////////////////
#include <Riostream.h>
#include <Riostream.h>
#include <TObject.h>

class TObjArray;
class TArrayI;
class TArrayF;


//______________________________________________________________________
class AliITSdigit: public TObject  {

 public:
    AliITSdigit();
    //Standard Constructor. Fills class from array digits
    AliITSdigit(const Int_t *digits);
    //Destructor
    virtual ~AliITSdigit() { }
    // returns the array size used to store Tracks and Hits
    // virtual Int_t GetNTracks() {return 0;}
    //returns pointer to array of tracks numbers
    virtual Int_t *GetTracks() {return 0;}
    // returns pointer to array of hits numbers for this module (as given by
    // AliITSmodule).
    virtual Int_t *GetHits() {return 0;}
    // returns track number kept in the array element i of fTracks 
    virtual Int_t GetTrack(Int_t) const {return 0;}
    // returns hit number kept in the array element i of fHits 
    virtual Int_t GetHit(Int_t) const {return 0;}
    virtual Int_t GetCoord1() const {return fCoord1;} // returns fCoord1
    virtual Int_t GetCoord2() const {return fCoord2;} // returns fCoord2
    virtual Int_t GetSignal() const {return fSignal;} // returns fSignal
    virtual Int_t GetCompressedSignal() const {return GetSignal();} // overloaded in AliITSdigitSDD
    virtual void SetCoord1(Int_t i){fCoord1 = i;} // Sets fCoord1 value
    virtual void SetCoord2(Int_t i){fCoord2 = i;} // Sets fCoord12value
    virtual void SetSignal(Int_t i){fSignal = i;} // Sets fSignal value

    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual void Print(ostream *os); // Class ascii print function
    virtual Int_t Read(const char *name) {return TObject::Read(name);}
    virtual void Read(istream *os);  // Class ascii read function

 protected:
    Int_t fCoord1; // Cell number on Z axis (SPD+SDD), flag for side type (SSD)
    Int_t fCoord2; // Cell number on X axis (SPD+SDD), strip number (SSD)
    Int_t fSignal; // Signal in ADC counts

    ClassDef(AliITSdigit,1)     // Real data digit object for set:ITS

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigit &source);
istream &operator>>(istream &os,AliITSdigit &source);

#endif
