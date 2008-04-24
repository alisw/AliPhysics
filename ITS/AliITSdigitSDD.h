#ifndef ALIITSDIGITSDD_H
#define ALIITSDIGITSDD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
////////////////////////////////////////////////////////
// Digit class for SDD                                //
////////////////////////////////////////////////////////
#include <AliITSdigit.h>

class AliITSCalibrationSDD;

//______________________________________________________________________
class AliITSdigitSDD: public AliITSdigit {

 public:
    AliITSdigitSDD(); //default creator
    //standard c.tor  with digits and "phys"
    AliITSdigitSDD(Float_t phys,const Int_t *digits);
    //standard c.tor with digits, tracls, hits, "phys", and charge
    AliITSdigitSDD( Float_t phys,const Int_t *digits,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges);
    //constructor setting also fSignalExpanded
    AliITSdigitSDD( Float_t phys,const Int_t *digits,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges, Int_t sige);
    virtual ~AliITSdigitSDD(){/* destructor*/}
    // returns the array size used to store Tracks and Hits
    virtual Int_t GetSignal() const {return fSignalExpanded;}
    virtual Int_t GetCompressedSignal() const {return fSignal;}
     // Return charge deposited by this track/hit
    virtual Float_t GetCharge(Int_t i) const {return fTcharges[i];}
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t,TArrayF &c);
    void SetSignalExpanded(Int_t sig){fSignalExpanded = sig;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Read(istream *os);  // Class ascii read function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    void InitObject(Float_t phys,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges);
    
                            // 3 hits temporarily - it will be only 1
    Float_t fTcharges[fgkSize];   //[fgkSize] charge per track making this digit 
    Float_t fPhysics;       // signal particles contribution to signal
    Int_t fSignalExpanded;  // 10 bit signal
    ClassDef(AliITSdigitSDD,4)   // Simulated digit object for SDD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSDD &source);
istream &operator>>(istream &os,AliITSdigitSDD &source);

#endif
