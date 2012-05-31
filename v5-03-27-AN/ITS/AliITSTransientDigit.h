#ifndef ALIITSTRANSIENTDIGIT_H
#define ALIITSTRANSIENTDIGIT_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include <AliITSdigitSDD.h>

//______________________________________________________________________
class AliITSTransientDigit : public AliITSdigitSDD {

 public:
    AliITSTransientDigit() : fTrackList(0) {}
    // Standard constructor with digits and "phys"
    AliITSTransientDigit(Float_t phys,const Int_t *digits);
    virtual ~AliITSTransientDigit(){/*destructor delets TObjArray fTracklist */
	delete fTrackList;}
    //copy constructor
    AliITSTransientDigit(const AliITSTransientDigit &source);
    //assignment operator
    AliITSTransientDigit& operator=(const AliITSTransientDigit &source);
    // returns pointer to the TObjArray of tracks and associated charges
    TObjArray  *TrackList() const {return fTrackList;}
    //returns element i of fTrackList
    TObject *TrackItem(Int_t i) const {return fTrackList->At(i);}
    //put TObject into fTrackList at location i
    void PutTrackItem(TObject *obj,Int_t i){fTrackList->AddAt(obj,i);}
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function
    virtual Int_t Read(const char *name) {return AliITSdigitSDD::Read(name);}
    virtual void Print(Option_t *option="") const {AliITSdigitSDD::Print(option);}
 protected:
    TObjArray *fTrackList;  // track list 

    ClassDef(AliITSTransientDigit,1)  // Transient digit for set: ITS

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTransientDigit &source);
istream &operator>>(istream &os,AliITSTransientDigit &source);

#endif
