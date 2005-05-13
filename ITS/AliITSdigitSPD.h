#ifndef ALIITSDIGITSPD_H
#define ALIITSDIGITSPD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include <AliITSdigit.h>

//______________________________________________________________________
class AliITSdigitSPD: public AliITSdigit {

 public:
    AliITSdigitSPD(); //default creator
    AliITSdigitSPD(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSdigitSPD(){/*destructor*/}
    // returns the array size used to store Tracks and Hits
    static Int_t GetNTracks() {return fgkSspd;}
    // returns the signal in electrons
    Int_t GetSignalSPD() const {return fSignalSPD;}
    // returns pointer to the array of tracks which make this digit
    virtual Int_t *GetTracks() {return &fTracks[0];}
     //returns the pointer to the array of hits which made this digit
    virtual Int_t *GetHits() {return &fHits[0];}
    // returns track number kept in the array element i of fTracks 
    virtual Int_t GetTrack(Int_t i) const {return fTracks[i];}
    // returns hit number kept in the array element i of fHits 
    virtual Int_t GetHit(Int_t i) const {return fHits[i];}
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t);
    //copy the array trks[fgkSspd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fgkSspd;i++) fTracks[i]=trks[i];}
    // set signal in electrons
    void SetSignalSPD(Int_t sig) {fSignalSPD = sig;}
    //copy the array hits[fgkSspd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fgkSspd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual void Read(istream *os);  // Class ascii read function
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    static const Int_t fgkSspd = 10; // size of fTracks and fHits arrays
    
    // debugging  -- goes to the dictionary
    Int_t fTracks[fgkSspd]; //[fgkSspd] tracks making this digit 
    Int_t fHits[fgkSspd];   //[fgkSspd] hits associated to the tracks
    Int_t fSignalSPD;   // Signal in electrons

    ClassDef(AliITSdigitSPD,2)   // Simulated digit object for SPD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSPD &source);
istream &operator>>(istream &os,AliITSdigitSPD &source);

#endif
