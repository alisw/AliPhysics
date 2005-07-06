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
    static Int_t GetNTracks() {return fgkSssd;}
    Int_t  GetSignal() const {/* returns signal*/return fSignal;}
    Int_t  GetStripNumber() const {/* returns strip number*/return fCoord2;}
    //returns 1  when side P and 0 when side N
    Int_t  IsSideP() const {return fCoord1;}
    // returns pointer to the array of tracks which make this digit
    virtual Int_t *GetTracks() {return &fTracks[0];}
    // returns the pointer to the array of hits which made this digit
    virtual Int_t *GetHits() {return &fHits[0];}
    // returns track number kept in the array element i of fTracks 
    virtual Int_t GetTrack(Int_t i) const {return fTracks[i];}
    // returns hit number kept in the array element i of fHits 
    virtual Int_t GetHit(Int_t i) const {return fHits[i];}
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t);
    //copy the array trks[fgkSssd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fgkSssd;i++) fTracks[i]=trks[i];}
    //copy the array hits[fgkSssd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fgkSssd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}


 protected:
    static const Int_t fgkSssd = 10; // size of fTracks and fHits arrays
    
    // debugging  -- goes to the dictionary
    Int_t fTracks[fgkSssd]; //[fgkSssd] tracks making this digit 
    Int_t fHits[fgkSssd];   //[fgkSssd] hits associated to the tracks
                        // 3 hits temporarily - it will be only 1
    
    ClassDef(AliITSdigitSSD,2)   // Simulated digit object for SSD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSSD &source);
istream &operator>>(istream &os,AliITSdigitSSD &source);

#endif
