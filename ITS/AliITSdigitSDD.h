#ifndef ALIITSDIGITSDD_H
#define ALIITSDIGITSDD_H
/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
////////////////////////////////////////////////////////
// Digit class for SDD                                //
////////////////////////////////////////////////////////
#include <AliITSdigit.h>

class AliITSresponseSDD;

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
    //constructor setting fSignalExpanded through AliITSresponseSDD
    AliITSdigitSDD( Float_t phys,const Int_t *digits,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges, 
                    AliITSresponseSDD *resp);
    virtual ~AliITSdigitSDD(){/* destructor*/}
    // returns the array size used to store Tracks and Hits
    static Int_t GetNTracks() {return fgkSsdd;}
    // returns pointer to the array of tracks which make this digit
    virtual Int_t *GetTracks()  {return &fTracks[0];}
    // returns the pointer to the array of hits which made this digit
    virtual Int_t *GetHits()  {return &fHits[0];}
    virtual Int_t GetSignal() const {return fSignalExpanded;}
    virtual Int_t GetCompressedSignal() const {return fSignal;}
    // returns track number kept in the array element i of fTracks 
    virtual Int_t GetTrack(Int_t i) const {return fTracks[i];}
    // returns hit number kept in the array element i of fHits 
    virtual Int_t GetHit(Int_t i) const {return fHits[i];}
    // Return charge deposited by this track/hit
    virtual Float_t GetCharge(Int_t i) const {return fTcharges[i];}
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t,TArrayF &c);
    //copy the array trks[fgkSsdd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fgkSsdd;i++) fTracks[i]=trks[i];}
    //copy the array hits[fgkSsdd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fgkSsdd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    void SetSignalExpanded(Int_t sig){fSignalExpanded = sig;}
    virtual void Print(ostream *os); // Class ascii print function
    virtual void Read(istream *os);  // Class ascii read function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 protected:
    void InitObject(Float_t phys,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges);
    static const Int_t fgkSsdd = 10; // size of fTracks and fHits arrays
    
    // debugging  -- goes to the dictionary
    Int_t   fTracks[fgkSsdd];   //[fgkSsdd] tracks making this digit 
    Int_t   fHits[fgkSsdd];     //[fgkSsdd] hits associated to the tracks
                            // 3 hits temporarily - it will be only 1
    Float_t fTcharges[fgkSsdd];   //[fgkSsdd] charge per track making this digit 
    Float_t fPhysics;       // signal particles contribution to signal
    Int_t fSignalExpanded;  // 10 bit signal
    ClassDef(AliITSdigitSDD,3)   // Simulated digit object for SDD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSDD &source);
istream &operator>>(istream &os,AliITSdigitSDD &source);

#endif
