#ifndef ALIITSDIGIT_H
#define ALIITSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
    AliITSdigit() {//default constructor. zero all values.
	fSignal=fCoord1=fCoord2=0;}
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
    virtual void SetCoord1(Int_t i){fCoord1 = i;} // Sets fCoord1 value
    virtual void SetCoord2(Int_t i){fCoord2 = i;} // Sets fCoord12value
    virtual void SetSignal(Int_t i){fSignal = i;} // Sets fSignal value
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function

 public:
    Int_t fCoord1; // Cell number on Z axis (SPD+SDD), flag for side type (SSD)
    Int_t fCoord2; // Cell number on X axis (SPD+SDD), strip number (SSD)
    Int_t fSignal; // Signal in ADC counts

    ClassDef(AliITSdigit,1)     // Real data digit object for set:ITS

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigit &source);
istream &operator>>(istream &os,AliITSdigit &source);
//______________________________________________________________________
class AliITSdigitSPD: public AliITSdigit {

 public:
    AliITSdigitSPD(); //default creator
    AliITSdigitSPD(const Int_t *digits);//standard creator digits only
    //standard creator with digits, tracks, and hits
    AliITSdigitSPD(const Int_t *digits,const Int_t *tracks,const Int_t *hits);
    virtual ~AliITSdigitSPD(){/*destructor*/}
    // returns the array size used to store Tracks and Hits
    static Int_t GetNTracks() {return fkSspd;}
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
    //copy the array trks[fkSspd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fkSspd;i++) fTracks[i]=trks[i];}
    //copy the array hits[fkSspd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fkSspd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function

 private:
    static const Int_t fkSspd = 10; // size of fTracks and fHits arrays
    
 public:  
    // debugging  -- goes to the dictionary
    Int_t fTracks[fkSspd]; //[fkSspd] tracks making this digit 
    Int_t fHits[fkSspd];   //[fkSspd] hits associated to the tracks
    Int_t fSignalSPD;   // Signal in electrons

    ClassDef(AliITSdigitSPD,2)   // Simulated digit object for SPD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSPD &source);
istream &operator>>(istream &os,AliITSdigitSPD &source);

//______________________________________________________________________
class AliITSdigitSDD: public AliITSdigit {

 public:
    AliITSdigitSDD(); //default creator
    //standard creator  with digits and "phys"
    AliITSdigitSDD(Float_t phys,const Int_t *digits);
    //standard creator with digits, tracls, hits, "phys", and charge
    AliITSdigitSDD( Float_t phys,const Int_t *digits,const Int_t *tracks,
		    const Int_t *hits,const Float_t *charges);
    virtual ~AliITSdigitSDD(){/* destructor*/}
    // returns the array size used to store Tracks and Hits
    static Int_t GetNTracks() {return fkSsdd;}
    // returns pointer to the array of tracks which make this digit
    virtual Int_t *GetTracks() {return &fTracks[0];}
    // returns the pointer to the array of hits which made this digit
    virtual Int_t *GetHits() {return &fHits[0];}
    // returns track number kept in the array element i of fTracks 
    virtual Int_t GetTrack(Int_t i) const {return fTracks[i];}
    // returns hit number kept in the array element i of fHits 
    virtual Int_t GetHit(Int_t i) const {return fHits[i];}
    // Return charge deposited by this track/hit
    virtual Float_t GetCharge(Int_t i){return fTcharges[i];}
    // returns TArrayI of unduplicated track numbers (summed over hits).
    virtual Int_t GetListOfTracks(TArrayI &t,TArrayF &c);
    //copy the array trks[fkSsdd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fkSsdd;i++) fTracks[i]=trks[i];}
    //copy the array hits[fkSsdd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fkSsdd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function

 private:
    static const Int_t fkSsdd = 10; // size of fTracks and fHits arrays
    
 public:
    // debugging  -- goes to the dictionary
    Int_t   fTracks[fkSsdd];   //[fkSsdd] tracks making this digit 
    Int_t   fHits[fkSsdd];     //[fkSsdd] hits associated to the tracks
                            // 3 hits temporarily - it will be only 1
    Float_t fTcharges[fkSsdd];   //[fkSsdd] charge per track making this digit 
    Float_t fPhysics;       // signal particles contribution to signal

    ClassDef(AliITSdigitSDD,2)   // Simulated digit object for SDD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSDD &source);
istream &operator>>(istream &os,AliITSdigitSDD &source);
//______________________________________________________________________
class AliITSTransientDigit : public AliITSdigitSDD {

 public:
    AliITSTransientDigit() {/*default constructor*/fTrackList=0;}
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

 public:
    TObjArray *fTrackList;  // track list 

    ClassDef(AliITSTransientDigit,1)  // Transient digit for set: ITS

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTransientDigit &source);
istream &operator>>(istream &os,AliITSTransientDigit &source);
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
    static Int_t GetNTracks() {return fkSssd;}
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
    //copy the array trks[fkSssd] into fTracks
    virtual void SetTracks(const Int_t *trks){
	for(Int_t i=0;i<fkSssd;i++) fTracks[i]=trks[i];}
    //copy the array hits[fkSssd] into fHits
    virtual void SetHits(const Int_t *hits){
	for(Int_t i=0;i<fkSssd;i++) fHits[i]=hits[i];}
    //set array element i of fTracks to trk.
    virtual void SetTrack(Int_t i,Int_t trk){fTracks[i]=trk;}
    //set array element i of fHits to hit.
    virtual void SetHit(Int_t i,Int_t hit){fHits[i]=hit;}
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function

 private:
    static const Int_t fkSssd = 10; // size of fTracks and fHits arrays
    
 public:
    // debugging  -- goes to the dictionary
    Int_t fTracks[fkSssd]; //[fkSssd] tracks making this digit 
    Int_t fHits[fkSssd];   //[fkSssd] hits associated to the tracks
                        // 3 hits temporarily - it will be only 1
    
    ClassDef(AliITSdigitSSD,2)   // Simulated digit object for SSD

};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSdigitSSD &source);
istream &operator>>(istream &os,AliITSdigitSSD &source);

#endif
