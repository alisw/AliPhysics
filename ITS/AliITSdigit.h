#ifndef ALIITSDIGIT_H
#define ALIITSDIGIT_H

////////////////////////////////////////////////
//  Digits classes for set:ITS                //
////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>

//___________________________________________
class AliITSdigit: public TObject  {
  
public:
  
  Int_t fCoord1;  // Cell number on Z axis (SPD+SDD) , flag for side type (SSD)
  Int_t fCoord2 ; // Cell number on X axis (SPD+SDD) , strip number (SSD)
  Int_t fSignal;  // Signal in ADC counts
  
public:
  AliITSdigit() {
    // constructor
    fSignal=fCoord1=fCoord2=0;
  }
  AliITSdigit(Int_t *digits);
  virtual ~AliITSdigit() {
    // destructor
  }
  
  ClassDef(AliITSdigit,1)     // Real data digit object for set:ITS
    };

//___________________________________________
class AliITSdigitSPD: public AliITSdigit {
  
public:
  
    // debugging  -- goes to the dictionary
  Int_t       fTracks[3];         // tracks making this digit 
  Int_t       fHits[3];           // hits associated to the tracks 
                                  // 3 hits temporarily - it will be only 1
  
public:
  AliITSdigitSPD() {
    // constructor
    fSignal=fCoord1=fCoord2=0;
    fTracks[0]=fTracks[1]=fTracks[2]=-3;
    fHits[0]=fHits[1]=fHits[2]=-1;
  }
  
  AliITSdigitSPD(Int_t *digits);
  AliITSdigitSPD(Int_t *digits, Int_t *tracks, Int_t *hits);
  
  virtual ~AliITSdigitSPD(){
    // destructor
  }
  virtual int *GetTracks() {
    // returns pointer to the array of tracks which make this digit
    return &fTracks[0];
  }
  
    ClassDef(AliITSdigitSPD,1)   // Simulated digit object for SPD
      };

//___________________________________________
class AliITSdigitSDD: public AliITSdigit {
  
public:
  
  // debugging  -- goes to the dictionary
  Int_t       fTracks[3];         // tracks making this digit 
  Int_t       fHits[3];           // hits associated to the tracks
                                  // 3 hits temporarily - it will be only 1
  Float_t     fTcharges[3];       // charge per track making this digit 
  Float_t     fPhysics;           // signal particles contribution to signal
  
public:
  AliITSdigitSDD() {
    // constructor
    fSignal=fCoord1=fCoord2=0;
    fTracks[0]=fTracks[1]=fTracks[2]=-3;
    fHits[0]=fHits[1]=fHits[2]=-1;
    fPhysics=0; fTcharges[0]=fTcharges[1]=fTcharges[2]=0;
  }
  
  AliITSdigitSDD(Float_t phys,Int_t *digits);
  AliITSdigitSDD( Float_t phys, Int_t *digits, Int_t *tracks, Int_t *hits, Float_t *charges);
  
  virtual ~AliITSdigitSDD(){
    // destructor
  }
  virtual int *GetTracks() {
    // returns pointer to the array of tracks which make this digit
    return &fTracks[0];
  }
  
  ClassDef(AliITSdigitSDD,1)   // Simulated digit object for SDD
    };

//_____________________________________________________________________________

class AliITSTransientDigit : public AliITSdigitSDD {
public:
  TObjArray     *fTrackList;  // track list 
public:
  AliITSTransientDigit() {
    // constructor
    fTrackList=0;
  }
  AliITSTransientDigit(Float_t phys,Int_t *digits);
  virtual ~AliITSTransientDigit() {
    // destructor
    delete fTrackList;
  }
  AliITSTransientDigit(const AliITSTransientDigit &source); // copy constructor
  AliITSTransientDigit& operator=(const AliITSTransientDigit &source); // ass. operator
  TObjArray  *TrackList()   {
    // returns pointer to the TObjArray of tracks and associated charges
    return fTrackList;
  }
  
  ClassDef(AliITSTransientDigit,1)  // Transient digit for set: ITS
    };

//___________________________________________
class AliITSdigitSSD: public AliITSdigit {
  
public:
  
  // debugging  -- goes to the dictionary
  Int_t       fTracks[3];         // tracks making this digit 
  Int_t       fHits[3];           // hits associated to the tracks
                                  // 3 hits temporarily - it will be only 1
  
public:
  AliITSdigitSSD() {
    // constructor
    fSignal=fCoord1=fCoord2=0;
    fTracks[0]=fTracks[1]=fTracks[2]=-3;
    fHits[0]=fHits[1]=fHits[2]=-1;
  }
  
  AliITSdigitSSD(Int_t *digits);
  AliITSdigitSSD(Int_t *digits, Int_t *tracks, Int_t *hits);
  
  virtual ~AliITSdigitSSD(){
    // destructor
  }
  
  Int_t  GetSignal() const {
    // returns signal
    return fSignal;
  }

  Int_t  GetStripNumber() const {
    // returns strip number
    return fCoord2;
  }
  
  Int_t  IsSideP() const {
    //returns 1  when side P and 0 when side N     
    return fCoord1;
  }
  
  virtual int *GetTracks() {
    // returns pointer to the array of tracks which make this digit
    return &fTracks[0];
    }
  
  ClassDef(AliITSdigitSSD,1)   // Simulated digit object for SSD
    };


#endif
