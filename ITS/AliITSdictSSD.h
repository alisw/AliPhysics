#ifndef ALIITSDICTSSD_H
#define ALIITSDICTSSD_H

#include <TArrayF.h>

class AliITSdictSSD  {

public:
  
  AliITSdictSSD() {
    // constructor
    fTracks = 0;
  }
  ~AliITSdictSSD() {
    // destructor
  };
  
  void   AddTrack(Int_t track); 
  Int_t  GetTrack(Int_t index);
  Int_t  GetNTracks() {
    // get num of tracks
    return fTracks;
  }
  void   ZeroTracks() { 
    // zero tracks
    for (Int_t i =0;i<10;i++) fTrack[i]=0; fTracks = 0;
  }
  
private:           
  Int_t fTrack[10];    // Track array
  Int_t fTracks;       // Tracks
};   


#endif
