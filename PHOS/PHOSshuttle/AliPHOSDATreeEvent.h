#ifndef AliPHOSDATreeEvent_H
#define AliPHOSDATreeEvent_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --

#include "Rtypes.h"
#include "AliPHOSDATreeCluster.h"

class AliPHOSDATreeCluster;

class AliPHOSDATreeEvent{

  friend std::ostream& operator<<(std::ostream& out,const AliPHOSDATreeEvent& event);

 public:
  
  AliPHOSDATreeEvent(): fTime(0), fNDigits(0), fDigits(0), fNClusters(0), fClusters(0){/* */};
  AliPHOSDATreeEvent(const AliPHOSDATreeEvent& evt);
  AliPHOSDATreeEvent& operator=(const AliPHOSDATreeEvent& evt);
  virtual ~AliPHOSDATreeEvent();
  time_t GetTime() const{return fTime;};
  void SetTime(time_t time){fTime=time;};
  int GetNDigits() const{ return fNDigits; };
  int GetNClusters() const{ return fNClusters; };
  AliPHOSDATreeCluster& GetCluster(int nclusters){
    return fClusters[nclusters];
  };
  AliPHOSDATreeDigit& GetDigit(int ndigits){
    return fDigits[ndigits];
  };
  bool Fill(float fenergy,int row,int col);
  bool Fill(AliPHOSDATreeDigit& digit);
  bool ExecuteClustering();
  void Reset(){
    delete[] fDigits;
    delete[] fClusters;
    fTime = 0;
    fNDigits = 0;
    fNClusters = 0;
  };
  void Print(Option_t *option="") const;

 private:
  bool Clusterize(AliPHOSDATreeDigit& digit);

  time_t fTime;                     // Time information
  int fNDigits;                     // Number of digits in event
  AliPHOSDATreeDigit* fDigits;      //[fNDigits]
  int fNClusters;                   // Number of clusters in event
  AliPHOSDATreeCluster* fClusters;  //[fNClusters]

  ClassDef(AliPHOSDATreeEvent,1) // Simple Event Structure for PHOS DA
};
#endif

