#ifndef ALIXISTARPP13TEVEVENTCOLLECTION_H
#define ALIXISTARPP13TEVEVENTCOLLECTION_H
//
// Class AliXiStarpp13TeVEventCollection, AliXiStarpp13TeVTrackStruct, AliXiStarpp13TeVEventStruct
//
// AliXiStarpp13TeVEventCollection, AliXiStarpp13TeVTrackStruct, AliXiStarpp13TeVEventStruct
// author:
//  (Original Code) Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//  (1st Modification) Jihye Song (jihye.song@cern.ch)
//  (last Modification) Bong-Hwi Lim (bong-hwi.lim@cern.ch)

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"

using namespace std;


class AliXiStarpp13TeVTrackStruct{
 public:

  AliXiStarpp13TeVTrackStruct();
  virtual ~AliXiStarpp13TeVTrackStruct();
  AliXiStarpp13TeVTrackStruct(const AliXiStarpp13TeVTrackStruct &obj ); 
  AliXiStarpp13TeVTrackStruct &operator=(const AliXiStarpp13TeVTrackStruct &obj );

  UInt_t fStatus;// track status
  UInt_t fFilterMap;// filter map for AOD filterbits
  Int_t fID;// track id
  Double_t fPhi;// track phi angle
  Double_t fPt;// track pt
  Float_t fMom;// track full momentum
  Double_t fP[3];// track 3d momentum
  Int_t fCharge;// track charge
  Double_t fEta;// track eta
  Double_t fMass;// track accepted mass
  Double_t fDCAXY;// track dca to PV in xy
  Double_t fDCAZ;// track dca to PV in z
  Double_t fDCA;// track full dca
  Double_t fX[3];// track x position
  Double_t fCov[21];// track covariance matrix
  Float_t fNSigmaPi;// track Nsigma pion
  Float_t fNSigmaK;// track Nsigma kaon
  Float_t fNSigmaPr;// track Nsigma proton
  Int_t fLabel;// track label for MC studies
  UShort_t fNclusTPC;// TPC N clusters

  ClassDef(AliXiStarpp13TeVTrackStruct, 1);
};

class AliXiStarpp13TeVEventStruct{
 public:

  AliXiStarpp13TeVEventStruct();
  virtual ~AliXiStarpp13TeVEventStruct();
  AliXiStarpp13TeVEventStruct(const AliXiStarpp13TeVEventStruct &obj ); 
  AliXiStarpp13TeVEventStruct &operator=(const AliXiStarpp13TeVEventStruct &obj );

  Int_t fNTracks;// Events track count
  AliXiStarpp13TeVTrackStruct *fTracks;// Events track structure

  ClassDef(AliXiStarpp13TeVEventStruct, 1);
};

class AliXiStarpp13TeVEventCollection {
 public:
  
  AliXiStarpp13TeVEventCollection();
  AliXiStarpp13TeVEventCollection(Short_t);
  virtual ~AliXiStarpp13TeVEventCollection();
  AliXiStarpp13TeVEventCollection(const AliXiStarpp13TeVEventCollection &obj ); 
  AliXiStarpp13TeVEventCollection &operator=(const AliXiStarpp13TeVEventCollection &obj );
  
  Short_t fFIFO; //Size of the Event Storage buffer. FIFO = first-in-first-out
  AliXiStarpp13TeVEventStruct *fEvtStr;// Event structure

  void FIFOShift();// remove event at end of buffer and add the new one
  void SetBuffSize(Short_t a){fFIFO = a;}// set size of event buffer (Nevents max to mix)
          
  ClassDef(AliXiStarpp13TeVEventCollection, 1);
};
#endif

















