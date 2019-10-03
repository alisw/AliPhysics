#ifndef AliXiStarPbPbEventCollection_H
#define AliXiStarPbPbEventCollection_H
//
// Class AliXiStarPbPbEventCollection, AliXiStarPbPbTrackStruct, AliXiStarPbPbEventStruct
//
// AliXiStarPbPbEventCollection, AliXiStarPbPbTrackStruct, AliXiStarPbPbEventStruct
// author:
//        Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//


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


class AliXiStarPbPbTrackStruct{
 public:

  AliXiStarPbPbTrackStruct();
  virtual ~AliXiStarPbPbTrackStruct();
  AliXiStarPbPbTrackStruct(const AliXiStarPbPbTrackStruct &obj); 
  AliXiStarPbPbTrackStruct &operator=(const AliXiStarPbPbTrackStruct &obj);

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

  ClassDef(AliXiStarPbPbTrackStruct, 1);
};

class AliXiStarPbPbEventStruct{
 public:

  AliXiStarPbPbEventStruct();
  virtual ~AliXiStarPbPbEventStruct();
  AliXiStarPbPbEventStruct(const AliXiStarPbPbEventStruct &obj); 
  AliXiStarPbPbEventStruct &operator=(const AliXiStarPbPbEventStruct &obj);

  Int_t fNTracks;// Events track count
  
  AliXiStarPbPbTrackStruct *fTracks;// Events track structure

  ClassDef(AliXiStarPbPbEventStruct, 1);
};

class AliXiStarPbPbEventCollection {
 public:
  
  AliXiStarPbPbEventCollection();
  AliXiStarPbPbEventCollection(Short_t);
  virtual ~AliXiStarPbPbEventCollection();
  AliXiStarPbPbEventCollection(const AliXiStarPbPbEventCollection &obj); 
  AliXiStarPbPbEventCollection &operator=(const AliXiStarPbPbEventCollection &obj);
  
  Short_t fFIFO; //Size of the Event Storage buffer. FIFO = first-in-first-out
  AliXiStarPbPbEventStruct *fEvtStr;// Event structure

  void FIFOShift();// remove event at end of buffer and add the new one
  void SetBuffSize(Short_t a){fFIFO = a;}// set size of event buffer (Nevents max to mix)
          
  ClassDef(AliXiStarPbPbEventCollection, 1);
};
#endif

















