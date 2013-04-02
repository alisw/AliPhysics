#ifndef ALICHAOTICITYEVENTCOLLECTION
#define ALICHAOTICITYEVENTCOLLECTION

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


class AliChaoticityTrackStruct{// detector tracks
  
 public:
  AliChaoticityTrackStruct();
  virtual ~AliChaoticityTrackStruct();
  AliChaoticityTrackStruct(const AliChaoticityTrackStruct &obj); 
  AliChaoticityTrackStruct &operator=(const AliChaoticityTrackStruct &obj);


  UInt_t fStatus;
  UInt_t fFiltermap;
  Int_t fId;
  Double_t fPhi;
  Float_t fPt;
  Float_t fMom;
  Double_t fP[3];
  Int_t fCharge;
  Float_t fEta;
  Float_t fMass;
  Float_t fDCAXY;
  Float_t fDCAZ;
  Float_t fDCA;
  Float_t fEaccepted;
  Short_t fKey;
  TBits fClusterMap;
  TBits fSharedMap;
  Double_t fX[3];
  Bool_t fTOFhit;
  Bool_t fElectron;
  Bool_t fMuon;
  Bool_t fPion;
  Bool_t fKaon;
  Bool_t fProton;
  Int_t fLabel;// MC

  ClassDef(AliChaoticityTrackStruct, 1);
};

class AliChaoticityPairStruct{// low Qinv pairs

 public:
  AliChaoticityPairStruct();
  virtual ~AliChaoticityPairStruct();
  AliChaoticityPairStruct(const AliChaoticityPairStruct &obj); 
  AliChaoticityPairStruct &operator=(const AliChaoticityPairStruct &obj);

  Float_t fP1[3];
  Float_t fP2[3];
  Float_t fE1;
  Float_t fE2;
  Short_t fCharge1;
  Short_t fCharge2;
  Int_t fIndex1;
  Int_t fIndex2;
  Float_t fQinv;
  Short_t fKey1;
  Short_t fKey2;
  Int_t fLabel1;
  Int_t fLabel2;
  Float_t fP1MC[3];
  Float_t fP2MC[3];

  ClassDef(AliChaoticityPairStruct, 1);
};

class AliChaoticityNormPairStruct{// Norm Qinv pairs

 public:
  AliChaoticityNormPairStruct();
  virtual ~AliChaoticityNormPairStruct();
  AliChaoticityNormPairStruct(const AliChaoticityNormPairStruct &obj); 
  AliChaoticityNormPairStruct &operator=(const AliChaoticityNormPairStruct &obj);

  Short_t fCharge1;
  Short_t fCharge2;
  Int_t fIndex1;
  Int_t fIndex2;
  Short_t fKey1;
  Short_t fKey2;

  ClassDef(AliChaoticityNormPairStruct, 1);
};

class AliChaoticityMCStruct{// MC info

 public:
  AliChaoticityMCStruct();
  virtual ~AliChaoticityMCStruct();
  AliChaoticityMCStruct(const AliChaoticityMCStruct &obj); 
  AliChaoticityMCStruct &operator=(const AliChaoticityMCStruct &obj);

  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Float_t fPtot;

  ClassDef(AliChaoticityMCStruct, 1);
};

class AliChaoticityEventStruct{// like particle_event
  
 public:
  AliChaoticityEventStruct();
  virtual ~AliChaoticityEventStruct();
  AliChaoticityEventStruct(const AliChaoticityEventStruct &obj); 
  AliChaoticityEventStruct &operator=(const AliChaoticityEventStruct &obj);


  Int_t fFillStatus;
  Int_t fNtracks;
  Int_t fNpairsSE;
  Int_t fNpairsME;
  Int_t fMCarraySize;
  AliChaoticityTrackStruct *fTracks;
  AliChaoticityPairStruct *fPairsSE;
  AliChaoticityPairStruct *fPairsME;
  AliChaoticityMCStruct *fMCtracks;

  ClassDef(AliChaoticityEventStruct, 1);
};



class AliChaoticityEventCollection {
  
  public:
    AliChaoticityEventCollection();
    AliChaoticityEventCollection(Short_t,Int_t,Int_t,Int_t,Bool_t);
    virtual ~AliChaoticityEventCollection();
    AliChaoticityEventCollection(const AliChaoticityEventCollection &obj); 
    AliChaoticityEventCollection &operator=(const AliChaoticityEventCollection &obj);
   
    void FIFOShift();
    void SetBuffSize(Short_t a){fFIFO = a;}
 
    Short_t fFIFO; //Size of the Event Storage buffer.
    Int_t fLimit; //Max number of tracks
    Int_t fPairLimit; //Max number of lowQ pairs
    Int_t fMCLimit; //Max number of MC tracks
    AliChaoticityEventStruct *fEvtStr;

   ClassDef(AliChaoticityEventCollection, 1);
};
#endif
